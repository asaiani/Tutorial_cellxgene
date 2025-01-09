import os
import wandb
import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
import numpy as np
import anndata
import sklearn.decomposition
#import seaborn as sns
import cytobench
import VAE
import score_VAE
from scipy.sparse import issparse
import json
import wandb
import census_function as cf
wandb.login(relogin=True, key = '7e568e5cae11c7b5987dfa18bf33457bc5a9f276') 
# I need to do this, because the terminal does not work properly


# Configure default plot style
#sns.set(style='whitegrid')
figsize = 6  # Default figure size
# --- 2. Dataset Loading and Preprocessing ---
# Load dataset
dataset_name = 'complete'
# read the dictonary 
import json
# Load the JSON file in which is stored the dictonary of the dataset
with open("../filtered_dataset_id_donor_id_mus_musculus.json", "r") as json_file:
    dataset_to_donor_map = json.load(json_file)

from itertools import islice

# Get the first two items of the dictionary
first_two_items = dict(islice(dataset_to_donor_map.items(), 2))

print(first_two_items)

import matplotlib.pyplot as plt
import pandas as pd




# Let's aggregate the dataset that has been downloaded and filtered to keep only the first 256 highly variance genes
dataset = cf.aggregate_across_datasets(dataset_to_donor_map, organism="Mus musculus", main_directory="dataset_small_split")
# Ensure X is always a dense array when needed
if issparse(dataset.X):
    dataset.X = dataset.X.toarray()

# Logarithmic transformation
import scanpy as sc
sc.pp.log1p(dataset)

# Extract count matrix and metadata
X = dataset.X
metadata = dataset.obs
var = dataset.var

dataset_train = dataset[dataset.obs['split']== 0] # train
dataset_val = dataset[dataset.obs['split']== 1] # validation
dataset_test = dataset[dataset.obs['split']== 2] # test
# keep in mind you need to normalize the dataset... 

# Assuming `dataset.obs` is a DataFrame-like object
'''# Count the number of samples for each donor_id
donor_sample_counts = dataset.obs.groupby('donor_id').size()

# Define a threshold for "small" donor IDs
threshold = 100  # Adjust this as needed

# Filter for donor IDs with fewer than `threshold` samples
small_donor_sample_counts = donor_sample_counts[donor_sample_counts < threshold]
# Plot the histogram for small donor IDs
plt.figure(figsize=(10, 6))
plt.hist(small_donor_sample_counts, bins=min(threshold, len(small_donor_sample_counts)), edgecolor='k', alpha=0.7)
plt.xlabel('Number of Samples per Donor ID')
plt.ylabel('Number of Donor IDs')
plt.title(f'Distribution of Small Samples (Donor IDs < {threshold} samples)')
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.savefig('small_distribution.png')
plt.savefig('distribution.png')'''




X_train = dataset_train.X
X_val = dataset_val.X
X_test = dataset_test.X

np.isnan(X)

sweep_config = {
    'method': 'random',
    'parameters': {
        'n_hidden_layers': {'values': [0, 1, 2, 3, 4, 5]},
        'bottleneck_units': {'values': [8, 16, 32, 64, 128, 256]},
        'kl_loss_factor': {'values': [ 1e-5, 1e-6]},
        'max_training_epochs': {'values': [1]},
        'batch_size':{'values' : [128]}, 
        'training': {
            'values': [
                {'learning_rate': 1e-3, 'learning_rate_decay': None},
                {'learning_rate': 1e-4, 'learning_rate_decay': None},
                {'learning_rate': 1e-3, 'learning_rate_decay': 'Exponential'}
                ]
                },
        'alpha': {'values': [1]} 
    }
}

# Initialize the sweep
sweep_id = wandb.sweep(sweep=sweep_config, project='PROVE_first_2_datasets' + dataset_name )
print(f"Sweep initialized with ID: {sweep_id}")


def train():
    wandb.init()  # Initialize a WandB run

    # Retrieve and debug parameters
    config = wandb.config
    print("Config received in train function:", dict(config))

    # Calculate `plateau_epochs` as a fraction of `max_training_epochs`
    plateau_fraction = 0.1  # Change this fraction as needed
    plateau_epochs = int(config['max_training_epochs'] * plateau_fraction)
    training_params = config['training']

    # Train model
    try:
        model = VAE.train_vae(
            X_train,
            X_val, 
            X_test,
            n_hidden_layers=config['n_hidden_layers'],
            bottleneck_units=config['bottleneck_units'],
            kl_loss_factor=config['kl_loss_factor'], 
            plateau_epochs=plateau_epochs,
            max_training_epochs=config['max_training_epochs'],
            batch_size=config['batch_size'], 
            learning_rate=training_params['learning_rate'], 
            learning_rate_decay=training_params['learning_rate_decay'], 
            alpha=config['alpha']
        )
        # Check if training returned None due to exploding/infinite loss
        if model is None:
            print("Training failed: Loss went to infinity or became NaN.")
            wandb.log({"Training Status": "Failed - Infinite Loss", "epoch": None})
            return
    except KeyError as e:
        print(f"KeyError: {e}. Check if the sweep config has this parameter.")
        return

    save_dir_coverage = './coverage_estimators/' + dataset_name

    # Prepare results dictionary
    results = {
        "train_mean_no_sampling_score": 0,
        "train_mean_sampling_score": 0,
        "val_mean_no_sampling_score": 0,
        "val_mean_sampling_score": 0,
        "test_mean_no_sampling_score": 0,
        "test_mean_sampling_score": 0
    }

    with open("splits.json", "r") as json_file:
        splits = json.load(json_file)

        # Minimum samples required to perform clustering
    MIN_SAMPLES = 1000  # or set it to the value of `knn`

    # Filter donor IDs with insufficient samples
    filtered_splits = {}
    for split_name, split_data in splits.items():
        if split_name == "categories":  # Skip the categories key
            continue

        filtered_splits[split_name] = {}
        for dataset_id, donor_ids in split_data.items():
            # Check the number of samples for each donor ID
            valid_donor_ids = [
                donor_id for donor_id in donor_ids
                if len(dataset.obs[
                    (dataset.obs['dataset_id'] == dataset_id) &
                    (dataset.obs['donor_id'] == donor_id)
                ]) >= MIN_SAMPLES
            ]
            if valid_donor_ids:
                filtered_splits[split_name][dataset_id] = valid_donor_ids

    # Replace splits with the filtered version
    splits = filtered_splits


    # Iterate over each split (train, val, test)
    for split_name, split_data in splits.items():
        if split_name == "categories":  # Skip the "categories" key
            continue

        print(f"Processing split: {split_name} ...")
        no_sampling_scores = []
        sampling_scores = []

        # Iterate over dataset IDs and their donor IDs
        for dataset_id, donor_ids in split_data.items():
            for donor_id in donor_ids:
                # Score for the current donor ID
                print(f"Scoring for donor ID: {donor_id} in {split_name} ...")
                score_no_sampling, score_with_sampling = score_VAE.score_autoencoder(
                    model,
                    metadata,
                    X,
                    column_name="donor_id",  # The column to split by
                    column_value=donor_id,
                    all_dataset=False,
                    save_dir=save_dir_coverage
                )
                # Append scores to the respective lists
                no_sampling_scores.append(score_no_sampling)
                sampling_scores.append(score_with_sampling)

        # Calculate mean scores for the current split
        results[f"{split_name}_mean_no_sampling_score"] = sum(no_sampling_scores) / len(no_sampling_scores)
        results[f"{split_name}_mean_sampling_score"] = sum(sampling_scores) / len(sampling_scores)

        print(
            f"Mean scores for {split_name}: "
            f"No sampling - {results[f'{split_name}_mean_no_sampling_score']}, "
            f"With sampling - {results[f'{split_name}_mean_sampling_score']}"
        )

    # Log results to WandB
    wandb.log(results)
    wandb.finish()  # End the WandB run


# Run the sweep
wandb.agent(sweep_id, train)
