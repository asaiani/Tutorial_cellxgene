import os
import json
import cellxgene_census
import numpy as np 
import pandas as pd
import anndata

# This function is to retrieve the donor IDs and can be used by itself 
def donor_id_information(dataset_id, organism, display_info=False):
    """
    Retrieves all available donor IDs for the specified dataset from the census.
    
    Parameters:
    - dataset_id (str): The ID of the dataset to retrieve information for.
    - organism can be 'Mus Musculus' or 'Homo Sapiens'
    - display_info (bool): Whether to display the dataset information as a DataFrame (default is False).
    
    Returns:
    - metadata_dict (dict): Dictionary containing dataset metadata including available donor IDs.
    """
    # Open the census and retrieve metadata
    census = cellxgene_census.open_soma()
    metadata = cellxgene_census.get_obs(census, organism = organism , value_filter=f"dataset_id == '{dataset_id}'")
    available_donor_ids = set(np.unique(metadata['donor_id']))

    # Read the census datasets table to get additional dataset information
    census_datasets = census["census_info"]["datasets"].read().concat().to_pandas()
    census_datasets = census_datasets.set_index("soma_joinid")

    # Filter for the specific dataset_id
    dataset_info = pd.DataFrame(census_datasets[census_datasets.dataset_id == dataset_id])
    if dataset_info.empty:
        raise ValueError(f"Dataset ID {dataset_id} not found in CellxGene Census.")
    else:
        if display_info:
            display(dataset_info)

        # Build the metadata dictionary
        metadata_dict = {
            "dataset_id": dataset_id,
            "available_donor_ids": list(available_donor_ids),  # All IDs from census
            "downloaded_donor_ids": [],  # To track IDs that have been downloaded
            "dataset_info": dataset_info.to_dict(orient="records")  # Additional dataset info
        }
        
        return metadata_dict


# These two function takes as input the dataset_id and the donor_id to retrieve
# The census is accessed only if necessary, if the file are already present in the directory then the census is not invoked 
# The file are stored in slices (one file for each donor_id) to garantee a low RAM necessities 
def update_metadata_file(dataset_id, donor_id, organism, main_directory):
    """
    Updates the metadata file for the specified dataset_id by appending the queried donor_id to downloaded_donor_ids.

    Parameters:
    - dataset_id (str): The ID of the dataset to update.
    - donor_id (str): The donor ID to append to downloaded_donor_ids.
    - organism can be 'Mus Musculus' or 'Homo Sapiens'
    """
    # Define the directory and path for the metadata file
    metadata_directory = os.path.join(main_directory, dataset_id)
    metadata_file_path = os.path.join(metadata_directory, f"{dataset_id}_metadata.json")
    
    # Load or create metadata structure
    if os.path.exists(metadata_file_path):
        # Load existing metadata
        with open(metadata_file_path, "r") as metadata_file:
            metadata_dict = json.load(metadata_file)
    else:
        # If file does not exist, retrieve available donor IDs and create metadata structure
        metadata_dict = donor_id_information(dataset_id, organism = organism)
        os.makedirs(metadata_directory, exist_ok=True)  # Ensure the directory exists

    # Ensure donor_id is not already in downloaded_donor_ids before appending
    if donor_id not in metadata_dict["downloaded_donor_ids"]:
        metadata_dict["downloaded_donor_ids"].append(donor_id)

        # Save the updated metadata back to the file
        with open(metadata_file_path, "w") as metadata_file:
            json.dump(metadata_dict, metadata_file, indent=4)

        print(f"Updated metadata for dataset_id '{dataset_id}': added donor_id '{donor_id}' to downloaded_donor_ids.")
    else:
        print(f"donor_id '{donor_id}' already exists in downloaded_donor_ids for dataset_id '{dataset_id}'.")



def save_adata_slices(dataset_id, donor_ids, organism, main_directory):
    """
    Retrieves and saves AnnData slices for each donor_id and dataset_id combination.

    Parameters:
    - donor_ids (list): List of donor IDs to retrieve data for.
    - dataset_ids (list): List of dataset IDs to retrieve data for.
    - organism (str): Organism name to use in the census query ('Mus musculus' or 'Homo Sapiens').
    """
    #main_directory = "my_data/adata_slices"
    os.makedirs(main_directory, exist_ok=True)

    
    dataset_directory = os.path.join(main_directory, dataset_id)
    metadata_file_path = os.path.join(dataset_directory, f"{dataset_id}_metadata.json")

    # Load metadata if it exists, or initialize by calling donor_id_information
    if os.path.exists(metadata_file_path):
        with open(metadata_file_path, "r") as metadata_file:
            metadata_dict = json.load(metadata_file)
        available_donor_ids = set(metadata_dict["available_donor_ids"])
        downloaded_donor_ids = set(metadata_dict["downloaded_donor_ids"])
    else:
        # Retrieve metadata and create metadata file if it doesn't exist
        metadata_dict = donor_id_information(dataset_id, organism)
        available_donor_ids = set(metadata_dict["available_donor_ids"])
        downloaded_donor_ids = set(metadata_dict["downloaded_donor_ids"])
        os.makedirs(dataset_directory, exist_ok=True)
        with open(metadata_file_path, "w") as metadata_file:
            json.dump(metadata_dict, metadata_file, indent=4)

    # Identify which donor_ids need to be downloaded
    donor_ids_to_download = set(donor_ids) - downloaded_donor_ids
    if not donor_ids_to_download:
        print(f"All requested donor IDs have already been downloaded for dataset_id '{dataset_id}'.")
    else: 
        print(f"Downloading {donor_ids_to_download} for dataset_id '{dataset_id}'.")
        # Open the census only if there are donor IDs to download
        census = cellxgene_census.open_soma()

    for donor_id in donor_ids_to_download:
        obs_value_filter = f"dataset_id == '{dataset_id}' and donor_id == '{donor_id}'"
        try:
            adata_slice = cellxgene_census.get_anndata(
                census=census,
                organism=organism,
                obs_value_filter=obs_value_filter
            )
            
            # Check if the slice is not empty
            if adata_slice.n_obs > 0:
                # Sanitize the donor_id to replace problematic characters
                sanitized_donor_id = donor_id.replace("/", "_")
                file_name = f"{dataset_id}_{sanitized_donor_id}.h5ad"
                file_path = os.path.join(dataset_directory, file_name)
                adata_slice.write(file_path)
                print(f"Saved AnnData slice for dataset_id '{dataset_id}' and donor_id '{donor_id}' in '{file_path}'")

                # Update metadata to include the newly downloaded donor_id
                update_metadata_file(dataset_id, donor_id, organism, main_directory)

            else:
                print(f"No data found for dataset_id '{dataset_id}' and donor_id '{donor_id}'")

        except Exception as e:
            print(f"Failed to retrieve data for dataset_id '{dataset_id}' and donor_id '{donor_id}': {e}")

    print(f"All AnnData slices have been saved in directory: {main_directory}")

'''
# This function instead serve to upload in memory the aggregated anndata objects
# Takes as input donor_ids, dataset_id and organism, gives as an output an anndata aggregated object
# This function calles save_adata_slices to prevent any missing file and concats the file with the specified donor ID 
def aggregate_donor_id(dataset_id, donor_ids, organism="Mus musculus", main_directory = "my_data/adata_slices"):
    """
    Aggregates the AnnData files for specified donor_ids and dataset_ids into a single AnnData object.

    Parameters:
    - donor_ids (list): List of donor IDs to include in the aggregation.
    - dataset_ids (list): List of dataset IDs to include in the aggregation.
    - organism (str): Organism name to use in the census query (default is 'Mus musculus').

    Returns:
    - AnnData: A single AnnData object containing all specified donor_id and dataset_id data.
    """
    # Ensure required data slices are downloaded
    save_adata_slices(dataset_id, donor_ids, organism, main_directory)

    # Initialize an empty list to collect AnnData slices
    adata_list = []

    # Iterate over each dataset_id and donor_id to load the saved AnnData slices
    
    dataset_directory = os.path.join(main_directory, dataset_id)
    
    for donor_id in donor_ids:
        # Construct the file path for the specific donor_id and dataset_id
        file_name = f"{dataset_id}_{donor_id}.h5ad"
        file_path = os.path.join(dataset_directory, file_name)

        # Check if the file exists and load it if it does
        if os.path.exists(file_path):
            # Load the AnnData slice and add to list
            adata_slice = anndata.read_h5ad(file_path)
            adata_list.append(adata_slice)
            print(f"Loaded data for dataset_id '{dataset_id}' and donor_id '{donor_id}' from '{file_path}'")
        else:
            print(f"Warning: Expected file '{file_path}' does not exist. Skipping this donor_id.")

    # Concatenate all the loaded AnnData objects if there are any
    if adata_list:
        aggregated_adata = anndata.concat(adata_list, join="outer", label="donor_id", index_unique="-")
        # grab all var DataFrames from our dictionary
        all_var = [x.var for x in adata_list]
        # concatenate them
        all_var = pd.concat(all_var, join="outer")
        # remove duplicates
        all_var = all_var[~all_var.index.duplicated()]
        # put all together
        aggregated_adata.var = all_var.loc[aggregated_adata.var_names]
        # This is for preventing a warning
        #aggregated_adata.obs_names_make_unique()
        print(f"Aggregated AnnData object created with {aggregated_adata.n_obs} observations.")
        return aggregated_adata
    else:
        print("No AnnData files were found to aggregate.")
        return None
'''

def aggregate_donor_id(dataset_id, donor_ids, organism="Mus musculus", main_directory="my_data/adata_slices"):
    """
    Aggregates the AnnData files for specified donor_ids and dataset_ids into a single AnnData object.

    Parameters:
    - donor_ids (list): List of donor IDs to include in the aggregation.
    - dataset_ids (list): List of dataset IDs to include in the aggregation.
    - organism (str): Organism name to use in the census query (default is 'Mus musculus').

    Returns:
    - AnnData: A single AnnData object containing all specified donor_id and dataset_id data.
    """
    # Ensure required data slices are downloaded
    save_adata_slices(dataset_id, donor_ids, organism, main_directory)

    # Initialize an empty list to collect AnnData slices
    adata_list = []

    # Iterate over each dataset_id and donor_id to load the saved AnnData slices
    dataset_directory = os.path.join(main_directory, dataset_id)
    
    for donor_id in donor_ids:
        # Construct the file path for the specific donor_id and dataset_id
        sanitized_donor_id = donor_id.replace("/", "_")
        file_name = f"{dataset_id}_{sanitized_donor_id}.h5ad"
        file_path = os.path.join(dataset_directory, file_name)

        # Check if the file exists and load it if it does
        if os.path.exists(file_path):
            # Load the AnnData slice and set donor_id as a new column in .obs
            adata_slice = anndata.read_h5ad(file_path)
            adata_slice.obs["donor_id"] = donor_id
            adata_list.append(adata_slice)
            print(f"Loaded data for dataset_id '{dataset_id}' and donor_id '{donor_id}' from '{file_path}'")
        else:
            print(f"Warning: Expected file '{file_path}' does not exist. Skipping this donor_id.")

    # Concatenate all the loaded AnnData objects if there are any
    if adata_list:
        aggregated_adata = anndata.concat(
            adata_list, 
            join="outer", 
            label="dataset_id", 
            keys=donor_ids,  # Use donor IDs as keys
            index_unique="-"
        )
        # Grab all var DataFrames from our list of AnnData objects
        all_var = [x.var for x in adata_list]
        # Concatenate them
        all_var = pd.concat(all_var, join="outer")
        # Remove duplicates
        all_var = all_var[~all_var.index.duplicated()]
        # Assign the combined var DataFrame back to aggregated_adata
        aggregated_adata.var = all_var.loc[aggregated_adata.var_names]
        
        print(f"Aggregated AnnData object created with {aggregated_adata.n_obs} observations and donor IDs preserved.")
        return aggregated_adata
    else:
        print("No AnnData files were found to aggregate.")
        return None

def aggregate_across_datasets(dataset_to_donor_map, organism="Mus musculus", main_directory="datasets"):
    """
    Aggregates the AnnData files for multiple dataset_ids and their respective donor_ids.

    Parameters:
    - dataset_to_donor_map (dict): A dictionary where keys are dataset IDs and values are lists of donor IDs.
    - organism (str): Organism name to use in the census query (default is 'Mus musculus').
    - main_directory (str): Directory where AnnData slices are stored.

    Returns:
    - AnnData: A single AnnData object containing aggregated data across all dataset_ids.
    """
    # Initialize a list to collect aggregated AnnData objects for each dataset
    dataset_adata_list = []

    for dataset_id, donor_ids in dataset_to_donor_map.items():
        print(f"Aggregating donors for dataset_id '{dataset_id}' with donors {donor_ids}...")
        
        # Aggregate all donor IDs for the current dataset
        aggregated_adata = aggregate_donor_id(dataset_id, donor_ids, organism=organism, main_directory=main_directory)
        
        if aggregated_adata is not None:
            # Add dataset_id information to the .obs for the concatenation
            aggregated_adata.obs["dataset_id"] = dataset_id
            dataset_adata_list.append(aggregated_adata)

    # Concatenate all dataset-level aggregated AnnData objects
    if dataset_adata_list:
        # Use `keys` to keep track of both donor_id and dataset_id
        aggregated_datasets = anndata.concat(dataset_adata_list, join="outer", index_unique="-")
        # Grab all var DataFrames from our list of AnnData objects
        all_var = [x.var for x in dataset_adata_list]
        # Concatenate them
        all_var = pd.concat(all_var, join="outer")
        # Remove duplicates
        all_var = all_var[~all_var.index.duplicated()]
        # Assign the combined var DataFrame back to aggregated_adata
        aggregated_datasets.var = all_var.loc[aggregated_adata.var_names]
        # Adding a new 'dataset_id' column to the combined dataset
        # This ensures that we preserve both donor_id and dataset_id
        print(f"Final aggregated AnnData object created across {len(dataset_to_donor_map)} datasets with {aggregated_datasets.n_obs} observations.")
        return aggregated_datasets
    else:
        print("No data found to aggregate across datasets.")
        return None


from anndata import AnnData

def save_anndata_sliced_h5ad(adata: AnnData, output_dir: str, organism = 'Mus musculus'):
    """
    Save an AnnData object in slices organized by dataset_id and donor_id as .h5ad files.
    
    Args:
        adata (AnnData): The input AnnData object with `dataset_id` and `donor_id` in `adata.obs`.
        output_dir (str): The directory to save the sliced data.
        
    Returns:
        None
    """
    # Check if required columns are in the AnnData object
    if 'dataset_id' not in adata.obs or 'donor_id' not in adata.obs:
        raise ValueError("The AnnData object must have 'dataset_id' and 'donor_id' in `adata.obs`.")
    
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Group the data by dataset_id
    for dataset_id in adata.obs['dataset_id'].unique():
        # Create a subfolder for the dataset_id
        dataset_dir = os.path.join(output_dir, dataset_id)
        os.makedirs(dataset_dir, exist_ok=True)
        
        # Filter the data for the current dataset_id
        dataset_data = adata[adata.obs['dataset_id'] == dataset_id]
        
        # Group by donor_id within the dataset_id
        for donor_id in dataset_data.obs['donor_id'].unique():
            # Filter the data for the current donor_id
            donor_data = dataset_data[dataset_data.obs['donor_id'] == donor_id].copy()
            
            # File name: dataset_id_donor_id.h5ad
            file_name = f"{dataset_id}_{donor_id}.h5ad"
            file_path = os.path.join(dataset_dir, file_name)
            update_metadata_file(dataset_id, donor_id, organism, output_dir)
            print('metadata_file')
            # Save the data as .h5ad
            donor_data.write_h5ad(file_path)
            
            print(f"Saved: {file_path}")
