import os
import requests


def download_alphafold_structures(uid, max_models=5, output_dir='data/mmcif/alphafold_structures'):
    # Ensure the output directory exists
    if not os.path.isdir(output_dir):
        print(f"Creating directory {output_dir}...")
        os.makedirs(output_dir)

    # Try to download each model
    for i in range(1, max_models + 1):
        # Define the URL for this model
        url = f"https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v{i}.cif"

        # Send a GET request to the URL
        print(f"Downloading AlphaFold structure from {url}...")
        response = requests.get(url)

        # If the request was successful, save the file
        if response.status_code == 200:
            filename = f"AF-{uid}-F1-model_v{i}.cif"
            filepath = os.path.join(output_dir, filename)

            with open(filepath, "wb") as f:
                f.write(response.content)

            print(f"Saved AlphaFold structure to {filepath}")
        else:
            print(
                f"Failed to download AlphaFold structure for Uniprot ID {uid}, model v{i}. Status code: {response.status_code}")
        print("\n")

