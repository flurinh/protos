import os
import requests
from pathlib import Path

# Import ProtosPaths for centralized path management
try:
    from protos.io.paths.path_config import ProtosPaths
    _HAS_PROTOS_PATHS = True
except ImportError:
    _HAS_PROTOS_PATHS = False


def download_alphafold_structures(uid, max_models=5, output_dir=None, processor=None):
    # Determine output path using ProtosPaths (preferred) or legacy fallback
    if processor is not None:
        # Use processor's structure directory (follows core philosophy)
        output_path = processor.path_structure_dir
    elif _HAS_PROTOS_PATHS and output_dir is None:
        # Use ProtosPaths directly if no processor provided
        paths = ProtosPaths()
        output_path = Path(paths.get_structure_subdir_path('structure_dir'))
    else:
        # Legacy fallback for backward compatibility
        output_path = Path(output_dir or 'data/mmcif/alphafold_structures')
        output_path.mkdir(parents=True, exist_ok=True)

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
            filepath = output_path / filename

            with open(filepath, "wb") as f:
                f.write(response.content)

            print(f"Saved AlphaFold structure to {filepath}")
        else:
            print(
                f"Failed to download AlphaFold structure for Uniprot ID {uid}, model v{i}. Status code: {response.status_code}")
        print("\n")


def download_alphafold_with_processor(uid, processor, max_models=5):
    """
    Download AlphaFold structures using a processor for ProtosPaths integration.
    
    This is the preferred method that fully follows the core philosophy.
    
    Args:
        uid: UniProt ID for AlphaFold structure
        processor: StructureProcessor instance (handles all path management)
        max_models: Maximum number of model versions to try
    """
    return download_alphafold_structures(
        uid=uid,
        processor=processor,
        max_models=max_models
    )

