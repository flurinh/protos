from protos.processing import *
from protos.io import *
from protos.visualization import *

# Public API for path management
from protos.io.paths import get_protos_paths


def set_data_path(path: str):
    """
    Set the global data path for all Protos operations.
    
    IMPORTANT: This must be called BEFORE creating any processors!
    Once a processor is initialized, the path is locked and cannot be changed.
    
    Args:
        path: The directory path where all Protos data will be stored.
        
    Raises:
        RuntimeError: If attempting to change path after initialization.
    """
    get_protos_paths(path)


def get_data_path() -> str:
    """
    Get the current global data path.
    
    Returns:
        The absolute path to the Protos data directory.
    """
    return get_protos_paths().data_root
