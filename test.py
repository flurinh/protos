
import protos

# Set custom data path (optional - defaults to ./data)
protos.set_data_path("/data/fast/projects/protos/data")

# Get current path
print(protos.get_data_path())

from protos.io.ingest.structure_loader import StructureLoader

loader = StructureLoader()

# Download from RCSB PDB
name = loader.download_and_register("1ubq")

# Download from AlphaFold
name = loader.download_and_register("P00533", source="alphafold")

# With custom name
name = loader.download_and_register("P00533", name="EGFR_HUMAN", source="alphafold")

# Import local file
name = loader.download_and_register("/path/to/structure.cif", source="local")