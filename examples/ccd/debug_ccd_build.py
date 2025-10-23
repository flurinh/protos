#!/usr/bin/env python3
"""
Debug CCD index building.
"""

import os
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from protos.io.paths.path_config import ProtosPaths
from protos.io.ingest import ccd_loader_unified as ccd_loader

# Set up paths
data_dir = Path("/mnt/c/Users/hidbe/PycharmProjects/protos/data")
os.environ["PROTOS_DATA_ROOT"] = str(data_dir.absolute())

# Initialize ProtosPaths
paths = ProtosPaths(data_root=str(data_dir.absolute()))

# Get CCD directory
ccd_dir = Path(paths.get_subdir_path("molecule", "cache_dir")) / "databases" / "ccd"

print(f"CCD directory: {ccd_dir}")

# Remove old index files
print("\nRemoving old index files...")
for f in ['ccd_index.json', 'ccd_ligands.pkl']:
    fp = ccd_dir / f
    if fp.exists():
        fp.unlink()
        print(f"  Removed {f}")

# Force rebuild with debug
print("\nRebuilding index...")
import logging
logging.basicConfig(level=logging.DEBUG)

success = ccd_loader.build_ccd_index(ccd_dir, force=True)

if success:
    print("\n✅ Build successful!")
    stats = ccd_loader.get_ccd_statistics(ccd_dir)
    print(f"  Total: {stats['total_components']}")
    print(f"  With SMILES: {stats['with_smiles']}")
else:
    print("\n✗ Build failed!")