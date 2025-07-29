#!/usr/bin/env python
"""
Patch to fix ChEMBL import issue when the service is down.
This moves the import inside the class methods where it's actually used.
"""

import os
import shutil
from pathlib import Path

# Find the chembl_loader.py file
chembl_loader_path = Path(__file__).parent / "src/protos/loaders/chembl_loader.py"

# Create backup
backup_path = chembl_loader_path.with_suffix('.py.bak')
shutil.copy2(chembl_loader_path, backup_path)
print(f"Created backup at {backup_path}")

# Read the file
with open(chembl_loader_path, 'r') as f:
    content = f.read()

# Replace the problematic import section
old_import = """# Try to import ChEMBL client
try:
    from chembl_webresource_client.new_client import new_client as chembl_client
    HAS_CHEMBL = True
except ImportError:
    logger.warning("chembl_webresource_client not available. ChEMBL functionality will be limited.")
    HAS_CHEMBL = False
    chembl_client = None"""

new_import = """# ChEMBL client will be imported on demand
HAS_CHEMBL = None
chembl_client = None

def _get_chembl_client():
    \"\"\"Lazy import of ChEMBL client.\"\"\"
    global HAS_CHEMBL, chembl_client
    
    if HAS_CHEMBL is None:
        try:
            from chembl_webresource_client.new_client import new_client as _chembl_client
            chembl_client = _chembl_client
            HAS_CHEMBL = True
        except (ImportError, Exception) as e:
            logger.warning(f"chembl_webresource_client not available: {e}")
            HAS_CHEMBL = False
            chembl_client = None
    
    return chembl_client"""

# Replace the import
content = content.replace(old_import, new_import)

# Also need to update references to chembl_client to use the getter
content = content.replace("if not HAS_CHEMBL:", "if not _get_chembl_client():")
content = content.replace("chembl_client.target", "_get_chembl_client().target")
content = content.replace("chembl_client.activity", "_get_chembl_client().activity")
content = content.replace("chembl_client.molecule", "_get_chembl_client().molecule")

# Write the patched file
with open(chembl_loader_path, 'w') as f:
    f.write(content)

print(f"Patched {chembl_loader_path}")
print("\nThe ChEMBL client will now be imported lazily only when needed.")
print("This prevents import errors when the ChEMBL service is down.")