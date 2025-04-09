#!/bin/bash

# Fix import statements for reorganized project structure
# This script systematically updates import paths according to the new protos package structure

set -e

PROJECT_DIR="/mnt/c/Users/hidbe/PycharmProjects/GPCR"
echo "Starting import fix script in $PROJECT_DIR"

# Find all Python files in the project
echo "Finding Python files..."
PYTHON_FILES=$(find "$PROJECT_DIR" -name "*.py" -type f -not -path "*/\.*" -not -path "*/venv/*" -not -path "*/env/*")

# 1. Fix general processing imports
echo "Fixing general processing imports..."
for file in $PYTHON_FILES; do
  echo "Processing file: $file"
  
  # from processing.X import Y -> from protos.processing.X import Y
  sed -i 's/from processing\./from protos.processing./g' "$file"
  
  # Fix module-specific imports based on new organization
  
  # Structure-related 
  sed -i 's/from protos.processing.struct_\([a-z_]*\) import/from protos.processing.structure.struct_\1 import/g' "$file"
  sed -i 's/from protos.processing.binding_domain_alignment import/from protos.processing.structure.binding_domain_alignment import/g' "$file"
  sed -i 's/from protos.processing.struct_heterogeneity import/from protos.processing.structure.struct_heterogeneity import/g' "$file"
  sed -i 's/from protos.processing.struct_alignment import/from protos.processing.structure.struct_alignment import/g' "$file"
  sed -i 's/from protos.processing.structure_bond_utils import/from protos.processing.structure.structure_bond_utils import/g' "$file"
  
  # Sequence-related
  sed -i 's/from protos.processing.seq_\([a-z_]*\) import/from protos.processing.sequence.seq_\1 import/g' "$file"
  
  # GRN-related
  sed -i 's/from protos.processing.grn_\([a-z_]*\) import/from protos.processing.grn.grn_\1 import/g' "$file"
  sed -i 's/from protos.processing.init_default_grn_table import/from protos.processing.grn.init_default_grn_table import/g' "$file"
  sed -i 's/from protos.processing.grn_dict import/from protos.processing.grn.grn_dict import/g' "$file"
  
  # Graph-related
  sed -i 's/from protos.processing.graph_\([a-z_]*\) import/from protos.processing.graph.graph_\1 import/g' "$file"
  
  # Embedding-related
  sed -i 's/from protos.processing.emb_\([a-z_]*\) import/from protos.embedding.emb_\1 import/g' "$file"
  sed -i 's/from protos.processing.embedding import/from protos.embedding.embedding import/g' "$file"
  
  # Table-related
  sed -i 's/from protos.processing.gene_table_utils import/from protos.processing.tables.gene_table_utils import/g' "$file"
  sed -i 's/from protos.processing.grn_table_\([a-z_]*\) import/from protos.processing.tables.grn_table_\1 import/g' "$file"
  sed -i 's/from protos.processing.grn_table_utils import/from protos.processing.tables.grn_table_utils import/g' "$file"
  
  # Visualization-related
  sed -i 's/from protos.processing.struct_visualization import/from protos.visualization.structure_vis import/g' "$file"
  sed -i 's/from protos.processing.graph_visualization import/from protos.visualization.graph_vis import/g' "$file"
  sed -i 's/from protos.processing.ligand_visualization import/from protos.visualization.ligand_vis import/g' "$file"
  
  # Property and IO related
  sed -i 's/from protos.processing.property_processor import/from protos.processing.property.property_processor import/g' "$file"
  sed -i 's/from protos.processing.mapping_processor import/from protos.processing.property.mapping_processor import/g' "$file"
  sed -i 's/from protos.processing.function_processor import/from protos.processing.property.function_processor import/g' "$file"
  sed -i 's/from protos.processing.fasta_utils import/from protos.io.fasta_utils import/g' "$file"
  
  # Ligand-related
  sed -i 's/from protos.processing.ligand_\([a-z_]*\) import/from protos.processing.ligand.ligand_\1 import/g' "$file"
  sed -i 's/from protos.processing.retinal import/from protos.processing.ligand.retinal import/g' "$file"
  
  # Opsin-specific
  sed -i 's/from protos.processing.opsin_\([a-z_]*\) import/from protos.processing.opsin.opsin_\1 import/g' "$file"
  sed -i 's/from protos.processing.retinal_mo_utils import/from protos.processing.opsin.retinal_mo_utils import/g' "$file"
done

# 2. Fix data_utils imports
echo "Fixing data_utils imports..."
for file in $PYTHON_FILES; do
  # from data_utils.uniprot.X import Y -> from protos.loaders.X import Y
  sed -i 's/from data_utils\.uniprot\./from protos.loaders./g' "$file"
  
  # from data_utils.structures.X import Y -> from protos.loaders.X import Y
  sed -i 's/from data_utils\.structures\./from protos.loaders./g' "$file"
  
  # from data_utils.alphafold.X import Y -> from protos.loaders.X import Y
  sed -i 's/from data_utils\.alphafold\./from protos.loaders./g' "$file"
  
  # from data_utils.gpcrdb.X import Y -> from protos.loaders.X import Y
  sed -i 's/from data_utils\.gpcrdb\./from protos.loaders./g' "$file"
  
  # from data_utils.X import Y -> from protos.loaders.X import Y (general case)
  sed -i 's/from data_utils\./from protos.loaders./g' "$file"
  
  # Fix 'import data_utils.X' cases
  sed -i 's/import data_utils\./import protos.loaders./g' "$file"
done

# 3. Fix wildcard imports (* imports)
echo "Fixing wildcard imports..."
for file in $PYTHON_FILES; do
  # Fix wildcard imports for processing module
  sed -i 's/from processing\.\([a-z_]*\) import \*/from protos.processing.\1 import */g' "$file"
  
  # Fix structure-related wildcard imports
  sed -i 's/from protos.processing.struct_\([a-z_]*\) import \*/from protos.processing.structure.struct_\1 import */g' "$file"
  
  # Fix sequence-related wildcard imports
  sed -i 's/from protos.processing.seq_\([a-z_]*\) import \*/from protos.processing.sequence.seq_\1 import */g' "$file"
  
  # Fix GRN-related wildcard imports
  sed -i 's/from protos.processing.grn_\([a-z_]*\) import \*/from protos.processing.grn.grn_\1 import */g' "$file"
done

# 4. Special handling for projects directory
echo "Special handling for projects directory..."
# Handle projects/opsin_analysis imports
OPSIN_FILES=$(find "$PROJECT_DIR/projects/opsin_analysis" -name "*.py" -type f)
for file in $OPSIN_FILES; do
  # Make sure opsin_analysis imports reference the protos package properly
  sed -i 's/from processing\.opsin_\([a-z_]*\) import/from protos.processing.opsin.opsin_\1 import/g' "$file"
  sed -i 's/from processing\.retinal\(_[a-z_]*\)\? import/from protos.processing.opsin.retinal\1 import/g' "$file"
  
  # Ensure other protos imports are correct
  sed -i 's/from processing\./from protos.processing./g' "$file"
done

# 5. Fix remaining special cases
echo "Fixing remaining special cases..."
for file in $PYTHON_FILES; do
  # Make sure visualization imports are consistently named
  sed -i 's/visualize_structure/structure_vis/g' "$file"
  
  # Make sure any remaining structure visualization imports are correct
  sed -i 's/from protos.processing.structure.struct_visualization import/from protos.visualization.structure_vis import/g' "$file"
  
  # Clean up any double protos references
  sed -i 's/from protos.protos\./from protos./g' "$file"
  
  # Fix any remaining direct imports
  sed -i 's/import processing\./import protos.processing./g' "$file"
  sed -i 's/import data_utils\./import protos.loaders./g' "$file"
done

echo "Import fixing complete! All Python files in $PROJECT_DIR have been processed."