# Protos Documentation

This document serves as the central reference point for all Protos documentation.

## Core Concepts

Protos is organized around several key concepts that work together to provide a comprehensive framework for protein structure analysis:

1. **[Processor Classes](PROCESSOR_CLASSES.md)** - The backbone of Protos functionality
   - GRNProcessor - Generic Residue Numbering
   - EMBProcessor - Embeddings
   - CifProcessor - Structure processing
   - PropertyProcessor - Property management

2. **[Path Management](PATH_CONFIG.md)** - How Protos handles file and directory paths
   - Environment variables
   - Programmatic configuration
   - Path resolution
   - Cross-platform compatibility

3. **[Directory Structure](PATH_STRUCTURE.md)** - The standardized directory hierarchy
   - Reference vs. user data separation
   - Registry system
   - Processor-specific organization

4. **[File Formats](FILE_FORMATS.md)** - Standard formats for Protos data
   - Sequence files (FASTA)
   - Structure files (PDB, mmCIF)
   - Table formats (GRN tables, property tables)
   - Embeddings, graphs, and temporary files

5. **[Dataset System](DATASET_SYSTEM.md)** - Working with collections of related data
   - Dataset class
   - DatasetManager
   - Creating, loading, and merging datasets

6. **[Generic Residue Numbering (GRN) System](GRN_SYSTEM.md)** - Standardized positional referencing
   - Notation formats
   - Assignment workflows
   - Utility functions
   - Implementation details

## Common Workflows

These core concepts integrate to support common workflows in Protos:

1. **Structure Analysis Pipeline**
   - Loading and processing structure files
   - Applying GRN assignments
   - Structure alignment and comparison

2. **Embedding Generation and Analysis**
   - Generating embeddings from sequences or structures
   - Analyzing and visualizing embedding spaces

3. **Property Management**
   - Assigning and tracking properties
   - Property-based filtering and analysis

4. **Dataset Manipulation**
   - Creating standardized datasets
   - Merging and transforming datasets

## Getting Started

For first-time users, we recommend reading the documentation in this order:
1. PATH_CONFIG.md - To understand how to set up your environment
2. PATH_STRUCTURE.md - To understand how data is organized
3. PROCESSOR_CLASSES.md - To understand the core functionality
4. FILE_FORMATS.md - To understand the data formats
5. DATASET_SYSTEM.md - To understand how to work with datasets
6. GRN_SYSTEM.md - To understand the specialized GRN system

For detailed API documentation, refer to the inline code documentation.