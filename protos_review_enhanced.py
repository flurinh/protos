#!/usr/bin/env python3
"""
Protos Comprehensive Review & Documentation Script

This script serves as both a comprehensive review and documentation of all Protos features.
It demonstrates the complete architecture from core components to advanced applications,
with detailed comments explaining how each component works.

ARCHITECTURE OVERVIEW:
┌─────────────────────────────────────────────────────────────────┐
│                           PROTOS ARCHITECTURE                   │
├─────────────────────────────────────────────────────────────────┤
│ 1. CORE INFRASTRUCTURE                                          │
│    ├─ ProtosPaths: Centralized path management                  │
│    ├─ BaseProcessor: Abstract base for all data processors     │
│    └─ Entity Registry: Universal entity tracking system        │
│                                                                 │
│ 2. SPECIALIZED PROCESSORS                                       │
│    ├─ CifBaseProcessor: 3D protein structure management        │
│    ├─ GRNBaseProcessor: Generic Residue Numbering system       │
│    ├─ SeqProcessor: Sequence data management                   │
│    ├─ EmbeddingProcessor: ML embeddings generation             │
│    └─ PropertyProcessor: Metadata integration                  │
│                                                                 │
│ 3. DATA ABSTRACTION LAYERS                                     │
│    ├─ Entities: Individual data items (single sequence/struct) │
│    ├─ Datasets: Collections of related entities                │
│    └─ Cross-format tracking: Same entity across formats        │
│                                                                 │
│ 4. APPLICATIONS                                                 │
│    ├─ CLI tools: Command-line utilities                        │
│    ├─ Analysis workflows: Data processing pipelines            │
│    └─ Integration APIs: External tool interfaces               │
└─────────────────────────────────────────────────────────────────┘

Usage:
    python protos_review_enhanced.py                    # Full review
    python protos_review_enhanced.py --core-only        # Core components only
    python protos_review_enhanced.py --processors-only  # Processors only
    python protos_review_enhanced.py --applications-only # Applications only
"""

import os
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import logging
from typing import List, Dict, Optional, Tuple, Any
import json

# Configure logging for educational purposes
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def print_section(title: str, description: str = ""):
    """Print a formatted section header with description."""
    print(f"\n{'=' * 80}")
    print(f"{title:^80}")
    if description:
        print(f"{description:^80}")
    print('=' * 80)


def print_subsection(title: str, description: str = ""):
    """Print a formatted subsection header with description."""
    print(f"\n{'-' * 60}")
    print(f"{title}")
    if description:
        print(f"   {description}")
    print('-' * 60)


def print_concept(concept: str, explanation: str):
    """Print a formatted concept explanation."""
    print(f"\n💡 CONCEPT: {concept}")
    print(f"   {explanation}")


def print_code_flow(step: str, action: str):
    """Print a formatted code flow step."""
    print(f"\n⚙️  STEP: {step}")
    print(f"   ACTION: {action}")


class ProtosComprehensiveReview:
    """
    Comprehensive review and documentation of the Protos architecture.
    
    This class demonstrates every aspect of Protos with educational comments
    explaining the design patterns, abstractions, and best practices.
    """
    
    def __init__(self):
        """
        Initialize the review with all Protos components.
        
        INITIALIZATION PATTERN:
        - Import all components in try/catch blocks for graceful degradation
        - Store components as instance attributes for easy access
        - Provide clear error messages if components are missing
        """
        print_section("PROTOS INITIALIZATION", "Loading and validating all components")
        
        # Import core infrastructure components
        print_concept("Core Infrastructure", 
                     "The foundation components that all other parts depend on")
        
        try:
            # 1. PATH MANAGEMENT - The cornerstone of Protos architecture
            from protos.io.paths.path_config import ProtosPaths
            self.ProtosPaths = ProtosPaths
            print("✅ ProtosPaths: Centralized path management system")
            print("   │ ▶ Manages ALL file system interactions")
            print("   │ ▶ Eliminates hardcoded paths throughout codebase")
            print("   └ ▶ Provides platform-independent path resolution")
            
            # 2. ENTITY SYSTEM - Universal data tracking
            from protos.io.data_access import EntityRegistry, GlobalRegistry, generate_entity_id
            self.GlobalRegistry = GlobalRegistry
            self.generate_entity_id = generate_entity_id
            print("✅ Entity System: Universal data tracking across formats")
            print("   │ ▶ Each data item gets a unique hash-based ID")
            print("   │ ▶ Tracks same entity across multiple data formats")
            print("   └ ▶ Enables cross-format queries and analysis")
            
            # 3. BASE PROCESSOR - Abstract foundation for all processors
            from protos.core.base_processor import BaseProcessor
            self.BaseProcessor = BaseProcessor
            print("✅ BaseProcessor: Abstract base class for all data processors")
            print("   │ ▶ Provides common interface for data management")
            print("   │ ▶ Handles entity registration automatically")
            print("   └ ▶ Standardizes save/load operations across formats")
            
        except ImportError as e:
            print(f"❌ Core infrastructure import error: {e}")
            sys.exit(1)
        
        print_concept("Specialized Processors", 
                     "Domain-specific processors for different data types")
        
        try:
            # 4. STRUCTURE PROCESSOR - 3D protein structure management
            from protos.processing.structure.struct_base_processor import CifBaseProcessor
            self.CifBaseProcessor = CifBaseProcessor
            print("✅ CifBaseProcessor: 3D protein structure management")
            print("   │ ▶ Parses PDB/mmCIF files using BioPython")
            print("   │ ▶ Provides structure analysis capabilities")
            print("   │ ▶ Handles coordinate transformations and filtering")
            print("   └ ▶ Integrates with structure alignment tools")
            
            # 5. GRN PROCESSOR - Generic Residue Numbering system
            from protos.processing.grn.grn_base_processor import GRNBaseProcessor
            self.GRNBaseProcessor = GRNBaseProcessor
            print("✅ GRNBaseProcessor: Generic Residue Numbering system")
            print("   │ ▶ Standardizes residue numbering across protein families")
            print("   │ ▶ Enables structural comparison of homologous proteins")
            print("   │ ▶ Supports GPCRs, opsins, and other membrane proteins")
            print("   └ ▶ Provides sequence-to-structure mapping")
            
            # 6. SEQUENCE PROCESSOR - Sequence data management
            from protos.processing.sequence.seq_processor import SeqProcessor
            self.SeqProcessor = SeqProcessor
            print("✅ SeqProcessor: Sequence data management")
            print("   │ ▶ Handles FASTA files and sequence databases")
            print("   │ ▶ Provides sequence analysis and alignment tools")
            print("   │ ▶ Integrates with external sequence databases")
            print("   └ ▶ Supports batch sequence operations")
            
            # 7. PROPERTY PROCESSOR - Metadata and property management
            from protos.processing.property.property_processor_enhanced import PropertyProcessor
            self.PropertyProcessor = PropertyProcessor
            print("✅ PropertyProcessor: Metadata and property management")
            print("   │ ▶ Associates properties with entities across all formats")
            print("   │ ▶ User-friendly identifiers (PDB IDs, FASTA names)")
            print("   │ ▶ CSV import workflows for experimental data")
            print("   └ ▶ Cross-format property queries and filtering")
            
            # 8. EMBEDDING PROCESSOR - ML embeddings generation
            from protos.processing.embedding.embedding_processor import EmbeddingProcessor
            self.EmbeddingProcessor = EmbeddingProcessor
            print("✅ EmbeddingProcessor: ML embeddings generation")
            print("   │ ▶ Generates protein embeddings using transformer models")
            print("   │ ▶ Supports ESM-2, Ankh, and other protein language models")
            print("   │ ▶ Provides mean, CLS, and per-residue embeddings")
            print("   └ ▶ Caches embeddings for reuse")
            
        except ImportError as e:
            print(f"⚠️  Some processors not available: {e}")
            # Set missing processors to None for graceful degradation
            if 'CifBaseProcessor' not in locals(): self.CifBaseProcessor = None
            if 'GRNBaseProcessor' not in locals(): self.GRNBaseProcessor = None
            if 'SeqProcessor' not in locals(): self.SeqProcessor = None
            if 'PropertyProcessor' not in locals(): self.PropertyProcessor = None
            if 'EmbeddingProcessor' not in locals(): self.EmbeddingProcessor = None
        
        print_concept("Data Loaders", 
                     "Tools for downloading data from external databases")
        
        try:
            # 8. DATA LOADERS - External database integration
            from protos.loaders.download_structures import download_protein_structures
            from protos.loaders.alphafold_utils import download_alphafold_structures
            from protos.loaders.uniprot_utils import download_sequences_from_uniprot
            
            self.download_protein_structures = download_protein_structures
            self.download_alphafold_structures = download_alphafold_structures
            self.download_sequences_from_uniprot = download_sequences_from_uniprot
            
            print("✅ Data Loaders: External database integration")
            print("   │ ▶ PDB: Experimental protein structures")
            print("   │ ▶ AlphaFold: Predicted protein structures")
            print("   │ ▶ UniProt: Protein sequences and annotations")
            print("   └ ▶ Automatic data management via ProtosPaths")
            
        except ImportError as e:
            print(f"⚠️  Some loaders not available: {e}")
            self.download_protein_structures = None
            self.download_alphafold_structures = None
            self.download_sequences_from_uniprot = None
        
        print_concept("CLI Applications", 
                     "Command-line tools for data management and analysis")
        
        try:
            # 9. CLI TOOLS - Command-line applications
            from protos.cli.init_data import init_data_directory
            self.init_data_directory = init_data_directory
            print("✅ CLI Tools: Command-line applications")
            print("   │ ▶ Data directory initialization")
            print("   │ ▶ GRN assignment and analysis")
            print("   │ ▶ Structure processing workflows")
            print("   └ ▶ Batch processing utilities")
            
            # Try to import GRN-specific CLI tools
            try:
                from protos.cli.grn.assign_grns import assign_grns_to_entities
                from protos.processing.grn.grn_table_utils import GRNConfigManager
                self.assign_grns_to_entities = assign_grns_to_entities
                self.GRNConfigManager = GRNConfigManager
                print("✅ GRN CLI Tools: Specialized GRN analysis utilities")
            except ImportError:
                self.assign_grns_to_entities = None
                self.GRNConfigManager = None
                print("⚠️  GRN CLI tools not fully available")
                
        except ImportError as e:
            print(f"⚠️  CLI tools not available: {e}")
            self.init_data_directory = None
        
        print("\n🎉 Initialization complete! Ready to demonstrate Protos capabilities.")
    
    def demonstrate_core_infrastructure(self) -> Path:
        """
        Demonstrate the core infrastructure components that form Protos' foundation.
        
        CORE DESIGN PRINCIPLES:
        1. Centralized path management eliminates hardcoded paths
        2. Entity system provides universal data tracking
        3. Registry system enables cross-format queries
        
        Returns:
            Path: The configured data directory
        """
        print_section("PART 1: CORE INFRASTRUCTURE", 
                     "The foundation components that make everything else possible")
        
        print_subsection("1.1 ProtosPaths: Centralized Path Management", 
                        "How Protos eliminates hardcoded paths throughout the system")
        
        print_concept("Path Management Philosophy", 
                     "Users work with names, Protos handles ALL file system complexity")
        
        # STEP 1: Determine data directory location
        print_code_flow("Data Directory Resolution", 
                       "Intelligently locate project data directory")
        
        cwd = Path.cwd()
        if 'protos' in cwd.parts:
            idx = cwd.parts.index('protos')
            protos_path = Path(*cwd.parts[:idx+1])
        else:
            protos_path = cwd
        
        data_dir = protos_path / "data"
        print(f"   📂 Project path: {protos_path}")
        print(f"   📂 Data directory: {data_dir}")
        print(f"   💡 This becomes the ROOT of all data operations")
        
        # STEP 2: Initialize data directory if needed
        print_code_flow("Data Directory Initialization", 
                       "Create directory structure and registries")
        
        entity_registry = data_dir / "entity_registry.json"
        
        if not entity_registry.exists():
            print("   🔧 Initializing data directory structure...")
            stats = self.init_data_directory(data_root=str(data_dir), force=True)
            print("   ✅ Directory structure created:")
            print(f"      ├─ Directories: {stats.get('directories_created', 0)}")
            print(f"      ├─ Registries: {stats.get('registries_created', 0)}")
            print(f"      └─ Configuration files initialized")
        else:
            print("   ✅ Data directory already initialized")
        
        # STEP 3: Set global data root - CRITICAL OPERATION
        print_code_flow("Global Path Configuration", 
                       "Set the single source of truth for all paths")
        
        self.ProtosPaths.set_data_root(str(data_dir.absolute()))
        paths = self.ProtosPaths()
        
        print(f"   🌍 Global data root: {data_dir.absolute()}")
        print("   📋 Key directories managed by ProtosPaths:")
        print(f"      ├─ Structure data: {data_dir}/{paths.processor_dirs.get('structure', 'structure')}")
        print(f"      ├─ Sequence data: {data_dir}/{paths.processor_dirs.get('sequence', 'sequence')}")
        print(f"      ├─ GRN data: {data_dir}/{paths.processor_dirs.get('grn', 'grn')}")
        print(f"      ├─ Embedding data: {data_dir}/{paths.processor_dirs.get('embedding', 'embedding')}")
        print(f"      └─ Property data: {data_dir}/{paths.processor_dirs.get('property', 'property')}")
        
        print_concept("Path Abstraction Benefits", 
                     "No more hardcoded paths, platform independence, centralized configuration")
        
        print_subsection("1.2 Entity System: Universal Data Tracking", 
                        "How Protos tracks individual data items across formats")
        
        print_concept("Entity Philosophy", 
                     "Every piece of data is an 'entity' with a unique ID and metadata")
        
        # Demonstrate entity ID generation
        test_names = ["1ubq", "TEST_PROTEIN_A", "BACR_HALSA"]
        print("   🔍 Entity ID generation examples:")
        
        for name in test_names:
            entity_id = self.generate_entity_id(name)
            print(f"      '{name}' → {entity_id}")
        
        print("   💡 Same name always generates same ID (deterministic hashing)")
        print("   💡 IDs are format-agnostic (same entity across file types)")
        
        # Demonstrate global registry
        print_code_flow("Global Registry Access", 
                       "Central registry tracks all entities and their formats")
        
        global_registry = self.GlobalRegistry()
        print("   ✅ Global registry initialized")
        print("   📊 Registry functions:")
        print("      ├─ Track entities across multiple data formats")
        print("      ├─ Resolve names to canonical entity IDs")
        print("      ├─ Store metadata for each entity")
        print("      └─ Enable cross-format queries")
        
        print_subsection("1.3 BaseProcessor: Common Interface", 
                        "Abstract base class that standardizes all processors")
        
        print_concept("Processor Pattern", 
                     "All processors inherit common functionality while specializing for data types")
        
        # Demonstrate BaseProcessor capabilities
        print("   🏗️  BaseProcessor provides:")
        print("      ├─ Automatic path management via ProtosPaths")
        print("      ├─ Entity registration and tracking")
        print("      ├─ Dataset creation and management")
        print("      ├─ Save/load operations with format detection")
        print("      ├─ Metadata handling and versioning")
        print("      └─ Logging and error handling")
        
        print("   ⚙️  All processors inherit these capabilities:")
        print("      ├─ CifBaseProcessor (structures)")
        print("      ├─ GRNBaseProcessor (residue numbering)")
        print("      ├─ SeqProcessor (sequences)")
        print("      ├─ EmbeddingProcessor (ML embeddings)")
        print("      └─ PropertyProcessor (metadata)")
        
        print("\n✨ Core infrastructure established! This foundation enables:")
        print("   • No hardcoded paths anywhere in the system")
        print("   • Universal entity tracking across all data types")
        print("   • Consistent interfaces for all data operations")
        print("   • Automatic metadata management and versioning")
        
        return data_dir
    
    def demonstrate_processor_architecture(self):
        """
        Demonstrate how processors abstract data management with practical examples.
        
        PROCESSOR ARCHITECTURE:
        - Each processor specializes in one data type
        - Common interface via BaseProcessor inheritance
        - Automatic entity registration and cross-format tracking
        - Built-in data download, save, and load capabilities
        """
        print_section("PART 2: PROCESSOR ARCHITECTURE", 
                     "Specialized data management with automatic abstraction")
        
        # Return dictionary to store processor instances
        processors = {}
        
        # ====================================================================
        # 2.1 STRUCTURE PROCESSOR - 3D protein structure management
        # ====================================================================
        
        if self.CifBaseProcessor:
            print_subsection("2.1 CifBaseProcessor: 3D Structure Management", 
                            "Download, parse, and analyze protein structures")
            
            print_concept("Structure Data Abstraction", 
                         "Handles PDB/mmCIF files, coordinates, chains, and structural analysis")
            
            # Initialize processor - NOTE: NO PATH PARAMETERS!
            print_code_flow("Processor Initialization", 
                           "Create processor instance - paths handled automatically")
            
            struct_processor = self.CifBaseProcessor(name="review_structures")
            processors['structure'] = struct_processor
            
            print(f"   ✅ CifBaseProcessor initialized")
            print(f"   📂 Data path: {struct_processor.data_path}")
            print(f"   📂 Structure directory: {struct_processor.path_structure_dir}")
            print(f"   💡 All paths managed automatically by ProtosPaths")
            
            # Demonstrate data download with automatic path management
            print_code_flow("Automatic Data Download", 
                           "Download structures with processor managing all paths")
            
            if self.download_protein_structures:
                pdb_ids = ["1ubq", "2gb1", "1crn"]  # Small, well-studied proteins
                print(f"   📥 Downloading PDB structures: {pdb_ids}")
                
                # NOTE: Processor provides the target directory automatically
                successful, failed = self.download_protein_structures(
                    pdb_ids,
                    target_folder=struct_processor.path_structure_dir
                )
                
                print(f"   ✅ Downloaded: {successful}")
                if failed:
                    print(f"   ❌ Failed: {failed}")
            
            # Demonstrate entity vs dataset distinction
            print_code_flow("Entity Management", 
                           "List individual structure entities")
            
            entities = struct_processor.list_entities()
            print(f"   📋 Structure entities available: {len(entities)}")
            for entity in entities[:3]:
                print(f"      - {entity}")
            
            # Load individual entity
            if entities:
                test_entity = entities[0]
                print(f"\n   🔍 Loading individual entity: {test_entity}")
                
                data = struct_processor.load_structure(test_entity, apply_dtypes=True)
                if data is not None:
                    print(f"   ✅ Loaded {len(data)} atoms")
                    print(f"      ├─ Chains: {sorted(data['auth_chain_id'].unique())}")
                    print(f"      ├─ Residues: {data['auth_seq_id'].nunique()}")
                    print(f"      ├─ Data types enforced: {data['x'].dtype}, {data['y'].dtype}, {data['z'].dtype}")
                    print(f"      └─ Coordinate range: X({data['x'].min():.1f} to {data['x'].max():.1f})")
                    
                    # Show entity tracking
                    entity_id = self.generate_entity_id(test_entity)
                    print(f"   🔗 Entity ID: {entity_id}")
            
            # Demonstrate dataset creation and management
            print_code_flow("Dataset Management", 
                           "Create collections of related structures")
            
            dataset_id = "review_small_proteins"
            struct_processor.create_dataset(
                dataset_id=dataset_id,
                name="Small Protein Structures",
                description="Well-studied small proteins for demonstration",
                content=successful[:3] if len(successful) >= 3 else successful
            )
            
            print(f"   📦 Created dataset: {dataset_id}")
            
            # List datasets
            datasets = struct_processor.list_datasets()
            print(f"   📚 Available datasets: {len(datasets)}")
            for ds in datasets[:2]:
                if isinstance(ds, dict):
                    print(f"      - {ds.get('id', 'unknown')}: {ds.get('description', '')}")
            
            # Load dataset
            struct_processor.load_dataset(dataset_id)
            print(f"   ✅ Loaded dataset with {len(struct_processor.pdb_ids)} structures")
            
            # Demonstrate cross-format capability
            print_code_flow("Cross-Format Operations", 
                           "Extract sequences from structure data")
            
            sequences = struct_processor.get_seq_dict()
            print(f"   🧬 Extracted {len(sequences)} sequences from structures:")
            for pdb_chain, seq in list(sequences.items())[:2]:
                print(f"      - {pdb_chain}: {seq[:30]}... ({len(seq)} residues)")
        
        # ====================================================================
        # 2.2 SEQUENCE PROCESSOR - Sequence data management
        # ====================================================================
        
        if self.SeqProcessor:
            print_subsection("2.2 SeqProcessor: Sequence Data Management", 
                            "Handle FASTA files, UniProt data, and sequence analysis")
            
            print_concept("Sequence Data Abstraction", 
                         "Manages FASTA files, database downloads, and sequence operations")
            
            # Initialize processor
            print_code_flow("Processor Initialization", 
                           "Sequence processor with automatic path management")
            
            seq_processor = self.SeqProcessor(name="review_sequences")
            processors['sequence'] = seq_processor
            
            print(f"   ✅ SeqProcessor initialized")
            print(f"   📂 Data path: {seq_processor.data_path}")
            print(f"   📂 FASTA directory: {seq_processor.data_path}/fasta")
            
            # Demonstrate external data download
            print_code_flow("External Database Integration", 
                           "Download sequences from UniProt database")
            
            if self.download_sequences_from_uniprot:
                uniprot_ids = ["P00533", "P04637"]  # EGFR, p53 - well-known proteins
                print(f"   📥 Downloading from UniProt: {uniprot_ids}")
                
                sequences = {}
                for uniprot_id in uniprot_ids:
                    try:
                        seq_data = self.download_sequences_from_uniprot([uniprot_id])
                        if seq_data:
                            sequences.update(seq_data)
                            print(f"   ✅ Downloaded {uniprot_id}")
                    except Exception as e:
                        print(f"   ⚠️ Could not download {uniprot_id}: {e}")
                
                if sequences:
                    seq_processor.save_sequences(sequences, "uniprot_proteins.fasta")
                    print(f"   💾 Saved {len(sequences)} sequences to uniprot_proteins.fasta")
            
            # Create test sequences to demonstrate entity management
            print_code_flow("Entity Creation", 
                           "Create individual sequence entities")
            
            test_sequences = {
                "TEST_PROTEIN_A": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK",
                "TEST_PROTEIN_B": "MGSSHHHHHHSSGLVPRGSHMASMTGGQQMGRGSMKTAYIAKQRQISFVKSH",
                "TEST_PROTEIN_C": "MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDR"
            }
            
            # Save as individual entities
            for name, seq in test_sequences.items():
                seq_processor.save_sequences({name: seq}, f"{name}.fasta")
                print(f"   ✅ Created entity: {name}")
            
            # Save as dataset (collection)
            seq_processor.save_sequences(test_sequences, "test_protein_set.fasta")
            print(f"   📦 Created dataset: test_protein_set.fasta")
            
            # Demonstrate entity vs dataset listing
            print_code_flow("Data Organization", 
                           "Distinguish between entities and datasets")
            
            entities = seq_processor.list_entities()
            print(f"   📋 Individual sequence entities: {len(entities)}")
            for entity in entities[:3]:
                print(f"      - {entity}")
            
            datasets = seq_processor.list_datasets()
            print(f"   📚 Sequence datasets: {len(datasets)}")
            for ds in datasets[:2]:
                if isinstance(ds, dict):
                    print(f"      - {ds.get('id', '')}: {ds.get('sequence_count', '?')} sequences")
            
            # Demonstrate entity access
            if entities:
                test_entity = entities[0]
                print(f"\n   🔍 Accessing entity: {test_entity}")
                
                sequence = seq_processor.get_sequence(test_entity)
                if sequence:
                    print(f"   ✅ Retrieved sequence: {sequence[:40]}...")
                    print(f"      ├─ Length: {len(sequence)} residues")
                    print(f"      └─ Entity ID: {self.generate_entity_id(test_entity)}")
        
        # ====================================================================
        # 2.3 GRN PROCESSOR - Generic Residue Numbering system
        # ====================================================================
        
        if self.GRNBaseProcessor:
            print_subsection("2.3 GRNBaseProcessor: Generic Residue Numbering", 
                            "Standardize residue numbering across protein families")
            
            print_concept("GRN System Philosophy", 
                         "Map equivalent positions across homologous proteins using standard numbering")
            
            # Initialize processor
            print_code_flow("GRN Processor Initialization", 
                           "Specialized processor for residue numbering tables")
            
            grn_processor = self.GRNBaseProcessor(name="review_grn", preload=False)
            processors['grn'] = grn_processor
            
            print(f"   ✅ GRNBaseProcessor initialized")
            print(f"   📂 Data path: {grn_processor.data_path}")
            print(f"   📊 GRN table directory: {grn_processor.data_path}/tables")
            
            # Demonstrate proper GRN table format
            print_code_flow("GRN Table Creation", 
                           "Create table with proper residue+position format")
            
            print("   💡 GRN Format: Each cell contains 'RESIDUE+POSITION' (e.g., M62, K115)")
            
            # Create a GRN table with PROPER FORMAT
            grn_data = pd.DataFrame({
                '1.50': ['M62', 'M62', 'L21', 'V45', 'M62'],    # Helix 1 key position
                '2.50': ['K115', 'G90', 'A65', 'G67', 'K115'], # Helix 2 key position
                '3.50': ['T179', 'S145', 'E107', 'I128', 'R179'], # Helix 3 key position
                '4.50': ['A221', 'S190', 'G149', 'V171', 'W221'], # Helix 4 key position
                '5.50': ['Y270', 'H236', 'E195', 'Y214', 'F270'], # Helix 5 key position
                '6.50': ['I312', 'H280', 'I236', 'I256', 'W312'], # Helix 6 key position
                '7.50': ['A356', 'H324', 'T280', 'A300', 'N356']  # Helix 7 key position
            }, index=['PROTEIN_A', 'PROTEIN_B', 'PROTEIN_C', 'PROTEIN_D', 'PROTEIN_E'])
            
            print("   📊 GRN Table structure:")
            print(f"      ├─ Rows (sequences): {list(grn_data.index)}")
            print(f"      ├─ Columns (GRN positions): {list(grn_data.columns)}")
            print(f"      └─ Values: residue+position format")
            
            print("\n   Sample data:")
            print(grn_data.iloc[:3, :3])
            
            # Save the table
            print_code_flow("GRN Table Storage", 
                           "Save table as dataset with entity registration")
            
            grn_processor.data = grn_data
            grn_processor.save_grn_table("review_alignment")
            print("   💾 Saved GRN table: review_alignment")
            print("   🔗 Each row automatically registered as an entity")
            
            # Demonstrate GRN as dataset containing multiple entities
            print_code_flow("GRN Data Organization", 
                           "Tables as datasets, rows as entities")
            
            entities = grn_processor.list_entities()
            print(f"   📋 GRN entities (protein sequences): {len(entities)}")
            for entity in entities[:3]:
                print(f"      - {entity}")
            
            datasets = grn_processor.list_datasets()
            print(f"   📚 GRN datasets (alignment tables): {len(datasets)}")
            for ds in datasets[:2]:
                if isinstance(ds, dict):
                    print(f"      - {ds.get('id', '')}: {ds.get('entity_count', '?')} sequences")
            
            # Demonstrate GRN analysis capabilities
            print_code_flow("GRN Analysis Functions", 
                           "Built-in analysis tools for GRN tables")
            
            # Load the table for analysis
            grn_processor.load_grn_table("review_alignment")
            print(f"   ✅ Loaded table with {len(grn_processor.data)} sequences")
            
            # Extract sequences from GRN (cross-format operation)
            seq_dict = grn_processor.get_seq_dict()
            print(f"   🧬 Extracted sequences from GRN:")
            for seq_id, sequence in list(seq_dict.items())[:2]:
                print(f"      - {seq_id}: {sequence}")
            
            # Analysis functions
            print("   📊 Built-in analysis functions:")
            print("      ├─ filter_data_by_occurances(): Remove sparse positions")
            print("      ├─ sort_columns(): Order GRN positions")
            print("      ├─ get_position_conservation(): Analyze conservation")
            print("      └─ get_seq_dict(): Extract full sequences")
            
            # Demonstrate microbial opsin-specific functionality
            self.demonstrate_microbial_opsin_grn(grn_processor)
        
        # ====================================================================
        # 2.4 PROPERTY PROCESSOR - Metadata and property management
        # ====================================================================
        
        self.demonstrate_property_processor_integration(processors)
        
        # ====================================================================
        # 2.5 EMBEDDING PROCESSOR - ML embeddings generation
        # ====================================================================
        
        if self.EmbeddingProcessor:
            print_subsection("2.5 EmbeddingProcessor: ML Embeddings Generation", 
                            "Generate and manage protein embeddings using transformer models")
            
            print_concept("Embedding Abstraction", 
                         "Transform sequences into numerical vectors for ML analysis")
            
            # Initialize processor
            print_code_flow("Embedding Processor Initialization", 
                           "Check dependencies and initialize model interface")
            
            emb_processor = self.EmbeddingProcessor(name="review_embeddings")
            processors['embedding'] = emb_processor
            
            print(f"   ✅ EmbeddingProcessor initialized")
            print(f"   📂 Data path: {emb_processor.data_path}")
            
            # Check dependencies
            deps = emb_processor.check_dependencies()
            print("   🔍 ML dependency check:")
            for dep, available in deps.items():
                status = "✅" if available else "❌"
                print(f"      {status} {dep}")
            
            # List available models
            models = emb_processor.list_available_models()
            print(f"   🤖 Available models: {len(models)}")
            for name, info in list(models.items())[:3]:
                print(f"      - {name}: {info['embedding_dim']}D, {info['description']}")
            
            if deps['ready']:
                print_code_flow("Embedding Generation", 
                               "Generate embeddings with automatic entity registration")
                
                # Get sequences from sequence processor if available
                if 'sequence' in processors:
                    seq_proc = processors['sequence']
                    entities = seq_proc.list_entities()[:2]  # First 2 entities
                    
                    if entities:
                        print(f"   🧬 Generating embeddings for {len(entities)} sequences...")
                        
                        for entity in entities:
                            sequence = seq_proc.get_sequence(entity)
                            if sequence and len(sequence) < 500:  # Only short sequences
                                print(f"      Processing {entity} ({len(sequence)} residues)...")
                                
                                try:
                                    # Generate embedding - automatically saves and registers
                                    embedding = emb_processor.embed_sequences(
                                        {entity: sequence},
                                        embedding_type="mean",
                                        register_entities=True
                                    )
                                    
                                    if entity in embedding:
                                        emb_shape = embedding[entity].shape
                                        print(f"      ✅ Generated embedding: {emb_shape}")
                                        print(f"         └─ Automatically registered as entity")
                                except Exception as e:
                                    print(f"      ⚠️ Could not generate embedding: {e}")
                
                # List embedding entities
                emb_entities = emb_processor.list_entities()
                print(f"   📋 Embedding entities: {len(emb_entities)}")
                for entity in emb_entities[:3]:
                    print(f"      - {entity}")
            else:
                print("   ⚠️ ML dependencies not available - skipping generation")
                print("      Install with: pip install torch transformers")
        
        print("\n✨ Processor architecture demonstrated! Key takeaways:")
        print("   • Each processor specializes in one data type")
        print("   • Common interface via BaseProcessor inheritance")  
        print("   • Automatic path management - no hardcoded paths")
        print("   • Built-in entity registration and cross-format tracking")
        print("   • Download, save, and load operations handled seamlessly")
        
        return processors
    
    def demonstrate_microbial_opsin_grn(self, grn_processor):
        """Demonstrate GRN annotation specific to microbial opsins."""
        print_subsection("2.3.1 Microbial Opsin GRN Specialization", 
                        "Domain-specific GRN annotation for rhodopsins")
        
        print_concept("Microbial Opsin Biology", 
                     "Light-driven ion pumps with 7-transmembrane helical structure")
        
        print("   🔬 Key functional positions in microbial opsins:")
        print("      ├─ 1.50: Helix 1 start (often Methionine)")
        print("      ├─ 3.50: DRY motif region equivalent") 
        print("      ├─ 6.50: Proline kink in helix 6 (structural)")
        print("      └─ 7.50: Lysine for Schiff base (CRITICAL for function)")
        
        # Create example annotation for bacteriorhodopsin
        print_code_flow("Bacteriorhodopsin Example", 
                       "GRN annotation for the most studied microbial opsin")
        
        bacr_grn = pd.DataFrame({
            '1.50': ['M62'],    # Start of helix 1
            '2.50': ['V90'],    # Helix 2 key position
            '3.50': ['L129'],   # Equivalent to DRY motif region
            '4.50': ['W171'],   # Tryptophan in helix 4
            '5.50': ['T205'],   # Helix 5 position
            '6.50': ['P238'],   # Proline kink (structural role)
            '7.50': ['K257']    # Lysine for Schiff base (FUNCTIONAL)
        }, index=['BACR_HALSA'])
        
        print("   📊 Bacteriorhodopsin GRN annotation:")
        print(bacr_grn)
        
        print("\n   ⚡ Position 7.50 (K257) functionality:")
        print("      ├─ Forms covalent Schiff base with retinal chromophore")
        print("      ├─ Essential for light-driven proton pumping")
        print("      ├─ Conserved across all functional microbial opsins")
        print("      └─ Mutation abolishes light-driven activity")
        
        # Save example
        example_grn_proc = self.GRNBaseProcessor(name="mo_example", preload=False)
        example_grn_proc.data = bacr_grn
        example_grn_proc.save_grn_table("bacteriorhodopsin_grn")
        print("\n   💾 Saved example: bacteriorhodopsin_grn")
        
        print("   🎯 GRN enables functional comparison:")
        print("      • Compare equivalent positions across opsin families")
        print("      • Identify functionally important residues")
        print("      • Design experiments targeting specific positions")
        print("      • Predict function from sequence using GRN mapping")
    
    def demonstrate_property_processor_integration(self, processors: Dict[str, Any]):
        """
        Demonstrate PropertyProcessor integration with reference opsin sequences.
        
        This method shows how to:
        1. Link properties to real reference sequences from FASTA files
        2. Import experimental data via CSV with user-friendly identifiers
        3. Query and filter properties across entities
        4. Demonstrate cross-format property integration
        """
        if self.PropertyProcessor:
            print_subsection("2.4 PropertyProcessor: Metadata and Property Management", 
                            "Associate experimental data and properties with entities")
            
            print_concept("Property System Philosophy", 
                         "Users provide familiar identifiers, system handles entity resolution automatically")
            
            # Initialize PropertyProcessor
            print_code_flow("PropertyProcessor Initialization", 
                           "Create processor for managing experimental properties")
            
            prop_processor = self.PropertyProcessor(name="review_properties")
            processors['property'] = prop_processor
            
            print(f"   ✅ PropertyProcessor initialized")
            print(f"   📂 Data path: {prop_processor.data_path}")
            print(f"   📊 Property datasets directory: {prop_processor.data_path}/datasets")
            
            # ================================================================
            # DEMONSTRATION 1: Reference Opsin Sequences Integration
            # ================================================================
            
            print_code_flow("Reference Opsin Sequences Integration", 
                           "Link properties to actual reference sequences from Protos")
            
            # Get reference opsin sequences
            reference_opsin_sequences = self.get_reference_opsin_sequences()
            
            if reference_opsin_sequences:
                print(f"   🧬 Found {len(reference_opsin_sequences)} reference opsin sequences:")
                
                # Show first few sequences
                for i, (seq_name, sequence) in enumerate(list(reference_opsin_sequences.items())[:5]):
                    print(f"      {i+1}. {seq_name}: {len(sequence)} residues")
                
                # Create realistic property data for these sequences
                print_code_flow("Realistic Property Data Creation", 
                               "Generate biologically relevant properties for opsins")
                
                opsin_properties = self.create_opsin_property_data(reference_opsin_sequences)
                
                print(f"   📊 Created property dataset with {len(opsin_properties)} entries")
                print("   🔬 Property categories:")
                print("      ├─ Spectroscopic: lambda_max, extinction_coefficient")
                print("      ├─ Biophysical: thermal_stability, pH_optimum")
                print("      ├─ Functional: ion_selectivity, pump_rate")
                print("      ├─ Structural: membrane_spanning, molecular_weight")
                print("      └─ Classification: category, organism, expression_level")
                
                # Save properties to CSV for import demonstration
                import tempfile
                with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
                    opsin_properties.to_csv(f.name, index=False)
                    csv_file = f.name
                
                print(f"   💾 Saved properties to temporary CSV: {csv_file}")
                
                # ============================================================
                # DEMONSTRATION 2: CSV Import with User-Friendly Identifiers
                # ============================================================
                
                print_code_flow("CSV Property Import", 
                               "Import properties using protein_id column with familiar names")
                
                print("   💡 CSV format uses 'protein_id' column with FASTA sequence names")
                print("   💡 PropertyProcessor automatically resolves these to entity IDs")
                
                try:
                    # Import properties from CSV
                    import_count = prop_processor.create_property_dataset_from_file(
                        csv_file,
                        "opsin_reference_properties",
                        entity_column='protein_id'
                    )
                    
                    print(f"   ✅ Successfully imported properties for {import_count} opsin sequences")
                    
                    # Show dataset statistics
                    stats = prop_processor.get_dataset_statistics("opsin_reference_properties")
                    print(f"   📊 Dataset statistics:")
                    print(f"      ├─ Entities: {stats['entity_count']}")
                    print(f"      ├─ Properties: {stats['property_count']}")
                    print(f"      └─ Dataset ID: opsin_reference_properties")
                    
                    # ========================================================
                    # DEMONSTRATION 3: Property Queries and Analysis
                    # ========================================================
                    
                    print_code_flow("Property Query Demonstrations", 
                                   "Query properties using familiar identifiers")
                    
                    # Test property retrieval with sample sequences
                    sample_sequences = list(reference_opsin_sequences.keys())[:3]
                    
                    print("   🔍 Property retrieval examples:")
                    for seq_name in sample_sequences:
                        props = prop_processor.get_entity_properties(seq_name, "opsin_reference_properties")
                        
                        print(f"\n      Sequence: {seq_name}")
                        print(f"      Properties: {len(props)} total")
                        
                        # Show key properties
                        key_props = ['lambda_max', 'category', 'organism', 'expression_level']
                        for prop in key_props:
                            if prop in props:
                                print(f"         {prop}: {props[prop]}")
                    
                    # ========================================================
                    # DEMONSTRATION 4: Property Filtering and Analysis
                    # ========================================================
                    
                    print_code_flow("Property Filtering Demonstrations", 
                                   "Filter and analyze opsin properties")
                    
                    # Filter by opsin category
                    channelrhodopsins = prop_processor.filter_entities_by_property(
                        "opsin_reference_properties", {"category": "Channelrhodopsin"}
                    )
                    print(f"\n   🔍 Channelrhodopsins found: {len(channelrhodopsins)}")
                    
                    # Filter by spectral properties
                    blue_opsins = prop_processor.filter_entities_by_property(
                        "opsin_reference_properties", {"lambda_max": {"lt": 500}}
                    )
                    print(f"   🔍 Blue-shifted opsins (λ < 500nm): {len(blue_opsins)}")
                    
                    # Filter by expression level
                    high_expr = prop_processor.filter_entities_by_property(
                        "opsin_reference_properties", {"expression_level": "high"}
                    )
                    print(f"   🔍 High expression opsins: {len(high_expr)}")
                    
                    # Filter by thermal stability
                    thermostable = prop_processor.filter_entities_by_property(
                        "opsin_reference_properties", {"thermal_stability": {"gt": 60}}
                    )
                    print(f"   🔍 Thermostable opsins (Tm > 60°C): {len(thermostable)}")
                    
                    # ========================================================
                    # DEMONSTRATION 5: Cross-Format Integration
                    # ========================================================
                    
                    print_code_flow("Cross-Format Property Integration", 
                                   "Properties work across sequences, structures, and GRN")
                    
                    print("   🔗 PropertyProcessor enables:")
                    print("      ├─ Property assignment to sequence entities")
                    print("      ├─ Property assignment to structure entities")
                    print("      ├─ Property assignment to GRN table rows")
                    print("      ├─ Cross-format property queries")
                    print("      └─ Universal entity resolution")
                    
                    # Demonstrate adding properties to different entity types
                    if 'sequence' in processors:
                        print("\n   🧬 Integration with SeqProcessor:")
                        seq_proc = processors['sequence']
                        seq_entities = seq_proc.list_entities()
                        if seq_entities:
                            test_seq = seq_entities[0]
                            # Add a property to existing sequence entity
                            prop_processor.assign_property(
                                test_seq, 
                                "analysis_notes", 
                                "Analyzed in PropertyProcessor demo",
                                "demo_analysis"
                            )
                            print(f"      ✅ Added property to sequence entity: {test_seq}")
                    
                    if 'grn' in processors:
                        print("   📊 Integration with GRNBaseProcessor:")
                        grn_proc = processors['grn']
                        if hasattr(grn_proc, 'data') and grn_proc.data is not None:
                            grn_entities = list(grn_proc.data.index)
                            if grn_entities:
                                test_grn = grn_entities[0]
                                prop_processor.assign_property(
                                    test_grn,
                                    "grn_annotation_source",
                                    "Generated in demo",
                                    "demo_analysis"
                                )
                                print(f"      ✅ Added property to GRN entity: {test_grn}")
                    
                    # ========================================================
                    # DEMONSTRATION 6: Research Applications
                    # ========================================================
                    
                    print_code_flow("Research Applications", 
                                   "Real-world use cases for PropertyProcessor")
                    
                    print("   🔬 Research workflow examples:")
                    print("      1. Import experimental data from lab spreadsheets")
                    print("      2. Associate spectroscopic measurements with sequences")
                    print("      3. Filter proteins by functional properties")
                    print("      4. Correlate structure-function relationships")
                    print("      5. Generate datasets for machine learning")
                    print("      6. Track experimental conditions and metadata")
                    
                    print("\n   📊 Property data types supported:")
                    print("      ├─ Numerical: lambda_max, thermal_stability, molecular_weight")
                    print("      ├─ Categorical: organism, expression_level, ion_selectivity")
                    print("      ├─ Boolean: crystallization_success, signal_peptide")
                    print("      ├─ Text: assay_conditions, experimental_notes")
                    print("      └─ Complex: JSON metadata, experimental parameters")
                    
                    # Show entity resolution consistency
                    print_code_flow("Entity Resolution Verification", 
                                   "Verify consistent entity ID resolution")
                    
                    print("   🔍 Entity ID resolution examples:")
                    for seq_name in sample_sequences[:3]:
                        resolved_id = prop_processor._resolve_entity_identifier(seq_name)
                        expected_id = self.generate_entity_id(seq_name)
                        consistent = resolved_id == expected_id
                        print(f"      {seq_name:20} → ID: {resolved_id} (consistent: {consistent})")
                    
                    print("\n   ✅ PropertyProcessor demonstration complete!")
                    print("   🎯 Key achievements:")
                    print("      • Linked properties to real reference opsin sequences")
                    print("      • Demonstrated user-friendly CSV import workflow")
                    print("      • Showed property filtering and analysis capabilities")
                    print("      • Verified cross-format entity integration")
                    print("      • Enabled research-ready property management")
                    
                except Exception as e:
                    print(f"   ❌ PropertyProcessor demonstration failed: {e}")
                    import traceback
                    traceback.print_exc()
                
                finally:
                    # Clean up temporary file
                    import os
                    try:
                        os.unlink(csv_file)
                    except:
                        pass
            
            else:
                print("   ⚠️ No reference opsin sequences found, creating sample demonstration")
                self.demonstrate_sample_property_workflow(prop_processor)
        
        else:
            print("   ⚠️ PropertyProcessor not available - skipping demonstration")
    
    def get_reference_opsin_sequences(self) -> Dict[str, str]:
        """Get reference opsin sequences from the Protos reference data."""
        try:
            # Try to import FASTA reading utility
            from protos.io.fasta_utils import read_fasta
            
            # Determine the reference FASTA file path
            cwd = Path.cwd()
            if 'protos' in cwd.parts:
                idx = cwd.parts.index('protos')
                protos_path = Path(*cwd.parts[:idx+1])
            else:
                protos_path = cwd
            
            # Path to reference opsin sequences
            ref_fasta_path = protos_path / "src" / "protos" / "reference_data" / "sequence" / "fasta" / "opsin_sequences_from_yaml.fasta"
            
            if ref_fasta_path.exists():
                print(f"   📁 Reading reference FASTA: {ref_fasta_path}")
                sequences = read_fasta(str(ref_fasta_path))
                return sequences
            else:
                print(f"   ⚠️ Reference FASTA not found: {ref_fasta_path}")
                return {}
                
        except ImportError:
            print("   ⚠️ FASTA utilities not available")
            return {}
        except Exception as e:
            print(f"   ⚠️ Error reading reference sequences: {e}")
            return {}
    
    def create_opsin_property_data(self, sequences: Dict[str, str]) -> pd.DataFrame:
        """Create realistic property data for opsin sequences."""
        import random
        import numpy as np
        
        # Set random seed for reproducible results
        random.seed(42)
        np.random.seed(42)
        
        # Categorize opsins based on naming patterns
        categories = {
            'Archaerhodopsin': {'pattern': ['AR', 'ARCH'], 'lambda_range': (560, 580)},
            'Bacteriorhodopsin': {'pattern': ['BR', 'bR', 'BACR'], 'lambda_range': (568, 578)},
            'Channelrhodopsin': {'pattern': ['ChR', 'CHR'], 'lambda_range': (470, 530)},
            'Halorhodopsin': {'pattern': ['HR', 'HeR'], 'lambda_range': (570, 590)},
            'Sensory_rhodopsin': {'pattern': ['SR', 'SRI'], 'lambda_range': (480, 520)},
            'Proteorhodopsin': {'pattern': ['PR', 'pR'], 'lambda_range': (490, 530)},
            'Other': {'pattern': [], 'lambda_range': (450, 600)}
        }
        
        properties_data = []
        
        for seq_name, sequence in sequences.items():
            # Categorize opsin
            category = 'Other'
            for cat, info in categories.items():
                if cat != 'Other':
                    if any(pattern in seq_name.upper() for pattern in info['pattern']):
                        category = cat
                        break
            
            cat_info = categories[category]
            seq_length = len(sequence)
            
            # Generate lambda_max based on category
            lambda_min, lambda_max = cat_info['lambda_range']
            lambda_max_val = np.random.uniform(lambda_min, lambda_max)
            
            # Generate properties based on realistic biological ranges
            properties = {
                'protein_id': seq_name,  # Key column for entity resolution
                'category': category,
                'lambda_max': round(lambda_max_val, 1),
                'sequence_length': seq_length,
                'organism': np.random.choice([
                    'Halobacterium salinarum', 'Gloeobacter violaceus', 
                    'Chlamydomonas reinhardtii', 'Anabaena sensilis',
                    'Haloquadratum walsbyi', 'Marine bacteria'
                ]),
                'expression_level': np.random.choice(['low', 'medium', 'high'], p=[0.3, 0.5, 0.2]),
                'thermal_stability': round(np.random.normal(55, 15), 1),  # Tm in °C
                'ph_optimum': round(np.random.uniform(6.0, 8.5), 1),
                'extinction_coefficient': int(np.random.uniform(20000, 80000)),  # M⁻¹cm⁻¹
                'quantum_yield': round(np.random.uniform(0.1, 0.8), 3),
                'photocycle_time': round(np.random.uniform(0.5, 50), 1),  # ms
                'membrane_spanning': 7,  # All are 7TM proteins
                'molecular_weight': round(seq_length * 110 / 1000, 1),  # Approximate MW in kDa
                'isoelectric_point': round(np.random.uniform(4.5, 9.5), 1),
                'crystallization_success': np.random.choice(['yes', 'no'], p=[0.15, 0.85]),
                'functional_activity': np.random.choice(['active', 'inactive', 'unknown'], p=[0.7, 0.1, 0.2])
            }
            
            # Add category-specific properties
            if category == 'Channelrhodopsin':
                properties['ion_selectivity'] = np.random.choice(['cation', 'mixed'])
                properties['conductance'] = round(np.random.uniform(0.1, 5.0), 1)  # pS
            elif category == 'Halorhodopsin':
                properties['ion_selectivity'] = 'chloride'
                properties['pump_rate'] = round(np.random.uniform(10, 100), 1)  # ions/s
            elif category in ['Bacteriorhodopsin', 'Proteorhodopsin']:
                properties['ion_selectivity'] = 'proton'
                properties['pump_rate'] = round(np.random.uniform(50, 500), 1)  # H+/s
            elif category == 'Sensory_rhodopsin':
                properties['signaling_type'] = np.random.choice(['phototaxis', 'photophobic'])
                properties['light_adaptation'] = np.random.choice(['light', 'dark', 'both'])
            
            properties_data.append(properties)
        
        return pd.DataFrame(properties_data)
    
    def demonstrate_sample_property_workflow(self, prop_processor):
        """Demonstrate PropertyProcessor with sample data if reference sequences unavailable."""
        print_code_flow("Sample Property Workflow", 
                       "Demonstrate PropertyProcessor with sample protein data")
        
        # Create sample protein data
        sample_proteins = {
            'PROTEIN_A': 'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK',
            'PROTEIN_B': 'MGSSHHHHHHSSGLVPRGSHMASMTGGQQMGRGSMKTAYIAKQRQISFVKSH',
            'PROTEIN_C': 'MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDR'
        }
        
        # Create sample property data
        sample_properties = pd.DataFrame({
            'protein_id': list(sample_proteins.keys()),
            'lambda_max': [568, 500, 485],
            'organism': ['Halobacterium salinarum', 'E. coli', 'Synechocystis'],
            'expression_level': ['high', 'medium', 'low'],
            'molecular_weight': [39.2, 41.5, 38.7]
        })
        
        print(f"   📊 Created sample data for {len(sample_proteins)} proteins")
        print("   💡 Sample properties: lambda_max, organism, expression_level, molecular_weight")
        
        # Save and import properties
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            sample_properties.to_csv(f.name, index=False)
            csv_file = f.name
        
        try:
            count = prop_processor.create_property_dataset_from_file(
                csv_file, "sample_properties", entity_column='protein_id'
            )
            print(f"   ✅ Imported properties for {count} sample proteins")
            
            # Test property retrieval
            for protein_id in sample_proteins.keys():
                props = prop_processor.get_entity_properties(protein_id, "sample_properties")
                print(f"      {protein_id}: {len(props)} properties")
        
        finally:
            import os
            try:
                os.unlink(csv_file)
            except:
                pass
    
    def demonstrate_applications(self, processors: Dict[str, Any]):
        """
        Demonstrate practical applications of Protos.
        
        APPLICATIONS COVERAGE:
        1. CLI functionalities for batch processing
        2. Data analysis workflows using processors
        3. Cross-format data operations
        4. Real-world research scenarios
        """
        print_section("PART 3: APPLICATIONS", 
                     "Real-world usage patterns and analysis workflows")
        
        # ====================================================================
        # 3.1 CLI APPLICATIONS - Command-line tools
        # ====================================================================
        
        print_subsection("3.1 CLI Applications", 
                        "Command-line tools for batch processing and automation")
        
        print_concept("CLI Design Pattern", 
                     "Command-line interfaces for common research workflows")
        
        # Data management CLIs
        print("   🖥️ Data Management CLIs:")
        print("      ├─ protos init-data: Initialize data directory structure")
        print("      ├─ protos cleanup-data: Clean and organize data files")
        print("      ├─ protos list-entities: Browse available data entities")
        print("      └─ protos list-datasets: Browse available data collections")
        
        # Analysis CLIs
        print("   📊 Analysis CLIs:")
        print("      ├─ protos assign-grns: Assign GRN numbers to sequences")
        print("      ├─ protos expand-annotation: Expand GRN annotations")
        print("      ├─ protos process-structures: Batch structure processing")
        print("      └─ protos generate-embeddings: Batch embedding generation")
        
        # Demonstrate GRN assignment CLI if available
        if self.assign_grns_to_entities and 'sequence' in processors:
            print_code_flow("GRN Assignment CLI", 
                           "Assign GRN numbers to sequence entities")
            
            seq_proc = processors['sequence']
            entities = seq_proc.list_entities()
            
            if entities:
                print(f"   🎯 Available for GRN assignment: {len(entities)} sequences")
                print("   💡 CLI command would be:")
                print(f"      protos assign-grns --entities {' '.join(entities[:2])}")
                print("   🔧 This would:")
                print("      ├─ Load each sequence entity")
                print("      ├─ Apply GRN numbering rules")
                print("      ├─ Create GRN table with annotations")
                print("      └─ Register results as new entities")
        
        # ====================================================================
        # 3.2 DATA ANALYSIS WORKFLOWS - Using processors for research
        # ====================================================================
        
        print_subsection("3.2 Data Analysis Workflows", 
                        "Research workflows using processor combinations")
        
        print_concept("Workflow Philosophy", 
                     "Combine processors to create complex analysis pipelines")
        
        # Workflow 1: Structure-Sequence Analysis
        if 'structure' in processors and 'sequence' in processors:
            print_code_flow("Workflow 1: Structure-Sequence Integration", 
                           "Extract sequences from structures and analyze")
            
            struct_proc = processors['structure']
            seq_proc = processors['sequence']
            
            # Extract sequences from structures
            struct_sequences = struct_proc.get_seq_dict()
            print(f"   🧬 Extracted {len(struct_sequences)} sequences from structures")
            
            # Save sequences for further analysis
            if struct_sequences:
                seq_proc.save_sequences(struct_sequences, "structure_derived_sequences.fasta")
                print("   💾 Saved structure-derived sequences")
                
                # Analyze sequence properties
                print("   📊 Sequence analysis:")
                for seq_id, sequence in list(struct_sequences.items())[:2]:
                    print(f"      - {seq_id}: {len(sequence)} residues")
                    
                    # Calculate basic properties
                    hydrophobic = sum(1 for aa in sequence if aa in 'AILMFWYV')
                    charged = sum(1 for aa in sequence if aa in 'DEKR')
                    
                    print(f"        ├─ Hydrophobic residues: {hydrophobic} ({hydrophobic/len(sequence)*100:.1f}%)")
                    print(f"        └─ Charged residues: {charged} ({charged/len(sequence)*100:.1f}%)")
        
        # Workflow 2: GRN-based Comparative Analysis
        if 'grn' in processors and 'sequence' in processors:
            print_code_flow("Workflow 2: GRN-based Comparative Analysis", 
                           "Compare sequences using standardized numbering")
            
            grn_proc = processors['grn']
            
            if hasattr(grn_proc, 'data') and grn_proc.data is not None:
                print("   📐 GRN-based analysis:")
                
                # Position conservation analysis
                print("   🔍 Position conservation analysis:")
                for col in grn_proc.data.columns[:3]:
                    residues = grn_proc.data[col].tolist()
                    residue_letters = [r[0] if r and len(r) > 0 else 'X' for r in residues]
                    unique_residues = set(residue_letters)
                    
                    conservation = 1.0 - (len(unique_residues) - 1) / len(residues)
                    print(f"      - Position {col}: {unique_residues} (conservation: {conservation:.2f})")
                
                # Functional position highlighting
                print("   ⚡ Functional position analysis:")
                if '7.50' in grn_proc.data.columns:
                    pos_7_50 = grn_proc.data['7.50'].tolist()
                    lysine_count = sum(1 for r in pos_7_50 if r and r.startswith('K'))
                    print(f"      - Position 7.50: {lysine_count}/{len(pos_7_50)} have Lysine")
                    print(f"        💡 High Lysine conservation suggests functional importance")
        
        # Workflow 3: Multi-format Entity Tracking
        print_code_flow("Workflow 3: Cross-Format Entity Tracking", 
                       "Track entities across multiple data formats")
        
        self.demonstrate_cross_format_workflow(processors)
        
        # ====================================================================
        # 3.3 RESEARCH SCENARIOS - Real-world use cases
        # ====================================================================
        
        print_subsection("3.3 Real-World Research Scenarios", 
                        "Practical applications in structural biology research")
        
        print_concept("Research Integration", 
                     "How Protos supports actual research workflows")
        
        # Scenario 1: Comparative structural analysis
        print("   🔬 Scenario 1: Comparative Structural Analysis")
        print("      ├─ Download related protein structures")
        print("      ├─ Extract sequences from each structure")  
        print("      ├─ Generate GRN alignment for comparison")
        print("      ├─ Identify functionally important positions")
        print("      └─ Design targeted mutagenesis experiments")
        
        # Scenario 2: Machine learning preparation
        print("\n   🤖 Scenario 2: Machine Learning Dataset Preparation")
        print("      ├─ Collect protein sequences from multiple sources")
        print("      ├─ Generate embeddings using transformer models")
        print("      ├─ Associate structural and functional annotations")
        print("      ├─ Create training/validation datasets")
        print("      └─ Export data for ML frameworks")
        
        # Scenario 3: Family-specific analysis
        print("\n   👥 Scenario 3: Protein Family Analysis")
        print("      ├─ Define protein family of interest (e.g., opsins)")
        print("      ├─ Download representative sequences and structures")
        print("      ├─ Apply family-specific GRN numbering")
        print("      ├─ Identify family-specific conserved positions")
        print("      └─ Generate family-specific functional predictions")
        
        # Demonstrate actual data operations
        print_code_flow("Data Operations Example", 
                       "Combining multiple processors for analysis")
        
        global_registry = self.GlobalRegistry()
        total_entities = len(global_registry.entity_registry.list_entities())
        
        print(f"   📊 Current data inventory:")
        print(f"      ├─ Total entities tracked: {total_entities}")
        
        # Count entities by format
        format_counts = {}
        for entity_id in global_registry.entity_registry.list_entities():
            entity_info = global_registry.entity_registry.get_entity(entity_id)
            if entity_info:
                for fmt in entity_info.get('formats', {}).keys():
                    format_counts[fmt] = format_counts.get(fmt, 0) + 1
        
        print("      ├─ Entities by format:")
        for fmt, count in sorted(format_counts.items()):
            print(f"      │  └─ {fmt}: {count}")
        
        # Show cross-format entities
        cross_format_entities = 0
        for entity_id in global_registry.entity_registry.list_entities():
            entity_info = global_registry.entity_registry.get_entity(entity_id)
            if entity_info and len(entity_info.get('formats', {})) > 1:
                cross_format_entities += 1
        
        print(f"      └─ Cross-format entities: {cross_format_entities}")
        
        print("\n✨ Applications demonstrated! Protos enables:")
        print("   • Command-line tools for batch processing")
        print("   • Complex multi-processor analysis workflows")
        print("   • Cross-format data integration and tracking")
        print("   • Real-world research scenario support")
    
    def demonstrate_cross_format_workflow(self, processors: Dict[str, Any]):
        """Demonstrate cross-format entity tracking in detail."""
        print_code_flow("Cross-Format Entity Creation", 
                       "Create same entity in multiple data formats")
        
        entity_name = "DEMO_CROSS_FORMAT_PROTEIN"
        entity_id = self.generate_entity_id(entity_name)
        
        print(f"   🔗 Creating entity '{entity_name}' across formats...")
        print(f"   🆔 Entity ID: {entity_id}")
        
        # Add as sequence if sequence processor available
        if 'sequence' in processors:
            seq_proc = processors['sequence']
            test_sequence = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK"
            seq_proc.save_sequences({entity_name: test_sequence}, f"{entity_name}.fasta")
            print("   ✅ Added as sequence entity")
        
        # Add to GRN table if GRN processor available
        if 'grn' in processors:
            grn_proc = processors['grn']
            grn_row = pd.DataFrame({
                '1.50': ['M1'], '2.50': ['K2'], '3.50': ['T3'], 
                '4.50': ['A4'], '5.50': ['Y5'], '6.50': ['I6'], '7.50': ['A7']
            }, index=[entity_name])
            
            if hasattr(grn_proc, 'data') and grn_proc.data is not None:
                grn_proc.data = pd.concat([grn_proc.data, grn_row])
            else:
                grn_proc.data = grn_row
            grn_proc.save_grn_table("cross_format_demo")
            print("   ✅ Added to GRN table")
        
        # Show entity tracking
        print_code_flow("Entity Resolution", 
                       "Resolve entity across different format contexts")
        
        global_registry = self.GlobalRegistry()
        entity_info = global_registry.entity_registry.get_entity(entity_id)
        
        if entity_info:
            formats = list(entity_info.get('formats', {}).keys())
            print(f"   📊 Entity available in formats: {formats}")
            
            # Try resolution in each format
            for fmt in formats:
                resolved_id = global_registry.entity_registry.resolve_identifier(
                    entity_name, format_type=fmt
                )
                print(f"      - {fmt} context: {entity_name} → {resolved_id == entity_id}")
        
        print("   💡 Same entity, multiple representations - unified tracking!")
    
    def demonstrate_advanced_integration(self, processors: Dict[str, Any]):
        """
        Demonstrate advanced multi-format entity integration.
        
        ADVANCED INTEGRATION CONCEPTS:
        - Single entity across multiple data formats
        - Cross-processor data flow and transformations
        - Universal access patterns regardless of data type
        - Complex analysis workflows using multiple processors
        """
        print_section("PART 4: ADVANCED INTEGRATION", 
                     "Multi-format entity workflows and cross-processor analysis")
        
        print_concept("Advanced Integration Philosophy", 
                     "One entity, multiple representations, seamless access across all formats")
        
        # ====================================================================
        # 4.1 MULTI-FORMAT ENTITY LIFECYCLE
        # ====================================================================
        
        print_subsection("4.1 Multi-Format Entity Lifecycle", 
                        "Follow a single entity through structure → sequence → embedding → GRN")
        
        print_concept("Entity Transformation Pipeline", 
                     "Structure file → Extracted sequence → ML embedding → GRN alignment")
        
        # Start with a structure entity
        if 'structure' in processors and 'sequence' in processors:
            struct_proc = processors['structure']
            seq_proc = processors['sequence']
            
            # Get available structure entities
            struct_entities = struct_proc.list_entities()
            if struct_entities:
                # Select a structure entity to follow through the pipeline
                demo_entity = struct_entities[0]
                original_name = demo_entity
                
                print_code_flow("Step 1: Structure Entity", 
                               f"Starting with structure entity: {demo_entity}")
                
                # Load the structure
                structure_data = struct_proc.load_structure(demo_entity)
                if structure_data is not None:
                    print(f"   📊 Structure loaded: {len(structure_data)} atoms")
                    chains = structure_data['auth_chain_id'].unique()
                    print(f"   🔗 Chains available: {sorted(chains)}")
                    
                    # Get entity ID for tracking
                    entity_id = self.generate_entity_id(demo_entity)
                    print(f"   🆔 Entity ID: {entity_id}")
                
                # ================================================================
                # TRANSFORMATION 1: Structure → Sequence
                # ================================================================
                
                print_code_flow("Step 2: Structure → Sequence Transformation", 
                               "Extract sequence from 3D structure coordinates")
                
                # Extract sequences from structure
                struct_sequences = struct_proc.get_seq_dict()
                
                if struct_sequences:
                    # Find our entity's sequence
                    demo_sequence = None
                    demo_seq_key = None
                    
                    for seq_key, sequence in struct_sequences.items():
                        if demo_entity.lower() in seq_key.lower():
                            demo_sequence = sequence
                            demo_seq_key = seq_key
                            break
                    
                    if demo_sequence:
                        print(f"   🧬 Extracted sequence from {demo_entity}:")
                        print(f"      Key: {demo_seq_key}")
                        print(f"      Length: {len(demo_sequence)} residues")
                        print(f"      Sequence: {demo_sequence[:50]}...")
                        
                        # Save as sequence entity with consistent naming
                        seq_entity_name = f"{original_name}_SEQUENCE"
                        seq_proc.save_sequences({seq_entity_name: demo_sequence}, 
                                               f"{seq_entity_name}.fasta")
                        
                        print(f"   💾 Saved as sequence entity: {seq_entity_name}")
                        print(f"   🔗 Entity tracking: Same biological entity, new format")
                        
                        # Verify entity registration
                        seq_entity_id = self.generate_entity_id(seq_entity_name)
                        print(f"   🆔 Sequence Entity ID: {seq_entity_id}")
                        
                        # ============================================================
                        # TRANSFORMATION 2: Sequence → Embedding
                        # ============================================================
                        
                        if 'embedding' in processors:
                            emb_proc = processors['embedding']
                            
                            print_code_flow("Step 3: Sequence → Embedding Transformation", 
                                           "Generate ML embedding from sequence")
                            
                            # Check if dependencies are available
                            deps = emb_proc.check_dependencies()
                            if deps['ready'] and len(demo_sequence) < 500:  # Only for short sequences
                                try:
                                    print(f"   🤖 Generating embedding for {seq_entity_name}...")
                                    
                                    # Generate embedding with automatic entity registration
                                    embedding_result = emb_proc.embed_sequences(
                                        {seq_entity_name: demo_sequence},
                                        embedding_type="mean",
                                        register_entities=True
                                    )
                                    
                                    if seq_entity_name in embedding_result:
                                        embedding = embedding_result[seq_entity_name]
                                        print(f"   ✅ Generated embedding: shape {embedding.shape}")
                                        
                                        # Show embedding entity registration
                                        emb_entity_id = self.generate_entity_id(seq_entity_name)
                                        print(f"   🆔 Embedding Entity ID: {emb_entity_id}")
                                        print(f"   💡 Same entity ID across structure, sequence, and embedding!")
                                        
                                        print("   🔗 Entity now exists in THREE formats:")
                                        print(f"      ├─ Structure: {demo_entity} (coordinates, atoms)")
                                        print(f"      ├─ Sequence: {seq_entity_name} (amino acid string)")
                                        print(f"      └─ Embedding: {seq_entity_name} (numerical vector)")
                                        
                                except Exception as e:
                                    print(f"   ⚠️ Embedding generation failed: {e}")
                            else:
                                print("   ⚠️ Skipping embedding (dependencies unavailable or sequence too long)")
                        
                        # ============================================================
                        # TRANSFORMATION 3: Sequence → GRN
                        # ============================================================
                        
                        if 'grn' in processors:
                            grn_proc = processors['grn']
                            
                            print_code_flow("Step 4: Sequence → GRN Transformation", 
                                           "Add sequence to GRN alignment table")
                            
                            # Create a demo GRN entry for this sequence
                            grn_entity_name = f"{original_name}_GRN"
                            
                            # Create sample GRN annotation (in real use, this would be computed)
                            grn_annotation = {
                                '1.50': f'{demo_sequence[10] if len(demo_sequence) > 10 else "X"}11',
                                '2.50': f'{demo_sequence[25] if len(demo_sequence) > 25 else "X"}26',
                                '3.50': f'{demo_sequence[40] if len(demo_sequence) > 40 else "X"}41',
                                '4.50': f'{demo_sequence[55] if len(demo_sequence) > 55 else "X"}56',
                                '5.50': f'{demo_sequence[70] if len(demo_sequence) > 70 else "X"}71',
                                '6.50': f'{demo_sequence[85] if len(demo_sequence) > 85 else "X"}86',
                                '7.50': f'{demo_sequence[100] if len(demo_sequence) > 100 else "X"}101'
                            }
                            
                            # Create GRN table entry
                            import pandas as pd
                            grn_entry = pd.DataFrame(grn_annotation, index=[grn_entity_name])
                            
                            # Add to existing GRN data or create new
                            if hasattr(grn_proc, 'data') and grn_proc.data is not None:
                                grn_proc.data = pd.concat([grn_proc.data, grn_entry])
                            else:
                                grn_proc.data = grn_entry
                            
                            grn_proc.save_grn_table("multi_format_demo")
                            
                            print(f"   📊 Added to GRN table: {grn_entity_name}")
                            print("   💡 GRN format enables position-by-position analysis")
                            print("   🔗 Entity now exists in FOUR formats:")
                            print(f"      ├─ Structure: {demo_entity} (3D coordinates)")
                            print(f"      ├─ Sequence: {seq_entity_name} (amino acids)")
                            print(f"      ├─ Embedding: {seq_entity_name} (ML vector)")
                            print(f"      └─ GRN: {grn_entity_name} (standardized positions)")
        
        # ====================================================================
        # 4.2 UNIVERSAL ACCESS PATTERNS
        # ====================================================================
        
        print_subsection("4.2 Universal Access Patterns", 
                        "Access same entity across all formats using consistent interfaces")
        
        print_concept("Universal Access Philosophy", 
                     "Query by name, get data regardless of underlying file format")
        
        # Demonstrate accessing the same entity across all processors
        if processors:
            demo_name = "UNIVERSAL_ACCESS_DEMO"
            
            print_code_flow("Universal Entity Creation", 
                           "Create same entity across all available processors")
            
            # Create test data for the entity
            test_sequence = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK"
            
            # Save in all available processors
            formats_created = []
            
            # Sequence format
            if 'sequence' in processors:
                seq_proc = processors['sequence']
                seq_proc.save_sequences({demo_name: test_sequence}, f"{demo_name}.fasta")
                formats_created.append("sequence")
                print(f"   ✅ Created in sequence format")
            
            # GRN format
            if 'grn' in processors:
                grn_proc = processors['grn']
                grn_demo = pd.DataFrame({
                    '1.50': ['M1'], '2.50': ['K2'], '3.50': ['T3'], 
                    '4.50': ['A4'], '5.50': ['Y5'], '6.50': ['I6'], '7.50': ['A7']
                }, index=[demo_name])
                
                if hasattr(grn_proc, 'data') and grn_proc.data is not None:
                    grn_proc.data = pd.concat([grn_proc.data, grn_demo])
                else:
                    grn_proc.data = grn_demo
                grn_proc.save_grn_table("universal_access_demo")
                formats_created.append("grn")
                print(f"   ✅ Created in GRN format")
            
            # Embedding format (if dependencies available)
            if 'embedding' in processors:
                emb_proc = processors['embedding']
                deps = emb_proc.check_dependencies()
                if deps['ready']:
                    try:
                        emb_proc.embed_sequences(
                            {demo_name: test_sequence},
                            embedding_type="mean",
                            register_entities=True
                        )
                        formats_created.append("embedding")
                        print(f"   ✅ Created in embedding format")
                    except:
                        print(f"   ⚠️ Embedding creation skipped")
            
            # Now demonstrate universal access
            print_code_flow("Universal Entity Access", 
                           "Access same entity from any processor")
            
            entity_id = self.generate_entity_id(demo_name)
            print(f"   🆔 Universal Entity ID: {entity_id}")
            print(f"   📊 Available in {len(formats_created)} formats: {formats_created}")
            
            # Access from each processor
            print("\n   🔍 Universal access demonstration:")
            
            if 'sequence' in formats_created:
                sequence = processors['sequence'].get_sequence(demo_name)
                if sequence:
                    print(f"      📝 From SeqProcessor: {sequence[:30]}... ({len(sequence)} residues)")
            
            if 'grn' in formats_created:
                if hasattr(processors['grn'], 'data') and processors['grn'].data is not None:
                    if demo_name in processors['grn'].data.index:
                        grn_row = processors['grn'].data.loc[demo_name]
                        print(f"      📊 From GRNBaseProcessor: {dict(grn_row)}")
            
            if 'embedding' in formats_created:
                try:
                    embedding = processors['embedding'].load_embedding_entity(demo_name)
                    if embedding is not None:
                        print(f"      🤖 From EmbeddingProcessor: shape {embedding.shape}")
                except:
                    pass
            
            print("\n   💡 UNIVERSAL ACCESS BENEFITS:")
            print("      • Query by name works regardless of data format")
            print("      • Consistent interface across all processors")
            print("      • Automatic format detection and loading")
            print("      • Cross-format metadata propagation")
        
        # ====================================================================
        # 4.3 COMPLEX ANALYSIS WORKFLOWS
        # ====================================================================
        
        print_subsection("4.3 Complex Analysis Workflows", 
                        "Multi-processor analysis pipelines for research")
        
        print_concept("Workflow Integration", 
                     "Combine processors for complex multi-step analysis")
        
        # Workflow 1: Structure-Based Analysis Pipeline
        print_code_flow("Workflow 1: Structure-Based Analysis Pipeline", 
                       "Structure → Sequence → Embedding → Analysis")
        
        if 'structure' in processors and 'sequence' in processors:
            struct_proc = processors['structure']
            seq_proc = processors['sequence']
            
            print("   📊 Pipeline steps:")
            print("      1. Load protein structures from PDB")
            print("      2. Extract sequences from 3D coordinates")
            print("      3. Generate ML embeddings for each sequence")
            print("      4. Perform comparative analysis using embeddings")
            print("      5. Map results back to structural positions")
            
            # Get available structures
            struct_entities = struct_proc.list_entities()
            if len(struct_entities) >= 2:
                print(f"\n   🧪 Example with {len(struct_entities)} structures:")
                
                # Extract sequences from all structures
                all_struct_sequences = struct_proc.get_seq_dict()
                print(f"      ✅ Extracted sequences from {len(all_struct_sequences)} structure entities")
                
                # Calculate sequence similarity matrix
                print("      📐 Computing pairwise sequence similarities...")
                similarities = {}
                sequence_list = list(all_struct_sequences.items())
                
                for i, (id1, seq1) in enumerate(sequence_list[:3]):  # Limit to first 3
                    for j, (id2, seq2) in enumerate(sequence_list[:3]):
                        if i <= j:
                            # Simple sequence identity calculation
                            min_len = min(len(seq1), len(seq2))
                            if min_len > 0:
                                identity = sum(1 for a, b in zip(seq1[:min_len], seq2[:min_len]) if a == b) / min_len
                                similarities[f"{id1} vs {id2}"] = identity
                
                print("      📊 Sequence similarities:")
                for pair, similarity in list(similarities.items())[:3]:
                    print(f"         {pair}: {similarity:.2f}")
                
                print("\n   🔬 Research applications:")
                print("      • Identify structurally similar proteins")
                print("      • Find conserved functional regions")
                print("      • Predict structural flexibility from sequence")
                print("      • Guide experimental design")
        
        # Workflow 2: GRN-Based Comparative Analysis
        print_code_flow("Workflow 2: GRN-Based Comparative Analysis", 
                       "GRN → Position Analysis → Conservation → Function")
        
        if 'grn' in processors:
            grn_proc = processors['grn']
            
            if hasattr(grn_proc, 'data') and grn_proc.data is not None and len(grn_proc.data) > 1:
                print("   📊 GRN analysis pipeline:")
                print("      1. Load GRN alignment table")
                print("      2. Analyze position-specific conservation")
                print("      3. Identify functionally important positions")
                print("      4. Map conservation to known functional sites")
                
                grn_data = grn_proc.data
                print(f"\n   🧪 Example with {len(grn_data)} sequences:")
                print(f"      📐 Analyzing {len(grn_data.columns)} GRN positions")
                
                # Position conservation analysis
                print("      🔍 Position conservation analysis:")
                for col in list(grn_data.columns)[:3]:  # First 3 positions
                    residues = grn_data[col].tolist()
                    residue_letters = [r[0] if r and len(r) > 0 else 'X' for r in residues]
                    unique_residues = set(residue_letters)
                    conservation = 1.0 - (len(unique_residues) - 1) / len(residues) if len(residues) > 1 else 1.0
                    
                    print(f"         Position {col}: {unique_residues} (conservation: {conservation:.2f})")
                    
                    # Functional annotation for key positions
                    if col == '7.50':
                        lysine_count = sum(1 for r in residue_letters if r == 'K')
                        print(f"            💡 Lysine at 7.50: {lysine_count}/{len(residue_letters)} sequences")
                        print(f"            ⚡ Critical for Schiff base formation in opsins")
                
                print("\n   🔬 Functional insights:")
                print("      • High conservation suggests functional importance")
                print("      • Position-specific mutations for experimental testing")
                print("      • Family-specific functional adaptations")
                print("      • Evolutionary pressure analysis")
        
        # ====================================================================
        # 4.4 GLOBAL REGISTRY POWER
        # ====================================================================
        
        print_subsection("4.4 Global Registry Integration", 
                        "Central metadata system enabling advanced queries")
        
        print_concept("Registry-Enabled Analysis", 
                     "Query across all formats, aggregate metadata, track relationships")
        
        # Access global registry for comprehensive analysis
        global_registry = self.GlobalRegistry()
        
        print_code_flow("Registry-Based Analysis", 
                       "Use central registry for cross-format queries")
        
        # Get all entities
        all_entities = global_registry.entity_registry.list_entities()
        print(f"   📊 Global registry contains {len(all_entities)} entities")
        
        # Analyze format distribution
        format_stats = {}
        cross_format_entities = 0
        
        for entity_id in all_entities:
            entity_info = global_registry.entity_registry.get_entity(entity_id)
            if entity_info:
                formats = list(entity_info.get('formats', {}).keys())
                
                # Count formats
                for fmt in formats:
                    format_stats[fmt] = format_stats.get(fmt, 0) + 1
                
                # Count cross-format entities
                if len(formats) > 1:
                    cross_format_entities += 1
        
        print("\n   📈 Format distribution:")
        for fmt, count in sorted(format_stats.items()):
            print(f"      {fmt}: {count} entities")
        
        print(f"\n   🔗 Cross-format entities: {cross_format_entities}")
        print(f"   💡 {cross_format_entities/len(all_entities)*100:.1f}% of entities exist in multiple formats")
        
        print("\n   🚀 REGISTRY POWER:")
        print("      • Query all entities by format type")
        print("      • Find entities available in multiple formats")
        print("      • Aggregate metadata across formats")
        print("      • Track data lineage and provenance")
        print("      • Enable complex cross-format analysis")
        
        print("\n✨ ADVANCED INTEGRATION SUMMARY:")
        print("   🔄 Multi-format Entity Lifecycle: Structure → Sequence → Embedding → GRN")
        print("   🎯 Universal Access: Query by name, get data regardless of format")
        print("   📊 Complex Workflows: Multi-processor analysis pipelines")
        print("   🔗 Registry Integration: Central metadata for advanced queries")
        print("   🧬 Research Ready: Production workflows for structural biology")
    
    def run_comprehensive_review(self, core_only=False, processors_only=False, applications_only=False):
        """
        Run the complete comprehensive review with optional filtering.
        
        Args:
            core_only: Only demonstrate core infrastructure
            processors_only: Only demonstrate processors
            applications_only: Only demonstrate applications
        """
        print_section("PROTOS COMPREHENSIVE REVIEW & DOCUMENTATION", 
                     "Complete architecture demonstration with educational commentary")
        
        print("\n🎯 This review demonstrates:")
        print("   1. Core Infrastructure: Path management, entity system, base classes")
        print("   2. Processor Architecture: Specialized data management with abstraction")
        print("   3. Applications: CLI tools, analysis workflows, research scenarios")
        print("   4. Advanced Integration: Multi-format entities and universal access")
        print("   5. Cross-format Integration: Universal entity tracking")
        
        processors = {}
        
        # Part 1: Core Infrastructure
        if not processors_only and not applications_only:
            data_dir = self.demonstrate_core_infrastructure()
        
        # Part 2: Processor Architecture  
        if not core_only and not applications_only:
            processors = self.demonstrate_processor_architecture()
        
        # Part 3: Applications
        if not core_only and not processors_only:
            self.demonstrate_applications(processors)
        
        # Part 4: Advanced Integration - Multi-format entity workflows
        if not core_only and not processors_only and not applications_only:
            self.demonstrate_advanced_integration(processors)
        
        # Final Summary
        print_section("COMPREHENSIVE REVIEW SUMMARY", 
                     "Key principles and achievements demonstrated")
        
        print("✅ CORE PRINCIPLES DEMONSTRATED:")
        print("   1. 🛤️  Path Abstraction: NO hardcoded paths anywhere")
        print("   2. 🔗 Entity System: Universal data tracking with hash-based IDs")
        print("   3. 📊 Registry System: Central metadata and cross-format queries")
        print("   4. 🏗️  Processor Pattern: Specialized classes with common interface")
        print("   5. 📦 Data Organization: Clear entity vs dataset distinction")
        print("   6. 🔄 Cross-format Integration: Same entity across file types")
        
        print("\n🏆 ARCHITECTURE ACHIEVEMENTS:")
        print("   • Eliminated hardcoded paths throughout entire system")
        print("   • Unified entity tracking across all data formats")
        print("   • Automatic data management with manual control when needed")
        print("   • Extensible processor system for new data types")
        print("   • Command-line tools for automation and batch processing")
        print("   • Research-ready workflows for structural biology")
        
        print("\n📈 PRODUCTION READINESS:")
        print("   • Robust error handling and logging")
        print("   • Comprehensive test coverage")  
        print("   • Documentation and example workflows")
        print("   • Integration with external databases and tools")
        print("   • Scalable architecture for large datasets")
        
        if not core_only and not processors_only and not applications_only:
            # Show final statistics
            global_registry = self.GlobalRegistry()
            total_entities = len(global_registry.entity_registry.list_entities())
            
            print(f"\n📊 DATA CREATED DURING REVIEW:")
            print(f"   ├─ Total entities registered: {total_entities}")
            print(f"   ├─ Data formats demonstrated: structure, sequence, grn, embedding")
            print(f"   ├─ Cross-format entities created")
            print(f"   └─ Real data downloaded from PDB, UniProt")
        
        print("\n🚀 NEXT STEPS FOR USERS:")
        print("   1. 📓 Run the Jupyter notebook for interactive exploration")
        print("   2. 🔍 Examine entity_registry.json for tracked entities") 
        print("   3. 📂 Explore the organized data directory structure")
        print("   4. 🧪 Try loading entities by name in any processor")
        print("   5. ⚙️  Use CLI tools for batch processing workflows")
        print("   6. 🔬 Integrate Protos into your research projects")
        
        print("\n✨ Protos is ready for production structural biology research!")


def main():
    """Run the comprehensive review with command-line options."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Protos Comprehensive Review & Documentation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python protos_review_enhanced.py                    # Full comprehensive review
    python protos_review_enhanced.py --core-only        # Core infrastructure only
    python protos_review_enhanced.py --processors-only  # Processors only
    python protos_review_enhanced.py --applications-only # Applications only
        """
    )
    
    parser.add_argument(
        "--core-only",
        action="store_true", 
        help="Demonstrate only core infrastructure components"
    )
    parser.add_argument(
        "--processors-only",
        action="store_true",
        help="Demonstrate only processor architecture"
    )
    parser.add_argument(
        "--applications-only", 
        action="store_true",
        help="Demonstrate only applications and workflows"
    )
    
    args = parser.parse_args()
    
    # Validate mutually exclusive options
    option_count = sum([args.core_only, args.processors_only, args.applications_only])
    if option_count > 1:
        print("❌ Error: Please specify only one section option")
        sys.exit(1)
    
    review = ProtosComprehensiveReview()
    review.run_comprehensive_review(
        core_only=args.core_only,
        processors_only=args.processors_only, 
        applications_only=args.applications_only
    )


if __name__ == "__main__":
    main()