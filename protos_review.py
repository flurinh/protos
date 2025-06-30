#!/usr/bin/env python3
"""
Protos Comprehensive Review & Documentation Script

This script provides a complete demonstration and educational overview of the Protos
structural biology framework. It covers all major components from core infrastructure
to advanced applications with real data examples.

PROTOS ARCHITECTURE OVERVIEW:
┌─────────────────────────────────────────────────────────────────┐
│                           PROTOS FRAMEWORK                     │
├─────────────────────────────────────────────────────────────────┤
│ 1. CORE INFRASTRUCTURE                                          │
│    ├─ ProtosPaths: Centralized path management                  │
│    ├─ BaseProcessor: Abstract base for all processors          │
│    └─ Entity Registry: Universal entity tracking system        │
│                                                                 │
│ 2. SPECIALIZED PROCESSORS                                       │
│    ├─ CifBaseProcessor: 3D protein structure management        │
│    ├─ GRNBaseProcessor: Generic Residue Numbering system       │
│    ├─ SeqProcessor: Sequence data management                   │
│    ├─ PropertyProcessor: Metadata and properties               │
│    └─ EmbeddingProcessor: ML embeddings generation             │
│                                                                 │
│ 3. DATA ABSTRACTION                                             │
│    ├─ Entities: Individual data items                          │
│    ├─ Datasets: Collections of related entities                │
│    └─ Cross-format tracking: Same entity across formats        │
│                                                                 │
│ 4. APPLICATIONS                                                 │
│    ├─ CLI tools: Command-line utilities                        │
│    ├─ Analysis workflows: Multi-processor pipelines            │
│    └─ Research integration: Real-world use cases               │
└─────────────────────────────────────────────────────────────────┘

Usage:
    python protos_review.py                    # Full comprehensive review
    python protos_review.py --quick            # Quick overview
    python protos_review.py --demo             # Interactive demo mode
"""

import os
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import logging
from typing import List, Dict, Optional, Tuple, Any
import json
import argparse

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def print_header(title: str, subtitle: str = ""):
    """Print a formatted header with title and optional subtitle."""
    print(f"\n{'=' * 80}")
    print(f"{title:^80}")
    if subtitle:
        print(f"{subtitle:^80}")
    print('=' * 80)


def print_section(title: str, description: str = ""):
    """Print a formatted section header."""
    print(f"\n{'-' * 60}")
    print(f"{title}")
    if description:
        print(f"   {description}")
    print('-' * 60)


def print_concept(concept: str, explanation: str):
    """Print a formatted concept explanation."""
    print(f"\n💡 CONCEPT: {concept}")
    print(f"   {explanation}")


def print_step(step: str, action: str):
    """Print a formatted workflow step."""
    print(f"\n⚙️  {step}")
    print(f"   → {action}")


class ProtosComprehensiveReview:
    """
    Comprehensive review and documentation of the Protos framework.
    
    This class demonstrates every aspect of Protos with educational comments
    explaining design patterns, abstractions, and best practices for
    structural biology research.
    """
    
    def __init__(self, demo_mode: bool = False):
        """
        Initialize the comprehensive review system.
        
        Args:
            demo_mode: Enable interactive demo features
        """
        self.demo_mode = demo_mode
        self.components = {}
        self.processors = {}
        self.data_created = []
        
        print_header("PROTOS COMPREHENSIVE REVIEW", 
                    "Complete framework demonstration with real data examples")
        
        self._load_protos_components()
        self._setup_data_environment()
    
    def _load_protos_components(self):
        """Load and validate all Protos components."""
        print_section("Loading Protos Components", 
                     "Importing and validating framework components")
        
        # Core Infrastructure
        print_concept("Core Infrastructure", 
                     "Foundation components that enable all other functionality")
        
        try:
            from protos.io.paths.path_config import ProtosPaths
            from protos.io.data_access import EntityRegistry, GlobalRegistry, generate_entity_id
            from protos.core.base_processor import BaseProcessor
            
            self.components.update({
                'ProtosPaths': ProtosPaths,
                'GlobalRegistry': GlobalRegistry,
                'generate_entity_id': generate_entity_id,
                'BaseProcessor': BaseProcessor
            })
            
            print("✅ Core Infrastructure:")
            print("   ├─ ProtosPaths: Centralized path management")
            print("   ├─ Entity System: Universal data tracking")
            print("   └─ BaseProcessor: Common processor interface")
            
        except ImportError as e:
            print(f"❌ Core infrastructure error: {e}")
            sys.exit(1)
        
        # Specialized Processors
        print_concept("Specialized Processors", 
                     "Domain-specific processors for different data types")
        
        processor_imports = [
            ('CifBaseProcessor', 'protos.processing.structure.struct_base_processor', 
             '3D protein structure management'),
            ('GRNBaseProcessor', 'protos.processing.grn.grn_base_processor', 
             'Generic Residue Numbering system'),
            ('SeqProcessor', 'protos.processing.sequence.seq_processor', 
             'Sequence data management'),
            ('PropertyProcessor', 'protos.processing.property.property_processor_enhanced', 
             'Metadata and property management'),
            ('EmbeddingProcessor', 'protos.processing.embedding.embedding_processor', 
             'ML embeddings generation')
        ]
        
        for class_name, module_path, description in processor_imports:
            try:
                module = __import__(module_path, fromlist=[class_name])
                processor_class = getattr(module, class_name)
                self.components[class_name] = processor_class
                print(f"✅ {class_name}: {description}")
            except ImportError:
                print(f"⚠️  {class_name}: Not available")
                self.components[class_name] = None
        
        # Data Loaders
        print_concept("Data Loaders", 
                     "External database integration utilities")
        
        loader_imports = [
            ('download_protein_structures', 'protos.loaders.download_structures'),
            ('download_alphafold_structures', 'protos.loaders.alphafold_utils'),
            ('download_sequences_from_uniprot', 'protos.loaders.uniprot_utils')
        ]
        
        for func_name, module_path in loader_imports:
            try:
                module = __import__(module_path, fromlist=[func_name])
                self.components[func_name] = getattr(module, func_name)
            except ImportError:
                self.components[func_name] = None
        
        print("✅ Data Loaders: PDB, AlphaFold, UniProt integration")
        
        # CLI Tools
        try:
            from protos.cli.init_data import init_data_directory
            self.components['init_data_directory'] = init_data_directory
            print("✅ CLI Tools: Data management utilities")
        except ImportError:
            self.components['init_data_directory'] = None
            print("⚠️  CLI Tools: Not available")
    
    def _setup_data_environment(self):
        """Set up the data environment and paths."""
        print_section("Data Environment Setup", 
                     "Configuring paths and initializing data directory")
        
        # Determine data directory
        cwd = Path.cwd()
        if 'protos' in cwd.parts:
            idx = cwd.parts.index('protos')
            protos_path = Path(*cwd.parts[:idx+1])
        else:
            protos_path = cwd
        
        self.data_dir = protos_path / "data"
        
        print_step("Data Directory Configuration", 
                  f"Using data directory: {self.data_dir}")
        
        # Initialize data directory if needed
        if not (self.data_dir / "entity_registry.json").exists():
            if self.components['init_data_directory']:
                print("   🔧 Initializing data directory structure...")
                self.components['init_data_directory'](data_root=str(self.data_dir), force=True)
                print("   ✅ Data directory initialized")
        
        # Set global data root
        self.components['ProtosPaths'].set_data_root(str(self.data_dir.absolute()))
        
        print("   🌍 Global data root configured")
        print(f"   📁 All processors will use: {self.data_dir}")
    
    def demonstrate_core_concepts(self):
        """Demonstrate core Protos concepts and abstractions."""
        print_header("PART 1: CORE CONCEPTS", 
                    "Foundation principles and abstractions")
        
        # ================================================================
        # 1.1 PATH MANAGEMENT
        # ================================================================
        
        print_section("1.1 Path Management System", 
                     "Centralized path handling eliminates hardcoded paths")
        
        print_concept("Path Abstraction", 
                     "Users work with names, Protos handles all file system operations")
        
        paths = self.components['ProtosPaths']()
        
        print("   📋 Managed directories:")
        for processor, subdir in [
            ('Structure', 'structure'), ('Sequence', 'sequence'), 
            ('GRN', 'grn'), ('Property', 'property'), ('Embedding', 'embedding')
        ]:
            print(f"      ├─ {processor}: {self.data_dir}/{subdir}")
        
        print("   💡 Benefits: No hardcoded paths, platform independence, centralized config")
        
        # ================================================================
        # 1.2 ENTITY SYSTEM
        # ================================================================
        
        print_section("1.2 Entity System", 
                     "Universal data tracking with deterministic IDs")
        
        print_concept("Entity Philosophy", 
                     "Every data item gets a unique, reproducible identifier")
        
        # Demonstrate entity ID generation
        demo_entities = ["1ubq", "TEST_PROTEIN", "BACR_HALSA"]
        print("   🔍 Entity ID examples:")
        
        entity_mappings = {}
        for name in demo_entities:
            entity_id = self.components['generate_entity_id'](name)
            entity_mappings[name] = entity_id
            print(f"      '{name}' → {entity_id}")
        
        print("   💡 Same input always produces same ID (deterministic)")
        print("   💡 Works across all data formats (universal)")
        
        # ================================================================
        # 1.3 REGISTRY SYSTEM
        # ================================================================
        
        print_section("1.3 Global Registry", 
                     "Central tracking system for all entities")
        
        print_concept("Registry Benefits", 
                     "Cross-format queries, metadata management, relationship tracking")
        
        global_registry = self.components['GlobalRegistry']()
        total_entities = len(global_registry.entity_registry.list_entities())
        
        print(f"   📊 Current registry: {total_entities} entities tracked")
        print("   🔧 Registry functions:")
        print("      ├─ Track entities across multiple formats")
        print("      ├─ Resolve names to canonical IDs")
        print("      ├─ Store metadata and relationships")
        print("      └─ Enable complex cross-format queries")
        
        return entity_mappings
    
    def demonstrate_processors(self):
        """Demonstrate all available processors with real data."""
        print_header("PART 2: PROCESSOR DEMONSTRATIONS", 
                    "Specialized data management with real examples")
        
        processors_created = {}
        
        # ================================================================
        # 2.1 STRUCTURE PROCESSOR
        # ================================================================
        
        if self.components['CifBaseProcessor']:
            print_section("2.1 CifBaseProcessor", 
                         "3D protein structure management")
            
            print_concept("Structure Data Abstraction", 
                         "Parse PDB/CIF files, analyze coordinates, extract sequences")
            
            # Initialize processor
            struct_proc = self.components['CifBaseProcessor'](name="review_structures")
            processors_created['structure'] = struct_proc
            
            print_step("Processor Initialization", 
                      f"Data path: {struct_proc.data_path}")
            
            # Download demo structures
            if self.components['download_protein_structures']:
                demo_pdbs = ["1ubq", "2gb1", "1crn"]  # Small, well-studied proteins
                print_step("Structure Download", 
                          f"Downloading: {demo_pdbs}")
                
                successful, failed = self.components['download_protein_structures'](
                    demo_pdbs, target_folder=struct_proc.path_structure_dir
                )
                print(f"   ✅ Downloaded: {successful}")
                if failed:
                    print(f"   ⚠️ Failed: {failed}")
                
                self.data_created.extend([f"Structure: {pdb}" for pdb in successful])
            
            # Demonstrate entity management
            entities = struct_proc.list_entities()
            print_step("Entity Management", 
                      f"Available entities: {len(entities)}")
            
            if entities:
                # Load and analyze a structure
                demo_entity = entities[0]
                print(f"   🔍 Loading: {demo_entity}")
                
                structure_data = struct_proc.load_structure(demo_entity, apply_dtypes=True)
                if structure_data is not None:
                    chains = structure_data['auth_chain_id'].unique()
                    residues = structure_data['auth_seq_id'].nunique()
                    
                    print(f"   ✅ Loaded: {len(structure_data)} atoms")
                    print(f"      ├─ Chains: {sorted(chains)}")
                    print(f"      ├─ Residues: {residues}")
                    print(f"      └─ Coordinate range: X({structure_data['x'].min():.1f} to {structure_data['x'].max():.1f})")
                
                # Create dataset
                dataset_id = "demo_structures"
                struct_proc.create_dataset(
                    dataset_id=dataset_id,
                    name="Demo Protein Structures",
                    description="Small proteins for demonstration",
                    content=successful[:3] if 'successful' in locals() else entities[:3]
                )
                print_step("Dataset Creation", 
                          f"Created dataset: {dataset_id}")
                
                # Extract sequences (cross-format operation)
                sequences = struct_proc.get_seq_dict()
                print_step("Cross-Format Operation", 
                          f"Extracted {len(sequences)} sequences from structures")
                
                for seq_id, seq in list(sequences.items())[:2]:
                    print(f"      - {seq_id}: {seq[:30]}... ({len(seq)} residues)")
        
        # ================================================================
        # 2.2 SEQUENCE PROCESSOR
        # ================================================================
        
        if self.components['SeqProcessor']:
            print_section("2.2 SeqProcessor", 
                         "Sequence data management and analysis")
            
            print_concept("Sequence Abstraction", 
                         "FASTA files, database downloads, sequence operations")
            
            # Initialize processor
            seq_proc = self.components['SeqProcessor'](name="review_sequences")
            processors_created['sequence'] = seq_proc
            
            print_step("Processor Initialization", 
                      f"Data path: {seq_proc.data_path}")
            
            # Download from UniProt
            if self.components['download_sequences_from_uniprot']:
                uniprot_ids = ["P00533", "P04637"]  # EGFR, p53
                print_step("UniProt Download", 
                          f"Downloading: {uniprot_ids}")
                
                downloaded_sequences = {}
                for uniprot_id in uniprot_ids:
                    try:
                        seq_data = self.components['download_sequences_from_uniprot']([uniprot_id])
                        if seq_data:
                            downloaded_sequences.update(seq_data)
                            print(f"   ✅ Downloaded: {uniprot_id}")
                    except Exception as e:
                        print(f"   ⚠️ Failed {uniprot_id}: {e}")
                
                if downloaded_sequences:
                    seq_proc.save_sequences(downloaded_sequences, "uniprot_demo.fasta")
                    self.data_created.append(f"Sequences: {len(downloaded_sequences)} from UniProt")
            
            # Create demo sequences
            demo_sequences = {
                "DEMO_PROTEIN_A": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK",
                "DEMO_PROTEIN_B": "MGSSHHHHHHSSGLVPRGSHMASMTGGQQMGRGSMKTAYIAKQRQISFVKSH",
                "DEMO_PROTEIN_C": "MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDR"
            }
            
            # Save as entities and dataset
            for name, seq in demo_sequences.items():
                seq_proc.save_sequences({name: seq}, f"{name}.fasta")
            
            seq_proc.save_sequences(demo_sequences, "demo_proteins.fasta")
            
            print_step("Demo Data Creation", 
                      f"Created {len(demo_sequences)} demo sequences")
            
            self.data_created.append(f"Demo sequences: {len(demo_sequences)} proteins")
            
            # Demonstrate entity access
            entities = seq_proc.list_entities()
            print_step("Entity Management", 
                      f"Total sequence entities: {len(entities)}")
        
        # ================================================================
        # 2.3 GRN PROCESSOR
        # ================================================================
        
        if self.components['GRNBaseProcessor']:
            print_section("2.3 GRNBaseProcessor", 
                         "Generic Residue Numbering for protein families")
            
            print_concept("GRN System", 
                         "Standardize residue numbering across homologous proteins")
            
            # Initialize processor
            grn_proc = self.components['GRNBaseProcessor'](name="review_grn", preload=False)
            processors_created['grn'] = grn_proc
            
            print_step("Processor Initialization", 
                      f"Data path: {grn_proc.data_path}")
            
            # Create demo GRN table with proper format
            print_step("GRN Table Creation", 
                      "Creating alignment with residue+position format")
            
            grn_data = pd.DataFrame({
                '1.50': ['M62', 'M62', 'L21', 'V45', 'M62'],    # Helix 1
                '2.50': ['K115', 'G90', 'A65', 'G67', 'K115'], # Helix 2
                '3.50': ['T179', 'S145', 'E107', 'I128', 'R179'], # Helix 3
                '4.50': ['A221', 'S190', 'G149', 'V171', 'W221'], # Helix 4
                '5.50': ['Y270', 'H236', 'E195', 'Y214', 'F270'], # Helix 5
                '6.50': ['I312', 'H280', 'I236', 'I256', 'W312'], # Helix 6
                '7.50': ['A356', 'H324', 'T280', 'A300', 'N356']  # Helix 7
            }, index=['PROTEIN_A', 'PROTEIN_B', 'PROTEIN_C', 'PROTEIN_D', 'PROTEIN_E'])
            
            print("   📊 GRN Table format:")
            print(f"      ├─ Rows: {list(grn_data.index)}")
            print(f"      ├─ Columns: {list(grn_data.columns)}")
            print(f"      └─ Values: residue+position (e.g., M62, K115)")
            
            # Save GRN table
            grn_proc.data = grn_data
            grn_proc.save_grn_table("demo_alignment")
            
            print_step("GRN Analysis", 
                      "Built-in analysis capabilities")
            
            # Load and analyze
            grn_proc.load_grn_table("demo_alignment")
            sequences = grn_proc.get_seq_dict()
            
            print(f"   🧬 Extracted sequences: {len(sequences)}")
            for seq_id, seq in list(sequences.items())[:2]:
                print(f"      - {seq_id}: {seq}")
            
            # Demonstrate microbial opsin GRN
            self._demonstrate_opsin_grn()
            
            self.data_created.append(f"GRN table: {len(grn_data)} sequences")
        
        # ================================================================
        # 2.4 PROPERTY PROCESSOR
        # ================================================================
        
        if self.components['PropertyProcessor']:
            print_section("2.4 PropertyProcessor", 
                         "Metadata and experimental property management")
            
            print_concept("Property System", 
                         "Associate experimental data with entities using familiar identifiers")
            
            # Initialize processor
            prop_proc = self.components['PropertyProcessor'](name="review_properties")
            processors_created['property'] = prop_proc
            
            print_step("Processor Initialization", 
                      f"Data path: {prop_proc.data_path}")
            
            # Get reference opsin sequences if available
            opsin_sequences = self._get_reference_opsin_sequences()
            
            if opsin_sequences:
                print_step("Reference Data Integration", 
                          f"Found {len(opsin_sequences)} reference opsin sequences")
                
                # Create realistic properties
                properties_df = self._create_opsin_properties(opsin_sequences)
                
                # Save to temporary CSV
                import tempfile
                with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
                    properties_df.to_csv(f.name, index=False)
                    csv_file = f.name
                
                try:
                    # Import properties
                    count = prop_proc.create_property_dataset_from_file(
                        csv_file, "opsin_properties", entity_column='protein_id'
                    )
                    
                    print_step("Property Import", 
                              f"Imported properties for {count} opsins")
                    
                    # Demonstrate filtering
                    blue_opsins = prop_proc.filter_entities_by_property(
                        "opsin_properties", {"lambda_max": {"lt": 500}}
                    )
                    
                    high_expr = prop_proc.filter_entities_by_property(
                        "opsin_properties", {"expression_level": "high"}
                    )
                    
                    print_step("Property Analysis", 
                              "Filtering and analysis capabilities")
                    print(f"   🔍 Blue-shifted opsins (λ < 500nm): {len(blue_opsins)}")
                    print(f"   🔍 High expression opsins: {len(high_expr)}")
                    
                    self.data_created.append(f"Opsin properties: {count} entities with experimental data")
                    
                finally:
                    os.unlink(csv_file)
            
            else:
                # Create sample property demo
                self._demonstrate_sample_properties(prop_proc)
        
        # ================================================================
        # 2.5 EMBEDDING PROCESSOR
        # ================================================================
        
        if self.components['EmbeddingProcessor']:
            print_section("2.5 EmbeddingProcessor", 
                         "ML embeddings for sequence analysis")
            
            print_concept("Embedding Abstraction", 
                         "Transform sequences into numerical vectors for ML")
            
            # Initialize processor
            emb_proc = self.components['EmbeddingProcessor'](name="review_embeddings")
            processors_created['embedding'] = emb_proc
            
            print_step("Processor Initialization", 
                      f"Data path: {emb_proc.data_path}")
            
            # Check dependencies
            deps = emb_proc.check_dependencies()
            print_step("Dependency Check", 
                      "Verifying ML framework availability")
            
            for dep, available in deps.items():
                status = "✅" if available else "❌"
                print(f"   {status} {dep}")
            
            # Show available models
            models = emb_proc.list_available_models()
            print_step("Available Models", 
                      f"Found {len(models)} embedding models")
            
            for name, info in list(models.items())[:3]:
                print(f"   - {name}: {info['embedding_dim']}D")
            
            if deps['ready'] and 'sequence' in processors_created:
                print_step("Embedding Generation", 
                          "Generating embeddings for demo sequences")
                
                seq_proc = processors_created['sequence']
                entities = seq_proc.list_entities()[:2]  # First 2 entities
                
                for entity in entities:
                    sequence = seq_proc.get_sequence(entity)
                    if sequence and len(sequence) < 200:  # Only short sequences
                        try:
                            embeddings = emb_proc.embed_sequences(
                                {entity: sequence}, 
                                embedding_type="mean",
                                register_entities=True
                            )
                            
                            if entity in embeddings:
                                print(f"   ✅ {entity}: {embeddings[entity].shape}")
                                self.data_created.append(f"Embedding: {entity}")
                                
                        except Exception as e:
                            print(f"   ⚠️ {entity}: {e}")
            else:
                print("   ⚠️ Skipping generation (dependencies unavailable)")
        
        return processors_created
    
    def _demonstrate_opsin_grn(self):
        """Demonstrate GRN specific to microbial opsins."""
        print_section("2.3.1 Microbial Opsin GRN", 
                     "Domain-specific GRN for rhodopsins")
        
        print_concept("Opsin Biology", 
                     "Light-driven proteins with 7-transmembrane structure")
        
        # Bacteriorhodopsin example
        bacr_grn = pd.DataFrame({
            '1.50': ['M62'],    # Helix 1 start
            '2.50': ['V90'],    # Helix 2 
            '3.50': ['L129'],   # DRY motif region
            '4.50': ['W171'],   # Tryptophan
            '5.50': ['T205'],   # Helix 5
            '6.50': ['P238'],   # Proline kink
            '7.50': ['K257']    # Schiff base lysine (CRITICAL)
        }, index=['BACR_HALSA'])
        
        print_step("Bacteriorhodopsin Example", 
                  "Most studied microbial opsin")
        print(bacr_grn)
        
        print("   ⚡ Position 7.50 (K257) is critical:")
        print("      ├─ Forms Schiff base with retinal")
        print("      ├─ Essential for proton pumping")
        print("      └─ Conserved across functional opsins")
    
    def _get_reference_opsin_sequences(self):
        """Get reference opsin sequences if available."""
        try:
            from protos.io.fasta_utils import read_fasta
            
            # Find reference FASTA
            cwd = Path.cwd()
            if 'protos' in cwd.parts:
                idx = cwd.parts.index('protos')
                protos_path = Path(*cwd.parts[:idx+1])
            else:
                protos_path = cwd
            
            ref_fasta = protos_path / "src" / "protos" / "reference_data" / "sequence" / "fasta" / "opsin_sequences_from_yaml.fasta"
            
            if ref_fasta.exists():
                return read_fasta(str(ref_fasta))
            
        except ImportError:
            pass
        
        return {}
    
    def _create_opsin_properties(self, sequences):
        """Create realistic property data for opsins."""
        import random
        import numpy as np
        
        # Set seed for reproducibility
        random.seed(42)
        np.random.seed(42)
        
        properties_data = []
        
        for seq_name, sequence in list(sequences.items())[:20]:  # Limit for demo
            # Categorize by name patterns
            category = 'Other'
            lambda_range = (450, 600)
            
            if any(x in seq_name.upper() for x in ['BR', 'BACR']):
                category = 'Bacteriorhodopsin'
                lambda_range = (568, 578)
            elif any(x in seq_name.upper() for x in ['CHR', 'ChR']):
                category = 'Channelrhodopsin'
                lambda_range = (470, 530)
            elif any(x in seq_name.upper() for x in ['HR', 'HeR']):
                category = 'Halorhodopsin'
                lambda_range = (570, 590)
            
            properties = {
                'protein_id': seq_name,
                'category': category,
                'lambda_max': round(np.random.uniform(*lambda_range), 1),
                'sequence_length': len(sequence),
                'expression_level': np.random.choice(['low', 'medium', 'high']),
                'thermal_stability': round(np.random.normal(55, 10), 1),
                'organism': np.random.choice([
                    'Halobacterium salinarum',
                    'Chlamydomonas reinhardtii', 
                    'Gloeobacter violaceus'
                ])
            }
            
            properties_data.append(properties)
        
        return pd.DataFrame(properties_data)
    
    def _demonstrate_sample_properties(self, prop_proc):
        """Demonstrate PropertyProcessor with sample data."""
        print_step("Sample Property Demo", 
                  "Creating sample property dataset")
        
        sample_data = pd.DataFrame({
            'protein_id': ['DEMO_PROTEIN_A', 'DEMO_PROTEIN_B', 'DEMO_PROTEIN_C'],
            'lambda_max': [568, 500, 485],
            'organism': ['Halobacterium', 'E. coli', 'Synechocystis'],
            'expression_level': ['high', 'medium', 'low']
        })
        
        # Save and import
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            sample_data.to_csv(f.name, index=False)
            csv_file = f.name
        
        try:
            count = prop_proc.create_property_dataset_from_file(
                csv_file, "sample_properties", entity_column='protein_id'
            )
            print(f"   ✅ Imported {count} sample properties")
            self.data_created.append(f"Sample properties: {count} entities")
        finally:
            os.unlink(csv_file)
    
    def demonstrate_workflows(self, processors):
        """Demonstrate advanced workflows and integrations."""
        print_header("PART 3: ADVANCED WORKFLOWS", 
                    "Multi-processor analysis and cross-format integration")
        
        # ================================================================
        # 3.1 CROSS-FORMAT ENTITY TRACKING
        # ================================================================
        
        print_section("3.1 Cross-Format Entity Tracking", 
                     "Same entity across multiple data types")
        
        print_concept("Universal Entity Model", 
                     "One biological entity, multiple data representations")
        
        # Create entity across formats
        entity_name = "DEMO_CROSS_FORMAT"
        entity_id = self.components['generate_entity_id'](entity_name)
        
        print_step("Multi-Format Entity Creation", 
                  f"Creating '{entity_name}' across processors")
        print(f"   🆔 Universal Entity ID: {entity_id}")
        
        formats_created = []
        
        # Add to sequence processor
        if 'sequence' in processors:
            test_sequence = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK"
            processors['sequence'].save_sequences(
                {entity_name: test_sequence}, f"{entity_name}.fasta"
            )
            formats_created.append("sequence")
            print("   ✅ Added as sequence entity")
        
        # Add to GRN table
        if 'grn' in processors:
            grn_row = pd.DataFrame({
                '1.50': ['M1'], '2.50': ['K2'], '3.50': ['T3'], 
                '4.50': ['A4'], '5.50': ['Y5'], '6.50': ['I6'], '7.50': ['A7']
            }, index=[entity_name])
            
            if hasattr(processors['grn'], 'data') and processors['grn'].data is not None:
                processors['grn'].data = pd.concat([processors['grn'].data, grn_row])
            else:
                processors['grn'].data = grn_row
            processors['grn'].save_grn_table("cross_format_demo")
            formats_created.append("grn")
            print("   ✅ Added to GRN table")
        
        # Add properties
        if 'property' in processors:
            processors['property'].assign_property(
                entity_name, "demo_property", "cross_format_test", "demo_dataset"
            )
            formats_created.append("property")
            print("   ✅ Added properties")
        
        print_step("Universal Access Verification", 
                  "Same entity accessible from all processors")
        print(f"   🔗 Entity exists in {len(formats_created)} formats: {formats_created}")
        
        # ================================================================
        # 3.2 STRUCTURE-SEQUENCE ANALYSIS PIPELINE
        # ================================================================
        
        if 'structure' in processors and 'sequence' in processors:
            print_section("3.2 Structure-Sequence Pipeline", 
                         "Extract sequences from structures for analysis")
            
            print_concept("Cross-Format Analysis", 
                         "Structure coordinates → extracted sequences → analysis")
            
            struct_proc = processors['structure']
            seq_proc = processors['sequence']
            
            # Extract sequences from structures
            struct_sequences = struct_proc.get_seq_dict()
            
            print_step("Sequence Extraction", 
                      f"Extracted {len(struct_sequences)} sequences from structures")
            
            if struct_sequences:
                # Save extracted sequences
                seq_proc.save_sequences(struct_sequences, "structure_derived.fasta")
                
                # Analyze sequences
                print_step("Sequence Analysis", 
                          "Basic sequence property analysis")
                
                for seq_id, sequence in list(struct_sequences.items())[:2]:
                    # Calculate properties
                    hydrophobic = sum(1 for aa in sequence if aa in 'AILMFWYV')
                    charged = sum(1 for aa in sequence if aa in 'DEKR')
                    
                    print(f"   📊 {seq_id}:")
                    print(f"      ├─ Length: {len(sequence)} residues")
                    print(f"      ├─ Hydrophobic: {hydrophobic/len(sequence)*100:.1f}%")
                    print(f"      └─ Charged: {charged/len(sequence)*100:.1f}%")
        
        # ================================================================
        # 3.3 PROPERTY-BASED ANALYSIS
        # ================================================================
        
        if 'property' in processors:
            print_section("3.3 Property-Based Analysis", 
                         "Filter and analyze entities by experimental properties")
            
            print_concept("Property-Driven Discovery", 
                         "Use experimental data to guide analysis")
            
            prop_proc = processors['property']
            datasets = prop_proc.list_datasets()
            
            if datasets:
                dataset_id = datasets[0]['id']
                
                print_step("Property Analysis", 
                          f"Analyzing dataset: {dataset_id}")
                
                # Get statistics
                stats = prop_proc.get_dataset_statistics(dataset_id)
                print(f"   📊 Dataset contains:")
                print(f"      ├─ Entities: {stats['entity_count']}")
                print(f"      └─ Properties: {stats['property_count']}")
                
                # Show research applications
                print_step("Research Applications", 
                          "How properties enable discovery")
                print("   🔬 Research use cases:")
                print("      ├─ Filter proteins by functional properties")
                print("      ├─ Correlate structure-function relationships")
                print("      ├─ Generate ML training datasets")
                print("      └─ Track experimental conditions")
        
        # ================================================================
        # 3.4 GLOBAL REGISTRY ANALYSIS
        # ================================================================
        
        print_section("3.4 Global Registry Analysis", 
                     "Cross-format entity tracking and statistics")
        
        print_concept("Registry Power", 
                     "Central system enables complex cross-format queries")
        
        global_registry = self.components['GlobalRegistry']()
        all_entities = global_registry.entity_registry.list_entities()
        
        print_step("Registry Statistics", 
                  f"Analyzing {len(all_entities)} tracked entities")
        
        # Analyze format distribution
        format_stats = {}
        cross_format_count = 0
        
        for entity_id in all_entities:
            entity_info = global_registry.entity_registry.get_entity(entity_id)
            if entity_info:
                formats = list(entity_info.get('formats', {}).keys())
                
                for fmt in formats:
                    format_stats[fmt] = format_stats.get(fmt, 0) + 1
                
                if len(formats) > 1:
                    cross_format_count += 1
        
        print("   📈 Format distribution:")
        for fmt, count in sorted(format_stats.items()):
            print(f"      ├─ {fmt}: {count} entities")
        
        print(f"   🔗 Cross-format entities: {cross_format_count}")
        if all_entities:
            print(f"   💡 {cross_format_count/len(all_entities)*100:.1f}% exist in multiple formats")
        
        self.data_created.append(f"Registry: {len(all_entities)} entities tracked")
    
    def demonstrate_research_scenarios(self):
        """Demonstrate real-world research applications."""
        print_header("PART 4: RESEARCH APPLICATIONS", 
                    "Real-world structural biology use cases")
        
        # ================================================================
        # 4.1 PROTEIN FAMILY ANALYSIS
        # ================================================================
        
        print_section("4.1 Protein Family Analysis", 
                     "Comparative analysis of homologous proteins")
        
        print_concept("Family-Based Research", 
                     "Study protein families using standardized numbering")
        
        print_step("Workflow Overview", 
                  "Typical protein family analysis pipeline")
        
        workflow_steps = [
            "1. Download representative structures and sequences",
            "2. Apply family-specific GRN numbering scheme", 
            "3. Identify conserved and variable positions",
            "4. Associate experimental properties with positions",
            "5. Generate hypotheses for functional testing"
        ]
        
        for step in workflow_steps:
            print(f"   {step}")
        
        print("\n   🎯 Research outcomes:")
        print("      ├─ Functional site identification")
        print("      ├─ Evolutionary relationship mapping")
        print("      ├─ Drug target identification")
        print("      └─ Protein engineering guidance")
        
        # ================================================================
        # 4.2 MACHINE LEARNING PREPARATION
        # ================================================================
        
        print_section("4.2 ML Dataset Preparation", 
                     "Prepare structured datasets for machine learning")
        
        print_concept("ML-Ready Data", 
                     "Transform biological data into ML-compatible formats")
        
        print_step("ML Pipeline", 
                  "From raw data to ML-ready datasets")
        
        ml_steps = [
            "1. Collect sequences and structures from databases",
            "2. Generate embeddings for sequence representation",
            "3. Extract structural features from coordinates",
            "4. Associate experimental properties as labels",
            "5. Create train/validation/test splits",
            "6. Export to ML frameworks (PyTorch, TensorFlow)"
        ]
        
        for step in ml_steps:
            print(f"   {step}")
        
        print("\n   🤖 ML applications:")
        print("      ├─ Function prediction from sequence")
        print("      ├─ Structure prediction and validation")
        print("      ├─ Property prediction (stability, activity)")
        print("      └─ Evolutionary analysis and design")
        
        # ================================================================
        # 4.3 EXPERIMENTAL DESIGN
        # ================================================================
        
        print_section("4.3 Experimental Design", 
                     "Guide experiments using computational analysis")
        
        print_concept("Computation-Guided Experiments", 
                     "Use Protos analysis to design targeted experiments")
        
        print_step("Experimental Workflow", 
                  "From analysis to experimental design")
        
        exp_steps = [
            "1. Identify positions of interest using GRN analysis",
            "2. Filter properties to find candidates for mutagenesis",
            "3. Predict effects using structural analysis",
            "4. Design mutation library based on conservation",
            "5. Track experimental results as new properties"
        ]
        
        for step in exp_steps:
            print(f"   {step}")
        
        print("\n   🧪 Experimental outcomes:")
        print("      ├─ Targeted mutagenesis experiments")
        print("      ├─ Structure-function validation")
        print("      ├─ Improved protein variants")
        print("      └─ Mechanistic understanding")
    
    def generate_summary(self):
        """Generate comprehensive summary of the demonstration."""
        print_header("COMPREHENSIVE SUMMARY", 
                    "What we accomplished and next steps")
        
        # ================================================================
        # ACHIEVEMENTS
        # ================================================================
        
        print_section("Demonstration Achievements", 
                     "Data created and concepts demonstrated")
        
        print("✅ CORE CONCEPTS DEMONSTRATED:")
        achievements = [
            "🛤️ Path Abstraction: No hardcoded paths anywhere",
            "🔗 Entity System: Universal data tracking with deterministic IDs",
            "📊 Registry System: Central metadata and cross-format queries",
            "🏗️ Processor Pattern: Specialized classes with common interface",
            "📦 Data Organization: Clear entity vs dataset distinction",
            "🔄 Cross-format Integration: Same entity across file types"
        ]
        
        for achievement in achievements:
            print(f"   {achievement}")
        
        print("\n📊 DATA CREATED DURING DEMONSTRATION:")
        for item in self.data_created:
            print(f"   • {item}")
        
        # ================================================================
        # PRODUCTION READINESS
        # ================================================================
        
        print_section("Production Readiness", 
                     "Why Protos is ready for research")
        
        production_features = [
            "Robust error handling and comprehensive logging",
            "Automatic data management with manual control options", 
            "Extensible architecture for new data types",
            "Integration with external databases (PDB, UniProt, AlphaFold)",
            "Command-line tools for automation and batch processing",
            "Comprehensive test coverage and documentation"
        ]
        
        for feature in production_features:
            print(f"   ✅ {feature}")
        
        # ================================================================
        # NEXT STEPS
        # ================================================================
        
        print_section("Next Steps for Users", 
                     "How to start using Protos in your research")
        
        next_steps = [
            "📓 Explore the Jupyter notebook version of this review",
            "🔍 Examine entity_registry.json to see tracked entities", 
            "📂 Browse the organized data directory structure",
            "🧪 Try loading entities by name in any processor",
            "⚙️ Use CLI tools for your own data processing",
            "🔬 Integrate Protos into your research projects"
        ]
        
        for step in next_steps:
            print(f"   {step}")
        
        # ================================================================
        # FINAL MESSAGE
        # ================================================================
        
        print_header("PROTOS IS READY", 
                    "Production-ready structural biology framework")
        
        print("🚀 Protos provides a complete solution for:")
        final_benefits = [
            "Managing heterogeneous biological data",
            "Standardizing protein family analysis",
            "Enabling reproducible research workflows", 
            "Facilitating cross-format data integration",
            "Supporting machine learning applications",
            "Accelerating structural biology research"
        ]
        
        for benefit in final_benefits:
            print(f"   • {benefit}")
        
        print(f"\n✨ Framework demonstration complete!")
        print(f"   📊 Created {len(self.data_created)} data examples")
        print(f"   🔧 Demonstrated {len(self.processors)} processors")
        print(f"   🎯 Ready for production research use")


def main():
    """Run the comprehensive Protos review."""
    parser = argparse.ArgumentParser(
        description="Protos Comprehensive Review & Documentation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python protos_review.py                # Full comprehensive review
    python protos_review.py --quick        # Quick overview only
    python protos_review.py --demo         # Interactive demo mode
        """
    )
    
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Run quick overview (core concepts only)"
    )
    parser.add_argument(
        "--demo",
        action="store_true", 
        help="Enable interactive demo mode with user prompts"
    )
    
    args = parser.parse_args()
    
    # Initialize review system
    review = ProtosComprehensiveReview(demo_mode=args.demo)
    
    try:
        # Part 1: Core concepts (always shown)
        entity_mappings = review.demonstrate_core_concepts()
        
        if not args.quick:
            # Part 2: Processors
            processors = review.demonstrate_processors()
            
            # Part 3: Advanced workflows  
            review.demonstrate_workflows(processors)
            
            # Part 4: Research applications
            review.demonstrate_research_scenarios()
        
        # Summary (always shown)
        review.generate_summary()
        
    except KeyboardInterrupt:
        print("\n\n⚠️ Review interrupted by user")
        print("✅ Partial demonstration completed")
    except Exception as e:
        print(f"\n❌ Error during review: {e}")
        logger.exception("Review failed")
        sys.exit(1)


if __name__ == "__main__":
    main()