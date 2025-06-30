#!/usr/bin/env python3
"""
Real-World PropertyProcessor Demo: User-Friendly Entity Integration

This demo shows how users can assign properties using familiar identifiers
(FASTA names, PDB IDs, file names) rather than entity IDs. The PropertyProcessor
automatically resolves these to the correct entity IDs behind the scenes.

Workflow:
1. User has a CSV with 'protein_id' column containing familiar names
2. PropertyProcessor resolves these names to entity IDs
3. Properties are correctly associated with entities across formats
4. Integration works seamlessly with existing processors

Example CSV structure:
protein_id,lambda_max,organism,expression_level
1ubq,568,Halobacterium salinarum,high
TEST_PROTEIN,500,Gloeobacter violaceus,medium
rhodopsin_sequence,485,Anabaena sensilis,low
"""

import os
import sys
from pathlib import Path
import pandas as pd
import tempfile

# Add protos to path
protos_path = Path(__file__).parent.parent
sys.path.insert(0, str(protos_path))

# Import Protos components
try:
    from protos.io.paths.path_config import ProtosPaths
    from protos.processing.property.property_processor_enhanced import PropertyProcessor
    from protos.processing.structure.struct_base_processor import CifBaseProcessor
    from protos.processing.sequence.seq_processor import SeqProcessor
    from protos.processing.grn.grn_base_processor import GRNBaseProcessor
    from protos.io.data_access import generate_entity_id, EntityRegistry
    print("✅ All imports successful!")
except ImportError as e:
    print(f"❌ Import error: {e}")
    print("Make sure protos is installed: pip install -e .")
    sys.exit(1)


def setup_demo_environment():
    """Set up demo environment with real data."""
    print("=== Setting Up Real-World Demo Environment ===")
    
    # Set up data directory
    demo_data_dir = protos_path / "data"
    ProtosPaths.set_data_root(str(demo_data_dir.absolute()))
    
    print(f"Demo data directory: {demo_data_dir}")
    return demo_data_dir


def create_sample_property_tables():
    """Create sample property tables that users would typically have."""
    print("\n=== Creating Sample Property Tables ===")
    
    # 1. Microbial Opsin Properties (typical experimental data)
    microbial_opsin_data = pd.DataFrame({
        'protein_id': [
            '1ubq',           # Structure identifier
            'TEST_PROTEIN',   # Sequence identifier  
            'OPSIN_001',      # Another sequence
            'rhodopsin_bovine', # Descriptive name
            'bacteriorhodopsin_1' # Another descriptive name
        ],
        'lambda_max': [568, 500, 485, 500, 568],
        'organism': [
            'Halobacterium salinarum',
            'Gloeobacter violaceus', 
            'Anabaena sensilis',
            'Bos taurus',
            'Halobacterium salinarum'
        ],
        'expression_system': ['native', 'E.coli', 'E.coli', 'native', 'native'],
        'expression_level': ['high', 'medium', 'low', 'high', 'high'],
        'temperature_optimum': [37.0, 25.0, 30.0, 37.0, 45.0],
        'pH_optimum': [7.0, 7.5, 7.2, 7.4, 6.8],
        'method': ['spectroscopy', 'spectroscopy', 'spectroscopy', 'literature', 'spectroscopy']
    })
    
    # 2. Structural Properties (from PDB/experimental data)
    structural_data = pd.DataFrame({
        'protein_id': [
            '1ubq',
            '2abc', 
            'rhodopsin_bovine'
        ],
        'resolution': [1.55, 2.1, 2.8],
        'space_group': ['P212121', 'P1', 'P321'],
        'crystal_system': ['orthorhombic', 'triclinic', 'hexagonal'],
        'r_factor': [0.198, 0.245, 0.223],
        'method': ['X-ray', 'X-ray', 'X-ray'],
        'ligand_bound': ['retinal', 'none', 'retinal']
    })
    
    # 3. Functional Properties (biochemical assays)
    functional_data = pd.DataFrame({
        'protein_id': [
            '1ubq',
            'TEST_PROTEIN',
            'rhodopsin_bovine',
            'bacteriorhodopsin_1'
        ],
        'kcat': [125.5, 89.2, 156.8, 134.2],  # turnover number
        'km': [2.5, 4.1, 1.8, 3.2],           # Michaelis constant
        'binding_affinity': [8.5, 7.2, 9.1, 8.8],  # pKd
        'thermal_stability': [65.5, 45.2, 58.9, 72.1], # Tm in °C
        'assay_conditions': ['pH 7.0, 25°C', 'pH 7.5, 25°C', 'pH 7.4, 37°C', 'pH 6.8, 45°C']
    })
    
    # 4. Evolutionary/Conservation Data
    evolutionary_data = pd.DataFrame({
        'protein_id': [
            '1ubq',
            'TEST_PROTEIN', 
            'OPSIN_001',
            'rhodopsin_bovine',
            'bacteriorhodopsin_1'
        ],
        'conservation_score': [0.89, 0.76, 0.82, 0.91, 0.87],
        'phylogenetic_group': ['archaeal', 'cyanobacterial', 'cyanobacterial', 'vertebrate', 'archaeal'],
        'domain_architecture': ['7TM', '7TM', '7TM', '7TM', '7TM'],
        'signal_peptide': [False, False, True, True, False],
        'transmembrane_domains': [7, 7, 7, 7, 7]
    })
    
    print(f"✅ Created microbial opsin data: {len(microbial_opsin_data)} entries")
    print(f"✅ Created structural data: {len(structural_data)} entries") 
    print(f"✅ Created functional data: {len(functional_data)} entries")
    print(f"✅ Created evolutionary data: {len(evolutionary_data)} entries")
    
    return {
        'microbial_opsins': microbial_opsin_data,
        'structural_properties': structural_data,
        'functional_assays': functional_data,
        'evolutionary_analysis': evolutionary_data
    }


def save_property_tables_to_files(property_tables, temp_dir):
    """Save property tables to CSV files."""
    print("\n=== Saving Property Tables to CSV Files ===")
    
    file_paths = {}
    for table_name, df in property_tables.items():
        file_path = temp_dir / f"{table_name}.csv"
        df.to_csv(file_path, index=False)
        file_paths[table_name] = file_path
        print(f"💾 Saved {table_name} to {file_path}")
        
        # Show sample of data
        print(f"   Sample data (first 3 rows):")
        print(f"   {df.head(3).to_string(index=False)}")
        print()
    
    return file_paths


def test_entity_id_resolution():
    """Test that protein_id values resolve to correct entity IDs."""
    print("\n=== Testing Entity ID Resolution ===")
    
    # Test identifiers that users would typically use
    test_identifiers = [
        '1ubq',           # PDB ID
        'TEST_PROTEIN',   # Sequence name from FASTA
        'rhodopsin_bovine', # Descriptive protein name
        'bacteriorhodopsin_1', # Another descriptive name
        'OPSIN_001'       # Generic sequence identifier
    ]
    
    print("🔍 Testing entity ID generation consistency:")
    entity_mappings = {}
    
    for identifier in test_identifiers:
        entity_id = generate_entity_id(identifier)
        entity_mappings[identifier] = entity_id
        print(f"   {identifier:20} → {entity_id}")
    
    # Test that same identifier always generates same entity ID
    print("\n🔍 Testing consistency (same input → same ID):")
    for identifier in test_identifiers:
        id1 = generate_entity_id(identifier)
        id2 = generate_entity_id(identifier) 
        id3 = generate_entity_id(identifier)
        consistent = (id1 == id2 == id3)
        print(f"   {identifier:20} → Consistent: {consistent}")
        assert consistent, f"Entity ID generation not consistent for {identifier}"
    
    return entity_mappings


def test_property_assignment_with_real_identifiers(property_processor, property_tables):
    """Test property assignment using real user identifiers."""
    print("\n=== Testing Property Assignment with Real Identifiers ===")
    
    # Test with microbial opsin data first
    microbial_data = property_tables['microbial_opsins']
    
    print("🧪 Assigning properties using protein_id column...")
    
    assignments = []
    for _, row in microbial_data.iterrows():
        protein_id = row['protein_id']
        
        # Create assignments for each property column (except protein_id)
        for column in microbial_data.columns:
            if column != 'protein_id':
                assignments.append({
                    'entity_identifier': protein_id,
                    'property_name': column,
                    'property_value': row[column],
                    'metadata': {'source': 'experimental_data', 'table': 'microbial_opsins'}
                })
    
    # Batch assign all properties
    entity_ids = property_processor.assign_properties_batch(
        assignments, "microbial_opsins_real"
    )
    
    print(f"✅ Assigned {len(assignments)} properties to {len(set(entity_ids))} entities")
    
    # Verify assignments worked
    print("\n🔍 Verifying property assignments:")
    for protein_id in microbial_data['protein_id'].unique():
        entity_props = property_processor.get_entity_properties(protein_id, "microbial_opsins_real")
        print(f"   {protein_id:20} → {len(entity_props)} properties")
        
        # Show sample properties
        if entity_props:
            sample_props = dict(list(entity_props.items())[:3])  # First 3 properties
            for prop_name, prop_value in sample_props.items():
                print(f"     {prop_name}: {prop_value}")
    
    return entity_ids


def test_cross_format_property_integration(property_processor, entity_mappings):
    """Test that properties work across different processor formats."""
    print("\n=== Testing Cross-Format Property Integration ===")
    
    # Simulate that some identifiers correspond to different format types
    format_examples = {
        '1ubq': 'structure',          # This would be in structure processor
        'TEST_PROTEIN': 'sequence',   # This would be in sequence processor  
        'rhodopsin_bovine': 'sequence', # This could be in sequence processor
        'OPSIN_001': 'grn'           # This might be in GRN processor
    }
    
    print("🔗 Testing property assignment across different format types:")
    
    for identifier, format_type in format_examples.items():
        # Assign format-specific properties
        format_specific_props = {
            'structure': [('resolution', 2.1), ('space_group', 'P21'), ('ligand_count', 1)],
            'sequence': [('length', 280), ('molecular_weight', 39000), ('isoelectric_point', 5.8)],
            'grn': [('conserved_positions', 15), ('motif_count', 3), ('alignment_score', 0.95)]
        }
        
        props_to_assign = format_specific_props.get(format_type, [])
        
        for prop_name, prop_value in props_to_assign:
            property_processor.assign_property(
                identifier, prop_name, prop_value, 
                f"{format_type}_properties",
                metadata={'format_type': format_type, 'source': 'processor_integration'}
            )
        
        print(f"   ✅ {identifier} ({format_type}): {len(props_to_assign)} properties assigned")
    
    # Verify cross-format integration
    print("\n🔍 Verifying cross-format property access:")
    
    test_entity = '1ubq'  # This entity should have properties from multiple datasets
    all_props = property_processor.get_entity_properties(test_entity)
    
    print(f"   Entity '{test_entity}' total properties: {len(all_props)}")
    print(f"   Properties across all datasets:")
    for prop_name, prop_value in all_props.items():
        print(f"     {prop_name}: {prop_value}")
    
    # Show which datasets this entity appears in
    datasets = property_processor.list_datasets()
    entity_datasets = []
    for dataset in datasets:
        entities = property_processor.list_entities(dataset['id'])
        entity_id = generate_entity_id(test_entity)
        if entity_id in entities:
            entity_datasets.append(dataset['id'])
    
    print(f"\n   Entity '{test_entity}' appears in datasets: {entity_datasets}")


def test_csv_import_workflow(property_processor, file_paths):
    """Test the complete CSV import workflow."""
    print("\n=== Testing CSV Import Workflow ===")
    
    # Import each CSV file as a dataset
    import_results = {}
    
    for table_name, file_path in file_paths.items():
        print(f"\n📥 Importing {table_name} from {file_path}")
        
        try:
            count = property_processor.create_property_dataset_from_file(
                str(file_path), 
                f"imported_{table_name}",
                entity_column='protein_id'
            )
            
            import_results[table_name] = count
            print(f"   ✅ Successfully imported {count} entities")
            
            # Show dataset statistics
            stats = property_processor.get_dataset_statistics(f"imported_{table_name}")
            print(f"   📊 Dataset stats: {stats['entity_count']} entities, {stats['property_count']} properties")
            
        except Exception as e:
            print(f"   ❌ Import failed: {e}")
            import_results[table_name] = 0
    
    return import_results


def test_property_filtering_with_real_data(property_processor):
    """Test property filtering with real-world queries."""
    print("\n=== Testing Property Filtering with Real Data ===")
    
    # Test filtering on imported datasets
    test_queries = [
        {
            'dataset': 'imported_microbial_opsins',
            'filters': {'lambda_max': {'gt': 520}},
            'description': 'Opsins with lambda_max > 520 nm'
        },
        {
            'dataset': 'imported_microbial_opsins', 
            'filters': {'organism': 'Halobacterium salinarum'},
            'description': 'Halobacterium proteins'
        },
        {
            'dataset': 'imported_structural_properties',
            'filters': {'resolution': {'lt': 2.0}},
            'description': 'High-resolution structures (< 2.0 Å)'
        },
        {
            'dataset': 'imported_functional_assays',
            'filters': {'thermal_stability': {'gt': 60.0}},
            'description': 'Thermostable proteins (Tm > 60°C)'
        }
    ]
    
    print("🔍 Running real-world property queries:")
    
    for query in test_queries:
        dataset_name = query['dataset']
        filters = query['filters']
        description = query['description']
        
        try:
            matching_entities = property_processor.filter_entities_by_property(
                dataset_name, filters
            )
            
            print(f"\n   Query: {description}")
            print(f"   Dataset: {dataset_name}")
            print(f"   Filters: {filters}")
            print(f"   Results: {len(matching_entities)} entities")
            
            if matching_entities:
                print(f"   Matching entities: {matching_entities}")
                
                # Show properties of first matching entity
                first_entity = matching_entities[0]
                props = property_processor.get_entity_properties(first_entity, dataset_name)
                print(f"   Sample properties for {first_entity}:")
                for prop_name, prop_value in list(props.items())[:3]:
                    print(f"     {prop_name}: {prop_value}")
            
        except Exception as e:
            print(f"   ❌ Query failed: {e}")


def test_entity_resolution_integration():
    """Test integration with actual processor entity resolution."""
    print("\n=== Testing Entity Resolution Integration ===")
    
    print("🔍 Testing resolution of common identifiers:")
    
    # Test identifiers that might exist in different processors
    test_cases = [
        ('1ubq', 'structure', 'PDB structure identifier'),
        ('TEST_PROTEIN', 'sequence', 'FASTA sequence name'),
        ('rhodopsin_bovine', None, 'Generic protein name'),
        ('OPSIN_001', 'sequence', 'Generic sequence identifier')
    ]
    
    prop_proc = PropertyProcessor(name="resolution_test")
    
    for identifier, expected_format, description in test_cases:
        print(f"\n   Testing: {identifier} ({description})")
        
        try:
            # Test resolution through PropertyProcessor
            resolved_id = prop_proc._resolve_entity_identifier(identifier, expected_format)
            expected_id = generate_entity_id(identifier)
            
            print(f"     Resolved ID: {resolved_id}")
            print(f"     Expected ID: {expected_id}")
            print(f"     Match: {resolved_id == expected_id}")
            
            # Test property assignment with resolution
            prop_proc.assign_property(
                identifier, 'test_property', 'test_value',
                'resolution_test',
                metadata={'identifier_type': description}
            )
            
            # Verify we can retrieve using both identifier and resolved ID
            value1 = prop_proc.get_entity_property(identifier, 'test_property', 'resolution_test')
            value2 = prop_proc.get_entity_property(resolved_id, 'test_property', 'resolution_test')
            
            print(f"     Property retrieval by identifier: {value1}")
            print(f"     Property retrieval by entity ID: {value2}")
            print(f"     Values match: {value1 == value2}")
            
        except Exception as e:
            print(f"     ❌ Resolution failed: {e}")


def demonstrate_complete_workflow():
    """Demonstrate the complete real-world workflow."""
    print("\n=== Complete Real-World Workflow Demonstration ===")
    
    # Initialize PropertyProcessor
    prop_proc = PropertyProcessor(name="real_world_demo")
    
    # Create temporary directory for CSV files
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        
        # 1. Create sample property tables (what users would have)
        property_tables = create_sample_property_tables()
        
        # 2. Save to CSV files (user's workflow)
        file_paths = save_property_tables_to_files(property_tables, temp_path)
        
        # 3. Import CSV files into PropertyProcessor
        import_results = test_csv_import_workflow(prop_proc, file_paths)
        
        print(f"\n📊 Import Summary:")
        total_imported = sum(import_results.values())
        print(f"   Total entities imported: {total_imported}")
        print(f"   Datasets created: {len(import_results)}")
        
        # 4. Test property filtering and analysis
        test_property_filtering_with_real_data(prop_proc)
        
        # 5. Show final statistics
        print(f"\n📈 Final Statistics:")
        all_datasets = prop_proc.list_datasets()
        all_entities = prop_proc.list_entities()
        
        print(f"   Total datasets: {len(all_datasets)}")
        print(f"   Total unique entities: {len(all_entities)}")
        
        for dataset in all_datasets:
            print(f"   - {dataset['id']}: {dataset['entity_count']} entities, {dataset['property_count']} properties")
        
        return prop_proc


def main():
    """Run the complete real-world demonstration."""
    print("=" * 80)
    print("PROPERTY PROCESSOR REAL-WORLD DEMO")
    print("User-Friendly Entity Integration")
    print("=" * 80)
    
    # Set up environment
    demo_data_dir = setup_demo_environment()
    
    # Test entity ID resolution
    entity_mappings = test_entity_id_resolution()
    
    # Create sample data
    property_tables = create_sample_property_tables()
    
    # Test property assignment with real identifiers
    prop_proc = PropertyProcessor(name="real_world_test")
    entity_ids = test_property_assignment_with_real_identifiers(prop_proc, property_tables)
    
    # Test cross-format integration
    test_cross_format_property_integration(prop_proc, entity_mappings)
    
    # Test entity resolution integration
    test_entity_resolution_integration()
    
    # Demonstrate complete workflow
    final_processor = demonstrate_complete_workflow()
    
    # Final summary
    print("\n" + "=" * 80)
    print("REAL-WORLD DEMO SUMMARY")
    print("=" * 80)
    
    print(f"\n✅ Successfully demonstrated real-world PropertyProcessor usage:")
    print(f"   🔑 User-friendly identifiers: PDB IDs, FASTA names, descriptive names")
    print(f"   🔧 Automatic entity ID resolution behind the scenes")
    print(f"   📊 CSV import workflow with 'protein_id' column")
    print(f"   🔗 Cross-format property integration")
    print(f"   🔍 Property filtering and analysis")
    print(f"   📈 Comprehensive statistics and reporting")
    
    print(f"\n🎯 Key User Benefits:")
    benefits = [
        "✅ Use familiar identifiers (no need to know entity IDs)",
        "✅ Import properties from standard CSV files", 
        "✅ Automatic resolution to correct entities",
        "✅ Properties work across all processor types",
        "✅ Powerful filtering and analysis capabilities",
        "✅ Seamless integration with existing Protos workflows"
    ]
    
    for benefit in benefits:
        print(f"   {benefit}")
    
    print(f"\n🚀 PropertyProcessor is ready for real-world use!")


if __name__ == "__main__":
    main()