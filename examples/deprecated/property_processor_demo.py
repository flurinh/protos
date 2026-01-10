#!/usr/bin/env python3
"""
PropertyProcessor Demonstration Script

This script demonstrates the comprehensive capabilities of the PropertyProcessor:

1. Basic property assignment to entities
2. Cross-format entity property management 
3. Property dataset creation and management
4. Advanced secondary selection (GRN positions, atom selections)
5. Data import/export functionality
6. Property filtering and analysis
7. Integration with entity registry system

The demo shows how properties can be assigned to any entity type (structure, 
sequence, GRN, embedding) and organized into meaningful datasets for analysis.
"""

import os
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import tempfile
import json

# Add protos to path
protos_path = Path(__file__).parent.parent
sys.path.insert(0, str(protos_path))

# Import Protos components
try:
    from protos.io.paths.path_config import ProtosPaths
    from protos.processing.property.property_processor_enhanced import PropertyProcessor
    from protos.io.data_access import generate_entity_id
    print("✅ All imports successful!")
except ImportError as e:
    print(f"❌ Import error: {e}")
    print("Make sure protos is installed: pip install -e .")
    sys.exit(1)


def setup_demo_environment():
    """Set up the demo environment with data directory."""
    print("=== Setting Up Demo Environment ===")
    
    # Set up data directory
    demo_data_dir = protos_path / "data"
    ProtosPaths.set_data_root(str(demo_data_dir.absolute()))
    
    print(f"Demo data directory: {demo_data_dir}")
    return demo_data_dir


def demo_basic_property_assignment():
    """Demonstrate basic property assignment functionality."""
    print("\n=== Basic Property Assignment Demo ===")
    
    # Initialize PropertyProcessor
    prop_proc = PropertyProcessor(name="demo_properties")
    print(f"✅ PropertyProcessor initialized at {prop_proc.data_path}")
    
    # Sample entities from different formats
    entities = {
        "36c2c0da93": {"original_id": "1ubq", "type": "structure"},
        "7e77394211": {"original_id": "TEST_PROTEIN", "type": "sequence"},
        "b3c4d5e6f7": {"original_id": "GRN_ENTRY_1", "type": "grn"},
        "a1b2c3d4e5": {"original_id": "OPSIN_001", "type": "sequence"}
    }
    
    print("\n📊 Assigning properties to entities...")
    
    # Assign properties to microbial opsins dataset
    microbial_opsin_properties = [
        ("36c2c0da93", "lambda_max", 568, {"method": "spectroscopy", "pH": 7.0}),
        ("7e77394211", "lambda_max", 500, {"method": "spectroscopy", "pH": 7.0}),
        ("a1b2c3d4e5", "lambda_max", 485, {"method": "spectroscopy", "pH": 7.0}),
        ("36c2c0da93", "organism", "Halobacterium salinarum", {"source": "NCBI"}),
        ("7e77394211", "organism", "Gloeobacter violaceus", {"source": "NCBI"}),
        ("a1b2c3d4e5", "organism", "Anabaena sensilis", {"source": "NCBI"}),
    ]
    
    for entity_id, prop_name, prop_value, metadata in microbial_opsin_properties:
        prop_proc.assign_property(
            entity_identifier=entity_id,
            property_name=prop_name,
            property_value=prop_value,
            dataset_name="microbial_opsins",
            metadata=metadata
        )
        print(f"   ✅ {entities[entity_id]['original_id']}: {prop_name} = {prop_value}")
    
    # Assign structural properties
    structural_properties = [
        ("36c2c0da93", "resolution", 1.55, {"method": "X-ray crystallography"}),
        ("36c2c0da93", "space_group", "P212121", {"method": "X-ray crystallography"}),
        ("36c2c0da93", "pdb_id", "1C3W", {"database": "PDB"}),
    ]
    
    for entity_id, prop_name, prop_value, metadata in structural_properties:
        prop_proc.assign_property(
            entity_identifier=entity_id,
            property_name=prop_name,
            property_value=prop_value,
            dataset_name="structural_data",
            metadata=metadata
        )
        print(f"   ✅ Structure properties: {prop_name} = {prop_value}")
    
    return prop_proc


def demo_property_retrieval(prop_proc):
    """Demonstrate property retrieval methods."""
    print("\n=== Property Retrieval Demo ===")
    
    # Get specific property
    lambda_max = prop_proc.get_entity_property("36c2c0da93", "lambda_max", "microbial_opsins")
    print(f"🔍 Lambda max for 1ubq: {lambda_max} nm")
    
    # Get all properties for an entity from specific dataset
    microbial_props = prop_proc.get_entity_properties("36c2c0da93", "microbial_opsins")
    print(f"🔍 Microbial opsin properties for 1ubq: {microbial_props}")
    
    # Get all properties for an entity across all datasets
    all_props = prop_proc.get_entity_properties("36c2c0da93")
    print(f"🔍 All properties for 1ubq: {all_props}")
    
    # Search for property across all datasets
    organism = prop_proc.get_entity_property("7e77394211", "organism")
    print(f"🔍 Organism for TEST_PROTEIN: {organism}")


def demo_dataset_management(prop_proc):
    """Demonstrate dataset management functionality."""
    print("\n=== Dataset Management Demo ===")
    
    # List available datasets
    datasets = prop_proc.list_datasets()
    print(f"📚 Available datasets ({len(datasets)}):")
    for dataset in datasets:
        print(f"   - {dataset['id']}: {dataset['entity_count']} entities, "
              f"{dataset['property_count']} properties")
    
    # Get dataset statistics
    for dataset_name in ["microbial_opsins", "structural_data"]:
        stats = prop_proc.get_dataset_statistics(dataset_name)
        print(f"\n📊 Statistics for '{dataset_name}':")
        print(f"   Entities: {stats['entity_count']}")
        print(f"   Properties: {stats['property_count']}")
        print(f"   Created: {stats['created_at']}")
        
        if 'properties' in stats:
            print("   Property details:")
            for prop_name, prop_stats in stats['properties'].items():
                print(f"     - {prop_name}: {prop_stats['type']}, "
                      f"non-null: {prop_stats['non_null_count']}")
                if prop_stats['type'] in ['int64', 'float64']:
                    print(f"       Mean: {prop_stats.get('mean', 'N/A'):.2f}, "
                          f"Range: [{prop_stats.get('min', 'N/A'):.2f}, "
                          f"{prop_stats.get('max', 'N/A'):.2f}]")
    
    # Get dataset as DataFrame
    print(f"\n📋 Microbial opsins dataset as DataFrame:")
    df = prop_proc.get_dataset_properties("microbial_opsins")
    print(df)
    
    return df


def demo_property_filtering(prop_proc):
    """Demonstrate property filtering functionality."""
    print("\n=== Property Filtering Demo ===")
    
    # Filter entities by lambda_max value
    high_lambda = prop_proc.filter_entities_by_property(
        "microbial_opsins", {"lambda_max": {"gt": 520}}
    )
    print(f"🔍 Entities with lambda_max > 520 nm: {high_lambda}")
    
    # Filter by specific organism
    halobacterium = prop_proc.filter_entities_by_property(
        "microbial_opsins", {"organism": "Halobacterium salinarum"}
    )
    print(f"🔍 Halobacterium entities: {halobacterium}")
    
    # Combined filters
    specific_opsins = prop_proc.filter_entities_by_property(
        "microbial_opsins", {
            "lambda_max": {"gt": 500},
            "organism": {"in": ["Halobacterium salinarum", "Gloeobacter violaceus"]}
        }
    )
    print(f"🔍 Specific filtered opsins: {specific_opsins}")


def demo_batch_operations(prop_proc):
    """Demonstrate batch property assignment."""
    print("\n=== Batch Operations Demo ===")
    
    # Prepare batch assignments for a new dataset
    batch_assignments = [
        {
            'entity_identifier': 'c1d2e3f4g5',
            'property_name': 'expression_level',
            'property_value': 'high',
            'metadata': {'system': 'E. coli', 'temperature': 37}
        },
        {
            'entity_identifier': 'c1d2e3f4g5',
            'property_name': 'stability_score',
            'property_value': 8.5,
            'metadata': {'assay': 'thermal_stability', 'units': 'score'}
        },
        {
            'entity_identifier': 'd2e3f4g5h6',
            'property_name': 'expression_level',
            'property_value': 'medium',
            'metadata': {'system': 'E. coli', 'temperature': 37}
        },
        {
            'entity_identifier': 'd2e3f4g5h6',
            'property_name': 'stability_score',
            'property_value': 6.2,
            'metadata': {'assay': 'thermal_stability', 'units': 'score'}
        }
    ]
    
    # Batch assign properties
    entity_ids = prop_proc.assign_properties_batch(
        batch_assignments, "expression_data"
    )
    
    print(f"✅ Batch assigned {len(entity_ids)} properties to {len(set(entity_ids))} unique entities")
    
    # Show results
    df = prop_proc.get_dataset_properties("expression_data")
    print("📋 Expression data dataset:")
    print(df)


def demo_secondary_selection(prop_proc):
    """Demonstrate advanced secondary selection features."""
    print("\n=== Secondary Selection Demo ===")
    
    # GRN position-specific properties
    print("🧬 Assigning GRN position-specific properties...")
    
    grn_entity = "b3c4d5e6f7"
    grn_properties = [
        ("3.50", "amino_acid", "Asp", {"conservation": "highly_conserved"}),
        ("3.50", "polarity", "polar", {"analysis": "hydrophobicity_scale"}),
        ("7.45", "amino_acid", "Lys", {"conservation": "moderately_conserved"}),
        ("7.45", "polarity", "polar", {"analysis": "hydrophobicity_scale"}),
        ("6.50", "amino_acid", "Trp", {"conservation": "highly_conserved"}),
        ("6.50", "polarity", "aromatic", {"analysis": "hydrophobicity_scale"}),
    ]
    
    for grn_pos, prop_name, prop_value, metadata in grn_properties:
        prop_proc.assign_grn_property(
            grn_entity, grn_pos, prop_name, prop_value, "grn_analysis"
        )
        print(f"   ✅ GRN {grn_pos}: {prop_name} = {prop_value}")
    
    # Atom-specific properties for structures
    print("\n🏗️ Assigning atom-specific properties...")
    
    structure_entity = "36c2c0da93"
    atom_properties = [
        ({"chain": "A", "residue": 100}, "b_factor", 25.5),
        ({"chain": "A", "residue": 100}, "occupancy", 1.0),
        ({"chain": "A", "atom_name": "CA"}, "coordination", "tetrahedral"),
        ({"chain": "B", "residue": 50}, "b_factor", 18.2),
    ]
    
    for atom_selector, prop_name, prop_value in atom_properties:
        prop_proc.assign_atom_property(
            structure_entity, atom_selector, prop_name, prop_value, "atomic_properties"
        )
        selector_str = "_".join([f"{k}{v}" for k, v in atom_selector.items()])
        print(f"   ✅ Atom {selector_str}: {prop_name} = {prop_value}")
    
    # Show secondary selection results
    grn_props = prop_proc.get_entity_properties(grn_entity, "grn_analysis")
    print(f"\n🔍 GRN-specific properties: {grn_props}")
    
    atom_props = prop_proc.get_entity_properties(structure_entity, "atomic_properties")
    print(f"🔍 Atom-specific properties: {atom_props}")


def demo_file_operations(prop_proc):
    """Demonstrate file import/export operations."""
    print("\n=== File Operations Demo ===")
    
    # Create sample CSV data for import
    sample_data = pd.DataFrame({
        'entity_id': ['e1f2g3h4i5', 'f2g3h4i5j6', 'g3h4i5j6k7'],
        'protein_family': ['rhodopsin', 'opsin', 'rhodopsin'],
        'membrane_spans': [7, 7, 7],
        'molecular_weight': [39000, 38500, 39200],
        'isoelectric_point': [5.8, 6.1, 5.9]
    })
    
    # Save to temporary CSV file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        sample_data.to_csv(f.name, index=False)
        csv_file = f.name
    
    print(f"📄 Created sample CSV file: {csv_file}")
    print("Sample data:")
    print(sample_data)
    
    # Import dataset from CSV
    count = prop_proc.create_property_dataset_from_file(
        csv_file, "imported_proteins", entity_column='entity_id'
    )
    
    print(f"✅ Imported {count} entities from CSV file")
    
    # Show imported dataset
    imported_df = prop_proc.get_dataset_properties("imported_proteins")
    print("📋 Imported dataset:")
    print(imported_df)
    
    # Save datasets to files
    print("\n💾 Saving datasets to files...")
    
    for dataset_name in ["microbial_opsins", "imported_proteins"]:
        prop_proc.save_property_dataset(dataset_name, 'both')  # Save as both JSON and CSV
        print(f"   ✅ Saved dataset '{dataset_name}' as JSON and CSV")
    
    # Clean up temporary file
    os.unlink(csv_file)
    
    return imported_df


def demo_cross_format_integration(prop_proc):
    """Demonstrate cross-format entity integration."""
    print("\n=== Cross-Format Integration Demo ===")
    
    # Simulate an entity that exists in multiple formats
    multi_format_entity = "36c2c0da93"  # This could be a structure that also has sequence, GRN, etc.
    
    print(f"🔗 Multi-format entity: {multi_format_entity}")
    
    # Assign properties from different "perspectives" or analyses
    cross_format_properties = [
        ("structural_analysis", "secondary_structure", "7TM_helical", {"method": "DSSP"}),
        ("sequence_analysis", "hydrophobicity", 0.65, {"method": "Kyte_Doolittle"}),
        ("functional_analysis", "ligand_binding", "retinal", {"assay": "spectroscopy"}),
        ("evolutionary_analysis", "conservation_score", 0.89, {"database": "ConSurf"}),
    ]
    
    for dataset, prop_name, prop_value, metadata in cross_format_properties:
        prop_proc.assign_property(
            multi_format_entity, prop_name, prop_value, dataset, metadata=metadata
        )
        print(f"   ✅ {dataset}: {prop_name} = {prop_value}")
    
    # Show all properties for this multi-format entity
    all_properties = prop_proc.get_entity_properties(multi_format_entity)
    print(f"\n🔍 All properties for multi-format entity {multi_format_entity}:")
    for prop_name, prop_value in all_properties.items():
        print(f"   {prop_name}: {prop_value}")
    
    # Show entity across different datasets
    print(f"\n📊 Entity {multi_format_entity} appears in datasets:")
    for dataset_name in prop_proc.property_datasets.keys():
        entities = prop_proc.list_entities(dataset_name)
        if multi_format_entity in entities:
            props = prop_proc.get_entity_properties(multi_format_entity, dataset_name)
            print(f"   - {dataset_name}: {list(props.keys())}")


def demo_advanced_analysis(prop_proc):
    """Demonstrate advanced analysis capabilities."""
    print("\n=== Advanced Analysis Demo ===")
    
    # Comprehensive analysis across all datasets
    all_datasets = prop_proc.list_datasets()
    print(f"📊 Analysis across {len(all_datasets)} datasets:")
    
    total_entities = len(prop_proc.list_entities())
    total_properties = sum(ds['property_count'] for ds in all_datasets)
    
    print(f"   Total unique entities: {total_entities}")
    print(f"   Total property assignments: {total_properties}")
    
    # Find entities with most properties
    entity_property_counts = {}
    for entity_id in prop_proc.list_entities():
        props = prop_proc.get_entity_properties(entity_id)
        entity_property_counts[entity_id] = len(props)
    
    # Sort by property count
    sorted_entities = sorted(entity_property_counts.items(), 
                           key=lambda x: x[1], reverse=True)
    
    print(f"\n🏆 Top entities by property count:")
    for i, (entity_id, count) in enumerate(sorted_entities[:5]):
        print(f"   {i+1}. {entity_id}: {count} properties")
    
    # Property usage statistics
    property_usage = {}
    for dataset_name in prop_proc.property_datasets:
        df = prop_proc.get_dataset_properties(dataset_name)
        for col in df.columns:
            if col not in property_usage:
                property_usage[col] = 0
            property_usage[col] += df[col].notna().sum()
    
    print(f"\n📈 Most commonly used properties:")
    sorted_props = sorted(property_usage.items(), key=lambda x: x[1], reverse=True)
    for prop_name, usage_count in sorted_props[:10]:
        print(f"   {prop_name}: {usage_count} assignments")


def main():
    """Run the complete PropertyProcessor demonstration."""
    print("=" * 80)
    print("PROPERTY PROCESSOR COMPREHENSIVE DEMONSTRATION")
    print("=" * 80)
    
    # Set up environment
    demo_data_dir = setup_demo_environment()
    
    # Run demonstrations
    prop_proc = demo_basic_property_assignment()
    demo_property_retrieval(prop_proc)
    dataset_df = demo_dataset_management(prop_proc)
    demo_property_filtering(prop_proc)
    demo_batch_operations(prop_proc)
    demo_secondary_selection(prop_proc)
    imported_df = demo_file_operations(prop_proc)
    demo_cross_format_integration(prop_proc)
    demo_advanced_analysis(prop_proc)
    
    # Final summary
    print("\n" + "=" * 80)
    print("DEMONSTRATION SUMMARY")
    print("=" * 80)
    
    print(f"\n✅ PropertyProcessor successfully demonstrated:")
    print(f"   📊 {len(prop_proc.list_datasets())} property datasets created")
    print(f"   🔗 {len(prop_proc.list_entities())} entities with properties")
    print(f"   🎯 Advanced features: GRN positions, atom selections")
    print(f"   💾 File I/O: CSV/JSON import/export")
    print(f"   🔍 Filtering and analysis capabilities")
    print(f"   🌐 Cross-format entity integration")
    
    print(f"\n🎯 Key Features Demonstrated:")
    features = [
        "✅ Entity-based property assignment across all format types",
        "✅ Dataset organization and management",
        "✅ Property filtering and querying",
        "✅ Batch operations for efficiency",
        "✅ Secondary selection (GRN positions, atom selections)",
        "✅ File import/export (CSV, JSON formats)",
        "✅ Cross-format entity integration",
        "✅ Comprehensive statistics and analysis",
        "✅ Metadata tracking and validation",
        "✅ Integration with Protos entity registry"
    ]
    
    for feature in features:
        print(f"   {feature}")
    
    print(f"\n📁 Data saved to: {prop_proc.data_path}")
    print(f"🚀 PropertyProcessor is ready for production use!")


if __name__ == "__main__":
    main()