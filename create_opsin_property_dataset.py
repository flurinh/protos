#!/usr/bin/env python3
"""
Create realistic opsin property dataset linked to reference sequences

This script:
1. Extracts sequence names from the reference opsin FASTA file
2. Creates realistic property data for these opsins based on known biology
3. Tests PropertyProcessor integration with real reference data
4. Demonstrates the complete workflow with actual Protos sequences
"""

import os
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import random
import re

# Add protos to path
protos_path = Path(__file__).parent
sys.path.insert(0, str(protos_path))

# Import Protos components
try:
    from protos.io.paths.path_config import ProtosPaths
    from protos.processing.property.property_processor_enhanced import PropertyProcessor
    from protos.processing.sequence.seq_processor import SeqProcessor
    from protos.io.data_access import generate_entity_id
    from protos.io.fasta_utils import read_fasta
    print("✅ All imports successful!")
except ImportError as e:
    print(f"❌ Import error: {e}")
    print("Make sure protos is installed: pip install -e .")
    sys.exit(1)


def extract_sequence_names_from_fasta():
    """Extract sequence names from the reference opsin FASTA file."""
    print("=== Extracting Sequence Names from Reference FASTA ===")
    
    # Path to reference opsin sequences
    fasta_path = protos_path / "src" / "protos" / "reference_data" / "sequence" / "fasta" / "opsin_sequences_from_yaml.fasta"
    
    if not fasta_path.exists():
        print(f"❌ Reference FASTA file not found: {fasta_path}")
        return []
    
    print(f"📁 Reading FASTA file: {fasta_path}")
    
    # Read FASTA sequences
    try:
        sequences = read_fasta(str(fasta_path))
        sequence_names = list(sequences.keys())
        
        print(f"✅ Found {len(sequence_names)} opsin sequences")
        print("First 10 sequences:")
        for i, name in enumerate(sequence_names[:10]):
            seq_length = len(sequences[name])
            print(f"   {i+1:2d}. {name:15} (length: {seq_length} AA)")
        
        return sequence_names, sequences
        
    except Exception as e:
        print(f"❌ Error reading FASTA file: {e}")
        return [], {}


def categorize_opsins(sequence_names):
    """Categorize opsins based on their names into functional groups."""
    print("\n=== Categorizing Opsins by Function ===")
    
    # Define opsin categories based on naming patterns
    categories = {
        'Archaerhodopsin': {
            'pattern': r'^AR\d*_',
            'description': 'Archaerhodopsin variants',
            'lambda_range': (560, 580),
            'organisms': ['Halobacterium', 'Haloquadratum', 'Halorubrum']
        },
        'Bacteriorhodopsin': {
            'pattern': r'^(BR|bR|ASR)',
            'description': 'Bacteriorhodopsin and related proton pumps',
            'lambda_range': (568, 578),
            'organisms': ['Halobacterium salinarum', 'Haloquadratum walsbyi']
        },
        'Channelrhodopsin': {
            'pattern': r'ChR|ChRmine|CnChR|VcChR|CoChR',
            'description': 'Light-gated ion channels',
            'lambda_range': (470, 530),
            'organisms': ['Chlamydomonas', 'Volvox', 'Chlorogonium']
        },
        'Halorhodopsin': {
            'pattern': r'HR|HeR|HrHR|HaHR',
            'description': 'Chloride pumps',
            'lambda_range': (570, 590),
            'organisms': ['Halobacterium', 'Haloarcula', 'Natronobacterium']
        },
        'Sensory_rhodopsin': {
            'pattern': r'SR|SRI|SRII',
            'description': 'Photosensory receptors',
            'lambda_range': (480, 520),
            'organisms': ['Halobacterium', 'Natronobacterium', 'Anabaena']
        },
        'Proteorhodopsin': {
            'pattern': r'PR|pR',
            'description': 'Marine proteorhodopsins',
            'lambda_range': (490, 530),
            'organisms': ['Marine bacteria', 'SAR11', 'Pelagibacter']
        },
        'Actinorhodopsin': {
            'pattern': r'ACR|ActR',
            'description': 'Anion channelrhodopsins',
            'lambda_range': (515, 545),
            'organisms': ['Mesostigma', 'Guillardia', 'Cryptophytes']
        },
        'Xanthorhodopsin': {
            'pattern': r'XR|XeR',
            'description': 'Light-driven proton pumps with antenna',
            'lambda_range': (540, 570),
            'organisms': ['Salinibacter', 'Gloeobacter']
        },
        'Other': {
            'pattern': r'.*',
            'description': 'Other rhodopsin variants',
            'lambda_range': (450, 600),
            'organisms': ['Various']
        }
    }
    
    categorized = {}
    category_counts = {}
    
    for name in sequence_names:
        assigned = False
        for category, info in categories.items():
            if category == 'Other':
                continue
            if re.search(info['pattern'], name, re.IGNORECASE):
                if category not in categorized:
                    categorized[category] = []
                categorized[category].append(name)
                category_counts[category] = category_counts.get(category, 0) + 1
                assigned = True
                break
        
        if not assigned:
            if 'Other' not in categorized:
                categorized['Other'] = []
            categorized['Other'].append(name)
            category_counts['Other'] = category_counts.get('Other', 0) + 1
    
    print("Opsin categorization results:")
    for category, count in category_counts.items():
        description = categories[category]['description']
        print(f"   {category:20} {count:3d} sequences - {description}")
    
    return categorized, categories


def generate_realistic_properties(sequence_names, sequences, categorized, categories):
    """Generate realistic property data for the opsin sequences."""
    print("\n=== Generating Realistic Property Data ===")
    
    properties_data = []
    
    # Set random seed for reproducible data
    random.seed(42)
    np.random.seed(42)
    
    for seq_name in sequence_names:
        # Determine category
        category = 'Other'
        for cat, seqs in categorized.items():
            if seq_name in seqs:
                category = cat
                break
        
        cat_info = categories[category]
        seq_length = len(sequences[seq_name])
        
        # Generate lambda_max based on category
        lambda_min, lambda_max = cat_info['lambda_range']
        lambda_max_val = np.random.uniform(lambda_min, lambda_max)
        
        # Generate other properties based on realistic biological ranges
        properties = {
            'protein_id': seq_name,
            'category': category,
            'lambda_max': round(lambda_max_val, 1),
            'sequence_length': seq_length,
            'organism': np.random.choice(cat_info['organisms']),
            'expression_level': np.random.choice(['low', 'medium', 'high'], p=[0.3, 0.5, 0.2]),
            'thermal_stability': round(np.random.normal(55, 15), 1),  # Tm in °C
            'ph_optimum': round(np.random.uniform(6.0, 8.5), 1),
            'extinction_coefficient': int(np.random.uniform(20000, 80000)),  # M⁻¹cm⁻¹
            'quantum_yield': round(np.random.uniform(0.1, 0.8), 3),
            'photocycle_time': round(np.random.uniform(0.5, 50), 1),  # ms
            'membrane_spanning': 7,  # All are 7TM proteins
            'molecular_weight': round(seq_length * 110 / 1000, 1),  # Approximate MW in kDa
            'isoelectric_point': round(np.random.uniform(4.5, 9.5), 1),
            'hydrophobicity_score': round(np.random.uniform(0.3, 0.8), 2),
            'conservation_score': round(np.random.uniform(0.6, 0.95), 2),
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
        elif category in ['Bacteriorhodopsin', 'Proteorhodopsin', 'Xanthorhodopsin']:
            properties['ion_selectivity'] = 'proton'
            properties['pump_rate'] = round(np.random.uniform(50, 500), 1)  # H+/s
        elif category == 'Sensory_rhodopsin':
            properties['signaling_type'] = np.random.choice(['phototaxis', 'photophobic'])
            properties['light_adaptation'] = np.random.choice(['light', 'dark', 'both'])
        
        properties_data.append(properties)
    
    # Convert to DataFrame
    df = pd.DataFrame(properties_data)
    
    print(f"✅ Generated properties for {len(df)} opsin sequences")
    print(f"Property columns: {list(df.columns)}")
    
    # Show sample data
    print("\nSample property data:")
    print(df.head(3).to_string(index=False))
    
    return df


def save_property_dataset(properties_df):
    """Save the property dataset to CSV file."""
    print("\n=== Saving Property Dataset ===")
    
    # Save to data directory
    data_dir = protos_path / "data"
    output_file = data_dir / "opsin_properties_reference.csv"
    
    # Ensure directory exists
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    # Save CSV
    properties_df.to_csv(output_file, index=False)
    
    print(f"💾 Saved property dataset to: {output_file}")
    print(f"   Rows: {len(properties_df)}")
    print(f"   Columns: {len(properties_df.columns)}")
    
    return output_file


def test_property_integration(properties_df, csv_file):
    """Test PropertyProcessor integration with the reference sequences."""
    print("\n=== Testing PropertyProcessor Integration ===")
    
    # Set up data directory
    data_dir = protos_path / "data"
    ProtosPaths.set_data_root(str(data_dir.absolute()))
    
    # Initialize PropertyProcessor
    prop_proc = PropertyProcessor(name="opsin_reference_test")
    print("✅ PropertyProcessor initialized")
    
    # Import properties from CSV
    try:
        count = prop_proc.create_property_dataset_from_file(
            str(csv_file),
            "opsin_reference_properties",
            entity_column='protein_id'
        )
        
        print(f"✅ Successfully imported properties for {count} opsin sequences")
        
        # Get dataset statistics
        stats = prop_proc.get_dataset_statistics("opsin_reference_properties")
        print(f"📊 Dataset statistics:")
        print(f"   Entities: {stats['entity_count']}")
        print(f"   Properties: {stats['property_count']}")
        
        # Test property retrieval
        print(f"\n🔍 Testing property retrieval:")
        sample_sequences = properties_df['protein_id'].head(5).tolist()
        
        for seq_name in sample_sequences:
            props = prop_proc.get_entity_properties(seq_name, "opsin_reference_properties")
            entity_id = generate_entity_id(seq_name)
            
            print(f"   {seq_name:15} (ID: {entity_id[:10]}): {len(props)} properties")
            
            # Show key properties
            key_props = ['lambda_max', 'category', 'organism', 'expression_level']
            for prop in key_props:
                if prop in props:
                    print(f"     {prop}: {props[prop]}")
        
        # Test filtering by opsin category
        print(f"\n🔍 Testing property filtering:")
        
        # Filter by category
        channelrhodopsins = prop_proc.filter_entities_by_property(
            "opsin_reference_properties", {"category": "Channelrhodopsin"}
        )
        print(f"   Channelrhodopsins: {len(channelrhodopsins)} found")
        
        # Filter by lambda_max range
        blue_opsins = prop_proc.filter_entities_by_property(
            "opsin_reference_properties", {"lambda_max": {"lt": 500}}
        )
        print(f"   Blue-shifted opsins (λ < 500nm): {len(blue_opsins)} found")
        
        # Filter by high expression
        high_expr = prop_proc.filter_entities_by_property(
            "opsin_reference_properties", {"expression_level": "high"}
        )
        print(f"   High expression opsins: {len(high_expr)} found")
        
        # Filter by thermal stability
        thermostable = prop_proc.filter_entities_by_property(
            "opsin_reference_properties", {"thermal_stability": {"gt": 60}}
        )
        print(f"   Thermostable opsins (Tm > 60°C): {len(thermostable)} found")
        
        return prop_proc
        
    except Exception as e:
        print(f"❌ PropertyProcessor integration failed: {e}")
        import traceback
        traceback.print_exc()
        return None


def test_sequence_processor_integration():
    """Test integration with SeqProcessor to verify sequences are available."""
    print("\n=== Testing SeqProcessor Integration ===")
    
    try:
        # Initialize SeqProcessor
        seq_proc = SeqProcessor(name="opsin_reference_test")
        print("✅ SeqProcessor initialized")
        
        # Check available sequences
        sequences = seq_proc.list_entities()
        print(f"📊 Available sequences in SeqProcessor: {len(sequences)}")
        
        # Check if our reference sequences are available
        reference_file = "opsin_sequences_from_yaml.fasta"
        
        # Look for sequences that might match our reference data
        opsin_sequences = [seq for seq in sequences if 'opsin' in seq.lower() or any(
            pattern in seq for pattern in ['AR', 'BR', 'ChR', 'HR', 'SR', 'PR']
        )]
        
        print(f"📊 Potential opsin sequences found: {len(opsin_sequences)}")
        if opsin_sequences:
            print("Sample opsin sequences:")
            for seq in opsin_sequences[:5]:
                print(f"   - {seq}")
        
        return len(opsin_sequences) > 0
        
    except Exception as e:
        print(f"❌ SeqProcessor integration test failed: {e}")
        return False


def analyze_property_statistics(properties_df):
    """Analyze the generated property statistics."""
    print("\n=== Property Dataset Analysis ===")
    
    # Category distribution
    print("📊 Opsin category distribution:")
    category_counts = properties_df['category'].value_counts()
    for category, count in category_counts.items():
        percentage = (count / len(properties_df)) * 100
        print(f"   {category:20} {count:3d} ({percentage:5.1f}%)")
    
    # Lambda max distribution
    print(f"\n📊 Lambda max statistics:")
    lambda_stats = properties_df['lambda_max'].describe()
    print(f"   Mean: {lambda_stats['mean']:.1f} nm")
    print(f"   Range: {lambda_stats['min']:.1f} - {lambda_stats['max']:.1f} nm")
    print(f"   Std: {lambda_stats['std']:.1f} nm")
    
    # Expression level distribution
    print(f"\n📊 Expression level distribution:")
    expr_counts = properties_df['expression_level'].value_counts()
    for level, count in expr_counts.items():
        percentage = (count / len(properties_df)) * 100
        print(f"   {level:10} {count:3d} ({percentage:5.1f}%)")
    
    # Organism distribution
    print(f"\n📊 Top organisms:")
    org_counts = properties_df['organism'].value_counts().head(10)
    for org, count in org_counts.items():
        print(f"   {org:25} {count:3d}")
    
    # Functional activity
    print(f"\n📊 Functional activity:")
    activity_counts = properties_df['functional_activity'].value_counts()
    for activity, count in activity_counts.items():
        percentage = (count / len(properties_df)) * 100
        print(f"   {activity:10} {count:3d} ({percentage:5.1f}%)")


def main():
    """Run the complete opsin property dataset creation and testing workflow."""
    print("=" * 80)
    print("OPSIN REFERENCE PROPERTY DATASET CREATION")
    print("=" * 80)
    
    # Extract sequence names from reference FASTA
    sequence_names, sequences = extract_sequence_names_from_fasta()
    
    if not sequence_names:
        print("❌ No sequences found, exiting")
        return
    
    # Categorize opsins by function
    categorized, categories = categorize_opsins(sequence_names)
    
    # Generate realistic properties
    properties_df = generate_realistic_properties(sequence_names, sequences, categorized, categories)
    
    # Analyze property statistics
    analyze_property_statistics(properties_df)
    
    # Save property dataset
    csv_file = save_property_dataset(properties_df)
    
    # Test PropertyProcessor integration
    prop_proc = test_property_integration(properties_df, csv_file)
    
    # Test SeqProcessor integration
    seq_integration = test_sequence_processor_integration()
    
    # Final summary
    print("\n" + "=" * 80)
    print("OPSIN REFERENCE DATASET SUMMARY")
    print("=" * 80)
    
    print(f"\n✅ Successfully created opsin property dataset:")
    print(f"   📊 {len(sequence_names)} opsin sequences from reference FASTA")
    print(f"   📊 {len(categorized)} functional categories identified")
    print(f"   📊 {len(properties_df.columns)-1} property columns generated")
    print(f"   💾 Dataset saved to: {csv_file}")
    
    print(f"\n🔗 Integration test results:")
    print(f"   PropertyProcessor: {'✅ SUCCESS' if prop_proc else '❌ FAILED'}")
    print(f"   SeqProcessor: {'✅ SUCCESS' if seq_integration else '❌ FAILED'}")
    
    print(f"\n🎯 Property categories include:")
    categories_list = [
        "✅ Spectroscopic properties (lambda_max, extinction_coefficient)",
        "✅ Biophysical properties (thermal_stability, pH_optimum)",
        "✅ Functional properties (expression_level, activity)",
        "✅ Structural properties (sequence_length, molecular_weight)",
        "✅ Photochemical properties (quantum_yield, photocycle_time)",
        "✅ Classification data (category, organism, ion_selectivity)"
    ]
    
    for category in categories_list:
        print(f"   {category}")
    
    print(f"\n🚀 Opsin reference property dataset ready for use!")
    
    return properties_df, prop_proc


if __name__ == "__main__":
    main()