# PropertyProcessor

The PropertyProcessor manages experimental properties, metadata, and annotations associated with biological entities. It provides a flexible system for storing, querying, and analyzing diverse property types linked to proteins, structures, and other biological data.

## Overview

PropertyProcessor handles:
- Experimental property storage (Kd, IC50, melting temperature, etc.)
- Metadata management for entities
- Property-based filtering and queries
- Integration with other data formats
- Statistical analysis of properties
- Import/export from various sources

## Basic Usage

### Initialization

```python
from protos.processing.property.property_processor import PropertyProcessor

# Create processor instance
pp = PropertyProcessor(name="properties")

# With options
pp = PropertyProcessor(
    name="experimental_data",
    validate_on_load=True,  # Validate property consistency
    enforce_types=True      # Enforce property data types
)
```

### Creating Property Datasets

```python
# Create from DataFrame
import pandas as pd

properties_df = pd.DataFrame({
    'protein_id': ['EGFR_HUMAN', 'BRAF_HUMAN', 'CDK2_HUMAN'],
    'Kd_nM': [0.5, 1.2, 3.4],
    'IC50_nM': [2.1, 5.6, 8.9],
    'expression_system': ['HEK293', 'E.coli', 'Insect'],
    'active': [True, True, False]
})

# Create property dataset
pp.create_property_dataset(
    dataset_id="kinase_inhibitors",
    properties_df=properties_df,
    entity_column='protein_id',
    description="Kinase inhibitor screening results"
)
```

### Loading Property Data

```python
# Load from CSV
pp.create_property_dataset_from_file(
    csv_file="experimental_data.csv",
    dataset_id="screening_results",
    entity_column='protein_name'
)

# Load from Excel
pp.create_property_dataset_from_file(
    excel_file="assay_results.xlsx",
    dataset_id="assay_data",
    entity_column='sample_id',
    sheet_name="Results"
)

# Load existing dataset
pp.load_property_dataset("kinase_inhibitors")
```

## Property Operations

### Querying Properties

```python
# Get properties for specific entity
properties = pp.get_entity_properties("EGFR_HUMAN")
# {'Kd_nM': 0.5, 'IC50_nM': 2.1, 'expression_system': 'HEK293', 'active': True}

# Get specific property across all entities
kd_values = pp.get_property_values("Kd_nM", dataset_id="kinase_inhibitors")
# {'EGFR_HUMAN': 0.5, 'BRAF_HUMAN': 1.2, 'CDK2_HUMAN': 3.4}

# Query multiple properties
multi_props = pp.get_properties(
    entities=["EGFR_HUMAN", "BRAF_HUMAN"],
    properties=["Kd_nM", "IC50_nM"]
)
```

### Filtering by Properties

```python
# Filter by single property
high_affinity = pp.filter_entities_by_property(
    dataset_id="kinase_inhibitors",
    filters={"Kd_nM": {"lt": 1.0}}  # Kd < 1.0 nM
)

# Complex filtering
active_potent = pp.filter_entities_by_property(
    dataset_id="kinase_inhibitors",
    filters={
        "active": True,
        "IC50_nM": {"lt": 10.0},
        "expression_system": {"in": ["HEK293", "CHO"]}
    }
)

# Range filtering
mid_range = pp.filter_entities_by_property(
    dataset_id="screening_results",
    filters={
        "pIC50": {"gte": 6.0, "lte": 8.0}  # 6 <= pIC50 <= 8
    }
)
```

### Property Assignment

```python
# Assign single property
pp.assign_property(
    entity_name="NEW_PROTEIN",
    property_name="melting_temp",
    value=65.5,
    dataset_id="thermal_stability"
)

# Batch assignment
new_properties = {
    "PROT1": {"Tm": 55.2, "pH_optimum": 7.4},
    "PROT2": {"Tm": 62.1, "pH_optimum": 8.0},
    "PROT3": {"Tm": 48.9, "pH_optimum": 6.5}
}

pp.batch_assign_properties(
    properties=new_properties,
    dataset_id="protein_stability"
)

# Update existing properties
pp.update_property(
    entity_name="EGFR_HUMAN",
    property_name="IC50_nM",
    new_value=1.8,
    dataset_id="kinase_inhibitors"
)
```

## Property Types and Validation

### Supported Property Types

```python
# Numeric properties
numeric_props = {
    "Kd": float,          # Dissociation constant
    "IC50": float,        # Half-maximal inhibitory concentration
    "EC50": float,        # Half-maximal effective concentration
    "Tm": float,          # Melting temperature
    "MW": float,          # Molecular weight
    "pI": float           # Isoelectric point
}

# Categorical properties
categorical_props = {
    "expression_system": ["E.coli", "HEK293", "CHO", "Insect"],
    "assay_type": ["FRET", "FP", "SPR", "ITC"],
    "activity_class": ["active", "inactive", "intermediate"]
}

# Boolean properties
boolean_props = {
    "is_active": bool,
    "passes_qc": bool,
    "is_control": bool
}

# Text properties
text_props = {
    "notes": str,
    "batch_id": str,
    "experimenter": str
}
```

### Property Validation

```python
# Define property schema
property_schema = {
    "Kd_nM": {
        "type": float,
        "min": 0.001,
        "max": 10000,
        "required": True
    },
    "expression_system": {
        "type": str,
        "allowed_values": ["E.coli", "HEK293", "CHO"],
        "required": True
    },
    "temperature": {
        "type": float,
        "min": 4,
        "max": 40,
        "required": False
    }
}

# Validate properties
pp.set_property_schema("kinase_assay", property_schema)
is_valid = pp.validate_properties(properties_df, schema="kinase_assay")
```

## Statistical Analysis

### Basic Statistics

```python
# Get property statistics
stats = pp.get_property_statistics(
    property_name="IC50_nM",
    dataset_id="screening_results"
)
# {
#     'count': 150,
#     'mean': 45.3,
#     'std': 23.1,
#     'min': 0.1,
#     'max': 500.0,
#     'median': 32.5,
#     'q25': 12.3,
#     'q75': 67.8
# }

# Dataset summary
summary = pp.get_dataset_statistics("screening_results")
# Returns summary statistics for all numeric properties
```

### Distribution Analysis

```python
# Analyze property distribution
distribution = pp.analyze_distribution(
    property_name="pIC50",
    dataset_id="screening_results"
)

# Plot distribution
pp.plot_property_distribution(
    property_name="pIC50",
    dataset_id="screening_results",
    bins=30,
    output="pIC50_distribution.png"
)

# Compare distributions
pp.compare_property_distributions(
    property_name="Tm",
    dataset_ids=["wild_type", "mutants"],
    output="tm_comparison.png"
)
```

### Correlation Analysis

```python
# Calculate property correlations
correlations = pp.calculate_property_correlations(
    dataset_id="screening_results",
    properties=["pIC50", "MW", "logP", "TPSA"]
)

# Plot correlation matrix
pp.plot_correlation_matrix(
    correlations,
    output="property_correlations.png"
)

# Find correlated properties
correlated = pp.find_correlated_properties(
    dataset_id="screening_results",
    threshold=0.7  # Correlation coefficient threshold
)
```

## Integration with Other Processors

### Sequence Properties

```python
# Link properties to sequences
sp = SeqProcessor()
sequences = sp.load_sequences(["PROT1", "PROT2", "PROT3"])

# Calculate sequence properties
seq_properties = {}
for seq_id, sequence in sequences.items():
    seq_properties[seq_id] = {
        "length": len(sequence),
        "MW": sp.calculate_molecular_weight(sequence),
        "pI": sp.calculate_isoelectric_point(sequence),
        "gravy": sp.calculate_gravy(sequence)
    }

# Create property dataset
pp.create_property_dataset(
    dataset_id="sequence_properties",
    properties_df=pd.DataFrame(seq_properties).T,
    entity_column=None  # Index contains entity names
)
```

### Structure Properties

```python
# Link properties to structures
cp = CifBaseProcessor()
structures = cp.list_entities()

# Extract structural properties
struct_properties = []
for pdb_id in structures:
    metadata = cp.get_structure_metadata(pdb_id)
    struct_properties.append({
        "protein_id": pdb_id,
        "resolution": metadata.get("resolution"),
        "r_factor": metadata.get("r_factor"),
        "method": metadata.get("method"),
        "chains": len(metadata.get("chains", []))
    })

# Add to property processor
pp.create_property_dataset(
    dataset_id="structure_quality",
    properties_df=pd.DataFrame(struct_properties),
    entity_column="protein_id"
)
```

### GRN Properties

```python
# Link properties to GRN positions
gp = GRNBaseProcessor()
gp.load_grn_table("gpcr_alignment")

# Analyze position-specific properties
position_properties = {}
for position in ["3.50", "6.50", "7.50"]:
    residues = gp.get_position_column(position)
    position_properties[position] = {
        "conservation": gp.calculate_position_conservation(position),
        "hydrophobicity": np.mean([hydrophobicity[r] for r in residues if r != '-']),
        "charge": sum([1 for r in residues if r in 'RKH']) - sum([1 for r in residues if r in 'DE'])
    }

pp.create_property_dataset(
    dataset_id="grn_position_properties",
    properties_df=pd.DataFrame(position_properties).T,
    entity_column=None
)
```

## Advanced Features

### Multi-Dataset Operations

```python
# Merge property datasets
pp.merge_datasets(
    dataset_ids=["screen_1", "screen_2", "screen_3"],
    output_dataset_id="combined_screens",
    merge_on="entity",
    how="outer"  # Include all entities
)

# Compare properties across datasets
comparison = pp.compare_datasets(
    dataset_ids=["wild_type", "mutant"],
    property_name="activity",
    paired=True  # Same entities in both datasets
)

# Find dataset intersections
common_entities = pp.find_common_entities(
    dataset_ids=["structural_data", "functional_data"]
)
```

### Property Transformations

```python
# Transform properties
pp.transform_property(
    dataset_id="screening_results",
    property_name="IC50_nM",
    transform_fn=lambda x: -np.log10(x * 1e-9),  # Convert to pIC50
    new_property_name="pIC50"
)

# Normalize properties
pp.normalize_properties(
    dataset_id="multi_assay",
    properties=["assay_1", "assay_2", "assay_3"],
    method="z_score"  # or "min_max", "robust"
)

# Bin continuous properties
pp.bin_property(
    dataset_id="screening_results",
    property_name="pIC50",
    bins=[0, 5, 6, 7, 8, 10],
    labels=["inactive", "weak", "moderate", "potent", "very_potent"],
    new_property_name="activity_class"
)
```

### Export and Reporting

```python
# Export to various formats
pp.export_dataset(
    dataset_id="final_results",
    format="excel",
    output="results.xlsx",
    include_metadata=True
)

# Generate report
pp.generate_property_report(
    dataset_id="screening_results",
    output="screening_report.html",
    include_plots=True,
    include_statistics=True
)

# Export for machine learning
X, y, feature_names = pp.export_for_ml(
    dataset_id="qsar_data",
    target_property="pIC50",
    feature_properties=["MW", "logP", "TPSA", "HBD", "HBA"]
)
```

## Best Practices

### 1. Entity Naming Consistency

```python
# Ensure consistent entity names across datasets
# Map different naming conventions
name_mapping = {
    "EGFR": "EGFR_HUMAN",
    "P00533": "EGFR_HUMAN",
    "Epidermal growth factor receptor": "EGFR_HUMAN"
}

pp.standardize_entity_names(
    dataset_id="mixed_names",
    mapping=name_mapping
)
```

### 2. Property Units

```python
# Always specify units in property names
good_names = {
    "Kd_nM": "Dissociation constant in nM",
    "Tm_C": "Melting temperature in Celsius",
    "t_half_min": "Half-life in minutes"
}

# Convert units if needed
pp.convert_property_units(
    dataset_id="assay_data",
    property_name="Kd_uM",
    conversion_factor=1000,
    new_property_name="Kd_nM"
)
```

### 3. Missing Data Handling

```python
# Check for missing data
missing_report = pp.check_missing_data("screening_results")

# Handle missing values
pp.handle_missing_data(
    dataset_id="screening_results",
    strategy="drop",  # or "mean", "median", "forward_fill"
    subset=["IC50_nM", "Kd_nM"]  # Only these properties
)

# Impute missing values
pp.impute_missing_values(
    dataset_id="incomplete_data",
    method="knn",  # k-nearest neighbors
    n_neighbors=5
)
```

### 4. Data Validation

```python
# Set up validation rules
validation_rules = {
    "pH": lambda x: 0 <= x <= 14,
    "temperature_C": lambda x: -80 <= x <= 100,
    "concentration_mM": lambda x: x > 0,
    "percent_inhibition": lambda x: 0 <= x <= 100
}

# Validate data
issues = pp.validate_dataset(
    dataset_id="experimental_data",
    rules=validation_rules
)

# Fix or flag invalid entries
pp.flag_invalid_entries(issues, action="remove")  # or "flag", "correct"
```

## Common Workflows

### Workflow 1: SAR Analysis

```python
def structure_activity_analysis(pp, compound_properties):
    """Analyze structure-activity relationships."""
    
    # Create property dataset
    pp.create_property_dataset(
        dataset_id="sar_data",
        properties_df=compound_properties,
        entity_column="compound_id"
    )
    
    # Calculate derived properties
    pp.transform_property(
        dataset_id="sar_data",
        property_name="IC50_nM",
        transform_fn=lambda x: -np.log10(x * 1e-9),
        new_property_name="pIC50"
    )
    
    # Find activity cliffs
    activity_cliffs = pp.find_activity_cliffs(
        dataset_id="sar_data",
        activity_property="pIC50",
        similarity_property="tanimoto_similarity",
        cliff_threshold=1.0  # 1 log unit difference
    )
    
    # Identify key properties
    important_props = pp.find_important_properties(
        dataset_id="sar_data",
        target="pIC50",
        method="random_forest"
    )
    
    return {
        "activity_cliffs": activity_cliffs,
        "important_properties": important_props
    }
```

### Workflow 2: Mutant Characterization

```python
def characterize_mutants(pp, wild_type_id, mutant_data):
    """Compare mutant properties to wild-type."""
    
    # Load wild-type properties
    wt_props = pp.get_entity_properties(wild_type_id)
    
    results = []
    for mutant_id, mutant_props in mutant_data.items():
        # Calculate changes
        changes = {}
        for prop, value in mutant_props.items():
            if prop in wt_props:
                change = value - wt_props[prop]
                percent_change = (change / wt_props[prop]) * 100
                changes[f"{prop}_change"] = change
                changes[f"{prop}_percent_change"] = percent_change
        
        # Classify effect
        if "Tm_C_change" in changes:
            if changes["Tm_C_change"] > 5:
                stability_effect = "stabilizing"
            elif changes["Tm_C_change"] < -5:
                stability_effect = "destabilizing"
            else:
                stability_effect = "neutral"
        
        results.append({
            "mutant": mutant_id,
            **changes,
            "stability_effect": stability_effect
        })
    
    # Create results dataset
    pp.create_property_dataset(
        dataset_id="mutant_effects",
        properties_df=pd.DataFrame(results),
        entity_column="mutant"
    )
    
    return results
```

### Workflow 3: High-Throughput Screening

```python
def process_hts_data(pp, screening_file, control_wells):
    """Process high-throughput screening data."""
    
    # Load raw data
    raw_data = pd.read_csv(screening_file)
    
    # Calculate normalized values
    pos_control = raw_data[raw_data['well'].isin(control_wells['positive'])]['signal'].mean()
    neg_control = raw_data[raw_data['well'].isin(control_wells['negative'])]['signal'].mean()
    
    raw_data['percent_inhibition'] = (
        (neg_control - raw_data['signal']) / (neg_control - pos_control) * 100
    )
    
    # Calculate Z-prime
    z_prime = pp.calculate_z_prime(
        positive_controls=raw_data[raw_data['well'].isin(control_wells['positive'])]['signal'],
        negative_controls=raw_data[raw_data['well'].isin(control_wells['negative'])]['signal']
    )
    
    # Filter hits
    hits = raw_data[raw_data['percent_inhibition'] > 50]
    
    # Create property dataset
    pp.create_property_dataset(
        dataset_id="hts_results",
        properties_df=hits[['compound_id', 'percent_inhibition', 'signal']],
        entity_column="compound_id",
        metadata={"z_prime": z_prime, "date": screening_file.split('_')[1]}
    )
    
    # Statistical analysis
    hit_rate = len(hits) / len(raw_data) * 100
    
    return {
        "hits": hits['compound_id'].tolist(),
        "hit_rate": hit_rate,
        "z_prime": z_prime
    }
```

## Troubleshooting

### Common Issues

1. **Inconsistent entity names**
```python
# Standardize names before import
df['protein_id'] = df['protein_id'].str.upper().str.replace(' ', '_')
```

2. **Mixed data types**
```python
# Ensure consistent types
df['IC50'] = pd.to_numeric(df['IC50'], errors='coerce')
```

3. **Memory issues with large datasets**
```python
# Process in chunks
for chunk in pd.read_csv('large_file.csv', chunksize=10000):
    pp.append_to_dataset("large_dataset", chunk)
```

## Summary

PropertyProcessor provides:
- Flexible property and metadata management
- Property-based filtering and queries
- Statistical analysis capabilities
- Integration with all data formats
- Import/export functionality
- Validation and quality control

The processor enables comprehensive tracking and analysis of experimental data linked to biological entities across the Protos framework.