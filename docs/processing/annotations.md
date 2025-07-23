# Annotations

Protos provides comprehensive annotation capabilities for biological data, with a focus on Generic Residue Numbering (GRN) and functional annotations. The annotation system enables consistent position referencing, functional site identification, and cross-species comparison.

## Overview

Annotation types supported:
- **GRN annotations**: Standardized residue numbering for protein families
- **Functional annotations**: Active sites, binding sites, motifs
- **Secondary structure annotations**: Helices, sheets, loops
- **Domain annotations**: Protein domain boundaries and types
- **Conservation annotations**: Evolutionary conservation scores

## GRN Annotations

### Understanding GRN

Generic Residue Numbering (GRN) provides a universal coordinate system for protein families. See the comprehensive [GRN Documentation](../grn.md) for detailed format specifications.

```python
from protos.processing.grn import GRNProcessor
grnp = GRNProcessor(name="my_grn_processor")

# GRN formats:
# - TM positions: "3.50" (helix 3, position 50)
# - N-terminal: "n.10" (10 residues before TM1)
# - C-terminal: "c.5" (5 residues after last helix)
# - Loops: "34.001" (1st residue from helix 3 in loop 3-4)
```

### Creating GRN Annotations

```python
# Load reference GRN table
grnp.load_grn_table("gpcrdb_ref")  # or "mo_ref" for microbial opsins

# Annotate a new sequence
result = grnp.annotate_sequence(
    name="my_receptor",
    sequence="MNGTEGPNFYVPF...",
    protein_family="gpcr_a"
)

# The result includes:
# - Complete GRN assignments (TM regions, loops, N/C-terminals)
# - Missing positions that couldn't be annotated
# - Alignment details
```

### Detailed Annotation Process

```python
from protos.processing.grn.grn_table_utils import annotate_grnp

# Full annotation with verbosity
grn_row = annotate_grnp(
    grnp=grnp,
    query_name="NEW_RECEPTOR",
    query_seq="MNGTEGPNFYVPFSNKTGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRT",
    protein_family="gpcr_a",
    verbose=2  # Show detailed progress
)

# Key positions in GPCRs:
# - x.50 positions: Most conserved in each helix
# - 3.50: R in DRY motif
# - 7.50: N in NPxxY motif
# - 6.48: W in CWxP motif
```

### Working with GRN Tables

```python
# Access GRN data
grn_table = grnp.data  # pandas DataFrame

# Get sequence from GRN table
sequence = grnp.get_sequence("5HT2B_HUMAN")

# Get specific positions
dry_motif = grnp.get_grn_positions("5HT2B_HUMAN", ["3.49", "3.50", "3.51"])

# Filter by GRN range
tm3_residues = grnp.filter_by_grn_range("3.20", "3.60")
```

### GRN-Based Analysis

```python
# Find conserved positions
conserved_positions = gp.find_conserved_positions(
    grn_table,
    conservation_threshold=0.9
)

# Analyze specific motifs
motifs = {
    "DRY": ["3.49", "3.50", "3.51"],
    "NPxxY": ["7.49", "7.50", "7.51", "7.52", "7.53"],
    "E/DRY": ["3.49", "3.50", "3.51"],
    "CWxP": ["6.47", "6.48", "6.50"]
}

motif_conservation = gp.analyze_motif_conservation(grn_table, motifs)

# Find co-evolving positions
coevolution = gp.find_coevolving_positions(
    grn_table,
    method="mutual_information",
    threshold=0.8
)
```

## Functional Annotations

### Active Site Annotation

```python
from protos.processing.property.property_processor import PropertyProcessor
pp = PropertyProcessor()

# Define active site residues
active_site_annotation = {
    "KINASE_A": {
        "catalytic_residues": ["K72", "E91", "D166"],
        "atp_binding": ["G50", "G52", "K72", "E127"],
        "substrate_binding": ["D166", "N171", "D184"]
    }
}

# Store as property
pp.assign_property(
    entity_name="KINASE_A",
    property_name="active_site",
    value=active_site_annotation,
    dataset_id="functional_annotations"
)

# Map to GRN positions
grn_active_site = gp.map_residues_to_grn(
    sequence_id="KINASE_A",
    residues=active_site_annotation["catalytic_residues"]
)
```

### Binding Site Annotation

```python
# Annotate ligand binding sites
binding_sites = {
    "retinal_binding": {
        "residues": ["K296"],  # Schiff base
        "grn_positions": ["7.43"],
        "type": "covalent",
        "ligand": "retinal"
    },
    "g_protein_binding": {
        "residues": ["R135", "T136", "Y223", "Y306"],
        "grn_positions": ["3.50", "3.51", "5.58", "7.53"],
        "type": "protein-protein",
        "partner": "G-protein"
    }
}

# Analyze binding site conservation
binding_conservation = gp.analyze_site_conservation(
    grn_table,
    sites=binding_sites
)

# Find potential binding sites in new sequences
potential_sites = gp.predict_binding_sites(
    sequence="NEW_PROTEIN",
    known_sites=binding_sites,
    similarity_threshold=0.8
)
```

### Motif Annotations

```python
# Find and annotate sequence motifs
from protos.processing.sequence.seq_processor import SeqProcessor
sp = SeqProcessor()

# Search for known motifs
motif_annotations = sp.annotate_motifs(
    sequence_id="PROTEIN_X",
    databases=["prosite", "pfam"],
    e_value=0.001
)

# Custom motif patterns
custom_motifs = {
    "zinc_finger": "C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H",
    "nuclear_localization": "K[KR].[KR]",
    "glycosylation": "N[^P][ST]"
}

# Find custom motifs
for motif_name, pattern in custom_motifs.items():
    matches = sp.find_pattern(sequence, pattern)
    
    # Annotate matches
    for match in matches:
        annotation = {
            "type": "motif",
            "name": motif_name,
            "start": match.start(),
            "end": match.end(),
            "sequence": match.group()
        }
        sp.add_annotation(sequence_id, annotation)
```

## Secondary Structure Annotations

### Structure-Based Annotations

```python
from protos.processing.structure.struct_base_processor import CifBaseProcessor
cp = CifBaseProcessor()

# Assign secondary structure
structure = cp.load_structure("1ubq")
ss_annotation = cp.assign_secondary_structure(
    structure,
    method="dssp"  # or "stride"
)

# Map to sequence
sequence_ss = cp.map_ss_to_sequence("1ubq")
# Returns: "CCCHHHHHHHHHHHCCCCEEEEEEECCCCHHHHHHHHHH..."

# Annotate helices and sheets
structural_elements = cp.extract_structural_elements("1ubq")
# {
#     "helices": [(5, 17), (23, 34), (40, 51)],
#     "sheets": [(19, 22), (35, 39)],
#     "loops": [(1, 4), (18, 18), (52, 76)]
# }

# Map to GRN
grn_ss_mapping = gp.map_secondary_structure_to_grn(
    sequence_id="1ubq",
    ss_annotation=sequence_ss
)
```

### Transmembrane Annotations

```python
# Predict transmembrane regions
tm_predictions = sp.predict_transmembrane(
    sequence,
    method="tmhmm"  # or "phobius", "topcons"
)

# Annotate with GRN
tm_grn_annotation = {}
for i, tm_region in enumerate(tm_predictions["tm_regions"], 1):
    start, end = tm_region
    
    # Map to GRN helix
    tm_grn_annotation[f"TM{i}"] = {
        "sequence_positions": (start, end),
        "grn_helix": i,
        "grn_positions": gp.assign_grn_to_region(
            sequence[start:end],
            helix_number=i
        )
    }
```

## Domain Annotations

### Domain Identification

```python
# Find domains using multiple methods
domains = sp.annotate_domains(
    sequence_id="MULTI_DOMAIN_PROTEIN",
    methods=["pfam", "smart", "cdd"]
)

# Consolidate domain annotations
consolidated_domains = sp.consolidate_domain_annotations(
    domains,
    overlap_threshold=0.8
)

# Create domain architecture string
architecture = sp.create_domain_architecture(
    consolidated_domains
)
# Returns: "SH3-KINASE-SH2"

# Map domains to structure
if cp.has_entity("MULTI_DOMAIN_PROTEIN"):
    structural_domains = cp.map_domains_to_structure(
        "MULTI_DOMAIN_PROTEIN",
        consolidated_domains
    )
```

### Domain-Based GRN

```python
# Create domain-specific GRN
domain_grn = {}

for domain in consolidated_domains:
    domain_seq = sequence[domain["start"]:domain["end"]]
    domain_name = domain["name"]
    
    # Get domain-specific numbering
    if domain_name in ["Kinase", "Kinase_Tyr", "Pkinase"]:
        # Use kinase-specific numbering
        domain_grn[domain_name] = gp.apply_kinase_numbering(
            domain_seq,
            domain_type=domain_name
        )
    elif domain_name in ["SH2", "SH3"]:
        # Use SH domain numbering
        domain_grn[domain_name] = gp.apply_sh_numbering(
            domain_seq,
            domain_type=domain_name
        )
```

## Conservation Annotations

### Sequence Conservation

```python
# Calculate positional conservation
msa = sp.multiple_sequence_alignment(sequences)
conservation_scores = sp.calculate_conservation(
    msa,
    method="shannon_entropy"  # or "variance", "sum_of_pairs"
)

# Annotate highly conserved positions
conserved_positions = []
for pos, score in enumerate(conservation_scores):
    if score > 0.9:
        annotation = {
            "position": pos,
            "conservation_score": score,
            "type": "highly_conserved"
        }
        conserved_positions.append(annotation)

# Map conservation to GRN
grn_conservation = gp.map_conservation_to_grn(
    conservation_scores,
    sequence_id="REFERENCE_SEQ"
)
```

### Structural Conservation

```python
# Calculate structural conservation
structures = ["1abc", "2def", "3ghi"]
structural_conservation = cp.calculate_structural_conservation(
    structures,
    method="rmsd",  # or "tm_score", "gdt"
    radius=8.0
)

# Annotate structurally conserved regions
struct_conserved_regions = cp.find_conserved_regions(
    structural_conservation,
    threshold=0.8,
    min_length=5
)

# Combine sequence and structure conservation
combined_conservation = {}
for position in range(len(sequence)):
    combined_conservation[position] = {
        "sequence_conservation": conservation_scores[position],
        "structural_conservation": structural_conservation.get(position, None),
        "combined_score": (conservation_scores[position] + 
                          structural_conservation.get(position, 0)) / 2
    }
```

## Annotation Integration

### Cross-Processor Annotations

```python
# Comprehensive annotation pipeline
def annotate_protein_comprehensively(protein_id):
    """Generate all available annotations for a protein."""
    
    annotations = {
        "protein_id": protein_id,
        "annotations": {}
    }
    
    # Sequence annotations
    if sp.has_entity(protein_id):
        sequence = sp.load_sequence(protein_id)
        
        annotations["annotations"]["sequence"] = {
            "length": len(sequence),
            "composition": sp.calculate_composition(sequence),
            "motifs": sp.annotate_motifs(protein_id),
            "domains": sp.annotate_domains(protein_id)
        }
    
    # Structure annotations
    if cp.has_entity(protein_id):
        annotations["annotations"]["structure"] = {
            "secondary_structure": cp.get_secondary_structure(protein_id),
            "disorder": cp.predict_disorder(protein_id),
            "surface_accessibility": cp.calculate_sasa(protein_id)
        }
    
    # GRN annotations
    if gp.has_annotation(protein_id):
        annotations["annotations"]["grn"] = {
            "positions": gp.get_grn_annotation(protein_id),
            "family": gp.get_protein_family(protein_id),
            "conserved_motifs": gp.get_conserved_motifs(protein_id)
        }
    
    # Property annotations
    if pp.has_entity(protein_id):
        annotations["annotations"]["properties"] = pp.get_entity_properties(protein_id)
    
    return annotations
```

### Annotation Storage

```python
# Store annotations in standardized format
def store_annotations(annotations, output_format="json"):
    """Store annotations in various formats."""
    
    if output_format == "json":
        with open(f"{annotations['protein_id']}_annotations.json", 'w') as f:
            json.dump(annotations, f, indent=2)
    
    elif output_format == "gff":
        # Convert to GFF format
        gff_lines = []
        for ann_type, ann_data in annotations["annotations"].items():
            if "positions" in ann_data:
                for pos_data in ann_data["positions"]:
                    gff_line = format_gff_line(
                        seqid=annotations["protein_id"],
                        source=ann_type,
                        feature=pos_data["type"],
                        start=pos_data["start"],
                        end=pos_data["end"],
                        attributes=pos_data.get("attributes", {})
                    )
                    gff_lines.append(gff_line)
        
        with open(f"{annotations['protein_id']}.gff", 'w') as f:
            f.write('\n'.join(gff_lines))
```

## Visualization of Annotations

### Sequence Annotation Plots

```python
# Create annotation track visualization
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def visualize_annotations(sequence_id, annotations):
    """Create visual representation of annotations."""
    
    fig, axes = plt.subplots(len(annotations), 1, 
                            figsize=(12, 2*len(annotations)),
                            sharex=True)
    
    seq_length = annotations[0]["sequence_length"]
    
    for ax, (ann_type, ann_data) in zip(axes, annotations.items()):
        ax.set_ylim(0, 1)
        ax.set_xlim(0, seq_length)
        ax.set_ylabel(ann_type)
        
        # Draw annotations
        for annotation in ann_data:
            if "start" in annotation and "end" in annotation:
                rect = patches.Rectangle(
                    (annotation["start"], 0.2),
                    annotation["end"] - annotation["start"],
                    0.6,
                    facecolor=get_color_for_type(annotation["type"]),
                    edgecolor='black'
                )
                ax.add_patch(rect)
                
                # Add label
                ax.text(
                    (annotation["start"] + annotation["end"]) / 2,
                    0.5,
                    annotation.get("name", ""),
                    ha='center',
                    va='center'
                )
    
    plt.xlabel("Sequence Position")
    plt.tight_layout()
    plt.savefig(f"{sequence_id}_annotations.png")
```

### Structure Annotation Visualization

```python
# Generate PyMOL script for structure annotations
def create_pymol_annotation_script(structure_id, annotations):
    """Create PyMOL visualization script for annotations."""
    
    script_lines = [
        f"load {structure_id}.pdb",
        "hide everything",
        "show cartoon",
        "color gray"
    ]
    
    # Color by annotation type
    for ann_type, positions in annotations.items():
        color = get_pymol_color(ann_type)
        selection_name = ann_type.replace(" ", "_")
        
        # Create selection
        residue_list = "+".join([str(p) for p in positions])
        script_lines.append(
            f"select {selection_name}, resi {residue_list}"
        )
        
        # Color selection
        script_lines.append(f"color {color}, {selection_name}")
        
        # Show as sticks for important residues
        if ann_type in ["active_site", "binding_site"]:
            script_lines.append(f"show sticks, {selection_name}")
    
    # Save script
    with open(f"{structure_id}_annotations.pml", 'w') as f:
        f.write('\n'.join(script_lines))
```

## Best Practices

### 1. Annotation Consistency

```python
# Ensure consistent annotation across formats
def validate_annotation_consistency(protein_id):
    """Check annotation consistency across different sources."""
    
    issues = []
    
    # Get annotations from different sources
    seq_length = len(sp.load_sequence(protein_id))
    struct_length = cp.get_sequence_length(protein_id)
    
    if seq_length != struct_length:
        issues.append(f"Length mismatch: seq={seq_length}, struct={struct_length}")
    
    # Check domain boundaries
    seq_domains = sp.get_domains(protein_id)
    struct_domains = cp.get_domains(protein_id)
    
    # Compare domain positions
    for domain in seq_domains:
        matching_struct = find_matching_domain(domain, struct_domains)
        if not matching_struct:
            issues.append(f"Domain {domain['name']} not found in structure")
    
    return issues
```

### 2. Annotation Versioning

```python
# Track annotation versions
annotation_metadata = {
    "version": "2.0",
    "date": "2024-01-15",
    "methods": {
        "motifs": "ProSite 2023.2",
        "domains": "Pfam 35.0",
        "disorder": "IUPred2A",
        "conservation": "ConSurf"
    },
    "parameters": {
        "conservation_threshold": 0.7,
        "disorder_threshold": 0.5
    }
}

# Store with annotations
annotations["metadata"] = annotation_metadata
```

### 3. Annotation Updates

```python
# Update annotations when new data available
def update_annotations(protein_id, new_data_source):
    """Update existing annotations with new information."""
    
    # Load existing annotations
    existing = load_annotations(protein_id)
    
    # Get new annotations
    new_annotations = generate_annotations(protein_id, new_data_source)
    
    # Merge annotations
    updated = merge_annotations(
        existing,
        new_annotations,
        conflict_resolution="newest"  # or "highest_confidence"
    )
    
    # Track changes
    updated["changelog"] = {
        "date": datetime.now().isoformat(),
        "changes": compare_annotations(existing, new_annotations)
    }
    
    return updated
```

## Summary

Protos annotation system provides:
- Comprehensive GRN annotation framework
- Multi-level functional annotations
- Integration across data types
- Standardized storage formats
- Visualization capabilities
- Version tracking and updates

The annotation system enables rich biological context to be associated with sequences, structures, and other data types, facilitating advanced analysis and interpretation.