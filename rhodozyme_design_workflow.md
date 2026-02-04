# Rhodozyme Design Workflow

A light-activated enzyme created by grafting a serine protease catalytic triad (SER-HIS-ASP) onto rhodopsin's intracellular domain.

## Overview

This workflow uses RFdiffusion2 to design loop regions connecting fixed TM helices while preserving:
- The 7TM rhodopsin scaffold (fixed backbone)
- A grafted catalytic triad (SER-HIS-ASP)
- Retinal chromophore (RET)
- Tetrapeptide substrate (SB1-SB4)

## 1. Placement Structure Setup

### Source Structures
Input templates originate from BoltzGen constrained placements:
```
/data/fast/projects/protos/data/models/boltzgen/20260203_rhodozyme_constrained/
├── placement_00/  → Triad at positions 230, 245, 250
├── placement_02/  → Triad at positions 230, 242, 246
├── placement_06/  → Triad at positions 139, 227, 230
└── placement_07/  → Triad at positions 243, 247, 249
```

### Input Preparation

Each placement CIF was converted to PDB with:

1. **Duplicate removal**: BoltzGen outputs contained alternate conformations
2. **Substrate conversion**: The tetrapeptide substrate (X5P-ALA-PRO-OAR from 2AGE trypsin) was converted to HETATM ligands with ligand-style atom names:
   - X5P → SB1 (modified N-terminus)
   - ALA → SB2 (alanine)
   - PRO → SB3 (proline)
   - OAR → SB4 (modified arginine)
3. **ORI token**: Added at the centroid of triad CA atoms for RFdiffusion2 centering

Prepared input files:
```
/data/fast/projects/protos/data/models/rfdiffusion2/input/
├── placement_00_triad_ori.pdb
├── placement_02_triad_ori.pdb
├── placement_06_triad_ori.pdb
└── placement_07_triad_ori.pdb
```

## 2. RFdiffusion2 Configuration

### Fixed Regions (TM Helices)
| Region | Residues | Description |
|--------|----------|-------------|
| TM1 | 1-53 | N-terminus + TM1 |
| TM2-TM3 | 81-129 | TM2-TM3 bundle |
| TM4-TM5 | 162-225 | Extended to include full TM5 |
| TM7 | 261-301 | TM7 helix |

### Designed Regions (Loops)
| Region | Length | Description |
|--------|--------|-------------|
| ICL1 | 27 residues | Connects TM1 to TM2 |
| ICL2 | 32 residues | Connects TM3 to TM4 |
| ICL3 | Variable | Contains catalytic triad |
| H8 | 25 residues | C-terminal helix |

### Catalytic Triad Constraints

The triad is embedded within the ICL3 loop as fixed points:

**Placement 00** (triad at 230, 245, 250):
```
contigs: A1-53,27,A81-129,32,A162-225,4,A230,14,A245,4,A250,10,A261-301,25
contig_atoms: {A230:[OG,CB,CA],A245:[NE2,ND1,CE1,CD2,CG,CB,CA],A250:[OD1,OD2,CG,CB,CA]}
```

This specifies:
- SER230: Nucleophile (OG, CB, CA atoms constrained)
- HIS245: General base (imidazole ring atoms constrained)
- ASP250: Charge relay (carboxyl atoms constrained)

### Ligands
```
inference.ligand='RET,SB1,SB2,SB3,SB4'
```
- **RET**: Retinal chromophore (20 atoms)
- **SB1-SB4**: Tetrapeptide substrate (30 atoms total)

### Diffusion Parameters
```yaml
diffuser.T: 100                          # Diffusion timesteps
diffuser.preserve_motif_sidechains: True # Keep triad sidechains
```

## 3. Design Pipeline Scripts

### Step 1: Backbone Design (RFdiffusion2)

**Script**: `run_placement_XX_production.sh`

Generates 50 backbone designs per placement with fixed TM helices and catalytic triad.

```bash
# Run for each placement
cd /data/fast/projects/protos/data/models/rfdiffusion2
bash run_placement_00_production.sh  # ~3.3 hours
bash run_placement_02_production.sh
bash run_placement_06_production.sh
bash run_placement_07_production.sh
```

**Output**: `20260203_placement_XX_production/rhodozyme_XX_*.pdb`

### Step 2: Sequence Design + Validation (LigandMPNN + Chai)

**Script**: `run_ligandmpnn_chai.sh`

For each backbone design:
1. **LigandMPNN**: Generates 8 sequences per design (400 total per placement)
2. **Chai**: Validates that designed sequences fold correctly with ligands

```bash
# Run after RFdiffusion2 completes
bash run_ligandmpnn_chai.sh 00
bash run_ligandmpnn_chai.sh 02
bash run_ligandmpnn_chai.sh 06
bash run_ligandmpnn_chai.sh 07
```

**Output**:
- `20260203_placement_XX_mpnn_chai/ligmpnn/seqs/` - FASTA sequences
- `20260203_placement_XX_mpnn_chai/ligmpnn/packed/` - Full-atom structures
- `20260203_placement_XX_mpnn_chai/chai/out/` - Chai predictions and scores

## 4. Expected Results

| Placement | Backbone Designs | Sequences | Chai Predictions |
|-----------|------------------|-----------|------------------|
| 00 | 50 | 400 | 2000 (5 models each) |
| 02 | 50 | 400 | 2000 |
| 06 | 50 | 400 | 2000 |
| 07 | 50 | 400 | 2000 |
| **Total** | **200** | **1600** | **8000** |

## 5. Selection Criteria

After Chai validation, select designs based on:
1. **pLDDT > 70**: Confident structure prediction
2. **pTM > 0.7**: Good template modeling score
3. **Triad geometry**: SER-OG to HIS-NE2 distance 2.8-3.2 Å
4. **Ligand contacts**: Substrate positioned near triad
5. **TM helix RMSD < 1 Å**: Preserved rhodopsin fold

## File Structure

```
/data/fast/projects/protos/data/models/rfdiffusion2/
├── input/
│   ├── placement_00_triad_ori.pdb
│   ├── placement_02_triad_ori.pdb
│   ├── placement_06_triad_ori.pdb
│   └── placement_07_triad_ori.pdb
├── run_placement_00_debug.sh
├── run_placement_00_production.sh
├── run_placement_02_production.sh
├── run_placement_06_production.sh
├── run_placement_07_production.sh
├── run_ligandmpnn_chai.sh
├── 20260203_placement_00_production/  # RFdiffusion2 outputs
├── 20260203_placement_02_production/
├── 20260203_placement_06_production/
├── 20260203_placement_07_production/
└── 20260203_placement_XX_mpnn_chai/   # LigandMPNN + Chai outputs
```

## References

- RFdiffusion2: All-atom protein structure diffusion
- LigandMPNN: Ligand-aware protein sequence design
- Chai: Structure prediction with ligand support
- BoltzGen: Boltzmann generator for protein conformations
