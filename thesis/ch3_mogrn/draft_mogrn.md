# MOGRN: A Generic Residue Numbering System for Microbial Rhodopsins

This chapter summarizes the research contribution presented in the accompanying publication: **"A Generic Residue-Numbering system for Microbial Rhodopsins"** (Hidber, Tajima, Kishi et al.). The full manuscript is included as part of this cumulative thesis.

Comparative analysis of microbial rhodopsins requires a way to identify equivalent positions across the superfamily. This is the problem that generic residue numbering systems solve. For GPCRs, including type II opsins, the Ballesteros-Weinstein system provides this foundation: position 3.50 refers to the conserved arginine of the DRY motif in any Class A receptor, and this common vocabulary enabled GPCRdb to integrate sequences, structures, ligands, and mutations into a unified resource. Microbial rhodopsins lacked an equivalent system. Each laboratory developed ad hoc conventions—some referenced bacteriorhodopsin positions, others channelrhodopsin-2, still others the C1C2 chimera—and comparing equivalent positions across publications required manual alignment for every analysis.

The difficulty is that sequence identity between major microbial rhodopsin classes falls as low as 10-12%. Ion pumps, channels, sensors, and enzyme-fused variants share the seven-transmembrane fold but have diverged extensively. Standard sequence alignment cannot reliably identify equivalent positions across such evolutionary distances. A numbering system for microbial rhodopsins must therefore be anchored to structure rather than sequence.

The MO (Microbial Opsin) numbering system uses the retinal-binding pocket as this structural anchor. Every microbial rhodopsin binds all-trans retinal via a Schiff base linkage to a conserved lysine on helix 7, and the geometry of this pocket remains consistent even when overall sequence identity falls below 20%. The Schiff base-forming lysine becomes position 7.50 by definition. Anchor positions for helices 1 through 6 are defined as the residue closest to retinal on each helix. Following the Ballesteros-Weinstein convention, each position receives a two-part number indicating helix and position relative to the anchor, with residues toward the N-terminus receiving lower numbers and those toward the C-terminus receiving higher numbers.

This system maps functionally important sites to standardized coordinates. The Schiff base counterions occupy positions 3.45 and 7.46. The TM3 motif at positions 3.45-3.49-3.56 determines ion pump specificity: outward proton pumps like bacteriorhodopsin carry DTD (aspartate-threonine-aspartate), inward proton pumps carry FSE, chloride pumps carry TSA or NTQ, and sodium pumps carry NDQ. Spectral tuning switches cluster at positions 3.53, 4.51, 6.54, and 7.49. The non-G rule controlling retinal planarity operates through position 4.54. Lateral fenestration sites for retinal uptake occur at positions 5.44 and 5.47. These patterns become immediately comparable in MO coordinates rather than requiring translation between protein-specific residue numbers.

The system accommodates non-canonical architectures. Heliorhodopsins have inverted membrane topology, yet the binding pocket geometry is preserved and MO coordinates apply without modification. Enzyme-fused rhodopsins have an additional TM0 helix, but the numbering annotates the seven-helix retinal-binding core. Gaps and insertions arising from local structural distortions receive systematic notation.

Validation against 129 structures spanning all major functional classes—proton pumps, chloride pumps, sodium pumps, cation channels, anion channels, sensory rhodopsins, and enzyme-coupled variants—demonstrated sub-angstrom accuracy (0.51-0.75 Å iRMSD) for the binding pocket region that defines the anchor positions. Application to approximately 40,000 non-redundant sequences from genomic and metagenomic sources revealed 31 sequence clusters, of which 14 contain characterized rhodopsins and 17 contain only uncharacterized sequences representing unexplored functional diversity.

Chapter 2 described how the GRNProcessor annotates Type II opsins: sequences align to GPCRdb reference tables, and the alignment transfers Ballesteros-Weinstein coordinates from structurally characterized references to new sequences. Type I opsins lacked equivalent infrastructure—no reference table existed, so no automated annotation was possible.

The MO numbering system provides this missing reference. The structural analysis of 129 microbial rhodopsins, spanning all major functional classes, produced a curated alignment where each column represents a structurally equivalent position and each row represents a characterized family member. This reference table is registered in ProtOS, enabling the same workflow for Type I that GPCRdb enables for Type II: query sequences align to the reference, and the GRNProcessor transfers MO coordinates based on alignment position.

Consider bacteriorhodopsin, the founding member of the microbial rhodopsin family. The crystal structure (PDB 1C3W) contains 248 residues forming seven transmembrane helices around a retinal chromophore. Key functional positions have been characterized through decades of mutagenesis, spectroscopy, and structural biology, yet the literature references these positions by sequence number—D85, T89, D96, D212, K216—which apply only to this specific protein. Given a channelrhodopsin or a heliorhodopsin, identifying equivalent positions requires manual structural alignment.

The GRNProcessor assigns MO coordinates automatically:

```python
from protos import ProtOS

protos = ProtOS()
registry = protos.registry
grn_proc = protos.grn_processor

# Fetch bacteriorhodopsin structure
entity = registry.fetch("1C3W")
protos.structure_processor.process(entity)

# Assign MO coordinates
grn_table = grn_proc.annotate_with_grn(
    entity_id="1C3W",
    grn_system="mogrn"
)
```

The output maps sequence positions to MO coordinates. Channelrhodopsin C1C2 (PDB 3UG9), despite sharing only 15% sequence identity with bacteriorhodopsin, can be annotated with the same system. The GRNProcessor assigns MO coordinates to both structures, enabling direct comparison:

```python
# Annotate both structures
grn_br = grn_proc.annotate_with_grn(entity_id="1C3W", grn_system="mogrn")
grn_c1c2 = grn_proc.annotate_with_grn(entity_id="3UG9", grn_system="mogrn")

# Compare equivalent positions
for mo_pos in ["3.45", "3.49", "3.56", "7.46", "7.50"]:
    br_res = grn_br.loc["1C3W", mo_pos]
    c1c2_res = grn_c1c2.loc["3UG9", mo_pos]
```

| MO Position | Function | HsBR (1C3W) | C1C2 (3UG9) |
|-------------|----------|-------------|-------------|
| 3.45 | Counterion | D85 | E162 |
| 3.49 | Proton transfer | T89 | T166 |
| 3.56 | Proton donor | D96 | H173 |
| 7.46 | Counterion | D212 | D292 |
| 7.50 | Schiff base | K216 | K296 |

The counterion at position 3.45 is aspartate in the proton pump (D85) and glutamate in the channel (E162)—both negatively charged, stabilizing the protonated Schiff base. The proton donor at 3.56 is aspartate in the pump (D96) but histidine in the channel (H173), reflecting their different mechanisms. The Schiff base lysine at 7.50 is conserved in both (K216, K296). These comparisons, which previously required manual structural alignment, become routine queries against the GRN table.

The TM3 motif at positions 3.45-3.49-3.56 distinguishes functional classes across the entire superfamily:

| Function | 3.45 | 3.49 | 3.56 | Example |
|----------|------|------|------|---------|
| Outward H⁺ pump | D | T | D | HsBR |
| Inward H⁺ pump (SzR) | F | S | E | SzR4 |
| Inward H⁺ pump (XeR) | D | T | A/L/S | NsXeR |
| Cl⁻ pump (archaeal) | T | S | A/D | NpHR |
| Cl⁻ pump (bacterial) | N | T | Q | NmClR |
| Na⁺ pump | N | D | Q | KR2 |

This pattern, described in the accompanying publication, becomes queryable once MO coordinates are assigned. Rather than manually extracting motifs from structural alignments, the question "what TM3 motif does this sequence have?" reduces to reading positions 3.45, 3.49, and 3.56 from the GRN table.

[FIGURE 3.1: MO numbering enables cross-family comparison. (A) Bacteriorhodopsin (PDB 1C3W, green) and channelrhodopsin C1C2 (PDB 3UG9, blue) superimposed on the retinal-binding pocket. Despite 15% sequence identity, both structures receive consistent MO coordinates. (B) Close-up comparing binding pocket residues at equivalent MO positions. The counterion (3.45) is D85 in HsBR and E162 in C1C2; the Schiff base lysine (7.50) is K216 and K296 respectively. Sequence positions differ by 77-80 residues, yet MO coordinates identify equivalent sites. (C) The binding pocket represented as a contact graph with MO-annotated nodes. This representation, with standardized labels across the superfamily, enables the cross-family learning described in Chapter 4.]

The critical feature of the reference table is its coverage. Aligning a channelrhodopsin sequence to bacteriorhodopsin alone produces unreliable coordinates because the sequences share only 15% identity—too divergent for confident alignment. But the MOGRN reference includes channelrhodopsin representatives, so alignment succeeds. The same holds for heliorhodopsins, enzyme rhodopsins, and other distant families: each is represented in the reference, so each can be reliably annotated. With standardized coordinates registered in the system, cross-family queries become routine: which residues occupy position 3.45 across all characterized channelrhodopsins, or how the TM3 motif varies across ion pump families, are questions answerable through single queries rather than manual compilation.

The binding pocket of bacteriorhodopsin contains 19 residues within 4Å of retinal, spanning helices 3 through 7. Once MO coordinates are assigned, these 19 positions form the nodes of a contact graph where edges connect spatially adjacent residues. Position 3.45 connects to 3.46, 7.46, and 7.50; position 7.50 connects to the retinal-contacting residues on helices 3, 6, and 7. This graph topology is conserved across microbial rhodopsins because the binding pocket geometry is conserved. What varies is the amino acid identity at each node—the substitutions that determine ion specificity, kinetics, and spectral properties.

This graph representation bridges MOGRN to LAMBDA. The standardized coordinates provide consistent node labels across the superfamily; the graph topology encodes the spatial relationships that determine function; protein language model embeddings capture the amino acid identity and sequence context at each position. Chapter 4 describes how a graph attention network learns from these representations to predict spectral properties across families that share less than 20% sequence identity.

The full methodology, detailed results, supplementary figures, and complete reference list are provided in the accompanying publication included in this thesis.
****