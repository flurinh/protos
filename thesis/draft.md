# A Brief History of Protein Data

For most of history, proteins were invisible. Biochemists knew they existed, knew they catalyzed reactions and transmitted signals, but could not see them. X-ray crystallography changed this. In 1958, John Kendrew solved the structure of myoglobin [@kendrew1958]. For the first time, researchers could see how a protein folded, a chain of amino acids wrapped into a compact globule, with a heme group nestled in a hydrophobic pocket where oxygen binds. The structure revealed why myoglobin worked.

Myoglobin was soluble, floating freely in the cytoplasm. Membrane proteins proved far harder. These proteins (the channels, pumps, and receptors embedded in cell membranes) have hydrophobic surfaces that contact lipid and hydrophilic surfaces exposed to water. Removed from the membrane, they aggregate, and aggregated proteins cannot form the ordered crystals that X-ray crystallography requires. Membrane protein structures remained out of reach for decades. Yet these proteins handle nearly every interaction between a cell and its environment, making them essential targets for understanding biology and disease.

Despite these challenges, the first view of a membrane protein came in 1975. Richard Henderson and Nigel Unwin used electron microscopy to image bacteriorhodopsin, a light-driven proton pump from *Halobacterium salinarum* [@henderson1975]. The resolution was 7 Å, too coarse to resolve individual atoms, but the map showed seven rod-like densities spanning the membrane. Bacteriorhodopsin was built from seven transmembrane helices.

Bacteriorhodopsin was a microbial protein, but the seven-helix architecture also appeared in an important class of animal membrane proteins. G protein-coupled receptors would prove especially important for drug discovery. In 1993, Schertler, Villa, and Henderson published a projection map of bovine rhodopsin, the light-sensitive protein of vertebrate eyes and the first GPCR amenable to structural analysis [@schertler1993]. The 9 Å map showed that rhodopsin, like bacteriorhodopsin, contained seven transmembrane helices. In 2000, Krzysztof Palczewski solved bovine rhodopsin at 2.8 Å, the first GPCR at atomic resolution [@palczewski2000], confirming the seven-helix architecture and showing how retinal, the light-absorbing chromophore, sits in the transmembrane core. In 2007, Brian Kobilka's group solved the β₂-adrenergic receptor [@rasmussen2007], showing that other GPCRs could be crystallized. Kobilka shared the 2012 Nobel Prize with Robert Lefkowitz for this work.

Crystallography required proteins to form ordered crystals, a barrier that excluded many membrane proteins. Cryo-electron microscopy removed this requirement. The technique flash-freezes proteins in thin ice, preserving their native conformations without the distortions that crystallization can introduce. For decades, cryo-EM produced only blurry images, useful for overall shape but not atomic detail. The breakthrough came around 2013 with direct electron detectors, which dramatically improved signal quality, and new algorithms that could align and average millions of particle images. Resolution improved from 10-20 Å to 3 Å, then 2 Å. Membrane protein complexes that would not crystallize (large, flexible, or transient assemblies) could now be imaged at near-atomic resolution. Henderson shared the 2017 Nobel Prize for this work.

These advances transformed structural biology. The Protein Data Bank grew from roughly 1,000 structures in 1993 to over 200,000 by 2024, with GPCR structures alone rising from fewer than 5 in 2007 to over 700 today.

Structure determination reveals individual proteins, while sequencing reveals how many different proteins exist. Sanger sequencing gave way to genomics, then metagenomics, the sequencing of DNA extracted directly from environmental samples. UniProt now contains over 200 million protein sequences, many from ocean surveys, soil samples, and gut microbiomes. For some protein families, metagenomics revealed unexpected abundance. Microbial rhodopsins offer a striking example. What had seemed a peculiarity of salt-loving archaea turned out to be globally abundant. Proteorhodopsins, discovered in 2000 from Pacific Ocean samples [@beja2000], are now recognized as among the most abundant proteins on Earth, and metagenomics has revealed thousands of rhodopsin sequences across bacteria, archaea, algae, and fungi, a diversity invisible to traditional microbiology.

Yet sequence diversity far exceeds structural coverage. For every structure in the PDB, roughly 1,000 sequences have no structural information.

In 2020, DeepMind's AlphaFold2 addressed this disparity [@jumper2021]. The system predicts protein structures from sequence with accuracy matching experiment, and within two years, DeepMind and EMBL-EBI released predicted structures for over 200 million proteins, covering the entire known sequence universe [@varadi2022]. Structure prediction continues to advance. Systems like Boltz [@wohlwend2024] now predict proteins bound to their ligands—small molecules, cofactors, and interaction partners—enabling computational analysis of binding pockets without requiring experimental structures of each complex.

Alongside structure prediction, protein language models have transformed how researchers extract information from sequences. These models train on millions of protein sequences, learning statistical patterns that reflect evolutionary constraints and structural relationships. Models like ESM-2 [@lin2023] and Ankh [@elnaggar2023] produce embeddings—vectors of hundreds or thousands of numbers that summarize what the model has learned about a sequence. Unlike hand-crafted features, embeddings emerge from the training process itself, capturing patterns that may not be obvious to human researchers. Two proteins with diverged sequences but similar functions often have similar embeddings, even when their sequences share less than 20% identity. This property makes embeddings valuable as input features for machine learning models that predict protein properties, as explored in the Embedding Processor section of Chapter 2.

With these advances, the field shifted from data scarcity to data abundance, and every known protein sequence now has a predicted structure.

Yet predicted structures lack context. AlphaFold predictions come with confidence scores and little else, lacking ligands, experimental conditions, and functional annotations. Comparing predictions to experimental structures requires careful handling of coordinate systems and numbering conventions. Annotating millions of structures needs automation that does not yet exist.

The bottleneck has shifted from generating data to understanding it.

For researchers studying membrane protein families like GPCRs and microbial rhodopsins, this shift has practical consequences. These families now have hundreds of structures and are targets for drug discovery and biotechnology, yet systematic comparison across their members remains difficult because it requires knowing which positions in one protein correspond to which positions in another.
# Comparing Proteins

Understanding protein function requires more than accumulating structures; it requires comparing them. Data abundance enables diverse approaches to this challenge. For protein superfamilies—large families sharing a common architectural fold—one particularly powerful approach is systematic comparison across members. Hundreds of related proteins, each with slightly different sequences and functions, create natural experiments. Comparing which residues vary and which are conserved reveals the molecular basis of functional differences.

Such comparison requires a common reference frame. Position 85 in one protein may be functionally equivalent to position 113 in another, or it may not be. Sequence alignments help, but proteins in the same superfamily can share less than 20% sequence identity, too little for reliable alignment. Structure-based alignment is more reliable but requires structures for both proteins. The central problem is determining which positions in different proteins correspond to one another.

G protein-coupled receptors solved this problem in 1995. Ballesteros and Weinstein introduced a numbering system that assigns each position a two-part identifier consisting of the helix number and position relative to the most conserved residue on that helix [@ballesteros1995]. The format is X.50 for the most conserved residue on helix X, with other positions numbered relative to this anchor. Position 3.32 indicates helix 3, at the position numbered 32 (eighteen residues before the conserved position 50). This scheme translates across the entire GPCR superfamily, roughly 800 members in the human genome.

GPCRdb built on this foundation [@kooistra2021]. The database integrates sequences, structures, ligands, and mutations using Ballesteros-Weinstein positions as the common coordinate system, enabling queries such as identifying the residue at position 3.32 across all aminergic receptors within seconds.

This numbering system shaped how an entire field communicates. When a researcher reports that a mutation at position 3.32 affects ligand binding, every GPCR researcher immediately understands where that position is and can relate the finding to their own work. Functional insights accumulate across laboratories because everyone uses the same coordinate system. The numbering became the language of the field.

This standardization also enabled discovery. A systematic analysis of GPCR structures found that despite less than 20% sequence identity, all GPCRs share a conserved network of residue contacts [@venkatakrishnan2013]. The contacts form the structural scaffold that defines the GPCR fold. Without standardized positions, this pattern would have been invisible, with each structure analyzed in isolation and contacts described in local numbering that could not be compared.

The impact on drug discovery was substantial. Approximately 30% of approved drugs target GPCRs. With standardized positions, binding sites can be compared across receptors. A mutation in one receptor can be mapped to the equivalent position in another. If a drug binds at position 3.32 in one receptor, researchers can immediately ask what residue sits at 3.32 in related receptors and whether side effects might arise from off-target binding. Structure-based design became routine.

The GPCR story illustrates what a common coordinate system enables. This thesis examines another seven-transmembrane protein family where such standardization does not yet exist.
# Opsins

Opsins are light-sensitive proteins that use retinal as their chromophore. They exist across all domains of life (archaea, bacteria, algae, fungi, and animals) where they perform diverse functions including pumping ions, opening channels, sensing light direction, and activating signaling cascades.

All opsins bind retinal, a vitamin A derivative, through a covalent Schiff base linkage (a bond formed between retinal's aldehyde group and the amino group of a lysine) to a conserved lysine residue. When light hits the retinal, it isomerizes, triggering a conformational change in the protein that serves as the signal. This bond, and the protein environment around it, determines how the opsin responds to light. All opsins share a seven-transmembrane helix architecture and the same basic mechanism of retinal isomerization coupled to conformational change.

Despite sharing a chromophore and architecture, opsins fall into two evolutionarily distinct families.

Type I opsins, the microbial rhodopsins, are found in archaea, bacteria, algae, and fungi. These are not GPCRs and share no sequence homology with animal opsins, adopting a distinct fold despite sharing the seven-transmembrane architecture. Retinal binds in the all-trans configuration, and light triggers isomerization to 13-cis. The protein responds directly by pumping ions, opening a channel, or activating an enzymatic domain, with no G protein involved.

The founding member is bacteriorhodopsin, discovered in 1971 in *Halobacterium salinarum* [@oesterhelt1971], which pumps protons across the membrane to generate ATP in the absence of oxygen. The family has since expanded to include ion pumps with different specificities, light-gated channels, sensory receptors, and light-activated enzymes. Channelrhodopsins, which open cation channels in response to light, enabled optogenetics—the use of light to control genetically targeted cells—and transformed neuroscience [@nagel2003]. Despite this functional diversity, all Type I opsins share the same architecture, and the residues surrounding the chromophore determine their specific activities.

Type II opsins, the animal opsins, are found exclusively in animals. These are GPCRs, G protein-coupled receptors that signal through heterotrimeric G proteins. Retinal binds in the 11-cis configuration, and light triggers isomerization to all-trans. The conformational change activates a bound G protein, initiating a signaling cascade; the protein itself does not pump ions or open channels. The founding member is bovine rhodopsin, the photoreceptor in rod cells that enables dim-light vision. Cone opsins, which mediate color vision, are closely related and exist in multiple spectral variants. The spectral tuning of these opsins determines which wavelengths each type detects.

Type I and Type II opsins are not homologous. They share no detectable sequence similarity outside the retinal-binding lysine, and their seven-transmembrane folds are structurally distinct. These are analogous structures, not homologous ones. Two separate protein lineages faced the same problem of detecting light and arrived at the same solution, binding retinal via a protonated Schiff base (one that carries a positive charge) and using light-triggered isomerization to drive conformational change.

This shared solution has an important consequence. Both families use the same chromophore in geometrically similar binding pockets. In bacteriorhodopsin, retinal binds to lysine 216 (K216) on helix G, while in bovine rhodopsin, retinal binds to lysine 296 (K296) on helix VII. In both cases, the chromophore sits in a pocket formed by transmembrane helices, the Schiff base is protonated, and a negatively charged residue (aspartate or glutamate, called the counterion) sits nearby to stabilize the positive charge (Figure 1.1). Despite independent origins, both families converged on shared geometry because that geometry solves a physical problem, namely tuning retinal's electronic properties to absorb light at specific wavelengths.

[FIGURE 1.1: Structural comparison of Type I and Type II opsins. (A) Bacteriorhodopsin (PDB 1C3W), a microbial rhodopsin (Type I), shown in ribbon representation with retinal in yellow. (B) Bovine rhodopsin (PDB 1U19), an animal opsin (Type II/GPCR), shown in the same orientation. (C) Close-up overlay of the retinal binding pockets, highlighting the conserved features despite independent evolutionary origins. The retinal chromophore (yellow) is covalently bound to a lysine residue via a protonated Schiff base, and the negatively charged counterion (D85 in bacteriorhodopsin, E113 in bovine rhodopsin, shown in red) is positioned to stabilize the positive charge. The geometric similarity of these binding pockets, arrived at independently, enables cross-family analysis of spectral tuning.]

If local environment determines function, and both families have similar local environments, then functional principles might transfer across families. This premise enables cross-family analysis. The spectral properties of an opsin, meaning the wavelengths of light it absorbs, depend on the residues surrounding the chromophore, and if those residues occupy similar positions in both families, the rules governing spectral tuning might be shared.

Such analysis requires consistent positional annotation, knowing which residues in one protein correspond to which in another. Type II opsins, as GPCRs, benefit from the Ballesteros-Weinstein numbering system described in the previous section, and researchers can compare binding pocket residues across animal opsins using standardized positions. Type I opsins have no equivalent system. This matters because microbial rhodopsins are both abundant and widely used. Over 100 structures exist in the PDB and thousands of sequences are known from metagenomic surveys. Without standardized positions, each laboratory uses its own conventions. Some reference bacteriorhodopsin, others channelrhodopsin. Papers describe "the position equivalent to D85 in bR" without a standard way to identify that position across the family, and the literature remains fragmented.
# Spectral Tuning

Every opsin has a characteristic absorption spectrum, and the peak of that spectrum, the wavelength at which the protein absorbs light most strongly, is called λmax. This parameter determines what color of light activates the opsin. An opsin with λmax at 480 nm responds best to blue light, one at 560 nm responds to green, and one at 620 nm responds to red. Across the opsin superfamily, λmax ranges from approximately 350 nm (ultraviolet) to 650 nm (far-red), a span of 300 nm achieved using the same chromophore (Figure 1.2).

[FIGURE 1.2: The spectral tuning problem. (A) The UV-visible spectrum (350-700 nm) with representative opsins positioned at their λmax values, illustrating the 300 nm range achieved by the opsin superfamily. UV-sensitive opsins (λmax ~350-400 nm) include SWS1 cone opsins with deprotonated Schiff bases; blue-absorbing opsins (λmax ~400-480 nm) include short-wavelength cone opsins; green-absorbing opsins (λmax ~500-550 nm) include rhodopsin and many microbial rhodopsins; red-shifted opsins (λmax ~560-650 nm) include long-wavelength cone opsins and engineered channelrhodopsin variants. (B) Chemical structure of retinal showing the conjugated polyene chain responsible for light absorption. The aldehyde group forms a Schiff base linkage with the ε-amino group of a conserved lysine. (C) Schematic of the key determinants of spectral tuning. The protonated Schiff base carries a positive charge that is stabilized by a nearby counterion (Asp or Glu). The distance between these charges, along with other binding pocket residues, determines the distribution of electron density along the conjugated chain and thus the energy gap for light absorption.]

Spectral properties matter for applications. In optogenetics, controlling which wavelengths activate which opsins enables precise manipulation. Red-shifted opsins (λmax > 600 nm) allow deeper tissue penetration because red light travels further through tissue than blue. Blue-shifted opsins (λmax < 450 nm) enable spectral separation, allowing multiple opsins with different λmax values to be activated independently in the same tissue. Engineering opsins with specific spectral properties requires understanding what determines λmax.

The chromophore itself provides only part of the answer. Retinal in solution absorbs at approximately 440 nm, yet when bound to different opsins, absorption shifts by up to 200 nm in either direction. The protein environment tunes the chromophore, and understanding how evolution achieves 300 nm of spectral range using a single chromophore requires examining the local environment around the retinal.

Retinal attaches to the protein through a Schiff base linkage to a conserved lysine (K216 in bacteriorhodopsin, K296 in bovine rhodopsin). In most opsins, this Schiff base is protonated, meaning it carries a positive charge. The protonation state is critical: protonated opsins absorb in the visible range (above 400 nm), while deprotonated opsins absorb in the ultraviolet (below 400 nm), as seen in UV-sensitive SWS1 cone opsins. The protonated Schiff base creates a positive charge in the binding pocket. How the protein environment stabilizes this charge determines where the electron density in retinal's conjugated system sits, and that determines λmax.

Retinal is a polyene, a chain of alternating single and double bonds. The electrons in these double bonds are delocalized across the chain, forming a π-electron system. Light absorption promotes an electron from a lower-energy state to a higher-energy state. The energy gap between these states determines the wavelength of light absorbed. A smaller gap means lower energy, which corresponds to longer wavelengths (red shift). A larger gap means higher energy, which corresponds to shorter wavelengths (blue shift).

A negatively charged residue (aspartate or glutamate) sits near the protonated Schiff base and stabilizes the positive charge. This residue is called the counterion. In bacteriorhodopsin, the counterion is D85. In bovine rhodopsin, it is E113. The distance between the counterion and the protonated Schiff base affects λmax. A closer counterion stabilizes the positive charge more strongly, localizing it near the Schiff base. This increases the energy gap and shifts absorption to shorter wavelengths (blue shift). A more distant counterion allows the positive charge to delocalize along the polyene chain, decreasing the energy gap and shifting absorption to longer wavelengths (red shift).

Other residues in the binding pocket contribute. Polar residues can form hydrogen bonds with retinal or with water molecules in the pocket. Aromatic residues provide steric constraints. Charged residues elsewhere in the pocket create an electrostatic field that influences how electrons distribute along the conjugated chain. The cumulative effect of all these interactions (counterion distance, hydrogen bonding, electrostatics, sterics) determines λmax.

Single mutations at key positions produce large spectral shifts. The D85N mutation in bacteriorhodopsin replaces the negatively charged aspartate counterion with neutral asparagine, producing a red shift of approximately 40 nm and changing the protein from purple to yellow. The equivalent mutation in bovine rhodopsin, E113Q, produces a similar effect. The counterion position is a major determinant of λmax, though not the only one; other binding pocket residues contribute individual shifts of 10-30 nm.

These examples illustrate that understanding spectral tuning requires identifying which residues matter and how they vary across proteins. Studying this systematically requires comparing binding pocket residues, which returns to the standardization problem identified in the previous section. Beyond comparison, researchers seek prediction. The ability to predict λmax from sequence alone would enable rapid screening of metagenomic discoveries and guide engineering efforts. Family-specific methods exist for this purpose. OPTICS [@frazer2025] predicts λmax for type II opsins using sequence features and phylogenetic information trained on the VPOD dataset [@davis2024], achieving approximately 5.5 nm mean absolute error, but working only within the GPCR family. For microbial opsins, Inoue et al. [@inoue2021] developed a LASSO regression model using binding pocket residue identities, achieving 7.8 nm MAE on their test set. RhoMax [@gerstenbruch2024], the most recent approach, uses a graph neural network on type I opsin binding pocket structures to achieve 6.8 nm MAE. However, no method crosses the Type I / Type II divide. Each is trained on one family and fails when applied to the other.

This limitation is notable. If local chromophore environment determines λmax, and both families have similar local environments, then a model trained on both families should perform at least as well as family-specific models, and potentially better with more training data. Yet no method currently predicts λmax across opsin families from sequence alone.

Two research questions emerge from this context.

The previous section identified the first: how can standardized positions be established across microbial rhodopsins? The systematic comparison that enabled GPCR research requires a common coordinate system. Microbial rhodopsins need the same.

The analysis above raises the second: can cross-family learning improve spectral property prediction? The shared binding pocket geometry suggests that rules governing spectral tuning might transfer across the Type I / Type II divide. Testing this hypothesis requires a model that learns from both families.

Addressing both questions requires capabilities that current tools do not provide, specifically linking sequences, structures, and functional properties across databases while maintaining consistent positional annotation.
# Fragmentation in Bioinformatics

Answering the questions raised in the previous section requires linking multiple data types. These include sequences from UniProt and metagenomic databases, structures from the PDB and AlphaFold DB, properties like measured λmax values scattered across literature, embeddings from protein language models, and annotations such as standardized positions and binding pocket definitions. No single database contains all of this. Researchers must download files from multiple sources, convert between formats, and track which identifiers refer to the same protein.

The tools that exist are refined individually but do not work together. Format proliferation presents one challenge, with structures stored in PDB or mmCIF format, sequences in FASTA, and properties in CSV or JSON, each tool expecting its own input format. Identity fragmentation presents another challenge, as the same protein has different identifiers in different databases. UniProt calls it P02945, the PDB calls the structure 1C3W, and AlphaFold calls its prediction AF-P02945-F1. Linking these requires manual lookup or custom scripts.

A typical analysis illustrates these challenges. Comparing spectral properties across 100 microbial rhodopsins requires downloading sequences from UniProt, structures from the PDB, and AlphaFold predictions for sequences without experimental structures, while also collecting λmax values from the literature. Determining which binding pocket residues correlate with spectral shifts then requires aligning all sequences, mapping sequence positions to structure positions, extracting the coordinates of binding pocket residues, and correlating those features with λmax. Each step needs a different tool and a different file format, connected by custom scripts. The glue code (parsing, format conversion, identifier mapping) dominates the work, and the same code gets written repeatedly across projects.

Two groups are limited by this fragmentation. Experimentalists generate data but cannot analyze it computationally. A researcher who measures λmax values for 50 opsins cannot easily compare those values to structural features without programming skills or a collaborator who has them. Bioinformaticians can analyze data but spend disproportionate time on data engineering. Loading a structure, aligning it, extracting binding pocket residues, and fetching embeddings are each achievable, but chaining them requires effort that does not contribute to scientific insight. Both communities are limited by tools, not by ideas.

No framework manages multi-modal protein data with persistent identity and composable operations. Establishing standardized positions across microbial rhodopsins requires assembling and aligning structures from multiple sources, a task currently requiring custom scripts for each analysis. Testing whether cross-family learning improves spectral prediction requires linking thousands of sequences to their measured λmax values, embeddings, and structural features. These are not extraordinary requirements, but no existing tool addresses them.

A further limitation cuts across the others. Accessibility requires programming. Even with good tools, researchers without coding skills cannot use them. An experimentalist who measures λmax values cannot computationally compare them to structural features without learning to program or finding a collaborator who already knows. This barrier excludes much of the experimental community from computational analysis of their own data.
# Thesis Contributions

The preceding sections identified three gaps. Microbial rhodopsins lack the standardized positional annotation that enabled systematic GPCR comparison. No method predicts spectral properties across the type I and type II divide. And existing tools do not compose into workflows that link sequences, structures, and functional properties while maintaining consistent identity. This thesis addresses all three through an integrated system.

**ProtOS** provides the data infrastructure that underlies the other contributions. Protein research requires linking sequences, structures, functional properties, and computational representations across databases that use different identifiers, formats, and conventions. ProtOS manages these relationships through a processor architecture that handles specific data types (sequences, structures, embeddings, properties, residue annotations) through a consistent interface. The framework enables workflows that span data modalities without the glue code that currently dominates computational protein analysis. ProtOS-MCP extends this through the Model Context Protocol, enabling natural language access to protein data operations—an experimentalist can ask which residues contact retinal and receive an answer without writing code.

**MOGRN** (Microbial Opsin Generic Residue Numbering) establishes standardized positional annotation for microbial rhodopsins. Overall sequence identity across major classes falls below 20%, yet the retinal-binding pocket geometry remains conserved. MOGRN introduces a structure-based numbering system anchored to the retinal-binding site, following the Ballesteros-Weinstein convention used for GPCRs. The Schiff base lysine becomes position 7.50 by definition; other positions are numbered relative to helix-specific anchors. Validated against 129 structures spanning the functional diversity of microbial rhodopsins, the system provides a common coordinate system that enables systematic comparison across the superfamily.

**LAMBDA** (Light Absorption Modeling through Binding Domain Analysis) tests whether shared binding pocket geometry enables cross-family learning for spectral prediction. Existing methods work within single opsin families but do not transfer across the type I / type II divide. LAMBDA represents binding pockets as graphs aligned on retinal coordinates, enabling a single model to learn from both families. This is possible because both now have standardized positional annotation: MOGRN for type I, Ballesteros-Weinstein for type II. Trained on 2,120 sequences spanning both opsin families and hCRBPII (a lipocalin fold), LAMBDA achieves state-of-the-art accuracy: 5.18 nm mean absolute error on type II opsins and 6.86 nm on type I, while providing cross-family capability that existing methods lack. The Opsin Atlas extends these predictions to 47,700 sequences, enabling systematic identification of candidates with desired spectral properties.

These contributions integrate into a system. ProtOS manages the data that MOGRN annotates and LAMBDA consumes. The Opsin Atlas exists because ProtOS could route 47,700 sequences through GRN annotation and spectral prediction. The following chapters present each contribution in detail, beginning with ProtOS.
# ProtOS

The previous chapter identified three gaps in opsin research: microbial rhodopsins lack the standardized positional annotation that enabled systematic GPCR comparison, no method predicts spectral properties across the Type I / Type II divide, and existing tools do not compose into workflows that link sequences, structures, and functional properties while maintaining consistent identity. This chapter addresses the third gap through ProtOS, a framework for protein data access and management that also provides the infrastructure underlying the research contributions that follow.

Protein research requires data from multiple sources, and each source uses its own conventions. Bacteriorhodopsin appears as P02945 in UniProt, 1C3W in the Protein Data Bank (one of over fifty experimental structures), and AF-P02945-F1 in AlphaFold DB. Assembling data for even a single protein means navigating separate databases, downloading files in different formats, and tracking which data belongs together. The burden grows with scope. Comparing binding pockets across 129 microbial rhodopsin structures, the validation set for MOGRN, requires hundreds of downloads. Collecting the 1,799 sequences with measured absorption maxima for LAMBDA training means reconciling identifiers across sources and linking measurements scattered through decades of literature. Each step is simple; collectively, they dominate the work.

ProtOS addresses this through two connected capabilities. Unified data access provides a single interface to major protein databases, routing requests by identifier format and caching results locally. The entity system maintains identity across data types, so that a protein loaded from the PDB, annotated with standardized positions, enriched with embeddings, and linked to measured properties remains a single coherent object rather than a collection of disconnected files.

Six processors handle the data types opsin research requires: sequences, structures, residue contact graphs, standardized position annotations, embeddings, and properties. Each processor manages one data modality. Outputs from one become inputs to another through persistent identity, enabling workflows that span modalities without the format conversion and identifier reconciliation that currently fragment computational biology.

The following sections present the framework through examples drawn from opsin research. These examples are not demonstrations of ProtOS capabilities in the abstract; they reveal pieces of the spectral tuning puzzle, building toward the research contributions in subsequent chapters.

![Figure 2.1: ProtOS architecture overview. (A) ProtosPath and the Entity/Registry/Dataset system provide persistent identity and organization. (B, D) Processors handle specific data modalities: GRN Processor for standardized positions, Sequence Processor for FASTA and alignments, Structure Processor for coordinates, Embedding Processor for pLM representations, Graph Processor for contact networks, and Property Processor for tabular data. (C) The Model Manager orchestrates external compute for structure prediction (Boltz), embeddings (pLM), and spectral prediction (LAMBDA).](data/thesis_overview_protos_v1.png)

# Unified Data Access

ProtOS provides a single interface for acquiring protein data from multiple sources. A request for "1C3W" goes to the Protein Data Bank. A request for "P02945" goes to UniProt. The system recognizes identifier formats and routes each request to the correct database, downloads the data, and stores it locally. The same call works for any supported database; the user specifies what they want, not where to find it. Additional loaders extend coverage to other sources as needed.

```python
struct = StructureProcessor()
seq = SequenceProcessor()

struct.load("1C3W")           # Routes to RCSB PDB
seq.load("P02945")            # Routes to UniProt
```

Downloaded data is cached locally and persists across sessions. A structure fetched once remains available without re-downloading, and code that runs on a fresh installation produces identical results to code running against an existing cache. This caching transforms data assembly from a manual process into iteration. Collecting the 129 microbial rhodopsin structures for MOGRN required only a loop over PDB identifiers; collecting sequences for LAMBDA training required batch downloads with automatic registration. The hours spent navigating databases, clicking download buttons, and organizing files collapse into minutes of execution.

Every load operation does more than fetch data. It registers the result as an entity—ProtOS's abstraction for a biological object that persists across sessions and links related data regardless of source. This registration bridges data access and data management. A structure loaded from the PDB is not just a file on disk; it becomes an entity that can be cross-referenced with sequences from UniProt, annotations from the literature, and computed properties from subsequent analysis. The entity system addresses what unified access alone cannot: the same protein exists under different identifiers in different databases, and without persistent identity, the relationship between a structure, its sequence, and its measured properties would remain implicit, maintained only through file naming conventions and researcher memory.
# Entity, Registry, and Datasets

When is a protein the same protein? Bacteriorhodopsin exists under dozens of identifiers: UniProt accession P02945, PDB entries 1C3W, 1FBB, 2NTU and over fifty other structures, AlphaFold prediction AF-P02945-F1, gene names *bop* and *BR*, and designations like "bR," "BR," "Bacteriorhodopsin," and "BACR_HALSA." Each identifier is valid within its database, but none links the protein across systems.

This fragmentation reflects how databases evolved independently. The Protein Data Bank assigns identifiers to individual experiments, so the same protein crystallized under different conditions receives different codes. UniProt assigns accessions to track evolutionary relationships. Gene nomenclature follows its own logic. The result is a many-to-many mapping where one protein corresponds to multiple identifiers and the same string may reference different proteins in different databases.

Without explicit identity management, researchers maintain correspondence through folder hierarchies, file naming conventions, and memory. PDB structures provide coordinates but lack functional annotations. UniProt entries contain curated metadata but reference consensus sequences that may differ from experimental structures. Property measurements reside in spreadsheets disconnected from structural data. Every integrative analysis requires manually tracking which files correspond to which biological object.

The entity abstraction in ProtOS addresses this problem. An entity represents a single biological object—a protein, a ligand, a variant—independent of how it is stored or referenced. The bacteriorhodopsin entity connects to PDB structures, FASTA sequences, embeddings, residue annotations, and property measurements. Requests for structure, sequence, or embedding by name resolve to the same underlying identity.

The entity registry tracks these relationships. While the file system organizes data by type, the registry organizes it by identity, persisting as a JSON file at the data directory root. The format is human-readable, portable, compatible with version control, and editable with any text editor. Each entry contains a stable identifier, a canonical name, aliases for flexible lookup, and pointers to processor-specific files. Users interact with entities through names; requests for "BACR_HALSA," "P02945," "bacteriorhodopsin," or "1c3w" all resolve to the same entity through case-insensitive alias matching. When ambiguity arises, the registry raises an error rather than silently choosing.

Individual entities suffice for exploratory work, but scientific analysis rarely operates on isolated proteins. Comparative studies require defined sets: all GPCRs with known structures, all rhodopsins with measured absorption spectra, all kinases in a particular conformation. A dataset is a named, persistent collection of entities that survives across sessions and integrates with the processor ecosystem. A collection assembled with the Sequence Processor is available to the Embedding Processor for batch computation and to the Property Processor for annotation lookup. The relationships persist; the same entities receive embeddings, annotations, and properties as analysis proceeds.

With identity and organization addressed, the remaining question is where files live and how the system remains usable without configuration.
# Zero-Configuration Data Management

Systems requiring database installation, environment variables, and configuration files impose barriers that prevent adoption. A researcher under deadline will not pause to configure infrastructure; they will create another folder and another spreadsheet. For data management to succeed, setup cost must not exceed the cost of ad hoc alternatives.

ProtOS requires no configuration. A default directory in the user's home folder provides storage, requiring no special permissions and persisting across sessions. The first operation that requires storage creates this directory automatically; importing ProtOS without using it leaves no trace. For shared storage, containers, or custom paths, a single call or environment variable sets the data root. Directory structure follows processor organization, with each processor type maintaining its own subtree for raw inputs and computed outputs. The global registry at the root links files into coherent biological objects through the entity system described above.

Portability complements zero-configuration. The entity registry stores file locations as relative paths rather than absolute paths, so copying a ProtOS data directory preserves all internal relationships. A researcher can archive their data folder, transfer it to a collaborator, and the collaborator can use it from any location without rewriting paths or reconfiguring software. For reproducibility, a published analysis can include its ProtOS directory as supplementary material, allowing reviewers to access the full entity graph with all annotations. The MOGRN validation structures and LAMBDA training corpus exist as portable ProtOS datasets distributed alongside manuscript code.

With data access, identity, organization, and portability established, the framework is ready for the processors that transform raw data into analysis-ready representations.
# Processors

Every action in ProtOS flows through a processor. Loading a structure, computing an embedding, assigning residue numbers, storing a measured property—each operation belongs to a processor that handles one data modality. Processors integrate with the entity system so that outputs from one become inputs to another without explicit format conversion or identifier reconciliation.

Six processors handle the data types that opsin research requires. The Sequence Processor manages FASTA files, sequence searches, and alignments. The Structure Processor handles PDB and mmCIF coordinates, structural alignment, and binding site extraction. The Graph Processor builds residue contact networks for machine learning applications. The GRN Processor assigns Generic Residue Numbers—the standardized position annotation that enables cross-protein comparison. The Embedding Processor generates representations from protein language models. The Property Processor stores data associated with any of the other data formats in a tabular format, such as experimental measurements, computed annotations, or functional classifications.

These processors form a hierarchy (Figure 2.1). Sequence and Structure Processors are foundational; they load data from external sources. The GRN Processor requires sequence or structure data to annotate. The Embedding Processor requires sequences. The Graph Processor requires structures. The Property Processor links annotations to any entity. Outputs flow upward: a structure yields a sequence, which receives GRN annotation and embeddings, which enrich the nodes of a contact graph, which links to measured properties. This composition enables workflows that span modalities while the entity system maintains identity throughout.

The following sections present each processor through an example drawn from opsin research. The examples follow the progression that underlies LAMBDA: sequences reveal diversity, structures show how different folds bind the same chromophore, graphs capture binding pocket topology, GRN provides position identity, embeddings add features, and properties provide training targets.
# Sequence Processor

Sequence data is the most abundant modality in protein science. UniProt contains over 200 million protein sequences, growing as sequencing projects deposit new entries. Sequences are cheaper than structures, faster to analyze, and serve as input for nearly every computational method: homology search, alignment, phylogenetics, structure prediction, and machine learning.

For opsins, metagenomics transformed understanding of family diversity. What seemed a peculiarity of halophilic archaea and animal eyes turned out to be globally abundant. Proteorhodopsins, discovered from Pacific Ocean samples in 2000, are among the most common proteins on Earth. Metagenomics has since revealed tens of thousands of rhodopsin sequences across bacteria, archaea, algae, fungi, and viruses. Distant members share less than 20% sequence identity, yet all bind retinal and respond to light.

Sequence alone does not explain function. Cone opsins provide a clear example: short-wave and long-wave variants share over 80% sequence identity yet differ by approximately 140 nm in absorption maximum. Conversely, proteins with less than 20% identity may share similar spectral properties if their binding pocket residues are conserved. Sequence is the starting point, not the answer. Determining function requires examining structure, binding pocket composition, and chromophore-protein interactions. Every analysis begins with knowing what sequences exist.

Cone opsins illustrate this principle. These proteins mediate color vision in vertebrates, with short-wave (SW) opsins absorbing around 420 nm and long-wave (LW) opsins absorbing around 560 nm. The spectral classes correspond to blue and red sensitivity, enabling dichromatic or trichromatic vision depending on species. Before comparing structures or predicting spectral properties, the sequence diversity within this subfamily must be characterized: how many sequences exist, how they cluster, and whether sequence similarity correlates with spectral type.

Starting from human short-wave opsin OPN1SW (UniProt P03999) and long-wave opsin OPN1LW (UniProt P04000), BLAST searches against SwissProt retrieve homologs across vertebrates. The resulting dataset contains 200 cone opsins annotated by spectral type.

```python
loader = SequenceLoader(processor=seq_proc)
loader.download_and_register("uniprot:P03999", name="OPSB_HUMAN")

ncbi_loader = NCBILoader(processor=seq_proc)
result = ncbi_loader.blast_search(sequence, database="swissprot", hitlist_size=100)
ncbi_loader.download_batch(accessions, dataset_name="cone_opsin_diversity")

alignment_df = mmseqs2_align(sequences, sequences)
```

![Figure 2.2: Cone opsin sequence diversity. (a) Distribution by spectral type showing SW (n=100) and LW (n=23) opsins, with boxplots of within-type sequence identity. LW opsins show higher identity (median ~82%) compared to SW opsins (median ~47%), reflecting greater sequence diversity in the short-wave subfamily. (b) Hierarchical clustering dendrogram colored by spectral type, showing clear phylogenetic separation between SW and LW lineages.](data/protos/thesis/processors/figures/sequence_opsin_diversity.png)

![](data/protos/thesis/processors/figures/sequence_phylogenetic_tree.png)

The sequences cluster by spectral type: within-type identity exceeds 70% for LW opsins, while SW opsins show greater diversity with identity dropping to 40-50%. SW opsins (teal) and LW opsins (yellow) form distinct evolutionary lineages within the cone opsin subfamily.

Cone opsins represent a fraction of opsin diversity. Microbial rhodopsins alone contain over 40,000 sequences in public databases, with less than 20% identity between distant members. At this scale, systematic tools replace manual curation.

Sequences reveal diversity but not mechanism. SW and LW cone opsins differ by approximately 140 nm in λmax despite sharing over 80% sequence identity. The difference lies not in the sequence as a whole but in specific positions surrounding the chromophore. Identifying those positions requires examining structure.
# Structure Processor

Protein function emerges from three-dimensional structure. The spatial arrangement of amino acids determines how enzymes catalyze reactions, how receptors recognize ligands, and how channels gate ion flow. For opsins, structure explains how the same chromophore, retinal, drives different cellular responses depending on which protein binds it. Bacteriorhodopsin pumps protons across membranes. Channelrhodopsin opens a cation pore. Bovine rhodopsin activates G protein signaling. All use retinal isomerization as the trigger; the structural context determines the output.

The Protein Data Bank contains over 200,000 experimental structures from X-ray crystallography, cryo-EM, and NMR. AlphaFold adds 200 million predicted structures covering nearly all known sequences. This abundance transformed structural biology from data-scarce to data-rich, shifting the challenge from acquiring structures to managing them. Structures come in multiple formats from multiple sources with different metadata conventions.

Type I (microbial) and Type II (animal) opsins both bind retinal through a protonated Schiff base linkage to a conserved lysine, and both families contain seven transmembrane helices arranged around a central binding pocket. Yet these families evolved independently and share no detectable sequence homology outside the retinal-binding site. The seven-helix architecture arose twice through convergent evolution. The structural question is whether both families bind retinal the same way, or whether different folds arrive at similar function through different arrangements.

Bacteriorhodopsin (PDB 1C3W) is a proton pump from halophilic archaea and the founding member of the Type I family. It binds retinal covalently to lysine 216 in transmembrane helix 7. The Structure Processor extracts coordinates, identifies the retinal ligand, and computes binding pocket residues based on proximity.

```python
loader = StructureLoader(processor=struct_proc)
loader.download("1C3W", dataset_name="opsin_structures")

interactions = struct_proc.get_ligand_interactions("1C3W", ligand_id="RET", cutoff=7.0)
pocket_residues = struct_proc.get_binding_pocket("1C3W", ligand_id="RET", cutoff=4.0)
```

The binding pocket contains 19 residues within 4Å of retinal, spanning helices 3 through 7. These residues determine the protein's spectral and functional properties—the counterion at position 85, the proton donor at position 96, and the spectral tuning sites distributed across the pocket.

![Figure 2.3: Bacteriorhodopsin binding pocket. The retinal chromophore (yellow) sits in a pocket formed by residues from TM3-TM7. Key functional positions are labeled: D85 (counterion), D96 (proton donor), and K216 (Schiff base). The Structure Processor extracts these interactions automatically from any PDB structure containing a bound ligand.](data/protos/thesis/processors/figures/br_pocket.png)

Type I (microbial) and Type II (animal) opsins both bind retinal through a protonated Schiff base, yet these families evolved independently and share no detectable sequence homology. The seven-helix architecture arose twice through convergent evolution.

The binding poses are not homologous. Type I and Type II opsins have different helix topologies. The helices are arranged differently, the loops follow different paths, and retinal sits in geometrically distinct pockets. Global structural alignment fails because the folds are not superimposable.

Yet both families tune spectral properties through the same physical mechanism—the electrostatic environment around the chromophore. The protein environment shifts electron density along retinal's conjugated chain, changing the energy gap for light absorption. Aligning on the chromophore succeeds where global alignment fails. The retinal position is conserved because it must be; everything else can vary.

This observation enables cross-family prediction. The folds differ but the functional principle is shared. If binding pockets can be represented comparably despite different global structures, cross-family learning becomes possible.

Structural comparison raises a practical problem. The conserved counterion is position 85 in bacteriorhodopsin and position 113 in bovine rhodopsin. These residues serve the same function—stabilizing the protonated Schiff base—but their sequence positions are arbitrary. Systematically identifying equivalent positions across structures requires representing which residues contact which, independent of their sequence numbering. A graph provides this representation: residues as nodes, spatial contacts as edges.
# Graph Processor

Proteins are networks of interacting residues, not linear sequences. Amino acids distant in sequence may be adjacent in space, forming contacts that determine stability, dynamics, and function. A contact graph represents this spatial organization—nodes correspond to residues, edges connect pairs within a distance threshold. The graph topology encodes which residues interact, information sequence alone cannot provide.

Graph representations have become central to machine learning on proteins. Graph neural networks operate directly on molecular structure, propagating information along contact edges to learn spatial patterns. For property prediction, the graph provides a natural substrate—nodes carry features describing each residue, edges encode three-dimensional relationships.

For opsins, the binding pocket graph is the relevant representation. Spectral properties depend on residues surrounding the chromophore, not distant loops or terminal regions. Extracting the binding pocket as a graph focuses on the functionally relevant subset.

Structural alignment showed that bacteriorhodopsin (Type I) and bovine rhodopsin (Type II) have different global folds but share the same chromophore position. The question is whether this similarity extends to binding pocket topology—whether both pockets can be represented as comparable graphs despite the fold difference.

Binding pocket graphs are extracted from both aligned structures. The pocket is defined as residues within 7 Å of retinal, capturing the spectral tuning sites. Edges connect residue pairs whose Cα atoms are within 6 Å, a standard threshold for residue contacts.

```python
graph_proc = GraphProcessor(default_cutoff=6.0, default_level="residue")

graph_bR = graph_proc.generate_graph(
    structure_id="1C3W",
    near_hetatm=("RET", 7.0),
    edge_cutoff=6.0
)

graph_bovine = graph_proc.generate_graph(
    structure_id="1U19",
    near_hetatm=("RET", 7.0),
    edge_cutoff=6.0
)
```

![Figure 2.4: Binding pocket graphs for Type I and Type II opsins. (a) Bacteriorhodopsin (1C3W) binding pocket graph with protein residues (blue) and retinal (orange) shown in spring layout. (b) Bovine rhodopsin (1U19) binding pocket graph with protein residues (red) and retinal (orange). (c) Quantitative comparison showing similar topology: 1C3W contains 46 residues and 162 contacts (15.7% density), while 1U19 contains 44 residues and 146 contacts (15.4% density).](data/protos/thesis/processors/figures/1c3w_graph_pocket.png)

![](data/protos/thesis/processors/figures/1u19_graph_pocket.png)

![](data/protos/thesis/processors/figures/graph_comparison.png)

The bacteriorhodopsin pocket contains 46 residues forming 162 edges at the 6 Å threshold; the bovine rhodopsin pocket contains 44 residues forming 146 edges. Node counts and edge densities fall in the same range, demonstrating that the pockets have similar topological complexity despite the fold difference.

This comparability enables cross-family learning. Both Type I and Type II opsins form contact networks of similar density around retinal. The graph representation captures which residues matter and how they relate spatially, abstracting away the global architecture that differs between families.

Two pieces remain missing from a complete representation. The first is node identity: residue 85 in bacteriorhodopsin has no defined equivalent in bovine rhodopsin, though both residues serve as the counterion stabilizing the protonated Schiff base. Sequence position numbers are arbitrary and do not transfer across proteins. The second is node features: graph topology alone encodes spatial relationships but not amino acid identity, conservation, or evolutionary context. A node knows its neighbors but not what amino acid it represents. Generic residue numbering addresses the first problem by providing standardized position identifiers within protein families; protein language model embeddings address the second by providing per-residue feature vectors.
# GRN Processor

The residue numbering problem underlies every cross-protein comparison. Position 85 in bacteriorhodopsin and position 113 in bovine rhodopsin are structurally equivalent—both are the counterion stabilizing the protonated Schiff base. Yet their sequence positions differ because insertions and deletions accumulate differently in different lineages. Sequence position is a poor indicator of structural or functional equivalence.

Generic Residue Numbering (GRN) solves this problem through two steps. First, structural analysis of representative family members establishes which positions are equivalent—the same location in three-dimensional space across different proteins. This structural work produces a reference table: a curated alignment where each column represents a structurally equivalent position, and each row represents a family member with known coordinates. Second, new sequences receive GRN annotation by aligning to this reference table. The alignment transfers coordinates from structurally characterized references to the query sequence. The structural analysis is done once; subsequent annotation requires only sequence alignment.

The standard convention for seven-transmembrane proteins numbers each helix 1 through 7, with the anchor residue on each helix designated X.50. Position 3.50 is the anchor on helix 3; positions 3.49 and 3.51 are its immediate neighbors. Under this system, position 3.50 refers to the same structural location in any protein using that GRN system, regardless of sequence length or insertion/deletion history.

For GPCRs, Ballesteros and Weinstein introduced this numbering in 1995, and it became the field standard. GPCRdb maintains reference tables covering the entire GPCR superfamily—thousands of sequences with structurally validated annotations. The GRNProcessor in ProtOS aligns query sequences to these reference tables and transfers the GRN coordinates. Findings reported in GRN positions are immediately interpretable across laboratories.

The binding pocket graphs from Type I and Type II opsins have comparable topology, but graph nodes carry arbitrary sequence position labels that do not transfer across proteins. Comparing binding pocket positions systematically—the counterion and spectral tuning sites that determine absorption—requires standardized numbering.

Animal opsins (Type II) are GPCRs, so Ballesteros-Weinstein numbering applies. Bovine rhodopsin (PDB 1U19), a vertebrate visual opsin, and squid rhodopsin (PDB 4ZWJ), an invertebrate visual opsin, both have Ballesteros-Weinstein annotation through GPCRdb. The test is whether GRN identifies functionally equivalent positions across these evolutionarily distant opsins.

```python
grn_table, summary = seq_proc.annotate_with_grn(
    dataset_name="opsin_comparison",
    reference_table="gpcrdb_ref",
    protein_family="gpcr_a",
    output_table="opsin_comparison_grn"
)

counterion_bovine = grn_table.loc["bovine_1U19", "3.28"]  # E113
counterion_squid = grn_table.loc["squid_4ZWJ", "3.28"]    # E180

for grn_pos in ["3.28", "3.32", "5.46", "6.48", "7.43"]:
    bovine_res = grn_table.loc["bovine_1U19", grn_pos]
    squid_res = grn_table.loc["squid_4ZWJ", grn_pos]
```

![Figure 2.5: GRN annotation enables cross-species comparison. Bovine rhodopsin (pink, PDB 1U19) and squid rhodopsin (yellow, PDB 4ZWJ) superimposed, showing how Ballesteros-Weinstein numbering identifies equivalent positions despite 67-residue offset in sequence numbering. The counterion (E113 in bovine, E180 in squid) maps to position 3.28 in both structures. Retinal and binding pocket residues occupy equivalent GRN positions across these evolutionarily distant visual opsins.](data/protos/thesis/processors/figures/bovine_squid_aligned.png)

Functionally equivalent positions align despite different sequence numbering. The counterion is E113 in bovine rhodopsin and E180 in squid rhodopsin—67 sequence positions apart, yet both map to GRN position 3.28. The spectral tuning site at 3.32 is A117 in bovine and A184 in squid. The rotamer switch at 6.48 is W265 in bovine and W332 in squid. The retinal-binding lysine at 7.43 is K296 in bovine and K363 in squid. Standardized numbering makes functional equivalence explicit without manual alignment for every comparison.

Type II opsins benefit from Ballesteros-Weinstein numbering because they are GPCRs. Bovine rhodopsin was one of the founding structures for this system, and the animal opsin subfamily inherits the annotation infrastructure that GPCRs developed over three decades.

Type I opsins (microbial rhodopsins) are not GPCRs. They have a different fold and different evolutionary origin, and Ballesteros-Weinstein does not apply because the helix topology differs. No equivalent standardized numbering system exists for microbial rhodopsins. Comparing position 85 in bacteriorhodopsin to position 156 in channelrhodopsin requires manual structural alignment, and the literature uses inconsistent conventions that impede systematic comparison.

MOGRN fills this gap. It provides for microbial rhodopsins what Ballesteros-Weinstein provides for GPCRs: a shared coordinate system anchored to the retinal-binding site, enabling systematic comparison across the superfamily.
# Embedding Processor

Protein language models have transformed how researchers extract information from sequences. Models like ESM-2 and Ankh learn evolutionary constraints from sequence patterns, producing embeddings where similar functions cluster together even when sequences have diverged beyond recognizable similarity. Representations trained on millions of sequences capture patterns that transfer to specific prediction tasks.

Embeddings are dense vectors, typically hundreds to thousands of dimensions, representing proteins in a learned space. They can be computed at the sequence level, producing a single vector summarizing an entire protein, or per-residue, producing one vector for each amino acid position. Sequence-level embeddings support clustering and classification tasks. Per-residue embeddings support position-specific predictions, including binding site identification and property mapping.

The binding pocket graphs now have standardized position identifiers through GRN annotation, but graph nodes lack features describing the amino acids they represent. Protein language model embeddings provide these features, encoding amino acid identity, sequence context, and evolutionary patterns.

The cone opsin diversity dataset assembled in the Sequence Processor section offers a test of whether embeddings capture functional relationships. Those 200 sequences are annotated by spectral type: short-wave (SW) opsins absorbing around 420 nm and long-wave (LW) opsins around 560 nm. If embeddings encode functional similarity, opsins of the same spectral type should cluster together even when sequence identity is moderate. Sequence-level embeddings from ESM2-t12-35M allow this comparison.

```python
emb_proc = EmbeddingProcessor(model_name="esm2_t12_35m")

embeddings = emb_proc.embed_sequences(
    sequences,
    embedding_type="mean",
    save_dataset="cone_opsin_diversity__esm2_t12_35m__mean"
)

emb_norm = F.normalize(emb_matrix, p=2, dim=1)
sim_matrix = torch.mm(emb_norm, emb_norm.T)
```

![Figure 2.6: Embeddings capture functional relationships. (a) t-SNE projection of ESM2 embeddings for 123 cone opsins, colored by spectral type. SW opsins (teal) and LW opsins (yellow) form distinct clusters despite moderate sequence identity. (b) Pairwise cosine similarity matrix sorted by spectral type, showing block structure with higher within-type similarity (darker blue) and lower between-type similarity.](data/protos/thesis/processors/figures/embedding_tsne.png)

![](data/protos/thesis/processors/figures/embedding_heatmap.png)

The protein language model learned something about spectral function from sequence alone, without access to absorption measurements during training. SW and LW opsins cluster separately in embedding space, indicating that sequence patterns associated with spectral class are captured by the representation.

Clustering by spectral type is not the same as predicting exact absorption wavelength. SW opsins absorb between 360 and 440 nm; LW opsins between 555 and 565 nm. The embeddings separate broad classes but do not predict where within each range a specific opsin absorbs. A SW opsin at 360 nm and one at 440 nm may have similar embeddings despite an 80 nm difference. Embeddings are features, not predictors. They encode sequence context and evolutionary patterns but do not output functional properties directly. For LAMBDA, per-residue embeddings serve as node features in the binding pocket graph, describing each position's amino acid and sequence context. The graph structure provides topology, the GRN provides position identity, and measured spectra provide training signal; together these enable learning the mapping from binding pocket composition to spectral properties.
# Property Processor

Beyond sequences, structures, and embeddings, protein research requires measurements and annotations that characterize what proteins do. Experimental properties include absorption maxima, binding affinities, kinetic rates, and thermostability measurements. Computed properties include predicted values from models, conservation scores from alignments, and functional annotations from databases. These properties exist independently of any particular sequence or structure file, yet they must link to the entities they describe.

The representational foundation—sequences, structures, binding pocket graphs, standardized positions, and embeddings—is now in place. Training a spectral prediction model requires ground truth. Measured absorption spectra from laboratory experiments provide target values for supervised learning, and these measurements must associate with registered sequence entities.

The challenge is coverage. Returning to the cone opsin diversity dataset: those 200 sequences now have embeddings, but only 9 have published absorption maxima. This sparse coverage is typical: most opsins have been sequenced but not characterized spectroscopically. The ratio of sequences to measurements defines both the problem LAMBDA addresses and the opportunity it creates. Linking the available literature measurements to sequence entities establishes ground truth for the subset where it exists.

```python
prop_proc = PropertyProcessor()

spectral_data = {
    "entity_id": ["P03999", "P04000", "P51491", ...],
    "lambda_max": [420, 560, 360, ...],
    "opsin_type": ["short_wave", "long_wave", "short_wave", ...],
    "species": ["Human", "Human", "Mouse", ...]
}
prop_proc.create_table("cone_opsin_spectra", spectral_data)

blue_opsins = prop_proc.query_table(
    "cone_opsin_spectra",
    filters={"lambda_max": {"<": 450}}
)
```

![Figure 2.7: Spectral sensitivity curves for characterized opsins. Gaussian approximations centered at measured λmax values, overlaid on the visible spectrum. The curves span from UV-sensitive opsins (~350 nm) through blue (~420 nm), green (~500 nm), and red-shifted variants (~560-630 nm). Each curve represents an opsin with experimentally measured absorption maximum stored as a property in ProtOS.](data/protos/thesis/processors/figures/property_sensitivity.png)

The characterized opsins span the visible spectrum from approximately 350 nm to 630 nm.

Of the 200 sequences in the cone opsin dataset, 9 have literature absorption values. This 4.5% coverage rate reflects the broader challenge—sequences accumulate faster than functional characterization. For microbial rhodopsins, the disparity is larger still.

Across all opsins, the literature contains approximately 1,800 measured absorption maxima assembled from decades of spectroscopy across hundreds of publications. These measurements are scattered across different formats, different databases, and different naming conventions. Unified access requires systematic curation that links measurements to sequence identifiers.

This curated dataset becomes the training data for LAMBDA. Each measurement provides a target value; each sequence provides input features through the processing pipeline described above. When LAMBDA predicts absorption maxima for the roughly 40,000 opsins in public databases, those predictions are stored as properties. Users can query for opsins with predicted absorption above 600 nm, identifying candidates for experimental characterization or engineering. Predictions become queryable alongside experimental measurements, integrated with sequence, structure, and embedding data.

The first six processors establish the components for LAMBDA: sequences of proteins that bind retinal, structural comparison revealing different folds with the same ligand, binding pocket graphs capturing local topology, GRN annotation providing position identity, embeddings encoding sequence context, and property measurements providing training targets.
# Model Manager

Processors transform data within ProtOS. Some analyses require external models that run on GPU clusters, consume specific input formats, and produce outputs that must integrate back into the entity system. The Model Manager handles this orchestration, packaging processor outputs into valid job submissions and registering returned predictions as entity properties.

Structure prediction illustrates the need. AlphaFold and Boltz2 predict three-dimensional structure from sequence, but running these models at scale requires managing input formats, GPU resources, job queues, and output organization.

```python
from protos.models import ModelManager

models = ModelManager()

job = models.prepare_input(
    model="boltz2",
    sequence=sequence,
    output_name="rhodozyme_TRP_01"
)

manifest = models.prepare_batch(
    model="boltz2",
    sequences=design_sequences,
    output_dir="rhodozyme_predictions"
)
```

The Model Manager supports models across several categories: structure prediction (Boltz2), protein design (RFdiffusion2, Boltzgen, LigandMPNN), sequence embeddings (ESM-2, Ankh), and spectral prediction (LAMBDA). Each model has a defined input schema; the Model Manager validates inputs and constructs correctly formatted job directories. When predictions complete, they register back into ProtOS as entity properties, queryable alongside the experimental data that processors assembled.

For this thesis, the Model Manager orchestrates LAMBDA predictions for spectral properties, Boltz2 predictions for designed protein structures, and RFdiffusion2 for theozyme scaffolding in the rhodozyme design workflow. Data flows through processors, reaches the Model Manager for external inference, and returns as queryable predictions integrated with the entity system.
# Limitations and the Accessibility Gap

ProtOS solves the data management problem for protein research. It does not solve the adoption problem.

The framework has technical boundaries. Database coverage includes UniProt, PDB, AlphaFold DB, and NCBI, but not every biological database; extending to additional sources requires implementing new loaders. GRN annotation depends on having a numbering system defined for the protein family of interest—GPCRs have Ballesteros-Weinstein, microbial rhodopsins now have MOGRN, but families without standardized numbering cannot use the GRN Processor until such systems are developed.

More fundamental are the barriers to use. ProtOS requires Python. This excludes the experimental scientists who generate the data that computational tools analyze—crystallographers, spectroscopists, biochemists who understand protein function deeply but have never written a for-loop. The barrier is not intellectual; it is technical. An experimentalist who asks "which residues contact retinal in this structure?" understands exactly what they want. They cannot get it without a programmer.

Bioinformaticians face a different barrier. They can write Python. But adopting a new framework requires learning its conventions, understanding its abstractions, restructuring existing workflows. For a single protein analysis, this overhead exceeds the benefit. If you need to check binding pocket residues in one structure, a quick script with BioPython suffices. ProtOS pays off when managing hundreds of structures, thousands of sequences, systematic comparisons across protein families—exactly the scale where ad hoc scripts become unmaintainable. But that scale is not every project, and for smaller tasks, bioinformaticians reasonably judge that learning a framework costs more than it saves.

These are the barriers that hinder adoption of any computational framework. Experimental scientists cannot use tools that require programming. Bioinformaticians will not invest in learning new abstractions unless the payoff is clear.

ProtOS-MCP addresses these barriers. Rather than requiring Python syntax, a researcher describes what they want in natural language. The system translates intent into tool calls and returns results. The experimentalist who asks about retinal contacts receives an answer without writing code. The bioinformatician who needs a quick analysis can ask rather than script. The framework's capabilities become accessible without the framework's learning curve. Chapter 6 describes this interface and evaluates its effectiveness.

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
# LAMBDA: Cross-Family Spectral Prediction from Sequence

The preceding chapter described MOGRN, a generic residue numbering system for microbial rhodopsins. This chapter presents LAMBDA (Light Absorption Modeling through Binding Domain Analysis), which uses GRN coordinates to predict opsin spectral properties from sequence. Existing methods treat spectral prediction as a family-specific regression problem. LAMBDA instead predicts how a binding pocket tunes retinal absorption across chromophore conformations. This reframing reflects the underlying biology—the same protein can bind retinal in multiple states, each with distinct spectral properties determined by the same binding pocket environment—and motivates a structure-based approach that extends to any retinal-binding protein.
Three capabilities follow from this formulation. LAMBDA predicts λmax for multiple chromophore states (11-cis, all-trans, protonated, deprotonated) from one model. Type I and type II opsins train jointly because their binding pockets align through the chromophore despite different global folds, testing whether spectral tuning principles transfer across families. And the approach extends to new folds given appropriate GRN mappings and training examples, as demonstrated with hCRBPII. These capabilities serve three applications: mapping the spectral landscape of all opsins (the Opsin Atlas), addressing spectral tuning as a protein engineering challenge, and understanding how retinal-binding pockets might be incorporated into synthetic folds to couple light activation to arbitrary protein functions. The chapter demonstrates that standardized positional annotation—MOGRN for type I opsins, Ballesteros-Weinstein for type II—enables cross-family learning about chromophore-protein interactions, providing a foundation for understanding spectral tuning across retinal-binding proteins.

# Introduction

Optogenetics depends on matching light wavelength to protein response. Type I (microbial) opsins—channelrhodopsins, bacteriorhodopsins, halorhodopsins—provide the ion channels and pumps that enable neural activation and silencing. Type II (animal) opsins—rhodopsins, melanopsins, cone opsins—underlie visual neuroscience and increasingly serve as light-activated GPCRs for controlling intracellular signaling. In both cases, the absorption maximum λmax determines which wavelengths activate the protein, and selecting or engineering opsins with specific λmax values is a core challenge. Red-shifted opsins (λmax > 580 nm) enable deeper tissue penetration because longer wavelengths scatter less in biological tissue. Spectrally separated pairs—one blue-absorbing, one red-absorbing—allow multiplexed control of distinct cell populations in the same preparation. Many optogenetic applications also require a distinct off-switch: bistable opsins can be both photoactivated and photoreversed, but this bidirectional control requires sufficient separation between dark-state and activated-state absorption maxima. If the two states absorb at similar wavelengths, a single light source drives both forward and reverse photoreactions, precluding selective control. Predicting λmax for both chromophore states of the same protein, and the spectral gap between them, is therefore as important as predicting a single absorption maximum.

Beyond opsins, the landscape of retinal-binding proteins is expanding. Engineered systems like human cellular retinol binding protein II (hCRBPII) demonstrate that retinal can be bound and spectrally tuned in protein folds entirely unrelated to opsins. As protein engineering capabilities grow—through directed evolution, computational design, and de novo fold generation—the prospect of incorporating retinal-binding pockets into novel protein scaffolds creates the possibility of coupling light activation to functions beyond the ion transport, channel gating, and GPCR signaling that natural opsins provide. Understanding how binding pocket architecture determines spectral properties, independent of global fold, is a prerequisite for this direction.

Light absorption follows the relationship E = hc/λ: a photon is absorbed when its energy matches the gap between electronic states of the chromophore. The protein environment tunes this gap through electrostatic interactions with the conjugated polyene chain of retinal. The counterion—a negatively charged residue 3–4 Å from the protonated Schiff base—is the primary determinant, but the entire network of binding pocket residues contributes, and single mutations at key positions shift λmax by 10–40 nm. All retinal-binding proteins use the same chromophore yet span λmax values from ~350 nm (UV) to ~650 nm (far-red), with the variation arising along two independent axes: geometric isomerism and protonation state. Type II opsins bind 11-cis retinal in the dark state; type I opsins and photoactivated type II opsins bind all-trans retinal. Independently, the Schiff base linkage to the conserved lysine can be protonated (the dominant state for most opsins, absorbing in the visible range) or deprotonated (as in UV-sensitive SWS1 opsins, which bind 11-cis retinal but absorb below 400 nm). Both properties—isomer state and protonation—are determined by the same binding pocket environment (Figure 4.1). Because energy is inversely proportional to wavelength, spectral tuning becomes increasingly difficult to resolve at long wavelengths: a 10 nm shift near 400 nm represents ~0.08 eV, while the same shift near 600 nm represents only ~0.03 eV.

[FIGURE 4.1: Spectral tuning in opsins. (A) The relationship between energy gap and absorption wavelength (E = hc/λ). Equal wavelength differences correspond to smaller energy differences at longer wavelengths—a 10 nm shift near 400 nm represents ~0.08 eV, while the same shift near 600 nm represents only ~0.03 eV. (B) Retinal Schiff base configurations: 11-cis protonated (dark state in type II opsins), all-trans protonated (type I opsins and activated type II), and deprotonated (UV-sensitive opsins).]

Tens of thousands of opsin sequences exist in databases, but fewer than 3,000 have measured λmax values. Several computational methods predict spectral properties: OPTICS predicts λmax for type II opsins using sequence features and phylogenetic information trained on the VPOD dataset; Inoue et al. trained a regression model for type I opsins from amino acid composition at key positions; RhoMax uses a graph neural network on type I opsin binding pocket structures, requiring experimental or predicted structures as input. Each achieves reasonable accuracy within its domain, but they share two fundamental limitations. They are fold-specific—OPTICS cannot predict type I opsins, the Inoue model cannot predict type II opsins, none can handle retinal-binding proteins with different folds. And they predict only a single λmax for their target family's native chromophore state, unable to predict activated states, deprotonated states, or the spectral gap between chromophore states that determines whether bidirectional control is feasible. At a deeper level, pure sequence-based methods learn statistical correlations—a charged residue at a particular sequence position correlates with a blue shift—which captures patterns within a family but does not generalize, because the same spectral effect arises from different sequence positions in different folds.

LAMBDA addresses these limitations by framing spectral prediction as a question about binding pockets rather than sequences. The model predicts λmax for multiple retinal conformations and protonation states simultaneously from a single shared representation. Predicting both dark-state (11-cis) and activated-state (all-trans) absorption maxima for the same protein yields the spectral shift Δλ between on and off states—a quantity that determines whether bidirectional optogenetic control is feasible and that no existing method can estimate.

Since spectral tuning is determined by the binding pocket, the natural representation for prediction is the pocket itself—the residues surrounding retinal and their spatial contacts—rather than the full protein sequence or the global fold. In LAMBDA, binding pockets are represented as graphs: nodes are binding pocket residues, edges encode spatial contacts between them. This topology captures which residues interact, independent of their sequence positions or the global architecture that places them there. Figure 4.2 illustrates this across three structurally unrelated protein families. Type I opsins (bacteriorhodopsin, 1C3W) and type II opsins (bovine rhodopsin, 1U19) are both seven-transmembrane alpha-helical bundles, but with independently evolved architectures—different helix tilts, loop connectivities, and retinal binding orientation—and no detectable sequence homology. hCRBPII (4QYP) adopts an entirely different fold—a beta-barrel—yet engineered variants bind retinal and span over 200 nm of spectral tuning through only nine point mutations. Despite these architectural differences, the graph representation reduces each binding pocket to the same format: nodes are residues within contact distance of retinal, edges encode spatial adjacency within 4 Å. The opsin graphs have similar sizes (68 and 69 nodes); the hCRBPII graph has comparable complexity (67 nodes) despite the entirely different fold. The same model architecture processes all three graph types, learning spectral tuning determinants from binding pocket composition regardless of the surrounding protein scaffold.

[FIGURE 4.2: Binding pocket structures and graph representations for three retinal-binding protein families. Top: crystal structures with pocket residues shown as sticks. Bottom: corresponding contact graphs used as LAMBDA input. (A, D) Bacteriorhodopsin, Type I opsin. (B, E) Bovine rhodopsin, Type II opsin. (C, F) hCRBPII, a beta-barrel protein. Opsin residues colored by transmembrane helix (TM1–TM7).]

A binding pocket graph derived from a crystal structure is the definitive representation, but structures exist for only a fraction of known opsins. Two components make the approach scalable. Generic Residue Numbering (GRN) systems provide standardized positional coordinates for transmembrane proteins: the Ballesteros-Weinstein system for type II opsins and MOGRN—described in the preceding chapter—for type I opsins. GRN positions are defined relative to conserved structural features, so a sequence annotated with GRN implicitly specifies which binding pocket positions are present. Because the seven-transmembrane fold is conserved within each opsin family, the contact pattern between GRN positions in the binding pocket is highly conserved—the same core edges connect the same positions across opsins within each family, though flexible loop regions may contribute variable contacts. Binding pocket graphs can therefore be constructed from sequence alone, using the reference topology as a template. Aligning reference structures from both families on retinal reveals that certain positions occupy equivalent spatial locations relative to the chromophore: position 7.43 (type II) and 7.50 (type I) are both the Schiff base lysine; 3.28 and 3.57 are both in the counterion region. Treating mapped positions as equivalent graph nodes enables cross-family learning despite the independent evolutionary origins of these folds. Protein language model embeddings from Ankh-large (Elnaggar et al. 2023), chosen for its balance of embedding quality and computational efficiency, provide 1536-dimensional per-residue feature vectors that serve as node features, encoding amino acid identity and evolutionary context. Both GRN annotation and pLM embeddings require only a protein sequence as input.

Two pathways construct LAMBDA input graphs, depending on whether a structure is available (Figure 4.3). The structure-based path processes any retinal-binding protein with a known or predicted structure: binding pocket residues are identified based on proximity to retinal, GRN coordinates are assigned, per-residue pLM embeddings are computed, and the contact graph is assembled. The sequence-only path annotates the sequence with GRN positions, applies the family-specific binding pocket topology directly, and uses pLM embeddings as node features—this is the path used for the Opsin Atlas predictions on 47,700 sequences. Both paths produce the same output: a binding pocket graph where each node carries a positional encoding and a pLM embedding, and edges encode spatial contacts. A graph attention network learns the mapping from binding pocket composition to spectral properties, trained on 2,120 proteins with measured absorption values across three protein folds.

![Figure 4.3: Data flow for the two binding pocket graph construction pathways. The structure-based path (top) processes any retinal-binding protein through Structure Processor → GRN annotation → Graph assembly. The sequence-only path (bottom) exploits conserved GRN systems to construct graphs directly from sequence via GRN annotation → graph topology lookup → embedding enrichment. Both paths produce the same graph format consumed by LAMBDA.](data/thesis_overview_protos_v1.png)

# Methods

## Datasets

Training data were compiled from three sources that together span both opsin families and multiple retinal conformations.

VPOD 1.3 (Visual Physiology Opsin Database) contains 1,253 type II opsin sequences spanning vertebrate and invertebrate species. The dataset includes both wild-type and experimentally characterized mutant opsins, with λmax values for 11-cis retinal ranging from 350–611 nm. VPOD provides the largest collection of animal opsin spectral data, spanning vertebrate visual opsins, invertebrate opsins, and non-visual opsins such as melanopsin.

The Inoue dataset expands coverage to type I (microbial) opsins, containing 785 sequences from 30 species. These sequences were collected from NCBI, metagenomic databases, and the Tara Oceans project, with experimentally determined λmax values for all-trans retinal. The dataset spans bacteriorhodopsins, halorhodopsins, channelrhodopsins, and sensory rhodopsins, with absorption maxima ranging from 436–644 nm.

To enable prediction of the all-trans retinal state for type II opsins, training data were augmented with activated-state λmax values collected from the literature. For 31 type II opsins with reported photoactivated or meta-state absorption maxima, the all-trans λmax was recorded alongside the dark-state (11-cis) value. This approach is necessarily approximate: "activated state" encompasses different photointermediates depending on the opsin's photocycle (meta II in visual opsins, the stable photoproduct in bistable opsins), and the reported values may reflect different experimental conditions. Nevertheless, these dual-state annotations provide the only available ground truth for training the all-trans prediction head on type II opsins. These 31 opsins were assigned exclusively to the training split to maximize the model's exposure to dual-state measurements.

hCRBPII (human cellular retinol binding protein II) represents an engineered system where Wang et al. (2012) created retinal-binding variants through targeted mutations. This dataset of 51 sequences demonstrates spectral tuning spanning 425–644 nm with all-trans retinal—a range of 219 nm achieved through only nine point mutations. hCRBPII adopts a lipocalin fold unrelated to opsins, providing a test case for cross-fold generalization. As with type I and type II opsins, the hCRBPII binding pocket was aligned to retinal and equivalent positions mapped to the nearest GRN coordinates, enabling the same graph-based representation.

| Dataset | Samples | Fold | Retinal | λmax range (nm) |
|---------|---------|------|---------|-----------------|
| VPOD 1.3 | 1,253 | Type II | 11-cis | 350–611 |
| Inoue | 785 | Type I | all-trans | 436–644 |
| Dual-state (this work) | 31 | Type II | both | varies |
| hCRBPII | 51 | Lipocalin | all-trans | 425–644 |
| **Total** | **2,120** | | | |

The distribution of measured λmax values across these datasets is shown in Figure 4.4. The diversity of these datasets—spanning different opsin types, species, and retinal conformations—necessitates a multi-output approach. Type II opsins bind 11-cis retinal in the dark state, while type I opsins bind all-trans retinal. The model addresses this by predicting λmax for multiple retinal states simultaneously, allowing it to learn shared spectral tuning principles while capturing conformation-specific effects. This creates an inherent asymmetry: the all-trans prediction head is trained on type I opsins (Inoue, n=785), hCRBPII (n=51), and the 31 dual-state type II opsins, while the 11-cis head is trained exclusively on type II data (VPOD, n=1,253). The model does produce 11-cis predictions for type I opsins, but these are extrapolations with no training signal—type I opsins natively bind all-trans retinal, and predicting their absorption with 11-cis or other isomers (9-cis, 13-cis) would require isomer-specific training data that does not currently exist at scale.

## Binding Pocket Graphs

Type I (microbial) and type II (animal) opsins share no detectable sequence homology, yet both evolved seven-transmembrane architectures that bind retinal. LAMBDA exploits this convergence by representing binding pockets as graphs aligned on the chromophore itself.

A single reference structure defines the binding pocket graph for each protein family: bovine rhodopsin (PDB: 1U19) for type II opsins and bacteriorhodopsin (PDB: 1C3W) for type I opsins.

Within each reference structure, binding pocket residues are identified as those with any sidechain atom within 7 Å of retinal. These residues become nodes in the graph. Edges connect residue pairs with atoms within 4 Å of each other, encoding the network of spatial contacts through the binding pocket.

To enable cross-family learning, the two reference structures are aligned on the retinal polyene chain (atoms C7–C15) using Kabsch superposition. After alignment, binding pocket residues from both families that occupy equivalent spatial locations are identified through sidechain volume overlap: for each pair of GRN positions (one from each family) whose Cα atoms fall within 5 Å, a voxelized sidechain overlap is computed, and pairs exceeding 15% overlap of the smaller sidechain volume are mapped in a greedy one-to-one matching. This procedure yields 25 cross-family position pairs out of ~70 binding pocket positions per family. Positions without a cross-family partner are retained as family-specific nodes. The Schiff base lysine (7x43 in type II, 7x50 in type I) is hardcoded as a mapped pair.

| Structural feature | Type II (animal) | Type I (microbial) |
|-------------------|------------------|-------------------|
| Schiff base Lys | 7x43 | 7x50 |
| Counterion region | 3x28 | 3x57 |
| Retinal contact 1 | 7x39 | 7x46 |
| Retinal contact 2 | ECL2 | 4x51 |
| Retinal contact 3 | 2x58 | 7x45 |

Mapped positions share graph node indices during training, enabling the model to learn that residues at equivalent spatial locations serve analogous roles in spectral tuning despite belonging to different protein folds. Unmapped positions contribute family-specific nodes that capture tuning determinants unique to each fold.

## Preprocessing

Generating binding pocket graphs for any opsin sequence without requiring experimental structures is achieved through Generic Residue Numbering (GRN) systems, which assign structural coordinates based on sequence alignment alone.

GRN systems label positions within transmembrane proteins using a standardized notation: the helix number followed by a position relative to a conserved reference residue (e.g., 3x50 denotes position 50 in helix 3). For type II opsins, the established Ballesteros-Weinstein numbering system is used. For type I opsins, MOGRN is used—the GRN system for microbial rhodopsins described in the preceding chapter, with reference positions (Xx50) chosen as the residue closest to retinal in each helix.

The preprocessing pipeline was implemented in ProtOS, providing standardized procedures for assigning GRN coordinates to opsin sequences. The pipeline operates as follows:

1. **Family classification**: Input sequences are classified as type I or type II based on sequence similarity to reference opsins using MMseqs2.
2. **Reference selection**: The most similar reference sequence with known GRN assignments is identified.
3. **Sequence alignment**: Pairwise alignment against the reference using Biopython's pairwise2 module.
4. **GRN transfer**: Positions are labeled according to the alignment, assigning GRN coordinates to each residue.
5. **Graph construction**: Binding pocket residues (those with GRN positions in the reference graph) are extracted, and the family-specific edge topology is applied.
6. **Feature generation**: Protein language model embeddings are computed for each node using Ankh-large.

This pipeline enables LAMBDA to accept any opsin sequence as input, requiring only the amino acid sequence. The GRN system provides the structural context—which positions are in the binding pocket and how they relate spatially—while the pLM embeddings capture the biochemical properties of the specific amino acids at each position.

## Model and Training

LAMBDA uses a graph neural network that processes binding pocket graphs to predict λmax for multiple retinal conformations.

**Input representation.** Each node in the graph corresponds to a binding pocket residue and carries two feature types: (1) a positional encoding derived from the GRN coordinate, capturing the residue's structural location relative to retinal, and (2) a 1536-dimensional embedding from Ankh-large, capturing evolutionary and biochemical context. The positional encoding is learned during training, allowing the model to discover which binding pocket locations matter most for spectral tuning.

**Message passing.** The encoder uses a Graph Convolutional Network (GCN) with residual connections to update node representations through neighborhood aggregation. Each GCN layer propagates information between connected residues, integrating local interactions across the binding pocket.

**Pooling.** A position-aware attention pooling mechanism with 32 heads generates a global representation from node-level features. This mechanism uses the GRN positional encodings to compute attention weights over residue contributions, enabling the model to learn multiple complementary views of position-dependent importance.

**Multi-task output.** The pooled representation feeds into five regression heads that predict λmax across retinal conformations and protonation states. For both 11-cis and all-trans retinal, two heads operate in parallel: one predicts λmax for protonated Schiff base states only (λmax > 400 nm), while the other predicts λmax across the full spectral range including deprotonated states. A fifth head predicts λmax for deprotonated Schiff base (UV-absorbing states, λmax < 400 nm). Protonation state is classified by a threshold at 400 nm, chosen based on the bimodal distribution of λmax values in the training data: opsins with λmax below this boundary are classified as having a deprotonated Schiff base, and the corresponding UV head prediction is used; otherwise the protonated-state head prediction is reported. All prediction heads share the same encoder—the same learned representation of binding pocket structure feeds all outputs. This shared-encoder architecture reflects the biological reality that a single binding pocket determines spectral properties across chromophore states, and enables generalization to additional retinal conformations (9-cis, 13-cis) as training data become available.

**Loss function.** Training uses a combined loss that sums mean squared error across all prediction tasks, with masking to handle missing labels (e.g., type I opsins lack 11-cis measurements). The full-range heads (which predict across both protonated and deprotonated states) receive half weight relative to the state-specific heads, as they serve primarily as a regularizer that encourages globally consistent predictions. The split at 400 nm into separate protonated and UV heads allows the model to learn distinct distributions for each protonation state rather than fitting a single non-Gaussian distribution.

**Optimization.** The model is trained with a learning rate of 2×10⁻⁴ using a ReduceLROnPlateau schedule that decays the rate to a minimum of 1×10⁻⁷ when validation loss plateaus. A batch size of 2 is used; larger batch sizes consistently degraded learning, likely because the small dataset and high variance across opsin families benefit from the noisier gradient estimates of very small batches.

**Regularization.** To prevent overfitting, dropout is applied within the GCN encoder and attention pooling layers, weight decay on all parameters, and early stopping based on validation loss with a patience of 50 epochs.

**Data splitting.** The combined dataset is divided into training (80%), validation (10%), and test (10%) sets using stratified sampling. Stratification is performed at two levels: by source dataset (VPOD, Inoue, hCRBPII) and by species within each dataset. The 31 dual-state opsins are assigned exclusively to the training split to maximize exposure to dual-state measurements. The remaining data are split to ensure representative sampling across opsin types and taxonomic groups. For species with fewer than 3 samples, random assignment is used to avoid stratification failures.

**Evaluation metrics.** Mean absolute error (MAE) in nanometers after denormalizing predictions to the original λmax scale, coefficient of determination (R²), and the fraction of predictions within 5 nm and 10 nm of experimental values are reported. For protonation state, classification accuracy is reported.

**Normalization.** Target values are normalized to [0, 1] using the observed ranges in training data: 11-cis λmax (401–611 nm), all-trans λmax (436–644 nm), and UV λmax (340–399 nm).

## Opsin Atlas

A comprehensive sequence collection was assembled by mining the NCBI non-redundant protein database to achieve broad coverage of opsin diversity across all major subfamilies in both type I and type II families.

A phylogenetically informed sampling strategy was designed using BLAST searches seeded from representative sequences of each major opsin subfamily. Each query was searched against NCBI nr using BLASTp, retrieving up to 5,000 hits per query. Because queries from related subfamilies often retrieve overlapping sets of sequences, each protein was assigned to the subfamily of the query to which it had the highest sequence identity. This similarity-based classification supersedes NCBI annotations, which frequently use "bacteriorhodopsin" as a generic label for diverse microbial rhodopsins regardless of their phylogenetic placement—for example, 46% of type I sequences carry "bacteriorhodopsin" in their NCBI annotation, yet similarity-based classification distributes these across eight distinct families, reflecting their actual evolutionary relationships. All per-query results were then merged into unified type I and type II tables.

For type I (microbial) opsins, 13 query sequences covering 10 functional subfamilies were used:

| Subfamily | Query | Accession |
|-----------|-------|-----------|
| Proton pump | Bacteriorhodopsin (*H. salinarum*) | P02945 |
| Proton pump | Archaerhodopsin-1 (*H. salinarum*) | P69051 |
| Proton pump | Archaerhodopsin-2 (*H. salinarum*) | P29563 |
| Chloride pump | Halorhodopsin (*H. salinarum*) | P0DMH7 |
| Sensory rhodopsin I | SRI (*H. salinarum*) | P0DMH8 |
| Sensory rhodopsin II | SRII (*N. pharaonis*) | P42196 |
| Cation channel | Channelrhodopsin-1 (*C. reinhardtii*) | A0A2K3CXC9 |
| Cation channel | Channelrhodopsin-2 (*C. reinhardtii*) | Q8RUT8 |
| Anion channel | GtACR1 (*G. theta*) | L1JRS2 |
| Green-absorbing PR | Proteorhodopsin (*SAR86*) | Q9F7P4 |
| Blue-absorbing PR | Proteorhodopsin (*HOT75*) | Q9AFF7 |
| Heliorhodopsin | HeR (metagenome) | A0A2P2C3K4 |
| Xanthorhodopsin | Xanthorhodopsin (*S. ruber*) | Q2S2F8 |

For type II (animal) opsins, 14 query sequences covering 12 subfamilies were used:

| Subfamily | Query | Accession |
|-----------|-------|-----------|
| Rod opsin | Rhodopsin (*Bos taurus*) | P02699 |
| Rod opsin | Rhodopsin (*H. sapiens*) | P08100 |
| Short-wavelength-sensitive cone opsin 1 (SWS1) | OPN1SW (*H. sapiens*) | P03999 |
| Medium-wavelength-sensitive cone opsin (MWS) | OPN1MW (*H. sapiens*) | P04001 |
| Long-wavelength-sensitive cone opsin (LWS) | OPN1LW (*H. sapiens*) | P04000 |
| Encephalopsin (OPN3) | OPN3 (*H. sapiens*) | Q9H1Y3 |
| Melanopsin (OPN4) | OPN4 (*H. sapiens*) | Q9UHM6 |
| Neuropsin (OPN5) | OPN5 (*H. sapiens*) | Q6U736 |
| Parapinopsin | Parapinopsin (*I. punctatus*) | O42266 |
| Pinopsin | Pinopsin (*G. gallus*) | P51475 |
| Peropsin | Peropsin (*H. sapiens*) | O14718 |
| Retinal G protein-coupled receptor (RGR) | RGR (*H. sapiens*) | P47804 |
| Teleost multiple tissue opsin (TMT) | TMT opsin (*D. rerio*) | R9R6C7 |
| Vertebrate ancient opsin (VA) | VA opsin (*O. masou*) | O13018 |

The retrieved sequences underwent quality filtering. Sequences exceeding 500 residues (likely multi-domain proteins or fusion constructs) were removed. To eliminate redundancy while preserving phylogenetic diversity, exact sequence duplicates were removed and cases where the same gene from the same species appeared under multiple accession numbers were filtered, while retaining orthologous sequences from different species. Most critically, each sequence was validated through the GRN annotation pipeline, retaining only sequences that (1) contained the conserved Schiff base lysine at the expected GRN position (7.43 for type II, 7.50 for type I), and (2) achieved sufficient GRN coverage to generate the input graphs for LAMBDA.

For prediction, LAMBDA was applied with output heads matched to each opsin type's native chromophore configuration. Type I opsins, which bind all-trans retinal, received predictions for λmax^AT. Type II opsins received predictions for their dark-state chromophore (λmax^11-cis), but λmax^AT was additionally predicted for all type II sequences. This dual prediction for type II opsins—though training data for their all-trans state is limited to the 31 dual-state sequences—enables estimation of spectral shifts upon photoactivation, a property critical for optogenetic applications. All data management, sequence processing, and GRN annotation were performed using ProtOS.

# Results

LAMBDA was evaluated on the held-out test set (10% of data, stratified by dataset and species). Performance is reported for each prediction target alongside published results from family-specific methods where available.

| Metric | LAMBDA | OPTICS | Inoue et al. | RhoMax |
|--------|--------|--------|--------------|--------|
| **MAE ± std (nm)** | | | | |
| Type II (11-cis) | **5.18 ± 0.82** | 5.49 | -- | -- |
| Type I (all-trans) | **6.86 ± 0.89** | -- | 7.80 | 6.83 |
| hCRBPII (all-trans, n=6) | **5.84 ± 1.13** | -- | -- | -- |
| **R²** | | | | |
| Type II (11-cis) | **0.972** | 0.964 | -- | -- |
| Type I (all-trans) | **0.894** | -- | -- | -- |
| hCRBPII (all-trans, n=6) | **0.978** | -- | -- | -- |
| **Other** | | | | |
| Protonation acc. | 98.3% | -- | -- | -- |
| Within 5 nm (11c) | 77.6% | -- | -- | -- |
| Within 5 nm (AT) | 52.7% | -- | -- | -- |
| Within 10 nm (11c) | 93.9% | -- | -- | -- |
| Within 10 nm (AT) | 83.9% | -- | -- | -- |

Evaluation protocols differ across methods. LAMBDA reports held-out test set performance (10% stratified split). OPTICS uses k-fold cross-validation on VPOD WDS (n=1,211) with amino acid property encoding (Frazer & Oakley 2025). Inoue et al. report within-family evaluation on 884 type I sequences (Karasuyama et al. 2018). RhoMax reports median absolute error across 4 family-aware cross-validation splits on 884 type I sequences (mean AE = 10.45 nm; Sela et al. 2024) and requires structures as input. Despite these protocol differences, all methods operate on comparable dataset sizes and the reported accuracies are broadly comparable. The hCRBPII test set contains only 6 samples, limiting statistical power.

LAMBDA achieved a mean absolute error of 5.18 ± 0.82 nm for type II opsins (11-cis retinal, n=119) and 6.86 ± 0.89 nm for type I opsins (all-trans retinal, n=87) on the held-out test set (Figure 4.5). These values reflect the model's actual output, using the protonation-state classification described in Methods to select between head predictions. The model correctly classified protonation state with 98.3% accuracy. For hCRBPII (42 training, 3 validation, 6 test samples), the model achieved an MAE of 5.84 ± 1.13 nm (R² = 0.978)—comparable to opsin performance despite the fundamentally different lipocalin fold. The small test set limits statistical confidence, but the result demonstrates that the binding domain graph framework accommodates non-opsin folds when training examples and appropriate GRN mappings are available.

These results place LAMBDA's accuracy in the same range as the best family-specific methods while predicting across chromophore conformations and protein families from a unified model. That joint training does not degrade within-family accuracy suggests that spectral tuning principles transfer across folds and that the cross-family position mapping captures genuine structural equivalences rather than introducing noise. A systematic underestimation of far-red absorption maxima, where the model saturates near 610 nm, represents the primary limitation and is addressed in the Discussion.

The validated model was applied to 47,700 opsin sequences assembled from NCBI to produce the Opsin Atlas, spanning 20,061 type I opsins across 8 families and 27,639 type II opsins across 10 families (12 subfamilies).

**Type I opsins** (Figure 4.6A). Predictions of λmax^AT for native all-trans retinal show a mean of 519.5 ± 28.1 nm (range: 436–644 nm) across 20,061 sequences. The aggregate distribution approximates a single Gaussian centered near 520 nm, reflecting the relatively homogeneous spectral tuning of microbial rhodopsins around their shared proton-pumping or sensory functions. Within this envelope, individual subfamilies occupy characteristic spectral ranges: proton pump: 532.4 ± 23.8 nm (n=5,108, 25.5%), heliorhodopsin: 498.5 ± 37.6 nm (n=4,993, 24.9%), green-absorbing proteorhodopsin: 520.7 ± 16.3 nm (n=4,942, 24.6%), xanthorhodopsin: 527.3 ± 12.5 nm (n=4,304, 21.5%), blue-absorbing proteorhodopsin: 537.0 ± 20.7 nm (n=310, 1.5%), cation channel: 486.1 ± 25.6 nm (n=120, 0.6%), chloride pump: 525.6 ± 16.1 nm (n=103, 0.5%), sensory rhodopsin II: 508.5 ± 25.3 nm (n=89, 0.4%), sensory rhodopsin I: 509.0 ± 24.5 nm (n=86, 0.4%), and anion channel: 478.3 ± 15.7 nm (n=6, <0.1%). Four subfamilies—proton pumps, heliorhodopsins, green-absorbing proteorhodopsins, and xanthorhodopsins—together account for 96.5% of all type I sequences. Proteorhodopsins show a characteristic bimodal distribution, with distinct blue-absorbing and green-absorbing sub-populations reflecting ecological adaptation to different ocean depths. Heliorhodopsins span the widest spectral range but are blue-shifted relative to the classical proton pumps.

**Type II opsins** (Figure 4.6B). In contrast to the unimodal type I distribution, the 11-cis dark-state predictions for 27,639 type II opsins (mean 467.0 ± 42.9 nm, range 344–583 nm) resolve into multiple distinct peaks corresponding to the spectral classes of animal vision. Per-subfamily statistics are: rod opsin: 481.7 ± 24.0 nm (n=5,321, 19.3%), neuropsin: 459.8 ± 36.7 nm (n=4,731, 17.1%), encephalopsin: 465.7 ± 21.5 nm (n=4,317, 15.6%), cone SWS1: 413.0 ± 42.3 nm (n=3,237, 11.7%), melanopsin: 495.2 ± 24.6 nm (n=2,584, 9.3%), cone MWS: 532.4 ± 31.3 nm (n=2,473, 8.9%), RGR: 456.0 ± 25.3 nm (n=2,041, 7.4%), peropsin: 448.3 ± 20.5 nm (n=1,449, 5.2%), parapinopsin: 430.1 ± 33.6 nm (n=1,255, 4.5%), VA opsin: 436.6 ± 32.4 nm (n=102, 0.4%), cone LWS: 486.8 ± 26.3 nm (n=64, 0.2%), and pinopsin: 453.1 ± 34.2 nm (n=64, 0.2%). The dominant peak near 500 nm is formed by rod opsins, the primary dim-light photoreceptors. A second major cluster in the UV-blue range below 430 nm comprises cone SWS1 opsins and parapinopsin, which detect short wavelengths for color discrimination and circadian entrainment. Cone MWS opsins form a distinct peak near 530 nm, while melanopsin clusters near 495 nm. Non-visual opsins—encephalopsin and neuropsin—populate the intermediate blue-green range. This multimodal structure reflects the evolutionary diversification of animal opsins into specialized spectral channels for color vision, scotopic sensitivity, and non-visual photoreception.

LAMBDA also predicts λmax^AT for all type II opsins—the absorption maximum with all-trans retinal, as occurs during photoactivation or stably in bistable opsins. The λmax^AT distribution is systematically red-shifted relative to 11-cis predictions and collapses from the multimodal 11-cis pattern into a bimodal shape: a primary peak near 520 nm dominated by rod and long-wavelength cone opsins, and a secondary shoulder near 490 nm comprising SWS1, neuropsin, RGR, and encephalopsin. This dual-state analysis applies only to type II opsins; the model produces 11-cis predictions for type I sequences, but these lack biological basis and training signal and are omitted from the atlas.

For type II opsins, predictions for both chromophore states enable estimation of the spectral shift upon activation. Δλ is computed using the protonated Schiff base predictions only (λmax > 400 nm), excluding UV-absorbing deprotonated states—the goal is to isolate the spectral shift arising from binding pocket interactions with retinal, not the large blue shift that accompanies Schiff base deprotonation (which moves any opsin into the UV regardless of pocket composition). These Δλ values are derived quantities: the all-trans prediction for type II opsins is trained on only 31 dual-state sequences and has not been independently validated, so the spectral shift estimates should be interpreted as indicative rather than quantitative. The signed difference Δλ = λmax^AT − λmax^11-cis ranges from −108 to +148 nm (mean +18.8 ± 33.0 nm; Figure 4.6C), with 23% of predictions showing a blue shift (negative Δλ). Per-subfamily means reveal consistent differences: VA opsin (+47.5 ± 28.4 nm), cone SWS1 (+41.7 ± 30.7 nm), RGR (+34.5 ± 22.1 nm), and peropsin (+24.8 ± 18.6 nm) show the largest red shifts. Rod opsin (+17.7 ± 20.0 nm), neuropsin (+17.5 ± 28.7 nm), and cone LWS (+5.2 ± 15.8 nm) show modest shifts. Cone MWS is the only subfamily with a negative mean (−37.9 ± 29.9 nm), indicating that these opsins are predicted to blue-shift upon isomerization to all-trans. Melanopsin (+21.6 ± 43.4 nm) shows the widest spread, spanning both directions.

The atlas is available as a supplementary resource, providing predicted λmax values for both chromophore states, spectral separation estimates, GRN annotations, and taxonomic metadata for each sequence.

# Discussion

LAMBDA demonstrates that binding domain analysis—representing the chromophore environment as a graph aligned on retinal—enables spectral prediction across protein families from a single model. Where OPTICS and the Inoue model extract features from amino acid sequences, and RhoMax uses structure-derived graphs limited to type I opsins, LAMBDA defines graph nodes by their spatial relationship to the chromophore rather than by sequence position. GRN systems make this a structure-based model that accepts sequence as input.

Accuracy comparable to the best family-specific methods (5.18 nm MAE for type II, 6.86 nm for type I) is achieved while simultaneously predicting multiple chromophore states from a unified model—a capability that sequence-based approaches, which are inherently fold- and conformation-specific, cannot provide. The improvement over OPTICS' sequence-based approach (5.49 nm cross-validated) suggests that binding pocket structure provides more robust features than lineage-specific sequence patterns. Within this representation, pLM embeddings (Ankh-large, 1536 dimensions) describe *what* occupies each binding pocket position while the graph describes *where* those positions are relative to retinal—the model depends on both, but the graph structure is what enables cross-family learning.

The most informative limitation is LAMBDA's systematic underestimation of far-red absorption maxima. The model saturates at approximately 610 nm, compressing the 590–650 nm range even when these samples are included in training. This reflects physics rather than architecture: the relationship E = hc/λ means that energy gaps shrink dramatically at long wavelengths, and at these small scales spectral tuning becomes dominated by subtle quantum mechanical effects—precise positioning of partial charges, fine details of hydrogen bonding, electronic polarization—that produce disproportionately large wavelength shifts. The binding pocket graph captures residue identity and spatial arrangement, but not the electrostatic precision required to resolve 0.03 eV differences in excitation energy. Additional training data are unlikely to help; the features that distinguish a 600 nm opsin from a 640 nm opsin may require explicit quantum mechanical modeling. For optogenetic engineering, LAMBDA can identify candidates as deeply red-shifted but cannot rank them within that regime—precise spectral positioning above 590 nm requires either experimental characterization or hybrid approaches combining machine learning with QM/MM refinement.

More generally, the model cannot generalize to novel protein folds without training examples—the relationship between residue positions and spectral effects is fold-specific and must be learned. Predicting properties beyond λmax—quantum yield, photocycle kinetics, ion selectivity—would require expanded training data and potentially modified architectures, and the model does not account for post-translational modifications or lipid environment effects that may influence spectral properties in vivo.

Despite these limitations, the atlas enables computational screening across the full spectral range—identifying red-shifted opsins for tissue penetration, spectrally separated pairs for multiplexed control, and opsins with large dark-to-activated spectral separation for bidirectional tools. LAMBDA is the first method to predict the spectral separation between dark and activated states systematically. Across type II opsins, most subfamilies red-shift upon activation, but cone MWS opsins are predicted to blue-shift (mean −37.9 nm)—the only subfamily with a negative mean. This pattern reflects a convergence: all-trans predictions cluster near 520 nm regardless of subfamily, while 11-cis dark-state tuning is subfamily-specific. Opsins absorbing at short wavelengths in the dark state accumulate large positive Δλ because the all-trans state is red-shifted relative to their blue dark state; cone MWS opsins, already at ~530 nm, shift in the opposite direction. Cone LWS opsins, at the longest wavelengths, show minimal shift—selective pressure acts on the 11-cis dark state for wavelength discrimination, not on the activated conformation.

Among subfamilies with the largest positive Δλ, VA opsin, melanopsin, and RGR are known bistable photopigments or photoisomerases; their larger predicted separations could indicate that maintaining spectrally distinct states is linked to bistable signaling. Rod opsins show modest shifts, consistent with their rapid meta-II decay. These correspondences are suggestive but not conclusive—Δλ is a derived quantity, not a training target, so if these trends hold under experimental scrutiny they would provide indirect validation that the model captures biologically meaningful spectral tuning. The atlas highlights specific candidates for optogenetic applications: a heliorhodopsin from *Candidatus Kerfeldbacteria* (PIS41810) is predicted as the most red-shifted type I opsin at 644 nm, deep into the optical window for tissue penetration; an SWS1 opsin from *Hylia prasina* (NWU40863) as the most blue-shifted type II opsin at 344 nm; and for bidirectional control, a peropsin from the beluga whale *Delphinapterus leucas* (XP_022422838) shows the largest predicted spectral separation at +148 nm (458 nm dark state, 607 nm activated), followed by melanopsins from *Sturnira hondurensis* (+146 nm) and *Phyllostomus discolor* (+143 nm). The ability to predict how mutations affect both chromophore states enables engineering strategies that optimize spectral separation, not just dark-state absorption—candidates that would be impractical to discover through experimental screening alone.

Beyond opsin families, the hCRBPII results demonstrate that LAMBDA's framework extends to entirely different protein folds. Despite the fundamentally different lipocalin architecture, the model achieves accuracy comparable to opsin predictions when training examples are available. The boundary is principled: extending to a new fold requires defining its binding pocket graph and providing representative training data, rather than building a separate model.

An open limitation is the role of solvent in spectral tuning. In opsins, the binding pocket is largely occluded from bulk water, and the model may implicitly learn solvation effects through correlations with surrounding residue identities. In solvent-exposed binding pockets, water molecules directly contact the chromophore and modulate its absorption in ways that depend on occupancy and orientation rather than protein sequence. pLM embeddings may capture some of this context—residues flanking a solvent-exposed pocket differ systematically from those in a buried one—but are unlikely to fully resolve the spectral effects of specific water configurations. For proteins with open binding sites, this remains a source of prediction error that neither additional sequence data nor larger embeddings are likely to eliminate without explicit solvation modeling.

When structures are available, binding pocket graphs can be derived directly. For the sequence-only path—which enables the Opsin Atlas and any large-scale application—LAMBDA depends on the standardized positional annotation developed in the preceding chapter: without MOGRN coordinates for microbial rhodopsins, the graph representation that enables cross-family learning from sequence would not exist; without the chromophore-centered mapping between MOGRN and Ballesteros-Weinstein positions, type I and type II opsins could not share graph nodes. This dependency is also the source of extensibility—any protein family with a retinal-binding pocket and an appropriate GRN system can be incorporated. The binding pocket representation further enables spectral prediction for engineered proteins that preserve the chromophore environment while modifying other regions. The following chapter explores this through rhodozyme design—engineering the intracellular regions of type II opsins to incorporate enzymatic active sites that become accessible only upon light activation, with spectral properties that remain predictable because the binding domain is preserved.

# Rhodozyme — Light-Activated Enzyme Design

Type II opsins undergo a conformational change upon photon absorption. In the dark state, the seven transmembrane helices pack tightly, with the intracellular face closed. Photoisomerization of retinal triggers a cascade of sidechain rearrangements at conserved positions—the microswitch residues—that propagate the conformational change from the binding pocket to the intracellular surface. These microswitches are well characterized: they include residues on TM3, TM5, TM6, and TM7 whose rotameric states differ between the inactive and active conformations. The net effect is that TM6 moves 10–14 Å outward at its cytoplasmic end, exposing a cavity that in native rhodopsin binds the G protein transducin. This cavity is transient, stereochemically defined, and gated by light. The retinal binding pocket—buried in the transmembrane core—determines the absorption wavelength, while the intracellular cavity is a separate surface. If an enzymatic active site could be placed on this intracellular face, the result would be a light-gated enzyme whose activation wavelength is set by the rhodopsin scaffold.

This chapter explores whether such a design can be assembled computationally. The core technique—scaffolding a protein backbone around a fixed arrangement of catalytic residues (a theozyme)—was established by the Baker lab in RFdiffusion2 (Ahern et al., 2025). Their protocol generates 100 backbone designs per theozyme, fits 8 sequences per design with LigandMPNN, validates all candidates with structure prediction, and picks the best by predicted confidence. We follow the same protocol. The rhodozyme application adds an earlier step and an additional difficulty: the scaffold is not a novel fold but an existing rhodopsin in its active conformation, and the mask that separates fixed from designable regions must preserve both the retinal binding pocket (for light sensitivity) and the theozyme positions (for catalysis). This means candidate selection happens at two stages rather than one—first at theozyme placement (is the geometry compatible with the rhodopsin intracellular face?) and again at validation (does the predicted structure maintain both the rhodopsin fold and the catalytic geometry?). This is where ProtOS contributes. The models are difficult; the data integration between them is not, provided the intermediate representations—structures, annotations, graphs, masks—are managed consistently. ProtOS handles structure loading, GRN annotation, geometric matching, mask construction, and job dispatch to the Model Manager. The theozyme extraction and placement code (Steps 2–3) required only a small, application-specific addition outside ProtOS's general processor framework.

The design proceeds in six steps. Each step can produce multiple outputs that fan out into the next: the geometric search yielded 8 candidate placements, of which we selected 2 based on geometric quality and helix involvement; for each we generated 50 backbone designs with RFdiffusion, sampled 8 sequences per design with LigandMPNN, and predicted all resulting candidates with Boltz2. For placement 00, this produced 50 backbone designs, 405 sequences (from 45 successful designs), and 307 Boltz2 structure predictions.

## Step 1 — Starting Structures

Three structures are required. The first is a rhodopsin in its active conformation—metarhodopsin II (PDB: 3PQR), which captures the open intracellular state with TM6 displaced outward. The dark-state structure of the same protein (PDB: 1U19, bovine rhodopsin at 2.2 Å) serves as a reference for the conformational change: superimposing the two reveals the TM5/TM6 displacement that creates the intracellular cavity. The second is a reference enzyme bound to its substrate in a catalytic intermediate. We use bovine trypsin (PDB: 2AGE), an acyl-enzyme intermediate at 1.15 Å resolution with succinyl-AAPR covalently bound to the catalytic serine. This structure captures the triad geometry mid-catalysis—Ser195, His57, and Asp102 in their active arrangement, with the substrate positioned at the reaction center. The enzyme provides the catalytic geometry to be transplanted; the rhodopsin provides the scaffold that gates access to it.

> **Figure 5.1** — Input structures and the design premise. (A) Dark-state bovine rhodopsin (1U19, gray) superimposed with active-state metarhodopsin II (3PQR, terracotta), showing TM5/TM6 displacement and the intracellular cavity that opens upon activation. Retinal in rust. (B) Bovine trypsin acyl-enzyme intermediate (2AGE) with the catalytic triad as sticks (Ser195, His57, Asp102) and hydrogen-bond distances shown as dashed lines. Covalently bound substrate (succinyl-AAPR) in ochre.

## Step 2 — Theozyme Extraction

The theozyme is the minimal catalytic unit: the three sidechain positions that perform chemistry. For a serine protease, these are the nucleophilic serine, the general base histidine, and the orienting aspartate. From the trypsin structure, we extract three quantities per residue: the Cα coordinate, the Cα→Cβ vector (sidechain direction), and the residue identity. The pairwise Cα–Cα distances define a triangle; the Cβ vectors define the orientation of each sidechain within that triangle. Together, these six quantities (three positions, three directions) specify the catalytic geometry that must be reproduced in the rhodopsin scaffold.

The catalytic triad in 2AGE shows the characteristic hydrogen-bond relay: the Ser195 hydroxyl is 3.04 Å from His57 Nε2, and His57 Nδ1 is 2.75 Å from Asp102 Oδ2. These functional-atom distances define the active geometry. The Cα triangle that supports this arrangement has sides of 8.3 Å (Ser195–His57), 6.5 Å (His57–Asp102), and 10.1 Å (Ser195–Asp102).

> **Figure 5.2** — Theozyme extraction. The catalytic triad shown in the trypsin active site (gray cartoon) with the geometric abstraction overlaid: Cα positions as spheres, Cα→Cβ direction vectors as arrows, pairwise distances as dashed lines. Distances: Ser195–His57 = 8.3 Å, His57–Asp102 = 6.5 Å, Ser195–Asp102 = 10.1 Å. This triangle and these vectors are the input to the placement search—everything else about trypsin is discarded.

## Step 3 — Theozyme Placement

The intracellular face of active rhodopsin is identified using GRN annotation. Residues on TM helix ends that face the cytoplasm (TM1 ≥ 1.60, TM3 ≥ 3.55, TM5 ≥ 5.68, TM7 ≥ 7.53) and intracellular loop residues (ICL1, ICL2, ICL3, H8) form the candidate region. An exhaustive search over all triplets in this region finds positions whose Cα triangle matches the theozyme triangle within 2 Å RMSD. Candidates passing the distance filter are then tested for sidechain direction: the source Cα triangle is Kabsch-aligned onto the candidate, and the rotated Cβ vectors are compared. Matches within 30° are retained.

An additional constraint requires at least one residue on TM5 or TM6—the helices that move during activation. This ensures that the catalytic geometry depends on the active conformation: in the dark state, with TM6 packed inward, the triad distances break. The enzyme turns on with light and off without it.

This step has no equivalent in the Baker lab protocol. Ahern et al. (2025) start from a theozyme and generate a novel scaffold around it; the placement is implicit in the diffusion process. Here, the scaffold is fixed—we must find positions in an existing structure that can accommodate the catalytic geometry. The search produced 8 candidates. We selected placements 0 and 2 based on triangle RMSD, Cβ vector alignment, and the involvement of TM6 residues. For placement 00, the theozyme maps to Ser-230, His-245, and Asp-250 on the intracellular face of 3PQR. The placement reproduces the trypsin Cα triangle exactly (8.3 / 6.5 / 10.1 Å) and preserves the hydrogen-bond distances (Ser–His 3.04 Å, His–Asp 2.75 Å).

> **Figure 5.3** — Theozyme placement on the rhodopsin scaffold. (A) The placement structure (3PQR with theozyme mutations) viewed from the intracellular face. TM helices in gray, theozyme residues (Ser-230, His-245, Asp-250) shown as sticks with Cα spheres in green. Retinal in rust. (B) Same structure from an alternate angle showing the theozyme sidechain arrangement relative to the transmembrane core.

## Step 4 — Backbone Design with RFdiffusion

Point mutations at the matched positions would introduce the correct residues but not the correct backbone geometry. The surrounding loops and helix termini must accommodate the theozyme. RFdiffusion generates backbone designs under constraints.

The mask defines what is fixed and what is designed. The locked regions comprise the TM helices (TM1–TM2, TM3, TM4–TM5, TM6–TM7), preserving the transmembrane core, the retinal binding pocket, and the microswitch residues that mediate the activation mechanism—210 residues in total. The three theozyme positions (Ser-230, His-245, Asp-250) are locked with their sidechain atoms explicitly constrained. Between these locked segments, 116 residues are free for RFdiffusion to redesign—these correspond to the intracellular loops (ICL1–3), the theozyme-surrounding loops, and the C-terminal region. Retinal and the tetrapeptide substrate are included as ligand context.

We generate 50 designs per placement. The smaller number (vs. 100 in Ahern et al.) reflects the additional constraint: the locked TM scaffold reduces the conformational search space relative to unconstrained scaffolding, but also makes each design harder to solve.

> **Figure 5.4** — RFdiffusion mask. The mask applied to the rhodopsin scaffold: TM helices and theozyme positions locked (gray), intracellular loops free for design (terracotta). Theozyme residues marked as green spheres. Retinal (rust) sits in the locked transmembrane core, unchanged by design. 210 residues are locked; 116 are designed.

## Step 5 — Sequence Design with LigandMPNN

Each RFdiffusion backbone specifies a fold but not a sequence. LigandMPNN generates amino acid sequences compatible with the designed backbone while accounting for the retinal cofactor and the substrate ligand. The theozyme residues (Ser, His, Asp at the fixed positions) are provided as constraints—LigandMPNN designs the rest.

Following Ahern et al. (2025), 8 sequences are sampled per backbone at temperature T=0.1. The retinal SMILES and substrate SMILES are included in the LigandMPNN input so that the designed sequence accounts for both the cofactor that enables light activation and the substrate that the enzyme should bind. For placement 00, this produced 405 sequences across 45 successful backbone designs. The top candidate has 72.7% sequence identity to wild-type rhodopsin (3PQR), with 89 mutations concentrated in the redesigned loop regions. The locked TM helices retain the native sequence almost entirely.

> **Figure 5.5** — Sequence design. Sequence alignment of the top candidate against wild-type rhodopsin (3PQR). Identical positions in gray, mutations highlighted. Mutations cluster in the redesigned loop regions; the TM helices are largely unchanged. Theozyme positions (Ser-230, His-245, Asp-250) marked.

## Step 6 — Structure Prediction with Boltz2

Each designed sequence is predicted with Boltz2 to evaluate whether the intended fold and catalytic geometry are maintained. The prediction includes the protein chain, retinal (as covalent cofactor), and substrate. This is the second filtering stage. Ahern et al. (2025) rank candidates by predicted confidence; we evaluate on two criteria specific to the rhodozyme constraint.

The first criterion is theozyme alignment. We superimpose each predicted structure onto the placement from Step 3, aligning only the three theozyme residues (all atoms, not just Cα), and measure the RMSD. This tests whether the catalytic geometry survived the design-predict cycle. Across all 307 predictions, we rank candidates by this theozyme all-atom RMSD. The second criterion is pLDDT of the Boltz2 prediction, reported per-residue. We examine pLDDT separately for the locked regions (which should score high, as they reproduce known structure) and the designed regions (where low confidence indicates the model is uncertain about the backbone).

The comparison between predicted and reference theozyme geometry requires more than a rigid-body superposition. The Cα positions may align well while the sidechain rotamers differ—a common outcome when Boltz2 resolves the local environment differently from the reference. We therefore allow sidechain rotation around the Cα–Cβ bond axis when assessing the match, comparing the functional-atom distances (Ser OG – His Nε2, His Nδ1 – Asp Oδ) rather than insisting on identical χ1 angles. This relaxed comparison reflects the physical reality: what matters for catalysis is the hydrogen-bond relay, not the exact dihedral.

An unexpected observation emerged from inspecting the validated designs. Because the theozyme was placed flatly on the ICL3 surface, the designed loops could fold in two ways: with the substrate-binding face directed inward (into the helical bundle interior, as originally intended) or outward (toward the cytosolic side of ICL3). The top-ranking prediction adopted the outward orientation—the catalytic face and substrate-binding site point toward the cytoplasm rather than into the transmembrane pocket. The theozyme geometry itself is preserved in both orientations; what differs is where the substrate approaches. The outward orientation may be more favorable, since substrate access is not occluded by the surrounding helices. Crucially, the light-gating mechanism is unaffected: the triad residues sit on positions that are only in register when TM5/TM6 adopt the active conformation, regardless of which face the binding site presents.

We selected the top candidate by balancing theozyme RMSD (2.44 Å all-atom, 0.51 Å Cα-only, rank 6 of 307) against the highest global pLDDT in the top 10 (91.7). The overall backbone RMSD to the parent rhodopsin is 0.84 Å, confirming that the fold is preserved. The pLDDT breakdown shows high confidence in the locked TM core (94.9 mean) and good confidence in the designed loops (85.6 mean), with the theozyme residues at 72.1—lower, as expected for residues at a designed interface.

The catalytic geometry in the predicted structure shows the Ser–His distance at 4.03 Å and the His–Asp distance at 3.18 Å, compared to 3.04 Å and 2.75 Å in the reference placement. These distances are longer than the ideal hydrogen-bond geometry but within the range that could be optimized by further design iterations or molecular dynamics relaxation.

Because the retinal binding pocket is preserved by the mask, the Schiff base linkage to Lys-296 is intact in the predicted structure. LAMBDA can predict the spectral properties of retained candidates—the rhodozyme absorbs at the same wavelength as the parent rhodopsin.

> **Figure 5.6** — Boltz2 evaluation of the top candidate. (A) Predicted structure overlaid on the parent rhodopsin. Gray: reference TM core (locked regions from 3PQR). Terracotta: designed loop regions (Boltz2 prediction). Retinal and Schiff base Lys-296 in rust. Overall backbone RMSD: 0.84 Å. (B) Catalytic geometry comparison: predicted theozyme (green) overlaid on reference placement (gray), with catalytic interaction distances shown. Reference: Ser–His 3.04 Å, His–Asp 2.75 Å. Predicted: Ser–His 4.03 Å, His–Asp 3.18 Å. (C) Per-residue pLDDT confidence. Locked TM core: 94.9 mean. Designed loops: 85.6 mean. Theozyme residues: 72.1 mean. Global: 91.7.

We consider the computational design a success: Boltz2 predicts a well-folded structure (pLDDT 91.7) that preserves the rhodopsin fold (backbone RMSD 0.84 Å) and maintains the theozyme geometry within optimizable range. The most difficult step in the workflow is not the AI-driven design or validation, but the original theozyme placement. This is a combinatorial and structural biology problem that precedes the Baker lab protocol entirely. Each placement generates a full cascade of 50 backbone designs × 8 sequences × structure predictions—for two placements, this already produces over 600 candidates to evaluate. Selecting many placements without careful geometric and structural reasoning would produce a vast screening space that is computationally expensive and difficult to interpret. The placement decision requires intuition about protein geometry, knowledge of the rhodopsin conformational cycle, and judgment about which helix positions can support catalytic function. The same kind of expert judgment applies at the end of the pipeline: interpreting predicted structures, assessing whether catalytic distances are close enough, and deciding which candidates merit experimental follow-up. The AI models automate the generative steps; the structural biology reasoning that frames and interprets them remains human.

## Integration

Every step in this workflow—structure fetching, annotation, geometric matching, mask construction, model submission, and result registration—runs through ProtOS. The PDB structures are fetched and loaded by the StructureProcessor. GRN annotation identifies the intracellular face. The theozyme placement search operates on the annotated coordinates. The mask for RFdiffusion is built from GRN regions and theozyme positions. RFdiffusion, LigandMPNN, and Boltz2 jobs are submitted and tracked by the Model Manager, which registers each output—backbone PDB, FASTA sequence, predicted CIF—as an entity in ProtOS. The evaluation scripts query these entities to compute theozyme RMSD and pLDDT statistics across all candidates.

RFdiffusion2, LigandMPNN, and Boltz2 are published tools; the Baker lab's protocol for combining them is established. ProtOS reproduces that protocol within a managed data framework, applied to a different and more constrained problem. Three AI models, each with its own input format and output format, are integrated into a single pipeline where intermediate results are traceable and the full provenance from PDB entry to predicted structure is recorded. When a candidate fails at validation, we can trace back to its placement, its backbone, its sequence, and ask why.

No experimental validation of the rhodozyme concept exists at this time. Protein design is not an exact science—it requires generating large numbers of candidates, screening them computationally, and selecting the most promising for synthesis and assay. The designs shown here are ongoing work. The contribution of this chapter is not a validated enzyme but a demonstration of what ProtOS can do: integrate multiple structure-generation models into a single pipeline with consistent data management, enabling the large screens from which candidates emerge. We show one such candidate. The rhodozyme concept itself remains speculative; the ability to explore it at scale, using the Baker lab's protocol within a managed data framework, is what this chapter highlights.

The rhodopsin scaffold offers one further possibility. Because the retinal binding pocket is separate from the intracellular catalytic face, different rhodopsins—with different absorption maxima—could carry different enzymes. In principle, a trypsin-rhodozyme built on a 500 nm rhodopsin would absorb at 500 nm; a papain-rhodozyme on a 550 nm scaffold would absorb at 550 nm. Different wavelengths could activate different enzymes, and sequential catalytic steps on the same substrate could become addressable by color. This is the domain of LAMBDA, the spectral tuning model described in the next chapter.

# Thesis Discussion

The introduction identified three gaps in opsin research. Microbial rhodopsins lacked the standardized positional annotation that enabled systematic GPCR comparison. No method predicted spectral properties across the Type I and Type II divide. And existing tools did not compose into workflows linking sequences, structures, and functional properties while maintaining consistent identity. This thesis addresses these gaps while establishing infrastructure that continues to evolve.

# Contributions

This thesis began with a practical problem: analyzing opsin sequences required too much glue code. Downloading a structure, extracting its sequence, aligning it to a reference, annotating binding pocket positions, computing embeddings, predicting spectral properties—each step required different tools, different formats, different conventions. The scientific question came last, after the data engineering.

The contributions emerged from solving that problem. ProtOS provides infrastructure for managing protein data across computational workflows. MOGRN provides standardized coordinates for microbial rhodopsins, enabling the same positional language that GPCRs have had for decades. LAMBDA uses those coordinates to predict spectral properties across opsin families—the first model to treat color tuning as a problem of binding pocket structure rather than sequence features specific to one family.

These components integrate into a system. ProtOS manages the data that MOGRN annotates and LAMBDA consumes. The Opsin Atlas—predictions for 47,700 sequences—exists because ProtOS could route that many sequences through GRN annotation and spectral prediction without custom scripting for each step. The Prey Vision workflow exists because ProtOS-MCP exposes the full stack through natural language.

The deeper contribution is showing that standardized coordinates enable cross-family learning. Type I and type II opsins diverged billions of years ago and share no sequence similarity outside the binding pocket. Yet they tune the same chromophore using the same physics. LAMBDA learns from both because MOGRN and Ballesteros-Weinstein provide a shared vocabulary for binding pocket positions. This principle—that structural annotation enables learning across evolutionary divides—extends beyond opsins.

# Limitations

Each contribution has boundaries that define where future work is needed.

MOGRN works because microbial rhodopsins share a conserved binding pocket despite low sequence identity. Other protein families may lack such obvious structural anchors. Extending generic residue numbering requires identifying appropriate invariants for each family—conserved ligand-binding sites, functional motifs, or structural features that persist across evolutionary distance. The GRNProcessor architecture accommodates new numbering systems; defining them remains manual work.

LAMBDA predicts well within the spectral range where training data is abundant, but struggles with deeply red-shifted opsins. The physics of spectral tuning at long wavelengths makes this a hard problem for any sequence-based approach. More fundamentally, LAMBDA has not been experimentally validated. The Opsin Atlas contains predictions, but predictions are not measurements. Until candidates from the atlas are characterized in the lab, LAMBDA remains a computational tool rather than a tested framework for understanding color tuning.

ProtOS-MCP depends on language model reasoning. The natural language interface works for straightforward queries, but complex multi-step analyses can fail when the model misunderstands the request. As language models improve, this limitation will ease. The workflow dataset provides a benchmark for tracking progress.

Infrastructure enables analysis but does not replace expertise. A researcher using ProtOS-MCP can run analyses they could not run before, but interpreting results still requires understanding the biology. The tools reduce barriers; they do not eliminate the need for scientific judgment.

# Future Work

ProtOS is evolving toward a broader vision: making advanced AI models accessible to every researcher through natural language.

The field of protein design is transforming. Tools like RFdiffusion2, Boltzgen, and LigandMPNN enable computational design of proteins with specified structures and functions. Structure prediction through AlphaFold and Boltz has become routine. Sequence embeddings from protein language models capture evolutionary and functional relationships that sequence identity misses. These capabilities exist, but accessing them requires programming skills that most experimental researchers do not have.

ProtOS-MCP bridges this gap. A crystallographer who has never written code can ask questions in natural language and receive results from state-of-the-art models. The vision is not replacing computational expertise but democratizing access to computational tools. When running an analysis requires only describing what you want, the barrier to exploration drops. Researchers can ask questions they would not have asked if asking required learning Python first.

For optogenetics specifically, the opportunity is significant. The Opsin Atlas identifies candidates with predicted spectral properties—red-shifted opsins for tissue penetration, spectrally separated pairs for multiplexed control, sequences with large spectral separation between dark and activated states for bistable tools. These are computational predictions awaiting experimental validation. Characterizing even a fraction of these candidates would test whether LAMBDA captures real spectral tuning principles and provide new tools for the optogenetics community.

Beyond screening natural variants, the tools developed here could support rational design of light-activated proteins. The rhodozyme concept explored in Chapter 5—coupling light activation to enzymatic function—represents one direction. More broadly, understanding how binding pocket structure determines spectral properties opens the possibility of designing retinal-binding pockets with specified absorption characteristics. This remains speculative, but the foundation is in place: standardized coordinates for describing binding pockets, a model that learns how pocket composition affects spectral tuning, and infrastructure for managing the data that design workflows require.

LAMBDA itself would benefit most from experimental validation. Predictions are useful for prioritizing candidates, but they do not advance understanding until tested. Validating atlas predictions would move LAMBDA from a screening tool to a framework for understanding color tuning in retinal-binding proteins. That transition—from prediction to understanding—is the next step.

# Conclusion

This PhD started with opsins and ended with infrastructure. The original goal was understanding spectral tuning—how the same retinal chromophore produces absorption from UV to far-red depending on its protein environment. That goal required standardized coordinates for microbial rhodopsins, which did not exist. It required linking sequences to structures to spectral measurements, which meant building data management tools. It required making predictions accessible to researchers who might use them, which meant building natural language interfaces.

MOGRN, LAMBDA, and ProtOS emerged from these requirements. Each solves a specific problem: positional annotation, spectral prediction, data management. Together they form a system where biological questions can flow from sequence to prediction without the glue code that usually intervenes.

The work is not finished. ProtOS continues to evolve toward integrating the latest AI models for protein research. LAMBDA awaits experimental validation that would transform it from a prediction tool to a tested framework for understanding color tuning. The natural language interface continues to develop as language models improve in reasoning capability.

What excites me most is the potential for optogenetic engineering. The Opsin Atlas contains tens of thousands of predictions—candidates for tissue-penetrating red-shifted tools, for multiplexed control with spectrally separated opsins, for bistable tools with large dark-to-activated separation. These are hypotheses waiting to be tested. Beyond natural opsins, the binding domain approach could support design of light-activated proteins with tailored spectral properties. The rhodozyme concept hints at coupling light activation to arbitrary protein functions.

Infrastructure shapes what research gets done. When analysis requires programming, only programmers analyze. When analysis requires natural language, anyone can ask. The tools presented here lower barriers. The science that emerges from broader access remains to be seen.

# Supplementary

This supplementary chapter provides reference materials supporting the main text. These include the complete ProtOS-MCP tool catalog, workflow benchmark specifications, and technical details that would interrupt the narrative flow of earlier chapters.

### ProtOS-MCP Tool Catalog

The following table lists all tools exposed through the ProtOS-MCP interface, organized by functional category. Each tool includes its name, required parameters, and a brief description of its function.

#### Loader Tools

| Tool Name | Parameters | Description |
|-----------|------------|-------------|
| `structure_load` | `pdb_id` | Download and load a structure from PDB or AlphaFold |
| `structure_load_batch` | `pdb_ids`, `dataset_name` | Download multiple structures and optionally create a dataset |
| `structure_load_file` | `file_path` | Load a structure from a local mmCIF or PDB file |
| `sequence_load` | `uniprot_id` | Retrieve a sequence from UniProt |
| `sequence_load_batch` | `uniprot_ids`, `dataset_name` | Retrieve multiple sequences and optionally create a dataset |
| `sequence_load_file` | `file_path` | Load sequences from a local FASTA file |
| `molecule_load_from_ccd` | `ccd_id` | Load a molecule definition from the Chemical Component Dictionary |
| `molecule_load_from_structure` | `structure_id`, `ligand_id` | Extract a ligand from a loaded structure |
| `molecule_load_file` | `file_path` | Load a molecule from a local SDF or MOL2 file |

#### Entity and Dataset Tools

| Tool Name | Parameters | Description |
|-----------|------------|-------------|
| `entity_register` | `entity_id`, `entity_type`, `data` | Register a new entity in the registry |
| `entity_lookup` | `query` | Find an entity by ID, alias, or partial match |
| `entity_list` | `entity_type` | List all registered entities of a given type |
| `entity_get_metadata` | `entity_id` | Retrieve metadata for a registered entity |
| `dataset_create` | `name`, `entity_ids`, `metadata` | Create a new dataset from existing entities |
| `dataset_add` | `dataset_name`, `entity_ids` | Add entities to an existing dataset |
| `dataset_remove` | `dataset_name`, `entity_ids` | Remove entities from a dataset |
| `dataset_list` | `entity_type` | List all datasets of a given type |
| `dataset_get_members` | `dataset_name` | Get all entity IDs in a dataset |

#### Structure Tools

| Tool Name | Parameters | Description |
|-----------|------------|-------------|
| `structure_summarize` | `structure_id` | Summarize structure contents (chains, residues, atoms) |
| `structure_summarize_ligands` | `structure_id`, `min_atoms` | Identify and describe all ligands in a structure |
| `structure_get_chains` | `structure_id` | List all chains with their types and lengths |
| `structure_get_sequence` | `structure_id`, `chain_id` | Extract amino acid sequence from a structure chain |
| `structure_get_metadata` | `structure_id` | Get experimental metadata (resolution, method, date) |
| `structure_filter_chain` | `structure_id`, `chain_ids`, `new_id` | Create a new structure containing only specified chains |
| `structure_filter_residues` | `structure_id`, `start`, `end`, `chain_id` | Extract a residue range from a structure |
| `structure_remove_hetatm` | `structure_id` | Remove all non-protein atoms (ligands, water) |
| `structure_find_contacts` | `structure_id`, `chain_a`, `chain_b`, `distance` | Find residue contacts within a distance threshold |
| `structure_align` | `mobile_id`, `reference_id`, `method`, `selection` | Superimpose one structure onto another |
| `structure_calculate_rmsd` | `structure_id_1`, `structure_id_2`, `selection` | Calculate RMSD between two structures |
| `structure_compute_water_networks` | `structure_ids`, `distance` | Analyze water-mediated interactions |
| `structure_export` | `structure_id`, `format`, `output_path` | Export structure to mmCIF or PDB file |

#### Sequence Tools

| Tool Name | Parameters | Description |
|-----------|------------|-------------|
| `sequence_get` | `sequence_id` | Retrieve a loaded sequence |
| `sequence_extract_from_structure` | `structure_id`, `chain_id` | Extract sequence from a structure chain |
| `sequence_align_pairwise` | `query_id`, `target_id`, `method` | Align two sequences (global or local) |
| `sequence_align_multiple` | `dataset_name`, `method` | Perform multiple sequence alignment |
| `sequence_build_database` | `dataset_name` | Build MMseqs2 database from a sequence dataset |
| `sequence_search` | `query`, `database`, `evalue`, `min_identity` | Search for homologs in a sequence database |
| `sequence_cluster` | `dataset_name`, `identity_threshold` | Cluster sequences by identity |
| `sequence_export` | `sequence_id`, `output_path` | Export sequence to FASTA file |
| `sequence_export_dataset` | `dataset_name`, `output_path` | Export dataset to multi-sequence FASTA |

#### GRN Tools

| Tool Name | Parameters | Description |
|-----------|------------|-------------|
| `grn_load_reference` | `family`, `config_path` | Load a GRN reference system for a protein family |
| `grn_list_references` | - | List all available GRN reference systems |
| `grn_annotate_sequence` | `sequence_id`, `reference` | Annotate a sequence with generic residue numbers |
| `grn_annotate_structure` | `structure_id`, `chain_id`, `reference` | Annotate a structure with generic residue numbers |
| `grn_annotate_dataset` | `dataset_name`, `reference`, `output_table` | Annotate all sequences in a dataset |
| `grn_get_position` | `sequence_id`, `position` | Get the residue at a specific GRN position |
| `grn_get_segment` | `sequence_id`, `start`, `end` | Extract a segment by GRN range |
| `grn_compare_positions` | `position`, `sequence_ids` | Compare residues at a position across sequences |
| `grn_get_annotation_table` | `table_name` | Load a GRN annotation table |

#### Embedding Tools

| Tool Name | Parameters | Description |
|-----------|------------|-------------|
| `embedding_list_models` | - | List available embedding models |
| `embedding_ensure_model` | `model_name` | Ensure an embedding model is available |
| `embedding_compute` | `sequence_ids`, `model_name` | Compute embeddings for sequences |
| `embedding_compute_dataset` | `dataset_name`, `model_name` | Compute embeddings for all sequences in a dataset |
| `embedding_load` | `entity_ids` | Load previously computed embeddings |
| `embedding_similarity` | `entity_id_1`, `entity_id_2` | Calculate cosine similarity between two embeddings |
| `embedding_similarity_matrix` | `dataset_name` | Compute pairwise similarity matrix for a dataset |
| `embedding_reduce_dimensions` | `dataset_name`, `method`, `n_components` | Apply dimensionality reduction (PCA, UMAP, t-SNE) |
| `embedding_cluster` | `dataset_name`, `method`, `n_clusters` | Cluster embeddings |

#### Property Tools

| Tool Name | Parameters | Description |
|-----------|------------|-------------|
| `property_create_table` | `name`, `data` | Create a new property table |
| `property_load_table` | `table_name` | Load an existing property table |
| `property_add_rows` | `table_name`, `data` | Add rows to a property table |
| `property_update` | `table_name`, `entity_id`, `column`, `value` | Update a value in a property table |
| `property_query` | `table_name`, `filters` | Query a property table with filters |
| `property_correlate` | `table_name`, `column_x`, `column_y` | Calculate correlation between columns |
| `property_export` | `table_name`, `output_path` | Export property table to CSV or parquet |

#### Model Tools

| Tool Name | Parameters | Description |
|-----------|------------|-------------|
| `model_list` | - | List available computational models |
| `model_ensure` | `model_name` | Ensure a model is available (download if needed) |
| `model_predict_structure` | `sequence_id`, `model` | Predict structure from sequence |
| `model_embed_sequence` | `sequence_id`, `model` | Generate embedding using specified model |
| `model_predict_properties` | `sequence_id`, `properties` | Predict sequence properties |

### Workflow Benchmark Specifications

The ProtOS-MCP workflow collection is stratified by difficulty based on the number of tools required, the number of processors involved, and the complexity of result interpretation.

#### Difficulty Criteria

| Level | Tools | Processors | Interpretation | Example |
|-------|-------|------------|----------------|---------|
| Beginner | 1-2 | 1 | Minimal | Load structure, list ligands |
| Intermediate | 3-5 | 2-3 | Moderate | Align structures, annotate with GRN |
| Advanced | 5+ | 3+ | Substantial | Multi-stage functional analysis |

#### Benchmark Categories

| Category | Description | Example Questions |
|----------|-------------|-------------------|
| Data Retrieval | Loading and inspecting data | "What chains are in structure 3SN6?" |
| Sequence Analysis | Alignment, search, clustering | "How similar are these two sequences?" |
| Structure Analysis | Alignment, contacts, geometry | "What residues contact the ligand?" |
| Cross-Processor | Combining multiple data types | "Annotate this structure with GRN" |
| Functional Analysis | Property correlation, prediction | "Which positions correlate with function?" |

#### Evaluation Metrics

Each benchmark workflow is evaluated on:

- **Tool Selection**: Did the model choose appropriate tools for the question?
- **Parameter Accuracy**: Were tool parameters correctly specified?
- **Execution Order**: Were tools called in a logical sequence?
- **Result Interpretation**: Did the model correctly interpret and explain results?
- **Biological Accuracy**: Was the final answer scientifically correct?

### Processor Data Storage Conventions

ProtOS organizes data in a standardized directory structure. The following table describes the storage locations for each processor.

| Processor | Base Directory | Subdirectories | File Formats |
|-----------|---------------|----------------|--------------|
| Structure | `structure/` | `mmcif/`, `cache/`, `datasets/` | `.cif`, `.pkl`, `.json` |
| Sequence | `sequence/` | `fasta/entities/`, `fasta/datasets/`, `alignments/`, `databases/` | `.fasta`, `.aln`, `.mmseqs` |
| GRN | `grn/` | `reference/`, `tables/`, `configs/` | `.yaml`, `.parquet` |
| Embedding | `embedding/` | `embeddings/`, `datasets/` | `.npy`, `.h5` |
| Property | `property/` | `tables/`, `datasets/` | `.parquet`, `.csv` |
| Molecule | `molecule/` | `records/`, `datasets/` | `.sdf`, `.mol2`, `.json` |

### Entity Registry Schema

The entity registry maintains consistent identity across all processors. Each entity record contains:

| Field | Type | Description |
|-------|------|-------------|
| `id` | string | Primary identifier (e.g., PDB ID, UniProt accession) |
| `type` | string | Entity type (structure, sequence, molecule, etc.) |
| `aliases` | list | Alternative identifiers and names |
| `created` | datetime | When the entity was first registered |
| `modified` | datetime | When the entity was last modified |
| `source` | string | Origin of the data (PDB, UniProt, local file, etc.) |
| `metadata` | dict | Type-specific metadata (resolution, length, etc.) |
| `datasets` | list | Datasets containing this entity |
| `relations` | dict | Links to related entities (e.g., structure→sequence) |

### GRN Reference System Format

GRN reference systems are defined in YAML configuration files with the following structure:

```yaml
family: microbial_rhodopsin
description: Generic residue numbering for microbial rhodopsins
version: 1.0

anchors:
  helix_1: {position: 50, residue: P, reference_idx: 29}
  helix_2: {position: 50, residue: L, reference_idx: 57}
  helix_3: {position: 50, residue: D, reference_idx: 85}
  helix_4: {position: 50, residue: W, reference_idx: 115}
  helix_5: {position: 50, residue: P, reference_idx: 146}
  helix_6: {position: 50, residue: W, reference_idx: 182}
  helix_7: {position: 50, residue: K, reference_idx: 216}

boundaries:
  helix_1: {start: 1.29, end: 1.59}
  helix_2: {start: 2.38, end: 2.62}
  helix_3: {start: 3.29, end: 3.55}
  helix_4: {start: 4.38, end: 4.60}
  helix_5: {start: 5.36, end: 5.63}
  helix_6: {start: 6.28, end: 6.56}
  helix_7: {start: 7.29, end: 7.55}

reference_sequence: >
  MLELLPTAVEGVSQAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSD
  PDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADW
  LFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFV
  WWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYP
  VVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFGEAEAP
  EPSAGDGAAATSD
```

### Supported External Tools and Services

ProtOS integrates with external tools and services through standardized interfaces:

| Tool | Purpose | Integration Method |
|------|---------|-------------------|
| MMseqs2 | Sequence search and clustering | Command-line wrapper |
| Clustal Omega | Multiple sequence alignment | Command-line wrapper |
| MUSCLE | Multiple sequence alignment | Command-line wrapper |
| MAFFT | Multiple sequence alignment | Command-line wrapper |
| Gemmi | mmCIF parsing | Python library |
| BioPython | Sequence manipulation | Python library |
| ESM-2 | Protein embeddings | Docker container |
| Ankh | Protein embeddings | Docker container |
| Boltz2 | Structure prediction | Docker container |
| AlphaFold | Structure prediction | API / Docker |

### File Format Specifications

#### Structure Cache Format

Parsed structures are cached as pickle files containing pandas DataFrames with the following columns:

| Column | Type | Description |
|--------|------|-------------|
| `atom_id` | int | Unique atom identifier |
| `atom_name` | str | Atom name (CA, CB, N, O, etc.) |
| `element` | str | Element symbol |
| `res_name` | str | Three-letter residue code |
| `res_id` | int | Residue sequence number |
| `auth_chain_id` | str | Author-assigned chain ID |
| `x`, `y`, `z` | float | Cartesian coordinates (Å) |
| `occupancy` | float | Occupancy factor |
| `b_factor` | float | Temperature factor |
| `group` | str | ATOM or HETATM |
| `grn` | str | Generic residue number (if annotated) |

#### Embedding Storage Format

Embeddings are stored as HDF5 files with the following structure:

```
embeddings.h5
├── metadata/
│   ├── model_name        # str: e.g., "ankh-large", "esm2-650m"
│   ├── embedding_dim     # int: dimensionality of embeddings
│   ├── created           # str: ISO timestamp
│   └── version           # str: model version
├── entity_ids/           # dataset of entity identifiers
└── embeddings/           # float32 array [n_entities, embedding_dim]
```

For per-residue embeddings (when stored):

```
per_residue_embeddings.h5
├── metadata/
│   └── ...               # same as above
├── {entity_id}/
│   ├── sequence_length   # int
│   └── embeddings        # float32 array [seq_length, embedding_dim]
```

#### Property Table Format

Property tables are stored as Parquet files with the following schema:

| Column | Type | Description |
|--------|------|-------------|
| `entity_id` | string | Primary key linking to registry |
| `*` | varies | User-defined property columns |

Metadata is stored in the Parquet file's schema metadata field as JSON:

```json
{
  "table_name": "opsin_properties",
  "created": "2026-01-15T10:30:00Z",
  "columns": {
    "lambda_max": {"unit": "nm", "source": "experimental"},
    "ion_type": {"values": ["H+", "Cl-", "Na+"], "source": "annotation"}
  }
}
```

#### Alignment File Format

Multiple sequence alignments are stored in FASTA format with gap characters:

```
>entity_id_1
MLELLPTAVEGVSQAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSD...
>entity_id_2
MLELLP-AVEGVSQAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSD...
```

Alignment metadata is stored in a companion JSON file:

```json
{
  "alignment_name": "rhodopsin_msa",
  "method": "clustalo",
  "n_sequences": 150,
  "alignment_length": 280,
  "created": "2026-01-15T10:30:00Z",
  "entity_ids": ["entity_id_1", "entity_id_2", ...]
}
```

