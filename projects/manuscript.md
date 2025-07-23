This is a manuscript I am working on:











\subsection{Datasets and Data Preprocessing}







To develop LAMBDA's unified predictive framework, we assembled comprehensive training data spanning both major opsin families—Type I and Type II—which together represent nature's diversity in retinal-based light sensing. Type I opsins contain all-trans retinal in their ground state and isomerize to 13-cis upon activation, while Type II opsins contain 11-cis retinal in darkness and isomerize to all-trans when illuminated. The $\lambda_{\text{max}}$ values in available databases typically report these dark-adapted states—all-trans for Type I and 11-cis for Type II—making direct spectral comparisons between opsin families impossible without accounting for retinal conformation. To address this challenge, LAMBDA predicts $\lambda_{\text{max}}$ for all four possible retinal states—11-cis, 13-cis, all-trans, and deprotonated—from a single protein sequence. This multi-output approach requires training data that captures how diverse protein environments modulate retinal absorption across all conformations, even hypothetical combinations not found in nature. We assembled datasets spanning 350-650 nm with sequences from over 200 species, totaling 1,799 retinal-binding proteins—sufficient diversity for deep neural networks to learn these complex tuning principles. Our primary data sources are the Visual Physiology Opsin Database (VPOD) for Type II opsins and the Inoue dataset for Type I microbial rhodopsins, which together provide complementary coverage of natural opsin diversity.







The Visual Physiology Opsin Database (VPOD), compiled by Frazer et al., provides 864 unique opsin sequences from 166 species, representing the most comprehensive collection of Type II opsins with experimentally validated $\lambda_{\text{max}}$ values. This dataset comprises 721 vertebrate visual opsins and 143 invertebrate opsins, with all $\lambda_{\text{max}}$ measurements taken using 11-cis retinal. The spectral coverage spans 350-611 nm with characteristic clustering around four major visual pigment families: UV-sensitive opsins (350-380 nm), blue-sensitive opsins (420-470 nm), green-sensitive opsins (480-530 nm), and red-sensitive opsins (530-570 nm). Notably, VPOD includes 318 wild-type sequences and 546 experimental mutants—the latter being particularly valuable as they demonstrate how engineered mutations tune absorption wavelengths, directly relevant to optogenetic tool development where researchers modify natural opsins for specific applications.







For Type I opsins, we incorporated the Inoue dataset comprising 884 microbial rhodopsin sequences from bacterial, archaeal and fungal sources, all with $\lambda_{\text{max}}$ measurements using all-trans retinal. These sequences, collected from NCBI non-redundant proteins, environmental metagenomes, and the Tara Oceans expedition, include diverse functional types—bacteriorhodopsins, halorhodopsins, sensory rhodopsins, and channelrhodopsins. In contrast to VPOD's clustered distribution, the Inoue dataset exhibits a continuous spectral distribution across 400-600 nm, reflecting the diverse ecological adaptations of microbial opsins unconstrained by specific visual functions.







Together, these two datasets provide 1,748 opsin sequences spanning both major families, with Type II clustering revealing functional constraints while Type I's continuous distribution fills critical spectral gaps.







To augment this opsin training data and test LAMBDA's ability to generalize beyond seven-transmembrane architectures, we included 51 engineered variants of human cellular retinol binding protein II (hCRBPII). hCRBPII belongs to the lipocalin family with a $\beta$-barrel structure completely distinct from the seven-transmembrane opsin fold, naturally functioning as a retinol (vitamin A) transport protein. Wang et al. engineered these proteins to bind all-trans retinal instead of retinol, creating artificial light-absorbing proteins with no biological function beyond spectroscopy. Remarkably, just nine mutation sites across these 51 variants achieve a 219 nm spectral range (425-644 nm), demonstrating that dramatic spectral tuning is possible even in non-native retinal-binding scaffolds. Including these engineered proteins serves three purposes: testing whether spatial tuning principles extend beyond seven-transmembrane architectures, providing training examples in the far-red region (>600 nm) where natural opsins are sparse, and demonstrating the potential for engineering light sensitivity into arbitrary protein scaffolds.







The complete retinal-binding protein (RBP) dataset thus combines 1,799 sequences: 884 Type I opsins (49.1\%), 864 Type II opsins (48.0\%), and 51 hCRBPII variants (2.9\%), spanning 350-644 nm across three distinct protein folds. Each protein sequence in our dataset has $\lambda_{\text{max}}$ measurements for its native retinal conformation—11-cis for Type II opsins, all-trans for Type I opsins and hCRBPII—yet LAMBDA must predict absorption wavelengths for all four possible retinal states. During training, the model computes loss only for available measurements, learning to generalize from predominantly single-conformation data to predict unmeasured states based on shared tuning principles.







We ensured data quality by removing sequences with ambiguous amino acids (X, B, Z) and excluding entries missing $\lambda_{\text{max}}$ values. We verified the presence of the conserved retinal-binding lysine in all sequences, as this residue forms the essential Schiff base linkage with retinal required for light absorption. To prevent phylogenetic bias while maintaining spectral diversity across training, validation, and test sets, we implemented a two-level stratified splitting strategy. First, we stratified by dataset to maintain the proportions of VPOD, Inoue, and hCRBPII sequences in each split, then stratified by species within each dataset to ensure evolutionary diversity. Species with ten or more sequences received proportional representation across splits, while species with fewer sequences were randomly assigned to ensure all taxonomic groups appeared in the dataset. We allocated 90\% of sequences (1,619) to training, 5\% (90) to validation, and 5\% (90) to test sets—a distribution chosen to provide sufficient training examples for our complex graph neural network architecture while maintaining adequate validation and test sizes. This stratified approach ensures that our model learns generalizable spectral tuning principles rather than overfitting to specific protein families or wavelength ranges, setting the foundation for accurate predictions across diverse retinal-binding proteins.



With 1,799 sequences spanning three distinct protein folds and ranging from 100 to 800 amino acids in length, direct sequence comparison becomes impossible due to different architectures, variable sequence lengths, and varying numbers of transmembrane helices.







\subsection{Generic Residue Numbering (GRN) Assignment}



Generic Residue Numbering (GRN) systems assign standardized labels to amino acids based on their structural position within conserved protein architectures, allowing researchers to refer to functionally equivalent residues across different sequences using a universal coordinate system. In the GRN system, each residue receives a label in the format X.YY, where X denotes the transmembrane helix (1-7) and YY indicates the position within that helix, with position 50 assigned to a functionally important residue chosen as a reference point.







For example, in opsins the retinal-binding lysine that forms the Schiff base is designated 7.43 in Type II opsins, while the primary counterion stabilizing this positive charge is typically found at position 3.28, allowing researchers to discuss these critical residues consistently across hundreds of diverse sequences.







For Type I opsins, we utilized a GRN system developed through comprehensive structural analysis of the microbial rhodopsin landscape, where reference positions were defined based on proximity to the retinal. This system identifies key functional positions across Type I opsins—the Schiff base lysine at 7.50, the primary counterion at 3.57, and conserved waters coordinating at positions like 7.46—creating a common language for discussing microbial rhodopsin structure-function relationships. For Type II opsins, we annotated sequences using the established Ballesteros-Weinstein numbering system, originally developed for class A GPCRs.







Unlike traditional multiple sequence alignments that align residues by sequence similarity, GRN systems function as structure-based alignments where positions are defined by their location in three-dimensional space and functional role, making them robust to insertions and deletions that would disrupt conventional alignments. To assign GRN coordinates to each sequence in our dataset, we used a ProtOS, which also handles GRN assignments in loop and tail regions. This alignment process successfully assigned GRN coordinates to over 98\% of sequences in our dataset. This structural encoding proves essential for our graph-based approach, as residues with the same GRN coordinates occupy equivalent spatial positions relative to retinal across all proteins of a specific type, regardless of their evolutionary origin or sequence identity.











\subsection{Mapping generic residue numberings}



While Type I and Type II opsins have distinct binding pocket architectures with different GRN systems, we can create a unified spatial representation by mapping positions based on their location relative to retinal, guided by the principle that spatial arrangement determines spectral properties. To establish this mapping, we structurally aligned representative Type I and Type II opsin structures by superimposing their retinal molecules, then identified pairs of GRN positions from each system that occupy similar three-dimensional locations within the aligned binding pockets. The spatial mapping process yielded correspondence between positions that occupy similar locations relative to retinal: while the Schiff base lysines at 7.43 (Type II) and 7.50 (Type I) represent an obvious functional match, the mapping also identifies less obvious spatial correspondences such as 3.28 (Type II) with 3.57 (Type I), representing different functional roles but similar spatial positions in the binding pocket.







The complete mapping, shown in Table 1, identifies five key position pairs that form the core of the unified binding pocket representation, though notably many positions in each GRN system have no counterpart in the other due to the fundamentally different architectures.



\newline







\begin{table}[h]



\centering



\caption{Examples of spatial mapping between Type II and Type I opsin GRN positions based on structural alignment relative to retinal}



\begin{tabular}{lll}



\toprule



{} & Type II (animal opsins) & Type I (microbial opsins)\\



\midrule



0 & 7x43 & 7x50 \\



1 & 3x28 & 3x57 \\



2 & 7x39 & 7x46 \\



3 & 54x01 & 4x51 \\



4 & 2x58 & 7x45 \\



\bottomrule



\end{tabular}



\end{table}











For hCRBPII variants, which lack an established GRN system but share an identical fold across all 51 sequences, we used sequential position numbers as proxy GRN coordinates, applying the same conceptual framework of structure-based position labeling.\subsection{Graph Construction}







Spectral tuning emerges from networks of interactions within the binding pocket—electrostatic effects propagate through hydrogen bond networks, charged residues influence neighbors, and the collective environment modulates retinal absorption. Graph representations naturally capture these interaction networks by encoding amino acids as nodes and their spatial relationships as edges, combining the accessibility of sequence-based methods with the mechanistic insights of structure-based approaches.







We construct a unified graph topology through a systematic process. First, we identified all amino acid positions within 7 Ångströms of retinal across multiple high-resolution reference structures spanning Type I opsins, Type II opsins, and hCRBPII. This 7 Å cutoff captures not only direct retinal contacts (3-4 Å) but also second-shell residues that influence spectral tuning through hydrogen bond networks or medium-range electrostatic fields. This comprehensive approach yielded 108 unique GRN positions forming our master graph template: approximately 50 positions mapped between opsin types, plus 58 family-specific positions reflecting architectural differences.







Individual proteins populate only subsets of this topology. Type II opsins typically contain 90-98 nodes (lacking microbial-specific positions), Type I opsins have 50-56 nodes (missing animal-specific loops and segments), and hCRBPII variants contain $\sim$38 nodes (only $\beta$-barrel positions). This variable occupancy allows our single framework to accommodate three distinct protein architectures while preserving their structural relationships.







We defined edges between nodes when residues lie within 4 Ångströms in any reference structure. This threshold captures direct molecular interactions: hydrogen bonds (2.5-3.5 Å), salt bridges (3-4 Å), and tight van der Waals contacts. By taking the union of contacts across all reference structures, we create a superset adjacency that accounts for conformational flexibility—if positions 3.28 and 7.43 contact in any opsin structure, all graphs containing both positions include this edge. This ensures our model considers all plausible interactions within opsin-like binding pockets.







Critically, we exclude retinal and water molecules as explicit nodes. Only amino acid positions are represented, forcing the model to infer chromophore and solvent effects from the protein environment alone. This design choice prevents the model from memorizing static molecular configurations and instead requires learning how amino acid arrangements create the electrostatic fields, pocket geometries, and hydrogen bond networks that stabilize different retinal states. For instance, when polar residues flank the Schiff base, the model must learn that waters likely bridge them, even without explicit water representation.







Each node combines structural identity with biochemical context through a dual encoding strategy. Simple amino acid identity fails to capture functional differences—lysine at position 7.50 (Schiff base) differs fundamentally from surface lysines despite identical chemistry. We employ Ankh, a protein language model trained on 45 million sequences, to generate 1280-dimensional embeddings for each residue. These embeddings encode conservation patterns, secondary structure propensities, and functional motifs learned from evolutionary data. By providing Ankh with complete protein sequences, we capture how each position's biochemical properties emerge from its protein context.







We augment these Ankh embeddings with explicit GRN position information, ensuring the model distinguishes structurally equivalent positions across proteins. This dual representation—biochemical features from Ankh plus structural coordinates from GRN—enables learning both universal principles (how any lysine at position 7.50 forms the Schiff base) and sequence-specific effects (how surrounding residues modulate its properties).







From this unified graph representation, LAMBDA predicts $\lambda_{\text{max}}$ for four retinal states simultaneously: 11-\textit{cis}, 13-\textit{cis}, all-\textit{trans}, and deprotonated. Rather than training separate models or requiring retinal conformation as input, our multi-output design reflects biological reality—the same protein scaffold accommodates multiple retinal configurations, each with distinct spectral properties. The model must discover which amino acid arrangements favor specific conformations, how pocket geometry constrains isomerization, and how electrostatics affect both protonation and spectral shifts.







During training, we compute loss only for experimentally measured states: Type II opsins provide 11-\textit{cis} data, Type I opsins and hCRBPII contribute all-\textit{trans} measurements, while sparse 13-\textit{cis} and deprotonated data exist for select proteins. Despite each protein supervising just one output, the shared graph representation enables cross-conformational learning. Knowledge from Type II opsins (11-\textit{cis}) informs predictions for Type I opsins in hypothetical 11-\textit{cis} states, while patterns from microbial opsins enhance all-\textit{trans} predictions for animal opsins. This multi-task approach encourages the model to learn fundamental tuning principles that transcend specific protein-retinal combinations.







\subsection{Model Architecture and Training}







LAMBDA employs a graph neural network architecture designed to process our heterogeneous protein graphs while capturing the complex, non-local interactions that determine spectral properties. The architecture comprises three main components: an input projection layer, multiple graph attention layers for message passing, and task-specific output heads for each retinal state.







The input projection transforms our node features—1280-dimensional Ankh embeddings concatenated with positional encodings—into a 512-dimensional hidden representation. This dimensionality reduction creates a more computationally tractable representation while preserving essential information. We apply layer normalization and dropout ($p=0.1$) to improve training stability and generalization.







The core of our architecture consists of 8 graph attention network (GAT) layers that iteratively refine node representations through neighborhood aggregation. Each GAT layer employs 8 attention heads with 64-dimensional features per head, allowing the model to attend to different aspects of node relationships simultaneously. The multi-head attention mechanism proves particularly valuable for modeling binding pockets, where a single residue might participate in multiple interaction types—hydrogen bonding through its backbone while contributing electrostatic effects through its sidechain.







Within each GAT layer, attention coefficients are computed as:



\begin{equation}



\alpha_{ij} = \text{softmax}_j(\text{LeakyReLU}(\mathbf{a}^T[\mathbf{W}h_i || \mathbf{W}h_j]))



\end{equation}



where $h_i$ and $h_j$ are node features, $\mathbf{W}$ is a learned transformation, $\mathbf{a}$ is an attention vector, and $||$ denotes concatenation. These coefficients weight the importance of each neighboring node when updating representations:



\begin{equation}



h_i' = \sigma\left(\sum_{j \in \mathcal{N}(i)} \alpha_{ij} \mathbf{W} h_j\right)



\end{equation}







We employ residual connections around each GAT layer, adding the input to the output after message passing. This design allows training deeper networks by providing gradient highways and preserving information from earlier layers. Between GAT layers, we apply batch normalization and dropout ($p=0.2$) to prevent overfitting on our relatively small dataset.







After message passing, we aggregate node features into a graph-level representation using both global mean and max pooling. This dual pooling strategy captures both average binding pocket properties and extreme values that might dominate spectral tuning. The concatenated pooled features pass through a shared fully-connected layer (512 → 256 dimensions) with ReLU activation.







Four task-specific output heads, one per retinal state, transform the shared representation into $\lambda_{\text{max}}$ predictions. Each head consists of two fully-connected layers (256 → 128 → 1) with ReLU activation and dropout ($p=0.3$) between layers. This branched architecture allows specialization for each retinal conformation while benefiting from shared feature learning in earlier layers.







We train LAMBDA using mean squared error loss, computing gradients only for available experimental measurements. The multi-task setup naturally handles missing data—Type II opsins contribute only to 11-\textit{cis} loss, Type I opsins to all-\textit{trans} loss—while the shared architecture enables knowledge transfer across tasks. We employ the AdamW optimizer with an initial learning rate of $10^{-4}$, weight decay of $10^{-5}$, and cosine annealing schedule over 500 epochs.







To prevent overfitting on our 1,619 training proteins, we implement several regularization strategies beyond dropout: early stopping based on validation loss with patience of 50 epochs, gradient clipping at norm 1.0 to prevent training instabilities, and data augmentation through node feature dropout that randomly masks 10\% of Ankh embeddings during training. This last technique improves robustness to missing or uncertain sequence regions.







The trained model processes new sequences through the same pipeline: sequence → GRN assignment → graph construction → predictions for all four retinal states. Inference takes approximately 0.5 seconds per protein on a standard GPU, enabling large-scale screening of genomic databases. By learning from both opsin families simultaneously, LAMBDA captures universal principles of spectral tuning while respecting family-specific architectural constraints, providing a unified framework for understanding and engineering light-sensitive proteins.











for reference this is the introduction:















Light-sensitive membrane proteins that bind retinal are found across all domains of life, mediating diverse photoreceptive processes \cite{terakita2005opsins}. These proteins fall into two classes: type I opsins (microbial opsins) found in bacteria, archaea, and some eukaryotes, and type II opsins found in animals. Despite having no evolutionary relationship, both classes evolved remarkably similar protein folds—seven transmembrane helix bundles that form a binding pocket for retinal, a vitamin A-derived chromophore. This convergent evolution produced the same activation principle: photon absorption causes retinal isomerization, triggering conformational changes in the protein scaffold.







However, the molecular mechanisms differ fundamentally. Type I opsins include ion pumps and channels but also sensory proteins and light-activated enzymes. Type II opsins operate through G protein-coupled receptor signaling, initiating intracellular cascades through effector proteins. Even the orientation of their binding pockets differs relative to the membrane, reflecting their independent evolutionary origins \cite{cite} (Figure~\ref{fig:opsin_pockets}). In opsins, nature has engineered proteins that harness the retinal light switch for diverse functions—from microbial phototaxis and energy generation to animal vision and circadian rhythms \cite{palczewski2006g, wang2009chromophore}.







Type I opsins typically bind all-trans retinal and isomerize to 13-cis upon light absorption. Many are bistable—they can be both photoactivated and photoreversed—and reset through thermal re-isomerization \cite{ernst2014microbial}. In contrast, most type II opsins bind 11-cis retinal that isomerizes to all-trans, requiring enzymatic recycling through the visual cycle for regeneration \cite{palczewski2000crystal}. Notable exceptions like melanopsin show bistable properties, highlighting the mechanistic diversity within type II opsins \cite{cite}.







All opsins bind the same chromophore, retinal, yet exhibit remarkably different properties—absorption spectra spanning UV to red, kinetics from microseconds to minutes, and diverse ion selectivities. The protein environment modulates these chromophore properties through specific amino acid interactions, including electrostatic effects, hydrogen bonding, and steric constraints.







This diversity of light-responsive proteins has made opsins central to optogenetics, where light-activated proteins control biological processes with spatiotemporal precision. Opsins offer distinct advantages: varied ion selectivities enable both neuronal excitation and inhibition, different signaling mechanisms span direct ion flux to GPCR cascades, and their spectral range allows for wavelength-specific control.







To expand the optogenetic toolkit, researchers pursue two strategies: discovering useful variants among the thousands of uncharacterized opsins in genomic databases, and engineering new variants with tailored properties. For discovery, the challenge lies in identifying opsins with desired spectral properties without testing each one experimentally. For engineering, goals include red-shifted variants for deeper tissue penetration, spectrally separated tools for multiplexed control, chimeric receptors, or modified ion selectivities. Both approaches require understanding how sequence and structure of opsins determine their function.







Among these goals, spectral tuning—particularly predicting $\lambda_{\text{max}}$—represents a critical first step. Current methods rely on labor-intensive experimental characterization, scaffolding-design and even random mutagenesis. Accurate sequence-based prediction would enable rapid screening of genomic databases for naturally occurring variants with desired colors, guide targeted engineering efforts, and provide insights applicable to predicting other opsin properties.







To develop such predictive methods, we must first understand the molecular basis of spectral tuning. Opsin spectral properties arise from the interaction between the protonated retinal Schiff base and surrounding amino acids \cite{nathans1987molecular}. The external point charge model explains how charged residues, especially the counterion located $\sim$3--4~\AA{} from the Schiff base, modulate absorption by stabilizing the positive charge \cite{point_charge_ref}. This primary electrostatic effect, combined with hydrogen bonding and steric interactions from nearby residues, determines the wavelength of maximum absorption. Even small changes in this electrostatic environment can shift $\lambda_{\text{max}}$ by tens of nanometers \cite{lin2016interplay, wang2012evolutionary, yokoyama2008evolution, fahmy2000resonance, hart2000visual}.







QM/MM calculations can model these interactions accurately but require high-resolution structures, available for only a few opsins and are extremely compute heavy\cite{katritch2013structure}. Simpler models focusing on key tuning sites work for similar sequences but miss the complex, non-additive effects that produce large spectral shifts \cite{lin2016interplay}. Machine learning approaches have shown promise but were limited to either type I or type II opsins, missing the opportunity to learn universal principles across protein families \cite{cite}.







These limitations prevent exploring the vast diversity of opsins in genomic databases—over 70,000 sequences have been identified across microbial and animal genomes, yet experimental $\lambda_{\text{max}}$ values exist for fewer than 3,000 (only around 4\%) and only for specific retinal conformations. We need methods that predict $\lambda_{\text{max}}$ from sequence alone while capturing the complex electrostatic interactions that determine color.



More importantly, we need approaches that generalize across opsin types, accounting for different retinal conformations (11-\textit{cis} vs all-\textit{trans}) and potentially revealing universal principles of spectral tuning that transcend protein fold boundaries \cite{lapidus2006complete}. By training on both families together, we can leverage their diverse evolutionary solutions to the same problem—spectral tuning—and create a comprehensive resource combining experimental data from both opsin types for multiple retinal conformations.















\input{figures_tex/fig_opsin_pockets}







Here we present LAMBDA, which addresses these challenges by combining three key insights: structural positions can be inferred from sequence through GRN systems, spatial alignment can unify different protein architectures, and graph neural networks can capture the complex interaction networks that determine spectral properties.







First, we leverage Generic Residue Numbering (GRN) systems that both opsin types use to label structurally equivalent positions within their seven-transmembrane architectures \cite{isberg2015gpcrdb}. These systems assign consistent coordinates based on structural landmarks—for example, the retinal-binding lysine is position 7.42 in type II opsins and 7.50 in type I opsins. This provides a structural coordinate framework that can be determined directly from sequence alignment, eliminating the need for experimental structures.







Second, we create a unified spatial representation by aligning the binding pockets of type I and type II opsins using retinal as the common reference point. Through structural superposition, we identify which GRN positions from each opsin type occupy equivalent locations relative to the chromophore (Figure~\ref{fig:opsin_pockets}C). This mapping reveals that positions like 7.43 (type II) and 7.50 (type I) serve the same functional role despite residing in different protein folds. By treating these mapped positions as equivalent coordinates, we enable knowledge transfer between opsin types.







Third, we represent binding pockets as graphs where nodes correspond to GRN positions and edges encode spatial relationships. This graph representation captures how molecular properties propagate through networks of interacting residues—just as retinal's charge distribution affects nearby amino acids, which in turn influence their neighbors. Each node combines its structural coordinate (GRN position) with features from protein language models that capture the evolutionary and biochemical properties of the amino acid at that position \cite{rives2021biological}. This separation of position and identity allows LAMBDA to learn how specific amino acids at particular structural locations influence spectral properties. Through message passing neural networks, the model captures how electrostatic, hydrogen bonding, and steric effects propagate through the binding pocket network to determine $\lambda_{\text{max}}$.







By training on both opsin types simultaneously with their mapped coordinates, LAMBDA learns patterns that generalize across protein families. The model predicts absorption wavelengths for multiple retinal conformations—11-\textit{cis}, all-\textit{trans}, and 13-\textit{cis}—revealing both shared and distinct determinants of spectral tuning. This unified approach enables systematic exploration of genomic databases and rational engineering of opsins with desired spectral properties.







\FloatBarrier







=======================





My Manuscript Writing Instructions: LAMBDA Methods



Core Writing Principles



1. Language Precision







Write declarative "what is" statements







❌ "The GRN system might allow researchers to potentially identify..."



✅ "The GRN system identifies functionally equivalent positions"











Eliminate redundancy ruthlessly







Before writing a sentence, ask: "Have I already said this?"



If explaining a concept twice, consolidate into ONE strong explanation



Track key concepts and ensure each appears exactly once











Compress without losing meaning







❌ "allowing researchers to refer to functionally equivalent residues across different sequences using a universal coordinate system"



✅ "enabling cross-sequence comparison of equivalent positions"















2. Narrative Architecture



Build a clear dependency chain:







Each section must create a need for the next



End sections with transition sentences that create anticipation



Never introduce concepts without their prerequisites







Transition formula:







Last sentence of section A: State the challenge/limitation



First sentence of section B: Present this section as the solution



Example: "With sequences spanning three protein architectures, we need a universal coordinate system. [NEW SECTION] Generic Residue Numbering provides this coordinate system..."







Information hierarchy:







Paragraph 1: Why this matters (motivation)



Paragraph 2: What we did (approach)



Paragraph 3+: How it works (details)



Final paragraph: What this enables (implications)







3. Technical Completeness



The "Algorithm Test":







Could a competent programmer implement this from your description?



If no, add missing steps







Required elements for each method:







Physical/biological justification for choices



Specific numerical parameters with reasoning



Computational details (tools, versions, settings)



Validation metrics or quality checks







Example structure for parameters:



We selected a 7Å distance cutoff based on crystallographic analysis



showing that 95% of spectral tuning effects occur within this radius



of retinal. This threshold captures both direct contacts (3-4Å) and



water-mediated interactions (5-7Å) while excluding bulk solvent effects.



4. Fact Preservation Checklist



Before revising, create a fact inventory:







All numerical values (1,799 sequences, 350-644 nm range, etc.)



Percentages (49.1% Type I, 48.0% Type II, 2.9% hCRBPII)



Thresholds (7Å, 4Å, position 50)



Counts (108 unique positions, 50 mapped positions, 5 examples)







Never alter these during revision.



5. Reader-Centric Organization



The "Naive Reader Test":







Would someone unfamiliar with the field understand why each step matters?



Do concepts appear in dependency order?







Consolidation rules:







If explaining the same tool/concept multiple times, create ONE definitive paragraph



Place this paragraph at first meaningful use



Reference back to it (not re-explain) in subsequent uses







Example of consolidation:



❌ Paragraph 2: "Ankh provides embeddings..."



❌ Paragraph 5: "We use Ankh, which encodes..."



❌ Paragraph 8: "The Ankh model generates..."







✅ Single paragraph: "Node features employ Ankh, a protein language model



trained on 45 million sequences, generating 1280-dimensional embeddings



that capture evolutionary and biochemical context. These embeddings



distinguish identical amino acids in different structural contexts—the



lysine forming the Schiff base differs fundamentally from surface lysines



in its embedding representation."



6. Section-Specific Templates



Dataset Section:



1. State what LAMBDA does (one sentence)



2. Describe primary dataset with numbers



3. Explain secondary dataset and why included



4. State preprocessing steps



5. Conclude with what this enables



Method Sections:



1. Define the problem this method solves



2. State our solution approach



3. Provide implementation details



4. Include validation/quality metrics



5. Connect to next section need



7. Writing Process



First Draft:







Write without editing



Include all technical details



Don't worry about redundancy yet







Revision Pass 1 - Content:







Highlight every repeated concept



Choose the BEST explanation and delete others



Verify all facts remain unchanged







Revision Pass 2 - Flow:







Read only first and last sentences of each paragraph



Do they form a logical story?



Add transition sentences where needed







Revision Pass 3 - Compression:







Find every sentence over 20 words



Can it be split or shortened without losing meaning?



Remove every unnecessary adjective and adverb







Final Check:







Can each section stand alone while contributing to the whole?



Does every sentence add new information?



Would a competent scientist be able to reproduce this work?







8. Forbidden Patterns



Never use:







"might", "could", "potentially" (unless discussing future work)



"allows researchers to" (just state what it does)



"remarkable", "elegant", "powerful" (let the work speak)



Multiple explanations of the same concept



Vague transitions like "To address this challenge"







Always use:







Active voice



Specific numbers and thresholds



Clear causal relationships



One explanation per concept



Transitions that preview content











=================



These are editor 1's notes:



1. Language Precision



Declarative Clarity:



You largely follow clear declarative statements, but occasionally drift into speculative language.



Example:



❌ "demonstrating the potential for engineering light sensitivity into arbitrary protein scaffolds."



✅ "demonstrating that engineered proteins can exhibit dramatic spectral tuning despite non-native retinal-binding scaffolds."



Elimination of Redundancy:



Major Issue: You repeatedly describe the retinal conformations (11-cis, all-trans, etc.). Consolidate these explanations into one paragraph and subsequently refer back to it succinctly.



Example: Merge early paragraphs describing Type I and Type II retinal conformations into a single authoritative statement and reference it consistently.



Compression without Meaning Loss:



The manuscript includes verbose phrases that could be significantly condensed.



Example:



❌ "providing complementary coverage of natural opsin diversity."



✅ "covering natural opsin diversity."



2. Narrative Architecture



Dependency Chain and Transitions:



Moderate Issue: Some sections end abruptly, lacking transitions to the next topic.



Example:



End of GRN assignment section could better set up why mapping across families is required.



Recommended transition: "While GRN systems standardize residues within families, comparing across opsin types requires spatial alignment."



Information Hierarchy:



Moderate Issue: Paragraph structure sometimes deviates from your prescribed hierarchy, mixing motivation, approach, and details.



Example: In "Datasets and Data Preprocessing," your motivations (why datasets matter) and methodological details (numbers, splits, preprocessing) are slightly intermixed. Consider a clear separation:



Motivation: diversity, spectral coverage, limitations of current datasets.



Dataset details (VPOD, Inoue, hCRBPII separately).



Preprocessing specifics (removal criteria, splits).



Implication for LAMBDA.



3. Technical Completeness



Algorithm Test:



Your descriptions of the graph construction and GRN assignment are generally excellent, but a few gaps remain:



You state "ProtOS" is used for GRN assignment without version or parameters. Provide specific software parameters or at least a citation with configuration details.



You say "over 98% of sequences were assigned"—what happens to the remaining ~2%? Clarify how these cases were handled.



Numerical and Physical Justifications:



Excellent detail on distance thresholds (7 Å and 4 Å) and node definitions.



Minor suggestion: briefly justify the choice of a 4 Å contact threshold explicitly by referencing typical bond distances directly in-text rather than implicitly (currently scattered).



Computational Details:



Missing explicit computational details (version of Ankh model, computational resources—GPU type used for timing estimation). Include this explicitly to ensure reproducibility.



Validation and Quality Checks:



You include basic validation (sequence verification steps). However, no explicit quality check on structural alignments (e.g., RMSD threshold) or graph quality checks (e.g., node/edge distribution) is described. Add these briefly.



4. Fact Preservation Checklist



Good Adherence:



You consistently use precise numbers (1,799 sequences, percentage splits, GRN position numbering), and numerical parameters (7 Å and 4 Å cutoffs).



Double-check the consistency of GRN positions—positions such as 7.43 (Type II) and 7.50 (Type I) are critical and repeatedly referenced. Ensure you never accidentally interchange these numbers.



5. Reader-Centric Organization



Naive Reader Test:



Generally good, but complex concepts like GRN numbering and spatial mapping may require clearer first explanations. Consider adding one definitive paragraph explicitly explaining what GRN numbering is, why it is crucial, and its basic rationale, then reference it briefly afterwards.



Example:



✅ "Generic Residue Numbering (GRN) assigns coordinates based on structural landmarks, enabling researchers to compare functionally equivalent residues across protein sequences despite evolutionary divergence." Then briefly refer back: "Using GRN coordinates (as previously defined), we mapped..."



Consolidation of repeated concepts:



You repeatedly explain the Ankh embedding concept. Consolidate fully into one paragraph, explicitly stating embedding dimensionality (1280) and its evolutionary/biochemical relevance clearly at first use, then only briefly reference afterwards.



6. Section-Specific Templates



Datasets Section:



Moderate Issue:



Clearly state upfront (in first paragraph) what LAMBDA specifically requires from the datasets (multiple retinal conformations, diverse spectra, structural variability). Currently, this is scattered.



Clearly separate each dataset description (VPOD, Inoue, hCRBPII). Each should clearly state size, spectral range, and purpose succinctly, one after the other.



Conclude the dataset section clearly summarizing: "These datasets provide the spectral and structural diversity necessary for training LAMBDA to generalize spectral tuning principles across distinct retinal-binding proteins."



Methods Section (GRN assignment, mapping, Graph construction):



GRN assignment is good but lacks explicit validation. Add a short sentence describing how accuracy of GRN assignment was verified.



Mapping section well-described, but explicitly mention whether all positions from each GRN were tested for alignment, and clarify criteria for "key position pairs."



Graph construction is strong but clarify explicitly why water and retinal are excluded as nodes (currently explained, but the justification is scattered). Consolidate justification explicitly: "Water and retinal molecules are excluded explicitly as nodes to prevent model overfitting on static molecular configurations and force inference of dynamic solvent and chromophore interactions."



7. Writing Process Suggestions



Revision Pass Recommendations:



Pass 1 (Content):



Aggressively merge repeated explanations (retinal states, Ankh embeddings, GRN numbering rationale).



Consolidate dataset details more clearly.



Pass 2 (Flow):



Add explicit transitions between sections (end of GRN numbering clearly transitioning to GRN mapping, end of graph section clearly transitioning to Model Architecture).



Pass 3 (Compression):



Shorten verbose phrases.



Example: "By taking the union of contacts across all reference structures..." could simply be "By unifying contacts across reference structures..."



8. Forbidden Patterns (Check)



Language Usage:



You frequently use "allow researchers" or "potentially." Remove or rewrite these instances as explicit declarative statements.



Example rewrite:



❌ "allowing researchers to refer to functionally equivalent residues across different sequences"



✅ "enabling cross-sequence comparison of functionally equivalent residues."



Adjectives/Adverbs:



Remove subjective adjectives: "Remarkably," "particularly valuable," "critical," unless strictly necessary.



Additional Detailed Notes:



Table 1:



Excellent clarity. Consider briefly stating explicitly why these specific pairs were chosen (structural proximity, functional equivalence, etc.).



Computational Performance Metrics:



"Inference takes approximately 0.5 seconds per protein on a standard GPU"—add explicitly which GPU model or capability for clarity.



Validation Split:



Clarify rationale behind your chosen training (90%)/validation/test (5% each) split. Explicitly state why this distribution is adequate given your dataset size and complexity.



Overall Summary:



Your methods section is generally strong, technically detailed, and largely aligns with your writing principles. The key improvements required are:



Aggressive redundancy removal (especially retinal conformations, Ankh embeddings, and GRN explanations).



Clearer narrative transitions between sections (GRN → Mapping → Graph → Model).



Explicit computational parameters and validation checks to improve reproducibility.



Compression of verbose phrasing and removal of speculative language.



Explicit statement of motivations and implications at start and end of each key subsection.



Addressing these points will elevate the precision, readability, and reproducibility of your manuscript substantially.















==========================











please take the authors position, adhereing to the writing principles and review the manuscript.











make minimal changes, to achieve the editor's wishes.







Provide the full text. Do not change logical statements (facts your coauthor has already written), but try to adhere to the editor's notes. Do not invent facts.







strictly adhere to the writing principles and suggestions.







Avoid large changes.







before starting create a detailed plan to make the required edits.