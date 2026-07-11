# Chapter 2: ProtOS — Story and Setup

## The argument

ProtOS is a data framework for multi-modal protein research. It manages sequences, structures, contact graphs, standardized position annotations, embeddings, and property tables through a unified entity system. The chapter presents ProtOS not as a software manual but as a scientific tool, demonstrated through a connected analysis of Type II opsins that builds progressively across sections.

By the end of the chapter, the reader has seen every data modality that LAMBDA requires — without LAMBDA being explained. When Chapter 4 introduces the model, the reader already understands every input. Equally important: the reader understands that the processors compose freely. LAMBDA is one workflow through ProtOS; a different project could start from any modality and move in any direction.

## Narrative structure

Each processor section follows the same pattern:

1. **General concept** — What is this data type? Why does it matter for protein research in general? (Brief for sequence/structure/property; deeper for GRN and embeddings, which are less familiar to structural biologists.)

2. **Minimal ProtOS example** — Working code showing how the processor handles this modality. Not a demo — the interface a researcher actually uses.

3. **Opsin application** — The same operation applied to opsin data, producing a concrete result (figure, table, observation). This is the bulk of each section. Critically: it references outputs from previous sections, so the analysis accumulates.

The balance: a + b ≤ c. The opsin application always dominates.

## Section flow and continuity

```
Sequence  →  "27,640 opsin sequences across 13 subfamilies. Identity distributions
              reveal evolutionary divergence."
    ↓
Structure →  "bR and bovine rhodopsin bind retinal in convergent pockets. Sequences
              and structures now live on the same entity."
    ↓
Graph     →  "Pocket graphs: 46 vs 44 nodes, ~15.5% density. Comparable topology —
              and the graph is already linked to the entity's sequence and structure."
    ↓
GRN       →  "18 positions across 9 subfamilies: conserved anchors vs variable tuning
              sites. SWS1 has G at 2.50 — loss of sodium coordination. Graph nodes
              now carry positional identity."
    ↓
Embedding →  "27,639 Ankh Large embeddings. UMAP: subfamilies cluster without
              supervision. Per-residue features complement the sparse GRN positions."
    ↓
Property  →  "27,640 entities with subfamily, identity, length. ~1,800 have measured
              λmax. All modalities linked — query any combination."
```

Each arrow is a transition: the next section adds a complementary perspective, and ProtOS integrates it with everything that came before. The LAMBDA workflow is one path through this space of processors — a different project could start from a graph and work toward a sequence (e.g., RFdiffusion), or from an embedding and work toward a structure. The processors compose freely.

## What each section is NOT

- Not a feature list ("the Sequence Processor supports BLAST, MMseqs2, batch downloads...")
- Not a software tutorial ("first install, then configure, then run...")
- Not explaining LAMBDA ("building LAMBDA required...")
- Not telling the reader what they already know ("protein sequences are stored in FASTA format")

## What each section IS

- A general concept relevant to any structural biologist
- A minimal code example showing ProtOS's interface
- An opsin-specific application that produces a scientific finding and connects to the next section

## Proportions

For sequence, structure, graph, property: the general concept is 1–2 sentences.
For GRN: full paragraph (numbering problem, x.50 convention, why it matters for cross-family comparison).
For embeddings: full paragraph (PLMs, what embeddings encode, why they complement GRN).
In all sections: opsin application (part c) longer than parts a + b combined.

## Opening framing

1. Chapter 1 gaps → this chapter addresses the third (no unified framework linking sequence, structure, and function)
2. ProtOS manages multiple data modalities through processors sharing a persistent entity system
3. This chapter introduces each processor through a connected opsin analysis — each section adds a modality, building on the previous

No separate architecture sections. No "Unified Data Access" or "Entity/Registry" sections. The entity concept is introduced in the opening paragraph and demonstrated through use in every section that follows.

## Sections

1. **ProtOS** (opening) — Framework purpose, entity system concept, processor overview. 1 paragraph.
2. **Sequence Processor** — Atlas construction, identity distributions (Fig 2.2). Transition: sequences are now registered entities; adding structural data to the same entities is one function call.
3. **Structure Processor** — Type I/II pocket comparison, bR vs bovine rhodopsin (Fig 2.3). Transition: the aligned binding pockets suggest comparable topology — graphs make this explicit.
4. **Graph Processor** — Pocket graph topology comparison (Fig 2.4). Transition: graphs give comparable topology; combining them with a cross-family numbering scheme gives every node a transferable label.
5. **GRN Processor** — 9 subfamilies × 18 positions, microswitches and spectral tuning (Fig 2.5a, 2.5b). Transition: GRN positions are sparse and curated; embeddings provide dense per-residue features that complement them.
6. **Embedding Processor** — Atlas embeddings, UMAP clustering (Fig 2.6). Transition: every modality introduced so far is now associated with the same entities — the property table makes this queryable.
7. **Property Processor** — Atlas properties, λmax coverage, closing integration. Every entity carries sequence, structure (where available), graph, GRN, embedding, and measured properties. Any combination can be queried or fed to a model.
8. **Model Manager** (brief) — ProtOS wraps external tools (Boltz2, RFdiffusion, LigandMPNN, LAMBDA) behind a unified interface. No opsin example needed — this is infrastructure for Chapters 4 and 5.
9. **Limitations and Accessibility** — Honest scope. What ProtOS does not do.

## Closing

After Property Processor: the atlas sequences now carry GRN annotations, per-residue embeddings, binding pocket graphs, and measured spectra — all linked through ProtOS entities. Each processor generalizes beyond opsins to any protein family. The chapters that follow build on this multi-modal representation.

## Tone

Scientific prose. Dense but clear. No adjectives that don't carry information. Every sentence advances either the general concept, the ProtOS interface, or the opsin analysis. The reader should finish the chapter understanding what ProtOS does and what the opsin atlas contains — without having been told either directly.
