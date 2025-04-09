# Protos Overview

## What is Protos?

Protos is a Python library designed to **standardize and execute complex computational workflows** essential for structural biology research. It provides integrated capabilities for handling diverse biological data types – including sequences, 3D structures, alignments, and associated properties – through a unified framework.

The core function of Protos is to manage multi-step analyses by breaking them down into defined tasks handled by modular components.

## Core Architecture: Processors & Interoperability

Protos utilizes a modular architecture built upon distinct Python components called **'Processors'**. Each Processor is specialized for a specific domain, such as:

*   **`CifProcessor`**: Manages 3D structure data.
*   **`SeqProcessor`**: Handles sequence data and alignments.
*   **`GRNProcessor`**: Implements Generic Residue Numbering systems.
*   **`LigProcessor`**: Deals with ligand information and interactions.
*   **`EMBProcessor`**: Manages protein embeddings.
*   **`PropertyProcessor`**: Integrates metadata and calculated properties.

A key feature is the **interoperability** between these Processors. Outputs from one (e.g., selected residues from `CifProcessor`) can directly serve as inputs for another (e.g., for GRN mapping by `GRNProcessor`, followed by sequence analysis by `SeqProcessor`), enabling flexible construction of sophisticated analysis pipelines.

The relationships and primary data flow between these core processors are visualized below:

![Protos Processor Overview Diagram](resources/protos_overview.png)
*(Diagram showing connections between CP, SP, GRNP, LP, EMBP, and PP)*

---

## Protos-MCP: The Vision & Approach

**1. The Need: Simplifying Complex Bioinformatics Tasks**

Computational analysis is fundamental to modern structural biology. Ideally, researchers should be able to ask complex questions directly:

*   "Fetch PDB 7ZOU and align it to my AlphaFold model X."
*   "Identify residues within 4Å of the ligand in these 5 structures."
*   "How conserved are the residues corresponding to GRN 3.50 and 6.48 across these aligned sequences?"

Currently, answering such questions involves multiple distinct steps: finding and downloading data, converting formats, running specific software tools, parsing results, and integrating information – a process requiring significant bioinformatics effort or programming skills.

**2. Protos: A Functional Backend Lacking Accessibility**

The Protos library contains the programmatic building blocks to perform many of these underlying steps. It wasn't designed abstractly but **emerged organically** from consolidating helper functions and processing scripts developed while performing repetitive, complex workflows across different research projects (e.g., implementing a novel structure-based GRN system for opsins [manuscript in preparation], developing graph neural networks using GRNs and embeddings for property prediction).

While functional and internally useful for streamlining these specific research tasks, Protos requires Python proficiency and familiarity with its specific framework. Experienced bioinformaticians might script these steps themselves for one-off tasks, while bench biologists typically cannot use the library directly. **This accessibility gap is the primary reason Protos remains an unreleased tool.**

**3. The Solution: Leveraging MCP for an AI Interface**

Recent advances offer a path to make Protos broadly useful. We propose **Protos-MCP**, integrating Protos with Large Language Models (LLMs) using the **Model Context Protocol (MCP)**.

*   **What is MCP?** MCP is a **novel, standardized protocol** designed specifically for LLMs (like ChatGPT, Claude) to interact securely and reliably with external software. It allows applications like Protos to define their capabilities as discrete **'Tools'**. The LLM can then understand these Tools and request their execution to accomplish tasks. This standardization is poised to become a **major staple** for enabling AI to leverage specialized scientific software.

*   **The Protos-MCP Workflow:** We will implement an **MCP Server** acting as a controller that wraps the Protos library. Key Protos functions become MCP 'Tools'. Users interact with an LLM in natural language. The LLM plans the required analysis steps and instructs the MCP Server (via MCP) to execute the corresponding Protos Tools. Results are passed back to the LLM for interpretation and presentation to the user. This interaction is visualized below:

![Protos-MCP Interaction Flow](resources/MCP_integration.png)
*(Diagram showing User -> LLM -> MCP Server -> Protos Library -> MCP Server -> LLM -> User)*

*   **Implementation Path & Vision:** This approach leverages the existing, functional Protos backend. Development focuses on creating the MCP Tool interfaces. We will follow a **bottom-up strategy**, initially implementing MCP Tools for robust, simpler Protos functions and progressively adding more complex, multi-step workflow capabilities. The ultimate goal is a system where the broader structural biology community can effectively "talk" to Protos, using its power without needing to code, thereby transforming Protos into a widely accessible and valuable resource.

---

## Protos-MCP Example Workflows

The following examples illustrate the **target functionality** of the proposed Protos-MCP system. They demonstrate how user requests in natural language could be translated by an LLM into a series of calls to Protos capabilities, orchestrated via the MCP server.

For clarity, each example includes:
1.  The user's natural language request.
2.  A conceptual plan the LLM might formulate.
3.  **Illustrative Python code snippets** demonstrating how the task *could currently be approached programmatically* using the underlying Protos library functions. This highlights the existing capabilities and the complexity abstracted by the Protos-MCP layer.
4.  The **TODO**: Defining corresponding MCP Tools that wrap these Protos functions.
5.  The expected user outcome via Protos-MCP.

*(Note: The Python snippets are illustrative. Actual Protos methods and the implementation of robust MCP Tools are part of the proposed work).*

# Protos-MCP Example Workflows

Computational analysis is increasingly vital for structural biology, enabling insights from large datasets. The Protos Python library provides a foundational framework for these analyses, integrating data handling (sequences, structures, alignments) and core processing steps through interoperable modules ('Processors'). Protos allows complex, multi-step workflows to be constructed programmatically. However, this requires Python scripting and familiarity with the library, limiting its direct use by many researchers.

This project proposes **Protos-MCP**, which aims to bridge this gap by integrating Protos with Large Language Models (LLMs) via the Model Context Protocol (MCP). The goal is to expose Protos' capabilities as MCP 'Tools', allowing LLMs to orchestrate complex analyses based on user requests in natural language.

The following examples illustrate how users might interact with the proposed Protos-MCP system. For each example, we show:
1.  The user's natural language request.
2.  A conceptual plan the LLM might formulate.
3.  **Illustrative Python code snippets** showing how the same task could currently be achieved *programmatically* using existing or planned Protos functionalities. This highlights the steps Protos handles and the complexity Protos-MCP aims to abstract.
4.  The implementation goal: wrapping these Protos functions as MCP Tools.
5.  The expected user outcome via Protos-MCP.

*(Note: While Protos provides the foundation, some specific functions or integration points shown in the code snippets might require refinement or further implementation as part of creating robust MCP tools).*

## Example Use Cases

### Example 1: Simple Data Fetching

**Goal:** Retrieve standard data from common bioinformatics databases using natural language.

**User Request (Natural Language):**
> "Get me the PDB structure for 7ZOU and the protein sequence for human Rhodopsin from UniProt (ID P08100)."

**Conceptual LLM Plan:**
1. Identify PDB ID: `7ZOU`.
2. Identify UniProt ID: `P08100`.
3. Plan to execute: (a) Fetch PDB structure, (b) Fetch UniProt sequence.

**Current Protos Programmatic Approach:**
```python
from protos.processing.structure.struct_processor import CifProcessor
from protos.processing.sequence.seq_processor import SeqProcessor
cp = CifProcessor()
sp = SeqProcessor()

# Fetch and load PDB structure
cp.load_structure('7ZOU')
print(f"Structure 7ZOU loaded/downloaded.")

# Fetch UniProt sequence
rhodopsin_seq = sp.fetch_uniprot_sequence('P08100')
print(f"Sequence for P08100 fetched:\n{rhodopsin_seq[:60]}...")

# Alternative approach to load sequences directly into the processor
sp.load_sequences_from_uniprot(['P08100'])
print(f"Sequence for P08100 loaded into SeqProcessor.")

# Results are now available programmatically within cp.data or sp.sequences
```

**Protos-MCP Implementation:**
- MCP Tool `fetch_rcsb_structure` calls the underlying Protos method for fetching PDB files.
- MCP Tool `fetch_uniprot_sequence` calls the underlying Protos method for fetching sequences.
- The LLM orchestrates calling these Tools via the MCP server based on the user request.

**Outcome (via Protos-MCP):**
- The system confirms file download/retrieval and provides the requested sequence or file location info directly to the user.

### Example 2: GRN Mapping & Sequence Conservation

**Goal:** Link structural positions (via GRN) to sequence conservation using natural language.

**User Request (Natural Language):**
> "For structure 6CMO chain A, which residues correspond to GPCRdb positions 3.50, 6.48, and 7.39? Also, how conserved are these positions in a Class A GPCR alignment?"

**Conceptual LLM Plan:**
1. Load structure 6CMO.
2. Get sequence of chain A.
3. Assign GPCRdb GRNs to sequence.
4. Map specific GRNs to 6CMO_A residues.
5. Load relevant Class A GPCR alignment.
6. Calculate conservation for target GRNs in alignment.
7. Present results.

**Current Protos Programmatic Approach:**
```python
from protos.processing.structure.struct_processor import CifProcessor
from protos.processing.grn.grn_processor import GRNProcessor
from protos.processing.sequence.seq_processor import SeqProcessor
cp = CifProcessor()
grnp = GRNProcessor()
sp = SeqProcessor()

# 1 & 2: Load structure and get sequence
cp.load_structure('6CMO')
seq_A = cp.get_sequence('6CMO', 'A')

# 3: Assign GRNs
grn_map_info = grnp.assign_grns_to_sequence(seq_A, pdb_id='6CMO_A', grn_scheme='GPCRdb_A')

# 4: Map specific GRNs
target_grns = ['3.50', '6.48', '7.39']
residue_mapping = grnp.get_residues_from_grns(target_grns, protein_id='6CMO_A')
print(f"Residue Mapping for 6CMO_A: {residue_mapping}")

# 5: Load Alignment
alignment_loaded = sp.load_alignment('ClassA_GPCRs_GRNaligned.fasta')

# 6: Calculate Conservation
conservation_scores = sp.calculate_conservation(target_grns, position_type='GRN')
print(f"Conservation Scores: {conservation_scores}")
```

**Protos-MCP Implementation:**
- MCP Tools for each required Protos method: `load_structure`, `get_sequence`, `assign_grn`, `get_residues_from_grns`, `load_alignment`, `calculate_conservation`.
- The LLM orchestrates the sequence of Tool calls via the MCP server.

**Outcome (via Protos-MCP):**
- The system directly provides the user with the residue mapping and conservation scores.

### Example 3: Multi-Structure Pocket Analysis with Alignment & Ranking

**Goal:** Perform a complex, multi-step comparative analysis via natural language.

**User Request (Natural Language):**
> "Align structures 6CMO, 7ZOU, and 5XEZ. For each, find residues within 4Å of the ligand RET. Map these pocket residues to GPCRdb GRNs. Then, in each aligned structure, calculate the average distance between the C-alpha of residue 3.50 and the C-alphas of all identified binding pocket residues. Rank the structures by this average distance."

**Conceptual LLM Plan:**
1. Load structures.
2. Align structures.
3. Loop through each structure:
   a. Identify pocket residues near RET.
   b. Assign GRNs.
   c. Map pocket residues to GRNs.
   d. Get coordinates for GRN 3.50 and pocket GRNs (using alignment context).
   e. Calculate average distance.
   f. Store result.
4. Rank results.
5. Present ranked list.

**Current Protos Programmatic Approach:**
```python
import numpy as np
from protos.processing.structure.struct_processor import CifProcessor
from protos.processing.grn.grn_processor import GRNProcessor
cp = CifProcessor()
grnp = GRNProcessor()

pdb_ids = ['6CMO', '7ZOU', '5XEZ']
ligand_id = 'RET'
distance_cutoff = 4.0
anchor_grn = '3.50'
grn_scheme = 'GPCRdb_A'
results = {}

# 1. Load structures
for pdb_id in pdb_ids:
    cp.load_structure(pdb_id)
print("Structures loaded.")

# 2. Align structures
alignment_info = cp.align_structures(pdb_ids, method='Cealign')
print("Structures aligned.")

# 3. Loop through structures
for pdb_id in pdb_ids:
    # 3a. Identify pocket residues
    pocket_residues = cp.extract_binding_pocket(pdb_id, ligand=ligand_id, distance=distance_cutoff)
    
    # 3b & 3c. Assign GRNs and Map pocket residues
    seq = cp.get_sequence(pdb_id)
    grn_map_info = grnp.assign_grns_to_sequence(seq, pdb_id=pdb_id, grn_scheme=grn_scheme)
    pocket_grns = grnp.get_grns_from_residues(pocket_residues, protein_id=pdb_id)
    
    # 3d & 3e. Get relevant coordinates (C-alphas, using alignment context)
    coord_350 = cp.get_ca_coordinate(pdb_id, identifier_type='GRN', identifier_value=anchor_grn, alignment_info=alignment_info)
    pocket_coords = cp.get_ca_coordinates_for_list(pdb_id, identifier_type='GRN', identifier_list=pocket_grns, alignment_info=alignment_info)
    
    # 3f & 3g. Calculate average distance
    distances = [np.linalg.norm(coord_350 - p_coord) for p_coord in pocket_coords]
    avg_dist = np.mean(distances)
    results[pdb_id] = avg_dist
    print(f"Processed {pdb_id}, Avg Dist: {avg_dist:.2f} Å")

# 4 & 5. Rank results
ranked_results = sorted(results.items(), key=lambda item: item[1])
print("\nRanked Structures by Average Distance:")
for pdb, dist in ranked_results:
    print(f"- {pdb}: {dist:.2f} Å")
```

**Protos-MCP Implementation:**
- MCP Tools for each complex step: `load_structures`, `align_structures`, `extract_binding_pocket`, `assign_grn`, `get_grns_from_residues`, `get_ca_coordinate`, `get_ca_coordinates_for_list`, and `calculate_average_distance`.
- The LLM manages the loop and orchestrates the sequence of Tool calls via the MCP server.

**Outcome (via Protos-MCP):**
- The system returns the final ranked list of structures based on the calculated average distance, directly answering the user's complex comparative question.

## Installation

This section describes how to install Protos and its dependencies.

### Prerequisites

- **Python**: Protos requires Python 3.9 or later.
- **Package Manager**: pip (or optionally uv) for installing Python packages.

### Docker

The easiest way to run Protos-MCP with all dependencies is using the Docker image. This approach isolates the environment and bundles the core library along with necessary external tools.

```bash
docker run -it --rm \
  -v $(pwd)/my_protos_data:/app/data \
  protos-mcp-image:latest
```

This command runs the Protos-MCP container:
- The `-v $(pwd)/my_protos_data:/app/data` flag mounts a local directory (`my_protos_data` in your current location) into the container's `/app/data` directory.
- Protos (via its ProtosPaths system) reads from and writes to this `/app/data` directory within the container, allowing seamless interaction with your local data files.

### Standard Installation via PyPI

You can install the core Protos Python library using pip:

```bash
pip install protos
```

### External Dependencies

Protos leverages powerful, established bioinformatics tools for specific tasks. These are automatically installed with the Docker image but need to be installed separately if using the PyPI method.

Key external dependencies include:

- **MMseqs2**: Used for fast, sensitive sequence searching and clustering.
  - Installation: Follow the official MMseqs2 installation guide: https://github.com/soedinglab/MMseqs2
- **GTalign**: A GPU-accelerated tool for fast protein structure alignment.
  - Installation: Follow the official GTalign installation guide: https://github.com/BioinfoMachineLearning/GTalign

### Local Development Setup

If you want to contribute to Protos development:

**Clone the Repository:**
```bash
git clone https://github.com/protos-team/protos.git
cd protos
```

**Create a Virtual Environment:**
```bash
python -m venv venv
source venv/bin/activate  # On Windows use `venv\Scripts\activate`
```

**Install Dependencies:**
```bash
pip install -r requirements.txt
pip install -r requirements-dev.txt
```

**Install Protos in Editable Mode:**
```bash
pip install -e .
```

## Configuration

Protos uses a path management system (ProtosPaths) to handle data directories. By default, it creates and uses a data directory relative to where it's run. You can configure these paths via environment variables or programmatically as needed.
