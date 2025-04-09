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

**Current Protos Programmatic Approach (Illustrative Code):**
```python
# Assumes processor instances are initialized
# from protos.processing.structure.struct_processor import CifProcessor
# from protos.processing.sequence.seq_processor import SeqProcessor # Assuming SeqProcessor handles UniProt fetch
# cp = CifProcessor()
# sp = SeqProcessor()

# Fetch and load PDB structure (downloads if not found locally)
try:
    # cp.load_structure('7ZOU') # Method would handle fetch/load logic
    print(f"Structure 7ZOU loaded/downloaded.") # Placeholder confirmation
except FileNotFoundError:
    print(f"Error: Could not find or download PDB 7ZOU.")

# Fetch UniProt sequence
try:
    # Assuming fetch_uniprot_sequence returns the sequence string or adds it internally
    # rhodopsin_seq = sp.fetch_uniprot_sequence('P08100') # Placeholder method
    rhodopsin_seq = "MNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTPHEETNNESFVIYMFVVHFIIPLIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQGSDFGPIFMTIPAFFAKTSAVYNPVIYIMMNKQFRNCMVTTLCCGKNPLGDDEASTTVSKTETSQVAPA" # Example sequence
    if rhodopsin_seq:
        print(f"Sequence for P08100 fetched:\n{rhodopsin_seq[:60]}...")
    # Or if it loads internally:
    # sp.load_sequences_from_uniprot(['P08100'])
    # print(f"Sequence for P08100 loaded into SeqProcessor.")
except Exception as e:
    print(f"Error fetching UniProt sequence P08100: {e}")

# Results are now available programmatically within cp.data or sp.sequences
```

**Proposed Protos-MCP Implementation:**
- Create an MCP Tool (e.g., `fetch_rcsb_structure`) that calls the underlying Protos method for fetching PDB files.
- Create an MCP Tool (e.g., `fetch_uniprot_sequence`) that calls the underlying Protos method for fetching sequences.
- The LLM orchestrates calling these Tools via the MCP server based on the user request.

**Expected Outcome (via Protos-MCP):**
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

**Current Protos Programmatic Approach (Illustrative Code):**
```python
# Assumes processor instances are initialized
# from protos.processing.structure.struct_processor import CifProcessor
# from protos.processing.grn.grn_processor import GRNProcessor
# from protos.processing.sequence.seq_processor import SeqProcessor
# cp = CifProcessor()
# grnp = GRNProcessor()
# sp = SeqProcessor()

try:
    # 1 & 2: Load structure and get sequence
    # cp.load_structure('6CMO') # Assumes structure file exists or is fetched
    # seq_A = cp.get_sequence('6CMO', 'A') # Assumes method extracts sequence from loaded structure data

    # Placeholder sequence for example continuity
    seq_A = "MNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTPHEETNNESFVIYMFVVHFIIPLIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQGSDFGPIFMTIPAFFAKTSAVYNPVIYIMMNKQFRNCMVTTLCCGKNPLGDDEASTTVSKTETSQVAPA"

    # 3: Assign GRNs
    # Assumes assign_grns_to_sequence returns a mapping or updates grnp instance
    # grn_map_info = grnp.assign_grns_to_sequence(seq_A, pdb_id='6CMO_A', grn_scheme='GPCRdb_A') # Actual GRN assignment logic

    # 4: Map specific GRNs (using placeholder results)
    target_grns = ['3.50', '6.48', '7.39']
    # residue_mapping = grnp.get_residues_from_grns(target_grns, protein_id='6CMO_A') # Requires method to use stored map
    residue_mapping = {'3.50': 'K296', '6.48': 'W265', '7.39': 'Y306'} # Example result
    print(f"Residue Mapping for 6CMO_A: {residue_mapping}")

    # 5: Load Alignment (assuming pre-computed and GRN-aligned)
    # Depending on implementation, this might be via SeqProcessor or GRNProcessor
    # alignment_loaded = sp.load_alignment('ClassA_GPCRs_GRNaligned.fasta') # Assumes method exists
    # or grnp.load_grn_table('ClassA_GPCRs_GRNaligned.csv')
    alignment_loaded = True # Placeholder status

    # 6: Calculate Conservation (using placeholder results)
    if alignment_loaded:
        # conservation_scores = sp.calculate_conservation(target_grns, position_type='GRN') # Assumes method exists
        conservation_scores = {'3.50': 0.98, '6.48': 0.85, '7.39': 0.70} # Example result
        print(f"Conservation Scores: {conservation_scores}")
    else:
        print("Alignment could not be loaded.")

except Exception as e:
    print(f"An error occurred: {e}")
```

**Proposed Protos-MCP Implementation:**
- Create MCP Tools corresponding to each required Protos method (e.g., `load_structure`, `get_sequence`, `assign_grn`, `get_residues_from_grns`, `load_alignment`, `calculate_conservation`).
- The LLM orchestrates the sequence of Tool calls via the MCP server.

**Expected Outcome (via Protos-MCP):**
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

**Current Protos Programmatic Approach (Illustrative Code):**
```python
import numpy as np
# Assumes processor instances are initialized
# from protos.processing.structure.struct_processor import CifProcessor
# from protos.processing.grn.grn_processor import GRNProcessor
# Assume alignment functions/wrappers exist within Protos (e.g., in CifProcessor or separate module)
# cp = CifProcessor()
# grnp = GRNProcessor()

pdb_ids = ['6CMO', '7ZOU', '5XEZ']
ligand_id = 'RET'
distance_cutoff = 4.0
anchor_grn = '3.50'
grn_scheme = 'GPCRdb_A'
results = {}

# Placeholder data structures and functions for demonstration
# Real implementation would use pandas DataFrames, BioPython structures, etc.
mock_structures = {
    '6CMO': {'seq': 'SEQ_A_6CMO...', 'atoms': {'A:K296': [1.0, 2.0, 3.0], 'A:W265': [4.0, 5.0, 6.0], 'A:Y306': [7.0, 8.0, 9.0], 'A:L111': [10.0, 11.0, 12.0]}}, # Example coords
    '7ZOU': {'seq': 'SEQ_A_7ZOU...', 'atoms': {'A:K296': [1.1, 2.1, 3.1], 'A:W265': [4.1, 5.1, 6.1], 'A:Y306': [7.1, 8.1, 9.1], 'A:M111': [10.1, 11.1, 12.1]}},
    '5XEZ': {'seq': 'SEQ_A_5XEZ...', 'atoms': {'A:K296': [1.2, 2.2, 3.2], 'A:W265': [4.2, 5.2, 6.2], 'A:Y306': [7.2, 8.2, 9.2], 'A:I111': [10.2, 11.2, 12.2]}}
}
mock_grn_maps = { # ProteinID -> {GRN -> ResID}
    '6CMO': {'3.50': 'A:K296', '6.48': 'A:W265', '7.39': 'A:Y306', '3.33': 'A:L111'},
    '7ZOU': {'3.50': 'A:K296', '6.48': 'A:W265', '7.39': 'A:Y306', '3.33': 'A:M111'},
    '5XEZ': {'3.50': 'A:K296', '6.48': 'A:W265', '7.39': 'A:Y306', '3.33': 'A:I111'}
}
mock_pockets = { # ProteinID -> List[ResID]
    '6CMO': ['A:K296', 'A:W265', 'A:L111'],
    '7ZOU': ['A:K296', 'A:W265', 'A:M111'],
    '5XEZ': ['A:K296', 'A:W265', 'A:I111']
}

def get_mock_coord(pdb_id, res_id):
    # Placeholder: Simulates getting CA coord; needs alignment context in real version
    return np.array(mock_structures[pdb_id]['atoms'].get(res_id, [0,0,0]))

try:
    # 1. Load structures (Simulated: data assumed loaded into mock_structures)
    print("Structures loaded (simulated).")

    # 2. Align structures (Simulated: assume alignment done, coords are comparable)
    # alignment_info = cp.align_structures(pdb_ids, method='Cealign') # Actual call
    print("Structures aligned (simulated).")
    alignment_info = "mock_alignment" # Placeholder

    # 3. Loop through structures
    for pdb_id in pdb_ids:
        # 3a. Identify pocket residues (Simulated)
        # pocket_residues = cp.extract_binding_pocket(pdb_id, ligand=ligand_id, distance=distance_cutoff)
        pocket_residues = mock_pockets[pdb_id]

        # 3b & 3c. Assign GRNs and Map pocket residues (Simulated)
        # seq = cp.get_sequence(pdb_id)
        # grn_map_info = grnp.assign_grns_to_sequence(seq, pdb_id=pdb_id, grn_scheme=grn_scheme)
        # pocket_grns = grnp.get_grns_from_residues(pocket_residues, protein_id=pdb_id)
        grn_map = mock_grn_maps[pdb_id]
        pocket_grns = [grn for grn, res_id in grn_map.items() if res_id in pocket_residues]

        # 3d & 3e. Get relevant coordinates (C-alphas, using alignment context) (Simulated)
        # coord_350 = cp.get_ca_coordinate(pdb_id, identifier_type='GRN', identifier_value=anchor_grn, alignment_info=alignment_info)
        # pocket_coords = cp.get_ca_coordinates_for_list(pdb_id, identifier_type='GRN', identifier_list=pocket_grns, alignment_info=alignment_info)
        res_id_350 = grn_map.get(anchor_grn)
        coord_350 = get_mock_coord(pdb_id, res_id_350) if res_id_350 else None

        pocket_coords = []
        for grn in pocket_grns:
            res_id = grn_map.get(grn)
            if res_id:
                coord = get_mock_coord(pdb_id, res_id)
                pocket_coords.append(coord)

        # 3f & 3g. Calculate average distance
        if coord_350 is not None and pocket_coords:
            distances = [np.linalg.norm(coord_350 - p_coord) for p_coord in pocket_coords]
            avg_dist = np.mean(distances) if distances else float('inf')
            results[pdb_id] = avg_dist
            print(f"Processed {pdb_id}, Avg Dist: {avg_dist:.2f} Å")
        else:
            print(f"Could not calculate distance for {pdb_id} (missing coords?).")
            results[pdb_id] = float('inf') # Assign infinite distance if calculation failed


    # 4 & 5. Rank results
    ranked_results = sorted(results.items(), key=lambda item: item[1])
    print("\nRanked Structures by Average Distance:")
    for pdb, dist in ranked_results:
        dist_str = f"{dist:.2f} Å" if dist != float('inf') else "N/A"
        print(f"- {pdb}: {dist_str}")

except Exception as e:
    print(f"An error occurred during the workflow: {e}")
```

**Proposed Protos-MCP Implementation:**
- Create MCP Tools for each complex step: `load_structures`, `align_structures`, `extract_binding_pocket`, `assign_grn`, `get_grns_from_residues`, `get_ca_coordinate` (handling alignment), `get_ca_coordinates_for_list` (handling alignment), potentially a dedicated `calculate_average_distance` tool or have the LLM request raw coordinates and perform calculation itself.
- The LLM manages the loop and orchestrates the sequence of Tool calls via the MCP server.

**Expected Outcome (via Protos-MCP):**
- The system returns the final ranked list of structures based on the calculated average distance, directly answering the user's complex comparative question.

## Installation

This section describes how to install Protos and its dependencies. Please note that while the core Protos library can be installed via pip, certain functionalities rely on external bioinformatics tools that must be installed separately. The long-term goal is to provide a Docker image for simplified deployment.

### Prerequisites

- **Python**: Protos requires Python 3.9 or later.
- **Package Manager**: pip (or optionally uv) is used for installing Python packages.

### Recommended Method for Users (Future Goal): Docker

The easiest way to run Protos-MCP with all dependencies correctly configured will be using a provided Docker image (under development). This approach isolates the environment and bundles the core library along with necessary external tools.

**Intended Docker Usage:**
```bash
# Example command (Image name TBD)
docker run -it --rm \
  -v $(pwd)/my_protos_data:/app/data \
  protos-mcp-image:latest
```

This command would run the Protos-MCP container:
- Crucially, the `-v $(pwd)/my_protos_data:/app/data` flag mounts a local directory (`my_protos_data` in your current location) into the container's `/app/data` directory.
- Protos (via its ProtosPaths system) is designed to read from and write to this `/app/data` directory within the container, allowing seamless interaction with your local data files.

(Please check back for updates on the official Docker image release).

### Standard Installation (via PyPI)

You can install the core Protos Python library using pip (Note: Protos is not yet released on PyPI):

```bash
# Command for future release
# pip install protos
```

**Important**: When available, this command will only install the Python library itself. Functionalities relying on external tools (see below) will not work unless those tools are also installed on your system.

### External Dependencies

Protos leverages powerful, established bioinformatics tools for specific tasks. These must be installed separately on your system and be accessible via the system's PATH environment variable for Protos to find and use them.

Key external dependencies currently include:

- **MMseqs2**: Used for fast, sensitive sequence searching and clustering (often required for SeqProcessor alignment tasks).
  - Installation: Please follow the official MMseqs2 installation guide: https://github.com/soedinglab/MMseqs2
- **GTalign**: (If used) A GPU-accelerated tool for fast protein structure alignment. Requires compatible hardware (NVIDIA GPU) and drivers.
  - Installation: Please follow the official GTalign installation guide: https://github.com/BioinfoMachineLearning/GTalign
- **(Other Tools)**: Depending on the specific functionalities enabled (e.g., Cealign, FoldMason), additional tools might be required. Refer to the documentation for specific Processor requirements.

Note: Installation procedures for these tools vary depending on your operating system (Linux, macOS, Windows/WSL).

### Local Development Setup

If you want to contribute to Protos development or run the latest unreleased version:

**Clone the Repository:**
```bash
git clone https://github.com/your-username/protos.git # Replace with actual repo URL
cd protos
```

**Create a Virtual Environment (Recommended):**
```bash
python -m venv venv
source venv/bin/activate  # On Windows use `venv\Scripts\activate`
```

**Install Dependencies**: Install the core requirements and specific development dependencies.
```bash
pip install -r requirements.txt
# Optionally install development dependencies if available
# pip install -r requirements-dev.txt
```

**Install Protos in Editable Mode**: This links the installed package to your local source code.
```bash
pip install -e .
```

**Install External Dependencies**: Remember to install the necessary external tools (MMseqs2, GTalign, etc.) separately as described above, ensuring they are in your system's PATH.

## Configuration

Protos uses a path management system (ProtosPaths) to handle data directories. By default, it creates and uses a data directory relative to where it's run. This location might be configurable via environment variables or programmatically if needed. Refer to the Protos path management documentation for details.
