<h1 align="center">ProtOS</h1>

<p align="center">
  <b>A zero-configuration Python framework for structural biology.</b><br>
  Structures, sequences, GRNs, properties, and ML embeddings — every piece of data is a named entity in one registry.
</p>

<p align="center"><img src="docs/architecture.jpg" alt="ProtOS architecture: ProtosPath, Registry/Entity/Dataset, processors and loaders, and a Model Manager with Boltz, Ankh and Lambda backends" width="760"></p>

<p align="center"><a href="https://flurinh.github.io/aboutme">◆ Portfolio</a></p>

<p align="center"><i>You may also be interested in</i></p>

<table align="center"><tr>
<td align="left">←&nbsp; <b>Previous work</b><br><a href="https://github.com/flurinh/mt">Master thesis — GPCRs as graphs</a></td>
<td width="56"></td>
<td align="right"><b>Continuation of this project</b> &nbsp;→<br><a href="https://github.com/flurinh/MOGRN">MOGRN — one numbering for every opsin</a></td>
</tr></table>

---

## What it is

ProtOS handles structural-biology data the way a good toolkit should: you work with
**human-readable names** (`1ubq`, `EGFR_HUMAN`) and the framework manages every file path,
format conversion, and dataset behind the scenes. A modular **processor** handles each data
type, a central **registry** tracks everything, **datasets** group entities for batch work,
and a **Model Manager** dispatches jobs to ML backends (Boltz, Ankh, Lambda).

It's the analysis engine the rest of my PhD work is built on.

## Quick start

```bash
git clone https://github.com/flurinh/protos.git
cd protos
python -m pip install -e .
```

```python
from protos.io.ingest.structure_loader import StructureLoader
from protos.processing.structure import StructureProcessor

sp = StructureProcessor()
loader = StructureLoader(processor=sp)

# Fetch an RCSB mmCIF file, register it, and store its canonical representation.
name = loader.download_and_register("1ubq")
if name is None:
    raise RuntimeError("RCSB download failed")

structure = sp.load_entity(name)

# Group entities once, operate on the whole set
sp.create_dataset("example_structures", [name])
structures = sp.load_dataset("example_structures", return_format="dict")
```

Structures, sequences, GRNs, embeddings — all one tracked registry.

## Core ideas

- **Zero configuration** — set a data path once (or use the default); never touch a file path again.
- **Entity registry** — every structure / sequence / property / embedding is an entity with a name, linkable across formats.
- **Datasets** — named collections for batch operations across any processor.
- **Processors** — Structure, Sequence, GRN, Embedding, Graph, and Property, each with its own loaders.
- **Model Manager** — submit and track ML jobs (structure prediction, embeddings, property prediction).

## What it powers

| Project | What it adds |
|---------|--------------|
| **[MOGRN](https://github.com/flurinh/MOGRN)** | a generic residue numbering for type-I opsins, built on the GRN processor |
| **[Lambda](https://github.com/flurinh/lambda)** | a model that predicts opsin absorption (λmax), registered as a Model-Manager backend |
| **[ProtOS-MCP](https://github.com/flurinh/Protos_MCP)** | exposes the whole toolkit to LLM agents over the Model Context Protocol |

More detail lives in [`docs/`](docs/) — unified data access, the entity registry, zero-configuration design, and the Model Manager.
