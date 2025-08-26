# Installation Guide

This guide covers the installation of Protos, including prerequisites, installation methods, and post-installation setup.

## Prerequisites

### System Requirements

- **Python**: 3.8 or higher
- **Operating System**: Linux, macOS, or Windows (with WSL recommended)
- **Memory**: 8GB RAM minimum (16GB+ recommended for large datasets)
- **Storage**: 10GB+ free space for software and data

### Required Dependencies

Core dependencies installed automatically:
- `pandas >= 1.3.0`
- `numpy >= 1.21.0`
- `biopython >= 1.79`
- `click >= 8.0.0`
- `pyyaml >= 5.4.0`
- `requests >= 2.26.0`

Optional dependencies for full functionality:
- `torch >= 1.9.0` (for embeddings)
- `transformers >= 4.20.0` (for protein language models)
- `fair-esm >= 2.0.0` (for ESM models)
- `pymol-open-source` (for visualization)

## Installation Methods

### Method 1: Install from PyPI (Recommended)

```bash
# Basic installation
pip install protos

# Full installation with all features
pip install protos[all]

# Specific feature sets
pip install protos[embeddings]  # Embedding generation
pip install protos[structure]   # Structure analysis tools
pip install protos[ml]          # Machine learning features
```

### Method 2: Install from Source

```bash
# Clone repository
git clone https://github.com/protos/protos.git
cd protos

# Install in development mode
pip install -e .

# Or install with all dependencies
pip install -e ".[all]"
```

### Method 3: Conda Installation

```bash
# Create conda environment
conda create -n protos python=3.9
conda activate protos

# Install protos
conda install -c conda-forge protos

# Or use pip within conda
pip install protos[all]
```

## Environment Setup

### 1. Data Directory Configuration

Protos needs a data directory for storing files:

```bash
# Option 1: Use default location
mkdir ~/protos_data

# Option 2: Custom location
export PROTOS_DATA_ROOT=/path/to/your/data
mkdir -p $PROTOS_DATA_ROOT
```

### 2. Initialize Data Structure

```bash
# Initialize directory structure
protos init-data

# Or specify custom location
protos init-data --data-root /path/to/data

# Initialize with sample data
protos init-data --with-examples
```

### 3. Configure Environment Variables

Add to your shell configuration file (`~/.bashrc`, `~/.zshrc`, etc.):

```bash
# Required: Data root directory
export PROTOS_DATA_ROOT=/path/to/protos_data

# Optional: Default processor settings
export PROTOS_DEFAULT_PROCESSOR=structure
export PROTOS_LOG_LEVEL=INFO

# Optional: External tool paths
export MUSCLE_PATH=/usr/local/bin/muscle
export DSSP_PATH=/usr/local/bin/dssp
```

## External Tools Installation

### Structure Analysis Tools

```bash
# DSSP (secondary structure)
conda install -c conda-forge dssp

# Or compile from source
wget https://github.com/cmbi/dssp/releases/download/4.0.0/dssp-4.0.0.tar.gz
tar xzf dssp-4.0.0.tar.gz
cd dssp-4.0.0
cmake .
make
sudo make install
```

### Sequence Alignment Tools

```bash
# MUSCLE
conda install -c bioconda muscle

# Clustal Omega
conda install -c bioconda clustalo

# MAFFT
conda install -c bioconda mafft

# Or using package manager (Ubuntu/Debian)
sudo apt-get install muscle clustalo mafft
```

### Structure Alignment Tools

```bash
# FoldMason
wget https://github.com/steineggerlab/foldmason/releases/download/1.0/foldmason-linux-x86_64.tar.gz
tar xzf foldmason-linux-x86_64.tar.gz
sudo mv foldmason/bin/foldmason /usr/local/bin/

# TM-align
wget https://zhanggroup.org/TM-align/TMalign.cpp
g++ -O3 -o TMalign TMalign.cpp
sudo mv TMalign /usr/local/bin/
```

## Platform-Specific Instructions

### Windows (WSL)

```bash
# Install WSL2
wsl --install

# Inside WSL, follow Linux instructions
sudo apt update
sudo apt install python3-pip python3-dev

# Install protos
pip install protos[all]
```

### macOS

```bash
# Using Homebrew
brew install python@3.9

# Install protos
pip3 install protos[all]

# Install external tools
brew install muscle clustalo mafft
```

### Linux

```bash
# Ubuntu/Debian
sudo apt update
sudo apt install python3-pip python3-dev

# CentOS/RHEL
sudo yum install python3-pip python3-devel

# Install protos
pip3 install protos[all]
```

## Verification

### 1. Check Installation

```bash
# Verify protos is installed
protos --version

# Check available commands
protos --help

# Test import
python -c "import protos; print(protos.__version__)"
```

### 2. Run Test Suite

```bash
# Run basic tests
protos test

# Run full test suite
pytest tests/

# Test specific component
pytest tests/test_processors/test_structure/
```

### 3. Verify External Tools

```bash
# Check tool availability
protos check-tools

# Should show:
# ✓ muscle: Found at /usr/local/bin/muscle
# ✓ clustalo: Found at /usr/local/bin/clustalo
# ✓ dssp: Found at /usr/local/bin/dssp
# ...
```

## GPU Setup (Optional)

For embedding generation with GPU acceleration:

### NVIDIA GPUs

```bash
# Install CUDA toolkit (check PyTorch compatibility)
# Visit: https://developer.nvidia.com/cuda-downloads

# Install PyTorch with CUDA
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118

# Verify GPU access
python -c "import torch; print(torch.cuda.is_available())"
```

### Apple Silicon (M1/M2)

```bash
# Install PyTorch with Metal Performance Shaders
pip install torch torchvision torchaudio

# Verify MPS access
python -c "import torch; print(torch.backends.mps.is_available())"
```

## Docker Installation

### Using Pre-built Image

```bash
# Pull official image
docker pull protos/protos:latest

# Run with data volume
docker run -it -v ~/protos_data:/data protos/protos:latest

# With GPU support
docker run --gpus all -it -v ~/protos_data:/data protos/protos:gpu
```

### Building Custom Image

```dockerfile
# Dockerfile
FROM python:3.9-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    muscle \
    clustalo \
    mafft \
    && rm -rf /var/lib/apt/lists/*

# Install protos
RUN pip install protos[all]

# Set up data directory
ENV PROTOS_DATA_ROOT=/data
RUN mkdir -p /data

WORKDIR /workspace
CMD ["bash"]
```

Build and run:
```bash
docker build -t my-protos .
docker run -it -v $(pwd):/workspace my-protos
```

## Troubleshooting

### Common Installation Issues

1. **Permission Errors**
```bash
# Use user installation
pip install --user protos

# Or use virtual environment
python -m venv protos_env
source protos_env/bin/activate
pip install protos
```

2. **Dependency Conflicts**
```bash
# Create clean environment
conda create -n protos_clean python=3.9
conda activate protos_clean
pip install protos[all]
```

3. **Missing System Libraries**
```bash
# Ubuntu/Debian
sudo apt-get install python3-dev build-essential

# CentOS/RHEL
sudo yum install python3-devel gcc gcc-c++
```

4. **SSL Certificate Errors**
```bash
# Temporary fix (not recommended for production)
pip install --trusted-host pypi.org --trusted-host files.pythonhosted.org protos

# Better: Update certificates
pip install --upgrade certifi
```

### Verifying Components

```python
# Check all components
import protos
from protos.core.base_processor import BaseProcessor
from protos.processing.structure.struct_base_processor import CifBaseProcessor
from protos.processing.sequence.seq_processor import SeqProcessor
from protos.processing.grn.grn_base_processor import GRNBaseProcessor
from protos.processing.embedding.embedding_processor import EmbeddingProcessor

print("All components imported successfully!")

# Check data paths
from protos.io.paths.path_config import ProtosPaths
paths = ProtosPaths()
print(f"Data root: {paths.get_data_root()}")
```

## Next Steps

1. Follow the [Quick Start Tutorial](quickstart.md)
2. Explore [example notebooks](https://github.com/protos/protos/tree/main/examples)
3. Read about [Core Concepts](../core_concepts.md)
4. Check [CLI Commands](commands.md)

## Updating Protos

```bash
# Update to latest version
pip install --upgrade protos

# Update to specific version
pip install protos==2.0.0

# Update from source
cd /path/to/protos
git pull
pip install -e .
```

## Uninstallation

```bash
# Remove protos
pip uninstall protos

# Clean up data (optional)
# WARNING: This removes all data!
rm -rf ~/protos_data

# Remove environment variables
# Edit ~/.bashrc or ~/.zshrc and remove PROTOS_* variables
```