# Protos Development & Testing Container
FROM python:3.10-slim

# Set environment variables
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    curl \
    build-essential \
    mmseqs2 \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy requirements first for better caching
COPY requirements.txt setup.py pyproject.toml MANIFEST.in README.md ./
COPY src/ ./src/

# Install the package with all dependencies
RUN pip install -e ".[all]"

# Install PyTorch Geometric dependencies for lambda model
# Install compatible versions for the installed PyTorch
RUN TORCH_VERSION=$(python -c "import torch; print(torch.__version__.split('+')[0])") && \
    pip install torch_geometric && \
    pip install torch_scatter torch_sparse torch_cluster torch_spline_conv \
    -f https://data.pyg.org/whl/torch-${TORCH_VERSION}+cpu.html || \
    pip install torch_scatter torch_sparse torch_cluster torch_spline_conv

# Copy the rest of the project (tests, examples, etc.)
COPY . .

# Install the lambda (lmda) model package with all dependencies
RUN pip install -e src/protos/models/lambda/ && \
    pip install lmdb kaleido sparsemax tensorboard

# Create data directory for protos
RUN mkdir -p /data/protos

# Set environment variable for protos data directory
ENV PROTOS_DATA_ROOT=/data/protos

# Default command: run pytest
CMD ["pytest", "-v"]