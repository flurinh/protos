# Lambda Model Service
# This models requires specific versions due to geometric deep learning dependencies
FROM pytorch/pytorch:1.13.1-cuda11.6-cudnn8-runtime

# Install system dependencies
RUN apt-get update && apt-get install -y \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies with specific versions
RUN pip install --no-cache-dir \
    torch==1.13.1+cu116 \
    torch-scatter==2.1.0+pt113cu116 \
    torch-sparse==0.6.16+pt113cu116 \
    torch-cluster==1.6.0+pt113cu116 \
    torch-spline-conv==1.2.1+pt113cu116 \
    torch-geometric==2.2.0 \
    -f https://data.pyg.org/whl/torch-1.13.0+cu116.html

# Install additional dependencies
RUN pip install --no-cache-dir \
    flask \
    gunicorn \
    numpy \
    pandas \
    scikit-learn \
    biopython

# Create directories
RUN mkdir -p /models /app

# Copy service code
COPY service.py /app/service.py
COPY model_utils.py /app/model_utils.py
WORKDIR /app

# Expose port
EXPOSE 8003

# Health check
HEALTHCHECK --interval=30s --timeout=10s --retries=3 \
    CMD curl -f http://localhost:8003/health || exit 1

# Run service
CMD ["gunicorn", "--bind", "0.0.0.0:8003", "--workers", "1", "--threads", "2", "--timeout", "300", "service:app"]