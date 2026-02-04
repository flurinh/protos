# Protos Job Server
# FastAPI server for managing model execution jobs

FROM python:3.11-slim

LABEL maintainer="protos"
LABEL description="Protos Job Server for model execution management"

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    docker.io \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy server code
COPY job_server.py .

# Create data directory
RUN mkdir -p /data/models

# Environment variables
ENV PROTOS_DATA_DIR=/data/models
ENV PYTHONUNBUFFERED=1

# Expose port
EXPOSE 8000

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:8000/health || exit 1

# Run server
CMD ["uvicorn", "job_server:app", "--host", "0.0.0.0", "--port", "8000"]
