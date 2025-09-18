# ESM-2 Model Service
FROM pytorch/pytorch:2.0.1-cuda11.7-cudnn8-runtime

# Install dependencies
RUN pip install --no-cache-dir \
    fair-esm \
    flask \
    gunicorn \
    numpy \
    pandas

# Create model directory
RUN mkdir -p /models

# Copy service code
COPY service.py /app/service.py
WORKDIR /app

# Expose port
EXPOSE 8001

# Health check
HEALTHCHECK --interval=30s --timeout=10s --retries=3 \
    CMD curl -f http://localhost:8001/health || exit 1

# Run service
CMD ["gunicorn", "--bind", "0.0.0.0:8001", "--workers", "1", "--threads", "2", "--timeout", "300", "service:app"]