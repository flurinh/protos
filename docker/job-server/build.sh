#!/bin/bash
# Build the Protos Job Server Docker image

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Copy job_server.py from source
cp "$PROJECT_ROOT/src/protos/models/job_server.py" "$SCRIPT_DIR/"

# Build the image
docker build -t protos/job-server:latest "$SCRIPT_DIR"

# Clean up
rm -f "$SCRIPT_DIR/job_server.py"

echo "Built protos/job-server:latest"
