#!/bin/bash
# Start an interactive shell in the protos container

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_DIR"

echo "Starting protos development container..."
docker-compose -f docker-compose.dev.yml run --rm protos