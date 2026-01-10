#!/bin/bash
# Build the protos development Docker image

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_DIR"

echo "Building protos development Docker image..."
docker-compose -f docker-compose.dev.yml build protos

echo ""
echo "Build complete! You can now run:"
echo "  ./scripts/docker-shell.sh    - Interactive shell"
echo "  ./scripts/docker-test.sh     - Run all tests"
echo "  ./scripts/docker-test.sh <test_file.py>  - Run specific test"