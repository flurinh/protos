#!/bin/bash
# Run tests in the protos Docker container
#
# Usage:
#   ./scripts/docker-test.sh              # Run all tests
#   ./scripts/docker-test.sh test_foo.py  # Run specific test file
#   ./scripts/docker-test.sh -k "test_name"  # Run tests matching pattern

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_DIR"

if [ $# -eq 0 ]; then
    echo "Running all tests..."
    docker-compose -f docker-compose.dev.yml run --rm test
else
    echo "Running: pytest $@"
    docker-compose -f docker-compose.dev.yml run --rm pytest "$@"
fi
