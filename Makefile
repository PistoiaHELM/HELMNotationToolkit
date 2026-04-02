# Developer scripts for HELM Notation Toolkit
# Usage: make <command>

.PHONY: test test-s test-coverage lint format build install install-dev clean-build help

# Test commands
test:
	uv run python -m pytest python/tests/ -v

test-s:
	uv run python -m pytest python/tests/ -v -s

test-coverage:
	uv run python -m pytest python/tests/ -v --cov=helm_notation --cov-report=html

# Code quality
lint:
	@echo "Linting (add ruff to dev dependencies to enable)"

format:
	@echo "Formatting (add ruff to dev dependencies to enable)"

# Build commands
build:
	uv build

install:
	uv pip install -e .

install-dev:
	uv pip install -e ".[dev]"

clean-build:
	rm -rf dist/ build/ python/*.egg-info/

# Default target
help:
	@echo "Available commands:"
	@echo "  make test          - Run all tests"
	@echo "  make test-s        - Run tests with stdout output"
	@echo "  make test-coverage - Run tests with coverage report"
	@echo "  make build         - Build package with uv"
	@echo "  make install       - Install package in editable mode"
	@echo "  make install-dev   - Install with dev dependencies"
	@echo "  make clean-build   - Remove build artifacts"