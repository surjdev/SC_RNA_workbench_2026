# ==============================================================================
# Makefile
# Downstream Single-Cell Transcriptomics Analysis Workbench
# Complies with scRNAseq_Workbench_Requirements.md (FR-7, FR-8, Deliverable 7)
# ==============================================================================

.PHONY: help setup test lint format run report validate-config clean

SHELL := /bin/bash
NB ?= 00_end_to_end
CONFIG ?= configs/default_analysis.yaml

help:
	@echo "Available targets:"
	@echo "  make setup                - Initialize directory scaffold and environment"
	@echo "  make test                 - Run pytest suite (unit tests)"
	@echo "  make lint                 - Run ruff linter"
	@echo "  make format               - Format code using ruff"
	@echo "  make run NB=<name>        - Run parameterized notebook via papermill (FR-7)"
	@echo "  make report NB=<name>     - Execute notebook and export HTML report to reports/ (FR-8)"
	@echo "  make validate-config      - Inspect and validate YAML parameters (CONFIG=<path>)"
	@echo "  make clean                - Clean transient caches and checkpoints"

setup:
	@echo "Initializing directory scaffold..."
	@bash scripts/init_workbench.sh
	@echo "Setting up environment via pixi..."
	@pixi install
	@pixi run install-editable
	@pixi run pre-commit install || true
	@echo "Setup complete."

test:
	@pixi run pytest tests/ -v

lint:
	@pixi run ruff check src/ tests/ notebooks/

format:
	@pixi run ruff format src/ tests/ notebooks/

validate-config:
	@pixi run python -m workbench_utils.config $(CONFIG)

run:
	@bash scripts/run_notebook.sh -n $(NB) -c $(CONFIG) --no-html

report:
	@bash scripts/run_notebook.sh -n $(NB) -c $(CONFIG)

clean:
	@rm -rf .pytest_cache .coverage htmlcov __pycache__ src/**/__pycache__ tests/__pycache__
	@find . -name "*.pyc" -delete
	@find . -name ".ipynb_checkpoints" -type d -exec rm -rf {} +
	@echo "Clean complete."
