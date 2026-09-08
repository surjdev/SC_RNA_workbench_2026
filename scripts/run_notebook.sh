#!/usr/bin/env bash
# ==============================================================================
# scripts/run_notebook.sh
# Parameterized Notebook Runner & Report Generator
# Complies with scRNAseq_Workbench_Requirements.md (FR-7, FR-8)
# ==============================================================================

set -euo pipefail

# Text formatting
BOLD="\033[1m"
GREEN="\033[0;32m"
BLUE="\033[0;34m"
YELLOW="\033[0;33m"
RED="\033[0;31m"
NC="\033[0m"

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Default parameters
NOTEBOOK=""
CONFIG="configs/default_analysis.yaml"
OUTPUT_DIR="reports"
GENERATE_HTML=true

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -n|--notebook)
            NOTEBOOK="$2"
            shift 2
            ;;
        -c|--config)
            CONFIG="$2"
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --no-html)
            GENERATE_HTML=false
            shift
            ;;
        -h|--help)
            echo "Usage: $0 -n <notebook_path_or_name> [-c <config_path>] [-o <output_dir>] [--no-html]"
            exit 0
            ;;
        *)
            if [ -z "$NOTEBOOK" ]; then
                NOTEBOOK="$1"
                shift
            else
                echo -e "${RED}Unknown argument: $1${NC}"
                exit 1
            fi
            ;;
    esac
done

if [ -z "$NOTEBOOK" ]; then
    echo -e "${RED}Error: Notebook path or name is required (-n <notebook>).${NC}"
    echo "Usage: $0 -n <notebook> [-c configs/default_analysis.yaml]"
    exit 1
fi

# Resolve notebook path
if [ -f "$NOTEBOOK" ]; then
    NB_PATH="$NOTEBOOK"
elif [ -f "${ROOT_DIR}/notebooks/${NOTEBOOK}" ]; then
    NB_PATH="${ROOT_DIR}/notebooks/${NOTEBOOK}"
elif [ -f "${ROOT_DIR}/notebooks/${NOTEBOOK}.ipynb" ]; then
    NB_PATH="${ROOT_DIR}/notebooks/${NOTEBOOK}.ipynb"
else
    echo -e "${RED}Error: Notebook not found at: ${NOTEBOOK}${NC}"
    exit 1
fi

NB_BASENAME="$(basename "$NB_PATH" .ipynb)"
mkdir -p "${ROOT_DIR}/${OUTPUT_DIR}"

EXECUTED_NB="${ROOT_DIR}/${OUTPUT_DIR}/executed_${NB_BASENAME}.ipynb"
HTML_REPORT="${ROOT_DIR}/${OUTPUT_DIR}/${NB_BASENAME}.html"

echo -e "${BOLD}${BLUE}=== Parameterized Notebook Execution (Papermill) ===${NC}"
echo -e "Notebook:    ${BOLD}${NB_PATH}${NC}"
echo -e "Config:      ${BOLD}${CONFIG}${NC}"
echo -e "Executed NB: ${BOLD}${EXECUTED_NB}${NC}"
echo -e "Report HTML: ${BOLD}${HTML_REPORT}${NC}\n"

# Run command runner helper (uses pixi if available, else python)
RUNNER="python"
if command -v pixi >/dev/null 2>&1 && [ -f "${ROOT_DIR}/pixi.toml" ]; then
    RUNNER="pixi run"
fi

echo -e "${BLUE}1. Executing notebook with Papermill (FR-7)...${NC}"
$RUNNER papermill \
    "$NB_PATH" \
    "$EXECUTED_NB" \
    -k python3 \
    --language python \
    -p config_path "$CONFIG" \
    --log-output

echo -e "${GREEN}✔ Notebook executed successfully.${NC}"

if [ "$GENERATE_HTML" = true ]; then
    echo -e "\n${BLUE}2. Generating HTML report with nbconvert (FR-8)...${NC}"
    $RUNNER jupyter nbconvert \
        --to html \
        --output "$HTML_REPORT" \
        "$EXECUTED_NB"
    echo -e "${GREEN}✔ HTML report generated: ${BOLD}${HTML_REPORT}${NC}"
fi

echo -e "\n${GREEN}${BOLD}✔ Pipeline execution completed successfully!${NC}\n"
