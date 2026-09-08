#!/usr/bin/env bash
# ==============================================================================
# scripts/init_workbench.sh
# Downstream Single-Cell Transcriptomics Analysis Workbench Scaffolding Script
# Complies with Requirement FR-1 (scRNAseq_Workbench_Requirements.md)
# ==============================================================================

set -euo pipefail

# Text formatting
BOLD="\033[1m"
GREEN="\033[0;32m"
BLUE="\033[0;34m"
YELLOW="\033[0;33m"
RED="\033[0;31m"
NC="\033[0m" # No Color

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TEST_MODE=false

# Parse flags
for arg in "$@"; do
    case $arg in
        --test|--dry-run)
            TEST_MODE=true
            shift
            ;;
    esac
done

echo -e "${BOLD}${BLUE}=== scRNA-seq Workbench Scaffolding Setup ===${NC}"
echo -e "Project Root: ${ROOT_DIR}\n"

REQUIRED_DIRS=(
    "data/raw"
    "data/processed"
    "notebooks"
    "src/workbench_utils"
    "configs"
    "reports"
    "tests"
    "scripts"
    "docs"
)

if [ "$TEST_MODE" = true ]; then
    echo -e "${YELLOW}[TEST MODE] Validating existing directory structure...${NC}"
    MISSING=0
    for dir in "${REQUIRED_DIRS[@]}"; do
        if [ -d "${ROOT_DIR}/${dir}" ]; then
            echo -e "  [✔] ${dir}"
        else
            echo -e "  [✘] ${RED}Missing: ${dir}${NC}"
            MISSING=$((MISSING + 1))
        fi
    done

    if [ $MISSING -eq 0 ]; then
        echo -e "\n${GREEN}${BOLD}Verification PASSED: All required directories exist.${NC}"
        exit 0
    else
        echo -e "\n${RED}${BOLD}Verification FAILED: ${MISSING} directory(s) missing.${NC}"
        exit 1
    fi
fi

# Create required directories
echo -e "${BLUE}1. Creating standard directory structure (FR-1)...${NC}"
for dir in "${REQUIRED_DIRS[@]}"; do
    TARGET="${ROOT_DIR}/${dir}"
    if [ ! -d "$TARGET" ]; then
        mkdir -p "$TARGET"
        echo -e "  [+] Created: ${dir}"
    else
        echo -e "  [=] Exists:  ${dir}"
    fi
done

# Ensure .gitkeep in empty tracking directories
echo -e "\n${BLUE}2. Ensuring .gitkeep files...${NC}"
for keep_dir in "data/raw" "data/processed" "configs" "reports"; do
    KEEP_FILE="${ROOT_DIR}/${keep_dir}/.gitkeep"
    if [ ! -f "$KEEP_FILE" ]; then
        touch "$KEEP_FILE"
        echo -e "  [+] Created: ${keep_dir}/.gitkeep"
    fi
done

# Check system tools
echo -e "\n${BLUE}3. Checking required system tools...${NC}"
for tool in git python3; do
    if command -v "$tool" >/dev/null 2>&1; then
        echo -e "  [✔] ${tool}: $(command -v "$tool")"
    else
        echo -e "  [!] ${YELLOW}${tool} is not in PATH${NC}"
    fi
done

for pkg_tool in pixi mamba conda; do
    if command -v "$pkg_tool" >/dev/null 2>&1; then
        echo -e "  [✔] ${pkg_tool}: $(command -v "$pkg_tool")"
    else
        echo -e "  [ ] ${pkg_tool}: not found (optional alternative)"
    fi
done

echo -e "\n${GREEN}${BOLD}✔ Workbench scaffolding initialized successfully!${NC}"
echo -e "Next steps: Run ${BOLD}make setup${NC} or see ${BOLD}HANDOFF.md${NC} for Part 2.\n"
