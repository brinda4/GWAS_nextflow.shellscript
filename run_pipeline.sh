#!/usr/bin/env bash
# =============================================================================
# GWAS Analysis Pipeline Runner
# =============================================================================
# Supports:
#   1. Nextflow workflow execution (nextflow run main.nf)
#   2. Direct Python / UV CLI execution fallback (python run_pipeline.py --all)
# =============================================================================

set -e

# Automatically discover uv or python in common Windows/WSL/Linux paths
if [ -d "/mnt/c/Users/tsbri/AppData/Local/Programs/Python/Python313/Scripts" ]; then
    export PATH="$PATH:/mnt/c/Users/tsbri/AppData/Local/Programs/Python/Python313/Scripts"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Text styles & colors
BOLD="\033[1m"
GREEN="\033[0;32m"
BLUE="\033[0;34m"
YELLOW="\033[1;33m"
RED="\033[0;31m"
NC="\033[0m" # No Color

echo -e "${BOLD}${BLUE}======================================================================${NC}"
echo -e "${BOLD}${BLUE}             GWAS ANALYSIS PIPELINE EXECUTION SCRIPT                  ${NC}"
echo -e "${BOLD}${BLUE}======================================================================${NC}"

# Default parameters
PROFILE="standard"
RUN_MODE="auto"
OUTDIR="results"
EXTRA_ARGS=""

show_help() {
    cat << 'EOF'
Usage: ./run_pipeline.sh [options]

Options:
  --nextflow          Force execution using Nextflow (main.nf)
  --python            Force direct execution using Python (run_pipeline.py)
  --profile <prof>    Nextflow profile: standard, conda, docker, slurm (default: standard)
  --outdir <dir>      Output directory (default: results)
  --install-nextflow  Download and install standalone Nextflow binary locally
  -h, --help          Display this help message

Examples:
  ./run_pipeline.sh                  # Automatically detects Nextflow or runs Python/UV
  ./run_pipeline.sh --nextflow       # Run using Nextflow
  ./run_pipeline.sh --python         # Run using Python directly
  ./run_pipeline.sh --profile conda  # Run Nextflow with conda environment
  ./run_pipeline.sh --profile docker # Run Nextflow with docker container
EOF
    exit 0
}

# Parse command-line arguments
while [ $# -gt 0 ]; do
    case "$1" in
        --nextflow)
            RUN_MODE="nextflow"
            shift
            ;;
        --python)
            RUN_MODE="python"
            shift
            ;;
        --profile)
            PROFILE="$2"
            shift 2
            ;;
        --outdir)
            OUTDIR="$2"
            shift 2
            ;;
        --install-nextflow)
            echo -e "${YELLOW}Installing Nextflow locally...${NC}"
            curl -s https://get.nextflow.io | bash
            chmod +x nextflow
            echo -e "${GREEN}Nextflow downloaded as ./nextflow${NC}"
            if ! command -v java &>/dev/null; then
                echo -e "${YELLOW}Note: Nextflow requires Java (version 11 or later). Please ensure Java is installed.${NC}"
            fi
            exit 0
            ;;
        -h|--help)
            show_help
            ;;
        *)
            EXTRA_ARGS="$EXTRA_ARGS $1"
            shift
            ;;
    esac
done

# Mode determination
HAS_JAVA=false
if command -v java &>/dev/null; then
    HAS_JAVA=true
fi

HAS_NEXTFLOW=false
if [ -f "./nextflow" ] || command -v nextflow &>/dev/null; then
    HAS_NEXTFLOW=true
fi

if [ "$RUN_MODE" = "auto" ]; then
    if [ "$HAS_NEXTFLOW" = true ] && [ "$HAS_JAVA" = true ]; then
        RUN_MODE="nextflow"
    else
        RUN_MODE="python"
    fi
fi

if [ "$RUN_MODE" = "nextflow" ]; then
    if [ "$HAS_JAVA" = false ]; then
        echo -e "${RED}[!] Error: Nextflow requires Java (version 11 or later), but 'java' was not found on PATH.${NC}"
        echo -e "${YELLOW}Tip: You can install OpenJDK (e.g., 'sudo apt install default-jre') or run directly with Python via:${NC}"
        echo -e "     ${GREEN}./run_pipeline.sh --python${NC}"
        exit 1
    fi

    NF_BIN="nextflow"
    if [ -f "./nextflow" ] && ! command -v nextflow &>/dev/null; then
        NF_BIN="./nextflow"
    fi
    
    echo -e "${BOLD}${GREEN}[►] Executing Nextflow Pipeline (main.nf)...${NC}"
    echo -e "Profile         : ${YELLOW}${PROFILE}${NC}"
    echo -e "Output Directory: ${YELLOW}${OUTDIR}${NC}"
    
    $NF_BIN run main.nf \
        -profile "$PROFILE" \
        --outdir "$OUTDIR" \
        $EXTRA_ARGS

    echo -e "${BOLD}${GREEN}======================================================================${NC}"
    echo -e "${BOLD}${GREEN}     Nextflow Pipeline Execution Completed Successfully!              ${NC}"
    echo -e "${BOLD}${GREEN}======================================================================${NC}"
    echo -e "Results are available in: ${YELLOW}${OUTDIR}${NC}"
    echo -e "  - Tables : ${OUTDIR}/tables"
    echo -e "  - Figures: ${OUTDIR}/figures"
    echo -e "  - Reports: ${OUTDIR}/reports"

else
    echo -e "${BOLD}${YELLOW}[►] Running Standalone Python Pipeline...${NC}"
    
    # Check if uv is available
    if command -v uv &>/dev/null; then
        echo -e "${GREEN}[✓] Using 'uv' fast package runner:${NC}"
        uv run --with pandas,numpy,matplotlib,seaborn,scipy python run_pipeline.py --all $EXTRA_ARGS
    else
        # Find standard python
        PYTHON_BIN=""
        if command -v python3 &>/dev/null; then
            PYTHON_BIN="python3"
        elif command -v python &>/dev/null; then
            PYTHON_BIN="python"
        elif command -v py &>/dev/null; then
            PYTHON_BIN="py"
        else
            echo -e "${RED}Error: Neither 'uv' nor 'python' was found on PATH.${NC}"
            exit 1
        fi
        
        PY_VER=$($PYTHON_BIN --version 2>&1 || true)
        echo -e "${GREEN}[✓] Using Python: ${PY_VER}${NC}"
        
        # Check dependencies
        $PYTHON_BIN -c "import pandas, numpy, matplotlib, seaborn, scipy" 2>/dev/null || {
            echo -e "${YELLOW}Installing required Python dependencies...${NC}"
            $PYTHON_BIN -m pip install -r requirements.txt
        }
        
        $PYTHON_BIN run_pipeline.py --all $EXTRA_ARGS
    fi
    
    echo -e "${BOLD}${GREEN}======================================================================${NC}"
    echo -e "${BOLD}${GREEN}     Pipeline Execution Completed Successfully!                        ${NC}"
    echo -e "${BOLD}${GREEN}======================================================================${NC}"
    echo -e "Processed tables: ${YELLOW}data/processed/${NC}"
    echo -e "Figures         : ${YELLOW}figures/${NC}"
fi
