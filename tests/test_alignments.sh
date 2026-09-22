#!/bin/bash

# End-to-end test script for the IPP alignment pipeline
# Tests the complete pipeline from genome download to alignment output

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Configuration
SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
SPECIES_FILE="$REPO_ROOT/data/species_test.txt"
TEST_OUTPUT_DIR="$REPO_ROOT/test_alignments_output"
ORCHESTRATOR_SCRIPT="$REPO_ROOT/compute_alignments/orchestrator.py"
NTHREADS="${NTHREADS:-2}"
MAX_MEMORY_GB="${MAX_MEMORY_GB:-8}"

# Required tools
REQUIRED_TOOLS=(
    "lastal"
    "lastdb"
    "axtChain"
    "chainMergeSort"
    "chainPreNet"
    "maf-convert"
    "faToTwoBit"
    "twoBitInfo"
    "wget"
    "gzip"
)

# Function to print colored messages
print_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to check if a command exists
check_command() {
    if command -v "$1" &> /dev/null; then
        return 0
    else
        return 1
    fi
}

# Check all required tools
print_info "Checking required tools..."
MISSING_TOOLS=()
for tool in "${REQUIRED_TOOLS[@]}"; do
    if ! check_command "$tool"; then
        MISSING_TOOLS+=("$tool")
        print_error "Missing: $tool"
    else
        print_info "Found: $tool"
    fi
done

if [ ${#MISSING_TOOLS[@]} -ne 0 ]; then
    print_error "Missing required tools: ${MISSING_TOOLS[*]}"
    print_error "Please install the missing tools before running the test."
    exit 1
fi

# Check if species file exists
if [ ! -f "$SPECIES_FILE" ]; then
    print_error "Species file not found: $SPECIES_FILE"
    exit 1
fi

# Read species from file (skip comments and empty lines)
SPECIES_LIST=$(awk '/^[^#]/ && NF>0 {print $1}' "$SPECIES_FILE" | tr '\n' ',' | sed 's/,$//')
if [ -z "$SPECIES_LIST" ]; then
    print_error "No species found in $SPECIES_FILE"
    exit 1
fi

print_info "Testing with species: $SPECIES_LIST"

# Parse species array
SPECIES_ARRAY=($(echo "$SPECIES_LIST" | tr ',' ' '))

# Copy committed test FASTA files into the test output directory
print_info "Copying test FASTA files (skipping UCSC download)..."

TEST_FASTA_DIR="$REPO_ROOT/data/test_fasta"
if [ ! -d "$TEST_FASTA_DIR" ]; then
    print_error "Test FASTA directory not found: $TEST_FASTA_DIR"
    exit 1
fi

# Clean up test output directory if it exists
if [ -d "$TEST_OUTPUT_DIR" ]; then
    print_warn "Removing existing test output directory: $TEST_OUTPUT_DIR"
    rm -rf "$TEST_OUTPUT_DIR"
fi

FASTA_DIR="$TEST_OUTPUT_DIR/fasta"
mkdir -p "$FASTA_DIR"
for species in "${SPECIES_ARRAY[@]}"; do
    SRC="$TEST_FASTA_DIR/${species}.fa"
    if [ ! -f "$SRC" ]; then
        print_error "Missing test FASTA for species '$species': $SRC"
        exit 1
    fi
    cp "$SRC" "$FASTA_DIR/${species}.fa"
    print_info "  Copied: ${species}.fa"
done

# Set resource limits
NICE_LEVEL="${NICE_LEVEL:-10}"
print_info "Running with: threads=$NTHREADS, nice=$NICE_LEVEL"

# Run the alignment pipeline
print_info "Running alignment pipeline..."
print_info "Output directory: $TEST_OUTPUT_DIR"
print_info "Threads: $NTHREADS"
print_info "Memory limit: ${MAX_MEMORY_GB}GB"

# Get first two species for pwaln collection test
FIRST_SPECIES=$(echo "$SPECIES_LIST" | cut -d',' -f1)
SECOND_SPECIES=$(echo "$SPECIES_LIST" | cut -d',' -f2)

if [ -z "$FIRST_SPECIES" ] || [ -z "$SECOND_SPECIES" ]; then
    print_error "Need at least 2 species for testing"
    exit 1
fi

# Run the pipeline with pwaln collection creation
# Use nice to lower priority and prevent system freeze
print_info "Running: nice -n $NICE_LEVEL uv run python $ORCHESTRATOR_SCRIPT -s $SPECIES_FILE -r $FIRST_SPECIES -q $SECOND_SPECIES -c -d $TEST_OUTPUT_DIR -@ $NTHREADS"

if ! nice -n $NICE_LEVEL uv run python "$ORCHESTRATOR_SCRIPT" \
    -s "$SPECIES_FILE" \
    -r "$FIRST_SPECIES" \
    -q "$SECOND_SPECIES" \
    -c \
    -d "$TEST_OUTPUT_DIR" \
    -@ "$NTHREADS"; then
    print_error "Pipeline failed!"
    exit 1
fi

# Verify outputs
print_info "Verifying outputs..."

# Check for FASTA files
print_info "Checking FASTA files..."
for species in "${SPECIES_ARRAY[@]}"; do
    FASTA_FILE="$TEST_OUTPUT_DIR/fasta/${species}.fa"
    if [ -f "$FASTA_FILE" ]; then
        SIZE=$(du -h "$FASTA_FILE" | cut -f1)
        print_info "✓ Found: $FASTA_FILE ($SIZE)"
    else
        print_error "✗ Missing: $FASTA_FILE"
        exit 1
    fi
done

# Check for LAST database files
print_info "Checking LAST database files..."
for species in "${SPECIES_ARRAY[@]}"; do
    LASTDB_FILE="$TEST_OUTPUT_DIR/alignment/lastdb/${species}.prj"
    if [ -f "$LASTDB_FILE" ]; then
        print_info "✓ Found: $LASTDB_FILE"
    else
        print_error "✗ Missing: $LASTDB_FILE"
        exit 1
    fi
done

# Check for MAF files (pairwise alignments)
print_info "Checking MAF alignment files..."
MAF_COUNT=0
for i in "${!SPECIES_ARRAY[@]}"; do
    for j in "${!SPECIES_ARRAY[@]}"; do
        if [ $i -ne $j ]; then
            S1="${SPECIES_ARRAY[$i]}"
            S2="${SPECIES_ARRAY[$j]}"
            MAF_FILE="$TEST_OUTPUT_DIR/alignment/maf/${S1}.${S2}.maf"
            if [ -f "$MAF_FILE" ]; then
                SIZE=$(du -h "$MAF_FILE" | cut -f1)
                print_info "✓ Found: $MAF_FILE ($SIZE)"
                MAF_COUNT=$((MAF_COUNT + 1))
            else
                print_error "✗ Missing: $MAF_FILE"
                exit 1
            fi
        fi
    done
done

# Check for chain files
print_info "Checking chain files..."
CHAIN_COUNT=0
for i in "${!SPECIES_ARRAY[@]}"; do
    for j in "${!SPECIES_ARRAY[@]}"; do
        if [ $i -ne $j ]; then
            S1="${SPECIES_ARRAY[$i]}"
            S2="${SPECIES_ARRAY[$j]}"
            CHAIN_FILE="$TEST_OUTPUT_DIR/alignment/chain/${S1}.${S2}.all.pre.chain"
            if [ -f "$CHAIN_FILE" ]; then
                SIZE=$(du -h "$CHAIN_FILE" | cut -f1)
                print_info "✓ Found: $CHAIN_FILE ($SIZE)"
                CHAIN_COUNT=$((CHAIN_COUNT + 1))
            else
                print_error "✗ Missing: $CHAIN_FILE"
                exit 1
            fi
        fi
    done
done

# Check for pwaln collection file
PWALN_FILE="$TEST_OUTPUT_DIR/alignment/pwaln/${FIRST_SPECIES}.${SECOND_SPECIES}.pwaln.bin"
if [ -f "$PWALN_FILE" ]; then
    SIZE=$(du -h "$PWALN_FILE" | cut -f1)
    print_info "✓ Found pwaln collection: $PWALN_FILE ($SIZE)"
else
    print_error "✗ Missing pwaln collection: $PWALN_FILE"
    exit 1
fi

# Run projection
print_info "Running projection from $FIRST_SPECIES to $SECOND_SPECIES..."
PROJ_OUTPUT_DIR="$TEST_OUTPUT_DIR/projection"
mkdir -p "$PROJ_OUTPUT_DIR"
REGIONS_FILE="$REPO_ROOT/data/test_regions.bed"

# Use distance thresholds scaled for the test genome (~20 KB).
# The default -ddc 300 is calibrated for ~3 GB mammalian genomes; for a 20 KB
# genome it collapses to a score threshold of ~0, classifying everything as DC.
# 0.0001 bp gives score_DC = score_IC = 0.999, which correctly distinguishes
# on-block (score 1.0) from off-block (score 0.0) positions.
if ! uv run python "$REPO_ROOT/src/ipp/project.py" \
    -o "$PROJ_OUTPUT_DIR" \
    -n 1 \
    -ddc 0.0001 -dic 0.0001 \
    "$REGIONS_FILE" \
    "$FIRST_SPECIES" \
    "$SECOND_SPECIES" \
    "$PWALN_FILE"; then
    print_error "Projection failed!"
    exit 1
fi

# Verify projection outputs
REGIONS_BASENAME=$(basename "$REGIONS_FILE" .bed)
PROJ_FILE="$PROJ_OUTPUT_DIR/${REGIONS_BASENAME}.${FIRST_SPECIES}-${SECOND_SPECIES}.proj"
UNMAPPED_FILE="$PROJ_OUTPUT_DIR/${REGIONS_BASENAME}.${FIRST_SPECIES}-${SECOND_SPECIES}.unmapped"

if [ -f "$PROJ_FILE" ]; then
    N_PROJECTED=$(( $(wc -l < "$PROJ_FILE") - 1 ))  # subtract header
    print_info "✓ Found projection results: $PROJ_FILE ($N_PROJECTED regions projected)"
else
    print_error "✗ Missing projection output: $PROJ_FILE"
    exit 1
fi

if [ -f "$UNMAPPED_FILE" ]; then
    N_UNMAPPED=$(grep -c . "$UNMAPPED_FILE" || true)
    print_info "✓ Found unmapped regions file: $UNMAPPED_FILE ($N_UNMAPPED unmapped)"
else
    print_error "✗ Missing unmapped regions file: $UNMAPPED_FILE"
    exit 1
fi

if [ "$N_PROJECTED" -eq 0 ]; then
    print_error "✗ No regions were projected — check alignment quality"
    exit 1
fi

# Verify that all three conservation classes are represented
N_DC=$(awk -F'\t' 'NR>1 && $7=="DC"' "$PROJ_FILE" | wc -l | tr -d ' ')
N_IC=$(awk -F'\t' 'NR>1 && $7=="IC"' "$PROJ_FILE" | wc -l | tr -d ' ')
N_NC=$(awk -F'\t' 'NR>1 && $7=="NC"' "$PROJ_FILE" | wc -l | tr -d ' ')
print_info "Conservation classes: DC=$N_DC IC=$N_IC NC=$N_NC"

if [ "$N_DC" -eq 0 ]; then
    print_error "✗ No DC regions found — check alignment coverage"
    exit 1
fi
if [ "$N_IC" -eq 0 ]; then
    print_error "✗ No IC regions found — check bridging species alignment"
    exit 1
fi
if [ "$N_NC" -eq 0 ]; then
    print_error "✗ No NC regions found — check gap positions in test_regions.bed"
    exit 1
fi
print_info "✓ All three conservation classes (DC/IC/NC) present"

# Summary
print_info "=========================================="
print_info "Test completed successfully!"
print_info "=========================================="
print_info "Output directory: $TEST_OUTPUT_DIR"
print_info "FASTA files: ${#SPECIES_ARRAY[@]}"
print_info "MAF files: $MAF_COUNT"
print_info "Chain files: $CHAIN_COUNT"
print_info "Pwaln collection: 1"
print_info "Projected regions: $N_PROJECTED / $(grep -c . "$REGIONS_FILE")"
print_info "=========================================="

