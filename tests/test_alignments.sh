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
COMPUTE_SCRIPT="$REPO_ROOT/compute_alignments/compute_pairwise_alignments.sh"
NTHREADS="${NTHREADS:-2}"
MAX_MEMORY_GB="${MAX_MEMORY_GB:-8}"
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
    "snakemake"
    "sem"
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

# Note: Using dummy fasta files for fast testing (created below)
# Skip UCSC verification since we're not downloading real genomes
print_info "Using dummy fasta files for fast testing (skipping UCSC download)"

# Set resource limits
print_info "Setting resource limits..."
print_info "  Max memory: ${MAX_MEMORY_GB}GB"
print_info "  Max threads: $NTHREADS"

# Set memory limit (ulimit -v is in KB, so convert GB to KB)
MEMORY_KB=$((MAX_MEMORY_GB * 1024 * 1024))
ulimit -v $MEMORY_KB 2>/dev/null || print_warn "Could not set memory limit (ulimit -v)"

# Limit CPU time (optional, in seconds - 2 hours default)
ulimit -t 7200 2>/dev/null || print_warn "Could not set CPU time limit"

# Use nice to lower process priority (prevents system freeze)
# nice values: -20 (highest) to 19 (lowest), default is 10
NICE_LEVEL="${NICE_LEVEL:-10}"
print_info "  Process priority: nice $NICE_LEVEL"

# Clean up test output directory if it exists
if [ -d "$TEST_OUTPUT_DIR" ]; then
    print_warn "Removing existing test output directory: $TEST_OUTPUT_DIR"
    rm -rf "$TEST_OUTPUT_DIR"
fi

# Create dummy fasta files for fast testing
create_dummy_fastas() {
    print_info "Creating dummy fasta files for fast testing..."
    FASTA_DIR="$TEST_OUTPUT_DIR/fasta"
    ASSEMBLY_DIR="$TEST_OUTPUT_DIR/assembly"
    mkdir -p "$FASTA_DIR" "$ASSEMBLY_DIR"
    
    # Create a base sequence that will be shared/modified across species
    # This ensures sequences are alignable
    BASE_SEQ=$(python3 -c "import random; random.seed(42); print(''.join(random.choices('ATCG', k=1000)))")
    
    for i in "${!SPECIES_ARRAY[@]}"; do
        species="${SPECIES_ARRAY[$i]}"
        DUMMY_FASTA="$FASTA_DIR/${species}.fa"
        
        # Create sequences with ~80% similarity to base sequence (mutations)
        # This ensures they will align but with some differences
        python3 <<PYTHON > "$DUMMY_FASTA"
import random
random.seed(42 + $i)  # Different seed per species

base_seq = "$BASE_SEQ"

def mutate_sequence(seq, mutation_rate=0.2):
    """Mutate sequence with given mutation rate"""
    result = list(seq)
    for i in range(len(result)):
        if random.random() < mutation_rate:
            # Mutate to a different base
            bases = ['A', 'T', 'C', 'G']
            bases.remove(result[i])
            result[i] = random.choice(bases)
    return ''.join(result)

# Create two chromosomes with shared similarity
chr1 = mutate_sequence(base_seq, 0.15)  # 85% similarity
chr2 = mutate_sequence(base_seq, 0.20)  # 80% similarity

print(">chr1")
print(chr1)
print(">chr2")
print(chr2)
PYTHON
        
        # Create corresponding .sizes file (needed by pipeline)
        SIZES_FILE="$ASSEMBLY_DIR/${species}.sizes"
        cat > "$SIZES_FILE" <<EOF
chr1	1000
chr2	1000
EOF
        print_info "  Created: $DUMMY_FASTA ($(wc -c < "$DUMMY_FASTA" | xargs) bytes)"
    done
}

create_dummy_fastas

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
print_info "Running: nice -n $NICE_LEVEL $COMPUTE_SCRIPT -s $SPECIES_FILE -r $FIRST_SPECIES -q $SECOND_SPECIES -c -d $TEST_OUTPUT_DIR -@ $NTHREADS"

if ! nice -n $NICE_LEVEL bash "$COMPUTE_SCRIPT" \
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

# Summary
print_info "=========================================="
print_info "Test completed successfully!"
print_info "=========================================="
print_info "Output directory: $TEST_OUTPUT_DIR"
print_info "FASTA files: ${#SPECIES_ARRAY[@]}"
print_info "MAF files: $MAF_COUNT"
print_info "Chain files: $CHAIN_COUNT"
print_info "Pwaln collection: 1"
print_info "=========================================="

