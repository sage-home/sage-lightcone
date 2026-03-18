#!/bin/bash

# Load utilities
source $(dirname $0)/utils/benchmark_utils.sh

# Exit on error
set -e

# Parse command line arguments
FORCE_REBUILD=0
for arg in "$@"; do
    case $arg in
        --rebuild)
            FORCE_REBUILD=1
            echo "Force rebuild requested"
            ;;
        --help|-h)
            echo "Usage: $0 [OPTIONS]"
            echo "Options:"
            echo "  --rebuild    Force rebuild of all executables (including SAGE)"
            echo "  --help, -h   Show this help message"
            exit 0
            ;;
    esac
done

# What's my code's root directory
SCRIPT="${BASH_SOURCE[0]}"
[ -z "$SCRIPT" ] && SCRIPT="$0"
export MY_SCRIPTS_DIRECTORY=$(cd "$(dirname "$SCRIPT")" && pwd)
export MY_ROOT=$(cd "${MY_SCRIPTS_DIRECTORY}/.." && pwd)

# Store MY_SCRIPTS_DIRECTORY to restore it after sourcing setup scripts
# (setup_mac.sh overwrites it based on its own location)
ORIGINAL_SCRIPTS_DIRECTORY=$MY_SCRIPTS_DIRECTORY

if [[ "$OSTYPE" == "darwin"* ]]; then
    echo "macOS detected - using Homebrew setup"
    source ${MY_ROOT}/setup_mac.sh
elif [[ -f "/etc/modules" ]] || command -v module >/dev/null 2>&1; then
    echo "HPC/Linux environment detected - using module setup"
    source ${MY_ROOT}/setup.sh
else
    echo "Unknown platform - falling back to basic setup"
    source ${MY_ROOT}/setup.sh
fi

# Restore the correct scripts directory for this script
export MY_SCRIPTS_DIRECTORY=$ORIGINAL_SCRIPTS_DIRECTORY

# Change to script directory - required for relative paths like "${MY_SCRIPTS_DIRECTORY}/first_run.sh"
cd "${MY_ROOT}/tests/sage-model-tests"

# Ensure sage-model repository exists (clone if missing)
SAGE_REPO="https://github.com/MBradley1985/SAGE26.git"
if [ ! -d "${MY_ROOT}/sage-model" ]; then
    echo "sage-model directory not found - cloning from ${SAGE_REPO}..."
    pushd ${MY_ROOT} > /dev/null
    git clone ${SAGE_REPO} sage-model
    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to clone SAGE repository"
        exit 1
    fi
    popd > /dev/null
    echo "✓ sage-model cloned successfully"
fi

# Check if required executables exist and are working
echo "Checking required executables..."
NEED_REBUILD=$FORCE_REBUILD

# Check sage executable exists
if [ ! -f "${MY_ROOT}/bin/sage" ]; then
    echo "sage not found - rebuild needed"
    NEED_REBUILD=1
fi

# Check sage2kdtree exists
if [ ! -f "${MY_ROOT}/bin/sage2kdtree" ]; then
    echo "sage2kdtree not found - rebuild needed"
    NEED_REBUILD=1
fi

# Check cli_lightcone can run (catches library loading issues)
if [ -f "${MY_ROOT}/bin/cli_lightcone" ]; then
    ${MY_ROOT}/bin/cli_lightcone --help > /dev/null 2>&1
    if [ $? -ne 0 ]; then
        echo "cli_lightcone has library issues - rebuild needed"
        NEED_REBUILD=1
    fi
fi

if [ $NEED_REBUILD -eq 1 ]; then
    echo "Building required executables..."
    pushd ${MY_ROOT} > /dev/null
    ./build_platform_aware.sh
    BUILD_STATUS=$?
    popd > /dev/null
    if [ $BUILD_STATUS -ne 0 ]; then
        echo "ERROR: Build failed!"
        exit 1
    fi
    echo "Build complete."
fi

export RAWNAME=myhdf5millennium
export OUTPUTDIR=output_sage_hdf5_one_step_benchmark
rm -rf ${OUTPUTDIR}
export BENCHMARK_CSV=${OUTPUTDIR}/benchmark_timings.csv
export DISK_IO_CSV=${OUTPUTDIR}/benchmark_disk_io.csv

# Clean up previous test outputs, but preserve downloaded tree files
rm -rf output
rm -rf ${OUTPUTDIR}
rm -f *${RAWNAME}*
# Only remove input/*.par files, preserve input/millennium/trees/
rm -f input/*.par 2>/dev/null

# Create output directory early for CSV files
mkdir -p ${OUTPUTDIR}

# Initialize CSV files
echo "phase,duration_sec,peak_mem_mb" > $BENCHMARK_CSV
echo "phase,disk_mb" > $DISK_IO_CSV

# Ensure input directory exists and copy base millennium.par from sage-model
mkdir -p input
if [ -f "${MY_ROOT}/sage-model/input/millennium.par" ]; then
    cp "${MY_ROOT}/sage-model/input/millennium.par" input/millennium.par
    echo "✓ Copied base millennium.par from sage-model"
else
    echo "ERROR: sage-model/input/millennium.par not found"
    echo "Please ensure the sage-model submodule is initialized: git submodule update --init --recursive"
    exit 1
fi

# Check if tree files already exist - skip first_run.sh if so
if [ -f "input/millennium/trees/trees_063.7" ] && [ -f "input/millennium/trees/millennium.a_list" ]; then
    echo "✓ Tree files already present - skipping download."
else
    echo "Tree files not found - running first_run.sh to download..."
    "${MY_SCRIPTS_DIRECTORY}/first_run.sh"
    FIRST_RUN_STATUS=$?
    echo ==== Have finished with first_run.sh ====

    # Check if first_run.sh succeeded
    if [ $FIRST_RUN_STATUS -ne 0 ]; then
        echo "ERROR: first_run.sh failed with exit code $FIRST_RUN_STATUS"
        exit 1
    fi

    # Verify tree files were downloaded
    if [ ! -f "input/millennium/trees/trees_063.7" ]; then
        echo "ERROR: Tree files not found. first_run.sh may have failed to download them."
        echo "Expected file: input/millennium/trees/trees_063.7"
        echo "Please check your internet connection and try again."
        exit 1
    fi

    # Verify scale factor list was downloaded
    if [ ! -f "input/millennium/trees/millennium.a_list" ]; then
        echo "ERROR: Scale factor list not found."
        echo "Expected file: input/millennium/trees/millennium.a_list"
        exit 1
    fi

    echo "✓ Tree files and scale factor list verified."
fi

# Update paths in millennium.par to use current directory
CURRENT_DIR=$(pwd)
sed -i'' -e "s|^OutputDir.*|OutputDir   ${CURRENT_DIR}/output/millennium/|" input/millennium.par
sed -i'' -e "s|^SimulationDir.*|SimulationDir   ${CURRENT_DIR}/input/millennium/trees/|" input/millennium.par
sed -i'' -e "s|^FileWithSnapList.*|FileWithSnapList ${CURRENT_DIR}/input/millennium/trees/millennium.a_list|" input/millennium.par
echo "✓ Updated paths in millennium.par"

# Extract settings from millennium.par
echo "Extracting settings from millennium.par..."
python3 "${MY_SCRIPTS_DIRECTORY}/utils/extract_settings.py"

# Generate .par files by concatenating headers with settings
cat mypar_files/millennium_sage_binary_header.txt mypar_files/millennium_settings.txt > input/millennium.par
cat mypar_files/millennium_sage_binary_kdtreeindex_header.txt mypar_files/millennium_settings.txt > input/millennium_minus1.par
cat mypar_files/millennium_sage_hdf5_header.txt mypar_files/millennium_settings.txt > input/millennium_sage_hdf5.par
echo "✓ Parameter files generated."

mkdir -p output/millennium/

# PHASE 1: Run sage
# NOTE: SAGE is run once by benchmark_workflows.sh and shared between workflows
# to ensure identical input data (see SAGE_REPRODUCIBILITY_ISSUE.md)
echo "========== PHASE 1: Running SAGE =========="

# Check if shared SAGE output exists (from benchmark_workflows.sh)

echo "Running SAGE standalone"
run_with_profiling "sage" "$BENCHMARK_CSV" \
    ${MY_ROOT}/bin/sage input/millennium_sage_hdf5.par > /dev/null
mv output/millennium ${OUTPUTDIR}/
mv output/log ${OUTPUTDIR}/ 2>/dev/null || true

measure_disk_usage "${OUTPUTDIR}/millennium" "sage_output" "$DISK_IO_CSV"

# Clean up
rm -f log.00000
rm -rf log

echo ""
echo "=========================================="
echo "Setup for verify kdtree complete."
echo "Results:"
echo "  Timings: ${BENCHMARK_CSV}"
echo "  Disk I/O: ${DISK_IO_CSV}"
echo "=========================================="