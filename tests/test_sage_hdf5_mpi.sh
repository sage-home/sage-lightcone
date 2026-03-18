#!/bin/bash

# Load utilities
source $(dirname $0)/utils/benchmark_utils.sh

# Exit on error
set -e

# Parse command line arguments
FORCE_REBUILD=0
MPI_TASKS=1

while [[ $# -gt 0 ]]; do
    case $1 in
        --rebuild)
            FORCE_REBUILD=1
            echo "Force rebuild requested"
            shift
            ;;
        --np)
            MPI_TASKS="$2"
            shift 2
            ;;
        --help|-h)
            echo "Usage: $0 [OPTIONS]"
            echo "Options:"
            echo "  --rebuild    Force rebuild of all executables (including SAGE)"
            echo "  --np N       Number of MPI tasks (default: 1)"
            echo "  --help, -h   Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
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

# Change to script directory - required for relative paths like ./first_run.sh
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
export OUTPUTDIR=output_sage_hdf5_mpi${MPI_TASKS}

# Clean up previous test outputs, but preserve downloaded tree files
rm -rf output
rm -rf ${OUTPUTDIR}
rm -f *${RAWNAME}*
# Only remove input/*.par files, preserve input/millennium/trees/
rm -f input/*.par 2>/dev/null

# Create output directory early to avoid issues with SAGE output later
mkdir -p ${OUTPUTDIR}

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


# Call preparation function from benchmark_utils.sh
prepare_parameter_file

mkdir -p output/millennium/

echo "========== Run sage with mpirun -np $MPI_TASKS =========="


echo "No shared SAGE output found - running SAGE standalone"
mpirun -np $MPI_TASKS ${MY_ROOT}/bin/sage input/millennium_sage_hdf5.par

# Verify expected output files were created
echo "Verifying parallel SAGE output..."
for ((i=0; i<MPI_TASKS; i++)); do
    if [ ! -f "output/millennium/model_${i}.hdf5" ]; then
        echo "ERROR: Expected output file output/millennium/model_${i}.hdf5 not found!"
        echo "Contents of output/millennium:"
        ls -la output/millennium/ || true
        exit 1
    fi
done
echo "SUCCESS: Found all $MPI_TASKS model_<N>.hdf5 output files."

mv output/millennium ${OUTPUTDIR}/
mv output/log ${OUTPUTDIR}/ 2>/dev/null || true

