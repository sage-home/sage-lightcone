#!/bin/bash

# Load utilities
source $(dirname $0)/utils/benchmark_utils.sh

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
export MY_ROOT=$(cd "${MY_SCRIPTS_DIRECTORY}/../.." && pwd)

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

# Check if tree files already exist - skip first_run.sh if so
if [ -f "input/millennium/trees/trees_063.7" ] && [ -f "input/millennium/trees/millennium.a_list" ]; then
    echo "✓ Tree files already present - skipping download."
else
    echo "Tree files not found - running first_run.sh to download..."
    ./first_run.sh
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

mkdir -p output/millennium/
cat ${MY_SCRIPTS_DIRECTORY}/mypar_files/millennium_sage_binary_header.txt ${MY_SCRIPTS_DIRECTORY}/mypar_files/millennium_settings.txt >> ${MY_SCRIPTS_DIRECTORY}/input/millennium.par
cat ${MY_SCRIPTS_DIRECTORY}/mypar_files/millennium_sage_binary_kdtreeindex_header.txt ${MY_SCRIPTS_DIRECTORY}/mypar_files/millennium_settings.txt >> ${MY_SCRIPTS_DIRECTORY}/input/millennium_minus1.par
cat ${MY_SCRIPTS_DIRECTORY}/mypar_files/millennium_sage_hdf5_header.txt ${MY_SCRIPTS_DIRECTORY}/mypar_files/millennium_settings.txt >> ${MY_SCRIPTS_DIRECTORY}/input/millennium_sage_hdf5.par

# PHASE 1: Use shared SAGE output (run by benchmark_workflows.sh)
# NOTE: SAGE is run once by benchmark_workflows.sh and shared between workflows
# to ensure identical input data (see SAGE_REPRODUCIBILITY_ISSUE.md)
echo "========== PHASE 1: Using shared SAGE output =========="

# Check if shared SAGE output exists (from benchmark_workflows.sh)
if [ -d "../shared_sage_output" ]; then
    echo "Using shared SAGE output from benchmark_workflows.sh"
    cp -r ../shared_sage_output ${OUTPUTDIR}/millennium
elif [ -d "shared_sage_output" ]; then
    echo "Using shared SAGE output"
    cp -r shared_sage_output ${OUTPUTDIR}/millennium
else
    # Fallback: run SAGE if not run by benchmark_workflows.sh (for standalone testing)
    echo "No shared SAGE output found - running SAGE standalone"
    run_with_profiling "sage" "$BENCHMARK_CSV" \
        ${MY_ROOT}/bin/sage input/millennium_sage_hdf5.par
    mv output/millennium ${OUTPUTDIR}/
    mv output/log ${OUTPUTDIR}/ 2>/dev/null || true
fi

measure_disk_usage "${OUTPUTDIR}/millennium" "sage_output" "$DISK_IO_CSV"

# PHASE 2: sage2kdtree (replaces sageh5toh5 + sageimport + dstreeinit)
echo "========== PHASE 2: Running sage2kdtree =========="
run_with_profiling "sage2kdtree" "$BENCHMARK_CSV" \
    ${MY_ROOT}/bin/sage2kdtree \
    -s ${OUTPUTDIR}/millennium \
    -p input/millennium_sage_hdf5.par \
    -a input/millennium/trees/millennium.a_list \
    -o ${OUTPUTDIR}/${RAWNAME}-kdtree-onestep.h5 \
    --ppc 1000 -v 2

measure_disk_usage "${OUTPUTDIR}" "after_sage2kdtree" "$DISK_IO_CSV"

# PHASE 3: cli_lightcone
echo "========== PHASE 3: Running cli_lightcone =========="
run_with_profiling "cli_lightcone" "$BENCHMARK_CSV" \
    ${MY_ROOT}/bin/cli_lightcone \
    --dataset ${OUTPUTDIR}/${RAWNAME}-kdtree-onestep.h5 \
    --decmin 0 --decmax 1 --ramin 0 --ramax 1 --zmin 0 --zmax 1 \
    --outdir ${OUTPUTDIR} --outfile $RAWNAME-lightcone.h5

measure_disk_usage "${OUTPUTDIR}" "final" "$DISK_IO_CSV"

# PHASE 4: Plot (not benchmarked - visualization only)
echo "========== PHASE 4: Plotting =========="
source ${MY_ROOT}/.venv/bin/activate
# Field names are case-insensitive in plot_lightcone.py (SnapNum = snapnum)
python3 ${MY_ROOT}/src/plot_lightcone.py ${OUTPUTDIR}/$RAWNAME-lightcone.h5 SnapNum

# Clean up
rm -f log.00000
rm -rf log

echo ""
echo "=========================================="
echo "NEW workflow benchmark complete."
echo "Results:"
echo "  Timings: ${BENCHMARK_CSV}"
echo "  Disk I/O: ${DISK_IO_CSV}"
echo "=========================================="
