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

# Change to test data directory - required for relative paths like input/, output/, mypar_files/
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

if [ ! -f "${MY_ROOT}/bin/sage" ]; then
    echo "sage not found - rebuild needed"
    NEED_REBUILD=1
fi

if [ ! -f "${MY_ROOT}/bin/sage2kdtree" ]; then
    echo "sage2kdtree not found - rebuild needed"
    NEED_REBUILD=1
fi

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

REPORT_DIR=benchmark_reports
mkdir -p ${REPORT_DIR}
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
REPORT_FILE=${REPORT_DIR}/benchmark_${TIMESTAMP}.md

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

    if [ $FIRST_RUN_STATUS -ne 0 ]; then
        echo "ERROR: first_run.sh failed with exit code $FIRST_RUN_STATUS"
        exit 1
    fi

    if [ ! -f "input/millennium/trees/trees_063.7" ]; then
        echo "ERROR: Tree files not found. first_run.sh may have failed to download them."
        echo "Expected file: input/millennium/trees/trees_063.7"
        echo "Please check your internet connection and try again."
        exit 1
    fi

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

# PHASE 1: Run SAGE
echo "========== PHASE 1: Running SAGE =========="
run_with_profiling "sage" "$BENCHMARK_CSV" \
    ${MY_ROOT}/bin/sage input/millennium_sage_hdf5.par
mv output/millennium ${OUTPUTDIR}/
mv output/log ${OUTPUTDIR}/ 2>/dev/null || true

measure_disk_usage "${OUTPUTDIR}/millennium" "sage_output" "$DISK_IO_CSV"

# PHASE 2: sage2kdtree
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
python3 ${MY_ROOT}/src/plot_lightcone.py ${OUTPUTDIR}/$RAWNAME-lightcone.h5 SnapNum

# Clean up
rm -f log.00000
rm -rf log

# Generate timestamped Markdown report
export OUTPUTDIR
export BENCHMARK_CSV
export DISK_IO_CSV
export REPORT_FILE
export TIMESTAMP

python3 - <<'EOF'
import csv
import os

def read_csv(filename):
    data = {}
    try:
        with open(filename, 'r') as f:
            for row in csv.DictReader(f):
                data[row['phase']] = {
                    'duration': float(row['duration_sec']),
                    'peak_mem': float(row['peak_mem_mb'])
                }
    except FileNotFoundError:
        print(f"Warning: {filename} not found")
    except Exception as e:
        print(f"Error reading {filename}: {e}")
    return data

def read_disk_csv(filename):
    data = {}
    try:
        with open(filename, 'r') as f:
            for row in csv.DictReader(f):
                if row['disk_mb']:
                    data[row['phase']] = float(row['disk_mb'])
    except FileNotFoundError:
        print(f"Warning: {filename} not found")
    except Exception as e:
        print(f"Error reading {filename}: {e}")
    return data

outputdir   = os.environ['OUTPUTDIR']
report_file = os.environ['REPORT_FILE']
timestamp   = os.environ['TIMESTAMP']

timings = read_csv(f"{outputdir}/benchmark_timings.csv")
disk    = read_disk_csv(f"{outputdir}/benchmark_disk_io.csv")

with open(report_file, 'w') as f:
    f.write("# Workflow Benchmark Report\n\n")
    f.write(f"Generated: {timestamp}\n\n")

    if timings:
        total = sum(v['duration'] for v in timings.values())
        peak  = max(v['peak_mem'] for v in timings.values())
        f.write("## Summary\n\n")
        f.write(f"- **Total time**: {total:.2f} seconds\n")
        f.write(f"- **Peak memory**: {peak:.2f} MB\n\n")

    f.write("## Phase Timings\n\n")
    f.write("| Phase | Time (s) | Peak Mem (MB) |\n")
    f.write("|-------|----------|---------------|\n")
    for phase, v in timings.items():
        f.write(f"| {phase} | {v['duration']:.2f} | {v['peak_mem']:.2f} |\n")

    if disk:
        f.write("\n## Disk Usage\n\n")
        f.write("| Phase | Disk (MB) |\n")
        f.write("|-------|-----------|\n")
        for phase, mb in disk.items():
            f.write(f"| {phase} | {mb:.2f} |\n")

print(f"Report written to {report_file}")
EOF

echo ""
echo "=========================================="
echo "Benchmark complete."
echo "Results:"
echo "  Timings:  ${BENCHMARK_CSV}"
echo "  Disk I/O: ${DISK_IO_CSV}"
echo "  Report:   ${REPORT_FILE}"
echo "=========================================="
