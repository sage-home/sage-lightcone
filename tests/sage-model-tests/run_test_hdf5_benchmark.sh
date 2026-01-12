#!/bin/bash

# Load utilities
source $(dirname $0)/utils/benchmark_utils.sh

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

export RAWNAME=myhdf5millennium
export OUTPUTDIR=output_sage_hdf5_benchmark
export BENCHMARK_CSV=${OUTPUTDIR}/benchmark_timings.csv
export DISK_IO_CSV=${OUTPUTDIR}/benchmark_disk_io.csv

# Clean and prepare
rm -rf input
rm -rf output
rm -rf ${OUTPUTDIR}
rm -f *${RAWNAME}*

# Create output directory early for CSV files
mkdir -p ${OUTPUTDIR}

# Initialize CSV files
echo "phase,duration_sec,peak_mem_mb" > $BENCHMARK_CSV
echo "phase,disk_mb" > $DISK_IO_CSV

# Get tree files
./first_run.sh
echo ==== Have finished with first_run.sh ====

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

# PHASE 2: sageh5toh5
echo "========== PHASE 2: Running sageh5toh5 =========="
run_with_profiling "sageh5toh5" "$BENCHMARK_CSV" \
    ${MY_ROOT}/bin/sageh5toh5 -m convert -v 2 \
    -s ${OUTPUTDIR}/millennium \
    -p input/millennium_sage_hdf5.par \
    -a input/millennium/trees/millennium.a_list \
    -o ${OUTPUTDIR}/${RAWNAME}-depthfirstordered.h5

measure_disk_usage "${OUTPUTDIR}" "after_sageh5toh5" "$DISK_IO_CSV"

# PHASE 3: sageimport
echo "========== PHASE 3: Running sageimport =========="
run_with_profiling "sageimport" "$BENCHMARK_CSV" \
    ${MY_ROOT}/bin/sageimport \
    --settings ${OUTPUTDIR}/${RAWNAME}-depthfirstordered_import_settings.xml

measure_disk_usage "${OUTPUTDIR}" "after_sageimport" "$DISK_IO_CSV"

# PHASE 4: dstreeinit
echo "========== PHASE 4: Running dstreeinit =========="
run_with_profiling "dstreeinit" "$BENCHMARK_CSV" \
    ${MY_ROOT}/bin/dstreeinit --mode kdtree \
    --tree ${OUTPUTDIR}/out_$RAWNAME-depthfirstordered.h5 \
    --sage ${OUTPUTDIR}/$RAWNAME-bysnap.h5 \
    --output ${OUTPUTDIR}/$RAWNAME-kdtree.h5

measure_disk_usage "${OUTPUTDIR}" "after_dstreeinit" "$DISK_IO_CSV"

# PHASE 5: cli_lightcone
echo "========== PHASE 5: Running cli_lightcone =========="
run_with_profiling "cli_lightcone" "$BENCHMARK_CSV" \
    ${MY_ROOT}/bin/cli_lightcone \
    --dataset ${OUTPUTDIR}/$RAWNAME-kdtree.h5 \
    --decmin 0 --decmax 1 --ramin 0 --ramax 1 --zmin 0 --zmax 1 \
    --outdir ${OUTPUTDIR} --outfile $RAWNAME-lightcone.h5

measure_disk_usage "${OUTPUTDIR}" "final" "$DISK_IO_CSV"

# PHASE 6: Plot (not benchmarked - visualization only)
echo "========== PHASE 6: Plotting =========="
source ${MY_ROOT}/.venv/bin/activate
python3 ${MY_ROOT}/src/plot_lightcone.py ${OUTPUTDIR}/$RAWNAME-lightcone.h5 snapnum

# Clean up
rm -f log.00000
rm -rf log

echo ""
echo "=========================================="
echo "OLD workflow benchmark complete."
echo "Results:"
echo "  Timings: ${BENCHMARK_CSV}"
echo "  Disk I/O: ${DISK_IO_CSV}"
echo "=========================================="
