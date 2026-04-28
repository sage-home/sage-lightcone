#!/bin/bash
# What's my code's root directory
SCRIPT="${BASH_SOURCE[0]}"
[ -z "$SCRIPT" ] && SCRIPT="$0"
export MY_SCRIPTS_DIRECTORY=$(cd "$(dirname "$SCRIPT")" && pwd)
if [ ${MY_ROOT} ];then
    echo Users Code ${MY_ROOT}
else
    MY_ROOT="$(cd "${_SCRIPTS_DIR}/.." && pwd)"
    echo Default Code ${MY_ROOT}
fi

if [ ${SIMULATION} ];then
    echo Users Simulation ${SIMULATION}
else
    SIMULATION=millennium
    echo Default Simulation ${SIMULATION}
fi


# Load utilities
source ${MY_ROOT}/tests/utils/benchmark_utils.sh

# Exit on error
set -e

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

export OUTPUTDIR=../DB/${SIMULATION}
export BENCHMARK_CSV=${OUTPUTDIR}/lightcone_benchmark_timings.csv
export DISK_IO_CSV=${OUTPUTDIR}/lightcone_benchmark_disk_io.csv

REPORT_DIR=../REPORTS/benchmark_reports
mkdir -p ${REPORT_DIR}
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
REPORT_FILE=${REPORT_DIR}/benchmark_${TIMESTAMP}.md

# Clean up previous test outputs and recreate them
#rm -rf output
#rm -rf ${OUTPUTDIR}
#rm -f ${BENCHMARK_CSV}
#rm -f ${DISK_IO_CSV}

mkdir -p ${OUTPUTDIR}
mkdir -p $(dirname "${BENCHMARK_CSV}")
mkdir -p $(dirname "${DISK_IO_CSV}")

# Initialize CSV files
echo "phase,duration_sec,peak_mem_mb" > $BENCHMARK_CSV
echo "phase,disk_mb" > $DISK_IO_CSV

mkdir -p input
if [ -f "par_by_simulation/${SIMULATION}.par" ]; then
    cp "par_by_simulation/${SIMULATION}.par" input/${SIMULATION}.par
    echo "✓ Copied base ${SIMULATION}.par from sage-model"
else
    echo "ERROR: par_by_simulation/${SIMULATION}.par not found"
    exit 1
fi

mkdir -p output/${SIMULATION}/

# PHASE 1: Run SAGE
echo "========== PHASE 1: Skip Running SAGE =========="
# PHASE 2: sage2kdtree
echo "========== PHASE 2: Skip Running sage2kdtree =========="
# PHASE 3: cli_lightcone
echo "========== PHASE 3: Running cli_lightcone =========="
run_with_profiling "cli_lightcone" "${BENCHMARK_CSV}" \
    ${MY_ROOT}/bin/cli_lightcone \
    --dataset ../DB/${SIMULATION}/${SIMULATION}-kdtree.h5 \
    --decmin 0 --decmax 1 --ramin 0 --ramax 1 --zmin 0 --zmax 1 \
    --outdir ${OUTPUTDIR} --outfile $SIMULATION-lightcone.h5

measure_disk_usage "${OUTPUTDIR}" "final" "${DISK_IO_CSV}"

# PHASE 4: Plot (not benchmarked - visualization only)
echo "========== PHASE 4: Plotting =========="
source ${MY_ROOT}/.venv/bin/activate
python3 ${MY_ROOT}/src/plot_lightcone.py ${OUTPUTDIR}/$SIMULATION-lightcone.h5 SnapNum
mv lightcone_3d_SnapNum.png ${OUTPUTDIR}/${SIMULATION}_3d_SnapNum.png

# Clean up
#rm -rf output
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
