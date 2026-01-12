#!/bin/bash

# What's my code's root directory
SCRIPT="${BASH_SOURCE[0]}"
[ -z "$SCRIPT" ] && SCRIPT="$0"
export MY_SCRIPTS_DIRECTORY=$(cd "$(dirname "$SCRIPT")" && pwd)
export MY_ROOT=$(cd "${MY_SCRIPTS_DIRECTORY}/../.." && pwd)

echo "=========================================="
echo "  Workflow Benchmarking Suite"
echo "=========================================="
echo ""

# Parse arguments
RUN_OLD=1
RUN_NEW=1
CLEANUP=1

for arg in "$@"; do
    case $arg in
        --old-only)
            RUN_NEW=0
            ;;
        --new-only)
            RUN_OLD=0
            ;;
        --no-cleanup)
            CLEANUP=0
            ;;
        --help|-h)
            echo "Usage: $0 [OPTIONS]"
            echo "Options:"
            echo "  --old-only     Only run OLD workflow benchmark"
            echo "  --new-only     Only run NEW workflow benchmark"
            echo "  --no-cleanup   Keep intermediate files for inspection"
            echo "  --help, -h     Show this help message"
            exit 0
            ;;
    esac
done

# Create output directory for reports
REPORT_DIR=${MY_SCRIPTS_DIRECTORY}/benchmark_reports
mkdir -p ${REPORT_DIR}
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
REPORT_FILE=${REPORT_DIR}/benchmark_${TIMESTAMP}.md

# Run SAGE once to ensure both workflows use IDENTICAL data
# (See SAGE_REPRODUCIBILITY_ISSUE.md for why this is necessary)
if [ $RUN_OLD -eq 1 ] || [ $RUN_NEW -eq 1 ]; then
    echo "=========================================="
    echo "Running SAGE (shared for both workflows)"
    echo "=========================================="

    # Clean and prepare shared SAGE output
    rm -rf ${MY_SCRIPTS_DIRECTORY}/shared_sage_output
    rm -rf ${MY_SCRIPTS_DIRECTORY}/input
    rm -rf ${MY_SCRIPTS_DIRECTORY}/output

    # Get tree files
    pushd ${MY_SCRIPTS_DIRECTORY} > /dev/null
    python3 utils/extract_settings.py
    ./first_run.sh
    popd > /dev/null

    # Copy parameter files
    mkdir -p ${MY_SCRIPTS_DIRECTORY}/input
    cat ${MY_SCRIPTS_DIRECTORY}/mypar_files/millennium_sage_binary_header.txt ${MY_SCRIPTS_DIRECTORY}/mypar_files/millennium_settings.txt >> ${MY_SCRIPTS_DIRECTORY}/input/millennium.par
    cat ${MY_SCRIPTS_DIRECTORY}/mypar_files/millennium_sage_binary_kdtreeindex_header.txt ${MY_SCRIPTS_DIRECTORY}/mypar_files/millennium_settings.txt >> ${MY_SCRIPTS_DIRECTORY}/input/millennium_minus1.par
    cat ${MY_SCRIPTS_DIRECTORY}/mypar_files/millennium_sage_hdf5_header.txt ${MY_SCRIPTS_DIRECTORY}/mypar_files/millennium_settings.txt >> ${MY_SCRIPTS_DIRECTORY}/input/millennium_sage_hdf5.par

    # Run SAGE once
    mkdir -p ${MY_SCRIPTS_DIRECTORY}/output/millennium
    ${MY_ROOT}/bin/sage ${MY_SCRIPTS_DIRECTORY}/input/millennium_sage_hdf5.par
    SAGE_STATUS=$?

    if [ $SAGE_STATUS -ne 0 ]; then
        echo "ERROR: SAGE execution failed!"
        exit 1
    fi

    # Save shared SAGE output
    cp -r ${MY_SCRIPTS_DIRECTORY}/output/millennium ${MY_SCRIPTS_DIRECTORY}/shared_sage_output
    echo "✓ SAGE complete - output saved to shared_sage_output/"
    echo ""
fi

# Run OLD workflow
if [ $RUN_OLD -eq 1 ]; then
    echo "Running OLD workflow benchmark..."
    pushd ${MY_SCRIPTS_DIRECTORY} > /dev/null
    ./run_test_hdf5_benchmark.sh
    OLD_STATUS=$?
    popd > /dev/null

    if [ $OLD_STATUS -ne 0 ]; then
        echo "ERROR: OLD workflow failed!"
        exit 1
    fi
    echo "✓ OLD workflow complete"
fi

# Run NEW workflow
if [ $RUN_NEW -eq 1 ]; then
    echo "Running NEW workflow benchmark..."
    pushd ${MY_SCRIPTS_DIRECTORY} > /dev/null
    ./run_test_hdf5_one_step_benchmark.sh
    NEW_STATUS=$?
    popd > /dev/null

    if [ $NEW_STATUS -ne 0 ]; then
        echo "ERROR: NEW workflow failed!"
        exit 1
    fi
    echo "✓ NEW workflow complete"
fi

# Validate outputs if both ran
if [ $RUN_OLD -eq 1 ] && [ $RUN_NEW -eq 1 ]; then
    echo "Validating workflow outputs..."
    
    # Activate virtual environment if available
    if [ -f "${MY_ROOT}/.venv/bin/activate" ]; then
        source "${MY_ROOT}/.venv/bin/activate"
    elif [ -f "${MY_ROOT}/setup_mac.sh" ] && [[ "$OSTYPE" == "darwin"* ]]; then
        # Try to source setup script which activates venv
        source "${MY_ROOT}/setup_mac.sh" > /dev/null 2>&1
    fi

    python3 ${MY_SCRIPTS_DIRECTORY}/utils/validate_outputs.py \
        --old-kdtree ${MY_SCRIPTS_DIRECTORY}/output_sage_hdf5_benchmark/myhdf5millennium-kdtree.h5 \
        --new-kdtree ${MY_SCRIPTS_DIRECTORY}/output_sage_hdf5_one_step_benchmark/myhdf5millennium-kdtree-onestep.h5 \
        --old-lightcone ${MY_SCRIPTS_DIRECTORY}/output_sage_hdf5_benchmark/myhdf5millennium-lightcone.h5 \
        --new-lightcone ${MY_SCRIPTS_DIRECTORY}/output_sage_hdf5_one_step_benchmark/myhdf5millennium-lightcone.h5 \
        --output ${REPORT_DIR}/validation_${TIMESTAMP}.txt

    VALIDATION_STATUS=$?
    if [ $VALIDATION_STATUS -eq 0 ]; then
        echo "✓ Validation passed"
    else
        echo "✗ Validation failed - see ${REPORT_DIR}/validation_${TIMESTAMP}.txt"
    fi
fi

# Generate comparison report
echo "Generating benchmark report..."

# Export variables for Python script
export MY_SCRIPTS_DIRECTORY
export REPORT_FILE
export TIMESTAMP
export RUN_OLD
export RUN_NEW

python3 - <<EOF
import csv
import sys
import os

def read_csv(filename):
    data = {}
    try:
        with open(filename, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
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
            reader = csv.DictReader(f)
            for row in reader:
                # Handle empty values
                if row['disk_mb']:
                    data[row['phase']] = float(row['disk_mb'])
    except FileNotFoundError:
        print(f"Warning: {filename} not found")
    except Exception as e:
        print(f"Error reading {filename}: {e}")
    return data

scripts_dir = os.environ.get('MY_SCRIPTS_DIRECTORY', '.')
report_file = os.environ.get('REPORT_FILE', 'benchmark_report.md')
timestamp = os.environ.get('TIMESTAMP', 'unknown')
run_old = int(os.environ.get('RUN_OLD', '0'))
run_new = int(os.environ.get('RUN_NEW', '0'))

old_timings = {}
new_timings = {}
old_disk = {}
new_disk = {}

if run_old:
    old_timings = read_csv(f'{scripts_dir}/output_sage_hdf5_benchmark/benchmark_timings.csv')
    old_disk = read_disk_csv(f'{scripts_dir}/output_sage_hdf5_benchmark/benchmark_disk_io.csv')

if run_new:
    new_timings = read_csv(f'{scripts_dir}/output_sage_hdf5_one_step_benchmark/benchmark_timings.csv')
    new_disk = read_disk_csv(f'{scripts_dir}/output_sage_hdf5_one_step_benchmark/benchmark_disk_io.csv')

with open(report_file, 'w') as f:
    f.write("# Workflow Benchmark Report\n\n")
    f.write(f"Generated: {timestamp}\n\n")

    if run_old and run_new:
        f.write("## Executive Summary\n\n")

        # Calculate total times
        old_total = sum(v['duration'] for v in old_timings.values())
        new_total = sum(v['duration'] for v in new_timings.values())
        speedup = (old_total - new_total) / old_total * 100 if old_total > 0 else 0

        f.write(f"- **OLD Workflow Total Time**: {old_total:.2f} seconds\n")
        f.write(f"- **NEW Workflow Total Time**: {new_total:.2f} seconds\n")
        f.write(f"- **Speedup**: {speedup:.1f}%\n\n")

        # Peak memory
        old_peak = max(v['peak_mem'] for v in old_timings.values()) if old_timings else 0
        new_peak = max(v['peak_mem'] for v in new_timings.values()) if new_timings else 0
        mem_reduction = (old_peak - new_peak) / old_peak * 100 if old_peak > 0 else 0

        f.write(f"- **OLD Workflow Peak Memory**: {old_peak:.2f} MB\n")
        f.write(f"- **NEW Workflow Peak Memory**: {new_peak:.2f} MB\n")
        f.write(f"- **Memory Reduction**: {mem_reduction:.1f}%\n\n")

        # Disk I/O (if available)
        old_disk_total = old_disk.get('final', 0)
        new_disk_total = new_disk.get('final', 0)
        if old_disk_total > 0 and new_disk_total > 0:
            io_reduction = (old_disk_total - new_disk_total) / old_disk_total * 100
            f.write(f"- **OLD Workflow Disk Usage**: {old_disk_total:.2f} MB\n")
            f.write(f"- **NEW Workflow Disk Usage**: {new_disk_total:.2f} MB\n")
            f.write(f"- **Disk I/O Reduction**: {io_reduction:.1f}%\n\n")

    f.write("## Detailed Timing Breakdown\n\n")
    f.write("| Phase | OLD Time (s) | NEW Time (s) | OLD Mem (MB) | NEW Mem (MB) |\n")
    f.write("|-------|--------------|--------------|--------------|---------------|\n")

    all_phases = sorted(set(old_timings.keys()) | set(new_timings.keys()))
    for phase in all_phases:
        old_time = old_timings.get(phase, {}).get('duration', 0)
        new_time = new_timings.get(phase, {}).get('duration', 0)
        old_mem = old_timings.get(phase, {}).get('peak_mem', 0)
        new_mem = new_timings.get(phase, {}).get('peak_mem', 0)

        old_str = f"{old_time:.2f}" if old_time > 0 else "N/A"
        new_str = f"{new_time:.2f}" if new_time > 0 else "N/A"
        old_mem_str = f"{old_mem:.2f}" if old_mem > 0 else "N/A"
        new_mem_str = f"{new_mem:.2f}" if new_mem > 0 else "N/A"

        f.write(f"| {phase} | {old_str} | {new_str} | {old_mem_str} | {new_mem_str} |\n")

    # Only show disk I/O table if we have data
    if old_disk or new_disk:
        f.write("\n## Disk I/O Breakdown\n\n")
        f.write("| Phase | OLD Disk (MB) | NEW Disk (MB) |\n")
        f.write("|-------|---------------|---------------|\n")

        all_disk_phases = sorted(set(old_disk.keys()) | set(new_disk.keys()))
        for phase in all_disk_phases:
            old_d = old_disk.get(phase, 0)
            new_d = new_disk.get(phase, 0)

            old_str = f"{old_d:.2f}" if old_d > 0 else "N/A"
            new_str = f"{new_d:.2f}" if new_d > 0 else "N/A"

            f.write(f"| {phase} | {old_str} | {new_str} |\n")

print(f"Report written to {report_file}")
EOF

echo ""
echo "=========================================="
echo "  Benchmarking Complete"
echo "=========================================="
echo ""
echo "Results:"
echo "  - Benchmark report: ${REPORT_FILE}"
if [ $RUN_OLD -eq 1 ] && [ $RUN_NEW -eq 1 ]; then
    echo "  - Validation report: ${REPORT_DIR}/validation_${TIMESTAMP}.txt"
fi
echo ""

# Cleanup if requested
if [ $CLEANUP -eq 1 ]; then
    echo "Note: Intermediate files preserved for analysis."
    echo "      Use --no-cleanup to keep all outputs."
fi
