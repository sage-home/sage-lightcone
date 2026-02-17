#!/bin/bash
# Benchmarking utilities for workflow comparison

# Start timing a phase
# Usage: start_timer "phase_name"
start_timer() {
    local phase_name=$1
    echo "$phase_name" > /tmp/benchmark_phase.txt
    date +%s.%N > /tmp/benchmark_start_${phase_name}.txt
}

# End timing and record results
# Usage: end_timer "phase_name" "output_file"
end_timer() {
    local phase_name=$1
    local output_file=$2
    local start_time=$(cat /tmp/benchmark_start_${phase_name}.txt)
    local end_time=$(date +%s.%N)
    
    echo "DEBUG: timer $phase_name start=$start_time end=$end_time"

    local duration="0"
    # Basic numeric check for start_time (simple regex for float/int)
    if [[ "$start_time" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
        duration=$(echo "$end_time - $start_time" | bc)
    fi

    echo "${phase_name},${duration}" >> "$output_file"
    rm /tmp/benchmark_start_${phase_name}.txt
}

# Run command with memory profiling
# Usage: run_with_profiling "phase_name" "output_csv" "command" "args..."
run_with_profiling() {
    local phase_name=$1
    local output_csv=$2
    shift 2
    # Do NOT capture cmd="$@" as string, use "$@" directly to preserve argument splitting

    echo "DEBUG: run_with_profiling phase='$phase_name' cmd='$@'"
    echo "DEBUG: current directory: $(pwd)"

    start_timer "$phase_name"

    # Use /usr/bin/time for memory profiling (works on macOS and Linux)
    if [[ "$OSTYPE" == "darwin"* ]]; then
        # macOS version
        /usr/bin/time -l "$@" 2> /tmp/benchmark_time_${phase_name}.txt
        local exit_code=$?

        # If command failed, print stderr
        if [ $exit_code -ne 0 ]; then
             echo "ERROR: Command failed with exit code $exit_code"
             echo "Captured stderr:"
             cat /tmp/benchmark_time_${phase_name}.txt
        fi

        # Extract peak memory (maxrss is in bytes on macOS)
        local peak_mem=$(grep "maximum resident set size" /tmp/benchmark_time_${phase_name}.txt | awk '{print $1}' | head -n 1)
        local peak_mem_mb="0"
        if [[ "$peak_mem" =~ ^[0-9]+$ ]]; then
            peak_mem_mb=$(echo "scale=2; $peak_mem / 1024 / 1024" | bc)
        fi
    else
        # Linux version
        /usr/bin/time -v "$@" 2> /tmp/benchmark_time_${phase_name}.txt
        local exit_code=$?

        # If command failed, print stderr
        if [ $exit_code -ne 0 ]; then
             echo "ERROR: Command failed with exit code $exit_code"
             echo "Captured stderr:"
             cat /tmp/benchmark_time_${phase_name}.txt
        fi

        echo "DEBUG: time output file content:"
        cat /tmp/benchmark_time_${phase_name}.txt

        # Extract peak memory (maxrss is in KB on Linux)
        local peak_mem=$(grep "Maximum resident set size" /tmp/benchmark_time_${phase_name}.txt | awk '{print $6}' | head -n 1)
        local peak_mem_mb="0"
        if [[ "$peak_mem" =~ ^[0-9]+$ ]]; then
            peak_mem_mb=$(echo "scale=2; $peak_mem / 1024" | bc)
        else
            echo "DEBUG: Could not parse peak_mem '$peak_mem'. defaulting to 0"
        fi
    fi

    end_timer "$phase_name" "$output_csv"

    # Append memory data to CSV
    local duration=$(tail -1 "$output_csv" | cut -d',' -f2)
    # Replace last line with phase, duration, peak_mem
    sed -i.bak "$ s/$/,${peak_mem_mb}/" "$output_csv"
    rm "${output_csv}.bak" 2>/dev/null

    rm /tmp/benchmark_time_${phase_name}.txt

    return $exit_code
}

# Measure disk usage of a directory
# Usage: measure_disk_usage "directory" "phase_name" "output_csv"
measure_disk_usage() {
    local dir=$1
    local phase_name=$2
    local output_csv=$3

    if [ -d "$dir" ]; then
        local size_mb=""
        
        if [[ "$OSTYPE" == "darwin"* ]]; then
            # macOS: du -k gives kilobytes
            local size_kb=$(du -sk "$dir" 2>/dev/null | awk '{print $1}')
            if [[ "$size_kb" =~ ^[0-9]+$ ]]; then
                size_mb=$(echo "scale=2; $size_kb / 1024" | bc)
            fi
        else
            # Linux: du -sb gives bytes
            local size_bytes=$(du -sb "$dir" 2>/dev/null | cut -f1)
            if [[ "$size_bytes" =~ ^[0-9]+$ ]]; then
                size_mb=$(echo "scale=2; $size_bytes / 1024 / 1024" | bc)
            fi
        fi

        if [ -n "$size_mb" ]; then
            echo "${phase_name},${size_mb}" >> "$output_csv"
        fi
    fi
}

# Count files in directory
# Usage: count_files "directory"
count_files() {
    local dir=$1
    if [ -d "$dir" ]; then
        find "$dir" -type f | wc -l
    else
        echo "0"
    fi
}

# Prepare parameter files for SAGE run (downloads trees, updates paths, generates par files)
# Usage: prepare_parameter_file
prepare_parameter_file() {
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

    # Update paths in millennium.par to use current directory
    CURRENT_DIR=$(pwd)
    sed -i'' -e "s|^OutputDir.*|OutputDir   ${CURRENT_DIR}/output/millennium/|" input/millennium.par
    sed -i'' -e "s|^SimulationDir.*|SimulationDir   ${CURRENT_DIR}/input/millennium/trees/|" input/millennium.par
    sed -i'' -e "s|^FileWithSnapList.*|FileWithSnapList ${CURRENT_DIR}/input/millennium/trees/millennium.a_list|" input/millennium.par
    echo "✓ Updated paths in millennium.par"

    # Extract settings from millennium.par
    echo "Extracting settings from millennium.par..."
    python3 utils/extract_settings.py

    # Generate .par files by concatenating headers with settings
    cat mypar_files/millennium_sage_binary_header.txt mypar_files/millennium_settings.txt > input/millennium.par
    cat mypar_files/millennium_sage_binary_kdtreeindex_header.txt mypar_files/millennium_settings.txt > input/millennium_minus1.par
    cat mypar_files/millennium_sage_hdf5_header.txt mypar_files/millennium_settings.txt > input/millennium_sage_hdf5.par
    echo "✓ Parameter files generated."
}
