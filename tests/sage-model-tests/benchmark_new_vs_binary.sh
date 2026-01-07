#!/bin/bash

# Benchmark NEW workflow against Binary workflow baseline
# This compares the NEW HDF5-based workflow against the traditional binary SAGE workflow

# What's my code's root directory
export MY_SCRIPT=${BASH_SOURCE:-$0}
export MY_SCRIPTS_DIRECTORY=$(dirname $MY_SCRIPT)
export MY_ROOT=$(cd ${MY_SCRIPTS_DIRECTORY}/../.. && pwd)

echo "=========================================="
echo "  NEW vs BINARY Workflow Benchmark"
echo "=========================================="
echo ""

# Parse arguments
RUN_BINARY=1
RUN_NEW=1
CLEANUP=1

for arg in "$@"; do
    case $arg in
        --binary-only)
            RUN_NEW=0
            ;;
        --new-only)
            RUN_BINARY=0
            ;;
        --no-cleanup)
            CLEANUP=0
            ;;
        --help|-h)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Compare NEW workflow (HDF5 SAGE → sage2kdtree) against"
            echo "BINARY workflow (Binary SAGE → sage2h5 → sageimport → dstreeinit)"
            echo ""
            echo "Options:"
            echo "  --binary-only  Only run BINARY workflow benchmark"
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
REPORT_FILE=${REPORT_DIR}/new_vs_binary_${TIMESTAMP}.md

# NOTE: BINARY and NEW workflows require different SAGE output formats:
# - BINARY workflow: sage2h5 requires binary SAGE output (millennium.par)
# - NEW workflow: sage2kdtree processes HDF5 SAGE output (millennium_sage_hdf5.par)
# Unfortunately, we cannot share SAGE output between these workflows.
# Due to SAGE non-determinism (see SAGE_REPRODUCIBILITY_ISSUE.md), field values
# may differ slightly between runs, but kdtree structure should match.

# Prepare environment
rm -rf ${MY_SCRIPTS_DIRECTORY}/input
rm -rf ${MY_SCRIPTS_DIRECTORY}/output

# Get tree files
pushd ${MY_SCRIPTS_DIRECTORY} > /dev/null
./first_run.sh
popd > /dev/null

# Copy parameter files
mkdir -p ${MY_SCRIPTS_DIRECTORY}/input
cp ${MY_SCRIPTS_DIRECTORY}/mypar_files/millennium.par ${MY_SCRIPTS_DIRECTORY}/input/
cp ${MY_SCRIPTS_DIRECTORY}/mypar_files/millennium_minus1.par ${MY_SCRIPTS_DIRECTORY}/input/
cp ${MY_SCRIPTS_DIRECTORY}/mypar_files/millennium_sage_hdf5.par ${MY_SCRIPTS_DIRECTORY}/input/

# Run BINARY workflow (for baseline comparison)
if [ $RUN_BINARY -eq 1 ]; then
    echo "=========================================="
    echo "Running BINARY workflow baseline"
    echo "=========================================="
    pushd ${MY_SCRIPTS_DIRECTORY} > /dev/null

    export RAWNAME=mybinarymillennium
    export OUTPUTDIR=output_sage_binary_benchmark

    # Clean output directory
    rm -rf ${OUTPUTDIR}
    mkdir -p ${OUTPUTDIR}

    # Run SAGE with BINARY output format
    # SAGE writes to ./output/millennium/ as specified in parameter file
    echo "Running SAGE with binary output format..."
    mkdir -p output/millennium
    ${MY_ROOT}/bin/sage input/millennium.par
    SAGE_STATUS=$?

    if [ $SAGE_STATUS -ne 0 ]; then
        echo "ERROR: SAGE (binary) execution failed!"
        exit 1
    fi

    # Move SAGE output to benchmark directory
    mv output/millennium ${OUTPUTDIR}/

    # Run the conversion pipeline
    echo "Running sage2h5..."
    ${MY_ROOT}/bin/sage2h5 -m convert -v 2 \
        -s ${OUTPUTDIR}/millennium \
        -p input/millennium_minus1.par \
        -a input/millennium/trees/millennium.a_list \
        -o ${OUTPUTDIR}/${RAWNAME}-depthfirstordered.h5

    echo "Running sageimport..."
    ${MY_ROOT}/bin/sageimport --settings ${OUTPUTDIR}/${RAWNAME}-depthfirstordered_import_settings.xml

    echo "Running dstreeinit..."
    ${MY_ROOT}/bin/dstreeinit --mode kdtree \
        --tree ${OUTPUTDIR}/out_${RAWNAME}-depthfirstordered.h5 \
        --sage ${OUTPUTDIR}/${RAWNAME}-bysnap.h5 \
        --output ${OUTPUTDIR}/${RAWNAME}-kdtree.h5

    BINARY_STATUS=$?
    popd > /dev/null

    if [ $BINARY_STATUS -ne 0 ]; then
        echo "ERROR: BINARY workflow failed!"
        exit 1
    fi
    echo "✓ BINARY workflow complete"
    echo ""
fi

# Run NEW workflow
if [ $RUN_NEW -eq 1 ]; then
    echo "=========================================="
    echo "Running NEW workflow"
    echo "=========================================="
    pushd ${MY_SCRIPTS_DIRECTORY} > /dev/null

    export RAWNAME=myhdf5millennium
    export OUTPUTDIR=output_sage_hdf5_new_benchmark

    # Clean output directory
    rm -rf ${OUTPUTDIR}
    mkdir -p ${OUTPUTDIR}

    # Run SAGE with HDF5 output format
    # SAGE writes to ./output/millennium/ as specified in parameter file
    echo "Running SAGE with HDF5 output format..."
    mkdir -p output/millennium
    ${MY_ROOT}/bin/sage input/millennium_sage_hdf5.par
    SAGE_STATUS=$?

    if [ $SAGE_STATUS -ne 0 ]; then
        echo "ERROR: SAGE (HDF5) execution failed!"
        exit 1
    fi

    # Move SAGE output to benchmark directory
    mv output/millennium ${OUTPUTDIR}/

    # Run sage2kdtree (single-step conversion)
    echo "Running sage2kdtree..."
    ${MY_ROOT}/bin/sage2kdtree \
        -s ${OUTPUTDIR}/millennium \
        -p input/millennium_sage_hdf5.par \
        -a input/millennium/trees/millennium.a_list \
        -o ${OUTPUTDIR}/${RAWNAME}-kdtree-new.h5 \
        --ppc 1000 -v 2

    NEW_STATUS=$?
    popd > /dev/null

    if [ $NEW_STATUS -ne 0 ]; then
        echo "ERROR: NEW workflow failed!"
        exit 1
    fi
    echo "✓ NEW workflow complete"
    echo ""
fi

# Compare outputs if both ran
if [ $RUN_BINARY -eq 1 ] && [ $RUN_NEW -eq 1 ]; then
    echo "=========================================="
    echo "Generating test lightcones for validation"
    echo "=========================================="

    # Generate lightcone from BINARY kdtree
    echo "Generating BINARY lightcone..."
    ${MY_ROOT}/bin/cli_lightcone \
        --dataset ${MY_SCRIPTS_DIRECTORY}/output_sage_binary_benchmark/mybinarymillennium-kdtree.h5 \
        --decmin 0 --decmax 1 --ramin 0 --ramax 1 --zmin 0 --zmax 1 \
        --outdir ${MY_SCRIPTS_DIRECTORY}/output_sage_binary_benchmark \
        --outfile mybinarymillennium-lightcone.h5

    # Generate lightcone from NEW kdtree
    echo "Generating NEW lightcone..."
    ${MY_ROOT}/bin/cli_lightcone \
        --dataset ${MY_SCRIPTS_DIRECTORY}/output_sage_hdf5_new_benchmark/myhdf5millennium-kdtree-new.h5 \
        --decmin 0 --decmax 1 --ramin 0 --ramax 1 --zmin 0 --zmax 1 \
        --outdir ${MY_SCRIPTS_DIRECTORY}/output_sage_hdf5_new_benchmark \
        --outfile myhdf5millennium-lightcone-new.h5

    echo ""
    echo "Comparing workflow outputs..."
    python3 ${MY_SCRIPTS_DIRECTORY}/utils/validate_outputs.py \
        --old-kdtree ${MY_SCRIPTS_DIRECTORY}/output_sage_binary_benchmark/mybinarymillennium-kdtree.h5 \
        --new-kdtree ${MY_SCRIPTS_DIRECTORY}/output_sage_hdf5_new_benchmark/myhdf5millennium-kdtree-new.h5 \
        --old-lightcone ${MY_SCRIPTS_DIRECTORY}/output_sage_binary_benchmark/mybinarymillennium-lightcone.h5 \
        --new-lightcone ${MY_SCRIPTS_DIRECTORY}/output_sage_hdf5_new_benchmark/myhdf5millennium-lightcone-new.h5 \
        --output ${REPORT_DIR}/new_vs_binary_validation_${TIMESTAMP}.txt

    VALIDATION_STATUS=$?
    if [ $VALIDATION_STATUS -eq 0 ]; then
        echo "✓ Validation passed"
    else
        echo "✗ Validation failed - see ${REPORT_DIR}/new_vs_binary_validation_${TIMESTAMP}.txt"
    fi
fi

# Generate comparison report
echo "Generating comparison report..."

cat > ${REPORT_FILE} << EOF
# NEW vs BINARY Workflow Comparison

Generated: ${TIMESTAMP}

## Workflow Descriptions

### BINARY Workflow (Baseline)
\`\`\`
SAGE (HDF5) → sage2h5 → sageimport → dstreeinit → kdtree
\`\`\`
- Traditional multi-step conversion
- Processes SAGE HDF5 through intermediate formats
- Uses compound datatypes then splits to columnar

### NEW Workflow
\`\`\`
SAGE (HDF5) → sage2kdtree → kdtree
\`\`\`
- Direct single-step conversion
- Skips intermediate file formats
- Writes directly to columnar kdtree format

## Output Files

- **BINARY kdtree**: \`output_sage_binary_benchmark/mybinarymillennium-kdtree.h5\`
- **NEW kdtree**: \`output_sage_hdf5_new_benchmark/myhdf5millennium-kdtree-new.h5\`

## Validation Results

See: \`benchmark_reports/new_vs_binary_validation_${TIMESTAMP}.txt\`

## Key Differences

### Field Names
Both workflows should now produce identical field names matching SAGE output:
- SAGE fields: CamelCase (e.g., \`Posx\`, \`StellarMass\`, \`Mvir\`)
- Computed fields: lowercase_with_underscores (e.g., \`global_index\`, \`subsize\`)

### Units Preservation
Both workflows now preserve SAGE units:
- \`dT\`: Myr (megayears) - no Gyr conversion
- All other fields: original SAGE units

### Data Consistency
**IMPORTANT**: SAGE runs twice (once with binary output, once with HDF5 output).
Due to SAGE non-determinism (see SAGE_REPRODUCIBILITY_ISSUE.md), field values
may differ slightly between the two SAGE runs, but kdtree structures should match.

Differences in field values may reflect:
1. SAGE non-determinism (different output formats, execution ordering)
2. Pipeline-specific field processing
3. Unit conversions (should be identical after recent fixes)

## Notes

- SAGE cannot be shared between workflows (requires different output formats)
- BINARY workflow: Uses binary SAGE output (millennium.par)
- NEW workflow: Uses HDF5 SAGE output (millennium_sage_hdf5.par)
- Field names should match SAGE output (CamelCase for SAGE fields)
- Units should be preserved from SAGE (no Myr→Gyr conversions)
EOF

echo "Report written to ${REPORT_FILE}"

echo ""
echo "=========================================="
echo "  Comparison Complete"
echo "=========================================="
echo ""
echo "Results:"
echo "  - Comparison report: ${REPORT_FILE}"
if [ $RUN_BINARY -eq 1 ] && [ $RUN_NEW -eq 1 ]; then
    echo "  - Validation report: ${REPORT_DIR}/new_vs_binary_validation_${TIMESTAMP}.txt"
fi
echo ""

# Cleanup if requested
if [ $CLEANUP -eq 1 ]; then
    echo "Note: Intermediate files preserved for analysis."
    echo "      Use --no-cleanup to keep all outputs."
fi
