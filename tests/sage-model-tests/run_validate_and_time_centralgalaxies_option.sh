#!/bin/bash
#
# run_validate_and_time_centralgalaxies_option.sh
#
# Validates the --centralgalaxies option by running the full
# sage -> sage2kdtree -> cli_lightcone pipeline three times:
#
#   Pass 1: --centralgalaxies passed explicitly to sage2kdtree and cli_lightcone
#   Pass 2: DEFAULT_MODE="--centralgalaxies" passed via $DEFAULT_MODE to sage2kdtree only;
#           cli_lightcone gets no flag (mirrors the else-branch in run_test_sage_hdf5.sh)
#   Pass 3: DEFAULT_MODE="" passed to sage2kdtree (expands to nothing); cli_lightcone no flag
#
# Expected galaxy count ordering: pass1 > pass2 > pass3
#   Pass 1 extra galaxies (vs pass 2): satellites outside strict cone boundaries whose
#     central galaxy is inside the cone, included by cli_lightcone --centralgalaxies.
#   Galaxies strictly inside the cone are identical for pass 1 and pass 2.
#   Pass 2 extra galaxies (vs pass 3): central galaxy index built into kdtree gives
#     a richer query result even without --centralgalaxies in cli_lightcone.
#
# SAGE is run once and shared across all three passes.
# sage2kdtree and cli_lightcone are run and timed independently per pass.
# A comparison report is printed to stdout at the end.

# ---------------------------------------------------------------------------
# Timing helpers
# ---------------------------------------------------------------------------

now_sec() {
    python3 -c "import time; print(f'{time.time():.3f}')"
}

elapsed_sec() {
    local start=$1 end=$2
    python3 -c "print(f'{$end - $start:.1f}')"
}

# Count galaxies in a lightcone HDF5 output file.
# Prints three space-separated values: total inside outside
# where inside/outside are relative to the strict cone bounds.
count_galaxies_stats() {
    local file=$1 decmin=$2 decmax=$3 ramin=$4 ramax=$5 zmin=$6 zmax=$7
    if [ ! -f "$file" ]; then
        echo "MISSING MISSING MISSING"
        return
    fi
    python3 -c "
import h5py, sys
try:
    with h5py.File('$file', 'r') as f:
        keys = list(f.keys())
        ra_key  = next((k for k in keys if k.lower() in ['ra', 'rightascension', 'right_ascension']), None)
        dec_key = next((k for k in keys if k.lower() in ['dec', 'declination']), None)
        z_key   = next((k for k in keys if k.lower() in ['redshift_cosmological', 'redshift', 'z_cos', 'zcos']), None)
        if ra_key and dec_key and z_key:
            ra  = f[ra_key][:]
            dec = f[dec_key][:]
            z   = f[z_key][:]
            inside = ((ra  >= $ramin)  & (ra  <= $ramax)  &
                      (dec >= $decmin) & (dec <= $decmax) &
                      (z   >= $zmin)   & (z   <= $zmax))
            total    = len(ra)
            n_inside  = int(inside.sum())
            n_outside = total - n_inside
            print(f'{total} {n_inside} {n_outside}')
        else:
            # fallback: total only
            for key in keys:
                ds = f[key]
                if hasattr(ds, 'shape') and len(ds.shape) == 1 and ds.shape[0] > 0:
                    print(f'{ds.shape[0]} N/A N/A')
                    sys.exit(0)
            print('0 N/A N/A')
except Exception:
    print('N/A N/A N/A')
" 2>/dev/null
}

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------
# Locate root directory and load platform environment
# ---------------------------------------------------------------------------

export MY_SCRIPT="${BASH_SOURCE[0]:-$0}"
export MY_SCRIPTS_DIRECTORY=$(cd "$(dirname "$MY_SCRIPT")" && pwd)
export MY_ROOT=$(cd "${MY_SCRIPTS_DIRECTORY}/../.." && pwd)
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

export MY_SCRIPTS_DIRECTORY=$ORIGINAL_SCRIPTS_DIRECTORY
cd "${MY_SCRIPTS_DIRECTORY}"

# ---------------------------------------------------------------------------
# Ensure sage-model submodule is present
# ---------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------
# Build if needed
# ---------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

export RAWNAME=myhdf5millennium
export OUTPUTDIR=output_centralgalaxies_validation
export VALIDATIONDIR=validate_centralgalaxies_validation
mkdir -p ${VALIDATIONDIR}

# Lightcone query parameters (identical for all three passes)
LC_DECMIN=0
LC_DECMAX=1
LC_RAMIN=0
LC_RAMAX=1
LC_ZMIN=0
LC_ZMAX=1

# ---------------------------------------------------------------------------
# Prepare output directory and tree/parameter files
# ---------------------------------------------------------------------------

rm -rf output
rm -rf ${OUTPUTDIR}
rm -f *${RAWNAME}*
rm -f input/*.par 2>/dev/null

mkdir -p ${OUTPUTDIR}
mkdir -p input

if [ -f "${MY_ROOT}/sage-model/input/millennium.par" ]; then
    cp "${MY_ROOT}/sage-model/input/millennium.par" input/millennium.par
    echo "✓ Copied base millennium.par from sage-model"
else
    echo "ERROR: sage-model/input/millennium.par not found"
    echo "Please ensure the sage-model submodule is initialized: git submodule update --init --recursive"
    exit 1
fi

if [ -f "input/millennium/trees/trees_063.7" ] && [ -f "input/millennium/trees/millennium.a_list" ]; then
    echo "✓ Tree files already present - skipping download."
else
    echo "Tree files not found - running first_run.sh to download..."
    ./first_run.sh
    FIRST_RUN_STATUS=$?
    echo "==== Have finished with first_run.sh ===="
    if [ $FIRST_RUN_STATUS -ne 0 ]; then
        echo "ERROR: first_run.sh failed with exit code $FIRST_RUN_STATUS"
        exit 1
    fi
    if [ ! -f "input/millennium/trees/trees_063.7" ]; then
        echo "ERROR: Tree files not found after first_run.sh"
        exit 1
    fi
    if [ ! -f "input/millennium/trees/millennium.a_list" ]; then
        echo "ERROR: Scale factor list not found after first_run.sh"
        exit 1
    fi
    echo "✓ Tree files and scale factor list verified."
fi

CURRENT_DIR=$(pwd)
sed -i'' -e "s|^OutputDir.*|OutputDir   ${CURRENT_DIR}/output/millennium/|" input/millennium.par
sed -i'' -e "s|^SimulationDir.*|SimulationDir   ${CURRENT_DIR}/input/millennium/trees/|" input/millennium.par
sed -i'' -e "s|^FileWithSnapList.*|FileWithSnapList ${CURRENT_DIR}/input/millennium/trees/millennium.a_list|" input/millennium.par
echo "✓ Updated paths in millennium.par"

echo "Extracting settings from millennium.par..."
python3 utils/extract_settings.py

cat mypar_files/millennium_sage_binary_header.txt mypar_files/millennium_settings.txt > input/millennium.par
cat mypar_files/millennium_sage_binary_kdtreeindex_header.txt mypar_files/millennium_settings.txt > input/millennium_minus1.par
cat mypar_files/millennium_sage_hdf5_header.txt mypar_files/millennium_settings.txt > input/millennium_sage_hdf5.par
echo "✓ Parameter files generated."

mkdir -p output/millennium/

source ${MY_ROOT}/.venv/bin/activate

# ---------------------------------------------------------------------------
# SAGE: run once, output shared across all passes
# ---------------------------------------------------------------------------

echo ""
echo "========== SAGE: running (output shared across all passes) =========="

T_SAGE_START=$(now_sec)
${MY_ROOT}/bin/sage input/millennium_sage_hdf5.par
SAGE_STATUS=$?
T_SAGE_END=$(now_sec)
T_SAGE=$(elapsed_sec $T_SAGE_START $T_SAGE_END)

if [ $SAGE_STATUS -ne 0 ]; then
    echo "ERROR: SAGE failed with exit code $SAGE_STATUS"
    exit 1
fi

mv output/millennium ${OUTPUTDIR}/sage_output
echo "✓ SAGE complete (${T_SAGE}s)"

# ---------------------------------------------------------------------------
# Pass 1: --centralgalaxies passed explicitly to sage2kdtree and cli_lightcone
# ---------------------------------------------------------------------------

echo ""
echo "######################################################################"
echo "# PASS 1: explicit --centralgalaxies flag"
echo "######################################################################"

echo "========== PASS 1: sage2kdtree --centralgalaxies =========="
T_S2K1_START=$(now_sec)
${MY_ROOT}/bin/sage2kdtree \
    --centralgalaxies \
    -s ${OUTPUTDIR}/sage_output \
    -p input/millennium_sage_hdf5.par \
    -a input/millennium/trees/millennium.a_list \
    -o ${OUTPUTDIR}/pass1-${RAWNAME}-kdtree.h5 \
    --ppc 1000 -v 2
S2K1_STATUS=$?
T_S2K1_END=$(now_sec)
T_S2K1=$(elapsed_sec $T_S2K1_START $T_S2K1_END)

if [ $S2K1_STATUS -ne 0 ] || [ ! -f "${OUTPUTDIR}/pass1-${RAWNAME}-kdtree.h5" ]; then
    echo "ERROR: sage2kdtree failed (pass 1) - exit code ${S2K1_STATUS}, output file exists: $([ -f "${OUTPUTDIR}/pass1-${RAWNAME}-kdtree.h5" ] && echo yes || echo no)"
    T_S2K1="${T_S2K1:-N/A}"; T_CLI1="SKIPPED"; GALAXIES_PASS1="N/A"
else
    echo "========== PASS 1: cli_lightcone --centralgalaxies =========="
    T_CLI1_START=$(now_sec)
    ${MY_ROOT}/bin/cli_lightcone \
        --centralgalaxies \
        --dataset ${OUTPUTDIR}/pass1-${RAWNAME}-kdtree.h5 \
        --decmin ${LC_DECMIN} --decmax ${LC_DECMAX} \
        --ramin  ${LC_RAMIN}  --ramax  ${LC_RAMAX} \
        --zmin   ${LC_ZMIN}   --zmax   ${LC_ZMAX} \
        --outdir ${OUTPUTDIR} \
        --outfile pass1-${RAWNAME}-lightcone.h5
    CLI1_STATUS=$?
    T_CLI1_END=$(now_sec)
    T_CLI1=$(elapsed_sec $T_CLI1_START $T_CLI1_END)

    if [ $CLI1_STATUS -ne 0 ]; then
        echo "ERROR: cli_lightcone failed (pass 1)"
        GALAXIES_PASS1_TOTAL="N/A"; GALAXIES_PASS1_INSIDE="N/A"; GALAXIES_PASS1_OUTSIDE="N/A"
    else
        read GALAXIES_PASS1_TOTAL GALAXIES_PASS1_INSIDE GALAXIES_PASS1_OUTSIDE \
            <<< $(count_galaxies_stats ${OUTPUTDIR}/pass1-${RAWNAME}-lightcone.h5 \
                  $LC_DECMIN $LC_DECMAX $LC_RAMIN $LC_RAMAX $LC_ZMIN $LC_ZMAX)
        cp ${OUTPUTDIR}/pass1-${RAWNAME}-lightcone.h5 ${VALIDATIONDIR}/
        python3 ${MY_ROOT}/src/plot_lightcone.py ${VALIDATIONDIR}/pass1-${RAWNAME}-lightcone.h5 SnapNum
        cp lightcone_3d_SnapNum.png ${VALIDATIONDIR}/pass1-${RAWNAME}-lightcone-snapnum.png 2>/dev/null || true
        python3 ${MY_ROOT}/src/plot_lightcone.py ${VALIDATIONDIR}/pass1-${RAWNAME}-lightcone.h5 redshift_cosmological
        cp lightcone_3d_redshift_cosmological.png ${VALIDATIONDIR}/pass1-${RAWNAME}-lightcone-redshift.png 2>/dev/null || true
    fi
fi

# ---------------------------------------------------------------------------
# Pass 2: DEFAULT_MODE="--centralgalaxies", flag passed via $DEFAULT_MODE
# ---------------------------------------------------------------------------

echo ""
echo "######################################################################"
echo "# PASS 2: DEFAULT_MODE=\"--centralgalaxies\" (flag via variable)"
echo "######################################################################"
export DEFAULT_MODE=--centralgalaxies

echo "========== PASS 2: sage2kdtree \$DEFAULT_MODE =========="
T_S2K2_START=$(now_sec)
${MY_ROOT}/bin/sage2kdtree \
    -s ${OUTPUTDIR}/sage_output \
    -p input/millennium_sage_hdf5.par \
    -a input/millennium/trees/millennium.a_list \
    -o ${OUTPUTDIR}/pass2-${RAWNAME}-kdtree.h5 \
    --ppc 1000 -v 2 ${DEFAULT_MODE}
S2K2_STATUS=$?
T_S2K2_END=$(now_sec)
T_S2K2=$(elapsed_sec $T_S2K2_START $T_S2K2_END)

if [ $S2K2_STATUS -ne 0 ] || [ ! -f "${OUTPUTDIR}/pass2-${RAWNAME}-kdtree.h5" ]; then
    echo "ERROR: sage2kdtree failed (pass 2) - exit code ${S2K2_STATUS}, output file exists: $([ -f "${OUTPUTDIR}/pass2-${RAWNAME}-kdtree.h5" ] && echo yes || echo no)"
    T_CLI2="SKIPPED"; GALAXIES_PASS2="N/A"
else
    # cli_lightcone does NOT get --centralgalaxies in pass 2: only sage2kdtree does.
    # This mirrors the else-branch in run_test_sage_hdf5.sh where $DEFAULT_MODE is
    # only used in the sage2kdtree call, not in cli_lightcone.
    echo "========== PASS 2: cli_lightcone (no flag) =========="
    T_CLI2_START=$(now_sec)
    ${MY_ROOT}/bin/cli_lightcone \
        --dataset ${OUTPUTDIR}/pass2-${RAWNAME}-kdtree.h5 \
        --decmin ${LC_DECMIN} --decmax ${LC_DECMAX} \
        --ramin  ${LC_RAMIN}  --ramax  ${LC_RAMAX} \
        --zmin   ${LC_ZMIN}   --zmax   ${LC_ZMAX} \
        --outdir ${OUTPUTDIR} \
        --outfile pass2-${RAWNAME}-lightcone.h5
    CLI2_STATUS=$?
    T_CLI2_END=$(now_sec)
    T_CLI2=$(elapsed_sec $T_CLI2_START $T_CLI2_END)

    if [ $CLI2_STATUS -ne 0 ]; then
        echo "ERROR: cli_lightcone failed (pass 2)"
        GALAXIES_PASS2_TOTAL="N/A"; GALAXIES_PASS2_INSIDE="N/A"; GALAXIES_PASS2_OUTSIDE="N/A"
    else
        read GALAXIES_PASS2_TOTAL GALAXIES_PASS2_INSIDE GALAXIES_PASS2_OUTSIDE \
            <<< $(count_galaxies_stats ${OUTPUTDIR}/pass2-${RAWNAME}-lightcone.h5 \
                  $LC_DECMIN $LC_DECMAX $LC_RAMIN $LC_RAMAX $LC_ZMIN $LC_ZMAX)
        cp ${OUTPUTDIR}/pass2-${RAWNAME}-lightcone.h5 ${VALIDATIONDIR}/
        python3 ${MY_ROOT}/src/plot_lightcone.py ${VALIDATIONDIR}/pass2-${RAWNAME}-lightcone.h5 SnapNum
        cp lightcone_3d_SnapNum.png ${VALIDATIONDIR}/pass2-${RAWNAME}-lightcone-snapnum.png 2>/dev/null || true
        python3 ${MY_ROOT}/src/plot_lightcone.py ${VALIDATIONDIR}/pass2-${RAWNAME}-lightcone.h5 redshift_cosmological
        cp lightcone_3d_redshift_cosmological.png ${VALIDATIONDIR}/pass2-${RAWNAME}-lightcone-redshift.png 2>/dev/null || true
    fi
fi

# ---------------------------------------------------------------------------
# Pass 3: DEFAULT_MODE="", $DEFAULT_MODE expands to nothing
# ---------------------------------------------------------------------------

echo ""
echo "######################################################################"
echo "# PASS 3: DEFAULT_MODE=\"\" (blank - no centralgalaxies flag)"
echo "######################################################################"
export DEFAULT_MODE=

echo "========== PASS 3: sage2kdtree \$DEFAULT_MODE =========="
T_S2K3_START=$(now_sec)
${MY_ROOT}/bin/sage2kdtree \
    -s ${OUTPUTDIR}/sage_output \
    -p input/millennium_sage_hdf5.par \
    -a input/millennium/trees/millennium.a_list \
    -o ${OUTPUTDIR}/pass3-${RAWNAME}-kdtree.h5 \
    --ppc 1000 -v 2 ${DEFAULT_MODE}
S2K3_STATUS=$?
T_S2K3_END=$(now_sec)
T_S2K3=$(elapsed_sec $T_S2K3_START $T_S2K3_END)

if [ $S2K3_STATUS -ne 0 ] || [ ! -f "${OUTPUTDIR}/pass3-${RAWNAME}-kdtree.h5" ]; then
    echo "ERROR: sage2kdtree failed (pass 3) - exit code ${S2K3_STATUS}, output file exists: $([ -f "${OUTPUTDIR}/pass3-${RAWNAME}-kdtree.h5" ] && echo yes || echo no)"
    T_CLI3="SKIPPED"; GALAXIES_PASS3="N/A"
else
    echo "========== PASS 3: cli_lightcone (no flag) =========="
    T_CLI3_START=$(now_sec)
    ${MY_ROOT}/bin/cli_lightcone \
        --dataset ${OUTPUTDIR}/pass3-${RAWNAME}-kdtree.h5 \
        --decmin ${LC_DECMIN} --decmax ${LC_DECMAX} \
        --ramin  ${LC_RAMIN}  --ramax  ${LC_RAMAX} \
        --zmin   ${LC_ZMIN}   --zmax   ${LC_ZMAX} \
        --outdir ${OUTPUTDIR} \
        --outfile pass3-${RAWNAME}-lightcone.h5 ${DEFAULT_MODE}
    CLI3_STATUS=$?
    T_CLI3_END=$(now_sec)
    T_CLI3=$(elapsed_sec $T_CLI3_START $T_CLI3_END)

    if [ $CLI3_STATUS -ne 0 ]; then
        echo "ERROR: cli_lightcone failed (pass 3)"
        GALAXIES_PASS3_TOTAL="N/A"; GALAXIES_PASS3_INSIDE="N/A"; GALAXIES_PASS3_OUTSIDE="N/A"
    else
        read GALAXIES_PASS3_TOTAL GALAXIES_PASS3_INSIDE GALAXIES_PASS3_OUTSIDE \
            <<< $(count_galaxies_stats ${OUTPUTDIR}/pass3-${RAWNAME}-lightcone.h5 \
                  $LC_DECMIN $LC_DECMAX $LC_RAMIN $LC_RAMAX $LC_ZMIN $LC_ZMAX)
        cp ${OUTPUTDIR}/pass3-${RAWNAME}-lightcone.h5 ${VALIDATIONDIR}/
        python3 ${MY_ROOT}/src/plot_lightcone.py ${VALIDATIONDIR}/pass3-${RAWNAME}-lightcone.h5 SnapNum
        cp lightcone_3d_SnapNum.png ${VALIDATIONDIR}/pass3-${RAWNAME}-lightcone-snapnum.png 2>/dev/null || true
        python3 ${MY_ROOT}/src/plot_lightcone.py ${VALIDATIONDIR}/pass3-${RAWNAME}-lightcone.h5 redshift_cosmological
        cp lightcone_3d_redshift_cosmological.png ${VALIDATIONDIR}/pass3-${RAWNAME}-lightcone-redshift.png 2>/dev/null || true
    fi
fi

# ---------------------------------------------------------------------------
# Clean up log artefacts
# ---------------------------------------------------------------------------

rm -f log.00000
rm -rf log

# ---------------------------------------------------------------------------
# Compute per-pass totals
# ---------------------------------------------------------------------------

total_pass() {
    local s2k=$1 cli=$2
    if [[ "$s2k" =~ ^[0-9.]+$ ]] && [[ "$cli" =~ ^[0-9.]+$ ]]; then
        python3 -c "print(f'{$s2k + $cli:.1f}')"
    else
        echo "N/A"
    fi
}

TOTAL1=$(total_pass "${T_S2K1:-N/A}" "${T_CLI1:-N/A}")
TOTAL2=$(total_pass "${T_S2K2:-N/A}" "${T_CLI2:-N/A}")
TOTAL3=$(total_pass "${T_S2K3:-N/A}" "${T_CLI3:-N/A}")

# ---------------------------------------------------------------------------
# Validation:
#   PASS if pass1_inside == pass2_inside AND pass1_outside > 0
# ---------------------------------------------------------------------------

VALIDATION_STATUS="N/A"

P1_INSIDE="${GALAXIES_PASS1_INSIDE:-N/A}"
P1_OUTSIDE="${GALAXIES_PASS1_OUTSIDE:-N/A}"
P2_INSIDE="${GALAXIES_PASS2_INSIDE:-N/A}"

if [[ "$P1_INSIDE"  =~ ^[0-9]+$ ]] && \
   [[ "$P1_OUTSIDE" =~ ^[0-9]+$ ]] && \
   [[ "$P2_INSIDE"  =~ ^[0-9]+$ ]]; then
    INSIDE_MATCH=0
    PASS1_HAS_OUTSIDE=0
    [ "$P1_INSIDE" -eq "$P2_INSIDE" ] && INSIDE_MATCH=1
    [ "$P1_OUTSIDE" -gt 0 ]           && PASS1_HAS_OUTSIDE=1

    if [ $INSIDE_MATCH -eq 1 ] && [ $PASS1_HAS_OUTSIDE -eq 1 ]; then
        VALIDATION_STATUS="PASS"
    else
        REASONS=""
        [ $INSIDE_MATCH -eq 0 ]      && REASONS="${REASONS} inside_pass1(${P1_INSIDE})!=inside_pass2(${P2_INSIDE})"
        [ $PASS1_HAS_OUTSIDE -eq 0 ] && REASONS="${REASONS} pass1_outside==0"
        VALIDATION_STATUS="FAIL [${REASONS# }]"
    fi
fi

ORDERING_STATUS="N/A"
P1_TOTAL="${GALAXIES_PASS1_TOTAL:-N/A}"
P2_TOTAL="${GALAXIES_PASS2_TOTAL:-N/A}"
P3_TOTAL="${GALAXIES_PASS3_TOTAL:-N/A}"

if [[ "$P1_TOTAL" =~ ^[0-9]+$ ]] && \
   [[ "$P2_TOTAL" =~ ^[0-9]+$ ]] && \
   [[ "$P3_TOTAL" =~ ^[0-9]+$ ]]; then
    if [ "$P1_TOTAL" -gt "$P2_TOTAL" ] && [ "$P2_TOTAL" -gt "$P3_TOTAL" ]; then
        ORDERING_STATUS="PASS"
    else
        ORDERING_STATUS="FAIL [expected pass1(${P1_TOTAL})>pass2(${P2_TOTAL})>pass3(${P3_TOTAL})]"
    fi
fi

# ---------------------------------------------------------------------------
# Final report
# ---------------------------------------------------------------------------

echo ""
echo "======================================================================"
echo "  CENTRALGALAXIES OPTION: VALIDATION AND TIMING REPORT"
echo "======================================================================"
echo ""
echo "  Lightcone query:  dec=[${LC_DECMIN},${LC_DECMAX}]  ra=[${LC_RAMIN},${LC_RAMAX}]  z=[${LC_ZMIN},${LC_ZMAX}]"
echo ""
echo "  SAGE (shared, run once):  ${T_SAGE}s"
echo ""
printf "  %-22s  %-18s  %-18s  %-18s\n" \
    "Phase" "Pass 1" "Pass 2" "Pass 3"
printf "  %-22s  %-18s  %-18s  %-18s\n" \
    "Configuration" "explicit flag" "s2k:DM=flag" "s2k:DM=blank"
printf "  %-22s  %-18s  %-18s  %-18s\n" \
    "----------------------" "------------------" "------------------" "------------------"
printf "  %-22s  %-16ss  %-16ss  %-16ss\n" \
    "sage2kdtree" \
    "${T_S2K1:-N/A}" "${T_S2K2:-N/A}" "${T_S2K3:-N/A}"
printf "  %-22s  %-16ss  %-16ss  %-16ss\n" \
    "cli_lightcone" \
    "${T_CLI1:-N/A}" "${T_CLI2:-N/A}" "${T_CLI3:-N/A}"
printf "  %-22s  %-18s  %-18s  %-18s\n" \
    "----------------------" "------------------" "------------------" "------------------"
printf "  %-22s  %-16ss  %-16ss  %-16ss\n" \
    "Total (excl. SAGE)" \
    "${TOTAL1:-N/A}" "${TOTAL2:-N/A}" "${TOTAL3:-N/A}"
echo ""
echo "  Galaxy counts:"
printf "  %-54s  %10s  %10s  %10s\n" "Pass" "Total" "Inside" "Outside"
printf "  %-54s  %10s  %10s  %10s\n" \
    "------------------------------------------------------" "----------" "----------" "----------"
printf "  %-54s  %10s  %10s  %10s\n" \
    "Pass 1  (sage2kdtree --centralgalaxies, cli --centralgalaxies)" \
    "${GALAXIES_PASS1_TOTAL:-N/A}" "${GALAXIES_PASS1_INSIDE:-N/A}" "${GALAXIES_PASS1_OUTSIDE:-N/A}"
printf "  %-54s  %10s  %10s  %10s\n" \
    "Pass 2  (sage2kdtree DEFAULT_MODE=--centralgalaxies, cli none)" \
    "${GALAXIES_PASS2_TOTAL:-N/A}" "${GALAXIES_PASS2_INSIDE:-N/A}" "${GALAXIES_PASS2_OUTSIDE:-N/A}"
printf "  %-54s  %10s  %10s  %10s\n" \
    "Pass 3  (sage2kdtree DEFAULT_MODE=blank, cli none)" \
    "${GALAXIES_PASS3_TOTAL:-N/A}" "${GALAXIES_PASS3_INSIDE:-N/A}" "${GALAXIES_PASS3_OUTSIDE:-N/A}"
echo ""
echo "  Validation:"
printf "    inside_pass1==inside_pass2 AND pass1_outside>0:  %s\n" "${VALIDATION_STATUS}"
printf "    total pass1>pass2>pass3:                         %s\n" "${ORDERING_STATUS}"
echo ""
echo "  Plots saved in: ${VALIDATIONDIR}/"
echo "  Output files in: ${OUTPUTDIR}/"
echo "======================================================================"
