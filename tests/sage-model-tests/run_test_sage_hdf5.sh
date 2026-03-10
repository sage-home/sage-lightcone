#!/bin/bash

# Parse command line arguments
BUILD_CENTRAL_GALAXY_INDEX=0
FORCE_REBUILD=0
for arg in "$@"; do
    case $arg in
        --rebuild)
            FORCE_REBUILD=1
            echo "Force rebuild requested"
            ;;
        --centralgalaxies)
            BUILD_CENTRAL_GALAXY_INDEX=1
            echo "Central galaxies index build requested"
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
export MY_SCRIPT=${BASH_SOURCE:-$0}
export MY_SCRIPTS_DIRECTORY=$(cd "$(dirname $MY_SCRIPT)" && pwd)
export MY_ROOT=$(cd ${MY_SCRIPTS_DIRECTORY}/../.. && pwd)
export MY_SCRIPT=

# Change to script directory - required for relative paths like ./first_run.sh
cd "${MY_SCRIPTS_DIRECTORY}"

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
export OUTPUTDIR=output_test_sage_hdf5
export VALIDATIONDIR=validate_test_sage_hdf5
# keep the validation directory separate from the output directory to avoid accidentally deleting it during cleanup
mkdir -p ${VALIDATIONDIR}

# Clean up previous test outputs, but preserve downloaded tree files
rm -rf output
rm -rf ${OUTPUTDIR}
rm -f *${RAWNAME}*
# Only remove input/*.par files, preserve input/millennium/trees/
rm -f input/*.par 2>/dev/null

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

# Check if tree files already exist - skip download if so
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

mkdir -p output/millennium/

${MY_ROOT}/bin/sage input/millennium_sage_hdf5.par
mv output ${OUTPUTDIR}

echo "Create kdtree in one step, then generate lightcone and plot SnapNum..."
source ${MY_ROOT}/.venv/bin/activate

if [ $BUILD_CENTRAL_GALAXY_INDEX -eq 1 ]; then
    echo ******* Building with central galaxy index included *******
    ${MY_ROOT}/bin/sage2kdtree -s ${OUTPUTDIR}/millennium -p input/millennium_sage_hdf5.par -a input/millennium/trees/millennium.a_list -o ${OUTPUTDIR}/${RAWNAME}-kdtree-onestep.h5 --ppc 1000 -v 2 --centralgalaxies
    echo "Test by creating a lightcone and plotting SnapNum..."
    ${MY_ROOT}/bin/cli_lightcone --dataset ${OUTPUTDIR}/${RAWNAME}-kdtree-onestep.h5 --decmin 0 --decmax 1 --ramin 0 --ramax 1 --zmin 0 --zmax 1 --outdir ${OUTPUTDIR} --outfile $RAWNAME-lightcone.h5 --centralgalaxies
    cp ${OUTPUTDIR}/${RAWNAME}-lightcone.h5 ${VALIDATIONDIR}/${RAWNAME}-lightcone-centralgalaxies.h5
    python3 ${MY_ROOT}/src/plot_lightcone.py ${VALIDATIONDIR}/${RAWNAME}-lightcone-centralgalaxies.h5 SnapNum
    cp lightcone_3d_SnapNum.png ${VALIDATIONDIR}/$RAWNAME-lightcone-snapnum-centralgalaxies.png
    python3 ${MY_ROOT}/src/plot_lightcone.py ${VALIDATIONDIR}/${RAWNAME}-lightcone-centralgalaxies.h5 redshift_cosmological
    cp lightcone_3d_redshift_cosmological.png ${VALIDATIONDIR}/$RAWNAME-lightcone-redshift-centralgalaxies.png
else
#   Build with central galaxy index
#   export DEFAULT_MODE=--centralgalaxies
#    ${MY_ROOT}/bin/sage2kdtree -s ${OUTPUTDIR}/millennium -p input/millennium_sage_hdf5.par -a input/millennium/trees/millennium.a_list -o ${OUTPUTDIR}/${RAWNAME}-kdtree-onestep.h5 --ppc 1000 -v 2
#   Build without central galaxy index (should be faster to build and query, but larger file size)
    export DEFAULT_MODE=
    ${MY_ROOT}/bin/sage2kdtree -s ${OUTPUTDIR}/millennium -p input/millennium_sage_hdf5.par -a input/millennium/trees/millennium.a_list -o ${OUTPUTDIR}/${RAWNAME}-kdtree-onestep.h5 --ppc 1000 -v 2 ${DEFAULT_MODE}
    echo "Test by creating a lightcone and plotting SnapNum..."
    ${MY_ROOT}/bin/cli_lightcone --dataset ${OUTPUTDIR}/${RAWNAME}-kdtree-onestep.h5 --decmin 0 --decmax 1 --ramin 0 --ramax 1 --zmin 0 --zmax 1 --outdir ${OUTPUTDIR} --outfile $RAWNAME-lightcone.h5
    cp ${OUTPUTDIR}/${RAWNAME}-lightcone.h5 ${VALIDATIONDIR}/${RAWNAME}-lightcone.h5
    python3 ${MY_ROOT}/src/plot_lightcone.py ${VALIDATIONDIR}/${RAWNAME}-lightcone.h5 SnapNum
    cp lightcone_3d_SnapNum.png ${VALIDATIONDIR}/$RAWNAME-lightcone-snapnum.png
    python3 ${MY_ROOT}/src/plot_lightcone.py ${VALIDATIONDIR}/${RAWNAME}-lightcone.h5 redshift_cosmological
    cp lightcone_3d_redshift_cosmological.png ${VALIDATIONDIR}/$RAWNAME-lightcone-redshift.png
fi

# Clean up
rm -f log.00000
rm -rf log
