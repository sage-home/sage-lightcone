#!/bin/bash

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
export MY_SCRIPT=${BASH_SOURCE:-$0}
export MY_SCRIPTS_DIRECTORY=$(dirname $MY_SCRIPT)
export MY_ROOT=$(cd ${MY_SCRIPTS_DIRECTORY}/../.. && pwd)
export MY_SCRIPT=
export MY_SCRIPTS_DIRECTORY=

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
export OUTPUTDIR=output_sage_hdf5_one_step

# Clean up previous test outputs, but preserve downloaded tree files
rm -rf output
rm -rf ${OUTPUTDIR}
rm -f *${RAWNAME}*
# Only remove input/*.par files, preserve input/millennium/trees/
rm -f input/*.par 2>/dev/null

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
cp mypar_files/millennium.par input
cp mypar_files/millennium_minus1.par input
cp mypar_files/millennium_sage_hdf5.par input

${MY_ROOT}/bin/sage input/millennium_sage_hdf5.par
mv output ${OUTPUTDIR}

echo "Create kdtree in one step"

${MY_ROOT}/bin/sage2kdtree -s ${OUTPUTDIR}/millennium -p input/millennium_sage_hdf5.par -a input/millennium/trees/millennium.a_list -o ${OUTPUTDIR}/${RAWNAME}-kdtree-onestep.h5 --ppc 1000 -v 2

# echo "Creating kdtree indexed HDF5..."
# ${MY_ROOT}/bin/sageh5toh5 -m convert -v 2 -s ${OUTPUTDIR}/millennium -p input/millennium_sage_hdf5.par -a input/millennium/trees/millennium.a_list -o ${OUTPUTDIR}/${RAWNAME}-depthfirstordered.h5
# ${MY_ROOT}/bin/sageimport --settings ${OUTPUTDIR}/${RAWNAME}-depthfirstordered_import_settings.xml
# ${MY_ROOT}/bin/dstreeinit --mode kdtree --tree ${OUTPUTDIR}/out_$RAWNAME-depthfirstordered.h5 --sage ${OUTPUTDIR}/$RAWNAME-bysnap.h5 --output ${OUTPUTDIR}/$RAWNAME-kdtree.h5

echo "Test by creating a lightcone and plotting snapnum..."
#${MY_ROOT}/bin/cli_lightcone --dataset ${OUTPUTDIR}/${RAWNAME}-kdtree-onestep.h5 --decmin 0 --decmax 10 --ramin 0 --ramax 30 --zmin 0 --zmax 0.5 --outdir ${OUTPUTDIR} --outfile $RAWNAME-lightcone.h5
${MY_ROOT}/bin/cli_lightcone --dataset ${OUTPUTDIR}/${RAWNAME}-kdtree-onestep.h5 --decmin 0 --decmax 0 --ramin 0 --ramax 1 --zmin 0 --zmax 1 --outdir ${OUTPUTDIR} --outfile $RAWNAME-lightcone.h5
source ${MY_ROOT}/.venv/bin/activate
python3 ${MY_ROOT}/src/plot_lightcone.py ${OUTPUTDIR}/$RAWNAME-lightcone.h5 snapnum

# Clean up
rm -f log.00000
rm -rf log
