#!/bin/bash

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

export RAWNAME=myhdf5millennium
export OUTPUTDIR=output_sage_hdf5

rm -rf input
rm -rf output
rm -rf ${OUTPUTDIR}
rm -f *${RAWNAME}*
./first_run.sh
echo ==== Have finished with first_run.sh ====
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
${MY_ROOT}/bin/cli_lightcone --dataset ${OUTPUTDIR}/${RAWNAME}-kdtree-onestep.h5 --decmin 0 --decmax 1 --ramin 0 --ramax 1 --zmin 0 --zmax 1 --outdir ${OUTPUTDIR} --outfile $RAWNAME-lightcone.h5
source ${MY_ROOT}/.venv/bin/activate
python3 ${MY_ROOT}/src/plot_lightcone.py ${OUTPUTDIR}/$RAWNAME-lightcone.h5 snapnum

# Clean up
rm -f log.00000
rm -rf log
