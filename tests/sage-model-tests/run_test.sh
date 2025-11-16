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

rm -rf input
rm -rf output
rm -f *mymillennium*
./first_run.sh
echo ==== Have finished with first_run.sh ====
mkdir -p output/millennium/
cp mypar_files/millennium.par input
cp mypar_files/millennium_minus1.par input
${MY_ROOT}/bin/sage input/millennium.par
${MY_ROOT}/bin/sage2h5 -m convert -v 2 -s output/millennium -p input/millennium_minus1.par -a input/millennium/trees/millennium.a_list -o mymillennium.h5

export RAWNAME=mymillennium

source ${MY_ROOT}/.venv/bin/activate
cat > MandatoryList.txt << !
posx
posy
posz
velx
vely
velz
snapnum
type
cold_gas
metals_cold_gas
sfr_disk
sfr_bulge
descendant
stellar_mass
global_index
!
rm -rf log
mkdir -p log
rm -f log.00000
python3 ${MY_ROOT}/src/sageimport_mpi_HDF2HDF/main.py ${RAWNAME}_import_settings.xml False
cp log/logfile0.log log_1_of_3.txt
${MY_ROOT}/bin/dstreeinit --mode tree2sage --tree out_$RAWNAME.h5 --output $RAWNAME-bysnap.h5
cp log.00000 log_2_of_3.txt
${MY_ROOT}/bin/dstreeinit --mode init --tree out_$RAWNAME.h5 --sage $RAWNAME-bysnap.h5 --output $RAWNAME-kdtree.h5
cp log.00000 log_3_of_3.txt
cp out_$RAWNAME.xml $RAWNAME-kdtree.xml
${MY_ROOT}/bin/cli_lightcone --dataset $RAWNAME-kdtree.h5 --decmin 0 --decmax 1 --ramin 0 --ramax 1 --zmin 0 --zmax 1 --outfile $RAWNAME-lightcone.h5
python3 ${MY_ROOT}/src/plot_lightcone.py output/$RAWNAME-lightcone.h5 snapnum
