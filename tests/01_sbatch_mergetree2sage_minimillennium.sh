#!/bin/bash
export MY_ROOT=$(realpath "$(dirname "$0")/../..")
export SIMULATION=minimillennium
export task=mergetree2sage
./download_minimillennium.sh
cat > ./${SIMULATION}_${task}.sh  <<EOF
#!/bin/bash
#SBATCH --job-name=${SIMULATION}-${task}
#SBATCH --time=8:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=../LOGS/${SIMULATION}_${task}.o
#SBATCH --error=../LOGS/${SIMULATION}_${task}.e
export MY_ROOT=${MY_ROOT}
export SIMULATION=${SIMULATION}
./mergetree2sage.sh
EOF
chmod +x ${SIMULATION}_${task}.sh
if [ "${1}" = "interactive" ]; then
    ./${SIMULATION}_${task}.sh
else
    sbatch ${SIMULATION}_${task}.sh
fi
