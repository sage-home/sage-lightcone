#!/bin/bash
# What's my code's root directory
export MY_SCRIPT=${BASH_SOURCE:-$0}
export MY_SCRIPTS_DIRECTORY=$(dirname $MY_SCRIPT)
export MY_ROOT=$(cd ${MY_SCRIPTS_DIRECTORY} && pwd)
export MY_SCRIPT=

if [ "$SYS_ARCH" = "milan" ];then
   module purge &> /dev/null
   module load gcc/12.2.0
   module load boost/1.81.0
   module load openmpi/4.1.4
   module load hdf5/1.14.0
   module load gcccore/12.2.0
   module load gsl/2.7
   module load cmake/3.24.3
fi
# Check for Python virtual environment and install if needed
if [[ ! -d "$MY_ROOT/.venv" ]]; then
    echo "Virtual environment not found. Running install.sh..."
    if [[ -f "$MY_ROOT/install.sh" ]]; then
        cd "$MY_ROOT"
        chmod +x install.sh
        ./install.sh
        if [[ $? -eq 0 ]]; then
            echo "✅ Virtual environment created successfully"
        else
            echo "❌ Failed to create virtual environment"
            return 1
        fi
    else
        echo "❌ install.sh not found. Cannot create virtual environment."
        return 1
    fi
else
    echo "✅ Virtual environment found at $MY_ROOT/.venv"
fi

# Activate virtual environment
source $MY_ROOT/.venv/bin/activate
echo "Python virtual environment activated"
