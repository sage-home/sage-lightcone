#!/bin/bash

# What's my code's root directory
export MY_SCRIPT=${BASH_SOURCE:-$0}
export MY_SCRIPTS_DIRECTORY=$(dirname $MY_SCRIPT)
export MY_ROOT=$(cd ${MY_SCRIPTS_DIRECTORY} && pwd)

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
pip3 list

# Mac-specific setup script for M3 Mac with Homebrew
# This replaces the module-based setup.sh for macOS development

echo "Setting up environment for M3 Mac..."

# Ensure Homebrew is in PATH
export PATH="/opt/homebrew/bin:/opt/homebrew/sbin:$PATH"

# Set compiler preferences for Apple Silicon
export CC="clang"
export CXX="clang++"

# HDF5-MPI paths - robust discovery
HDF5_CANDIDATES=(
    "/opt/homebrew/Cellar/hdf5-mpi"         # Homebrew hdf5-mpi (any version)
    "/opt/homebrew"                         # Homebrew standard location
    "/usr/local"                            # Alternative location
)

HDF5_ROOT=""
for candidate in "${HDF5_CANDIDATES[@]}"; do
    if [[ "$candidate" == "/opt/homebrew/Cellar/hdf5-mpi" ]]; then
        # Find the latest version in the Cellar
        if [[ -d "$candidate" ]]; then
            latest_version=$(ls "$candidate" | sort -V | tail -1)
            if [[ -n "$latest_version" && -d "$candidate/$latest_version" ]]; then
                HDF5_ROOT="$candidate/$latest_version"
                echo "Found HDF5-MPI at: $HDF5_ROOT"
                break
            fi
        fi
    elif [[ -f "$candidate/bin/h5cc" ]]; then
        HDF5_ROOT="$candidate"
        echo "Found HDF5-MPI at: $HDF5_ROOT"
        break
    fi
done

if [[ -n "$HDF5_ROOT" ]]; then
    export HDF5_ROOT
    export HDF5_DIR="$HDF5_ROOT"
    export PATH="$HDF5_ROOT/bin:$PATH"
else
    echo "Warning: HDF5-MPI not found. Using system PATH."
    if ! command -v h5cc &> /dev/null; then
        echo "Error: h5cc not found. Please install HDF5 with MPI support:"
        echo "  brew install hdf5-mpi"
        return 1
    fi
fi

# OpenMPI paths  
export MPI_ROOT="/opt/homebrew"
export PATH="$MPI_ROOT/bin:$PATH"

# Boost paths
export BOOST_ROOT="/opt/homebrew"

# GSL paths
export GSL_ROOT="/opt/homebrew"

# CMake paths - robust discovery
CMAKE_CANDIDATES=(
    "/opt/homebrew"                          # Homebrew on Apple Silicon
    "/usr/local"                            # Homebrew on Intel
    "/Applications/CMake.app/Contents"      # CMake.app installation
    "/usr"                                  # System installation
)

CMAKE_ROOT=""
for candidate in "${CMAKE_CANDIDATES[@]}"; do
    if [[ -x "$candidate/bin/cmake" ]]; then
        CMAKE_ROOT="$candidate"
        echo "Found CMake at: $CMAKE_ROOT"
        break
    fi
done

if [[ -n "$CMAKE_ROOT" ]]; then
    export CMAKE_ROOT
    export PATH="$CMAKE_ROOT/bin:$PATH"
else
    echo "Warning: CMake not found in expected locations. Using system PATH."
    # Check if cmake is available in PATH
    if ! command -v cmake &> /dev/null; then
        echo "Error: cmake not found in PATH either. Please install CMake."
        echo "You can install it via:"
        echo "  - Homebrew: brew install cmake"
        echo "  - Download from: https://cmake.org/download/"
        return 1
    fi
fi

# Set PKG_CONFIG_PATH to help CMake find libraries
export PKG_CONFIG_PATH="/opt/homebrew/lib/pkgconfig:$PKG_CONFIG_PATH"

# Set library and include paths
export LIBRARY_PATH="/opt/homebrew/lib:$LIBRARY_PATH"
export CPATH="/opt/homebrew/include:$CPATH"

echo "Environment configured for:"
echo "  - HDF5-MPI: $(h5cc -showconfig 2>/dev/null | grep 'HDF5 Version' || echo 'Not found')"
echo "  - OpenMPI: $(mpicc --version 2>/dev/null | head -1 || echo 'Not found')"
echo "  - Boost: $(brew list boost --versions 2>/dev/null || echo 'Not found')"
echo "  - GSL: $(brew list gsl --versions 2>/dev/null || echo 'Not found')"
echo "  - CMake: $(cmake --version 2>/dev/null | head -1 || echo 'Not found')"

# Validate critical dependencies
validation_failed=false

if ! command -v cmake &> /dev/null; then
    echo "❌ CMake validation failed"
    validation_failed=true
else
    echo "✅ CMake: $(cmake --version | head -1)"
fi

if ! command -v h5cc &> /dev/null; then
    echo "❌ HDF5 validation failed"
    validation_failed=true
else
    # Check if HDF5 has parallel support
    if h5cc -showconfig 2>/dev/null | grep -i "Parallel HDF5.*ON"; then
        echo "✅ HDF5-MPI: Parallel support enabled"
    else
        echo "⚠️  HDF5: Found but parallel support may not be enabled"
    fi
fi

if ! command -v mpicc &> /dev/null; then
    echo "❌ MPI validation failed"
    validation_failed=true
else
    echo "✅ MPI: $(mpicc --version 2>/dev/null | head -1)"
fi

if [[ "$validation_failed" == "true" ]]; then
    echo ""
    echo "❌ Some dependencies are missing. Build may fail."
    echo "Consider installing missing dependencies with Homebrew."
    return 1
else
    echo ""
    echo "✅ All critical dependencies found. Ready to build!"
fi
