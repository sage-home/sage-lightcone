#!/bin/bash

# Clean up any old CMake cache files
find . -name CMakeCache.txt -exec rm -f {} \;

echo "Building main project..."
# CMAKE_BUILD_TYPE = Debug or Release
cmake -B="./bin" -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=./runtime \
.
cd bin
make clean
make -j20

echo "Build sage-model for testing..."
cd ..
./build_sage.sh
