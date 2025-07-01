#!/bin/bash

# Clean up any old CMake cache files
find . -name CMakeCache.txt -exec rm -f {} \;

# Build PugiXML first
echo "Building PugiXML..."
cd dep/pugixml-1.7/scripts
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=.. ..
make -j20
make install
cd ../../../..

echo "Building main project..."
# CMAKE_BUILD_TYPE = Debug or Release
cmake -B="./bin" -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=./runtime \
-DPUGIXML_HOME=./dep/pugixml-1.7 \
.
cd bin
make clean
make -j20
