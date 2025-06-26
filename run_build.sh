#!/bin/bash
find . -name CMakeCache.txt -exec rm -f {} \;
# CMAKE_BUILD_TYPE = Debug or Release
cmake -B="./bin" -DCMAKE_BUILD_TYPE=Debug \
-DPUGIXML_HOME=./dep/pugixml-1.7 \
.
cd bin
make clean
make -j20
