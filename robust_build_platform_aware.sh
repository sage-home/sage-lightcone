#!/bin/bash

# Platform-aware build script that automatically selects the right setup
# This preserves Linux HPC compatibility while enabling macOS builds

# Ensure we are running from the script directory
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd "$DIR" || exit 1

rm -rf sage-model
rm -rf bin

./build_platform_aware.sh
