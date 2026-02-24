#!/bin/bash
# Check all executables in bin/ for MPI dependencies using otool

FOUND_MPI=0

echo "Checking bin/ executables for MPI dependencies..."
echo ""

for exe in bin/*; do
  if [ -f "$exe" ] && [ -x "$exe" ]; then
    name=$(basename "$exe")
    if otool -L "$exe" 2>/dev/null | grep -qi mpi; then
      echo "FAIL: $name has MPI dependencies:"
      otool -L "$exe" | grep -i mpi
      FOUND_MPI=1
    else
      echo "OK:   $name - no MPI dependencies"
    fi
  fi
done

echo ""
if [ $FOUND_MPI -eq 0 ]; then
  echo "SUCCESS: All executables are MPI-free"
  exit 0
else
  echo "FAILURE: Some executables have MPI dependencies"
  exit 1
fi
