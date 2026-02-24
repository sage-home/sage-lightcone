# macOS Serial Build (USE_MPI=no) Status

**Date:** 24 February 2026
**Status:** WORKING (macOS), CI/CD FIXED

## Objective
Enable a "serial" build of the project on macOS (`USE_MPI=no`) that completely avoids linking to MPI libraries or including system MPI headers (`<mpi.h>`). This is useful for development and debugging on machines without complex MPI setups.

## The Problem
The project was originally designed with heavy MPI integration. Even when `USE_MPI=no` was requested, the build system was:
1.  **Header Conflict:** HDF5 headers (`H5public.h`) from the parallel installation include `<mpi.h>`. This conflicted with our internal `src/libhpc/mpi/mpi_stub.hh`, which defines mock MPI types for serial builds.
2.  **Include Path Order:** Boost's CMake configuration added `/opt/homebrew/include` (which contains symlinks to parallel HDF5) BEFORE the serial HDF5 include path.
3.  **Missing MPI Stubs:** `libhpc` relies on MPI functions that were not implemented in `mpi_stub.hh`.

## Fixes Applied

### 1. HDF5 Include Path Priority (`.cmake/3rd_party.cmake`)
*   **Problem:** CMake's `include_directories()` was adding HDF5 paths after Boost paths, causing `/opt/homebrew/include/hdf5.h` (parallel symlink) to be found first.
*   **Fix:** Changed `include_directories(${HDF5_INCLUDE_DIRS})` to `include_directories(BEFORE ${HDF5_INCLUDE_DIRS})` to ensure serial HDF5 headers take precedence.

### 2. Environment Isolation (`setup_mac.sh`)
*   **Fix:** Script now skips HDF5-MPI discovery when `USE_MPI=no` and avoids adding generic `/opt/homebrew/include` to paths.

### 3. CMake Configuration (`build_platform_aware.sh`)
*   **Fix:** Passes `-DUSE_MPI=OFF` to CMake and forces `HDF5_ROOT` to the serial HDF5 path.

### 4. MPI Stubs (`src/libhpc/mpi/mpi_stub.hh`)
*   **Fix:** Expanded the stub file with all MPI symbols required by `libhpc`:
    *   **Constants:** `MPI_MAX_PROCESSOR_NAME`, `MPI_STATUSES_IGNORE`, `MPI_STATUS_IGNORE`
    *   **Communicator:** `MPI_Comm_split`
    *   **Async Operations:** `MPI_Isend`, `MPI_Issend`, `MPI_Irecv`
    *   **Request Handling:** `MPI_Waitall`, `MPI_Testall`
    *   **Probing:** `MPI_Iprobe`
    *   **Collectives:** `MPI_Alltoallw`
    *   **Utility:** `MPI_Wtime`, `MPI_Get_processor_name`

## How to Build Serial Version

```bash
USE_MPI=no ./build_platform_aware.sh
```

## Verification

```bash
# Check executables exist
ls -la bin/cli_lightcone bin/sage2kdtree bin/sage2kdtree_four

# Test help output
./bin/cli_lightcone --help
```

## Current Status
*   **Environment:** Correctly detecting Serial HDF5 and ignoring MPI
*   **Compilation:** SUCCESSFUL
*   **Executables Built:**
    - `cli_lightcone`
    - `sage2kdtree`
    - `sage2kdtree_four`
    - `sage`
*   **Warnings:** Minor warning about archiver/compiler toolchain mismatch (can be ignored)

## CI/CD Fix (Feb 2026)

### Problem
`test_verify_kdtree` CI job was failing with:
```
H5public.h:65:10: fatal error: mpi.h: No such file or directory
```
The CI caches a parallel HDF5 build (`--enable-parallel`, `h5pcc` only, no `h5cc`).
`build_platform_aware.sh` with `USE_MPI=no` was supposed to detect this and use `mpicc`
as the compiler (so `mpi.h` is found via mpicc's include path), but `USING_PARALLEL_HDF5`
was evaluating to 0.

### Root Causes
1. **Detection logic bug (`build_platform_aware.sh`)**: The outer `if/elif` structure
   meant that when `HDF5_ROOT` was set but the file check `[ -f "$HDF5_ROOT/bin/h5pcc" ]`
   failed (for whatever reason - variable expansion issues, etc.), the fallback
   `command -v h5pcc` was never reached (it was in an `elif`, only executed when
   `HDF5_ROOT` was NOT set).

2. **Cache invalidation bug (`.gitlab-ci.yml`)**: The HDF5 cache check used
   `[ ! -f $HDF5_PREFIX/bin/h5cc ]`, but parallel HDF5 creates `h5pcc`, not `h5cc`.
   This caused HDF5 to be rebuilt on every CI run.

### Fixes Applied
1. **`build_platform_aware.sh`**: Simplified `USING_PARALLEL_HDF5` to a flat OR
   across all indicators — file checks under `$HDF5_ROOT`, `h5pcc` in PATH, or
   explicit CMake flags. Now any single indicator is sufficient.

2. **`.gitlab-ci.yml`**: Changed cache check from `h5cc` to `h5pcc` so the
   parallel HDF5 cache is correctly reused rather than rebuilt each time.

## Second CI/CD Fix (Feb 2026)

### Problem
After the first fix was committed, the `test_cli_validation` CI job still failed with
`cli_lightcone not found in artifacts`. The upstream `build` job was silently using
`gcc` instead of `mpicc` despite `cmake` being invoked with `-DCMAKE_C_COMPILER=mpicc`.

### Root Cause
**CMake `set()` before `project()` overrides `-D` cache entries** (`CMakeLists.txt`):
```cmake
# BROKEN: plain set() before project() shadows -DCMAKE_C_COMPILER=mpicc
set(CMAKE_C_COMPILER   "gcc")
set(CMAKE_CXX_COMPILER "g++")
```
Per CMake docs, a **normal variable** (no `CACHE` keyword) set before `project()` is
processed before the cache is consulted and **overrides** any `-DVAR=value` passed on
the command line. So even though the `build` CI job ran:
```
cmake -DCMAKE_C_COMPILER=mpicc ...
```
CMakeLists.txt silently forced `gcc`, which then failed to find `mpi.h` from the
parallel HDF5 headers.

### Fixes Applied
1. **`CMakeLists.txt`**: Changed compiler defaults from plain `set()` to
   `set(... CACHE STRING "...")` **without `FORCE`**. Without FORCE, if the cache
   already contains `-DCMAKE_C_COMPILER=mpicc` (from the command line), the `set()`
   is a no-op — the command-line value wins:
   ```cmake
   set(CMAKE_C_COMPILER   "gcc" CACHE STRING "C compiler")
   set(CMAKE_CXX_COMPILER "g++" CACHE STRING "CXX compiler")
   ```

2. **`build_platform_aware.sh`**: When parallel HDF5 is detected and `mpicc` is
   available, now also adds `-DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx`
   to `CMAKE_EXTRA_FLAGS`. This ensures the compiler reaches cmake as a **cache
   entry** (not just a `CC` env var, which the old `set()` would shadow):
   ```bash
   export CMAKE_EXTRA_FLAGS="${CMAKE_EXTRA_FLAGS} -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx"
   ```

## Notes
- The serial build uses stub implementations for all MPI functions
- Serial mode is single-process only (no parallel HDF5 I/O)
- Suitable for development, debugging, and small-scale testing
- In CI, `USE_MPI=no` still works with parallel HDF5: `mpicc` is used as compiler
  (satisfying `#include <mpi.h>` from HDF5 headers) while `-DUSE_MPI=0` disables
  MPI code paths in our project
