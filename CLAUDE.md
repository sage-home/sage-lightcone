# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

TAO Lightcone CLI is a C++17/MPI/HDF5 workflow for converting SAGE semi-analytic galaxy formation model outputs into KD-tree indexed HDF5 files suitable for extracting cosmological lightcones. The pipeline handles both legacy SAGE binary format and modern SAGE HDF5 format.

## Build System

### Platform-Aware Building

Use the platform-aware build script that automatically detects macOS vs HPC/Linux:

```bash
./build_platform_aware.sh
```

This script:
- Detects platform and loads appropriate environment (setup_mac.sh or setup.sh)
- Builds main CMake project in Debug mode
- Builds SAGE executable from sage-model submodule
- Places executables in `bin/`

### CMake Configuration

- **Default C++ Standard**: C++17 (configurable via `-DCMAKE_CXX_STANDARD=XX`)
- **Build Type**: Debug (default) or Release
- **Backend**: KD-tree indexing (always enabled via `-DBACKEND_KDTREE`)
- **Important Flags**:
  - `-DNO_PROGRESS`: Progress tracking disabled (currently broken)
  - `-DNO_CHECKPOINT`: Checkpointing disabled (currently broken)
  - `-DPARALLELHDF5`: Parallel HDF5 enabled

**Note**: Release mode with `-DNDEBUG` is known to break at runtime. Use Debug builds.

## Key Executables

### Workflow Executables (in bin/)

1. **sage** - SAGE semi-analytic galaxy formation model (from sage-model submodule)
2. **sage2h5** - Converts SAGE binary output to intermediate HDF5 format
3. **sageh5toh5** - Converts SAGE HDF5 output to intermediate HDF5 format (modern workflow)
4. **sageimport** - Python-based import tool (processes intermediate HDF5)
5. **dstreeinit** - Creates KD-tree indexed HDF5 from depth-first ordered galaxy data
6. **cli_lightcone** - Extracts lightcones from KD-tree indexed HDF5 files
7. **sageh5tokdtree** - Direct SAGE HDF5 to KD-tree conversion (future/experimental)

### Workflow Comparison

**Legacy Binary Workflow**:
```
SAGE (binary) → sage2h5 → sageimport → dstreeinit → cli_lightcone
```

**Modern HDF5 Workflow**:
```
SAGE (HDF5) → sageh5toh5 → sageimport → dstreeinit → cli_lightcone
```

**Target Workflow** (experimental):
```
SAGE (HDF5) → sageh5tokdtree → cli_lightcone
```

## Testing

### End-to-End Testing

```bash
cd tests/sage-model-tests
./run_test_hdf5.sh        # Test HDF5 workflow
./run_test_binary.sh      # Test binary workflow (baseline)
```

These scripts:
- Run SAGE model with Mini-Millennium simulation
- Convert outputs to KD-tree indexed format
- Extract test lightcone
- Plot snapnum verification

### Validation Testing

```bash
./tests/test_cli_validation.sh
```

Tests command-line argument validation for cli_lightcone.

## Architecture

### HDF5 File Formats

#### SAGE HDF5 Output (Input to sageh5toh5)
```
model_X.hdf5
├── Header/              # Cosmology parameters
├── Snap_0/              # Columnar datasets per snapshot
│   ├── BlackHoleMass    # 1D scalar fields: shape (n_galaxies,)
│   ├── BulgeMass
│   ├── Posx, Posy, Posz
│   ├── SAGETreeIndex
│   ├── SfrBulgeSTEPS    # 2D array fields: shape (n_galaxies, n_timesteps) - NOT YET SUPPORTED
│   ├── SfrDiskSTEPS     # 2D array fields: shape (n_galaxies, n_timesteps) - NOT YET SUPPORTED
│   └── ... (50+ fields)
├── Snap_1/
└── TreeInfo/
```

**Field Types**:
- **1D scalar fields** (shape `(n_galaxies,)`): One value per galaxy. Examples: `StellarMass`, `Posx`, `ColdGas`. **Currently supported**.
- **2D array fields** (shape `(n_galaxies, n_timesteps)`): Time-series data with multiple values per galaxy. Examples: `SfrBulgeSTEPS`, `SfrDiskSTEPS`. **Currently skipped** - see Known Issues.

#### Intermediate HDF5 (Output of sage2h5/sageh5toh5)
```
mybinarymillennium.h5
├── cosmology/           # box_size, hubble_constant, omega_l, omega_m
├── galaxies             # Compound dataset with all galaxy properties
├── snapshot_redshifts   # Redshift per snapshot
├── tree_counts          # Galaxy count per tree
└── tree_displs          # Cumulative displacement per tree
```

#### KD-Tree Indexed HDF5 (Output of dstreeinit, Input to cli_lightcone)
```
mymillennium-kdtree.h5
├── cosmology/           # box_size, hubble_constant, omega_l, omega_m
├── data/                # All fields spatially reordered within each snapshot
│   ├── Posx, Posy, Posz
│   ├── StellarMass
│   └── ... (50+ fields)
├── lightcone/
│   ├── snapshot000/     # KD-tree index for snapshot 0
│   │   ├── bounds       # Spatial bounds
│   │   ├── cell_counts  # Galaxies per KD-tree cell
│   │   ├── cell_offs    # Offset into data arrays (relative to snapshot)
│   │   └── splits       # KD-tree split structure
│   ├── snapshot001/
│   └── ...
├── snapshot_counts      # Galaxies per snapshot
├── snapshot_displs      # Cumulative offset per snapshot
└── snapshot_redshifts   # Redshift per snapshot
```

### Key Data Organization Principles

1. **Snapshot Ordering**: All data fields are organized in contiguous blocks by snapshot, indexed via `snapshot_displs` and `snapshot_counts`.

2. **KD-Tree Spatial Ordering**: Within each snapshot, galaxies are reordered spatially using KD-tree partitioning (x, then y, then z). This enables efficient spatial queries.

3. **Field Alignment**: All field arrays (Posx, Posy, StellarMass, etc.) maintain identical ordering and length. Index `i` refers to the same galaxy across all fields.

4. **Relative Indexing**: KD-tree `cell_offs` are relative to the snapshot start (add `snapshot_displs[N]` for absolute index).

### Core Library (libhpc)

Custom C++ library under `src/libhpc/` providing:
- **algorithm/kdtree.hh**: KD-tree indexing implementation
- **h5/**: HDF5 C++ wrappers (file, dataset, dataspace, group)
- **mpi/**: MPI communication abstractions
- **containers/**: Custom containers (CSR, fibre, range_map)
- **logging/**: Debug logging with conditional compilation
- **debug/**: Assertions, stack traces, instrumentation

## Field Naming Conventions

The pipeline preserves SAGE's original field names throughout the conversion process:

- **SAGE fields** use **CamelCase**: `StellarMass`, `Posx`, `Velx`, `SnapNum`, `ColdGas`, etc.
- **Computed fields** use **lowercase**: `globaltreeid`, `snapnum`, `breadthfirst_traversalorder`, etc.
- **Case-insensitive matching** in user-facing tools (e.g., `plot_lightcone.py`)
- **Mandatory field validation** at pipeline entry points (`sage2kdtree`, `kdtree_backend`)

See **[FIELD_NAMING.md](FIELD_NAMING.md)** for the complete list of:
- Mandatory SAGE fields required for the pipeline
- Computed fields created by `sage2kdtree`
- Field naming conventions and best practices
- Validation implementation references

**Important**: Field names are preserved exactly as they appear in SAGE HDF5 output. Do not normalize or convert field names - the pipeline handles both CamelCase (SAGE) and lowercase (computed) fields natively.

## Common Development Commands

### Run a single test
```bash
cd tests/sage-model-tests
source ../../setup_mac.sh  # or setup.sh for HPC
../../bin/cli_lightcone --dataset output_sage_hdf5/myhdf5millennium-kdtree.h5 \
  --decmin 0 --decmax 1 --ramin 0 --ramax 1 --zmin 0 --zmax 1 \
  --outdir output_sage_hdf5 --outfile test-lightcone.h5
```

### Inspect HDF5 files
```bash
h5ls -r output_sage_hdf5/myhdf5millennium-kdtree.h5  # Recursive listing
h5dump -n output_sage_hdf5/myhdf5millennium-kdtree.h5  # Dataset names only
h5dump -d /cosmology/box_size output.h5  # Specific dataset
```

### Build specific executables
```bash
cd bin
make cli_lightcone    # Build single executable
make clean && make -j  # Clean rebuild with parallel jobs
```

## Important Implementation Notes

### Platform Differences

- **macOS (Apple Silicon)**: Uses Homebrew dependencies, no `-march=native` flag, HDF5 from brew
- **HPC (Linux)**: Uses module system, MPI via mpicc, can use `-march=native`

### Submodule: sage-model

SAGE is included as a Git submodule. Initialize with:
```bash
git submodule update --init --recursive
```

The build system automatically compiles SAGE when `sage-model/Makefile` exists.

### HDF5 Parallel Requirements

HDF5 must be compiled with MPI support for parallel I/O. On macOS, install:
```bash
brew install hdf5-mpi
```

### Descendant Field Calculation

When converting SAGE output, the `Descendant` field must be calculated by tracking merger trees. See `sage2h5.cc` for reference implementation.

### Subsize Field for Lightcones

The `subsize` field in `lightcone/data` represents the total number of galaxies in a cluster/tree rooted at a galaxy. This is computed during tree traversal, not directly available in SAGE HDF5 output.

## Current Development Status

### Completed
- `sageh5toh5` tool for HDF5 input workflow
- Platform-aware build system
- End-to-end test scripts for both binary and HDF5 workflows

### In Progress (see SAGE_PLAN_SAGEH5_TO_KDTREE.md)
- Direct `sageh5tokdtree` conversion (bypassing intermediate steps)
- Proper subsize calculation from SAGE HDF5
- Validation against binary workflow baseline
- **Array field support**: Currently only 1D scalar fields are processed. Multi-dimensional array fields (e.g., `SfrBulgeSTEPS`, `SfrDiskSTEPS`) with shape `(n_galaxies, n_timesteps)` are skipped. Future work needed to flatten or properly handle these time-series fields.

### Known Issues
- Release mode (`-DNDEBUG`) breaks at runtime - use Debug builds only
- Progress tracking disabled due to MPI thread joining issues
- Checkpointing currently disabled
- **Array fields not supported**: Multi-dimensional fields (ndims > 1) are currently skipped during conversion. This includes time-series fields like `SfrBulgeSTEPS` and `SfrDiskSTEPS` that store values at multiple timesteps. Only 1D scalar fields (one value per galaxy) are processed.

## Code Style & Conventions

- C++17 standard with some C++14 compatibility
- CamelCase for dataset/field names in HDF5 (e.g., `StellarMass`, not `stellar_mass`)
- MPI rank 0 typically handles I/O coordination
- Extensive use of HDF5 parallel I/O with collective operations
- Debug logging controlled by compile-time flags (`-DNLOGDEBUG`, etc.)
