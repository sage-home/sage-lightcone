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
2. **sage2kdtree** - Direct SAGE HDF5 to KD-tree conversion (**primary workflow** — 4-phase pure C++)
3. **cli_lightcone** - Extracts lightcones from KD-tree indexed HDF5 files
4. **sage2h5** - Converts SAGE binary output to intermediate HDF5 format (legacy)
5. **sageh5toh5** - Converts SAGE HDF5 output to intermediate HDF5 format (legacy, superseded)
6. **sageimport** - Python-based import tool (legacy, superseded)
7. **dstreeinit** - Creates KD-tree indexed HDF5 from depth-first ordered galaxy data (legacy, superseded)

### Workflow Comparison

**Legacy Binary Workflow** (baseline reference only):
```
SAGE (binary) → sage2h5 → sageimport → dstreeinit → cli_lightcone
```

**Old HDF5 Workflow** (superseded):
```
SAGE (HDF5) → sageh5toh5 → sageimport → dstreeinit → cli_lightcone
```

**Current Workflow** (operational):
```
SAGE (HDF5) → sage2kdtree → cli_lightcone
```

**Central Galaxies Workflow** (operational, feature branch):
```
SAGE (HDF5) → sage2kdtree --centralgalaxies → cli_lightcone --centralgalaxies
```

## Testing

### End-to-End Testing

```bash
cd tests/sage-model-tests
./run_test_sage_hdf5.sh                        # Primary: sage2kdtree workflow
./run_test_sage_hdf5.sh --centralgalaxies      # With central galaxies feature
./run_test_hdf5.sh                             # Legacy HDF5 workflow (reference)
./run_test_binary.sh                           # Legacy binary workflow (baseline)
```

These scripts:
- Run SAGE model with Mini-Millennium simulation
- Convert outputs to KD-tree indexed format via `sage2kdtree`
- Extract test lightcone via `cli_lightcone`
- Plot SnapNum and redshift verification plots

### Validation Testing

```bash
./tests/test_cli_validation.sh
```

Tests command-line argument validation for `cli_lightcone`.

### Satellite Validation

Full central/satellite analysis for the `--centralgalaxies` feature is documented in:
```
plans/VALIDATE_SATELLITE.md
```

Key findings (Mini-Millennium, ra/dec/z ∈ [0,1]):

| Run | Total | Inside | Outside | T0 | T1+ |
|-----|-------|--------|---------|----|----|
| baseline (default) | 211,344 | 211,344 | 0 | 185,107 | 26,237 |
| `--centralgalaxies` | 211,283 | 211,278 | 5 | 185,107 | 26,176 |
| `--cg --includeorphansatellites` | 211,698 | 211,344 | 354 | 185,107 | 26,237 |

- Tile transformation is correct: same tile applied to satellites as their central
- Max within-tile central–satellite distance: 0.83 Mpc/h — physically plausible
- Baseline and `--cg --orphans` have identical inside sets (185,107 T0 + 26,237 T1)
- Baseline produces zero outside-cone galaxies; `--cg --orphans` allows 354 CSR
  satellites outside (by design — those satellites' centrals are in-cone)

**Note on central–satellite separation analysis**: The output `Posx/Posy/Posz` fields
contain observer-frame coordinates (overwritten by `_calc_fields()`), not raw SAGE box
coordinates. Also, the same galaxy appears multiple times in the output (once per
lightcone tile). When matching satellites to centrals by `(SnapNum, CentralGalaxyIndex)`,
always pick the **nearest** candidate in 3D observer space — the `/lightcone-sat-distance`
skill does this correctly. Naively keeping the last match per key produces false large
separations (~150 Mpc/h) from cross-tile confusion.

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

# Standard lightcone
../../bin/cli_lightcone --dataset output_sage_hdf5_one_step/myhdf5millennium-kdtree-onestep.h5 \
  --decmin 0 --decmax 1 --ramin 0 --ramax 1 --zmin 0 --zmax 1 \
  --outdir output_sage_hdf5_one_step --outfile test-lightcone.h5

# With central galaxies feature
../../bin/cli_lightcone --dataset output_sage_hdf5_one_step/myhdf5millennium-kdtree-onestep.h5 \
  --decmin 0 --decmax 1 --ramin 0 --ramax 1 --zmin 0 --zmax 1 \
  --centralgalaxies \
  --outdir output_sage_hdf5_one_step --outfile test-lightcone-centralgalaxies.h5
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

### Completed (see `plans/SAGE_PLAN_SAGEH5_TO_KDTREE.md`)
- `sage2kdtree` consolidated 4-phase pure C++ pipeline (operational, ~95% of vision)
  - Phase 1: SAGE HDF5 → depth-first ordered
  - Phase 2: traversal metadata (BFS/DFS orders, globaltreeid, subtree_count)
  - Phase 3: tree order → snapshot order (columnar per-snapshot groups)
  - Phase 4: KD-tree spatial indexing → final `*-kdtree.h5`
- Columnar storage throughout all intermediate phases
- Dynamic field discovery from SAGE HDF5 at runtime
- Input validation with clear error messages
- Lustre performance fixes (`H5D_FILL_TIME_NEVER`, batch I/O)
- Field name preservation: SAGE CamelCase flows through unchanged
- Platform-aware build system (`build_platform_aware.sh`)
- End-to-end test scripts for all workflows

### Completed (see `plans/CENTRAL_GALAXIES_INDEXING.md`)
- `--centralgalaxies` flag for `sage2kdtree`: builds CSR satellite lookup index
  alongside the normal KD-tree output (`centralgalaxies/snapshotNNN/` groups)
- `--centralgalaxies` flag for `cli_lightcone`: emits satellites automatically for
  every central that enters the cone, with same tile transformation applied
- `--includeorphansatellites` flag: also emits satellites whose central is not in cone
- `central_spatial_index` output field (absolute KD-tree index of host central; -1 for centrals)
- Graceful degradation when index absent; non-regression when flag not set
- Bug fix: `central_spatial_index` no longer included in output field list when
  `--centralgalaxies` is not active (was causing crash)
- **Baseline epoch coupling** (always-on): baseline mode now always loads the CSR index
  and couples each satellite's epoch to its central's epoch. Satellites whose central is
  in-cone are suppressed from the normal stream and re-emitted via CSR; others pass
  through unchanged. This ensures central–satellite pairs always share the same snapshot.
- **Baseline cone filter**: CSR satellites emitted in baseline mode are post-filtered
  against the cone bounds (ra/dec/distance) so no galaxy lands outside the query volume.
  `--centralgalaxies` mode intentionally bypasses this filter.
- Validated end-to-end with Mini-Millennium dataset (see `plans/VALIDATE_SATELLITE.md`)

### Completed (branch `feature/minimum-image-convention`)
- **Minimum-image convention (MIC) fix** in `src/libtao/base/kdtree_backend/kdtree_backend.hh`:
  Satellites emitted via CSR are now placed in the periodic image nearest to their central
  before observer coordinates are computed. Without this fix, a satellite and central on
  opposite sides of the simulation box (e.g., box-x≈0 vs box-x≈60 for box_size=62.5 Mpc/h)
  can appear ~60 Mpc/h apart in observer space depending on the tile translation.
  - **Key insight**: MIC must be applied AFTER the tile rotation+translation in `_calc_fields()`,
    not before. Applying it to raw box coordinates first is wrong because `_calc_fields()`
    re-applies modulo wrapping which undoes the correction.
  - **Implementation**: `_collect_satellites()` stores the central's fully-transformed
    (post-rotation, post-translation, post-origin-shift) position. `_calc_fields()` in
    satellite mode applies MIC in that transformed frame, after translating the satellite
    but before computing observer distance/ra/dec.
  - **Validated**: Mini-Millennium test (ra/dec/z ∈ [0,1]) gives max separation 3.4 Mpc/h,
    median 0.29 Mpc/h (using correct nearest-tile matching in analysis). The fix is needed
    for correctness in general — the test data happens to use tile translations that avoid
    the worst case, but other query geometries will trigger it.

### Remaining (plans/CENTRAL_GALAXIES_INDEXING.md — To Validate)
- [ ] Verify `central_spatial_index` values are correct for all satellites in output
- [ ] Test edge cases: snapshots with zero satellites, orphan galaxies
- [ ] Confirm MPI behaviour: `_write_central_galaxy_index()` is rank-0-only; verify
  `_sat_offs`/`_sat_list` loading is consistent across ranks in parallel runs

### Remaining (plans/SAGE_PLAN_SAGEH5_TO_KDTREE.md — Optional Polish)
- [ ] Unit tests for critical functions (depth-first ordering, subsize, KD-tree partitioning)
- [ ] Memory profiling at phase boundaries

### In Progress (see SAGE_PLAN_SAGEH5_TO_KDTREE.md)
- Direct `sageh5tokdtree` conversion (bypassing intermediate steps)
- Proper subsize calculation from SAGE HDF5
- Validation against binary workflow baseline
- **Array field support**: Currently only 1D scalar fields are processed. Multi-dimensional array fields (e.g., `SfrBulgeSTEPS`, `SfrDiskSTEPS`) with shape `(n_galaxies, n_timesteps)` are skipped. Future work needed to flatten or properly handle these time-series fields.

### Known Issues
- Release mode (`-DNDEBUG`) breaks at runtime — use Debug builds only
- Progress tracking disabled due to MPI thread joining issues
- Checkpointing currently disabled
- `_fetch_satellites()` issues individual point reads per field per satellite (performance
  only — correctness unaffected); future optimisation: sort and batch contiguous ranges
- **Array fields not supported**: Multi-dimensional fields (ndims > 1) are currently skipped during conversion. This includes time-series fields like `SfrBulgeSTEPS` and `SfrDiskSTEPS` that store values at multiple timesteps. Only 1D scalar fields (one value per galaxy) are processed.

## Code Style & Conventions

- C++17 standard with some C++14 compatibility
- CamelCase for dataset/field names in HDF5 (e.g., `StellarMass`, not `stellar_mass`)
- MPI rank 0 typically handles I/O coordination
- Extensive use of HDF5 parallel I/O with collective operations
- Debug logging controlled by compile-time flags (`-DNLOGDEBUG`, etc.)
