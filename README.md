# TAO LIGHTCONE CLI #

## Description

These are modules designed to carry out the workflow of converting sage binary output to a kd-tree indexed format suitable for input to the cli_lightcone executable to extract lightcones in a simple HDF5 format.

## Dependencies

 * C++17 compiler or later (GCC 7+, Clang 5+, Apple Clang 10+)
 * CMake 3.10 or later
 * Boost - A set of efficient C++ extensions (filesystem, program options).
 * HDF5 - Hierarchical Data Format libraries. Currently tested with version 1.14.0.
 * MPI - Message Passing Interface implementation. Tested with OpenMPI. (Optional)

## Set up: Clone repository with submodules

```bash
git clone --recurse-submodules https://gitlab.com/CAS-eResearch/external/tao-managed/tao-lightcone-cli.git
cd tao-lightcone-cli
```

## Building

On the HPC system at nt.swin.edu.au or on MacOS the build should be done as follows:

```bash
./build_platform_aware.sh
```

after which the executables will be under the bin director.


### To run the end to end test (From sage to kd-indexing to extracting a lightcone)

```bash
./verify_kdtree_output.sh
```
This end to end test runs the following processes based on the mini-millennium test dataset.
1. sage with the "sage_hdf5" option. The "sage_hdf5" output option is mandatory.
2. sage2kdtree which converts the sage_hdf5 output from sage to a KD-tree indexed for lightcone extraction.
3. cli_lightcone which extracts lightcones from the KD-tree indexed hdf5 dataset.
 

## Building on other environments

To build in other environments refer to the .gitlab-ci.yml on how this is done within a CI/CD framework. The science modules use cmake as a build system. As with any package
that has numerous dependencies, we recommend installing as many of them
as possible using your system's package management software.

## Command Line Reference

The two main executables in the pipeline each accept a set of command line arguments described below. Both can be run with `-h` / `--help` to print a summary. `sage2kdtree` converts SAGE HDF5 output into a KD-tree indexed file ready for lightcone extraction. `cli_lightcone` queries that file to produce a flat HDF5 lightcone catalogue for a given sky volume and redshift range.

### sage2kdtree

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `--sage` | `-s` | *(required)* | Path to the directory containing SAGE HDF5 output files (`model_N.hdf5`) |
| `--param` | `-p` | *(required)* | SAGE parameter file (`.par`) used to read cosmology and simulation settings |
| `--alist` | `-a` | *(required)* | SAGE expansion factor list file (maps snapshot numbers to scale factors) |
| `--output` | `-o` | *(required)* | Output KD-tree HDF5 file path |
| `--ppc` | | `1000` | Particles per KD-tree cell — controls spatial index granularity |
| `--centralgalaxies` | | off | Build the central galaxy satellite lookup index (CSR index) alongside the KD-tree, required for `cli_lightcone --centralgalaxies` |
| `--noarrays` | | off | Skip 2D array fields (e.g. `SFHMassBulge`, `SFHMassDisk`) from the SAGE output |
| `--verbose` | `-v` | `1` | Verbosity level: 0 = quiet, 1 = progress, 2 = info, 3 = debug |
| `--help` | `-h` | | Print help and exit |

### cli_lightcone

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `--dataset` | `-d` | *(required)* | Path to the KD-tree indexed HDF5 file produced by `sage2kdtree` |
| `--decmin` | | `0` | Minimum declination of the output cone (degrees) |
| `--decmax` | | `10` | Maximum declination of the output cone (degrees) |
| `--ramin` | | `0` | Minimum right ascension of the output cone (degrees) |
| `--ramax` | | `10` | Maximum right ascension of the output cone (degrees) |
| `--zmin` | | `0` | Minimum redshift |
| `--zmax` | | `1` | Maximum redshift |
| `--outfile` | `-o` | `output.hdf5` | Output lightcone HDF5 filename |
| `--outdir` | | `output` | Directory for output files |
| `--outfields` | | *(all fields)* | Space-separated list of fields to include in the output; omit to include all fields plus all calculated fields |
| `--centralgalaxies` | | off | Include all satellites whose central galaxy falls in the cone, regardless of the satellite's own sky position |
| `--unique` | | `0` | If `1`, produce a unique lightcone (no repeated galaxies); if `0`, produce a random lightcone |
| `--seed` | | `0` | RNG seed for random lightcone mode |
| `--filterfield` | `-f` | `stellar_mass` | Galaxy property field to apply a value filter on |
| `--filtermin` | | `0` | Minimum value for the filter field |
| `--filtermax` | | *(none)* | Maximum value for the filter field |
| `--verbose` | `-v` | off | Enable verbose output |
| `--debug` | | off | Enable debug output |
| `--version` | | | Print version number and exit |
| `--help` | `-h` | | Print help and exit |

## Testing

All test and validation scripts live in the `tests/` directory and can be run from any working directory.

| Script | Purpose |
|--------|---------|
| `first_run.sh` | Downloads the Mini-Millennium tree files and scale-factor list into `tests/sage-model-tests/input/` if not already present. Called automatically by the test scripts when input data is missing. |
| `test_sage_hdf5.sh` | **Primary end-to-end test.** Runs the full pipeline from scratch: builds executables if needed, downloads tree data if missing, runs SAGE, `sage2kdtree`, and `cli_lightcone`, then plots SnapNum and redshift verification plots. Supports `--centralgalaxies` to exercise that code path. |
| `test_sage_hdf5_mpi.sh` | Same as `test_sage_hdf5.sh` but runs `sage2kdtree` and `cli_lightcone` under MPI with `--np N` tasks. Confirms the parallel code path produces correct output. |
| `test_cli_validation.sh` | Tests command-line argument validation for `cli_lightcone`: exercises invalid ra/dec/redshift ranges and confirms the executable reports errors correctly. |
| `validate_kdtree_output.sh` | Runs `sage2kdtree` and `cli_lightcone` against pre-existing SAGE output (does not re-run SAGE) and writes a full log to `validate_kdtree_output.txt` at the project root. Use this to check correctness of the conversion and lightcone steps in isolation. |
| `validate_kdtree_output_mpi.sh` | Same as `validate_kdtree_output.sh` but runs under MPI with `--np N` tasks. |
| `setup_validate_kdtree_output.sh` | Prepares the environment for `validate_kdtree_output.sh`: builds executables if needed and ensures SAGE output exists. Called automatically by `validate_kdtree_output.sh`. |
| `validate_and_time_workflow.sh` | Runs the full pipeline end-to-end (SAGE → `sage2kdtree` → `cli_lightcone`), timing and profiling each phase with peak memory and disk usage recorded to CSV files. Generates a timestamped Markdown report in `benchmark_reports/`. |
| `validate_and_time_centralgalaxies_option.sh` | Validates the `--centralgalaxies` option by running the full pipeline three times with different flag combinations and confirming the expected galaxy count ordering between passes. |
| `profile_xctrace_cli_lightcone.sh` | macOS only. Profiles `cli_lightcone` using Instruments `xctrace` (Time Profiler). Produces a `.trace` file openable in Instruments.app. |
| `profile_xctrace_sage2kdtree.sh` | macOS only. Profiles `sage2kdtree` from both the `main` and `before` branches using `xctrace` and the `sample` command, producing `.trace` and plain-text call-tree files in `profile_output/`. |

# SAGE HDF5 Output Format Requirements

`sage2kdtree` reads SAGE HDF5 output files directly. This section defines what a conformant SAGE HDF5 output must contain.

## File Layout

Each SAGE output file (e.g. `model_0.hdf5`, `model_1.hdf5`, …) must have this top-level structure:

```
model_N.hdf5
├── Header/
│   └── Simulation/          # Cosmology attributes (see below)
├── Snap_0/                  # One group per snapshot
│   ├── Posx                 # Columnar datasets, shape (n_galaxies,)
│   ├── Posy
│   ├── ...
├── Snap_1/
├── ...
└── TreeInfo/                # (optional but expected by SAGE)
```

- Snapshot groups must be named `Snap_0`, `Snap_1`, … (consecutive integers starting at 0).
- Each snapshot group contains one HDF5 dataset per galaxy field.
- **1D scalar fields** have shape `(n_galaxies,)`.
- **2D array fields** (e.g. star-formation history bins) have shape `(n_galaxies, n_cols)` and are fully supported throughout the pipeline.

## Header/Simulation Attributes

The `Header/Simulation` group must contain the following scalar attributes (case-insensitive matching is applied):

| Attribute | Description |
|-----------|-------------|
| `BoxSize` | Simulation box size in Mpc/h |
| `Hubble_h` | Dimensionless Hubble parameter *h* |
| `Omega_m` | Matter density parameter Ω_m |
| `Omega_lambda` | Dark energy density parameter Ω_Λ |

If `Header/Simulation` is absent or any of these attributes is missing, `sage2kdtree` will emit a warning and attempt to read cosmology from the parameter file (`-p` argument) instead.

## Mandatory Galaxy Fields

The following fields must be present in every `Snap_N` group. `sage2kdtree` validates their presence at startup and aborts with a clear error message listing any that are missing. Field name matching is case-insensitive.

### Spatial coordinates (required for KD-tree indexing)
| Field | Description |
|-------|-------------|
| `Posx` | x position (Mpc/h, comoving) |
| `Posy` | y position (Mpc/h, comoving) |
| `Posz` | z position (Mpc/h, comoving) |

### Snapshot identification
| Field | Description |
|-------|-------------|
| `SnapNum` | Snapshot number this galaxy belongs to |

### Galaxy and tree identifiers
| Field | Description |
|-------|-------------|
| `GalaxyIndex` | Unique galaxy index within the tree |
| `CentralGalaxyIndex` | Index of this galaxy's central (= `GalaxyIndex` for centrals) |

### Velocities (required for lightcone redshift calculations)
| Field | Description |
|-------|-------------|
| `Velx` | x peculiar velocity (km/s) |
| `Vely` | y peculiar velocity (km/s) |
| `Velz` | z peculiar velocity (km/s) |

All other fields present in the snapshot groups are passed through the pipeline unchanged and will appear in the final KD-tree HDF5 output.

## Output Format Requirement

SAGE must be run with `OutputFormat = sage_hdf5`. Binary (`sage_binary`) output is not accepted by `sage2kdtree`.

# Structure of HDF5 format simulations with KD-TREE indexing

Using the h5ls tool we find the following (for the mini-millennium-public-sage-sed-kdtree.h5)

```
cosmology                Group
data                     Group
lightcone                Group
snapshot_counts          Dataset {64}
snapshot_displs          Dataset {65}
snapshot_redshifts       Dataset {64/Inf}
```

1. cosmology holds the meta data for the simulation.  These are box_size, hubble constant, omega_l, and omega_m.  These are mandatory and used to calculate comoving distances from galaxy redshifts during the light-cone extraction.
2. data contains the data for all the fields.  Each field has its own dataset. The order of the gaflaxiy entries are the same for all datasets.  All the data belonging to the same snapshot will be contiguous (As described by snapshot_counts and snapshot_displ).  The order within a snapshot is such that they are in KD-TREE order using the x, y, and z coordinates to order them.  This means for example that all entries with an x less than half the the distance from the minimum x to the maximum x will be contiguous.
3. lightcone contains the KD-TREE indices for each snapshot.
4. snapshot_counts, snapshot_displ, and snapshot_redshifts contains the indices to each snapshot and the associated redshift for that snapshot.

## HDF5 File Structure Diagram

```
mini-millennium-public-sage-sed-kdtree.h5
|-- SageOutputHeader -> Copied from sage output's Header group
├── cosmology/
│   ├── box_size
│   ├── hubble_constant  
│   ├── omega_l
│   └── omega_m
├── data/
│   ├── field1 (e.g., x_pos)     [N total galaxies]
│   ├── field2 (e.g., y_pos)     [N total galaxies]  
│   ├── field3 (e.g., z_pos)     [N total galaxies]
│   ├── field4 (e.g., mass)      [N total galaxies]
│   └── ...
├── lightcone/
│   ├── snapshot_0/
│   │   ├── kdtree_index_0       [start_idx, length] -> KD-TREE element 0
│   │   ├── kdtree_index_1       [start_idx, length] -> KD-TREE element 1
│   │   └── ...
│   ├── snapshot_1/
│   │   ├── kdtree_index_0       [start_idx, length] -> KD-TREE element 0
│   │   └── ...
│   └── ...
├── snapshot_counts     [64] -> Number of galaxies in each snapshot
├── snapshot_displs     [65] -> Cumulative displacement for each snapshot
└── snapshot_redshifts  [64] -> Redshift for each snapshot
```

### Data Organization by Snapshot

The data arrays are organized as contiguous blocks by snapshot, with each field having identically-length arrays:

```
Each Field Array (ALL have identical structure and length):

x_pos:   ┌─────────────────┬─────────────────┬─────────────────┬─────┐
         │   Snapshot 0    │   Snapshot 1    │   Snapshot 2    │ ... │
         │ [0 ... count0-1]│[count0...count1]│[count1...count2]│     │
         └─────────────────┴─────────────────┴─────────────────┴─────┘

y_pos:   ┌─────────────────┬─────────────────┬─────────────────┬─────┐
         │   Snapshot 0    │   Snapshot 1    │   Snapshot 2    │ ... │
         │ [0 ... count0-1]│[count0...count1]│[count1...count2]│     │
         └─────────────────┴─────────────────┴─────────────────┴─────┘

z_pos:   ┌─────────────────┬─────────────────┬─────────────────┬─────┐
         │   Snapshot 0    │   Snapshot 1    │   Snapshot 2    │ ... │
         │ [0 ... count0-1]│[count0...count1]│[count1...count2]│     │
         └─────────────────┴─────────────────┴─────────────────┴─────┘

mass:    ┌─────────────────┬─────────────────┬─────────────────┬─────┐
         │   Snapshot 0    │   Snapshot 1    │   Snapshot 2    │ ... │
         │ [0 ... count0-1]│[count0...count1]│[count1...count2]│     │
         └─────────────────┴─────────────────┴─────────────────┴─────┘

... (all other fields follow the same pattern)

         ALL arrays use the SAME KD-TREE ordering within each snapshot
>>>>>>> main

Snapshot Displacement Array (snapshot_displs):
Index:  0    1    2    3   ...
Value: [0, 1000, 2500, 4200, ...]
       ↑     ↑     ↑     ↑
   Start  Start  Start  Start
   Snap0  Snap1  Snap2  Snap3

Snapshot Count Array (snapshot_counts):
Index:  0    1     2     3   ...  
Value: [1000, 1500, 1700, 800, ...]
        ↑     ↑     ↑     ↑
     Galaxies in each snapshot
```

### KD-TREE Indexing Within Snapshots

Within each snapshot, ALL field arrays use the SAME KD-TREE spatial ordering:

```
For Snapshot N (e.g., 1000 galaxies from index 0-999):

All Field Arrays (x_pos, y_pos, z_pos, mass, etc.) have IDENTICAL ordering:
┌─────────────────────────────────────────────────────────────┐
│  Galaxies 0-499       │  Galaxies 500-999                   │
│  (x < x_mid)          │  (x >= x_mid)                       │
│ ┌───────┬───────┐     │ ┌───────┬───────┐                   │
│ │ 0-249 │250-499│     │ │500-749│750-999│                   │
│ │y<y_mid│y>=y_mid     │ │y<y_mid│y>=y_mid                   │
│ └───────┴───────┘     │ └───────┴───────┘                   │
└─────────────────────────────────────────────────────────────┘

ONE KD-TREE Index Set for Snapshot N (applies to ALL fields):
├── kdtree_index_0: [start=0, length=500]      -> Left half (x < x_mid)
├── kdtree_index_1: [start=500, length=500]    -> Right half (x >= x_mid)  
├── kdtree_index_2: [start=0, length=250]      -> LL quarter (x<x_mid, y<y_mid)
├── kdtree_index_3: [start=250, length=250]    -> LR quarter (x<x_mid, y>=y_mid)
├── kdtree_index_4: [start=500, length=250]    -> RL quarter (x>=x_mid, y<y_mid)
├── kdtree_index_5: [start=750, length=250]    -> RR quarter (x>=x_mid, y>=y_mid)
└── ... (continues recursively)

Example Query: "Get all galaxies in LL quarter of Snapshot N"
1. snapshot_displs[N] = 1000 (starting index for this snapshot)
2. kdtree_index_2 = [start=0, length=250] (relative to snapshot start)
3. Absolute indices: 1000 + 0 to 1000 + 249 = indices 1000-1249
4. Extract from ALL fields: x_pos[1000:1249], y_pos[1000:1249], mass[1000:1249], etc.
```

This structure allows efficient spatial queries by:
1. Using `snapshot_displs` to find the data range for a specific redshift/snapshot
2. Using KD-TREE indices within that snapshot to quickly locate galaxies in specific spatial regions
3. All field arrays (x_pos, y_pos, mass, etc.) maintain the same ordering for consistent indexing

## Output lightcone HDF5 format

`cli_lightcone` writes a single flat HDF5 file containing all galaxies that fall within the requested sky volume. There are no snapshot or spatial sub-groups — every dataset lives directly at the file root.

### File layout

```
output-lightcone.h5
├── SageOutputHeader/          # Copied verbatim from the input KD-tree file
│   └── ...                    # (which copied it from the original SAGE output)
├── LightconeOutputHeader/     # Attributes recording the query parameters
│   ├── dataset                # Path to the input KD-tree file
│   ├── ramin, ramax           # Right ascension bounds (degrees)
│   ├── decmin, decmax         # Declination bounds (degrees)
│   ├── zmin, zmax             # Redshift bounds
│   ├── outfile, outdir        # Output path components
│   ├── output_fields          # Comma-separated list of fields written
│   ├── filter_field, filter_min, filter_max  # Optional galaxy filter
│   ├── unique                 # Whether duplicate suppression was applied
│   └── rng_seed               # RNG seed used
├── Posx                       # Observer-frame x coordinate (Mpc/h) [N galaxies]
├── Posy                       # Observer-frame y coordinate (Mpc/h) [N galaxies]
├── Posz                       # Observer-frame z coordinate (Mpc/h) [N galaxies]
├── ra                         # Right ascension (degrees)            [N galaxies]
├── dec                        # Declination (degrees)                [N galaxies]
├── distance                   # Comoving distance from observer (Mpc/h)
├── redshift_cosmological      # Redshift from comoving distance
├── redshift_observed          # Redshift including peculiar velocity
├── sfr                        # Total SFR = SfrDisk + SfrBulge (Msun/yr)
├── SnapNum                    # Snapshot the galaxy is drawn from
├── StellarMass                # (and all other SAGE fields)
├── SFHMassDisk                # 2D array fields shape (N, n_cols) if present
└── ...                        # All other fields from the KD-tree /data/ group
```

`central_spatial_index` (absolute KD-tree index of the host central; -1 for centrals) is also written when `--centralgalaxies` mode is active.

### Relationship to the SAGE and KD-tree formats

| Aspect | SAGE HDF5 | KD-tree HDF5 | Lightcone output |
|--------|-----------|--------------|------------------|
| Galaxy grouping | Per snapshot (`Snap_N/`) | Per snapshot in `data/`, spatially ordered | **Flat** — one list of all matched galaxies |
| Spatial ordering | None | KD-tree order within each snapshot | Extraction order (tile × snapshot) |
| Field location | Inside `Snap_N/` groups | Inside `data/` group | **Root level** |
| Field names | CamelCase SAGE names | SAGE names + lowercase computed names | Same as KD-tree |
| 2D array fields | `(n_galaxies, n_cols)` | `(n_galaxies, n_cols)` in `data/` | `(N, n_cols)` at root |
| Cosmology | `Header/Simulation` attrs | `cosmology/` group | Via `SageOutputHeader` |
| Position fields | Box coordinates | Box coordinates | **Observer-frame** coordinates (overwritten by `_calc_fields()`) |
| Extra metadata | — | — | `LightconeOutputHeader` with query parameters |

### Dataset attributes

Each calculated field dataset (`ra`, `dec`, `distance`, `redshift_cosmological`, `redshift_observed`, `sfr`, `central_spatial_index`) carries `Description` and `Units` string attributes. SAGE pass-through fields carry these attributes only if they were present in the KD-tree file.

