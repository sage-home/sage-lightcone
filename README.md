# TAO LIGHTCONE CLI #

## Description

These are the science modules, designed to carry out the HPC operation of extracting light-cones from "KD TREE" indexed simulations.

## Dependencies

 * Boost - A set of efficient C++ extensions.
 * MPI - Message Passing Interface implementation. Tested with OpenMPI.
 * HDF5 - Hierarchical Data Format libraries. Currently tested with version 1.14.0.
    (Must be compiled with OpenMPI or compilation errors will be encountered)
 * GSL - The GNU Scientific Library.
 * PugiXML - A C++ XML parsing library.
 

## Building

To build on Swinburne's nt.swin.edu.au after cloning from the repository just run
```bash
source ./setup.sh
./run_build.sh
```

To build in other environments refer to the .gitlab-ci.yml on how this is done within a CI/CDD framework. The science modules use cmake as a build system. As with any package
that has numerous dependencies, we recommend installing as many of them
as possible using your system's package management software. Note that the source to a compatible puxixml library is already included under the dep directory.

## Testing

The script test_cli_validation.sh will test that the build has succeeded and that a suite of validation tests correctly inform the user of errors in the command line arguments.

## Installing

```bash
cd bin
make install
```

The script lightcone.sh is installed under the runtime directory and allows a user to easily run the cli-lightcone tool on nt.swin.edu.au (for members of the oz114 project) ensuring that the right modules are loaded first.

# Structure of HDF5 format simulations with KD-TREE indexing

Using the h5ls tool we find the following (for the mini-millennium-public-sage-sed-kdtree.h5)

```
cosmology                Group
data                     Group
lightcone                Group
sed                      Group
snapshot_counts          Dataset {64}
snapshot_displs          Dataset {65}
snapshot_redshifts       Dataset {64/Inf}
```

1. cosmology holds the meta data for the simulation.  These are box_size, hubble constant, omega_l, and omega_m.
2. data contains the data for all the fields.  Each field has its own dataset and the order of the entries are the same for all datasets.  All the data belonging to the same snapshot will be contiguous (As described by snapshot_counts and snapshot_displ).  The order within a snapshot is such that they are in KD-TREE order using the x, y, and z coordinates to order them.  This means for example that all entries with an x less than half the the distance from the minimum x to the maximum x will be contiguous.
3. lightcone contains the KD-TREE indices for each snapshot.
4. sed contains (for each galaxy) the information required to compute the sed.
5. snapshot_counts, snapshot_displ, and snapshot_redshifts contains the indices to each snapshot and the associated redshift for that snapshot.

## HDF5 File Structure Diagram

```
mini-millennium-public-sage-sed-kdtree.h5
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
├── sed/
│   └── [SED data for each galaxy]
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
│  Galaxies 0-499       │  Galaxies 500-999                  │
│  (x < x_mid)          │  (x >= x_mid)                       │
│ ┌───────┬───────┐     │ ┌───────┬───────┐                  │
│ │ 0-249 │250-499│     │ │500-749│750-999│                  │
│ │y<y_mid│y>=y_mid     │ │y<y_mid│y>=y_mid                  │
│ └───────┴───────┘     │ └───────┴───────┘                  │
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

