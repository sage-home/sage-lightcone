# A plan to create a version of sage2h5 that takes input from sage_hdf5 format
This document outlines the plan to create a version of sage2h5.cc (called sageh5toh5.cc) that takes input from the sage_hdf5 format and produces the same outputs as sage2h5.cc does when the sage_binary format is used as input.
The idea is to replicate the run_test_binary.sh workflow but starting with the sage_hdf5 option for sage.  We shall name this script run_test_sage_hdf5_input.sh.
The only difference in the sequence of events between run_test_binary.sh and run_test_sage_hdf5_input.sh is  

## Goal
Replace the complex `sage2h5` (binary) -> `sageimport` (Python) -> `dstreeinit` (tree2sage) workflow with the more standard 'sageh5toh5' (hdf5) -> `sageimport` (Python) -> `dstreeinit` (tree2sage) workflow.

## Workflow Steps

### 0. Status Update (2025-12-07)
- **Completed:** `sageh5toh5` tool created and partially verified.
- **Completed:** `run_test_sage_hdf5_input.sh` created and runs successfully.
- **Verified:** `sageh5toh5` output (`myhdf5millennium.h5`) matches `sage2h5` output (`mybinarymillennium.h5`) in terms of `snapshot_redshifts` and galaxy counts.
- **Next:** Final cleanup and merging.

### 0. Continue with current state
sageh5toh5.cc runs and has been verified to produce identical `snapshot_redshifts` as the baseline `mybinarymillennium.h5`.
The integration test `run_test_sage_hdf5_input.sh` runs the full pipeline successfully.
My root directory is /Users/rseikel/code/tao-lightcone-cli/sage_integration.
The root has a bin directory where I can run make.
The tests run_test_binary.sh and run_test_sageh5.sh are in /Users/rseikel/code/tao-lightcone-cli/sage_integration/tests/sage-model-tests

### 1. Configure SAGE for HDF5 Output
Modify the SAGE parameter file (`.par`) to output in HDF5 format instead of binary.
- **Setting:** `OutputFormat sage_hdf5`
- **Output:** SAGE produces `model_zX.XXX_Y.hdf5` files containing `Snap_N` groups with columnar datasets.

### 2. Create a new tool `sageh5toh5`
Create `sageh5toh5` to accept the SAGE HDF5 format directly and output the intermediate HDF5 format used by `sageimport`.
This tool will be based on `sage2h5.cc` but will read from HDF5 instead of binary files.

#### Sample SAGE HDF5 file structure (Input)
Knowing that the SAGE HDF5 file has the following sample structure:
```
h5ls model_0.hdf5
Header                   Group
Snap_0                   Group
Snap_1                   Group
...
Snap_N                   Group
TreeInfo                 Group
```

And the Snap_N groups contain columnar datasets for each property:
```
h5ls model_0.hdf5/Snap_0
BlackHoleMass            Dataset {0/Inf}
BulgeMass                Dataset {0/Inf}
GalaxyIndex              Dataset {0/Inf}
...
SAGETreeIndex            Dataset {0/Inf}
...
```

#### Sample expected HDF5 file structure (Output)
The output of `sageh5toh5` must match the output of `sage2h5`.
Structure:
```
h5ls mybinarymillennium.h5
cosmology                Group
galaxies                 Dataset {N_Galaxies} (Compound Type)
snapshot_redshifts       Dataset {N_Snapshots}
tree_counts              Dataset {N_Trees}
tree_displs              Dataset {N_Trees + 1}
```

The `galaxies` dataset contains the same fields as the binary output, but in a single compound dataset.

### 3.  Verify some assumptions and establish unit tests for them
The `run_test_binary.sh` script produces a valid HDF5 file from the SAGE binary workflow. This result is the baseline we want to replicate.

Assumptions are:
1. We can map the columnar data from SAGE HDF5 (`Snap_N/Property`) to the compound `galaxies` dataset in the output.
2. We can reconstruct the tree structure (`tree_counts`, `tree_displs`) from the `SAGETreeIndex` in the input.
3. The resulting file should be identical (or effectively identical) to the one produced by `sage2h5`.

### 4. Integration
Create `run_test_sage_hdf5_input.sh` which:
1. Runs SAGE in HDF5 mode.
2. Runs `sageh5toh5` to convert the output.
3. Runs `sageimport` and `dstreeinit` as usual.
4. Compares the final result with the binary workflow.


#### Sample expected HDF5 KDTREE file structure
Knowing that the kdtree indexedd HDF5 file will have the following sample structure
Sample KDTREE HDF5 format:
h5ls mymillennium-kdtree.h5           
cosmology                Group
data                     Group
lightcone                Group
snapshot_counts          Dataset {64}
snapshot_displs          Dataset {65}
snapshot_redshifts       Dataset {64}

where snapshot_counts is the total number of galaxies in each snapshot
and snapshot_displs is an offset into where the galaxies for each snapshot start
and snapshot_redshifts is the redshift for each snapshot

and the data dataset:
h5ls mymillennium-kdtree.h5/data
BlackHoleMass            Dataset {0/Inf}
BulgeMass                Dataset {0/Inf}
CentralGalaxyIndex       Dataset {0/Inf}
CentralMvir              Dataset {0/Inf}
ColdGas                  Dataset {0/Inf}
Cooling                  Dataset {0/Inf}
DiskRadius               Dataset {0/Inf}
EjectedMass              Dataset {0/Inf}
GalaxyIndex              Dataset {0/Inf}
Heating                  Dataset {0/Inf}
HotGas                   Dataset {0/Inf}
IntraClusterStars        Dataset {0/Inf}
Len                      Dataset {0/Inf}
MetalsBulgeMass          Dataset {0/Inf}
MetalsColdGas            Dataset {0/Inf}
MetalsEjectedMass        Dataset {0/Inf}
MetalsHotGas             Dataset {0/Inf}
MetalsIntraClusterStars  Dataset {0/Inf}
MetalsStellarMass        Dataset {0/Inf}
Mvir                     Dataset {0/Inf}
OutflowRate              Dataset {0/Inf}
Posx                     Dataset {0/Inf}
Posy                     Dataset {0/Inf}
Posz                     Dataset {0/Inf}
QuasarModeBHaccretionMass Dataset {0/Inf}
Rvir                     Dataset {0/Inf}
SAGEHaloIndex            Dataset {0/Inf}
SAGETreeIndex            Dataset {0/Inf}
SfrBulge                 Dataset {0/Inf}
SfrBulgeZ                Dataset {0/Inf}
SfrDisk                  Dataset {0/Inf}
SfrDiskZ                 Dataset {0/Inf}
SimulationHaloIndex      Dataset {0/Inf}
SnapNum                  Dataset {0/Inf}
Spinx                    Dataset {0/Inf}
Spiny                    Dataset {0/Inf}
Spinz                    Dataset {0/Inf}
StellarMass              Dataset {0/Inf}
TimeOfLastMajorMerger    Dataset {0/Inf}
TimeOfLastMinorMerger    Dataset {0/Inf}
Type                     Dataset {0/Inf}
VelDisp                  Dataset {0/Inf}
Velx                     Dataset {0/Inf}
Vely                     Dataset {0/Inf}
Velz                     Dataset {0/Inf}
Vmax                     Dataset {0/Inf}
Vvir                     Dataset {0/Inf}
dT                       Dataset {0/Inf}
infallMvir               Dataset {0/Inf}
infallVmax               Dataset {0/Inf}
infallVvir               Dataset {0/Inf}
mergeIntoID              Dataset {0/Inf}
mergeIntoSnapNum         Dataset {0/Inf}
mergeType                Dataset {0/Inf}

and the lightcone dataset:
h5ls millennium-public-sage-ozstar-sed-kdtree.h5/lightcone
data                     Dataset {753745268}
snapshot000              Group
snapshot001              Group
snapshot002              Group
snapshot003              Group
snapshot004              Group
snapshot005              Group
snapshot006              Group
snapshot007              Group
...
snapshot063              Group

and the lightcone/snapshotNNN:
h5ls millennium-public-sage-ozstar-sed-kdtree.h5/lightcone/snapshot063
bounds                   Dataset {3}
cell_counts              Dataset {35849}
cell_offs                Dataset {35849}
splits                   Dataset {17924}
which represent the kdtree index - cell_offs is a relative offset (to the snapshot_displs)
cell_counts is the number of galaxies
splits is the balanced tree info

and some lightcone/data sample:
data                     Dataset {753745268}
    Data:
         {328.946746826172, 419.109222412109, 114.894584655762, 204815448, 1},
         {78.9146423339844, 104.139892578125, 166.457153320312, 17229474, 1},
         {28.453540802002, 313.152557373047, 157.291946411133, 129785797, 1},
         {134.044082641602, 292.632507324219, 19.2555103302002, 182868279, 1},
         {329.011596679688, 419.080413818359, 115.026695251465, 204815449, 2},
         {337.436309814453, 380.859527587891, 149.605285644531, 214979719, 1},
where the tuples are (I believe) {posx, posy, posz, some_index, subsize}
and subsize is the total number of galaxies that are in the galaxy cluster given by global_galaxy_index.
**Note:** `Len` in SAGE HDF5 is likely unrelated to `subsize`. We must investigate the `sage_binary` workflow (specifically `sage2h5` or `dstreeinit` legacy modes) to determine how `subsize` is calculated.

### 3.  Verify some assumptions and establish unit tests for them
The run_test_binary.sh script will produce a valid HDF5 KDTREE indexed file from the SAGE binary workflow.  This result is the baseline we want to replicate.

Assumptions are:
1. The lightcone/data tuples have an entry for each non-child galaxy
2. We can kdtree index our SAGE HDF5 file such that for snapshot000 we end up with an identical lightcone/snapshot000 entry as we do with the run_test_binary.sh example.
3.  We assume the field names of the input in the case of sage_binary are mappable to field names of the input in the case of the sage_hdf5.  The output hdf5 in both cases should be identical.


