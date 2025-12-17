# SAGE HDF5 to KDTree Workflow Plan
This document outlines the plan to convert SAGE HDF5 results directly to the KDTree indexed format required by the TAO Lightcone CLI, replacing the legacy binary-based workflow.

## Goal
Replace the complex `sage2h5` (binary) -> `sageimport` (Python) -> `dstreeinit` (tree2sage) workflow with a more direct path: **SAGE HDF5 ->sageh5tokdtree**.

## Workflow Steps

### 1. Configure SAGE for HDF5 Output
Modify the SAGE parameter file (`.par`) to output in HDF5 format instead of binary.
- **Setting:** `OutputFormat sage_hdf5`
- **Output:** SAGE produces `model_zX.XXX_Y.hdf5` files containing `Snap_N` groups with columnar datasets.

### 2. Create a new process and investigate the SAGE binary workflow
Create `sageh5tokdtree` to accept the SAGE HDF5 format directly.
Use C++ code in dstreeinit.cc to investigate (but never change dstreeinit.cc but instead copy to sageh5tokdtree.cc if needed):

#### Sample SAGE HDF5 file structure
Knowing that the SAGE HDF5 file has the following sample structure:
h5ls model_0.hdf5
h5ls   model_0.hdf5
Header                   Group
Snap_0                   Group
Snap_1                   Group
...
Snap_9                   Group
TreeInfo                 Group

and the Snap_N where N=0 to number of snaphots
h5ls model_0.hdf5/Snap_0
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

1. How to create the cosmology dataset directly from the input SAGE HDF5 Header
2. How to process galaxies in one snapshot at a time and reorder the galaxy entries such that they are "kdtree indexed".  Each field in the data dataset is reorder in the same way according to the "kdtree index"
3.  It is key to understand how and where dsdtreeinit.cc does the kdtree indexing as we want to copy this approach for sageh5tokdtree.cc
4. How to populate the lightcone/data array of tuples giving the position of each galaxy (or only those that are not children - I'm not sure - a galaxy that aes a central_galaxy_id of -1 is in fact a central galaxy I believe)

### 3.  Verify some assumptions and establish unit tests for them
The run_test_binary.sh script will produce a valid HDF5 KDTREE indexed file from the SAGE binary workflow

Assumptions are:
1. The lightcone/data tuples have an entry for each non-child galaxy
2. We can kdtree index our SAGE HDF5 file such that for snapshot000 we end up with an identical lightcone/snapshot000 entry as we do with the run_test_binary.sh example.

### 4.  Implement the direct SAGE HDF5 to KDTREE indexed conversion

Create a new tool `sageh5tokdtree` (based on `dstreeinit.cc` logic but adapted for SAGE HDF5 input).

**Key understanding of the sage_binary workflow which forms the baseline.**
1.  sage creates a binary output
2.  sage2h5 takes the sage output
    a. Reorders the galaxies into depth first recreating gals[ii].global_index = first_global_index + ii; where ii starts at 0 and increments by 1.
    b. and creates a hdf5 output with a calculated extra field: descendant
3. src/sageimport_mpi_HDF2HDF/main.py createGalaxies shows the following fields are calculated
dynamic_dtype_description.append(('globaltreeid', '<i8'))
        dynamic_dtype_description.append(('breadthfirst_traversalorder', '<i8'))
        dynamic_dtype_description.append(('depthfirst_traversalorder', '<i8'))
        dynamic_dtype_description.append(('subtree_count', '<i8'))
        dynamic_dtype_description.append(('localgalaxyid', '<i4'))


**Key Implementation Details:**
1.  **Input Reading:** Read SAGE HDF5 `Snap_N` groups and CamelCase datasets.
2.  **Descendant calculation
'descendant', look at sage2h5.cc on how it is calculated.
3.  **Subsize Calculation:** Since SAGE HDF5 likely lacks `subsize`, implement logic to calculate it. `subsize` is the total number of galaxies in the galaxy cluster (tree) rooted at the current galaxy. This will likely require traversing the tree structure (using `Descendant`/`Progenitor` links or `SAGETreeIndex`).
4.  **KDTree Indexing:** Use the `hpc::kdtree` and `data_permuter` logic from `dstreeinit.cc` to reorder galaxies spatially.
5.  **Attribute Copying:** Copy attributes from input to output, preserving CamelCase names (e.g., `StellarMass`, `Vvir`) to avoid mapping errors.
6.  **Output Structure:** Generate the `/cosmology`, `/data` (with reordered CamelCase datasets), and `/lightcone` groups as per the KDTree format.

### 5. Validation
- Compare the output of `sageh5tokdtree` with the output of the legacy binary workflow (`run_test_binary.sh`).
- Verify `subsize` values match for equivalent trees.
- Ensure `lightcone/data` tuples are correct.
