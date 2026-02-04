# Implementation Summary: Simplify sage2kdtree Process

**Date:** 2026-02-04
**Goal:** Simplify the `sage2kdtree` pipeline by removing Phase 2 (Traversal Metadata) and associated fields (`subsize`, `bfs_idx`, `dfs_idx`) and SED logic.

## Iteration 1: Remove Traversal Metadata & SED (Completed)

**File:** `src/apps/sage2kdtree.cc`

### 1. Phase 1 (SAGE → Depth-First Ordered)
*   **Removed** `subsize` logic from `TreeNode` struct.
*   **Removed** `subsize` summation logic from `_depthfirst_recurse`.
*   **Added** direct calculation and writing of `globaltreeid` and `localgalaxyid` in Phase 1 output (copying `global_index` and `local_index`).
*   **Registered** these new fields in `field_list` and the `computed_fields` exclusion list.

### 2. Phase 2 (Add Traversal Metadata)
*   **Removed** the entire `phase2_add_traversal_metadata` function.
*   **Removed** helper function `_dfs_visit`.
*   **Removed** the call to Phase 2 in the main execution flow.
*   **Outcome:** Pipeline no longer performs a full read/write pass to calculate traversal indices.

### 3. Phase 3 (Tree → Snapshot)
*   **Updated** input source to read directly from Phase 1 output (`_depthfirst_fn`) instead of the now-removed Phase 2 output (`_enhanced_fn`).
*   **Outcome:** Phase 1 output flows directly into Phase 3.

### 4. Phase 4 (Build KD-Tree)
*   **Removed** `_write_sed_data` function and its call.
*   **Updated** `computed_fields` exclusion list to remove references to `subtree_count`, `breadthfirst_traversalorder`, and `depthfirst_traversalorder` (since they are no longer generated).
*   **Outcome:** SED logic eliminated; field exclusions consistent with new pipeline data.

## Iteration 2: Remove Depth-First Intermediate File (Completed)

**Date:** 2026-02-04
**Goal:** Further optimize IO by removing the "Depth-First Ordered" intermediate file. The legacy pipeline read SAGE data into memory, sorted it by tree (depth-first), wrote it to disk, then read it back to sort by snapshot. The new approach skips the tree sort and writes directly to snapshot-organized groups.

### 1. New "Direct" Phase 1
*   **Implemented** `phase1_direct_sage_to_snapshot`:
    *   Reads `model_X.hdf5` files directly.
    *   Using `MPI_Allreduce` and `MPI_Exscan` to calculate global galaxy counts and offsets without a full pre-read.
    *   Iterates through snapshots available in the input files (e.g. `Snap_63` -> `snapshot063`).
    *   Writes directly to the intermediate `*-bysnap.h5` file.
    *   **Optimization**: Uses `std::vector<std::shared_ptr<hpc::h5::file>>` to manage open file handles, preventing double-free errors during vector resizing.
    *   **Fix**: Explicitly selects file dataspaces for `write` operations to ensure hyperslab selections are respected (fixing "Dataspace sizes don't match" errors).

### 2. Header & Cosmology Handling
*   **Manual Attribute Writing**: The direct phase now manually creates the `cosmology` group and writes required attributes (`SimulationBoxSize`, `HubbleParam`, `OmegaMatter`, `OmegaLambda`) using the HDF5 C API.
*   **Rationale**: This replaces the `H5Ocopy` from the intermediate file which previously carried this metadata. Explicit writing ensures downstream tools (Phase 4 / Lightcone) have the necessary context.
*   **Robustness**: The input validator now accepts lowercase attribute names (e.g., `omega_matter`), commonly found in some SAGE outputs.

### 3. Removed Phases
*   **Removed** `_depthfirst_ordering` and `_depthfirst_recurse`.
*   **Removed** legacy Phase 3 (which read the depth-first file).
*   The pipeline flow is now: `SAGE HDF5` -> `Snapshot Organized HDF5` -> `KD-Tree HDF5`.

### 4. Verification
*   **Script**: `verify_kdtree_output.sh`
*   **Result**: Validated against `sage2kdtree_four` (a copy of the legacy code).
*   **Outcome**: The new direct pipeline produces bit-exact content matches and identical KD-tree structures compared to the legacy approach, but with significantly fewer IO passes.
*   **Performance vs Legacy (Single Core, Millennium test)**:
    *   **Legacy**: 205s
    *   **New**: 14s
    *   **Speedup**: ~14.6x
