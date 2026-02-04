# Plan: Simplify sage2kdtree Process

**Goal**: Simplify the `sage2kdtree` pipeline by removing Phase 2 (Traversal Metadata) and associated fields (`subsize`, `bfs_idx`, `dfs_idx`) and SED logic, as they are no longer required.

## Status: COMPLETE

## Steps

### 1. Modify Phase 1 (SAGE → Depth-First Ordered) [DONE]
*   **File**: [src/apps/sage2kdtree.cc](src/apps/sage2kdtree.cc)
*   **Action**:
    *   Removed `subsize` calculation from `_depthfirst_recurse`.
    *   Ensured output format remains compatible with what Phase 3 expects (columnar).

### 2. Remove Phase 2 (Add Traversal Metadata) [DONE]
*   **File**: [src/apps/sage2kdtree.cc](src/apps/sage2kdtree.cc)
*   **Action**:
    *   Removed the `phase2_add_traversal_metadata` function entirely.
    *   Removed the call to this function in `main`.

### 3. Update Phase 3 (Tree → Snapshot) [DONE]
*   **File**: [src/apps/sage2kdtree.cc](src/apps/sage2kdtree.cc)
*   **Action**:
    *   Changed input file source: Use Phase 1 output (`_depthfirst_fn`) instead of Phase 2 output (`_enhanced_fn`).

### 4. Simplify Phase 4 (Build KD-Tree) [DONE]
*   **File**: [src/apps/sage2kdtree.cc](src/apps/sage2kdtree.cc)
*   **Action**:
    *   Removed `_write_sed_data` function and its call.
    *   Updated `computed_fields` exclusion set: Removed `subtree_count`, `breadthfirst_traversalorder`, `depthfirst_traversalorder`.

## Validation [DONE]
*   **Strategy**: "Dual Binary" Verification.
    *   Created `src/apps/sage2kdtree_four.cc` (restored from git) to act as a legacy baseline.
    *   Compiled both binaries (`bin/sage2kdtree`, `bin/sage2kdtree_four`).
*   **Execution**:
    *   Ran `verify_kdtree_output.sh` (wraps `sage2kdtree` and `sage2kdtree_four` on `millennium` dataset).
    *   Compared outputs using schema-aware Python script (`verify_kdtree.py`).
*   **Results**:
    *   Coordinates (`x`, `y`, `z`) match exactly.
    *   Refactor confirmed: `global_index` and `subsize` are successfully absent in the new output.

## Vision Alignment Review

The plan strongly aligns with and actively advances several core principles in `VISION.md`.

### Alignment Analysis

**1. Unified Processing Model (Principle 5)**
*   **Vision:** "Clear separation between tree traversal logic and physics calculations."
*   **Plan:** By removing "Phase 2," you are stripping away legacy traversal logic (`bfs_idx`, `dfs_idx`, `subsize`) that is no longer needed for the core "SAGE to KD-Tree" transformation. This simplifies the processing model significantly, leaving only the essential steps: Depth-First Ordering -> Snapshot Grouping -> KD-Tree Construction.

**2. Metadata-Driven Architecture (Principle 3) & Single Source of Truth (Principle 4)**
*   **Vision:** "Maintain all the meta data available" and "Single Source of Truth".
*   **Plan:** The plan actively reduces legitimate "noise". By removing derived/calculated fields (`subsize`, traversal orders) that were added in the middle of the pipeline but not part of the SAGE source of truth, you move closer to a pure translation of the source data. The removal of `_write_sed_data` (which appeared to be hard-coded or specific logic) also aligns with relying on the metadata present in the input files rather than hard-coded structures.

**3. Memory Efficiency (Principle 6)**
*   **Vision:** "Memory usage is bounded, predictable... efficient".
*   **Plan:** Eliminating an entire processing phase (Phase 2), which involved reading and re-writing the entire dataset to add metadata, is a massive win for I/O efficiency and reduces total runtime memory pressure. It removes a "full copy" step.

**4. SAGE output hdf5 as the single point of truth (Principle 1)**
*   **Vision:** "...eventual light-cone extraction maintains all the meta data available in the sage output hdf5 file."
*   **Plan:** The plan focuses on passing through `globaltreeid` and original data while dropping intermediate artifacts that were specific to previous implementations (SED, specific traversal indices). This respects the input SAGE file as the primary source of truth rather than augmenting it with derived data that isn't strictly necessary for the KD-tree output.