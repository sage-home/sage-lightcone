# SAGE MPI Reproducibility Issue

## Problem Description
Discrepancies were observed in SAGE model outputs (specifically `GalaxyIndex` counts and galaxy properties) when running with different numbers of MPI tasks (e.g., np=1 vs np=4). The physical results of the simulation should be independent of the parallelization strategy (number of MPI ranks), but they were diverging.

## Investigation

### 1. Verification Tooling
To isolate the issue, we enhanced the testing infrastructure:
*   **`verify_kdtree_bysnap.py`**: A new Python script was created to compare HDF5 outputs. Unlike simple file comparison, this tool sorts datasets by `GalaxyIndex` before comparison. This eliminates false positives caused by non-deterministic row ordering (which is expected in parallel file writing) and highlights actual content differences.
*   **Shell Script Refactoring**: `run_test_sage_hdf5_mpi4.sh` and `benchmark_utils.sh` were updated to support flexible argument passing (like `--np`) and robust validation.

### 2. Source Code Analysis
We analyzed the C source code in `sage-model/src` to identify sources of non-determinism:
*   **Logic Check**: `distribute_weighted_forests_over_ntasks` in `forest_utils.c` appeared deterministic in how it assigns work.
*   **RNG Usage**: We searched for random number generation. We found usage of the standard C library `rand()` function in `sage-model/src/model_misc.c`, specifically inside `init_galaxy` for assigning the `FFBRegime` (Feedback Free Burst Regime).
*   **Missing Seeding**: We found no evidence of `srand()` or GSL RNG initialization occurring at the start of each forest/tree processing.

## Root Cause
The standard `rand()` function maintains a global state. The sequence of random numbers returned depends on the order of calls.

*   **Scenario A (Serial)**: The processor handles Forest 1, then Forest 2. The random numbers drawn for Forest 2 depend on the final state of the RNG after Forest 1 finishes.
*   **Scenario B (Parallel)**: Task 0 handles Forest 1. Task 1 handles Forest 2. Task 1 starts with a fresh (or default) RNG seed. The sequence of random numbers drawn for Forest 2 is now different than in the serial case.

This led to different physical decisions (e.g., FFB regime assignment) for the same galaxies depending on the MPI configuration, ultimately causing divergence in the output populations.

## Fix Implementation

We implemented a fix to enforce deterministic seeding at the **forest level**. This ensures that regardless of which MPI rank processes a forest, or in what order, the random number stream for that forest is identical.

**File**: `sage-model/src/sage.c`
**Function**: `sage_per_forest`

**Code Change**:
We injected the following initialization code at the very beginning of the `sage_per_forest` function:

```c
    // Seed the random number generator based on the forest number.
    // This allows for reproducibility across different MPI counts.
    int64_t global_tree_id = forest_info->original_treenr[forestnr];
    int32_t file_nr = forest_info->FileNr[forestnr];
    unsigned int seed = global_tree_id + file_nr * 100000;
    srand(seed);
```

## Next Steps
1.  The `sage` model has been recompiled with this change.
2.  Run the comparison tests (`run_test_sage_hdf5_mpi4.sh` or similar) to verify that `np=1` and `np=4` now produce identical science results (using `verify_kdtree_bysnap.py`).
