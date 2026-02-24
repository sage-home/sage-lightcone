# SAGE Reproducibility Issue

## Problem

SAGE produces non-deterministic outputs when run multiple times with identical inputs and parameters.

## Evidence

Tested with Mini-Millennium simulation using identical conditions:
- Same parameter file: `millennium_sage_hdf5.par`
- Same input tree files: `input/millennium/trees/trees_063.*`
- Same SAGE binary (commit: 177dd42)
- Runs executed 1 minute apart

### Observed Differences

Comparing `tdeplete` field in Snap_41 between two runs:
- **Max absolute difference**: 999.9
- **Mean absolute difference**: 6.57
- **Files are NOT bitwise identical** despite same inputs

## Root Cause (Suspected)

SAGE exhibits non-deterministic behavior likely due to:

1. **Internal parallelization** with non-deterministic task ordering (OpenMP/MPI)
2. **Floating-point operation reordering** in parallel loops
3. **No fixed random seed** - may use system time/process ID implicitly
4. **Numerical precision** - different execution paths can accumulate different rounding errors

## Investigation Results

✅ **Parameter files**: Identical (`millennium_sage_hdf5.par`)
✅ **Input tree files**: Identical (shared from `input/millennium/trees/`)
✅ **SAGE binary**: Identical (compiled once, used for both runs)
✅ **No explicit random seed**: Not found in parameter file or source code
❌ **Outputs differ**: Significant differences in physical quantities

## Current Workaround

For benchmark comparisons between OLD and NEW workflows:
- Run SAGE **once** and share the output between both workflow benchmarks
- This ensures both workflows process **identical** SAGE data
- Eliminates SAGE non-determinism from workflow comparison

## Implementation

Modified benchmark scripts (`benchmark_workflows.sh`) to:
1. Run SAGE only once in a shared location
2. Copy/link SAGE output to both workflow directories
3. Ensures OLD and NEW workflows process identical data

## Future Work

To achieve reproducible SAGE runs:
- [ ] Investigate SAGE parallelization settings (OpenMP threads, MPI ranks)
- [ ] Check if SAGE supports a `--seed` parameter or reproducibility mode
- [ ] Test single-threaded execution: `export OMP_NUM_THREADS=1`
- [ ] Report issue to SAGE developers if reproducibility is critical

## Date

Issue identified: 2026-01-06
Documented by: Claude Code
