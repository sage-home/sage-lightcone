# Benchmark Script Changes for SAGE Output Sharing

## Summary

Modified benchmark scripts to ensure both OLD and NEW workflows use **identical** SAGE output data for fair comparison.

## Problem

SAGE produces non-deterministic outputs when run multiple times with identical parameters (see `SAGE_REPRODUCIBILITY_ISSUE.md`). This makes it impossible to do fair comparisons between workflows since they start with different input data.

## Solution

**Run SAGE once, share output between both workflows:**

1. `benchmark_workflows.sh` runs SAGE once before any workflow execution
2. SAGE output saved to `shared_sage_output/`
3. Individual workflow scripts use the shared output instead of running SAGE

## Modified Files

### 1. `benchmark_workflows.sh`
- **Added**: SAGE execution block (lines 47-84)
- Runs SAGE once before workflows
- Saves output to `shared_sage_output/`
- Both workflows now process identical SAGE data

### 2. `run_test_hdf5_benchmark.sh` (OLD workflow)
- **Modified**: PHASE 1 (lines 51-72)
- Changed from running SAGE to copying shared output
- Falls back to standalone SAGE run if shared output not found

### 3. `run_test_hdf5_one_step_benchmark.sh` (NEW workflow)
- **Modified**: PHASE 1 (lines 151-172)
- Changed from running SAGE to copying shared output
- Falls back to standalone SAGE run if shared output not found

## Usage

### Running both workflows (recommended):
```bash
./benchmark_workflows.sh
```
- SAGE runs once
- Both workflows use identical data
- Valid comparison results

### Running individual workflows:
```bash
./run_test_hdf5_benchmark.sh          # OLD workflow
./run_test_hdf5_one_step_benchmark.sh # NEW workflow
```
- Falls back to running SAGE standalone
- Each workflow may have different SAGE outputs
- Comparison less meaningful

## Benefits

✅ **Guaranteed identical inputs** - Both workflows process same SAGE data
✅ **Fair comparisons** - Differences reflect workflow, not SAGE randomness
✅ **Faster execution** - SAGE runs once instead of twice
✅ **Backward compatible** - Individual scripts still work standalone

## Validation

After this change, the kdtree outputs from OLD and NEW workflows should differ **only** due to:
- Conversion pipeline differences (OLD: sageh5toh5→sageimport→dstreeinit vs NEW: sage2kdtree)
- Field name formatting (now standardized to match SAGE output)
- NOT due to different SAGE input data

## Date

Changes implemented: 2026-01-06
