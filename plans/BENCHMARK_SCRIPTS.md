# Benchmark Scripts Guide

## Overview

Three benchmark scripts are available for comparing different workflow implementations:

## 1. `benchmark_workflows.sh` - OLD vs NEW Comparison

Compares the OLD and NEW HDF5-based workflows.

### Workflows Compared

**OLD Workflow (HDF5):**
```
SAGE (HDF5) → sageh5toh5 → sageimport → dstreeinit → kdtree
```

**NEW Workflow (HDF5):**
```
SAGE (HDF5) → sage2kdtree → kdtree
```

### Usage

```bash
# Run both workflows
./benchmark_workflows.sh

# Run only OLD workflow
./benchmark_workflows.sh --old-only

# Run only NEW workflow
./benchmark_workflows.sh --new-only

# Keep intermediate files
./benchmark_workflows.sh --no-cleanup
```

### Output Files

- **OLD kdtree**: `output_sage_hdf5_benchmark/myhdf5millennium-kdtree.h5`
- **NEW kdtree**: `output_sage_hdf5_one_step_benchmark/myhdf5millennium-kdtree-onestep.h5`
- **Reports**: `benchmark_reports/benchmark_YYYYMMDD_HHMMSS.md`

### Purpose

- Compare OLD multi-step HDF5 workflow against NEW single-step workflow
- Measure performance improvements (speed, memory, disk I/O)
- Validate output consistency

---

## 2. `benchmark_new_vs_binary.sh` - NEW vs BINARY Comparison

Compares the NEW workflow against the traditional binary SAGE baseline.

**IMPORTANT**: Unlike the OLD vs NEW comparison, SAGE must run twice with different output formats:
- BINARY workflow requires binary SAGE output
- NEW workflow requires HDF5 SAGE output

Due to SAGE non-determinism, field values may differ slightly, but kdtree structures should match.

### Workflows Compared

**BINARY Workflow (Baseline):**
```
SAGE (binary) → sage2h5 → sageimport → dstreeinit → kdtree
```

**NEW Workflow:**
```
SAGE (HDF5) → sage2kdtree → kdtree
```

### Usage

```bash
# Run both workflows
./benchmark_new_vs_binary.sh

# Run only BINARY workflow
./benchmark_new_vs_binary.sh --binary-only

# Run only NEW workflow
./benchmark_new_vs_binary.sh --new-only

# Keep intermediate files
./benchmark_new_vs_binary.sh --no-cleanup
```

### Output Files

- **BINARY kdtree**: `output_sage_binary_benchmark/mybinarymillennium-kdtree.h5`
- **NEW kdtree**: `output_sage_hdf5_new_benchmark/myhdf5millennium-kdtree-new.h5`
- **Reports**: `benchmark_reports/new_vs_binary_YYYYMMDD_HHMMSS.md`

### Purpose

- Validate NEW workflow against traditional baseline
- Ensure NEW workflow produces compatible kdtree output
- Verify no regressions compared to established workflow

---

## 3. Individual Test Scripts

### `run_test_binary.sh`
Runs the traditional binary SAGE workflow (baseline).

```bash
./run_test_binary.sh
```

Output: `output_sage_binary/mybinarymillennium-kdtree.h5`

### `run_test_hdf5_benchmark.sh`
Runs the OLD HDF5 workflow.

```bash
./run_test_hdf5_benchmark.sh
```

Output: `output_sage_hdf5_benchmark/myhdf5millennium-kdtree.h5`

### `run_validate_and_time_workflow.sh`
Runs the NEW single-step workflow.

```bash
./run_validate_and_time_workflow.sh
```

Output: `output_sage_hdf5_one_step_benchmark/myhdf5millennium-kdtree-onestep.h5`

---

## Key Features

### Shared SAGE Output

**`benchmark_workflows.sh` (OLD vs NEW)** uses **shared SAGE output** to ensure fair comparisons:
- SAGE runs once before workflows
- Output saved to `shared_sage_output/`
- Both workflows process identical HDF5 data
- Eliminates SAGE non-determinism (see `SAGE_REPRODUCIBILITY_ISSUE.md`)

**`benchmark_new_vs_binary.sh`** cannot share SAGE output:
- BINARY workflow requires binary SAGE output
- NEW workflow requires HDF5 SAGE output
- SAGE runs twice (once for each format)
- Field values may differ due to SAGE non-determinism, but structures should match

### Validation

The `utils/validate_outputs.py` script compares kdtree files:
- Field name consistency
- Data value comparison
- Structure validation
- Generates detailed validation reports

### Benchmark Reports

Generated in `benchmark_reports/` with:
- Timing breakdown per phase
- Peak memory usage
- Disk I/O statistics
- Speedup calculations
- Validation results

---

## Workflow Comparison Matrix

| Feature | OLD (HDF5) | NEW | BINARY |
|---------|------------|-----|--------|
| **Input** | SAGE HDF5 | SAGE HDF5 | SAGE binary* |
| **Steps** | 4 (sageh5toh5 → sageimport → dstreeinit) | 1 (sage2kdtree) | 3 (sage2h5 → sageimport → dstreeinit) |
| **Intermediate files** | 3 files | 0 files | 2 files |
| **Field names** | SAGE CamelCase | SAGE CamelCase | SAGE CamelCase |
| **Units** | SAGE units (Myr) | SAGE units (Myr) | SAGE units (Myr) |
| **Output format** | Columnar HDF5 | Columnar HDF5 | Columnar HDF5 |

*Note: BINARY workflow uses binary SAGE output (`millennium.par`). `sage2h5` requires binary format.

---

## Expected Results

After recent fixes, all workflows should produce:

1. **Identical field names** matching SAGE output
2. **Identical units** (dT in Myr, not Gyr)
3. **Identical data values** (when processing same SAGE output)
4. **Compatible kdtree structures**

Differences should only reflect:
- Performance (speed, memory, disk I/O)
- Intermediate file generation
- Pipeline implementation details

---

## Troubleshooting

### Different SAGE outputs between runs

See `SAGE_REPRODUCIBILITY_ISSUE.md`. Solution: Use benchmark scripts which share SAGE output.

### Field name mismatches

Ensure latest code with field name mapping fixes. Rebuild executables.

### Unit conversion issues

Ensure dT unit conversions removed from sageh5toh5.cc and sage2h5.cc.

### Validation failures

Check validation report in `benchmark_reports/` for specific differences.

---

## Date

Documentation created: 2026-01-06
