# Plan: Batch Satellite Reads in `_fetch_satellites()`

## Problem

`_fetch_satellites()` in `kdtree_backend.hh` reads satellite data with individual
point reads — one HDF5 read call per field per satellite:

```
for each satellite i:
    for each field:
        ds.read(buf + i, ..., 1, sat_abs_idx)   // one read per field per satellite
```

For a batch of N satellites and F fields this is O(N × F) HDF5 read calls. Each call
has fixed HDF5 overhead (dataset open is cached, but hyperslab selection and transfer
is not). With ~50 fields and batches of thousands of satellites this adds up on both
local and parallel filesystems.

## Root Cause

The absolute indices in `_pending_sats` are not sorted. Satellites belonging to
different centrals from different snapshots are interleaved in the order centrals were
encountered. Without sorting, contiguous reads are not possible without rearranging
the output buffer.

## Proposed Solution: Sort → Batch Read → Scatter

Restructure `_fetch_satellites()` as three phases per field:

### Phase 1 — Sort by absolute index
Build a `sort_order` permutation that orders the current batch's satellite indices
ascending:

```cpp
std::vector<unsigned> sort_order(count);
std::iota(sort_order.begin(), sort_order.end(), 0);
std::sort(sort_order.begin(), sort_order.end(), [&](unsigned a, unsigned b) {
    return _pending_sats[_pending_sat_pos + a] < _pending_sats[_pending_sat_pos + b];
});
```

### Phase 2 — Coalesce into contiguous runs and multi-read
Walk the sorted indices and group them into contiguous runs (consecutive absolute
indices). For each run, issue a single hyperslab read into a temporary buffer.
For non-contiguous gaps between runs, use HDF5 point selection
(`H5Sselect_elements`) to read all sorted indices in a single call per field.

The simplest correct implementation uses `H5Sselect_elements` for the full sorted
index list — one read per field for the whole batch regardless of contiguity:

```cpp
// sorted_indices: sat_abs_idx values in ascending order
std::vector<hsize_t> sorted_indices(count);
for (unsigned i = 0; i < count; ++i)
    sorted_indices[i] = _pending_sats[_pending_sat_pos + sort_order[i]];

// One read per field:
hpc::h5::dataspace file_space = ds.dataspace();
H5Sselect_elements(file_space.hid(), H5S_SELECT_SET, count, sorted_indices.data());
hpc::h5::dataspace mem_space(count);  // simple 1D
ds.read(sorted_buf.data(), native_type, mem_space, file_space);
```

### Phase 3 — Scatter back to batch order
The read fills `sorted_buf` in sorted order. Scatter to the batch buffer using the
inverse permutation:

```cpp
for (unsigned i = 0; i < count; ++i)
    field_buf[sort_order[i]] = sorted_buf[i];
```

## Scope

- **File**: `src/libtao/base/kdtree_backend/kdtree_backend.hh`
- **Function**: `_fetch_satellites()` (lines ~1098–1222)
- **Change**: inner loop body only — the outer structure (batch sizing, `_calc_fields`,
  cone filter, `_pending_sat_pos` advance) is unchanged
- **Both field types**: 1D scalar fields and 2D array fields need the same treatment.
  For 2D, `H5Sselect_elements` selects rows; the sorted buffer holds `count × n_cols`
  values and the scatter copies whole rows.

## Scratch Buffers

Reuse the same approach as the Phase 4 scratch buffers in `sage2kdtree.cc`:
add two member vectors grown on demand, never shrunk:

```cpp
std::vector<hsize_t>  _sat_sorted_indices;   // sorted abs indices for H5Sselect_elements
std::vector<unsigned> _sat_sort_order;       // permutation: sorted → batch order
std::vector<char>     _sat_sorted_buf;       // read target (sized to max field × batch)
```

## Expected Benefit

| Scenario | Before | After |
|----------|--------|-------|
| N sats, F fields | N × F reads | F reads (one per field, all sats) |
| 1000 sats, 50 fields | 50,000 reads | 50 reads |

On a parallel Lustre filesystem (OzStar) each HDF5 read has ~0.1–1 ms latency, so
50,000 → 50 reads is a significant wall-time saving when lightcones include many
satellites.

## Correctness Constraints

- `_calc_fields(true)` is called **after** the read loop and uses the batch buffer
  directly — scatter must write into the same positions expected by `_calc_fields`.
- `central_spatial_index` is set separately from the field read loop (line 1183) and
  uses `_pending_sat_pos + i` indexing — this is unaffected.
- The cone filter after `_calc_fields` iterates `0..count` — also unaffected.
- Sort and scatter are local to `_fetch_satellites()` with no MPI interaction.

## Steps

- [ ] Add scratch member vectors `_sat_sorted_indices`, `_sat_sort_order`,
      `_sat_sorted_buf` to `kdtree_backend` class
- [ ] Build `sort_order` permutation at the top of `_fetch_satellites()` once per call
- [ ] Replace the per-satellite per-field read loop with per-field
      `H5Sselect_elements` + single read + scatter for 1D fields
- [ ] Same for 2D array fields (row selection + scatter of whole rows)
- [ ] Validate: Mini-Millennium end-to-end test output identical to baseline
- [ ] Benchmark: compare wall time for a satellite-heavy lightcone before/after
