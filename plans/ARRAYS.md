# 2D Array Field Support

## Summary

SAGE HDF5 output contains two classes of fields:

- **1D scalar fields** — shape `(n_galaxies,)`: one value per galaxy. e.g. `StellarMass`, `Posx`.
- **2D array fields** — shape `(n_galaxies, n_cols)`: multiple values per galaxy. The Mini-Millennium SAGE output contains `SFHMassBulge` and `SFHMassDisk` (64 star-formation history bins each). These were silently skipped by the `ndims != 1` guard in Phase 1 of `sage2kdtree`.

This plan documents the implementation that carries 2D fields through the full pipeline:

```
SAGE HDF5 → sage2kdtree → KD-tree HDF5 → cli_lightcone → lightcone output HDF5
```

## Status: Completed (branch `feature/minimum-image-convention`)

Validated end-to-end with Mini-Millennium dataset:
- `SFHMassBulge` and `SFHMassDisk` both present in kdtree and lightcone output
- Shape `(1519086, 64)` in kdtree; `(211262, 64)` in lightcone output
- Non-zero values present in 34 k of 211 k lightcone rows — physically expected

## Implementation

### sage2kdtree.cc — Phase 1 (SAGE HDF5 → snapshot-organised HDF5)

**File**: `src/apps/sage2kdtree.cc`

Previously skipped fields with `ndims != 1`. Changed to:
- Skip only `ndims > 2` (genuinely unsupported ranks)
- Broadcast `ndims_val` (int) and `n_cols` (unsigned long long) from `global_info_rank` to all MPI ranks so collective dataset creation is consistent
- For `ndims == 2`: create `(total_gals, n_cols)` dataset; write via 2D hyperslab with `start = {current_offset, 0}`, `count = {n_rows, n_cols}`
- For `ndims == 1`: existing 1D path unchanged

### sage2kdtree.cc — Phase 4 (snapshot HDF5 → KD-tree HDF5)

Dataset discovery now detects `ndims` from the intermediate file and creates:
- `(total_gals, n_cols)` 2D datasets for array fields
- `(total_gals,)` 1D datasets for scalar fields

`_write_attributes()` uses `row_size = n_cols * elem_size` for row-wise KD-tree permutation (memcpy of full rows, not individual elements).

### kdtree_backend.hh / kdtree_backend.cc — reading

**File**: `src/libtao/base/kdtree_backend/kdtree_backend.hh` and `.cc`

- `_array_field_ncols` map populated in `connect()` by inspecting dataset dataspace dims
- `field_binding` struct extended with `n_cols` field (0 for scalar, >0 for 2D)
- `_ensure_field_cache()` detects 2D datasets and stores `n_cols` in the binding
- `_fetch()` branches on `binding.n_cols > 0`: reads via 2D hyperslab into `hpc::matrix<T>` from `_bat->vector<T>(field)`
- `_fetch_satellites()` reads one row at `sat_abs_idx` for 2D fields
- `init_batch` shadow (non-virtual, shadows `backend::init_batch`): pre-registers 2D fields as `batch::VECTOR` before the base class registers them, preventing the base class `set_scalar` from overwriting the VECTOR rank

### hdf5.hh — writing (cli_lightcone output)

**File**: `src/libtao/modules/hdf5.hh`

Three new/modified sections:

**First batch (dataset creation + write)**:
- Detect `batch::VECTOR` rank via `std::get<1>(bat.field(field)) == VECTOR`
- Create 2D chunked extendable HDF5 dataset via raw HDF5 API (`H5Screate_simple`, `H5Pcreate`, `H5Pset_chunk`, `H5Dcreate2`) with `maxdims = {H5S_UNLIMITED, n_cols}`
- Write via `_write_2d_batch(dset, bat, field, row_offset=0, fvtype)`
- Store `field_ncols` and `field_name` in `_dset_ncols` / `_dset_names` parallel lists

**Subsequent batches (`else if (alt)`)** — iterates `_dset_ncols`/`_dset_names` in parallel with `_dsets`:
- For 2D fields: open dataset by name via `H5Dopen2(_file.id(), ...)`, call `H5Dset_extent` with `{old_rows + thistime, n_cols}`, close, then call `_write_2d_batch`
- For scalar fields: existing `dset->set_extent()` + switch-on-type path unchanged

**`_write_2d_batch` helper** (new protected method):
```cpp
void _write_2d_batch(h5::dataset* dset, const tao::batch<real_type>& bat,
                     const std::string& field, hsize_t row_offset,
                     typename tao::batch<real_type>::field_value_type fvtype)
```
- Collects filtered galaxy rows from `bat.vector<T>(field)` (a `hpc::matrix<T>`) into a flat row-major buffer
- Selects 2D file hyperslab at `{row_offset, 0}` with count `{nrows, n_cols}`
- Creates 2D memory space `{nrows, n_cols}`
- Writes via `dset->write(buf.data(), h5::datatype::native_double/llong, mem_space, file_space)`
- Handles DOUBLE and LONG_LONG element types; no-ops if filter produces 0 rows

## Key Design Decisions

**Why `H5Dset_extent` via raw HDF5 in the alt branch?**
`hpc::h5::dataset::set_extent(hsize_t)` asserts `ndims == 1`. There is no 2D variant. Since `dataset::_id` is protected, the workaround is to re-open the dataset by name (`H5Dopen2`), resize via `H5Dset_extent`, and close. The already-open `dset` pointer reflects the new extent on the next `dataspace()` call because HDF5 tracks extent state at the file level.

**Why `batch::VECTOR` pre-registration in `init_batch`?**
`backend::init_batch` iterates `qry.output_fields()` and calls `bat.set_scalar(field, type)`. `batch::set_scalar` and `batch::set_vector` both no-op if the field is already registered. So `kdtree_backend::init_batch` must pre-register 2D fields as VECTOR *before* calling the base class, or the base class would register them as SCALAR first and the VECTOR registration would be a no-op.

**Why MPI broadcast of `ndims` and `n_cols` in Phase 1?**
All ranks must create the HDF5 dataset collectively. Only `global_info_rank` reads the SAGE input and knows `ndims`/`n_cols`. These are broadcast before the collective dataset creation call.

**Row-wise permutation in Phase 4 KD-tree reordering**
`_write_attributes()` uses `row_size = n_cols * elem_size` and `memcpy` of full rows so that the KD-tree spatial reordering treats each galaxy's full array as an atomic unit — the `n_cols` values for a galaxy stay together.

## Files Changed

| File | Change |
|------|--------|
| `src/apps/sage2kdtree.cc` | Phase 1: support ndims==2 datasets; broadcast ndims/n_cols; 2D write. Phase 4: detect ndims, create 2D datasets, row-wise permutation. |
| `src/libtao/base/kdtree_backend/kdtree_backend.hh` | `_array_field_ncols` map; `field_binding.n_cols`; `init_batch` shadow; 2D hyperslab read in `_fetch()`; 2D row read in `_fetch_satellites()`. |
| `src/libtao/base/kdtree_backend/kdtree_backend.cc` | Populate `_array_field_ncols` in `connect()`; implement `init_batch`. |
| `src/libtao/modules/hdf5.hh` | `_write_2d_batch` helper; `_dset_ncols`/`_dset_names` members; first-batch 2D dataset creation; alt-batch 2D extent+write. |

## Performance after the feature/arrays changes

| Version | Option     | time sage2kdtree | time cli | time total |
|---------|------------|------------------|----------|------------|
| before. | default.   | 11.3s            | 6.7s     | 18.0s      |
| array.  | default.   | 31.6s            | 10.2s    | 41.8s      |
| array.  | --noarrays | 15.1s            | 9.9s     | 23.8s      |

Confirmed in `main` branch as of March 2026 (Pass 1/2/3 across three repeated runs):

| Branch | Option       | sage2kdtree (P1/P2/P3) | total excl. SAGE (P1/P2/P3) |
|--------|--------------|------------------------|------------------------------|
| before | (none)       | 11.3 / 11.9 / 12.5 s  | 18.0 / 18.7 / 16.7 s        |
| main   | --noarrays   | 15.2 / 15.7 / 14.9 s  | 23.5 / 24.1 / 19.8 s        |

Goals:
- **G1**: `main --noarrays` ≈ `before` speed (currently ~3–4 s slower per pass)
- **G2**: `main` (default, includes 2D array fields) faster than 31.6 s

---

## Performance Analysis

### Observed regression: `--noarrays` vs `before`

With `--noarrays` the two branches process identical 1D fields, so they should be
equivalent. The ~3–4 s regression is consistent across all three passes and across
both sage2kdtree and cli\_lightcone, which points to structural overhead rather than
a single operation.

**Likely causes (from diff analysis):**

1. **Per-snapshot redundant HDF5 metadata queries in Phase 1**
   For each field in each snapshot, Phase 1 opens the source dataset twice:
   once to check `ndims` and once (twice!) to read the field type:
   ```cpp
   // ndims check – per field per snapshot:
   g.open(*open_files[0], group_in);        // H5Gopen
   hpc::h5::dataset ds = g.dataset(fname);  // H5Dopen1
   hpc::h5::dataspace sp = ds.dataspace();  // H5Dget_space

   // type check – same field, same snapshot, separately:
   g.open(*open_files[0], group_in);              // H5Gopen (second open!)
   H5T_class_t cls = g.dataset(fname).type_class(); // H5Dopen1 + close
   size_t sz = H5Tget_size(g.dataset(fname).datatype().id()); // H5Dopen1 again!
   ```
   These metadata are **constant across snapshots** — `ndims` and the field type
   never change. With 64 snapshots × ~52 fields this is ~9,984 redundant HDF5
   open/query operations that could be reduced to ~156 (one pass over all fields
   from the first non-empty snapshot).

   On a warm SSD each H5Dopen/H5Gopen call is ≤5 µs, so 10k calls ≈ 50 ms —
   negligible. On a cold filesystem or with HDF5 metadata cache pressure (as seen
   during first-time runs on the Mini-Millennium test data) each call can take
   100–500 µs, pushing total overhead to 1–5 s. This matches the observed
   regression magnitude and explains why the gap is consistent across passes
   (each pass starts with a warm cache for subsequent fields but cold for the
   first pass).

2. **Redundant `dset.dataspace()` call inside the write loop (hot path)**
   Inside the per-file write loop, for each 1D field write the code calls
   `dset.dataspace()` twice: once for the bounds check and once for the hyperslab
   write. The bounds check invokes `H5Dget_space` on the *output* dataset:
   ```cpp
   hpc::h5::dataspace check_space = dset.dataspace();  // H5Dget_space – bounds check
   hsize_t actual_dataset_size = check_space.size();
   // ... then ...
   hpc::h5::dataspace file_space = dset.dataspace();   // H5Dget_space – write
   ```
   The bounds check is logically guaranteed by construction (`total_gals` matches
   the dataset size), so the first call is pure overhead. With 64 × 52 × n_files
   iterations this adds up.

3. **Two extra `MPI_Bcast` calls per non-skipped field per snapshot**
   Added for `ndims_val` and `n_cols` after the `skip_field` check. On a single
   MPI rank these are no-ops but still involve a function call and memory access.
   Impact is small (~few ms total) but adds up.

4. **Phase 4 `_write_attributes`: extra `simple_extent_num_dims()` per field per snapshot**
   `main` calls `field_ds.dataspace()` + `simple_extent_num_dims()` for every
   field in every snapshot even when `--noarrays` guarantees all fields are 1D.
   With 64 × 52 calls at ~1–10 µs each, this adds 3–30 ms.

5. **cli\_lightcone regression: `init_batch` overhead**
   The new `kdtree_backend::init_batch` shadow iterates `qry.output_fields()` and
   does two `std::unordered_map::find` calls per field (once for the resolved name,
   once for lowercase). With many output fields this adds O(n\_fields) work per
   batch setup. More importantly, `connect()` now calls `dsp.simple_extent_num_dims()`
   for every field in the KD-tree file to populate `_array_field_ncols`. Even with
   `--noarrays` the 1D fields are all inspected.

### Why `before` is faster despite identical HDF5 call counts per snapshot

The `before` branch has the same per-snapshot metadata call pattern. The difference
is that those numbers were measured with the test dataset already cached in the OS
filesystem cache from prior SAGE runs, while the `main` timings include a cold-start
run of SAGE (7.4 s vs 4.0 s for SAGE) suggesting a cooler filesystem cache state
throughout. Profiling on a fully warm or RAM-backed filesystem would clarify the
true code-level contribution. Regardless, the fixes below eliminate the redundant
calls and are clearly correct improvements.

---

## Performance Optimization Plan

### Fix 1 — Cache field metadata before the snapshot loop (highest impact)

**Target**: Phase 1 in `sage2kdtree.cc` (`phase1_direct_sage_to_snapshot`).

**Change**: Before the `for (snap ...)` loop, scan the **first non-empty snapshot**
once to build a metadata map:

```cpp
struct FieldMeta {
    int      type_code;  // 0=float, 1=int, 2=llong
    int      ndims;      // 1 or 2
    hsize_t  n_cols;     // 0 for 1D fields
};
std::unordered_map<std::string, FieldMeta> field_meta;

// --- pre-scan (rank == global_info_rank) ---
// Open first non-empty snapshot, iterate fields, populate field_meta.
// Broadcast the map to all ranks (same MPI pattern as field list broadcast).
```

Inside the snapshot loop replace the two separate `if (rank == global_info_rank)`
blocks (ndims detection + type detection) with a single map lookup:

```cpp
const auto& meta = field_meta.at(fname);
if (meta.ndims > 2 || (meta.ndims == 2 && _no_arrays)) {
    // skip
    MPI_Bcast(&skip_field, ...);
    continue;
}
// ndims_val, n_cols, type_code all come from meta — no HDF5 calls needed
```

**Reduction**: from `n_snapshots × n_fields × 3` H5Dopen calls to `n_fields × 3`
(one-time pre-scan). For 64 snapshots × 52 fields: 9,984 → 156.

**Expected savings**: 1–4 s depending on filesystem cache state; makes `--noarrays`
match `before` even on a cold filesystem.

### Fix 2 — Remove redundant bounds-check `dset.dataspace()` from write hot path

**Target**: Phase 1, 1D field write loop in `sage2kdtree.cc`.

**Change**: Remove the `check_space = dset.dataspace(); actual_dataset_size = ...`
block. The constraint `current_offset + count ≤ total_gals` is guaranteed by
construction: `total_gals` was computed by summing galaxy counts across all files
for this snapshot, and the loop iterates those same files in the same order.

```cpp
// Remove these lines from the 1D inner write path:
// hpc::h5::dataspace check_space = dset.dataspace();
// hsize_t actual_dataset_size = check_space.size();
// if (current_offset + count > actual_dataset_size) { throw ...; }
```

**Reduction**: Eliminates one `H5Dget_space` call per file per field per snapshot in
the hot write path (64 × 52 × n\_files calls).

### Fix 3 — Skip ndims check in Phase 4 when `--noarrays` is active

**Target**: `_write_attributes()` in `sage2kdtree.cc`.

**Change**: When `_no_arrays` is true the intermediate file contains only 1D fields
(array fields were filtered in Phase 1). Bypass `simple_extent_num_dims()`:

```cpp
int field_ndims = 1;
hsize_t n_cols = 1;
if (!_no_arrays) {
    hpc::h5::dataspace src_space = field_ds.dataspace();
    field_ndims = src_space.simple_extent_num_dims();
    if (field_ndims == 2) { /* get n_cols */ }
}
size_t row_size = n_cols * elem_size;
```

**Reduction**: Eliminates `n_snapshots × n_fields` `simple_extent_num_dims()` calls
when `--noarrays` is set (64 × 52 ≈ 3,328 calls).

### Fix 4 — Open field dataset once for both ndims and type detection (cleanup)

**Target**: Phase 1 metadata block in `sage2kdtree.cc` (applies even after Fix 1
for the pre-scan).

**Change**: The current code opens `g.dataset(fname)` three times per field
(ndims block + two opens in the type block). Consolidate into one open:

```cpp
if (rank == global_info_rank) {
    hpc::h5::group g;
    g.open(*open_files[0], group_in);
    hpc::h5::dataset ds = g.dataset(fname);   // single open
    hpc::h5::dataspace sp = ds.dataspace();
    int ndims = sp.simple_extent_num_dims();
    H5T_class_t cls = ds.type_class();         // reuse ds
    size_t sz = H5Tget_size(ds.datatype().id()); // reuse ds
    // ... set ndims_val, n_cols, type_code from above ...
}
```

### Fix 5 — Array mode: use HDF5 chunked datasets for 2D fields

**Target**: Phase 1 and Phase 4 dataset creation for 2D fields.

**Change**: Create 2D array-field datasets with chunking enabled:
- chunk shape `(ppc, n_cols)` where `ppc` is the points-per-cell parameter,
  matching the KD-tree cell size for aligned I/O in Phase 4.
- This avoids small partial-row writes that thrash the HDF5 metadata cache.

**Expected savings**: Significant for default mode (31.6 s → potentially < 20 s)
since SFHMassBulge and SFHMassDisk are 64-column datasets and the current code
writes each snapshot chunk with individual 2D hyperslab calls.

### Fix 6 — Array mode: read/write all snapshots of a 2D field together

**Target**: Phase 4 `_write_attributes()`.

**Change** (larger refactor): Instead of iterating field-by-field within each
snapshot call, refactor to iterate snapshot-by-snapshot within each field. This
allows reading the full column of a 2D field across all snapshots in one large
contiguous read (`total_gals × n_cols`), permuting all rows at once, and writing
in one pass. This trading of latency for throughput is the dominant win for large
2D fields.

---

## Implementation Priority

| Fix | Target | Difficulty | Expected impact |
|-----|--------|------------|----------------|
| 1   | `--noarrays` parity + array mode | Medium | High (1–4 s) |
| 2   | `--noarrays` parity | Easy | Medium (~0.5 s) |
| 3   | `--noarrays` parity | Easy | Low–medium |
| 4   | Both modes (cleanup) | Easy | Low |
| 5   | Array mode | Medium | High (estimated 5–10 s) |
| 6   | Array mode | Hard | Potentially highest |

Fixes 1–4 together should bring `main --noarrays` to parity with `before`.
Fix 5 is the primary lever for making default array mode competitive.
Fix 6 is a stretch goal.
