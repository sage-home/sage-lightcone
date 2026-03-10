# Central Galaxies Indexing Plan

## Overview

This document describes the design and implementation of an additional indexing scheme for
central galaxy-based lightcone extraction. The feature extends the existing KD-tree pipeline
without disrupting any existing functionality.

**Key idea:** Central galaxies are processed normally through the KD-tree spatial index. When
`--centralgalaxies` is active, satellite galaxies are not placed in the lightcone by the
spatial query at all. Instead, whenever a central galaxy is found inside the lightcone cone,
all of its satellites are automatically included as well, regardless of their own positions.

---

## Definitions

- **Central galaxy**: a galaxy for which `GalaxyIndex == CentralGalaxyIndex`. The galaxy is its
  own group host.
- **Satellite galaxy**: a galaxy for which `GalaxyIndex != CentralGalaxyIndex`. Its host
  central is the galaxy whose `GalaxyIndex == satellite.CentralGalaxyIndex`.

**Important:** `GalaxyIndex` and `CentralGalaxyIndex` are **opaque identifiers** from the SAGE
merger tree, **not** array indices. A separate map must be built to translate them to array
positions.

---

## Flags and Options

| Option | Tool | Effect |
|--------|------|--------|
| `--centralgalaxies` | `sage2kdtree` | Build the CSR satellite lookup index alongside the normal KD-tree output |
| `--centralgalaxies` | `cli_lightcone` | Use the CSR satellite lookup index; include satellites automatically for every central that enters the cone |

Both flags are independent and additive — the existing KD-tree functionality is untouched when
neither flag is set.

---

## New HDF5 Datasets (written by `sage2kdtree --centralgalaxies`)

These datasets are written in addition to the existing `/data/` and `/lightcone/` groups.

```
<output>.h5
├── centralgalaxies/              # New group, present only when --centralgalaxies was used
│   ├── snapshot000/
│   │   ├── satellite_offsets     # uint64[n_gals+1], CSR row pointers
│   │   └── satellite_list        # uint64[n_satellites], CSR column values
│   ├── snapshot001/
│   └── ...
└── has_central_galaxy_index      # Root attribute (int), flag indicating index was built
```

### CSR Index Layout

For snapshot `s` with `N` galaxies in KD-tree spatial order:

- `satellite_offsets[i]` = start of galaxy `i`'s satellite list in `satellite_list`
- `satellite_offsets[i+1]` = end (exclusive)
- `satellite_list[satellite_offsets[i] .. satellite_offsets[i+1]]` = **snapshot-relative**
  spatial indices of all satellites belonging to galaxy `i`

For non-central galaxies (satellites), `satellite_offsets[i] == satellite_offsets[i+1]`
(empty range).

To convert a snapshot-relative satellite index to an absolute file index:
```
abs_idx = snapshot_displs[snap] + satellite_list[k]
```

---

## Index Construction (`sage2kdtree --centralgalaxies`)

Implemented in `sage2kdtree.cc::_write_central_galaxy_index()` (line 2145), called from
Phase 4 of the conversion pipeline after the KD-tree spatial ordering is established.

### Algorithm

Given:
- `idxs[spatial_idx]` = original (pre-reorder) array index of galaxy at spatial position
  `spatial_idx` (i.e., the permutation array from KD-tree construction)
- `GalaxyIndex[orig]` read from intermediate HDF5 snapshot group
- `CentralGalaxyIndex[orig]` read from intermediate HDF5 snapshot group

Steps:

1. **Build inverse permutation**: `original_to_spatial[orig] = spatial_idx`

2. **Build opaque-ID-to-row map**: `gi_to_row[GalaxyIndex[orig]] = orig`

3. **Classify and bucket satellites**: for each galaxy `j` in spatial order:
   - `orig_j = idxs[j]`
   - If `GalaxyIndex[orig_j] != CentralGalaxyIndex[orig_j]` → satellite
   - Look up `central_orig = gi_to_row[CentralGalaxyIndex[orig_j]]`
   - `central_spatial = original_to_spatial[central_orig]`
   - Append `j` to `satellite_bucket[central_spatial]`

4. **Convert buckets to CSR** and write to output HDF5.

### Trigger

`_write_central_galaxy_index()` is called only when `_build_central_galaxies_index == true`,
controlled by the `--centralgalaxies` flag on `sage2kdtree`.

After all snapshots are processed, the root attribute `has_central_galaxy_index` is written
on the output file so that `cli_lightcone` can detect the index.

---

## Lightcone Query Algorithm (`cli_lightcone --centralgalaxies`)

Implemented across `kdtree_backend.hh` and `kdtree_backend.cc`. The query iterator class
`kdtree_galaxy_iterator` is extended with central-galaxies state.

### Activation

In `kdtree_backend::connect()` (`kdtree_backend.cc` line 79):
```cpp
_central_galaxies_mode = global_cli_dict._central_galaxies;
```

If the mode is active but the file has no `has_central_galaxy_index` attribute, a warning is
printed and the feature degrades gracefully.

An extra output field `central_spatial_index` (type `LONG_LONG`) is added to the batch when
the mode is active.

### Per-Snapshot Initialization

In `kdtree_backend::load_snapshot()` (`kdtree_backend.cc` line 519):
- Reads `centralgalaxies/snapshot###/satellite_offsets` into `_sat_offs`
- Reads `centralgalaxies/snapshot###/satellite_list` into `_sat_list`

### Normal KD-tree Traversal (unchanged)

The KD-tree spatial traversal proceeds exactly as before. Galaxies whose 3D positions fall
outside the lightcone cone are masked. This stage operates on all galaxies (centrals and
satellites alike) during cone testing.

### `_collect_satellites()` — Post-cell processing

Called after each KD-tree cell batch is loaded and cone-filtered
(`kdtree_backend.hh` line 780):

1. Reads the `type` field for each galaxy in the batch:
   - `type == 0`: central galaxy
   - `type != 0`: satellite or orphan
2. Masks out (`_bat->mask(i)`) all satellites from normal output — they will not appear in
   the regular batch stream.
3. For each **unmasked central** (i.e., central inside the cone):
   - Computes `snap_rel = (abs_idx) - snapshot_displs[snap]`
   - Looks up satellites from CSR: `sat_offs[snap_rel]` to `sat_offs[snap_rel+1]`
   - Appends absolute satellite indices to `_pending_sats`
   - Appends central's absolute index to `_pending_sat_centrals`
4. Sets `central_spatial_index = -1` for all centrals in the batch.
5. Sets `_in_sat_flush = true` if pending satellites exist.

### `_fetch_satellites()` — Satellite emission

Called in place of normal batch loading while `_in_sat_flush == true`
(`kdtree_backend.hh` line 829):

1. Emits the next batch (up to `max_size`) of pending satellites.
2. For each satellite: reads all fields individually from `/data/` at the satellite's
   absolute index.
3. Sets `central_spatial_index` to the corresponding central's absolute index.
4. Calls `_calc_fields(satellite_mode=true)` to apply the same position transformation
   (tile rotation, redshift, ra/dec) that was applied to the central.
5. Advances `_pending_sat_pos`; clears pending lists when all satellites are flushed.

### Output Field: `central_spatial_index`

| Galaxy type | Value |
|-------------|-------|
| Central | `-1` |
| Satellite | Absolute KD-tree file index of its host central |

This field is only present in the output when `--centralgalaxies` is active.

---

## Data Flow Summary

```
SAGE HDF5 (Snap_N/)
    │  GalaxyIndex, CentralGalaxyIndex, Posx/y/z, ...
    ▼
sage2kdtree --centralgalaxies
    │
    ├── Phase 1-3: normal aggregation & spatial reordering
    │
    └── Phase 4: KD-tree construction + satellite index
        ├── Write /data/, /lightcone/  (unchanged)
        ├── _write_central_galaxy_index()
        │   ├── Read GalaxyIndex, CentralGalaxyIndex from bysnap file
        │   ├── Build gi_to_row (opaque ID → original row)
        │   ├── Build satellite_bucket[central_spatial] per snapshot
        │   └── Write centralgalaxies/snapshot###/{satellite_offsets, satellite_list}
        └── Write has_central_galaxy_index root attribute
    ▼
mykdtree.h5
    ▼
cli_lightcone --centralgalaxies
    │
    ├── load_snapshot(): load sat_offs, sat_list per snapshot
    │
    └── For each KD-tree cell:
        ├── Load batch, cone-filter positions
        ├── _collect_satellites():
        │   ├── Mask out all type!=0 (satellites) from main stream
        │   ├── For each unmasked central: queue its satellites
        │   └── Emit central batch with central_spatial_index=-1
        └── While _in_sat_flush:
            └── _fetch_satellites():
                ├── Read satellite fields from /data/ by absolute index
                ├── Set central_spatial_index = host central's abs idx
                └── _calc_fields(satellite_mode=true)
    ▼
lightcone.h5  (centrals + their satellites, all positioned correctly)
```

---

## Key Implementation Files

| File | Location | Role |
|------|----------|------|
| `types.hh` | `src/libtao/base/` | `cli_dict._central_galaxies` flag |
| `kdapplication.cc` | `src/apps/` | `--centralgalaxies` CLI option for `cli_lightcone` |
| `sage2kdtree.cc` | `src/apps/` | `--centralgalaxies` CLI option; `_write_central_galaxy_index()` (line 2145) |
| `kdtree_backend.cc` | `src/libtao/base/kdtree_backend/` | `connect()`, `open()`, `load_snapshot()` |
| `kdtree_backend.hh` | `src/libtao/base/kdtree_backend/` | `_collect_satellites()` (line 780), `_fetch_satellites()` (line 829) |

---

## Implementation Status

### Completed

- [x] `cli_dict._central_galaxies` field in `types.hh`
- [x] `--centralgalaxies` CLI option in `kdapplication.cc` (cli_lightcone)
- [x] `--centralgalaxies` CLI option in `sage2kdtree.cc`
- [x] `_write_central_galaxy_index()` in `sage2kdtree.cc` — builds and writes CSR satellite index
- [x] `has_central_galaxy_index` root attribute written on output file
- [x] `kdtree_backend::open()` — detects `has_central_galaxy_index` attribute
- [x] `kdtree_backend::connect()` — enables mode and adds `central_spatial_index` field
- [x] `kdtree_backend::load_snapshot()` — loads `sat_offs`, `sat_list` per snapshot
- [x] `_collect_satellites()` — masks satellites, queues pending satellite lookups
- [x] `_fetch_satellites()` — emits satellite batches with field reads and position calc

### Remaining / To Validate

- [ ] End-to-end test with Mini-Millennium dataset using `--centralgalaxies` on both tools
- [ ] Verify satellite counts in output match expected values from SAGE input
- [ ] Verify `central_spatial_index` values are correct in output lightcone
- [ ] Verify that running without `--centralgalaxies` produces identical output to baseline
  (non-regression)
- [ ] Test edge cases: snapshots with zero satellites, orphan galaxies (no central found in
  `gi_to_row`), empty satellite_list datasets
- [ ] Confirm MPI behaviour: `_write_central_galaxy_index()` is rank-0-only; verify
  `_sat_offs`/`_sat_list` loading is consistent across ranks in parallel runs

---

## Known Limitations / Future Optimizations

### `_fetch_satellites()` reads one galaxy at a time

`_fetch_satellites()` (`kdtree_backend.hh:829`) issues individual point reads from `/data/`
for each satellite — one HDF5 read per field per satellite. For halos with large satellite
counts this can generate many small, non-contiguous reads, which is expensive on both local
and parallel (Lustre) filesystems.

A better approach would be to sort the pending satellite indices, group contiguous runs into
range reads, and issue one read per contiguous run per field. This is a pure performance
optimization with no correctness implications and can be deferred until profiling shows it is
a bottleneck.

---

## Non-Disruption Guarantees

- All new code paths are gated behind `_build_central_galaxies_index` (sage2kdtree) or
  `_central_galaxies_mode` (cli_lightcone).
- When neither flag is set the pipeline behaves exactly as before.
- The new `centralgalaxies/` HDF5 group is absent from files built without the flag,
  so existing files remain fully compatible.
- The `has_central_galaxy_index` root attribute is only written when the index is built,
  allowing cli_lightcone to detect gracefully whether the index is available.
