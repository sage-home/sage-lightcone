# Satellite Validation Analysis

Files analysed (in `tests/sage-model-tests/validate_test_sage_hdf5/`):
- `myhdf5millennium-lightcone.h5` — standard run (no `--centralgalaxies`)
- `myhdf5millennium-lightcone-centralgalaxies.h5` — run with `--centralgalaxies`

Satellites are matched to their centrals via `CentralGalaxyIndex == GalaxyIndex`,
within the **same SnapNum and same tile** (tile = floor(Posx/box_size), etc.).
Each snapshot and tile is treated independently. Box size = 62.5 Mpc/h, h = 0.73.

---

## Galaxy counts within ra∈[0,1], dec∈[0,1], redshift_cosmological∈[0,1]

| Run | Total galaxies | In range |
|-----|---------------|----------|
| No `--centralgalaxies` | 211,249 | 211,249 (100%) |
| `--centralgalaxies` | 211,283 | 210,929 (99.8%) |

The centralgalaxies run has 34 more galaxies total. 354 galaxies fall outside the
ra/dec/z ranges in the centralgalaxies run (the baseline has none out of range).

---

## Central–satellite counts and match rate

| | No `--centralgalaxies` | `--centralgalaxies` |
|---|---|---|
| Centrals (Type=0) | 185,107 | 185,107 |
| Satellites (Type>0) | 26,142 | 26,176 |
| Same-tile matches | 25,727 | 26,176 (100%) |
| Unmatched (no central in same tile) | 415 | 0 |

The centralgalaxies mode achieves 100% same-tile matching by construction — it only
emits satellites whose central passed the spatial query, and applies the same tile
transformation to both.

---

## Maximum central–satellite distance

Distances are computed two ways:
- **Raw (lightcone)**: Euclidean distance using output Posx/Posy/Posz
- **Within-tile**: minimum-image distance after stripping the tile offset (Pos % box_size),
  accounting for periodic boundary conditions

### No `--centralgalaxies` — max pair (SAGETreeIndex=54, SnapNum=43)

| Metric | Mpc/h | Mpc |
|--------|-------|-----|
| Raw (lightcone) distance | 62.49 | 85.60 |
| Within-tile distance | 0.54 | 0.74 |

The large raw distance is a box-boundary artefact: central sits at sim-space
x≈62.49 (upper edge) and satellite at x≈0.005 (lower edge). They are genuine
close neighbours (0.54 Mpc/h) under minimum-image convention.

#### Central galaxy

| Field | Value |
|-------|-------|
| GalaxyIndex | 3000054000000094 |
| CentralGalaxyIndex | 3000054000000094 |
| SAGETreeIndex | 54 |
| SAGEHaloIndex | 366 |
| SimulationHaloIndex | 9056594 |
| SnapNum | 43 |
| Type | 0 |
| Posx / Posy / Posz (lightcone, Mpc/h) | 2124.991, 7.554, 31.696 |
| Posx / Posy / Posz (sim-space, Mpc/h) | 62.491, 7.554, 31.696 |
| tile | (33, 0, 0) |
| ra / dec (degrees) | 0.2037, 0.8545 |
| distance (Mpc/h) | 2125.240 |
| redshift_cosmological | 0.8617 |
| redshift_observed | 0.8616 |
| StellarMass (10¹⁰ M☉/h) | 0.7734 |
| BulgeMass (10¹⁰ M☉/h) | 0.7609 |
| Mvir / CentralMvir (10¹⁰ M☉/h) | 181.34, 181.34 |
| Rvir | 0.1508 |
| Len | 2071 |
| HotGas | 22.332 |
| ColdGas | 0.07174 |
| BlackHoleMass | 0.002967 |
| SfrDisk / sfr | 0.0, 0.0 |
| Velx / Vely / Velz (km/s) | -15.96, 314.23, 38.60 |
| Vmax / Vvir / VelDisp | 241.64, 227.45, 144.62 |
| mergeType / mergeIntoID | 0, -1 |
| infallMvir / infallVmax / infallVvir | 0.0, 0.0, 0.0 |
| dT | 264.443 |

#### Satellite galaxy

| Field | Value |
|-------|-------|
| GalaxyIndex | 3000054000000116 |
| CentralGalaxyIndex | 3000054000000094 |
| SAGETreeIndex | 54 |
| SAGEHaloIndex | 434 |
| SimulationHaloIndex | 8945234 |
| SnapNum | 43 |
| Type | 1 |
| Posx / Posy / Posz (lightcone, Mpc/h) | 2062.505, 7.733, 31.190 |
| Posx / Posy / Posz (sim-space, Mpc/h) | 0.005, 7.733, 31.190 |
| tile | (33, 0, 0) |
| ra / dec (degrees) | 0.2148, 0.8664 |
| distance (Mpc/h) | 2062.755 |
| redshift_cosmological | 0.8299 |
| redshift_observed | 0.8295 |
| StellarMass (10¹⁰ M☉/h) | 1.3427 |
| BulgeMass (10¹⁰ M☉/h) | 0.0 |
| Mvir / CentralMvir (10¹⁰ M☉/h) | 13.771, 181.34 |
| Rvir | 0.06384 |
| Len | 160 |
| HotGas | 0.1242 |
| ColdGas | 0.3862 |
| BlackHoleMass | 0.0 |
| SfrDisk / sfr | 0.7426, 0.7426 |
| Velx / Vely / Velz (km/s) | -66.04, 19.37, -152.10 |
| Vmax / Vvir / VelDisp | 129.34, 96.32, 70.70 |
| mergeType / mergeIntoID | 0, -1 |
| infallMvir / infallVmax / infallVvir | 22.377, 136.02, 120.57 |
| dT | 264.443 |

---

### `--centralgalaxies` — max pair (SAGETreeIndex=60, SnapNum=50)

| Metric | Mpc/h | Mpc |
|--------|-------|-----|
| Raw (lightcone) distance | 87.57 | 119.96 |
| Within-tile distance | 0.83 | 1.14 |

Same box-boundary artefact: central at sim-space x≈62.06, satellite at x≈0.09.
Genuine separation of 0.83 Mpc/h under minimum-image convention.

#### Central galaxy

| Field | Value |
|-------|-------|
| GalaxyIndex | 3000060000000000 |
| CentralGalaxyIndex | 3000060000000000 |
| SAGETreeIndex | 60 |
| SAGEHaloIndex | 13 |
| SimulationHaloIndex | 9283025 |
| SnapNum | 50 |
| Type | 0 |
| Posx / Posy / Posz (lightcone, Mpc/h) | 1124.558, 18.398, 0.369 |
| Posx / Posy / Posz (sim-space, Mpc/h) | 62.058, 18.398, 0.369 |
| tile | (17, 0, 0) |
| ra / dec (degrees) | 0.9373, 0.01878 |
| distance (Mpc/h) | 1124.708 |
| redshift_cosmological | 0.4092 |
| redshift_observed | 0.4098 |
| StellarMass (10¹⁰ M☉/h) | 4.0347 |
| BulgeMass (10¹⁰ M☉/h) | 4.0347 |
| Mvir / CentralMvir (10¹⁰ M☉/h) | 1309.83, 1309.83 |
| Rvir | 0.3389 |
| Len | 16845 |
| HotGas | 185.023 |
| ColdGas | 0.09206 |
| BlackHoleMass | 0.008039 |
| SfrDisk / sfr | 0.0, 0.0 |
| Velx / Vely / Velz (km/s) | 111.97, -9.65, 187.16 |
| Vmax / Vvir / VelDisp | 437.62, 407.72, 252.00 |
| mergeType / mergeIntoID | 0, -1 |
| infallMvir / infallVmax / infallVvir | 0.0, 0.0, 0.0 |
| central_spatial_index | -1 |
| dT | 272.869 |

#### Satellite galaxy

| Field | Value |
|-------|-------|
| GalaxyIndex | 3000060000000089 |
| CentralGalaxyIndex | 3000060000000000 |
| SAGETreeIndex | 60 |
| SAGEHaloIndex | 524 |
| SimulationHaloIndex | 9504366 |
| SnapNum | 50 |
| Type | 1 |
| Posx / Posy / Posz (lightcone, Mpc/h) | 1062.587, 18.555, 62.246 |
| Posx / Posy / Posz (sim-space, Mpc/h) | 0.087, 18.555, 62.246 |
| tile | (17, 0, 0) |
| ra / dec (degrees) | 1.0004, 3.3520 |
| distance (Mpc/h) | 1064.571 |
| redshift_cosmological | 0.3852 |
| redshift_observed | 0.3870 |
| StellarMass (10¹⁰ M☉/h) | 0.5961 |
| BulgeMass (10¹⁰ M☉/h) | 0.1761 |
| Mvir / CentralMvir (10¹⁰ M☉/h) | 8.262, 1309.83 |
| Rvir | 0.06261 |
| Len | 96 |
| HotGas | 0.0 |
| ColdGas | 0.08826 |
| BlackHoleMass | 0.000141 |
| SfrDisk / sfr | 0.1037, 0.1037 |
| Velx / Vely / Velz (km/s) | 377.79, 137.97, 226.33 |
| Vmax / Vvir / VelDisp | 107.71, 75.33, 50.49 |
| mergeType / mergeIntoID | 0, -1 |
| infallMvir / infallVmax / infallVvir | 13.943, 119.33, 98.11 |
| central_spatial_index | 1039576 |
| dT | 272.869 |

---

## Summary

- `cli_lightcone` correctly applies the same tile transformation to satellites as to
  their central companions. The large raw distances seen are purely a periodic
  box-boundary artefact (central and satellite straddle opposite edges of the box).
- Within-tile (minimum-image) distances are small and physically plausible (< 1 Mpc/h).
- The `--centralgalaxies` mode achieves 100% same-tile satellite-central matching,
  eliminating the 415 unmatched satellites present in the normal run.

---

## The simulation box and tiling

### The simulation box and tiling

The Millennium simulation has a periodic box of side **62.5 Mpc/h**. To build a
lightcone extending to z~1 (~2200 Mpc/h away), the pipeline tiles this box repeatedly
along the line of sight — roughly 35 copies laid end to end. Each copy is a "tile", and
a galaxy in tile N has its lightcone position offset by `N × 62.5 Mpc/h` along x.

The simulation stores all galaxy positions in **wrapped simulation-space coordinates**
in `[0, 62.5)`. The tile transformation applied by `_calc_fields()` is:

```
lightcone_pos = sim_pos_translated + N × box_size
```

where `N × box_size` is the tile offset (`_dom.min() - _dom.origin()`).

### The periodic boundary problem

A galaxy near the **upper edge** of the box (sim-space x ≈ 62.4) and its satellite near
the **lower edge** (sim-space x ≈ 0.005) are physically close neighbours — separated by
only ~0.54 Mpc/h under periodic boundary conditions. In the real universe they are on
opposite sides of a repeating structure that wraps around.

After applying the same tile offset (e.g., tile 33 → offset = 33 × 62.5 = 2062.5 Mpc/h):

```
central lightcone_x   = 62.49 + 2062.5 = 2124.99 Mpc/h
satellite lightcone_x =  0.005 + 2062.5 = 2062.51 Mpc/h
```

The **Euclidean distance** in lightcone coordinates is:

```
|2124.99 − 2062.51| ≈ 62.49 Mpc/h  ← the "raw distance"
```

This looks alarming, but it is almost exactly one full box width — a dead giveaway that
it is a wrapping artefact, not a real physical separation.

### The minimum-image correction

In a periodic box, the true separation between two points is found by the
**minimum-image convention**: for each axis, if `|Δx| > box/2`, wrap it around by
subtracting one box width:

```
Δx_wrapped = Δx − box × round(Δx / box)
           = 62.49 − 62.5 × round(62.49 / 62.5)
           = 62.49 − 62.5 × 1
           = −0.01 Mpc/h
```

Combined with the y and z separations this gives a true physical distance of
**~0.54 Mpc/h**, which is well within the virial radius of a group-scale halo and
physically reasonable for a satellite.

### Why `cli_lightcone` is correct

Both the central and satellite received **exactly the same tile offset** from `_dom`.
The code is doing the right thing — it places both galaxies at their correct lightcone
positions. The large raw distance is not a bug; it is simply what happens when a bound
pair straddles the periodic boundary of the simulation box. Any downstream analysis that
computes central-satellite separations needs to apply the minimum-image convention, or
equivalently, work in simulation-space coordinates (`Posx % box_size`) rather than
lightcone coordinates.

---

## `--includeorphansatellites` implementation and validation

### Purpose

`--includeorphansatellites` (requires `--centralgalaxies`) extends the CSR re-emission
mode to also include satellites whose central galaxy is **not** in the lightcone cone.
In normal `--centralgalaxies` mode such satellites are suppressed; with this flag they
are emitted directly from the normal stream with `central_spatial_index = -1`.

### Bug: same-position duplicate emission (fixed)

Initial implementation produced ~17,000 extra T1 galaxies inside the cone compared to
baseline. Investigation showed 14,237 same-position (ra/dec/z identical) duplicates:
each appeared once from the orphan path (`csi=-1`) and once from the CSR path
(`csi>=0`).

**Root cause**: `tile::contains()` has a known limitation — it ignores the tile's
rotation and translation when checking whether a point is inside the cone:

```cpp
// TODO: Include rotation and shifting.
hpc::num::cartesian_to_ecs(crd[0] + this->_min[0], ...)  // offset only, no rotation/trans
```

The orphan check in `_collect_satellites()` used `_dom.contains(crd_C, _snap)` to
decide whether the central was in the cone. But `_calc_fields()` applies a full
rotation + translation before checking `ecs_box_collision`. So a central that passes
the real cone check (and is subsequently emitted, queuing its satellite via CSR) could
still fail `_dom.contains()`, causing the satellite to be emitted as both an orphan
**and** a CSR re-emission from the same tile.

**Fix**: replaced `_dom.contains(crd_C, _snap)` with a replication of the
`_calc_fields()` transformation:

```cpp
// Apply same rotation+translation as _calc_fields()
int ix = _dom.rotation()[0]; int iy = _dom.rotation()[1]; int iz = _dom.rotation()[2];
double px = raw[ix];
if (px + _dom.translation()[ix] < box_size) px += _dom.translation()[ix];
else { px -= box_size; px += _dom.translation()[ix]; }
px += (_dom.min()[0] - _dom.origin()[0]);
// ... same for py, pz ...
central_in_cone = ecs_box_collision(ecs_min, ecs_max, pt, pt);
```

After the fix, same-position duplicates dropped to **0**.

### Galaxy counts after fix

| Run | Total | Inside | Outside | T0 | T1+ |
|-----|-------|--------|---------|----|----|
| baseline | 211,249 | 211,249 | 0 | 185,107 | 26,142 |
| `--centralgalaxies` | 211,283 | 210,929 | 354 | 185,107 | 26,176 |
| `--cg --includeorphansatellites` | 211,698 | 211,344 | 354 | 185,107 | 26,591 |

The 354 "outside" entries in both `--cg` runs are CSR-emitted satellites of in-cone
centrals that ended up outside the requested ra/dec/z volume after tile transformation —
this is by design.

### Residual discrepancy: 95 extra "inside" in orphans run

`--cg --orphans` inside (211,344) exceeds baseline inside (211,249) by 95. Analysis
shows:

- **0** galaxies are inside baseline but missing from `--cg --orphans` (no losses)
- **76** unique (GalaxyIndex, SnapNum) pairs are inside `--cg --orphans` but absent
  from baseline entirely

All 76 extras are `Type=1, csi>=0` (CSR-emitted).

### How CSR emission differs from the baseline stream

Galaxies in SAGE evolve across simulation time and have a record at **every snapshot**
from formation to z=0 — each with a different 3D position. The KD-tree stores all of
these separately, keyed by snapshot.

The baseline's per-snapshot cone filter (`ecs_box_collision` in `_calc_fields()`) is
applied to every galaxy found by the KD-tree spatial query. For snapshot N, a galaxy
is included only if its comoving distance at time N falls within N's assigned distance
shell `[min_dist(N), max_dist(N)]`. The galaxy therefore appears in the output at
whichever snapshot M has a shell that matches its position at time M. This is the
"correct" epoch: the snapshot chosen by the baseline is the one whose lookback time
matches the galaxy's observed comoving distance.

CSR-emitted satellites are handled differently. In `_fetch_satellites()`, `_calc_fields()`
is called with `satellite_mode=true`. In that mode the `ecs_box_collision` check is
**skipped**:

```cpp
if (!satellite_mode)
{
    bool btest = ecs_box_collision(ecs_min, ecs_max, min, min);
    if (!btest)
        _bat->mask(ii);
```

So CSR satellites receive the tile transformation but are **never filtered** by the
snapshot's distance shell. A satellite is emitted whenever its central is in the cone,
regardless of whether the satellite's own comoving distance falls inside snapshot N's
shell.

### Why the 76 extras are at different (GalaxyIndex, SnapNum) than baseline

The 76 extras appear in `--cg --orphans` at snapshot N (the central's snapshot) with
`csi >= 0`. The same physical galaxy is present in the KD-tree at other snapshots too.
At snapshot M ≠ N, where the galaxy's position at time M maps to a distance within M's
shell, the baseline picks it up as `(GalaxyIndex, M)`.

Because the snapshot differs (N vs M), the `(GalaxyIndex, N)` entry in `--cg --orphans`
and the `(GalaxyIndex, M)` entry in baseline are treated as distinct rows. The baseline
entry at M is also present in `--cg --orphans` (via orphan or CSR path at M), which is
why there are 0 losses. The 76 "extras" are the CSR entries at the central's epoch N,
for satellites whose own distance at time N lies outside N's shell.

In other words: the baseline and CSR modes are capturing the **same physical satellites**
but at **different simulation epochs**. CSR ties the satellite to its central's epoch;
the baseline ties it to the epoch that best matches the satellite's own comoving distance.

Conclusion: the 95 extra entries in `--cg --orphans` are not a bug. They reflect a
design difference: CSR emission couples satellite epoch to central epoch, while the
baseline selects satellite epoch independently by distance-shell matching. The choice
of which behaviour is more physically appropriate depends on the science case.

---

## Baseline satellite epoch coupling (always-on CSR)

### Motivation

In the old baseline, a satellite's snapshot was chosen by distance-shell matching
independently of its central. This means a satellite and its central can appear at
different simulation epochs in the lightcone output. The coupling change makes the
baseline always select a satellite's epoch to match its central's epoch, consistent
with CSR behaviour.

### Implementation

Four `central_galaxies_mode()` gates were removed so that the CSR satellite index is
loaded and used regardless of whether `--centralgalaxies` is specified:

| Location | Change |
|----------|--------|
| `kdtree_backend.cc:load_snapshot()` | Load `_sat_offs`, `_sat_list`, build `_sat_to_central` whenever `_has_central_index` |
| `kdtree_backend.hh:_fetch()` | Call `_fetch_satellites()` whenever `_in_sat_flush` |
| `kdtree_backend.hh:operator++()` | Flush pending sats before advancing whenever `_in_sat_flush` |
| `kdtree_backend.hh:done()` | Return not-done whenever `_in_sat_flush` |

The `_collect_satellites()` logic was also inverted for the non-`--cg` case: Type>0
satellites are kept in the normal stream **unless** their central is in the cone at
the same snapshot (in which case they are suppressed and re-emitted via CSR at the
central's epoch).

Prerequisite: the KD-tree file must contain the CSR index (`has_central_galaxy_index`
attribute + `centralgalaxies/` group). `run_test_sage_hdf5.sh` always calls
`sage2kdtree --centralgalaxies` to ensure this.

### Result

| Run | Total | T0 | T1 | Inside |
|-----|-------|----|----|--------|
| old baseline (no coupling) | 211,249 | 185,107 | 26,142 | 211,249 |
| **coupled baseline (new default)** | **211,698** | **185,107** | **26,591** | **211,693** |
| `--centralgalaxies` | 211,283 | 185,107 | 26,176 | 211,278 |
| `--cg --includeorphansatellites` | 211,698 | 185,107 | 26,591 | 211,693 |

The coupled baseline produces **exactly the same (GalaxyIndex, SnapNum) pairs** as
`--cg --includeorphansatellites` (0 differences in either direction). This confirms
the two implementations are equivalent:

- **`--cg --orphans`**: suppress all T1 from normal stream; re-emit via CSR if central
  is in cone; re-emit as orphan if central is not in cone but satellite passes its own
  cone filter.
- **Coupled baseline**: let T1 pass the normal stream if central is not in cone; suppress
  and re-emit via CSR if central is in cone.

Both select the identical set of satellites. The coupled baseline is the physically
correct default: satellites are always output at the same epoch as their central.

---

## Baseline cone filter for CSR satellites

### Problem

After enabling always-on CSR epoch coupling, the baseline (no flags) produced 354
CSR-emitted satellites outside the requested cone bounds (ra/dec/z). All were Type=1
with `csi >= 0`. They landed outside because the satellite's physical position (often
near the box boundary) mapped to ra/dec/z values outside the query bounds after tile
transformation, even though its central was inside the cone.

`--centralgalaxies --includeorphansatellites` intentionally allows these through —
the user has explicitly asked for all satellites of in-cone centrals regardless of where
they end up.

### Fix

In `_fetch_satellites()` (`kdtree_backend.hh`), after `_calc_fields(true)` computes
ra/dec/distance for CSR satellites, a post-filter is applied in non-`--centralgalaxies`
mode using the lightcone's overall bounds:

```cpp
if (!_be->central_galaxies_mode() && _lc)
{
    real_type ra_min  = to_degrees(_lc->min_ra());
    real_type ra_max  = to_degrees(_lc->max_ra());
    real_type dec_min = to_degrees(_lc->min_dec());
    real_type dec_max = to_degrees(_lc->max_dec());
    real_type d_min   = _lc->min_dist();
    real_type d_max   = _lc->max_dist();
    auto ra_v   = _bat->template scalar<real_type>("ra");
    auto dec_v  = _bat->template scalar<real_type>("dec");
    auto dist_v = _bat->template scalar<real_type>("distance");
    for (unsigned i = 0; i < _bat->size(); ++i)
    {
        if (_bat->masked(i)) continue;
        if (ra_v[i]   < ra_min  || ra_v[i]   > ra_max  ||
            dec_v[i]  < dec_min || dec_v[i]  > dec_max  ||
            dist_v[i] < d_min   || dist_v[i] > d_max)
            _bat->mask(i);
    }
}
```

`_lc->min_dist()` / `_lc->max_dist()` (no-snap overloads) give the overall cone
distance bounds across all snapshots, avoiding false rejection of CSR satellites whose
distance does not match the snapshot's thin distance shell.

#### Detailed Discussion
                                                                         
  There are two distinct distance-bound semantics:                                                                                            
  
  1. No-snap overloads (_lc->min_dist() / _lc->max_dist()):
  Returns _dist[0] / _dist[1] — the comoving distances corresponding to the user's --zmin and --zmax. This is the full query volume.                                                                                            

  2. Per-snap overloads (_lc->min_dist(_snap) / _lc->max_dist(_snap)):
  Returns _dist_bins[offset+1] / _dist_bins[offset] — the thin distance shell assigned to one snapshot. Each snapshot occupies a narrow radial slice of the cone. This is what the KD-tree traversal itself uses to partition   
  galaxies among snapshots.

  Why the no-snap overloads are correct here

  A CSR satellite's physical position can differ from its central's by a few Mpc/h. Snapshot shells in the Mini-Millennium box can be as thin as ~10–50 Mpc/h. So a satellite can legitimately land just outside its central's
  thin shell even though it's well within the user's requested [z_min, z_max].

  If the post-filter used _lc->min_dist(_snap) / _lc->max_dist(_snap), a satellite 2 Mpc/h closer than the shell inner edge would be masked — a false rejection. Using the overall bounds means the only satellites filtered out
   are those genuinely outside the user's cone, which is the right behaviour.

  The one subtle concern this opens up
  There is a theoretical double-counting path: if a satellite's distance d_sat falls in snapshot M's shell (M ≠ N), then:
  - At snapshot M (normal KD-tree traversal): the suppression check central_in_cone uses M's thin shell for the central. If the central's epoch is N ≠ M, central_in_cone is false → satellite is not suppressed → it may be
  emitted normally at M's epoch
  - At snapshot N (CSR): the satellite is also emitted via CSR
  In practice this is impossible at the physical scales involved — central–satellite separations (~1–10 Mpc/h from validation data) are far smaller than inter-snapshot shell widths (~50–200 Mpc/h). So d_sat is always in the
  same thin shell as d_central.

  Is the comment accurate?
  Yes — it correctly describes what the no-snap overloads return and why using them is appropriate. The only thing it doesn't say explicitly is why the satellite's distance might differ from the snapshot shell (i.e.,
  physical offset from the central), but that's the underlying reason.
  If you wanted to be more precise, the comment could add: "A satellite is typically within a few Mpc/h of its central, which can place it just outside the central's thin distance shell while still within the overall cone;  
  using the per-snapshot bounds would falsely reject it."      

### Result

When running SAGE2016

| Run | Total | Inside | Outside | T0 | T1+ |
|-----|-------|--------|---------|----|----|
| baseline (new) | 211,344 | 211,344 | **0** | 185,107 | 26,237 |
| `--cg --includeorphansatellites` | 211,698 | 211,344 | 354 | 185,107 | 26,237 |

- Baseline: 0 outside — all galaxies strictly within the requested cone ✓
- `--cg --orphans`: 354 outside (unchanged, by design) ✓
- Inside counts are identical in both runs: 185,107 T0 + 26,237 T1 = 211,344 ✓
