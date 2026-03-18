# Field Naming Conventions

This document describes the field naming conventions used throughout the TAO Lightcone CLI pipeline, specifically for the `sage2kdtree` → `cli_lightcone` workflow.

## Overview

The pipeline preserves SAGE's original field names throughout the conversion process, ensuring consistency from input to output. Field names use two distinct conventions depending on their source.

## Naming Conventions

### SAGE Fields (CamelCase)

Fields that originate from SAGE HDF5 output use **CamelCase** naming:

**Position Fields** (critical for KD-tree spatial indexing):
- `Posx`, `Posy`, `Posz` - Galaxy position coordinates

**Velocity Fields** (critical for redshift calculations):
- `Velx`, `Vely`, `Velz` - Galaxy velocity components

**Mass Fields**:
- `StellarMass` - Total stellar mass
- `ColdGas` - Cold gas mass
- `HotGas` - Hot gas mass
- `BulgeMass` - Bulge mass
- `BlackHoleMass` - Black hole mass
- `EjectedMass` - Ejected gas mass

**Star Formation Fields**:
- `SfrDisk` - Star formation rate in disk
- `SfrBulge` - Star formation rate in bulge
- `SfrDiskZ` - Metallicity-weighted SFR in disk
- `SfrBulgeZ` - Metallicity-weighted SFR in bulge

**Metallicity Fields**:
- `MetalsColdGas` - Metals in cold gas
- `MetalsStellarMass` - Metals in stars
- `MetalsHotGas` - Metals in hot gas
- `MetalsBulgeMass` - Metals in bulge
- `MetalsEjectedMass` - Metals in ejected gas

**Halo Properties**:
- `Mvir` - Virial mass
- `Rvir` - Virial radius
- `Vvir` - Virial velocity
- `Vmax` - Maximum circular velocity

**Mixed Case Fields**:
- `infallMvir` - Virial mass at infall
- `infallVvir` - Virial velocity at infall
- `infallVmax` - Maximum velocity at infall

**Tree Structure Fields**:
- `SnapNum` - Snapshot number
- `GalaxyIndex` - Unique galaxy identifier within tree
- `CentralGalaxyIndex` - Index of central galaxy
- `SAGETreeIndex` - SAGE tree identifier (passes through; not mandatory)
- `mergeType` - Type of merger (0=none, 1=minor, 2=major)
- `mergeIntoID` - Galaxy index this galaxy merges into
- `mergeIntoSnapNum` - Snapshot number of merger

Note: `mergeType`, `mergeIntoID`, and `mergeIntoSnapNum` are present in SAGE HDF5 output and pass through the pipeline unchanged, but are **not mandatory** — the pipeline does not read them.

### Pipeline-Computed Fields (lowercase)

Fields created by the `sage2kdtree` pipeline use **lowercase** naming:

**Tree Traversal Fields**:
- `globaltreeid` - Global tree identifier across all files
- `breadthfirst_traversalorder` - Breadth-first traversal index
- `depthfirst_traversalorder` - Depth-first traversal index
- `localgalaxyid` - Local galaxy ID within tree
- `subtree_count` - Number of galaxies in subtree

**Index Fields**:
- `local_index` - Local index within snapshot
- `global_index` - Global index across all galaxies
- `descendant` - Local index of descendant galaxy
- `global_descendant` - Global index of descendant galaxy
- `subsize` - Size of subtree rooted at this galaxy

**Backward Compatibility Alias**:
- `snapnum` - Lowercase alias for `SnapNum` (supported for queries)

**Lightcone-Derived Fields** (computed during extraction):
- `distance` - Comoving distance
- `ra` - Right ascension
- `dec` - Declination
- `redshift_cosmological` - Cosmological redshift
- `redshift_observed` - Observed redshift (includes peculiar velocity)
- `sfr` - Total star formation rate (SfrDisk + SfrBulge)

## Mandatory Fields

The following SAGE fields **must** be present in the input HDF5 files for the pipeline to function:

### Critical for KD-tree Operation:
```
Posx, Posy, Posz       # Spatial coordinates
SnapNum                 # Snapshot number
```

### Required for Tree Structure:
```
GalaxyIndex            # Galaxy identifier
CentralGalaxyIndex     # Central galaxy reference
```

### Required for Lightcone Extraction:
```
Velx, Vely, Velz       # Velocity components for redshift calculation
```

These mandatory fields are validated at pipeline entry points. If any are missing, the pipeline will abort with a clear error message listing the missing fields.

## Case-Insensitive Field Lookup

User-facing tools (e.g., `plot_lightcone.py`) support **case-insensitive** field matching for flexibility:

```bash
# All of these are equivalent when querying:
python plot_lightcone.py output.h5 StellarMass
python plot_lightcone.py output.h5 stellarmass
python plot_lightcone.py output.h5 stellar_mass
```

The tool will find the field regardless of case, making scripts more robust to field name variations.

## HDF5 File Structure

### SAGE HDF5 Output (Input to `sage2kdtree`)
```
model_X.hdf5
├── Header/                    # Cosmology and simulation parameters
├── Snap_0/                    # Columnar datasets per snapshot
│   ├── Posx, Posy, Posz      # CamelCase field names
│   ├── StellarMass
│   ├── ColdGas
│   └── ... (50+ SAGE fields)
├── Snap_1/
└── TreeInfo/
```

### KD-tree Indexed HDF5 (Output of `sage2kdtree`, Input to `cli_lightcone`)
```
kdtree-output.h5
├── cosmology/                 # Cosmology parameters
├── data/                      # All fields spatially reordered by snapshot
│   ├── Posx, Posy, Posz      # SAGE CamelCase fields preserved
│   ├── StellarMass
│   ├── ColdGas
│   ├── globaltreeid          # Computed lowercase fields
│   ├── snapnum               # Lowercase alias
│   └── ... (SAGE + computed fields)
├── lightcone/
│   ├── snapshot000/          # KD-tree index per snapshot
│   │   ├── bounds
│   │   ├── cell_counts
│   │   ├── cell_offs
│   │   └── splits
│   └── ...
├── snapshot_counts
├── snapshot_displs
└── snapshot_redshifts
```

## Validation

Field validation occurs at two points:

1. **`sage2kdtree` (Phase 1 - Input Validation)**:
   - Validates all mandatory SAGE fields are present in input HDF5
   - Aborts with clear error message if any are missing
   - See `src/libtao/base/mandatory_fields.hh` for the complete list

2. **`kdtree_backend::open()` (KD-tree File Validation)**:
   - Validates KD-tree file has required spatial fields (Posx, Posy, Posz, SnapNum)
   - Uses case-insensitive matching for flexibility
   - Throws error if fields are missing with suggestion to regenerate file

## Migration Notes

**No backward compatibility** for existing KD-tree files with old (snake_case) field names.

If you have existing KD-tree files created before this field naming standardization:
1. Delete the old KD-tree files
2. Regenerate from original SAGE HDF5 input using `sage2kdtree`

Example:
```bash
cd tests/sage-model-tests
rm -rf output_sage_hdf5_one_step/*
./run_test_hdf5_one_step.sh
```

## Implementation References

- **Field name constants**: `src/libtao/base/mandatory_fields.hh`
- **Validation logic**: `src/libtao/base/mandatory_fields.cc`
- **Phase 1 validation**: `src/apps/sage2kdtree.cc` (lines 1013-1034)
- **KD-tree validation**: `src/libtao/base/kdtree_backend/kdtree_backend.cc` (lines 208-234)
- **Case-insensitive plotting**: `src/plot_lightcone.py` (function `find_field_case_insensitive`)

## Best Practices

1. **Use canonical SAGE names** in scripts and documentation (e.g., `StellarMass`, not `stellar_mass`)
2. **Let tools handle case-insensitivity** - don't manually normalize field names
3. **Check mandatory_fields.hh** before modifying required field lists
4. **Validate early** - field validation at pipeline entry prevents late failures
5. **Document new computed fields** - add them to `computed_fields` namespace in `mandatory_fields.hh`
