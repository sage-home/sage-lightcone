# Computed Fields in sage2kdtree Pipeline

This document lists fields that are computed by the `sage2kdtree` pipeline and do not exist in the original SAGE HDF5 output.

---

## Do we need to keep these?

### In `data/` group (columnar):

| Field | Description | Added in Phase | Required? |
|-------|-------------|----------------|-----------|
| `breadthfirst_traversalorder` | BFS traversal index within tree | Phase 2 | ? |
| `depthfirst_traversalorder` | DFS traversal index within tree | Phase 2 | ? |
| `descendant` | Local index of descendant galaxy | Phase 1 | ? |
| `global_descendant` | Global index of descendant galaxy | Phase 1 | ? |
| `global_index` | Global index across all galaxies | Phase 1 | ? |
| `globaltreeid` | Global tree identifier | Phase 2 | ? |
| `local_index` | Local index within tree | Phase 1 | ? |
| `localgalaxyid` | Local galaxy ID within tree | Phase 2 | ? |
| `subtree_count` | Number of descendants in subtree | Phase 2 | Yes (used for `subsize`) |

### In `lightcone/data` compound only:

| Field | Description | Added in Phase | Required? |
|-------|-------------|----------------|-----------|
| `subsize` | Subtree size (copied from `subtree_count`) | Phase 4 | Yes (KD-tree queries) |

---

## Notes

- `subsize` in `lightcone/data` is required for KD-tree spatial queries
- `subtree_count` provides the source values for `subsize`
- The other traversal fields may be useful for debugging/analysis but may not be required for lightcone extraction
- Consider whether `cli_lightcone` or downstream tools depend on any of these fields

---

## Comparison with OLD workflow

The OLD workflow (`sageh5toh5 → sageimport → dstreeinit`) produces the same computed fields plus:
- `data/subsize` (separate columnar field, int32)

The NEW workflow (`sage2kdtree`) uses `data/subtree_count` (int64) instead, with identical values.


# Computed fields for cli_lightcone

Computed Fields Created by cli_lightcone                                                                                                                                                                       
  ┌───────────────────────┬─────────────────────────────────────────────┬───────────────────────────────────┐
  │         Field         │                 Computation                 │         Saved to Output?          │
  ├───────────────────────┼─────────────────────────────────────────────┼───────────────────────────────────┤
  │ distance              │ sqrt(posx² + posy² + posz²)                 │ Yes, if requested via --outfields │
  ├───────────────────────┼─────────────────────────────────────────────┼───────────────────────────────────┤
  │ ra                    │ Cartesian → spherical conversion (degrees)  │ Yes, if requested                 │
  ├───────────────────────┼─────────────────────────────────────────────┼───────────────────────────────────┤
  │ dec                   │ Cartesian → spherical conversion (degrees)  │ Yes, if requested                 │
  ├───────────────────────┼─────────────────────────────────────────────┼───────────────────────────────────┤
  │ cosmological_redshift │ Distance → redshift via cosmology           │ Yes, if requested                 │
  ├───────────────────────┼─────────────────────────────────────────────┼───────────────────────────────────┤                                    
  │ observed_redshift     │ Cosmological + peculiar velocity correction │ Yes, if requested                 │                       
  ├───────────────────────┼─────────────────────────────────────────────┼───────────────────────────────────┤                                    
  │ sfr                   │ sfrdisk + sfrbulge                          │ Yes, if requested                 │
  └───────────────────────┴─────────────────────────────────────────────┴───────────────────────────────────┘                                                                                                                                                                                                     
  Output Behavior                                                                                                                                                                                                
                                                                                                                                                                                                                 
  The output HDF5 lightcone file contains:                                                                                                                                                                       
  1. If --outfields specified: Only the requested fields (can include computed fields)                                                                                                                           
  2. If no --outfields: All fields from the input KD-tree /data/* group                                                                                                                                          
                                                                                                                                                                                                                 
  Key implementation locations:                                                                                                                                                                                  
  - Computed fields: src/libtao/base/kdtree_backend/kdtree_backend.hh:600-761                                                                                                                                    
  - Base output fields: src/libtao/base/query.hh:70-90                                                                                                                                                           
  - HDF5 output: src/libtao/modules/hdf5.hh:400-665                                                                                                                                                              
                                                                                                                                                                                                                 
  So unlike sage2kdtree, cli_lightcone creates new computed fields (ra, dec, redshifts, distance, sfr) that are useful for astronomical analysis and are written to output if requested.   
