# Plan to Remove Compound HDF5 Datatypes

## Executive Summary

**Goal**: Replace compound HDF5 datatypes with columnar storage (separate datasets per field) throughout the entire pipeline.

**Alignment with VISION.md**:
- ✅ **Principle 1**: SAGE HDF5 already uses columnar storage - we should match it
- ✅ **Principle 3**: Metadata-driven architecture - dynamic field discovery
- ✅ **Principle 4**: Single source of truth - SAGE HDF5 structure is canonical
- ✅ **Principle 6**: Memory efficiency - read only fields you need

**Impact**: High value, medium effort refactoring that modernizes the codebase.

---

## Current State Analysis

### Where Compound Types Are Used

#### 1. **Intermediate Files (Phases 1-3)**
**Files affected:**
- `*-depthfirstordered.h5` (Phase 1 output)
- `out_*-depthfirstordered.h5` (Phase 2 output)
- `*-bysnap.h5` (Phase 3 output)

**Structure:**
```
depthfirstordered.h5/
  ├── galaxies              [COMPOUND DATASET - 60+ fields]
  ├── tree_displs           [separate dataset]
  └── tree_counts           [separate dataset]
```

#### 2. **Code Using Compound Types**

| File | Function | Purpose |
|------|----------|---------|
| `src/libtao/base/sage.cc` | `make_hdf5_types()` | Defines compound type for galaxy struct |
| `src/apps/sageh5toh5.cc` | Phase 1 | Writes compound dataset |
| `src/apps/sage2kdtree.cc` | Phase 1 | Writes compound dataset |
| `src/apps/sageimport.cc` | Phase 2 | Reads/writes compound dataset |
| `src/apps/dstreeinit.cc` | Phase 3→4 | Reads compound, writes separate |
| `src/apps/sage2kdtree.cc` | Phase 3→4 | Reads compound, writes separate |

#### 3. **What Already Works Correctly (Columnar)**

✅ **SAGE HDF5 Input** (already columnar):
```
model_0.hdf5/
  ├── Snap_0/
  │   ├── Posx              [separate dataset]
  │   ├── Posy              [separate dataset]
  │   ├── StellarMass       [separate dataset]
  │   └── ... (50+ fields)
```

✅ **Final KD-tree Output** (already columnar):
```
kdtree.h5/
  ├── data/
  │   ├── posx              [separate dataset]
  │   ├── posy              [separate dataset]
  │   ├── stellar_mass      [separate dataset]
  │   └── ... (60+ fields)
```

### The Problem

**Current flow:**
```
SAGE (columnar)
  → Phase 1: Convert to COMPOUND
  → Phase 2: Keep COMPOUND
  → Phase 3: Keep COMPOUND
  → Phase 4: Convert back to COLUMNAR
```

**Issues:**
1. **Unnecessary conversions**: Columnar → Compound → Columnar
2. **Inflexibility**: Must read all 60+ fields even if you only need positions
3. **Maintenance burden**: Keep `sage::galaxy` struct in sync with fields
4. **Not metadata-driven**: Hardcoded struct definitions
5. **Memory inefficient**: Loading full records when only need subset
6. **Against HDF5 best practices**: Modern HDF5 encourages columnar storage

---

## Proposed Solution

### Target Architecture

**Simplified flow:**
```
SAGE (columnar)
  → Phase 1: Keep COLUMNAR
  → Phase 2: Keep COLUMNAR
  → Phase 3: Keep COLUMNAR
  → Phase 4: Keep COLUMNAR
```

### New File Structures

#### Phase 1 Output: Depth-First Ordered (Columnar)
```
depthfirstordered.h5/
  ├── cosmology/
  │   ├── box_size
  │   ├── hubble
  │   ├── omega_l
  │   └── omega_m
  ├── fields/              [NEW: separate dataset per field]
  │   ├── snapnum
  │   ├── type
  │   ├── galaxy_index
  │   ├── posx, posy, posz
  │   ├── stellar_mass
  │   └── ... (all 60+ fields)
  ├── tree_displs
  ├── tree_counts
  └── snapshot_redshifts
```

**Benefits:**
- Match SAGE HDF5 input structure
- Can read individual fields on demand
- Easy to add new fields without struct changes
- Metadata-driven field discovery

---

## Implementation Plan

### Phase A: Refactor Core Libraries (Foundation)

#### A.1: Remove `make_hdf5_types()` dependency

**File**: `src/libtao/base/sage.cc`

**Current approach:**
```cpp
void make_hdf5_types(h5::datatype& mem_type, h5::datatype& file_type) {
    h5::derive der(sizeof(galaxy));
    der.add(h5::datatype::native_int, HOFFSET(galaxy, snapshot),
            h5::datatype::native_int, "snapnum");
    // ... 60+ more fields
    der.commit(mem_type, file_type);
}
```

**New approach:**
```cpp
// Field metadata structure
struct FieldInfo {
    std::string name;
    hpc::h5::datatype type;
    size_t offset;
};

// Generate field list from SAGE HDF5 dynamically
std::vector<FieldInfo> discover_fields_from_sage_hdf5(hpc::h5::file& sage_file);

// Or use static list (simpler migration)
std::vector<FieldInfo> get_galaxy_field_list();
```

**Action Items:**
- [ ] Create `FieldInfo` struct in `sage.hh`
- [ ] Implement `get_galaxy_field_list()` with all 60+ fields
- [ ] Mark `make_hdf5_types()` as deprecated

**Effort**: Small (1-2 hours)
**Risk**: Low

---

### Phase B: Refactor Phase 1 (SAGE → Depth-First)

#### B.1: Update `sage2kdtree.cc` Phase 1

**Current**: Writes compound `galaxies` dataset
**New**: Writes separate dataset per field

**Changes needed:**

1. **Remove compound type creation**
```cpp
// DELETE: _make_output_hdf5_types()
// DELETE: _mem_type, _file_type member variables
```

2. **Create datasets per field**
```cpp
void phase1_sageh5_to_depthfirst() {
    // ... existing tree counting logic ...

    // Create separate dataset for each field
    auto fields = get_galaxy_field_list();
    std::map<std::string, hpc::h5::dataset> field_datasets;
    std::map<std::string, hpc::h5::buffer<T>> field_buffers;

    for (const auto& field : fields) {
        field_datasets[field.name] = hpc::h5::dataset(
            out_file, "fields/" + field.name, field.type,
            hpc::h5::dataspace(total_gals)
        );
        // Create buffer for this field
        field_buffers[field.name].create(
            field_datasets[field.name], field.type,
            hpc::h5::buffer_default_size, global_gal_offset
        );
    }

    // Process galaxies - write field by field
    for (auto& gal : galaxies_in_tree) {
        field_buffers["snapnum"].append(gal.snapshot);
        field_buffers["posx"].append(gal.pos[0]);
        field_buffers["stellar_mass"].append(gal.stellar_mass);
        // ... for each field
    }

    // Flush all buffers
    for (auto& [name, buffer] : field_buffers) {
        buffer.flush();
    }
}
```

**Effort**: Medium (4-6 hours)
**Risk**: Medium (need to test thoroughly)

#### B.2: Update `sageh5toh5.cc`

Similar changes to sage2kdtree Phase 1.

**Effort**: Medium (4-6 hours)
**Risk**: Medium

---

### Phase C: Refactor Phase 2 (Add Traversal Metadata)

#### C.1: Update `sageimport.cc`

**Current**: Reads compound dataset, adds BFS/DFS order fields, writes compound
**New**: Reads columnar, adds new fields, writes columnar

**Changes:**

1. **Read fields individually**
```cpp
// Instead of reading compound "galaxies" dataset:
auto snapshot_ds = in_file.dataset("fields/snapnum");
auto tree_idx_ds = in_file.dataset("fields/sage_tree_index");
// etc.

std::vector<int> snapnums(n_gals);
std::vector<int> tree_indices(n_gals);
snapshot_ds.read(snapnums);
tree_idx_ds.read(tree_indices);
```

2. **Write new traversal fields**
```cpp
// Create new fields for BFS/DFS order
hpc::h5::dataset bfs_ds(out_file, "fields/breadthfirst_traversalorder",
                        hpc::h5::datatype::native_llong,
                        hpc::h5::dataspace(n_gals));
hpc::h5::dataset dfs_ds(out_file, "fields/depthfirst_traversalorder",
                        hpc::h5::datatype::native_llong,
                        hpc::h5::dataspace(n_gals));

// Write computed orders
bfs_ds.write(bfs_order);
dfs_ds.write(dfs_order);
```

3. **Copy other fields through**
```cpp
// Copy all existing fields from input to output
auto field_names = get_all_field_names(in_file, "fields/");
for (const auto& name : field_names) {
    if (name != "breadthfirst_traversalorder" &&
        name != "depthfirst_traversalorder") {
        copy_dataset(in_file, out_file, "fields/" + name);
    }
}
```

**Effort**: Medium (6-8 hours) - sageimport is Python, needs careful handling
**Risk**: Medium

---

### Phase D: Refactor Phase 3 (Tree → Snapshot Order)

#### D.1: Update `sage2kdtree.cc` Phase 3

**Current**: Reads compound, reorders, writes compound per snapshot
**New**: Reads columnar, reorders, writes columnar per snapshot

**Key insight**: Reordering is **per-field**, not per-record!

```cpp
void phase3_tree_to_snapshot_order() {
    auto field_names = get_all_field_names(input_file, "fields/");

    // For each field, reorder independently
    for (const auto& field_name : field_names) {
        auto in_ds = input_file.dataset("fields/" + field_name);
        auto field_type = in_ds.datatype();

        // Read entire field (tree-ordered)
        std::vector<T> tree_ordered_data(total_gals);
        in_ds.read(tree_ordered_data);

        // Reorder by snapshot
        std::vector<T> snapshot_ordered_data(total_gals);
        for (size_t i = 0; i < total_gals; ++i) {
            int snap = snapshot_indices[i];  // from snapnum field
            size_t new_idx = snapshot_offsets[snap] + snapshot_positions[snap]++;
            snapshot_ordered_data[new_idx] = tree_ordered_data[i];
        }

        // Write snapshot-organized field
        // Option A: One dataset per snapshot per field (snapshot000/posx, snapshot000/posy, ...)
        // Option B: All snapshots in one field, indexed by snapshot_displs (simpler!)

        auto out_ds = out_file.dataset("fields/" + field_name, field_type,
                                       hpc::h5::dataspace(total_gals));
        out_ds.write(snapshot_ordered_data);
    }
}
```

**Effort**: Medium (4-6 hours)
**Risk**: Medium

---

### Phase E: Simplify Phase 4 (Already Mostly Columnar)

Phase 4 already writes columnar output under `data/`. Main changes:

1. **Read from columnar input** (instead of compound snapshots)
2. **Simplify field extraction** (no compound type introspection)

```cpp
void phase4_build_kdtree_index() {
    // Read columnar fields directly
    auto posx_ds = snap_file.dataset("fields/posx");
    auto posy_ds = snap_file.dataset("fields/posy");
    auto posz_ds = snap_file.dataset("fields/posz");

    std::vector<float> posx(n_gals), posy(n_gals), posz(n_gals);
    posx_ds.read(posx);
    posy_ds.read(posy);
    posz_ds.read(posz);

    // Build KD-tree from positions (existing logic)
    // ...

    // Copy all fields to output (already columnar)
    auto field_names = get_all_field_names(snap_file, "fields/");
    for (const auto& name : field_names) {
        copy_dataset(snap_file, out_file, "fields/" + name, "data/" + name);
    }
}
```

**Effort**: Small (2-3 hours)
**Risk**: Low

---

## Benefits Summary

### 1. **Flexibility** ⭐⭐⭐⭐⭐
- Read only fields you need (positions for visualization, stellar_mass for analysis)
- Add new fields without changing struct definitions
- Drop deprecated fields easily

### 2. **Performance** ⭐⭐⭐⭐
- Selective I/O: Read 3 fields instead of 60+ for spatial queries
- Better HDF5 chunking/compression per field
- Parallel I/O easier (different MPI ranks read different fields)

### 3. **Metadata-Driven** ⭐⭐⭐⭐⭐ (VISION.md Principle 3)
- Dynamic field discovery from SAGE HDF5
- No hardcoded struct definitions
- Automatically handle new SAGE fields

### 4. **Maintainability** ⭐⭐⭐⭐⭐
- No `sage::galaxy` struct to keep in sync
- No compound type definitions (60+ insert() calls)
- Simpler code - field operations explicit

### 5. **Memory Efficiency** ⭐⭐⭐⭐ (VISION.md Principle 6)
- Load only needed fields into memory
- For 1M galaxies x 60 fields x 4 bytes = 240MB
- Load only positions: 1M x 3 x 4 bytes = 12MB (20x reduction!)

### 6. **Consistency** ⭐⭐⭐⭐⭐ (VISION.md Principle 4)
- Match SAGE HDF5 input structure (single source of truth)
- Match final KD-tree output structure
- No format conversions in pipeline

---

## Implementation Strategy

### Recommended Sequence

**Iteration 1: Proof of Concept** (1 week)
1. Implement Phase B.1 (sage2kdtree Phase 1 columnar output)
2. Modify Phase 3 to accept columnar input
3. Test end-to-end with subset of fields (posx, posy, posz, stellar_mass)
4. Validate output matches current approach

**Iteration 2: Complete sage2kdtree** (1 week)
5. Extend to all 60+ fields
6. Update Phase 2 (traversal metadata)
7. Update Phase 4 (KD-tree building)
8. Full validation against baseline

**Iteration 3: Update sageh5toh5** (3-5 days)
9. Apply same changes to multi-step workflow
10. Update sageimport (Python)
11. Update dstreeinit

**Iteration 4: Polish & Documentation** (2-3 days)
12. Add field metadata (units, descriptions) as HDF5 attributes
13. Create field registry/catalog
14. Update documentation
15. Add tests for individual field I/O

---

## Risks & Mitigations

| Risk | Impact | Probability | Mitigation |
|------|--------|-------------|------------|
| Breaking existing workflows | High | Low | Keep compound code in separate branch initially |
| Performance regression | Medium | Low | Benchmark before/after, HDF5 chunking should be similar |
| Increased code complexity | Medium | Medium | Use helper functions for common operations |
| Field name inconsistencies | Medium | Medium | Create canonical field name mapping |
| Backward compatibility | High | High | NOT NEEDED - intermediate files are transient |

**Key insight**: Intermediate HDF5 files are **transient** (deleted after pipeline completes), so backward compatibility is not a concern!

---

## Alternative Approaches

### Option 1: Hybrid Approach ❌
Keep compound for intermediate files, columnar for final output.
**Rejected**: Still has conversion overhead, doesn't solve flexibility issues.

### Option 2: Keep Compound, Add Field Subsets ❌
Allow reading individual fields from compound dataset.
**Rejected**: HDF5 doesn't efficiently support this - must read full record.

### Option 3: Columnar with Index Dataset ✅ (RECOMMENDED)
Store fields separately + optional compound "index" for record-oriented access.
**Accepted**: Best of both worlds, but index only if needed.

---

## Success Metrics

1. ✅ **Correctness**: Output validates against current baseline (byte-for-byte for KD-tree)
2. ✅ **Performance**: No regression in end-to-end time (should improve)
3. ✅ **Memory**: Reduced peak memory usage (measure with profiling)
4. ✅ **Flexibility**: Can extract single field in < 1 second (currently needs full read)
5. ✅ **Code quality**: Reduced lines of code in type definitions (delete ~200 lines)

---

## Open Questions

1. **Field naming**: Keep SAGE naming (Posx) or normalize (posx)?
   - **Recommendation**: Normalize to lowercase for consistency

2. **Group structure**: `fields/*` or flat structure?
   - **Recommendation**: Use `fields/` group to separate from metadata

3. **Snapshot organization**: Separate groups (`Snap_000/posx`) or indices?
   - **Recommendation**: Use indices with `snapshot_displs` (simpler)

4. **Traversal fields**: Keep in same `fields/` or separate `metadata/`?
   - **Recommendation**: Keep in `fields/` (they're just more fields)

5. **Type preservation**: Keep SAGE types or normalize to native?
   - **Recommendation**: Native types (already decided)

---

## Alignment with Existing Plans

This plan directly implements **Priority 2** from `SAGE2KDTREE_IMPROVEMENTS.md`:

### Priority 2: Metadata-Driven Architecture ✅

**2.1 Dynamic Field Discovery** → Implemented by columnar approach
- Introspect SAGE HDF5 to find fields
- No hardcoded struct definitions

**2.2 Attribute Preservation** → Enhanced by columnar approach
- Each field dataset can have own attributes (units, description)
- Preserve SAGE attributes per-field

---

## Next Steps

1. **Review & Approve Plan** - Get team consensus
2. **Create Feature Branch** - `feature/columnar-storage`
3. **Start Iteration 1** - Proof of concept
4. **Validate Approach** - Compare outputs
5. **Iterate** - Complete remaining phases

---

## Conclusion

Removing compound datatypes is a **high-value refactoring** that:
- ✅ Aligns perfectly with VISION.md principles
- ✅ Modernizes codebase to HDF5 best practices
- ✅ Improves flexibility and performance
- ✅ Reduces maintenance burden
- ✅ Enables metadata-driven architecture

**Recommendation**: **PROCEED** with phased implementation starting with Iteration 1 proof of concept.
