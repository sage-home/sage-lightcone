# Phase 2 Lustre Performance Analysis

## Problem Summary

`phase2_add_traversal_metadata()` is significantly slower on Lustre (HPC) compared to local filesystems (MacOS) due to I/O patterns that are poorly suited for parallel distributed filesystems.

## Root Causes

### 1. Tree-by-Tree I/O Pattern (src/apps/sage2kdtree.cc:1911-2003)

**Current Implementation:**
- Processes each tree individually in a loop
- For each tree: 1 read + 5 writes = 6 I/O operations
- For N trees: 6N total I/O operations with hyperslab selections

**Lustre Impact:**
- Each hyperslab + I/O operation involves metadata operations (~1-10ms latency)
- For 10,000 trees: 60,000 operations × 5ms = **5 minutes of metadata overhead alone**
- Lock contention when multiple MPI ranks access the same file
- Defeats Lustre's large sequential I/O optimization

### 2. No MPI-IO Hints

**Current Code (src/libhpc/h5/property_list.cc:92):**
```cpp
H5Pset_fapl_mpio(_id, comm.mpi_comm(), MPI_INFO_NULL);
```

**Missing Optimizations:**
- `striping_factor`: Number of OSTs to stripe across
- `striping_unit`: Size of each stripe
- `cb_nodes`: Number of collective buffering aggregators
- `romio_cb_write/read`: Enable collective buffering
- `romio_ds_write/read`: Enable data sieving

### 3. No HDF5 Dataset Chunking

**Current Code (line 1811-1813):**
```cpp
out_datasets[field.name] = std::make_shared<hpc::h5::dataset>(
    out_file, "fields/" + field.name, field.type, hpc::h5::dataspace(total_gals)
);
```

**Issue:** Contiguous datasets without chunking perform poorly for non-contiguous access patterns.

### 4. Independent I/O Operations

**Current Code (lines 1998-2002):**
```cpp
out_datasets["globaltreeid"]->write(globaltreeid.data(), ...);
out_datasets["breadthfirst_traversalorder"]->write(bfs_order.data(), ...);
// ... more writes
```

**Issue:** No explicit collective I/O hints, defaults to independent I/O which can't leverage:
- Collective buffering
- Aggregation across ranks
- Coordinated I/O scheduling

## Performance Comparison: MacOS vs Lustre

| Characteristic | MacOS (Local SSD) | Lustre (HPC) |
|----------------|-------------------|--------------|
| Metadata latency | <1ms | 1-10ms |
| Small I/O overhead | Minimal | Very high |
| Network latency | None | 1-5ms |
| Caching | Aggressive | Conservative |
| Optimal access | Any | Large sequential |

## Recommended Solutions (Ordered by Impact)

### Priority 1: Batch Tree I/O Operations

**Instead of:** Read/write each tree individually
**Do:** Batch multiple trees into larger I/O operations

**Implementation Strategy:**
```cpp
// Read all descendant data in larger chunks (e.g., 100 trees at a time)
const size_t TREE_BATCH_SIZE = 100;
for (size_t batch_start = start_tree; batch_start < end_tree; batch_start += TREE_BATCH_SIZE) {
    size_t batch_end = std::min(batch_start + TREE_BATCH_SIZE, end_tree);

    // Calculate total size for this batch
    long long batch_size = tree_offsets[batch_end] - tree_offsets[batch_start];
    long long batch_offset = tree_offsets[batch_start];

    // Read entire batch in one operation
    std::vector<int> descendant_batch(batch_size);
    hpc::h5::dataspace mem_space(batch_size);
    hpc::h5::dataspace file_space(total_gals);
    file_space.select_hyperslab(H5S_SELECT_SET, batch_size, batch_offset);
    get_dataset("descendant").read(descendant_batch.data(), ...);

    // Process trees within batch
    // Write results in one operation per field (not per tree)
}
```

**Expected Impact:** 10-100x reduction in I/O operations

### Priority 2: Add MPI-IO Hints for Lustre

**File:** `src/libhpc/h5/property_list.cc`

**Add before H5Pset_fapl_mpio:**
```cpp
void property_list::set_parallel(mpi::comm const &comm) {
    MPI_Info info;
    MPI_Info_create(&info);

    // Lustre optimizations
    MPI_Info_set(info, "romio_cb_write", "enable");        // Collective buffering for writes
    MPI_Info_set(info, "romio_cb_read", "enable");         // Collective buffering for reads
    MPI_Info_set(info, "romio_ds_write", "enable");        // Data sieving for writes
    MPI_Info_set(info, "romio_ds_read", "enable");         // Data sieving for reads

    // Set aggregators based on rank count
    int cb_nodes = std::min(comm.size() / 4, 16);  // 1 aggregator per 4 ranks, max 16
    char cb_nodes_str[32];
    sprintf(cb_nodes_str, "%d", cb_nodes);
    MPI_Info_set(info, "cb_nodes", cb_nodes_str);

    // Collective buffer size (8MB per aggregator)
    MPI_Info_set(info, "cb_buffer_size", "8388608");

    // Optional: Set striping if known
    // MPI_Info_set(info, "striping_factor", "4");
    // MPI_Info_set(info, "striping_unit", "1048576");  // 1MB

    INSIST(H5Pset_fapl_mpio(_id, comm.mpi_comm(), info), >= 0);
    MPI_Info_free(&info);
}
```

**Expected Impact:** 2-5x improvement for collective operations

### Priority 3: Enable HDF5 Dataset Chunking

**Add chunking to dataset creation (line 1811):**
```cpp
hpc::h5::property_list props;
props.create(H5P_DATASET_CREATE);

// Set chunk size (tune based on typical tree size)
hsize_t chunk_size = 100000;  // 100K galaxies per chunk
props.set_chunk_size(chunk_size);

out_datasets[field.name] = std::make_shared<hpc::h5::dataset>(
    out_file, "fields/" + field.name, field.type,
    hpc::h5::dataspace(total_gals), props
);
```

**Expected Impact:** 1.5-3x improvement for non-contiguous access

### Priority 4: Use Collective I/O Calls

**Modify write operations to be explicitly collective:**
```cpp
// Set collective I/O mode
hpc::h5::property_list xfer_props;
xfer_props.create(H5P_DATASET_XFER);
H5Pset_dxpl_mpio(xfer_props.id(), H5FD_MPIO_COLLECTIVE);

// All ranks participate in collective write (even if writing 0 bytes)
out_datasets["globaltreeid"]->write(globaltreeid.data(), dtype,
                                     mem_space, file_space, xfer_props);
```

**Expected Impact:** 1.5-2x improvement when combined with batching

## Diagnostic Commands

### Check Lustre Striping on HPC
```bash
# Check file striping
lfs getstripe /path/to/your/output.h5

# Recommended striping for large files (>1GB)
lfs setstripe -c 4 -S 1M /path/to/output_directory/
```

### Monitor I/O Performance
```bash
# Use Darshan for I/O profiling
module load darshan
export LD_PRELOAD=$DARSHAN_LIB
# Run your application

# Analyze Darshan log
darshan-parser your-log.darshan
```

### HDF5 Performance Metrics
Add to your code:
```cpp
// Before phase 2
double t_start = MPI_Wtime();

// After phase 2
double t_end = MPI_Wtime();
if (rank == 0) {
    printf("Phase 2 time: %.2f seconds\n", t_end - t_start);
    printf("Number of trees: %zu\n", n_trees);
    printf("Time per tree: %.2f ms\n", (t_end - t_start) * 1000.0 / n_trees);
}
```

## Implementation Priority

1. **Quick Win (1-2 hours):** Add MPI-IO hints (Priority 2)
2. **Medium Effort (4-6 hours):** Implement tree batching (Priority 1)
3. **Polish (2-3 hours):** Add chunking + collective I/O (Priorities 3 & 4)

## Expected Overall Improvement

With all optimizations: **20-100x speedup on Lustre** depending on:
- Number of trees
- Tree size distribution
- Number of MPI ranks
- Lustre configuration

## References

- [HDF5 Parallel I/O Best Practices](https://portal.hdfgroup.org/display/HDF5/Parallel+HDF5)
- [ROMIO MPI-IO Hints](https://www.mcs.anl.gov/research/projects/romio/hints/)
- [Lustre I/O Best Practices](https://wiki.lustre.org/Optimizing_Lustre_IO)
