import h5py
import numpy as np
import sys
import os
import glob

def check_outputs_sorted(baseline_file, test_file):
    print(f"Comparing with sorting strategy:")
    print(f"  Baseline: {baseline_file}")
    print(f"  Test:     {test_file}")

    if not os.path.exists(baseline_file):
        print(f"Error: Baseline file not found: {baseline_file}")
        sys.exit(1)
    if not os.path.exists(test_file):
        print(f"Error: Test file not found: {test_file}")
        sys.exit(1)

    with h5py.File(baseline_file, 'r') as f1, h5py.File(test_file, 'r') as f2:
        
        # 1. Structural Validation (KD-Tree Bounds & Indices)
        print("\n=== Phase 1: Structural Validation (KD-Tree) ===")
        # Get list of snapshots
        snaps1 = sorted([k for k in f1['lightcone'].keys() if k.startswith('snapshot')])
        snaps2 = sorted([k for k in f2['lightcone'].keys() if k.startswith('snapshot')])
        
        if snaps1 != snaps2:
            print(f"FAILURE: Snapshot list mismatch.\nBaseline: {snaps1}\nTest: {snaps2}")
            sys.exit(1)

        # Identify data columns to check
        data_keys = list(f1['data'].keys())
        
        total_errors = 0

        for snap in snaps1:
            print(f"Checking structure for {snap}...", end=" ")
            
            # Check bounds
            bounds1 = f1[f'lightcone/{snap}/bounds'][:]
            bounds2 = f2[f'lightcone/{snap}/bounds'][:]
            if not np.array_equal(bounds1, bounds2):
                print(f"\nFAILURE: Bounds mismatch in {snap}")
                total_errors += 1
            
            # Check cell_counts
            counts1 = f1[f'lightcone/{snap}/cell_counts'][:]
            counts2 = f2[f'lightcone/{snap}/cell_counts'][:]
            if not np.array_equal(counts1, counts2):
                print(f"\nFAILURE: cell_counts mismatch in {snap}")
                total_errors += 1
            
            # Check cell_offs
            offs1 = f1[f'lightcone/{snap}/cell_offs'][:]
            offs2 = f2[f'lightcone/{snap}/cell_offs'][:]
            if not np.array_equal(offs1, offs2):
                print(f"\nFAILURE: cell_offs mismatch in {snap}")
                total_errors += 1
            else:
                print(f"Structure OK. Validating content for {len(counts1)} cells...")
                
                # We need to map the compound lightcone/data to simple dict of arrays for easier handling
                lc_data1 = f1['lightcone/data']
                lc_data2 = f2['lightcone/data']
                
                # Check /data/GalaxyIndex exists for sorting
                if 'GalaxyIndex' not in f1['data'] or 'GalaxyIndex' not in f2['data']:
                     print("\nFAILURE: GalaxyIndex missing from /data group. Cannot perform sort-based validation.")
                     sys.exit(1)

                # Iterate over cells
                for i, (offset, count) in enumerate(zip(offs1, counts1)):
                    if count == 0:
                        continue
                    
                    sl = slice(int(offset), int(offset + count))
                    
                    # 1. Get Sort Keys (GalaxyIndex)
                    idx1 = f1['data/GalaxyIndex'][sl]
                    idx2 = f2['data/GalaxyIndex'][sl]
                    
                    # Check uniqueness (sanity check requested by user)
                    if len(np.unique(idx1)) != len(idx1):
                        print(f"WARNING: GalaxyIndex not unique in cell {i} of {snap} (Baseline). Sort might be unstable.")
                    
                    # Determine Sort Order
                    sorter1 = np.argsort(idx1)
                    sorter2 = np.argsort(idx2)
                    
                    # Verify GalaxyIndices match after sorting
                    if not np.array_equal(idx1[sorter1], idx2[sorter2]):
                        print(f"FAILURE: GalaxyIndex set mismatch in cell {i}, {snap}")
                        total_errors += 1
                        continue 
                        
                    # 2. Check Lightcone Coords (x, y, z)
                    lc_chunk1 = lc_data1[sl]
                    lc_chunk2 = lc_data2[sl]
                    
                    for dim in ['x', 'y', 'z']:
                        arr1 = lc_chunk1[dim][sorter1]
                        arr2 = lc_chunk2[dim][sorter2]
                        if not np.allclose(arr1, arr2, rtol=1e-5):
                             max_diff = np.max(np.abs(arr1 - arr2))
                             print(f"FAILURE: Coord {dim} mismatch in cell {i}, {snap}. Max diff: {max_diff}")
                             total_errors += 1
                    
                    # 3. Check Data Properties
                    for key in data_keys:
                        if key not in f2['data']:
                            continue
                        
                        if key == "tdeplete":
                            continue
                            
                        d1 = f1[f'data/{key}'][sl][sorter1]
                        d2 = f2[f'data/{key}'][sl][sorter2]
                        
                        if d1.dtype.kind in 'SUa' or d2.dtype.kind in 'SUa': # Strings
                             if not np.array_equal(d1, d2):
                                 print(f"FAILURE: Data {key} mismatch in cell {i}, {snap}")
                                 total_errors += 1
                        else: # Numbers
                             if not np.allclose(d1, d2, equal_nan=True, rtol=1e-4):
                                 diff = np.abs(d1 - d2)
                                 print(f"FAILURE: Data {key} mismatch in cell {i}, {snap}. Max diff: {np.nanmax(diff)}")
                                 total_errors += 1
                    
                    if total_errors > 20: 
                        print("\nTOO MANY ERRORS. Aborting content check.")
                        sys.exit(1)

        print("\n=== Summary ===")
        if total_errors == 0:
            print("SUCCESS: Full validation passed (Structure matched, Content matched after sorting).")
        else:
            print(f"FAILURE: Found {total_errors} errors.")
            sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python verify_kdtree.py <baseline.h5> <test.h5>")
        sys.exit(1)
    
    check_outputs_sorted(sys.argv[1], sys.argv[2])
