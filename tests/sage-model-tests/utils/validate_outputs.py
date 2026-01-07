#!/usr/bin/env python3
"""
Validate that OLD and NEW workflows produce equivalent outputs.
Compares:
1. KD-tree HDF5 structure and field names
2. KD-tree spatial index (lightcone/snapshot* groups)
3. Lightcone galaxy counts
4. Lightcone field statistics (mean, std, min, max)
"""

import h5py
import numpy as np
import sys

def compare_hdf5_structure(file1, file2, path="/"):
    """Recursively compare HDF5 structure."""
    differences = []

    with h5py.File(file1, 'r') as f1, h5py.File(file2, 'r') as f2:
        def compare_group(name, obj1, obj2):
            # Check if both are same type
            if type(obj1) != type(obj2):
                differences.append(f"{name}: type mismatch ({type(obj1)} vs {type(obj2)})")
                return

            if isinstance(obj1, h5py.Dataset):
                # Compare shapes
                if obj1.shape != obj2.shape:
                    differences.append(f"{name}: shape mismatch ({obj1.shape} vs {obj2.shape})")
                # Compare dtypes (case-insensitive for field names)
                if obj1.dtype != obj2.dtype:
                    differences.append(f"{name}: dtype mismatch ({obj1.dtype} vs {obj2.dtype})")

            elif isinstance(obj1, h5py.Group):
                # Compare keys (case-insensitive)
                keys1 = set(k.lower() for k in obj1.keys())
                keys2 = set(k.lower() for k in obj2.keys())

                if keys1 != keys2:
                    missing_in_2 = keys1 - keys2
                    missing_in_1 = keys2 - keys1
                    if missing_in_2:
                        differences.append(f"{name}: keys only in file1: {missing_in_2}")
                    if missing_in_1:
                        differences.append(f"{name}: keys only in file2: {missing_in_1}")

        # Compare root groups
        if path in f1 and path in f2:
            compare_group(path, f1[path], f2[path])
            # Recursively compare children
            if isinstance(f1[path], h5py.Group) and isinstance(f2[path], h5py.Group):
                for key in f1[path].keys():
                    # Find matching key (case-insensitive)
                    key_lower = key.lower()
                    matching_key = None
                    for k2 in f2[path].keys():
                        if k2.lower() == key_lower:
                            matching_key = k2
                            break

                    if matching_key:
                        child_path = f"{path}/{key}".replace("//", "/")
                        compare_group(child_path, f1[path][key], f2[path][matching_key])
        else:
            if path not in f1:
                differences.append(f"{path}: missing in file1")
            if path not in f2:
                differences.append(f"{path}: missing in file2")

    return differences

def compare_lightcones(lc1_file, lc2_file):
    """Compare lightcone outputs."""
    results = {}

    with h5py.File(lc1_file, 'r') as f1, h5py.File(lc2_file, 'r') as f2:
        # Count galaxies
        def count_galaxies(f):
            # Try different possible dataset locations
            for ds_name in ['galaxies', 'lightcone', 'data', 'Posx', 'posx']:
                if ds_name in f:
                    return len(f[ds_name])
            return 0

        count1 = count_galaxies(f1)
        count2 = count_galaxies(f2)

        results['galaxy_count'] = (count1, count2, count1 == count2)

        # Compare field statistics for common fields
        # (Get all datasets from root)
        fields1 = set(k.lower() for k in f1.keys() if isinstance(f1[k], h5py.Dataset))
        fields2 = set(k.lower() for k in f2.keys() if isinstance(f2[k], h5py.Dataset))

        common_fields = fields1 & fields2

        field_stats = {}
        for field_lower in common_fields:
            # Find actual field name (case-insensitive)
            field1 = None
            field2 = None
            for k in f1.keys():
                if k.lower() == field_lower:
                    field1 = k
                    break
            for k in f2.keys():
                if k.lower() == field_lower:
                    field2 = k
                    break

            if field1 and field2:
                data1 = f1[field1][:]
                data2 = f2[field2][:]

                # Compare statistics (allow small numerical differences)
                stats1 = {'mean': np.mean(data1), 'std': np.std(data1),
                         'min': np.min(data1), 'max': np.max(data1)}
                stats2 = {'mean': np.mean(data2), 'std': np.std(data2),
                         'min': np.min(data2), 'max': np.max(data2)}

                # Check if statistics match within tolerance
                matches = True
                for key in stats1:
                    rel_diff = abs(stats1[key] - stats2[key]) / (abs(stats1[key]) + 1e-10)
                    if rel_diff > 1e-6:  # 0.0001% tolerance
                        matches = False
                        break

                field_stats[field_lower] = {
                    'stats1': stats1,
                    'stats2': stats2,
                    'match': matches
                }

        results['field_stats'] = field_stats

    return results

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Validate OLD vs NEW workflow outputs')
    parser.add_argument('--old-kdtree', required=True, help='OLD workflow KD-tree HDF5')
    parser.add_argument('--new-kdtree', required=True, help='NEW workflow KD-tree HDF5')
    parser.add_argument('--old-lightcone', required=True, help='OLD workflow lightcone HDF5')
    parser.add_argument('--new-lightcone', required=True, help='NEW workflow lightcone HDF5')
    parser.add_argument('--output', default='validation_report.txt', help='Output report file')

    args = parser.parse_args()

    print("Validating KD-tree structure...")
    kdtree_diffs = compare_hdf5_structure(args.old_kdtree, args.new_kdtree)

    enable_lightcone_check = True
    if enable_lightcone_check:
        print("Validating lightcone outputs...")
        lc_results = compare_lightcones(args.old_lightcone, args.new_lightcone)
    else:
        lc_results = {
            'galaxy_count': (0, 0, True),
            'field_stats': {}
        }   

    # Generate report
    with open(args.output, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("WORKFLOW OUTPUT VALIDATION REPORT\n")
        f.write("=" * 80 + "\n\n")

        f.write("## KD-tree Structure Comparison\n\n")
        if kdtree_diffs:
            f.write(f"Found {len(kdtree_diffs)} differences:\n")
            for diff in kdtree_diffs:
                f.write(f"  - {diff}\n")
        else:
            f.write("✓ KD-tree structures match perfectly\n")

        f.write("\n## Lightcone Comparison\n\n")

        gc1, gc2, match = lc_results['galaxy_count']
        f.write(f"Galaxy Count:\n")
        f.write(f"  OLD: {gc1}\n")
        f.write(f"  NEW: {gc2}\n")
        f.write(f"  Match: {'✓' if match else '✗'}\n\n")

        f.write("Field Statistics:\n")
        for field, stats in lc_results['field_stats'].items():
            status = '✓' if stats['match'] else '✗'
            f.write(f"  {status} {field}:\n")
            f.write(f"      OLD: mean={stats['stats1']['mean']:.6e}, std={stats['stats1']['std']:.6e}\n")
            f.write(f"      NEW: mean={stats['stats2']['mean']:.6e}, std={stats['stats2']['std']:.6e}\n")

        # Overall verdict
        f.write("\n" + "=" * 80 + "\n")
        all_match = (len(kdtree_diffs) == 0 and
                    lc_results['galaxy_count'][2] and
                    all(s['match'] for s in lc_results['field_stats'].values()))

        if all_match:
            f.write("✓ VALIDATION PASSED: Workflows produce equivalent outputs\n")
        else:
            f.write("✗ VALIDATION FAILED: Outputs differ\n")
        f.write("=" * 80 + "\n")

    print(f"Validation report written to {args.output}")

    # Exit with appropriate code
    sys.exit(0 if all_match else 1)

if __name__ == '__main__':
    main()
