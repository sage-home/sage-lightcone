#!/usr/bin/env python3
"""
Plot 3D lightcone data from HDF5 file
Reads Posx, Posy, Posz fields and creates a 3D scatter plot

Usage:
    python plot_lightcone.py <hdf5_file> [color_field]

Arguments:
    hdf5_file   : Path to the HDF5 file
    color_field : Field name to color by (default: StellarMass, case-insensitive)

Note: For headless environments (HPC clusters), plots are automatically saved to PNG files.
      Field names are matched case-insensitively (e.g., 'StellarMass' = 'stellarmass' = 'stellar_mass').
"""

import h5py
import numpy as np
# Use non-interactive backend for headless environments (HPC clusters)
import matplotlib
matplotlib.use('Agg')  # Force non-interactive backend
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys

def find_field_case_insensitive(h5_obj, field_name):
    """
    Find field using case-insensitive search.

    Args:
        h5_obj: HDF5 file, group, or dataset object
        field_name: Field name to search for

    Returns:
        Actual field name as it appears in HDF5, or None if not found
    """
    # For compound datasets, check dtype.names
    if hasattr(h5_obj, 'dtype') and h5_obj.dtype.names:
        # Try exact match first
        if field_name in h5_obj.dtype.names:
            return field_name

        # Try case-insensitive
        field_lower = field_name.lower()
        for key in h5_obj.dtype.names:
            if key.lower() == field_lower:
                return key

    # For groups/files, check keys
    elif hasattr(h5_obj, 'keys'):
        # Try exact match first
        if field_name in h5_obj:
            return field_name

        # Try case-insensitive
        field_lower = field_name.lower()
        for key in h5_obj.keys():
            if key.lower() == field_lower:
                return key

    return None

def plot_lightcone_3d(hdf5_file, color_field="StellarMass"):
    """
    Read lightcone data from HDF5 file and create a 3D scatter plot

    Args:
        hdf5_file (str): Path to the HDF5 file
        color_field (str): Field name to color the points by (case-insensitive, default: StellarMass)

    Note:
        Field names are matched case-insensitively. For example, 'StellarMass',
        'stellarmass', and 'stellar_mass' will all match the same field.
    """
    # File path is now passed as argument
    # hdf5_file is already set by the function parameter
    
    try:
        # Open the HDF5 file
        with h5py.File(hdf5_file, 'r') as f:
            print(f"Successfully opened {hdf5_file}")
            
            # Print available datasets for debugging
            print("Available datasets:")
            def print_datasets(name, obj):
                if isinstance(obj, h5py.Dataset):
                    print(f"  {name}: shape={obj.shape}, dtype={obj.dtype}")
            f.visititems(print_datasets)
            
            # Try to read position data and color field - check different possible dataset names
            dataset_names = ['galaxies', 'lightcone', 'data', '/']
            
            posx = posy = posz = color_data = None
            
            for dataset_name in dataset_names:
                try:
                    if dataset_name in f or dataset_name == '/':
                        if dataset_name == '/':
                            # Check if fields are directly in root (case-insensitive)
                            posx_field = find_field_case_insensitive(f, 'posx')
                            if posx_field:
                                posy_field = find_field_case_insensitive(f, 'posy')
                                posz_field = find_field_case_insensitive(f, 'posz')
                                color_field_actual = find_field_case_insensitive(f, color_field)

                                posx = f[posx_field][:]
                                posy = f[posy_field][:]
                                posz = f[posz_field][:]
                                color_data = f[color_field_actual][:] if color_field_actual else None

                                print(f"Found position data in root directory")
                                if color_field_actual and color_field_actual != color_field:
                                    print(f"Using field '{color_field_actual}' (requested: '{color_field}')")
                                break
                        else:
                            # Check if fields are in a dataset (case-insensitive)
                            ds = f[dataset_name]
                            posx_field = find_field_case_insensitive(ds, 'posx')
                            if hasattr(ds, 'dtype') and ds.dtype.names and posx_field:
                                data = ds[:]
                                posy_field = find_field_case_insensitive(ds, 'posy')
                                posz_field = find_field_case_insensitive(ds, 'posz')
                                color_field_actual = find_field_case_insensitive(ds, color_field)

                                posx = data[posx_field]
                                posy = data[posy_field]
                                posz = data[posz_field]
                                color_data = data[color_field_actual] if color_field_actual else None

                                print(f"Found position data in dataset: {dataset_name}")
                                if color_field_actual and color_field_actual != color_field:
                                    print(f"Using field '{color_field_actual}' (requested: '{color_field}')")
                                break
                except (KeyError, ValueError) as e:
                    continue
            
            if posx is None:
                print("Could not find posx, posy, posz fields. Available fields:")
                # Try to find any dataset with structured data
                for key in f.keys():
                    try:
                        dataset = f[key]
                        if hasattr(dataset, 'dtype') and dataset.dtype.names:
                            print(f"Dataset '{key}' has fields: {list(dataset.dtype.names)}")
                    except:
                        pass
                return
            
            if color_data is None:
                print(f"Warning: {color_field} field not found, will color by z-position instead")
                color_data = posz
                color_label = 'Z Position'
            else:
                color_label = color_field.replace('_', ' ').title()
            
            print(f"Loaded {len(posx)} points")
            print(f"X range: [{np.min(posx):.2f}, {np.max(posx):.2f}]")
            print(f"Y range: [{np.min(posy):.2f}, {np.max(posy):.2f}]")
            print(f"Z range: [{np.min(posz):.2f}, {np.max(posz):.2f}]")
            if color_data is not None and not np.array_equal(color_data, posz):
                print(f"{color_label} range: [{np.min(color_data):.2e}, {np.max(color_data):.2e}]")
            
    except FileNotFoundError:
        print(f"Error: File {hdf5_file} not found")
        return
    except Exception as e:
        print(f"Error reading HDF5 file: {e}")
        return
    
    # Create 3D scatter plot
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111, projection='3d')
    
    # Sample data if too many points (for performance)
    max_points = 10000
    if len(posx) > max_points:
        print(f"Sampling {max_points} points from {len(posx)} total points")
        indices = np.random.choice(len(posx), max_points, replace=False)
        posx = posx[indices]
        posy = posy[indices]
        posz = posz[indices]
        color_data = color_data[indices]
    
    # Create scatter plot colored by the specified field
    scatter = ax.scatter(posx, posy, posz, 
                        c=color_data,  # Color by specified field
                        cmap='plasma',   # Good colormap for most data
                        alpha=0.6, 
                        s=1)
    
    # Set labels and title
    ax.set_xlabel('X Position')
    ax.set_ylabel('Y Position')
    ax.set_zlabel('Z Position')
    ax.set_title(f'3D Lightcone Galaxy Distribution (Colored by {color_label})')
    
    # Set equal aspect ratio for all axes
    # Get the range of each axis
    x_range = np.max(posx) - np.min(posx)
    y_range = np.max(posy) - np.min(posy)
    z_range = np.max(posz) - np.min(posz)
    
    # Find the maximum range and set all axes to this range
    max_range = max(x_range, y_range, z_range)
    
    # Set the same scale for all axes
    x_center = (np.max(posx) + np.min(posx)) / 2
    y_center = (np.max(posy) + np.min(posy)) / 2
    z_center = (np.max(posz) + np.min(posz)) / 2
    
    ax.set_xlim(x_center - max_range/2, x_center + max_range/2)
    ax.set_ylim(y_center - max_range/2, y_center + max_range/2)
    ax.set_zlim(z_center - max_range/2, z_center + max_range/2)
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.5, aspect=20, label=color_label)
    
    # Use log scale for colorbar if the range is large
    if color_label != 'Z Position':
        positive_data = color_data[color_data > 0]
        if len(positive_data) > 0 and np.max(color_data) / np.min(positive_data) > 100:
            import matplotlib.ticker as ticker
            cbar.formatter = ticker.LogFormatterSciNotation()
    
    plt.tight_layout()
    
    # Save the plot to a file (HPC environment - no display available)
    output_filename = f"lightcone_3d_{color_field}.png"
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"Plot saved as: {output_filename}")
    
    # Close the figure to free memory
    plt.close()

def main():
    """Main function to handle command-line arguments."""
    import sys

    # Check for correct number of arguments
    if len(sys.argv) < 2:
        print("Usage: python plot_lightcone.py <path_to_hdf5_file> [color_field]")
        print("Example: python plot_lightcone.py output/mymillennium_lightcone.h5 StellarMass")
        print("         python plot_lightcone.py output/mymillennium_lightcone.h5 SnapNum")
        print("Default color field is 'StellarMass' if not specified")
        print("Note: Field names are case-insensitive")
        sys.exit(1)

    # Parse command-line arguments
    hdf5_file = sys.argv[1]
    color_field = sys.argv[2] if len(sys.argv) > 2 else "StellarMass"
    
    # Check if file exists
    import os
    if not os.path.exists(hdf5_file):
        print(f"Error: File '{hdf5_file}' does not exist")
        sys.exit(1)
    
    try:
        # Call the plotting function
        plot_lightcone_3d(hdf5_file, color_field)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

# Usage example
if __name__ == "__main__":
    main()