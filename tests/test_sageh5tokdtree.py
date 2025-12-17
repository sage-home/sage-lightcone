import h5py
import numpy as np
import subprocess
import os
import sys

def create_dummy_sage_h5(filename):
    with h5py.File(filename, 'w') as f:
        # Create Header
        sim = f.create_group("Header/Simulation")
        for attr in ["Hubble_h", "Omega_b", "Omega_lambda", "Omega_m", "PartMass", "BoxSize"]:
            sim.attrs[attr] = 0.7  # Dummy values

        # Snapshot 0 (Early time)
        # 2 Galaxies. Gal 0 merges into Gal 0 of Snap 1. Gal 1 merges into Gal 0 of Snap 1.
        s0 = f.create_group("Snap_0")
        s0.create_dataset("Posx", data=np.array([0.1, 0.2]))
        s0.create_dataset("Posy", data=np.array([0.1, 0.2]))
        s0.create_dataset("Posz", data=np.array([0.1, 0.2]))
        s0.create_dataset("Descendant", data=np.array([0, 0])) # Both merge to index 0
        s0.create_dataset("GalaxyIndex", data=np.array([10, 11]))

        # Snapshot 1 (Late time)
        # 1 Galaxy. It is the descendant of the two above.
        s1 = f.create_group("Snap_1")
        s1.create_dataset("Posx", data=np.array([0.15]))
        s1.create_dataset("Posy", data=np.array([0.15]))
        s1.create_dataset("Posz", data=np.array([0.15]))
        s1.create_dataset("Descendant", data=np.array([-1])) # No descendant
        s1.create_dataset("GalaxyIndex", data=np.array([20]))

def test_sageh5tokdtree():
    input_file = "test_input.h5"
    output_file = "test_output.h5"
    
    if os.path.exists(input_file):
        os.remove(input_file)
    if os.path.exists(output_file):
        os.remove(output_file)

    create_dummy_sage_h5(input_file)
    
    # Run the tool
    # Assuming bin/sageh5tokdtree is relative to current working directory
    cmd = [
        "./bin/sageh5tokdtree", 
        "--sage", input_file, 
        "--output", output_file, 
        "--ppc", "10"
    ]
    print(f"Running: {' '.join(cmd)}")
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {e}")
        sys.exit(1)
    
    if not os.path.exists(output_file):
        print(f"Error: Output file {output_file} was not created.")
        sys.exit(1)

    with h5py.File(output_file, 'r') as f:
        # Check Groups
        if "cosmology" not in f:
            print("Error: 'cosmology' group missing.")
            sys.exit(1)
        if "lightcone" not in f:
            print("Error: 'lightcone' group missing.")
            sys.exit(1)
        if "data" not in f:
            print("Error: 'data' group missing.")
            sys.exit(1)
        
        # Check for KDTree index datasets
        for dset in ["snapshot_counts", "snapshot_displs", "snapshot_redshifts"]:
            if dset not in f:
                print(f"Error: '{dset}' dataset missing.")
                sys.exit(1)
        
        # Verify values
        counts = f["snapshot_counts"][:]
        displs = f["snapshot_displs"][:]
        
        print(f"Counts: {counts}")
        print(f"Displs: {displs}")
        
        if not np.array_equal(counts, [2, 1]):
             print("Error: snapshot_counts mismatch. Expected [2, 1]")
             sys.exit(1)
             
        if not np.array_equal(displs, [0, 2]):
             print("Error: snapshot_displs mismatch. Expected [0, 2]")
             sys.exit(1)

        # Check lightcone structure
        if "lightcone/snapshot000" not in f:
            print("Error: 'lightcone/snapshot000' missing.")
            sys.exit(1)
            
        lc0 = f["lightcone/snapshot000"]
        for dset in ["bounds", "splits", "cell_counts", "cell_offs"]:
            if dset not in lc0:
                print(f"Error: 'lightcone/snapshot000/{dset}' missing.")
                sys.exit(1)

        # Check data group content
        data = f["data"]
        for dset in ["Posx", "Posy", "Posz", "subsize", "GalaxyIndex"]:
            if dset not in data:
                print(f"Error: 'data/{dset}' missing.")
                sys.exit(1)
                
        # Verify subsize values
        subsize = data["subsize"][:]
        print(f"Subsize: {subsize}")
        
        # Expected: Snap 0 (2 gals) -> 1, 1. Snap 1 (1 gal) -> 3.
        # Total 3 values.
        if len(subsize) != 3:
             print(f"Error: subsize length mismatch. Expected 3, got {len(subsize)}")
             sys.exit(1)
             
        # Since order is permuted, we check if values are present.
        # We expect two 1s and one 3.
        unique, counts = np.unique(subsize, return_counts=True)
        subsize_counts = dict(zip(unique, counts))
        
        if subsize_counts.get(1, 0) != 2:
             print(f"Error: Expected two galaxies with subsize 1. Got {subsize_counts.get(1, 0)}")
             sys.exit(1)
             
        if subsize_counts.get(3, 0) != 1:
             print(f"Error: Expected one galaxy with subsize 3. Got {subsize_counts.get(3, 0)}")
             sys.exit(1)

        print("Output file structure verified (Groups and Datasets exist).")
        
        # Check Subsize Logic
        # Snap 1, Gal 0 should have subsize = 1 (itself) + 1 (Snap0_Gal0) + 1 (Snap0_Gal1) = 3
        # Wait, standard subsize is usually inclusive of self + all progenitors.
        # Your C++ logic: subsize[k] = 1 + incoming_subsize[k]
        # Snap 0: incoming is 0. subsize = 1.
        # Snap 1: incoming is sum of progenitors (1+1=2). subsize = 1 + 2 = 3.
        
        # We need to read the reordered data to verify. 
        # Since we don't write 'subsize' to disk yet in your C++ code (it was commented out/missing),
        # this test will fail until we add the write logic.
        pass

    # Cleanup
    # if os.path.exists(input_file):
    #     os.remove(input_file)
    # if os.path.exists(output_file):
    #     os.remove(output_file)
    
    print("Test passed successfully!")

if __name__ == "__main__":
    test_sageh5tokdtree()
