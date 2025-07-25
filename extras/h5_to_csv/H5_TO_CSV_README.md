# HDF5 to CSV Converter

## Description

The `h5_to_csv` tool converts HDF5 datasets to CSV format, making it easy to analyze simulation output data in spreadsheet applications or with data analysis tools.

## Building

The tool is built automatically as part of the main project:

```bash
source ./setup.sh
./run_build.sh
cd bin
make install
```

After installation, the executable will be available in `./runtime/bin/h5_to_csv`.

## Usage

```bash
./runtime/bin/h5_to_csv <input_hdf5_file> <output_csv_file> [dataset1] [dataset2] ...
```

### Parameters

- `input_hdf5_file`: Path to the input HDF5 file
- `output_csv_file`: Path to the output CSV file  
- `dataset1, dataset2, ...`: Names of the datasets to extract (optional - if not provided, you'll need to modify the code for automatic discovery)

### Example

```bash
# Extract specific datasets
./runtime/bin/h5_to_csv simulation_output.h5 results.csv galaxy_id x_pos y_pos z_pos mass

# This will create a CSV file with columns: galaxy_id, x_pos, y_pos, z_pos, mass
```

## Requirements

- All specified datasets must have the same number of elements (same length)
- Datasets must be 1-dimensional arrays
- Supported data types: integers, floats, doubles

## Output Format

The tool creates a CSV file with:
- Header row containing dataset names
- One row per data element
- Proper CSV escaping for special characters
- Scientific notation for floating-point numbers

## Implementation Notes

The tool uses the project's `libhpc` HDF5 wrapper library for reading data. Currently, dataset names must be provided as command-line arguments. To enable automatic dataset discovery, you would need to extend the implementation with HDF5 iteration functions.

## Compilation as Standalone Tool

If you need to compile the tool independently of the main project:

```bash
# Make sure HDF5 development libraries are installed
# Then compile with h5c++ wrapper
./compile_h5_to_csv.sh for ozstar
or perhaps
h5c++ -std=c++17 -O2 -o h5_to_csv h5_to_csv.cpp
```

