#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <memory>
#include <iomanip>
#include <hdf5.h>
#include <mpi.h>

// Structure to hold dataset information
struct DatasetInfo {
    std::string name;
    hid_t type_id;
    H5T_class_t type_class;
    size_t type_size;
    hsize_t num_elements;
    std::vector<char> data;
};

// Callback function for H5Literate to collect dataset names
herr_t dataset_iterator(hid_t group_id, const char* name, const H5L_info2_t* info, void* operator_data) {
    std::vector<std::string>* dataset_names = static_cast<std::vector<std::string>*>(operator_data);
    
    // Check if this is a dataset
    H5O_info2_t obj_info;
    if (H5Oget_info_by_name3(group_id, name, &obj_info, H5O_INFO_BASIC, H5P_DEFAULT) >= 0) {
        if (obj_info.type == H5O_TYPE_DATASET) {
            dataset_names->push_back(std::string(name));
        }
    }
    
    return 0; // Continue iteration
}

// Function to get the string representation of a value based on its type
std::string value_to_string(const void* data, hid_t type_id, size_t element_index) {
    H5T_class_t type_class = H5Tget_class(type_id);
    size_t type_size = H5Tget_size(type_id);
    const char* byte_data = static_cast<const char*>(data) + (element_index * type_size);
    
    std::ostringstream oss;
    
    switch (type_class) {
        case H5T_INTEGER: {
            if (H5Tequal(type_id, H5T_NATIVE_INT) > 0) {
                oss << *reinterpret_cast<const int*>(byte_data);
            } else if (H5Tequal(type_id, H5T_NATIVE_LONG) > 0) {
                oss << *reinterpret_cast<const long*>(byte_data);
            } else if (H5Tequal(type_id, H5T_NATIVE_LLONG) > 0) {
                oss << *reinterpret_cast<const long long*>(byte_data);
            } else if (H5Tequal(type_id, H5T_NATIVE_SHORT) > 0) {
                oss << *reinterpret_cast<const short*>(byte_data);
            } else if (H5Tequal(type_id, H5T_NATIVE_CHAR) > 0) {
                oss << static_cast<int>(*reinterpret_cast<const char*>(byte_data));
            } else {
                // Default to int for unknown integer types
                oss << *reinterpret_cast<const int*>(byte_data);
            }
            break;
        }
        case H5T_FLOAT: {
            if (H5Tequal(type_id, H5T_NATIVE_FLOAT) > 0) {
                oss << std::scientific << std::setprecision(6) << *reinterpret_cast<const float*>(byte_data);
            } else if (H5Tequal(type_id, H5T_NATIVE_DOUBLE) > 0) {
                oss << std::scientific << std::setprecision(15) << *reinterpret_cast<const double*>(byte_data);
            } else {
                // Default to double for unknown float types
                oss << std::scientific << std::setprecision(15) << *reinterpret_cast<const double*>(byte_data);
            }
            break;
        }
        case H5T_STRING: {
            // Handle string data
            if (H5Tis_variable_str(type_id)) {
                // Variable length string
                char** str_data = const_cast<char**>(reinterpret_cast<const char* const*>(byte_data));
                if (str_data && *str_data) {
                    oss << *str_data;
                }
            } else {
                // Fixed length string
                std::string str(byte_data, type_size);
                // Remove null terminators
                size_t null_pos = str.find('\0');
                if (null_pos != std::string::npos) {
                    str = str.substr(0, null_pos);
                }
                oss << str;
            }
            break;
        }
        default:
            oss << "[unsupported_type]";
            break;
    }
    
    return oss.str();
}

int main(int argc, char* argv[]) {
    // Initialize MPI (required for parallel HDF5)
    int mpi_initialized;
    MPI_Initialized(&mpi_initialized);
    if (!mpi_initialized) {
        MPI_Init(&argc, &argv);
    }
    
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input_hdf5_file> <output_csv_file>" << std::endl;
        if (!mpi_initialized) MPI_Finalize();
        return 1;
    }
    
    const std::string input_file = argv[1];
    const std::string output_file = argv[2];
    
    // Open HDF5 file
    hid_t file_id = H5Fopen(input_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        std::cerr << "Error: Could not open HDF5 file: " << input_file << std::endl;
        if (!mpi_initialized) MPI_Finalize();
        return 1;
    }
    
    try {
        // Get all dataset names
        std::vector<std::string> dataset_names;
        if (H5Literate2(file_id, H5_INDEX_NAME, H5_ITER_INC, nullptr, dataset_iterator, &dataset_names) < 0) {
            throw std::runtime_error("Failed to iterate over datasets");
        }
        
        if (dataset_names.empty()) {
            std::cerr << "Warning: No datasets found in the HDF5 file" << std::endl;
            H5Fclose(file_id);
            if (!mpi_initialized) MPI_Finalize();
            return 0;
        }
        
        std::cout << "Found " << dataset_names.size() << " datasets:" << std::endl;
        for (const auto& name : dataset_names) {
            std::cout << "  - " << name << std::endl;
        }
        
        // Load all datasets
        std::vector<DatasetInfo> datasets;
        hsize_t num_rows = 0;
        
        for (const std::string& name : dataset_names) {
            DatasetInfo info;
            info.name = name;
            
            // Open dataset
            hid_t dataset_id = H5Dopen2(file_id, name.c_str(), H5P_DEFAULT);
            if (dataset_id < 0) {
                std::cerr << "Warning: Could not open dataset: " << name << std::endl;
                continue;
            }
            
            // Get datatype and dataspace
            info.type_id = H5Dget_type(dataset_id);
            info.type_class = H5Tget_class(info.type_id);
            info.type_size = H5Tget_size(info.type_id);
            
            hid_t dataspace_id = H5Dget_space(dataset_id);
            int ndims = H5Sget_simple_extent_ndims(dataspace_id);
            
            if (ndims != 1) {
                std::cerr << "Warning: Dataset " << name << " is not 1D, skipping" << std::endl;
                H5Sclose(dataspace_id);
                H5Tclose(info.type_id);
                H5Dclose(dataset_id);
                continue;
            }
            
            hsize_t dims[1];
            H5Sget_simple_extent_dims(dataspace_id, dims, nullptr);
            info.num_elements = dims[0];
            
            // Check if all datasets have the same length
            if (num_rows == 0) {
                num_rows = info.num_elements;
            } else if (num_rows != info.num_elements) {
                std::cerr << "Error: Dataset " << name << " has " << info.num_elements 
                         << " elements, but expected " << num_rows << std::endl;
                H5Sclose(dataspace_id);
                H5Tclose(info.type_id);
                H5Dclose(dataset_id);
                H5Fclose(file_id);
                if (!mpi_initialized) MPI_Finalize();
                return 1;
            }
            
            // Read data
            info.data.resize(info.num_elements * info.type_size);
            if (H5Dread(dataset_id, info.type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, info.data.data()) < 0) {
                std::cerr << "Error: Failed to read dataset: " << name << std::endl;
                H5Sclose(dataspace_id);
                H5Tclose(info.type_id);
                H5Dclose(dataset_id);
                H5Fclose(file_id);
                if (!mpi_initialized) MPI_Finalize();
                return 1;
            }
            
            datasets.push_back(std::move(info));
            
            // Clean up
            H5Sclose(dataspace_id);
            H5Dclose(dataset_id);
        }
        
        // Write CSV file
        std::ofstream csv_file(output_file);
        if (!csv_file.is_open()) {
            H5Fclose(file_id);
            if (!mpi_initialized) MPI_Finalize();
            throw std::runtime_error("Could not open output CSV file: " + output_file);
        }
        
        // Write header
        for (size_t i = 0; i < datasets.size(); ++i) {
            if (i > 0) csv_file << ",";
            csv_file << datasets[i].name;
        }
        csv_file << std::endl;
        
        // Write data rows
        for (hsize_t row = 0; row < num_rows; ++row) {
            for (size_t col = 0; col < datasets.size(); ++col) {
                if (col > 0) csv_file << ",";
                
                std::string value = value_to_string(datasets[col].data.data(), 
                                                  datasets[col].type_id, row);
                
                // Escape commas and quotes in CSV
                if (value.find(',') != std::string::npos || 
                    value.find('"') != std::string::npos ||
                    value.find('\n') != std::string::npos) {
                    // Replace quotes with double quotes and wrap in quotes
                    std::string escaped;
                    escaped.reserve(value.size() + 10);
                    escaped += '"';
                    for (char c : value) {
                        if (c == '"') {
                            escaped += "\"\"";
                        } else {
                            escaped += c;
                        }
                    }
                    escaped += '"';
                    csv_file << escaped;
                } else {
                    csv_file << value;
                }
            }
            csv_file << std::endl;
        }
        
        // Clean up dataset types
        for (auto& dataset : datasets) {
            H5Tclose(dataset.type_id);
        }
        
        std::cout << "Successfully wrote " << num_rows << " rows and " 
                  << datasets.size() << " columns to " << output_file << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        H5Fclose(file_id);
        if (!mpi_initialized) MPI_Finalize();
        return 1;
    }
    
    H5Fclose(file_id);
    if (!mpi_initialized) MPI_Finalize();
    return 0;
}
