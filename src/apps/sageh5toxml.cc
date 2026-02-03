/*!
 * @file sageh5toxml.cc
 * @brief Utility to read SAGE HDF5 files and extract Header/Simulation and
 * Snap_ group information
 * @author Generated for tao-lightcone-cli project
 */

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <iostream>
#include <libhpc/libhpc.hh>
#include <map>
#include <string>
#include <vector>

namespace po = boost::program_options;

class SageH5Reader {
private:
  hpc::h5::file _file;
  std::string _filename;
  std::string _output_filename;

public:
  SageH5Reader(const std::string &filename, const std::string &output_filename)
      : _filename(filename), _output_filename(output_filename) {
    try {
      _file.open(filename, H5F_ACC_RDONLY, hpc::mpi::comm::self);
    } catch (const std::exception &e) {
      throw std::runtime_error("Failed to open HDF5 file: " +
                               std::string(e.what()));
    }
  }

  ~SageH5Reader() {
    if (_file.is_open()) {
      _file.close();
    }
  }

  void readHeaderSimulation() {
    std::cout << "\n=== Reading Header/Simulation Information ===" << std::endl;

    try {
      // Try to open Header group
      hpc::h5::group header_group;
      header_group.open(_file, "Header");
      std::cout << "Found Header group" << std::endl;

      // Read Header group itself
      readGroupInfo(header_group, "Header");

      try {
        // Try to open Simulation subgroup
        hpc::h5::group sim_group;
        sim_group.open(header_group, "Simulation");
        std::cout << "Found Header/Simulation group" << std::endl;
        readGroupInfo(sim_group, "Header/Simulation");
      } catch (const std::exception &e) {
        std::cout << "No Simulation subgroup found in Header" << std::endl;
      }

      try {
        // Try to open Misc subgroup (where attributes might be)
        hpc::h5::group misc_group;
        misc_group.open(header_group, "Misc");
        std::cout << "Found Header/Misc group" << std::endl;
        readGroupInfo(misc_group, "Header/Misc");
      } catch (const std::exception &e) {
        std::cout << "No Misc subgroup found in Header" << std::endl;
      }

    } catch (const std::exception &e) {
      std::cout << "No Header group found in file: " << e.what() << std::endl;
    }
  }

  void readFirstSnapGroup() {
    std::cout << "\n=== Reading First Snap_ Group ===" << std::endl;

    try {
      // Get links (groups/datasets) in the root
      std::vector<std::string> links = _file.links();

      std::string first_snap_group;
      for (const auto &link : links) {
        if (link.substr(0, 5) == "Snap_") {
          first_snap_group = link;
          break;
        }
      }

      if (!first_snap_group.empty()) {
        std::cout << "Found first Snap_ group: " << first_snap_group
                  << std::endl;

        try {
          hpc::h5::group snap_group;
          snap_group.open(_file, first_snap_group);
          readGroupInfo(snap_group, first_snap_group);
        } catch (const std::exception &e) {
          std::cout << "Error opening Snap_ group " << first_snap_group << ": "
                    << e.what() << std::endl;
        }

      } else {
        std::cout << "No Snap_ groups found in file" << std::endl;
      }

    } catch (const std::exception &e) {
      std::cerr << "Error reading Snap_ groups: " << e.what() << std::endl;
    }
  }

private:
  void readGroupInfo(hpc::h5::group &group, const std::string &group_path) {
    std::cout << "\n--- Contents of " << group_path << " ---" << std::endl;

    try {
      // Get size of group (number of links)
      hsize_t size = group.size();
      std::cout << "Number of items in group: " << size << std::endl;

      // Read group attributes using raw HDF5 API
      readGroupAttributes(group, group_path);

      // Read dataset attributes if this is a Snap_ group
      if (group_path.find("Snap_") != std::string::npos) {
        readDatasetAttributes(group, group_path);
      }

    } catch (const std::exception &e) {
      std::cerr << "Error reading group contents: " << e.what() << std::endl;
    }
  }

  void readGroupAttributes(hpc::h5::group &group,
                           const std::string &group_path) {
    std::cout << "\n--- Attributes of " << group_path << " ---" << std::endl;

    try {
      hid_t group_id = group.id();

      // Get number of attributes using deprecated but more reliable function
      int num_attrs = H5Aget_num_attrs(group_id);
      if (num_attrs < 0) {
        std::cerr << "Failed to get number of attributes" << std::endl;
        return;
      }

      std::cout << "Number of attributes: " << num_attrs << std::endl;

      // Iterate through all attributes
      for (int i = 0; i < num_attrs; i++) {
        readSingleAttributeByIndex(group_id, i);
      }

    } catch (const std::exception &e) {
      std::cerr << "Error reading attributes: " << e.what() << std::endl;
    }
  }

  void readSingleAttributeByIndex(hid_t group_id, int attr_index) {
    try {
      // Open attribute by index using older, more compatible API
      hid_t attr_id = H5Aopen_by_idx(group_id, ".", H5_INDEX_NAME, H5_ITER_INC,
                                     attr_index, H5P_DEFAULT, H5P_DEFAULT);
      if (attr_id < 0) {
        // Try the older deprecated method
        attr_id = H5Aopen_idx(group_id, attr_index);
        if (attr_id < 0) {
          std::cerr << "Failed to open attribute " << attr_index << std::endl;
          return;
        }
      }

      // Get attribute name
      ssize_t name_size = H5Aget_name(attr_id, 0, nullptr);
      if (name_size <= 0) {
        H5Aclose(attr_id);
        return;
      }

      std::vector<char> name_buffer(name_size + 1);
      H5Aget_name(attr_id, name_size + 1, name_buffer.data());
      std::string attr_name(name_buffer.data());

      // Get attribute type and space info
      hid_t type_id = H5Aget_type(attr_id);
      hid_t space_id = H5Aget_space(attr_id);

      if (type_id >= 0 && space_id >= 0) {
        H5T_class_t type_class = H5Tget_class(type_id);
        size_t type_size = H5Tget_size(type_id);

        std::cout << "  " << attr_name << ": ";

        // Read the attribute value based on its type
        readAttributeValue(attr_id, type_class, type_size);

        std::cout << std::endl;

        H5Sclose(space_id);
        H5Tclose(type_id);
      }

      H5Aclose(attr_id);

    } catch (const std::exception &e) {
      std::cerr << "Error reading attribute " << attr_index << ": " << e.what()
                << std::endl;
    }
  }

  std::string readAttributeValue(hid_t attr_id, H5T_class_t type_class,
                                 size_t type_size) {
    try {
      switch (type_class) {
      case H5T_INTEGER: {
        if (type_size == sizeof(int)) {
          int value;
          if (H5Aread(attr_id, H5T_NATIVE_INT, &value) >= 0) {
            return std::to_string(value);
          }
        } else if (type_size == sizeof(long)) {
          long value;
          if (H5Aread(attr_id, H5T_NATIVE_LONG, &value) >= 0) {
            return std::to_string(value);
          }
        } else {
          std::cout << "[integer, " << type_size << " bytes]";
        }
        break;
      }
      case H5T_FLOAT: {
        if (type_size == sizeof(float)) {
          float value;
          if (H5Aread(attr_id, H5T_NATIVE_FLOAT, &value) >= 0) {
            return std::to_string(value);
          }
        } else if (type_size == sizeof(double)) {
          double value;
          if (H5Aread(attr_id, H5T_NATIVE_DOUBLE, &value) >= 0) {
            return std::to_string(value);
          }
        } else {
          return "[float, " + std::to_string(type_size);
        }
        break;
      }
      case H5T_STRING: {
        hid_t type_id = H5Aget_type(attr_id);
        if (H5Tis_variable_str(type_id)) {
          char *str_data;
          if (H5Aread(attr_id, type_id, &str_data) >= 0 &&
              str_data != nullptr) {
            std::string result = "\"" + std::string(str_data);
            H5free_memory(str_data);
          } else {
            return "[variable string]";
          }
        } else {
          size_t str_size = H5Tget_size(type_id);
          std::vector<char> str_data(str_size + 1);
          if (H5Aread(attr_id, type_id, str_data.data()) >= 0) {
            str_data[str_size] = '\0';
            return std::string(str_data.data());
          } else {
            return "[fixed string, " + std::to_string(str_size) + " chars]";
          }
        }
        H5Tclose(type_id);
        break;
      }
      default:
        return "[unknown type: " + std::to_string(type_class) + "]";

        break;
      }
    } catch (const std::exception &e) {
      return "[error reading value: " + std::string(e.what()) + "]";
    }
    return "[unreadable value]";
  }

  void readDatasetAttributes(hpc::h5::group &group,
                             const std::string &group_path) {
    std::ofstream results_file(_output_filename, std::ios::out);
    if (!results_file.is_open()) {
      std::cerr << "Error: Could not open " << _output_filename
                << " for writing" << std::endl;
      return;
    }
    std::cout << "\n--- Dataset Attributes in " << group_path << " ---"
              << std::endl;

    try {
      hid_t group_id = group.id();

      // Get number of objects in the group
      hsize_t num_objects = group.size();

      // Iterate through all objects in the group
      int order = 1;
      for (hsize_t i = 0; i < num_objects; i++) {
        // Get object name
        ssize_t name_size =
            H5Lget_name_by_idx(group_id, ".", H5_INDEX_NAME, H5_ITER_INC, i,
                               nullptr, 0, H5P_DEFAULT);
        if (name_size <= 0)
          continue;

        std::vector<char> name_buffer(name_size + 1);
        if (H5Lget_name_by_idx(group_id, ".", H5_INDEX_NAME, H5_ITER_INC, i,
                               name_buffer.data(), name_size + 1,
                               H5P_DEFAULT) < 0) {
          continue;
        }

        std::string obj_name(name_buffer.data());

        // Check if it's a dataset
        H5O_info2_t obj_info;
        if (H5Oget_info_by_name3(group_id, obj_name.c_str(), &obj_info,
                                 H5O_INFO_BASIC, H5P_DEFAULT) >= 0) {
          if (obj_info.type == H5O_TYPE_DATASET) {
            readSingleDatasetAttributes(group_id, obj_name, results_file,
                                        order);
            order++;
          }
        }
      }

    } catch (const std::exception &e) {
      std::cerr << "Error reading dataset attributes: " << e.what()
                << std::endl;
    }
    results_file.close();
  }

  void readSingleDatasetAttributes(hid_t group_id,
                                   const std::string &dataset_name,
                                   std::ofstream &results_file, int order) {
    try {
      // Open the dataset
      hid_t dataset_id = H5Dopen2(group_id, dataset_name.c_str(), H5P_DEFAULT);
      if (dataset_id < 0) {
        std::cerr << "Failed to open dataset: " << dataset_name << std::endl;
        return;
      }

      // Get dataset datatype information
      hid_t dtype_id = H5Dget_type(dataset_id);
      std::string datatype_info = getDataTypeInfo(dtype_id);

      // Get dataset space information
      hid_t space_id = H5Dget_space(dataset_id);
      hsize_t dataset_size = 0;
      if (space_id >= 0) {
        dataset_size = H5Sget_simple_extent_npoints(space_id);
        H5Sclose(space_id);
      }

      // Get number of attributes on this dataset
      int num_attrs = H5Aget_num_attrs(dataset_id);
      if (num_attrs < 0) {
        H5Tclose(dtype_id);
        H5Dclose(dataset_id);
        return;
      }
      results_file << "<Field type=\"" << datatype_info << "\"" << std::endl;
      results_file << " label=\"" << dataset_name << "\"" << std::endl;

      std::cout << "\nDataset '" << dataset_name << "':" << std::endl;
      std::cout << "  DATATYPE: " << datatype_info << std::endl;
      std::cout << "  SIZE: " << dataset_size << " elements" << std::endl;

      if (num_attrs > 0) {
        std::cout << "  ATTRIBUTES (" << num_attrs << "):" << std::endl;

        // Read all attributes of this dataset
        for (int i = 0; i < num_attrs; i++) {
          std::cout << "    ";
          readSingleDatasetAttribute(dataset_id, i, results_file);
        }
      } else {
        std::cout << "  ATTRIBUTES: None" << std::endl;
      }
      results_file << " order=\"" << order << "\"" << std::endl;
      results_file << " group=\"" << "general\"" << ">" << dataset_name
                   << "</Field>" << std::endl;

      H5Tclose(dtype_id);
      H5Dclose(dataset_id);

    } catch (const std::exception &e) {
      std::cerr << "Error reading attributes for dataset " << dataset_name
                << ": " << e.what() << std::endl;
    }
  }

  std::string getDataTypeInfo(hid_t dtype_id) {
    if (dtype_id < 0)
      return "Unknown";

    H5T_class_t type_class = H5Tget_class(dtype_id);
    size_t type_size = H5Tget_size(dtype_id);

    std::string result;

    switch (type_class) {
    case H5T_INTEGER: {
      H5T_sign_t sign = H5Tget_sign(dtype_id);
      H5T_order_t order = H5Tget_order(dtype_id);

      result = (sign == H5T_SGN_NONE) ? "Unsigned " : "Signed ";
      result += "Integer, " + std::to_string(type_size * 8) + "-bit";
      result += (order == H5T_ORDER_LE) ? " (Little Endian)" : " (Big Endian)";
      break;
    }
    case H5T_FLOAT: {
      H5T_order_t order = H5Tget_order(dtype_id);

      if (type_size == 4) {
        result = "float";
      } else if (type_size == 8) {
        result = "double";
      } else {
        result = "Float, " + std::to_string(type_size * 8) + "-bit";
      }
      result += (order == H5T_ORDER_LE) ? " (Little Endian)" : " (Big Endian)";
      break;
    }
    case H5T_STRING: {
      result = "String";
      if (H5Tis_variable_str(dtype_id)) {
        result += " (Variable Length)";
      } else {
        result += " (Fixed Length: " + std::to_string(type_size) + " chars)";
      }

      H5T_cset_t cset = H5Tget_cset(dtype_id);
      if (cset == H5T_CSET_ASCII) {
        result += " ASCII";
      } else if (cset == H5T_CSET_UTF8) {
        result += " UTF-8";
      }
      break;
    }
    case H5T_COMPOUND:
      result = "Compound Type (" + std::to_string(H5Tget_nmembers(dtype_id)) +
               " members)";
      break;
    case H5T_ARRAY:
      result = "Array Type";
      break;
    case H5T_ENUM:
      result = "Enumeration Type";
      break;
    default:
      result = "Unknown Type (Class: " + std::to_string(type_class) + ")";
      break;
    }

    return result;
  }

  void readSingleDatasetAttribute(hid_t dataset_id, int attr_index,
                                  std::ofstream &results_file) {
    try {
      // Open attribute by index using the older, more reliable method
      hid_t attr_id = H5Aopen_idx(dataset_id, attr_index);
      if (attr_id < 0) {
        std::cerr << "Failed to open attribute " << attr_index << std::endl;
        return;
      }

      // Get attribute name
      ssize_t name_size = H5Aget_name(attr_id, 0, nullptr);
      if (name_size <= 0) {
        H5Aclose(attr_id);
        return;
      }

      std::vector<char> name_buffer(name_size + 1);
      H5Aget_name(attr_id, name_size + 1, name_buffer.data());
      std::string attr_name(name_buffer.data());

      // Get attribute type and space info
      hid_t type_id = H5Aget_type(attr_id);
      hid_t space_id = H5Aget_space(attr_id);

      if (type_id >= 0 && space_id >= 0) {
        H5T_class_t type_class = H5Tget_class(type_id);
        size_t type_size = H5Tget_size(type_id);

        std::cout << attr_name << ": ";

        // Read the attribute value based on its type
        std::string attr_value =
            readAttributeValue(attr_id, type_class, type_size);

        std::cout << std::endl;
        results_file << " " << attr_name << "=\"" << attr_value << "\""
                     << std::endl;

        H5Sclose(space_id);
        H5Tclose(type_id);
      }

      H5Aclose(attr_id);

    } catch (const std::exception &e) {
      std::cerr << "Error reading dataset attribute " << attr_index << ": "
                << e.what() << std::endl;
    }
  }
};

void printUsage(const char *program_name) {
  std::cout << "Usage: " << program_name << " [options] <hdf5_file>"
            << std::endl;
  std::cout << std::endl;
  std::cout << "Reads SAGE HDF5 files and extracts information from:"
            << std::endl;
  std::cout << "  - Header/Simulation group" << std::endl;
  std::cout << "  - First Snap_ group found" << std::endl;
  std::cout << std::endl;
}

int main(int argc, char *argv[]) {
  // Initialize MPI via libhpc
  hpc::mpi::initialise(argc, argv);

  try {
    // Program options
    po::options_description desc("Options");
    desc.add_options()("help,h", "Show this help message")(
        "file,f", po::value<std::string>(), "HDF5 file to read")(
        "output,o", po::value<std::string>()->default_value("results.xml"),
        "Output XML file")("verbose,v", "Verbose output");

    po::positional_options_description pos_desc;
    pos_desc.add("file", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv)
                  .options(desc)
                  .positional(pos_desc)
                  .run(),
              vm);
    po::notify(vm);

    if (vm.count("help") || !vm.count("file")) {
      printUsage(argv[0]);
      std::cout << desc << std::endl;
      hpc::mpi::finalise();
      return vm.count("help") ? 0 : 1;
    }

    std::string filename = vm["file"].as<std::string>();
    std::string output_filename = vm["output"].as<std::string>();
    bool verbose = vm.count("verbose") > 0;

    std::cout << "SAGE HDF5 to XML Utility" << std::endl;
    std::cout << "========================" << std::endl;
    std::cout << "Input file: " << filename << std::endl;
    std::cout << "Output file: " << output_filename << std::endl;

    if (verbose) {
      std::cout << "Verbose mode enabled" << std::endl;
    }

    // Check if file exists
    if (!boost::filesystem::exists(filename)) {
      std::cerr << "Error: File does not exist: " << filename << std::endl;
      hpc::mpi::finalise();
      return 1;
    }

    // Create reader and process file
    SageH5Reader reader(filename, output_filename);

    // Read Header/Simulation group
    reader.readHeaderSimulation();

    // Read first Snap_ group
    reader.readFirstSnapGroup();

    std::cout << "\n=== Processing Complete ===" << std::endl;

  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    hpc::mpi::finalise();
    return 1;
  }

  hpc::mpi::finalise();
  return 0;
}
