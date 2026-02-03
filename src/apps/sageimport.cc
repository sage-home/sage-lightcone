#include <deque>
#include <fstream>
#include <iostream>
#include <libhpc/libhpc.hh>
#include <libtao/tao.hh>
#include <map>
#include <pugixml.hpp>
#include <string>
#include <vector>

struct SageField {
  std::string name;
  std::string type;
  std::string label;
  std::string description;
  std::string units;
  std::string group;
  int order;
};

struct Config {
  std::string input_file;
  std::string output_file;
  std::string tree_mapping[3]; // 0: GlobalIndex, 1: Descendant, 2: SnapNum
  std::vector<SageField> fields;
};

struct Node {
  size_t id;
  long long original_index;
  Node *descendant = nullptr;
  std::vector<Node *> progenitors;

  long long bfs_idx = 0;
  long long dfs_idx = 0;
  long long subtree_count = 0;
};

void dfs_visit(Node *u, long long &counter) {
  u->dfs_idx = counter++;
  u->subtree_count = 1;

  // Sort progenitors descending to match Python's stack-based DFS (LIFO)
  std::sort(u->progenitors.begin(), u->progenitors.end(), [](Node *a, Node *b) {
    return a->original_index > b->original_index;
  });

  for (auto *v : u->progenitors) {
    dfs_visit(v, counter);
    u->subtree_count += v->subtree_count;
  }
}

class SageImportApp : public hpc::mpi::application {
public:
  SageImportApp(int argc, char *argv[]) : hpc::mpi::application(argc, argv) {
    options().add_options()("settings",
                            hpc::po::value<std::string>(&_settings_file),
                            "XML Settings File")(
        "resume", hpc::po::value<bool>(&_resume)->default_value(false),
        "Resume Processing");

    // hpc::po::positional_options_description p;
    // p.add("settings", 1);
    // p.add("resume", 1);

    // parse_options(argc, argv, p);
    parse_options(argc, argv);
  }

  void operator()() {
    if (_settings_file.empty()) {
      std::cerr << "Error: No settings file provided." << std::endl;
      return;
    }

    parse_settings();
    process();
  }

private:
  std::string _settings_file;
  bool _resume;
  Config _config;

  // Offsets for key fields in input buffer
  size_t _offset_global_index = 0;
  size_t _offset_descendant = 0;
  size_t _offset_snap_num = 0;

  // Types for key fields (to handle different integer sizes if needed)
  hid_t _type_global_index = H5T_NATIVE_LLONG;
  hid_t _type_descendant = H5T_NATIVE_LLONG;
  hid_t _type_snap_num = H5T_NATIVE_INT;

  // Offsets for new fields in output buffer
  size_t _out_offset_global_tree_id = 0;
  size_t _out_offset_bfs = 0;
  size_t _out_offset_dfs = 0;
  size_t _out_offset_subtree = 0;
  size_t _out_offset_local_id = 0;
  bool _add_local_id = false;

  size_t _input_element_size = 0;
  size_t _output_element_size = 0;

  void parse_settings() {
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(_settings_file.c_str());
    if (!result) {
      throw std::runtime_error("Failed to parse settings file: " +
                               std::string(result.description()));
    }

    pugi::xml_node root = doc.first_child(); // Settings

    // 1. Parse Sage Fields (Index 0)
    pugi::xml_node sage_fields_node =
        root.child("sageFields"); // Or by index if tag names vary
    if (!sage_fields_node)
      sage_fields_node = *root.begin(); // Fallback to first child

    int order = 1;
    for (pugi::xml_node field : sage_fields_node.children()) {
      SageField sf;
      sf.name = field.text().as_string();
      sf.type = field.attribute("Type").as_string();
      sf.label = field.attribute("label").as_string();
      sf.description = field.attribute("description").as_string();
      sf.units = field.attribute("units").as_string();
      sf.group = field.attribute("group").as_string();
      sf.order = field.attribute("order").as_int(order++);
      _config.fields.push_back(sf);
    }

    // 2. Parse Running Settings (Index 2 usually, but let's search by tag)
    pugi::xml_node running_settings = root.child("RunningSettings");
    if (!running_settings) {
      // Fallback by index if needed, but tag search is safer if tags are
      // consistent
      auto it = root.begin();
      std::advance(it, 2);
      running_settings = *it;
    }

    _config.input_file = running_settings.child("InputFile").text().as_string();
    _config.output_file =
        running_settings.child("OutputFile").text().as_string();

    // 3. Parse Tree Mappings (Index 3)
    pugi::xml_node tree_mappings = root.child("TreeMappings");
    if (!tree_mappings) {
      tree_mappings = root.child("TreeTraversal");
    }
    if (!tree_mappings) {
      auto it = root.begin();
      std::advance(it, 3);
      tree_mappings = *it;
    }

    int i = 0;
    for (pugi::xml_node mapping : tree_mappings.children()) {
      if (i < 3) {
        _config.tree_mapping[i] = mapping.text().as_string();
        i++;
      }
    }

    LOGBLOCK("Configuration");
    LOGDLN("Input File: ", _config.input_file);
    LOGDLN("Output File: ", _config.output_file);
    LOGDLN("Tree Mapping 0 (GlobalIndex): ", _config.tree_mapping[0]);
    LOGDLN("Tree Mapping 1 (Descendant): ", _config.tree_mapping[1]);
    LOGDLN("Tree Mapping 2 (SnapNum): ", _config.tree_mapping[2]);
    LOGDLN("Number of Fields: ", _config.fields.size());
  }

  void process() {
    LOGBLOCK("Processing");

    int my_rank = rank();
    int n_ranks = size();

    // Rank 0 prepares the output file
    if (my_rank == 0) {
      hpc::h5::file in_file(_config.input_file, H5F_ACC_RDONLY);
      if (!in_file.is_open()) {
        std::cerr << "Rank 0: Failed to open input file: " << _config.input_file
                  << std::endl;
        throw std::runtime_error("Rank 0: Failed to open input file");
      }
      hpc::h5::file out_file(_config.output_file, H5F_ACC_TRUNC);
      if (!out_file.is_open()) {
        std::cerr << "Rank 0: Failed to open output file: "
                  << _config.output_file << std::endl;
        throw std::runtime_error("Rank 0: Failed to open output file");
      }

      LOGDLN("Copying groups...");
      // Use try-catch to handle missing groups gracefully if needed
      try {
        hpc::h5::copy(in_file, "cosmology", out_file);
      } catch (...) {
      }
      try {
        hpc::h5::copy(in_file, "snapshot_redshifts", out_file);
      } catch (...) {
      }
      try {
        hpc::h5::copy(in_file, "tree_counts", out_file);
      } catch (...) {
      }
      try {
        hpc::h5::copy(in_file, "tree_displs", out_file);
      } catch (...) {
      }

      LOGDLN("Groups copied.");

      create_sidecar_xml();
      in_file.close();
      out_file.close();
    }

    // Wait for rank 0
    hpc::mpi::comm::world.barrier();

    // Process Trees
    process_trees(_config.input_file, _config.output_file, my_rank, n_ranks);
  }

  void create_sidecar_xml() {
    std::string filename =
        _config.output_file.substr(0, _config.output_file.length() - 3) +
        ".xml";
    std::ofstream xml_file(filename);
    xml_file << "  <sageinput>\n";
    int count = 1;
    for (const auto &field : _config.fields) {
      xml_file << "    <Field Type=\"" << field.type << "\"\n";
      xml_file << "      label=\"" << field.label << "\"\n";
      xml_file << "      description=\"" << field.description << "\"\n";
      xml_file << "      order=\"" << field.order << "\"\n";
      xml_file << "      units=\"" << field.units << "\"\n";
      xml_file << "      group=\"" << field.group << "\">" << field.name
               << "</Field>\n";
      count++;
    }

    auto add_field = [&](const std::string &name, const std::string &type) {
      xml_file << "    <Field Type=\"" << type << "\"\n";
      xml_file << "      label=\"" << name << "\"\n";
      xml_file << "      description=\"" << name << "\"\n";
      xml_file << "      order=\"" << count++ << "\"\n";
      xml_file << "      units=\"\"\n";
      xml_file << "      group=\"internal\">" << name << "</Field>\n";
    };

    add_field("globaltreeid", "long long");
    add_field("breadthfirst_traversalorder", "long long");
    add_field("depthfirst_traversalorder", "long long");
    add_field("subtree_count", "long long");
    add_field("localgalaxyid", "int");

    xml_file << "  </sageinput>\n";
    xml_file.close();
  }

  void process_trees(const std::string &in_fname, const std::string &out_fname,
                     int rank, int n_ranks) {
    hpc::h5::file in_file(in_fname, H5F_ACC_RDONLY, hpc::mpi::comm::world);
    if (!in_file.is_open()) {
      std::cerr << "Rank " << rank
                << ": Failed to open input file (MPI): " << in_fname
                << std::endl;
      throw std::runtime_error("Failed to open input file (MPI)");
    }
    hpc::h5::file out_file(out_fname, H5F_ACC_RDWR, hpc::mpi::comm::world);
    if (!out_file.is_open()) {
      std::cerr << "Rank " << rank
                << ": Failed to open output file (MPI): " << out_fname
                << std::endl;
      throw std::runtime_error("Failed to open output file (MPI)");
    }

    hpc::h5::dataset in_ds(in_file, "galaxies");
    hpc::h5::dataset out_ds;

    // Read tree counts
    hpc::h5::dataset counts_ds(in_file, "tree_counts");
    hsize_t n_trees = counts_ds.dataspace().size();

    std::vector<int> tree_counts(n_trees);
    counts_ds.read(tree_counts.data(), hpc::h5::datatype::native_int);

    // We need to know the offset of each tree in the galaxies dataset
    std::vector<long long> tree_offsets(n_trees + 1, 0);
    for (size_t i = 0; i < n_trees; ++i) {
      tree_offsets[i + 1] = tree_offsets[i] + tree_counts[i];
    }

    // Determine which trees this rank processes
    size_t trees_per_rank = n_trees / n_ranks;
    size_t start_tree = rank * trees_per_rank;
    size_t end_tree =
        (rank == n_ranks - 1) ? n_trees : (rank + 1) * trees_per_rank;

    // Prepare memory types for reading/writing
    hpc::h5::datatype in_file_type = in_ds.datatype();

    // Create a native memory type for reading to handle endianness conversion
    hid_t in_mem_type_id =
        H5Tget_native_type(in_file_type.id(), H5T_DIR_DEFAULT);
    hpc::h5::datatype in_mem_type(in_mem_type_id);
    size_t in_struct_size = in_mem_type.size();

    // We need to find the offset of "Descendant" field in the NATIVE type
    size_t descendant_offset = 0;
    bool found_descendant = false;
    for (unsigned i = 0; i < in_mem_type.n_members(); ++i) {
      if (in_mem_type.member_name(i) == _config.tree_mapping[1]) {
        descendant_offset = in_mem_type.member_offset(i);
        found_descendant = true;
        break;
      }
    }

    if (!found_descendant) {
      throw std::runtime_error("Descendant field '" + _config.tree_mapping[1] +
                               "' not found in input file");
    }

    // Output memory type
    size_t out_struct_size =
        in_struct_size + 4 * sizeof(long long) + sizeof(int);
    hid_t out_mem_type_id = H5Tcreate(H5T_COMPOUND, out_struct_size);
    // Copy members from input type (using native type)
    for (unsigned i = 0; i < in_mem_type.n_members(); ++i) {
      std::string name = in_mem_type.member_name(i);
      size_t offset = in_mem_type.member_offset(i);
      hid_t member_type = H5Tget_member_type(in_mem_type.id(), i);
      H5Tinsert(out_mem_type_id, name.c_str(), offset, member_type);
      H5Tclose(member_type);
    }
    H5Tinsert(out_mem_type_id, "globaltreeid", in_struct_size,
              H5T_NATIVE_LLONG);
    H5Tinsert(out_mem_type_id, "breadthfirst_traversalorder",
              in_struct_size + sizeof(long long), H5T_NATIVE_LLONG);
    H5Tinsert(out_mem_type_id, "depthfirst_traversalorder",
              in_struct_size + 2 * sizeof(long long), H5T_NATIVE_LLONG);
    H5Tinsert(out_mem_type_id, "subtree_count",
              in_struct_size + 3 * sizeof(long long), H5T_NATIVE_LLONG);
    H5Tinsert(out_mem_type_id, "localgalaxyid",
              in_struct_size + 4 * sizeof(long long), H5T_NATIVE_INT);

    hpc::h5::datatype out_mem_type(out_mem_type_id);

    out_ds.create(out_file, "galaxies", out_mem_type,
                  (hsize_t)tree_offsets[n_trees]);

    // Buffer for reading one tree
    std::vector<char> in_buffer;
    std::vector<char> out_buffer;

    for (size_t t = start_tree; t < end_tree; ++t) {
      long long count = tree_counts[t];
      if (count == 0)
        continue;

      long long start = tree_offsets[t];

      in_buffer.resize(count * in_struct_size);
      out_buffer.resize(count * out_struct_size);

      // Read galaxies for this tree
      hpc::h5::dataspace mem_space(count);
      hpc::h5::dataspace file_space = in_ds.dataspace();

      file_space.select_hyperslab(H5S_SELECT_SET, (hsize_t)count,
                                  (hsize_t)start);

      in_ds.read(in_buffer.data(), in_mem_type, mem_space, file_space);

      // Build tree structure
      std::vector<Node> nodes(count);
      std::vector<Node *> node_ptrs(count);

      // Map from galaxy index (in this tree) to Node*
      for (int i = 0; i < count; ++i) {
        nodes[i].id = i;
        nodes[i].original_index = start + i;

        // Read Descendant
        int descendant_global = *reinterpret_cast<int *>(
            in_buffer.data() + i * in_struct_size + descendant_offset);

        if (descendant_global != -1) {
          // Descendant in the input file is a LOCAL index (relative to the
          // start of the tree)
          long long descendant_local = descendant_global;
          if (descendant_local >= 0 && descendant_local < count) {
            nodes[i].descendant = &nodes[descendant_local];
          }
        }
        node_ptrs[i] = &nodes[i];
      }

      // Link progenitors
      for (int i = 0; i < count; ++i) {
        if (nodes[i].descendant) {
          nodes[i].descendant->progenitors.push_back(&nodes[i]);
        }
      }

      // Find roots (nodes with no descendant)
      std::vector<Node *> roots;
      for (int i = 0; i < count; ++i) {
        if (!nodes[i].descendant) {
          roots.push_back(&nodes[i]);
        }
      }

      // BFS
      long long bfs_counter = 0;
      std::deque<Node *> queue;
      for (auto *root : roots)
        queue.push_back(root);

      while (!queue.empty()) {
        Node *u = queue.front();
        queue.pop_front();
        u->bfs_idx = bfs_counter++;

        // Add progenitors to queue
        // Let's sort progenitors by original index to be deterministic
        std::sort(u->progenitors.begin(), u->progenitors.end(),
                  [](Node *a, Node *b) {
                    return a->original_index < b->original_index;
                  });

        for (auto *v : u->progenitors) {
          queue.push_back(v);
        }
      }

      // DFS
      long long dfs_counter = 0;
      // Python implementation uses a stack and pushes children in increasing
      // order, which means they are popped in decreasing order. To match this,
      // we should visit progenitors in reverse order (descending index).

      // ALSO: Python pushes all roots to the stack in increasing order, so they
      // are popped in decreasing order. So we must iterate roots in reverse for
      // DFS.
      for (auto it = roots.rbegin(); it != roots.rend(); ++it) {
        dfs_visit(*it, dfs_counter);
      }

      // Pack data into output buffer
      for (int i = 0; i < count; ++i) {
        char *in_ptr = in_buffer.data() + i * in_struct_size;
        char *out_ptr = out_buffer.data() + i * out_struct_size;

        // Copy original data
        std::memcpy(out_ptr, in_ptr, in_struct_size);

        // Append new fields
        char *extra_ptr = out_ptr + in_struct_size;
        long long *extra_ll = reinterpret_cast<long long *>(extra_ptr);
        extra_ll[0] = start + i; // globaltreeid (Actually Global Galaxy Index
                                 // based on Python code)
        extra_ll[1] = nodes[i].bfs_idx; // breadthfirst_traversalorder
        extra_ll[2] = nodes[i].dfs_idx; // depthfirst_traversalorder
        extra_ll[3] = nodes[i].subtree_count;

        int *extra_int =
            reinterpret_cast<int *>(extra_ptr + 4 * sizeof(long long));
        extra_int[0] = i; // localgalaxyid
      }

      // Write back
      out_ds.write(out_buffer.data(), out_mem_type, mem_space, file_space);
    }
  }
};

#define HPC_APP_CLASS SageImportApp
#include <libhpc/mpi/main.hh>
