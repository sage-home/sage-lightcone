/**
 * @file sage2kdtree.cc
 * @brief Consolidated SAGE HDF5 to KD-tree HDF5 converter
 *
 * This tool consolidates the functionality of three separate executables:
 * - sageh5toh5:  SAGE HDF5 → Depth-first ordered HDF5 (Phase 1)
 * - sageimport:  Add traversal metadata (Phase 2)
 * - dstreeinit:  Tree→Snapshot + KD-tree indexing (Phases 3 & 4)
 *
 * Workflow:
 *   SAGE HDF5 → Phase 1 → Phase 2 → Phase 3 → Phase 4 → KD-tree HDF5
 *
 * Intermediate files are written for debugging purposes.
 */

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <map>
#include <set>
#include <list>
#include <deque>
#include <algorithm>
#include <fstream>
#include <boost/range/algorithm/fill.hpp>
#include <pugixml.hpp>
#include <libhpc/system/type_traits.hh>
#include <libhpc/libhpc.hh>
#include <libtao/base/utils.hh>
#include <libtao/base/sage.hh>
#include <libtao/base/xml_dict.hh>
#include <libtao/base/data_dict.hh>
#include <libtao/base/batch.hh>

using namespace ::tao;
using namespace ::hpc;
using namespace ::sage;

static unsigned long long const chunk_size = 10000;

//==============================================================================
// Supporting Structures
//==============================================================================

// From sageimport.cc - Tree graph node for BFS/DFS traversal
struct Node {
    size_t id;
    long long original_index;
    Node* descendant = nullptr;
    std::vector<Node*> progenitors;

    long long bfs_idx = 0;
    long long dfs_idx = 0;
    long long subtree_count = 0;
};

// From dstreeinit.cc - Functor for KD-tree spatial reordering
struct data_permuter {
    data_permuter() : crds(nullptr), tree_idxs(nullptr), subsize(nullptr), idxs(nullptr) {
    }

    data_permuter(std::array<std::vector<double>, 3> &crds,
                  std::vector<unsigned long long> &   tree_idxs,
                  std::vector<unsigned> &             subsize,
                  std::vector<unsigned long long> &   idxs) :
        crds(&crds),
        tree_idxs(&tree_idxs), subsize(&subsize), idxs(&idxs) {
    }

    void operator()(hpc::mpi::balanced_partition const &part) {
        for(unsigned ii = 0; ii < 3; ++ii)
            part.transfer((*crds)[ii]);
        part.transfer(*tree_idxs);
        part.transfer(*subsize);
        part.transfer(*idxs);
    }

    std::array<std::vector<double>, 3> *crds;
    std::vector<unsigned long long> *   tree_idxs;
    std::vector<unsigned> *             subsize;
    std::vector<unsigned long long> *   idxs;
};

// From dstreeinit.cc - SED data structure
struct sed_data_t {
    int    descendant;
    int    snapshot;
    int    local_index;
    int    merge_type;
    double dt;
    double disk_sfr;
    double bulge_sfr;
    double disk_sfr_z;
    double bulge_sfr_z;
};

// From sageimport.cc - Field metadata
struct SageField {
    std::string name;
    std::string type;
    std::string label;
    std::string description;
    std::string units;
    std::string group;
    int order;
};

//==============================================================================
// Exception class for pipeline errors
//==============================================================================

class PipelineException : public std::runtime_error {
public:
    PipelineException(int phase, const std::string& msg)
        : std::runtime_error("Phase " + std::to_string(phase) + ": " + msg),
          _phase(phase) {}
    int phase() const { return _phase; }
private:
    int _phase;
};

//==============================================================================
// Main Application Class
//==============================================================================

class sage2kdtree_application : public hpc::mpi::application {
public:
    sage2kdtree_application(int argc, char *argv[]);
    void operator()();

private:
    //==========================================================================
    // Phase 1: SAGE HDF5 → Depth-First Ordered (from sageh5toh5.cc)
    //==========================================================================
    void phase1_sageh5_to_depthfirst();
    void _load_param(hpc::fs::path const &fn);
    void _load_redshifts(hpc::fs::path const &fn);
    void _depthfirst_ordering(hpc::view<std::vector<sage::galaxy>> gals);
    unsigned _depthfirst_recurse(unsigned idx,
                                 hpc::view<std::vector<sage::galaxy>> &gals,
                                 std::vector<unsigned> &map,
                                 unsigned &dfi,
                                 std::multimap<unsigned, unsigned> const &parents);

    //==========================================================================
    // Phase 2: Add Traversal Metadata (from sageimport.cc)
    //==========================================================================
    void phase2_add_traversal_metadata();
    void _dfs_visit(Node* u, long long& counter);
    void _create_import_settings_sidecar();

    //==========================================================================
    // Phase 3: Tree → Snapshot (from dstreeinit.cc tree2sage)
    //==========================================================================
    void phase3_tree_to_snapshot();
    std::vector<unsigned long long> _count_galaxies_by_snapshot(hpc::h5::file &file,
                                                                 int nfields,
                                                                 int snapnum_index,
                                                                 std::vector<tao::data_dict_field> const &fields);

    //==========================================================================
    // Phase 4: Build KD-Tree Index (from dstreeinit.cc init)
    //==========================================================================
    void phase4_build_kdtree_index();
    void _process_snapshot(hpc::h5::file &file,
                          std::string const &name,
                          hpc::h5::file &out_file,
                          hpc::h5::dataset &data,
                          unsigned long long &displ);
    void _read_coords(hpc::h5::dataset const &dset,
                     std::array<std::vector<double>, 3> &crds,
                     std::vector<unsigned long long> &tree_idxs,
                     std::vector<unsigned> &subsize);
    void _write_attributes(hpc::h5::dataset &in_ds,
                          hpc::h5::file &out_file,
                          std::vector<unsigned long long> const &idxs,
                          unsigned long long displ);
    void _write_kdtree(hpc::h5::file &file,
                      std::string const &snap_name,
                      hpc::kdtree<> const &kdt,
                      hpc::binary_partitioner<std::array<std::vector<double>, 3>::iterator, data_permuter> const &bp,
                      std::array<std::vector<double>, 3> const &crds,
                      unsigned long long displ);
    void _write_empty_kdtree(hpc::h5::file &file, int snap_num);
    void _write_sed_data(hpc::h5::file &out_file, hpc::fs::path const &tree_fn);

    //==========================================================================
    // Utility Methods
    //==========================================================================
    void _validate_inputs();         // Validate all inputs before processing (Priority 1.1)
    void _setup_intermediate_paths();
    void _make_output_hdf5_types();  // Create HDF5 types for output (big-endian, correct precision)

    //==========================================================================
    // Member Variables
    //==========================================================================

    // Command-line arguments
    hpc::fs::path _sage_dir;          // Input: SAGE HDF5 directory
    hpc::fs::path _param_fn;          // Input: SAGE parameter file
    hpc::fs::path _alist_fn;          // Input: expansion factor list
    hpc::fs::path _output_fn;         // Output: final KD-tree HDF5

    // Intermediate file paths (for debugging)
    hpc::fs::path _depthfirst_fn;     // Phase 1 output
    hpc::fs::path _enhanced_fn;       // Phase 2 output
    hpc::fs::path _bysnap_fn;         // Phase 3 output

    unsigned _ppc;                    // Particles per cell for KD-tree
    int _verb;                        // Verbosity level

    // Cosmology parameters (loaded in Phase 1)
    double _box_size;
    double _hubble;
    double _omega_l;
    double _omega_m;
    std::vector<double> _redshifts;   // Redshift per snapshot

    // HDF5 datatypes (created in Phases 1 & 4)
    hpc::h5::datatype _mem_type;
    hpc::h5::datatype _file_type;
    hpc::h5::datatype _lc_mem_type;
    hpc::h5::datatype _lc_file_type;

    // Field metadata (used in Phases 2-4)
    std::vector<SageField> _fields;

    // MPI communicator
    hpc::mpi::comm *_comm;
};

//==============================================================================
// Constructor: Parse command-line arguments
//==============================================================================

sage2kdtree_application::sage2kdtree_application(int argc, char *argv[])
    : hpc::mpi::application(argc, argv),
      _comm(const_cast<hpc::mpi::comm*>(&hpc::mpi::comm::world)),
      _box_size(0), _hubble(0), _omega_l(0), _omega_m(0), _ppc(1000), _verb(1) {

    // Setup command-line options
    options().add_options()
        ("sage,s", hpc::po::value<hpc::fs::path>(&_sage_dir)->required(),
         "SAGE HDF5 output directory")
        ("param,p", hpc::po::value<hpc::fs::path>(&_param_fn)->required(),
         "SAGE parameter file")
        ("alist,a", hpc::po::value<hpc::fs::path>(&_alist_fn)->required(),
         "SAGE expansion factor list file")
        ("output,o", hpc::po::value<hpc::fs::path>(&_output_fn)->required(),
         "Output KD-tree HDF5 file")
        ("ppc", hpc::po::value<unsigned>(&_ppc)->default_value(1000),
         "Particles per cell for KD-tree")
        ("verbose,v", hpc::po::value<int>(&_verb)->default_value(1),
         "Verbosity level (0=quiet, 1=progress, 2=info, 3=debug)");

    // Allow positional arguments
    positional_options().add("sage", 1);
    positional_options().add("param", 2);
    positional_options().add("alist", 3);
    positional_options().add("output", 4);

    // Parse options
    parse_options(argc, argv);

    // Validate required arguments
    EXCEPT(!_sage_dir.empty(), "No SAGE output directory given.");
    EXCEPT(!_param_fn.empty(), "No SAGE parameter file given.");
    EXCEPT(!_alist_fn.empty(), "No SAGE expansion file given.");
    EXCEPT(!_output_fn.empty(), "No output file given.");

    // Setup logging based on verbosity
    if(_verb > 0) {
        hpc::log::levels_type lvl;
        if(_verb == 1) lvl = hpc::log::info;
        else if(_verb == 2) lvl = hpc::log::debug;
        else if(_verb == 3) lvl = hpc::log::trivial;

        if(_comm->size() > 1)
            LOG_PUSH(new hpc::mpi::logger("sage2kdtree.log.", lvl));
        else
            LOG_PUSH(new hpc::log::stdout(lvl));
    }
}

//==============================================================================
// Main Operator: Execute four-phase pipeline
//==============================================================================

void sage2kdtree_application::operator()() {
    try {
        // Validate inputs before processing
        if (_verb && _comm->rank() == 0) {
            std::cout << "=== Validating Inputs ===" << std::endl;
        }
        _validate_inputs();
        if (_verb && _comm->rank() == 0) {
            std::cout << "  ✓ All inputs validated successfully\n" << std::endl;
        }

        // Setup intermediate file paths
        _setup_intermediate_paths();

        // Phase 1: SAGE HDF5 → Depth-First Ordered
        if (_verb && _comm->rank() == 0) {
            std::cout << "=== Phase 1: SAGE HDF5 → Depth-first ordered ===" << std::endl;
        }
        phase1_sageh5_to_depthfirst();
        _comm->barrier();
        if (_verb && _comm->rank() == 0) {
            std::cout << "  → Written: " << _depthfirst_fn << std::endl;
        }

        // Phase 2: Add Traversal Metadata
        if (_verb && _comm->rank() == 0) {
            std::cout << "\n=== Phase 2: Adding traversal metadata ===" << std::endl;
        }
        phase2_add_traversal_metadata();
        _comm->barrier();
        if (_verb && _comm->rank() == 0) {
            std::cout << "  → Written: " << _enhanced_fn << std::endl;
        }

        // Phase 3: Tree Order → Snapshot Order
        if (_verb && _comm->rank() == 0) {
            std::cout << "\n=== Phase 3: Tree order → Snapshot order ===" << std::endl;
        }
        phase3_tree_to_snapshot();
        _comm->barrier();
        if (_verb && _comm->rank() == 0) {
            std::cout << "  → Written: " << _bysnap_fn << std::endl;
        }

        // Phase 4: Build KD-Tree Index
        if (_verb && _comm->rank() == 0) {
            std::cout << "\n=== Phase 4: Building KD-tree spatial index ===" << std::endl;
        }
        phase4_build_kdtree_index();
        _comm->barrier();

        // Success summary
        if (_verb && _comm->rank() == 0) {
            std::cout << "\n=== Conversion Complete ===" << std::endl;
            std::cout << "Final output: " << _output_fn << std::endl;
            std::cout << "\nIntermediate files (for debugging):" << std::endl;
            std::cout << "  - " << _depthfirst_fn << std::endl;
            std::cout << "  - " << _enhanced_fn << std::endl;
            std::cout << "  - " << _bysnap_fn << std::endl;
        }

    } catch (PipelineException &e) {
        if (_comm->rank() == 0) {
            std::cerr << "\n=== PIPELINE FAILED ===" << std::endl;
            std::cerr << "Error: " << e.what() << std::endl;
            std::cerr << "\nIntermediate files preserved for debugging:" << std::endl;
            if (e.phase() >= 1 && hpc::fs::exists(_depthfirst_fn))
                std::cerr << "  Phase 1 output: " << _depthfirst_fn << std::endl;
            if (e.phase() >= 2 && hpc::fs::exists(_enhanced_fn))
                std::cerr << "  Phase 2 output: " << _enhanced_fn << std::endl;
            if (e.phase() >= 3 && hpc::fs::exists(_bysnap_fn))
                std::cerr << "  Phase 3 output: " << _bysnap_fn << std::endl;
        }
        _comm->abort(1);
    } catch (std::exception &e) {
        if (_comm->rank() == 0) {
            std::cerr << "\n=== UNEXPECTED ERROR ===" << std::endl;
            std::cerr << "Error: " << e.what() << std::endl;
        }
        _comm->abort(1);
    }
}

//==============================================================================
// Utility: Setup intermediate file paths
//==============================================================================

void sage2kdtree_application::_setup_intermediate_paths() {
    hpc::fs::path output_dir = _output_fn.parent_path();
    std::string base_stem = _output_fn.stem().string();

    // Remove -kdtree suffix if present
    if (base_stem.size() > 7 &&
        base_stem.substr(base_stem.size() - 7) == "-kdtree") {
        base_stem = base_stem.substr(0, base_stem.size() - 7);
    }

    _depthfirst_fn = output_dir / (base_stem + "-depthfirstordered.h5");
    _enhanced_fn = output_dir / ("out_" + base_stem + "-depthfirstordered.h5");
    _bysnap_fn = output_dir / (base_stem + "-bysnap.h5");
}

//==============================================================================
// Validate Inputs (Priority 1.1)
// Comprehensive input validation with clear error messages
//==============================================================================

void sage2kdtree_application::_validate_inputs() {
    // Only rank 0 performs validation to avoid redundant checks
    if (_comm->rank() != 0) return;

    std::vector<std::string> errors;
    std::vector<std::string> warnings;

    // 1. Check file existence
    if (_verb >= 2) std::cout << "  → Checking input file existence..." << std::endl;

    if (!hpc::fs::exists(_param_fn)) {
        errors.push_back("Parameter file not found: " + _param_fn.string() +
                        "\n      Suggestion: Check the path to your SAGE parameter file (.par)");
    }

    if (!hpc::fs::exists(_alist_fn)) {
        errors.push_back("Expansion factor list not found: " + _alist_fn.string() +
                        "\n      Suggestion: This file should contain scale factors (a_list) from your simulation");
    }

    if (!hpc::fs::exists(_sage_dir)) {
        errors.push_back("SAGE output directory not found: " + _sage_dir.string() +
                        "\n      Suggestion: Check the path to your SAGE HDF5 output directory");
    } else if (!hpc::fs::is_directory(_sage_dir)) {
        errors.push_back("SAGE output path is not a directory: " + _sage_dir.string());
    }

    // Early exit if files don't exist
    if (!errors.empty()) {
        std::string msg = "\n=== INPUT VALIDATION FAILED ===\n";
        for (const auto& err : errors) {
            msg += "  ✗ " + err + "\n";
        }
        throw PipelineException(0, msg);
    }

    // 2. Validate SAGE HDF5 directory has .hdf5 or .h5 files
    if (_verb >= 2) std::cout << "  → Checking for SAGE HDF5 files..." << std::endl;

    bool found_hdf5_files = false;
    hpc::fs::directory_iterator end_iter;
    for(hpc::fs::directory_iterator dir_itr(_sage_dir); dir_itr != end_iter; ++dir_itr) {
        if(hpc::fs::is_regular_file(dir_itr->status())) {
            std::string ext = dir_itr->path().extension().string();
            if (ext == ".hdf5" || ext == ".h5") {
                found_hdf5_files = true;
                break;
            }
        }
    }

    if (!found_hdf5_files) {
        errors.push_back("No HDF5 files found in SAGE directory: " + _sage_dir.string() +
                        "\n      Suggestion: SAGE HDF5 output should have .hdf5 or .h5 extension");
    }

    // 3. Validate first SAGE HDF5 file structure
    if (_verb >= 2) std::cout << "  → Validating SAGE HDF5 structure..." << std::endl;

    // Find first substantial HDF5 file (skip small master files)
    // SAGE typically creates model_0.hdf5, model_1.hdf5, etc. (data files)
    // and model.hdf5 (small master file with just header)
    hpc::fs::path first_hdf5;
    const size_t min_file_size = 100000; // 100KB - data files are typically much larger

    for(hpc::fs::directory_iterator dir_itr(_sage_dir); dir_itr != end_iter; ++dir_itr) {
        if(hpc::fs::is_regular_file(dir_itr->status())) {
            std::string ext = dir_itr->path().extension().string();
            std::string stem = dir_itr->path().stem().string();

            // Look for .hdf5 or .h5 files
            if (ext == ".hdf5" || ext == ".h5") {
                // Skip small files (likely master files without data)
                size_t file_size = hpc::fs::file_size(dir_itr->path());
                if (file_size < min_file_size) {
                    if (_verb >= 2) {
                        std::cout << "    Skipping small file (likely master): "
                                  << dir_itr->path().filename() << " (" << file_size << " bytes)" << std::endl;
                    }
                    continue;
                }

                // Prefer files matching model_N.hdf5 pattern (actual data files)
                if (stem.find("model_") == 0 && stem.length() > 6) {
                    first_hdf5 = dir_itr->path();
                    break;
                }

                // Otherwise accept any substantial file
                if (first_hdf5.empty()) {
                    first_hdf5 = dir_itr->path();
                }
            }
        }
    }

    if (!first_hdf5.empty()) {
        if (_verb >= 2) {
            std::cout << "    Validating file: " << first_hdf5.filename() << std::endl;
        }

        try {
            hpc::h5::file test_file(first_hdf5.string(), H5F_ACC_RDONLY);

            // Check for Header/Simulation group
            bool has_header_simulation = test_file.has_link("Header/Simulation");
            if (!has_header_simulation) {
                warnings.push_back("Header/Simulation group not found in " + first_hdf5.filename().string() +
                                 "\n      Cosmology parameters may not be available");
            } else {
                // Validate cosmology attributes
                auto sim_group = test_file.group("Header/Simulation");
                std::vector<std::string> required_attrs = {"BoxSize", "Hubble_h", "Omega_m", "Omega_lambda"};
                std::vector<std::string> missing_attrs;

                for (const auto& attr : required_attrs) {
                    if (H5Aexists(sim_group.id(), attr.c_str()) <= 0) {
                        missing_attrs.push_back(attr);
                    }
                }

                if (!missing_attrs.empty()) {
                    std::string missing_list;
                    for (const auto& attr : missing_attrs) {
                        missing_list += attr + " ";
                    }
                    warnings.push_back("Missing cosmology attributes in Header/Simulation: " + missing_list +
                                     "\n      These will be read from parameter file instead");
                }
            }

            // Check for at least one Snap_N group
            int snap_count = 0;
            while (test_file.has_link("Snap_" + std::to_string(snap_count))) {
                snap_count++;
            }

            if (snap_count == 0) {
                errors.push_back("No Snap_N groups found in " + first_hdf5.filename().string() +
                                "\n      Suggestion: SAGE HDF5 output should contain Snap_0, Snap_1, etc.");
            } else {
                if (_verb >= 2) std::cout << "    Found " << snap_count << " snapshots" << std::endl;

                // Validate required fields in Snap_0
                auto snap0 = test_file.group("Snap_0");
                std::vector<std::string> required_fields = {"Posx", "Posy", "Posz", "SAGETreeIndex"};
                std::vector<std::string> missing_fields;

                for (const auto& field : required_fields) {
                    if (!snap0.has_link(field)) {
                        missing_fields.push_back(field);
                    }
                }

                if (!missing_fields.empty()) {
                    std::string missing_list;
                    for (const auto& field : missing_fields) {
                        missing_list += field + " ";
                    }
                    errors.push_back("Missing required fields in Snap_0: " + missing_list +
                                    "\n      Suggestion: SAGE HDF5 output must contain position (Posx/y/z) and tree index fields");
                }
            }

        } catch (std::exception& e) {
            errors.push_back("Failed to open/validate SAGE HDF5 file " + first_hdf5.filename().string() +
                            ": " + std::string(e.what()) +
                            "\n      Suggestion: Ensure SAGE output is in HDF5 format (OutputFormat sage_hdf5)");
        }
    }

    // 4. Validate parameter file has required keys
    if (_verb >= 2) std::cout << "  → Validating parameter file..." << std::endl;

    std::ifstream param_file(_param_fn.native());
    if (param_file.good()) {
        std::set<std::string> found_params;
        std::string line;
        while(std::getline(param_file, line)) {
            if(line.empty() || line[0] == '%') continue;
            std::stringstream ss(line);
            std::string key;
            ss >> key;
            found_params.insert(key);
        }

        std::vector<std::string> required_params = {"BoxSize", "Hubble_h"};
        std::vector<std::string> recommended_params = {"Omega_Lambda", "OmegaLambda", "Omega_m", "Omega"};

        for (const auto& param : required_params) {
            if (found_params.find(param) == found_params.end()) {
                warnings.push_back("Parameter '" + param + "' not found in " + _param_fn.filename().string());
            }
        }

        bool has_omega_lambda = (found_params.find("Omega_Lambda") != found_params.end() ||
                                 found_params.find("OmegaLambda") != found_params.end());
        bool has_omega_m = (found_params.find("Omega_m") != found_params.end() ||
                           found_params.find("Omega") != found_params.end());

        if (!has_omega_lambda || !has_omega_m) {
            warnings.push_back("Cosmology parameters (Omega_Lambda, Omega_m) incomplete in parameter file");
        }
    }

    // 5. Validate expansion factor list
    if (_verb >= 2) std::cout << "  → Validating expansion factor list..." << std::endl;

    std::ifstream alist_file(_alist_fn.native());
    if (alist_file.good()) {
        int line_count = 0;
        std::string line;
        while(std::getline(alist_file, line)) {
            if(line.empty() || line[0] == '%') continue;
            line_count++;
        }

        if (line_count == 0) {
            errors.push_back("Expansion factor list is empty: " + _alist_fn.filename().string() +
                            "\n      Suggestion: File should contain one scale factor per snapshot");
        } else {
            if (_verb >= 2) std::cout << "    Found " << line_count << " scale factors" << std::endl;
        }
    }

    // Report results
    if (!warnings.empty() && _verb >= 1) {
        std::cout << "\n  Warnings:" << std::endl;
        for (const auto& warn : warnings) {
            std::cout << "  ⚠  " << warn << std::endl;
        }
        std::cout << std::endl;
    }

    if (!errors.empty()) {
        std::string msg = "\n=== INPUT VALIDATION FAILED ===\n";
        for (const auto& err : errors) {
            msg += "  ✗ " + err + "\n";
        }
        throw PipelineException(0, msg);
    }
}

//==============================================================================
// Create HDF5 Types for Output
// Uses native endian for both memory and file (simpler, more efficient)
// Changed from big-endian to native endian to avoid conversion overhead
//==============================================================================

void sage2kdtree_application::_make_output_hdf5_types() {
    using namespace hpc::h5;

    // Create memory type (native byte order for in-memory representation)
    _mem_type.compound(sizeof(sage::galaxy));
    _mem_type.insert(datatype::native_int,   "snapnum",              offsetof(sage::galaxy, snapshot));
    _mem_type.insert(datatype::native_int,   "type",                  offsetof(sage::galaxy, type));
    _mem_type.insert(datatype::native_llong, "galaxy_index",            offsetof(sage::galaxy, galaxy_idx));
    _mem_type.insert(datatype::native_llong, "central_galaxy_index",    offsetof(sage::galaxy, central_galaxy_idx));
    _mem_type.insert(datatype::native_int,   "sage_halo_index",         offsetof(sage::galaxy, sage_halo_idx));
    _mem_type.insert(datatype::native_int,   "sage_tree_index",         offsetof(sage::galaxy, sage_tree_idx));
    _mem_type.insert(datatype::native_llong, "simulation_halo_index",   offsetof(sage::galaxy, simulation_halo_idx));
    _mem_type.insert(datatype::native_int,   "local_index",           offsetof(sage::galaxy, local_index));
    _mem_type.insert(datatype::native_llong, "global_index",          offsetof(sage::galaxy, global_index));
    _mem_type.insert(datatype::native_int,   "descendant",            offsetof(sage::galaxy, descendant));
    _mem_type.insert(datatype::native_llong, "global_descendant",     offsetof(sage::galaxy, global_descendant));
    _mem_type.insert(datatype::native_int,   "subsize",               offsetof(sage::galaxy, subsize));
    _mem_type.insert(datatype::native_int,   "merge_type",            offsetof(sage::galaxy, merge_type));
    _mem_type.insert(datatype::native_int,   "merge_into_id",         offsetof(sage::galaxy, merge_into_id));
    _mem_type.insert(datatype::native_int,   "merge_into_snapshot",   offsetof(sage::galaxy, merge_into_snapshot));
    _mem_type.insert(datatype::native_float, "dt",                    offsetof(sage::galaxy, dt));
    _mem_type.insert(datatype::native_float, "posx",                offsetof(sage::galaxy, pos[0]));
    _mem_type.insert(datatype::native_float, "posy",                offsetof(sage::galaxy, pos[1]));
    _mem_type.insert(datatype::native_float, "posz",                offsetof(sage::galaxy, pos[2]));
    _mem_type.insert(datatype::native_float, "velx",                offsetof(sage::galaxy, vel[0]));
    _mem_type.insert(datatype::native_float, "vely",                offsetof(sage::galaxy, vel[1]));
    _mem_type.insert(datatype::native_float, "velz",                offsetof(sage::galaxy, vel[2]));
    _mem_type.insert(datatype::native_float, "spin_x",               offsetof(sage::galaxy, spin[0]));
    _mem_type.insert(datatype::native_float, "spin_y",               offsetof(sage::galaxy, spin[1]));
    _mem_type.insert(datatype::native_float, "spin_z",               offsetof(sage::galaxy, spin[2]));
    _mem_type.insert(datatype::native_int,   "n_darkmatter_particles",         offsetof(sage::galaxy, num_particles));
    _mem_type.insert(datatype::native_float, "virial_mass",                  offsetof(sage::galaxy, mvir));
    _mem_type.insert(datatype::native_float, "central_galaxy_mvir",          offsetof(sage::galaxy, central_mvir));
    _mem_type.insert(datatype::native_float, "virial_radius",                  offsetof(sage::galaxy, rvir));
    _mem_type.insert(datatype::native_float, "virial_velocity",                  offsetof(sage::galaxy, vvir));
    _mem_type.insert(datatype::native_float, "max_velocity",                  offsetof(sage::galaxy, vmax));
    _mem_type.insert(datatype::native_float, "velocity_dispersion",              offsetof(sage::galaxy, vel_disp));
    _mem_type.insert(datatype::native_float, "cold_gas",              offsetof(sage::galaxy, cold_gas));
    _mem_type.insert(datatype::native_float, "stellar_mass",          offsetof(sage::galaxy, stellar_mass));
    _mem_type.insert(datatype::native_float, "bulge_mass",            offsetof(sage::galaxy, bulge_mass));
    _mem_type.insert(datatype::native_float, "hot_gas",               offsetof(sage::galaxy, hot_gas));
    _mem_type.insert(datatype::native_float, "ejected_mass",          offsetof(sage::galaxy, ejected_mass));
    _mem_type.insert(datatype::native_float, "blackhole_mass",        offsetof(sage::galaxy, blackhole_mass));
    _mem_type.insert(datatype::native_float, "ics",                   offsetof(sage::galaxy, ics));
    _mem_type.insert(datatype::native_float, "metals_cold_gas",       offsetof(sage::galaxy, metals_cold_gas));
    _mem_type.insert(datatype::native_float, "metals_stellar_mass",   offsetof(sage::galaxy, metals_stellar_mass));
    _mem_type.insert(datatype::native_float, "metals_bulge_mass",     offsetof(sage::galaxy, metals_bulge_mass));
    _mem_type.insert(datatype::native_float, "metals_hot_gas",        offsetof(sage::galaxy, metals_hot_gas));
    _mem_type.insert(datatype::native_float, "metals_ejected_mass",   offsetof(sage::galaxy, metals_ejected_mass));
    _mem_type.insert(datatype::native_float, "metals_ics",            offsetof(sage::galaxy, metals_ics));
    _mem_type.insert(datatype::native_float, "sfr_disk",              offsetof(sage::galaxy, sfr_disk));
    _mem_type.insert(datatype::native_float, "sfr_bulge",             offsetof(sage::galaxy, sfr_bulge));
    _mem_type.insert(datatype::native_float, "sfr_disk_z",            offsetof(sage::galaxy, sfr_disk_z));
    _mem_type.insert(datatype::native_float, "sfr_bulge_z",           offsetof(sage::galaxy, sfr_bulge_z));
    _mem_type.insert(datatype::native_float, "disk_scale_radius",     offsetof(sage::galaxy, disk_scale_radius));
    _mem_type.insert(datatype::native_float, "cooling",               offsetof(sage::galaxy, cooling));
    _mem_type.insert(datatype::native_float, "heating",               offsetof(sage::galaxy, heating));
    _mem_type.insert(datatype::native_float, "quasar_mode_bh_accretion_mass", offsetof(sage::galaxy, quasar_mode_bh_accretion_mass));
    _mem_type.insert(datatype::native_float, "time_of_last_major_merger", offsetof(sage::galaxy, time_of_last_major_merger));
    _mem_type.insert(datatype::native_float, "time_of_last_minor_merger", offsetof(sage::galaxy, time_of_last_minor_merger));
    _mem_type.insert(datatype::native_float, "outflow_rate",          offsetof(sage::galaxy, outflow_rate));
    _mem_type.insert(datatype::native_float, "infall_mvir",           offsetof(sage::galaxy, infall_mvir));
    _mem_type.insert(datatype::native_float, "infall_vvir",           offsetof(sage::galaxy, infall_vvir));
    _mem_type.insert(datatype::native_float, "infall_vmax",           offsetof(sage::galaxy, infall_vmax));

    // Create file type (native byte order - same as memory type)
    // Simpler and more efficient than big-endian conversion
    _file_type.compound(sizeof(sage::galaxy));
    _file_type.insert(datatype::native_int,   "snapnum",              offsetof(sage::galaxy, snapshot));
    _file_type.insert(datatype::native_int,   "type",                  offsetof(sage::galaxy, type));
    _file_type.insert(datatype::native_llong, "galaxy_index",          offsetof(sage::galaxy, galaxy_idx));
    _file_type.insert(datatype::native_llong, "central_galaxy_index",  offsetof(sage::galaxy, central_galaxy_idx));
    _file_type.insert(datatype::native_int,   "sage_halo_index",       offsetof(sage::galaxy, sage_halo_idx));
    _file_type.insert(datatype::native_int,   "sage_tree_index",       offsetof(sage::galaxy, sage_tree_idx));
    _file_type.insert(datatype::native_llong, "simulation_halo_index", offsetof(sage::galaxy, simulation_halo_idx));
    _file_type.insert(datatype::native_int,   "local_index",           offsetof(sage::galaxy, local_index));
    _file_type.insert(datatype::native_llong, "global_index",          offsetof(sage::galaxy, global_index));
    _file_type.insert(datatype::native_int,   "descendant",            offsetof(sage::galaxy, descendant));
    _file_type.insert(datatype::native_llong, "global_descendant",     offsetof(sage::galaxy, global_descendant));
    _file_type.insert(datatype::native_int,   "subsize",               offsetof(sage::galaxy, subsize));
    _file_type.insert(datatype::native_int,   "merge_type",            offsetof(sage::galaxy, merge_type));
    _file_type.insert(datatype::native_int,   "merge_into_id",         offsetof(sage::galaxy, merge_into_id));
    _file_type.insert(datatype::native_int,   "merge_into_snapshot",   offsetof(sage::galaxy, merge_into_snapshot));
    _file_type.insert(datatype::native_float, "dt",                    offsetof(sage::galaxy, dt));
    _file_type.insert(datatype::native_float, "posx",                  offsetof(sage::galaxy, pos[0]));
    _file_type.insert(datatype::native_float, "posy",                  offsetof(sage::galaxy, pos[1]));
    _file_type.insert(datatype::native_float, "posz",                  offsetof(sage::galaxy, pos[2]));
    _file_type.insert(datatype::native_float, "velx",                  offsetof(sage::galaxy, vel[0]));
    _file_type.insert(datatype::native_float, "vely",                  offsetof(sage::galaxy, vel[1]));
    _file_type.insert(datatype::native_float, "velz",                  offsetof(sage::galaxy, vel[2]));
    _file_type.insert(datatype::native_float, "spin_x",                offsetof(sage::galaxy, spin[0]));
    _file_type.insert(datatype::native_float, "spin_y",                offsetof(sage::galaxy, spin[1]));
    _file_type.insert(datatype::native_float, "spin_z",                offsetof(sage::galaxy, spin[2]));
    _file_type.insert(datatype::native_int,   "n_darkmatter_particles", offsetof(sage::galaxy, num_particles));
    _file_type.insert(datatype::native_float, "virial_mass",           offsetof(sage::galaxy, mvir));
    _file_type.insert(datatype::native_float, "central_galaxy_mvir",   offsetof(sage::galaxy, central_mvir));
    _file_type.insert(datatype::native_float, "virial_radius",         offsetof(sage::galaxy, rvir));
    _file_type.insert(datatype::native_float, "virial_velocity",       offsetof(sage::galaxy, vvir));
    _file_type.insert(datatype::native_float, "max_velocity",          offsetof(sage::galaxy, vmax));
    _file_type.insert(datatype::native_float, "velocity_dispersion",   offsetof(sage::galaxy, vel_disp));
    _file_type.insert(datatype::native_float, "cold_gas",              offsetof(sage::galaxy, cold_gas));
    _file_type.insert(datatype::native_float, "stellar_mass",          offsetof(sage::galaxy, stellar_mass));
    _file_type.insert(datatype::native_float, "bulge_mass",            offsetof(sage::galaxy, bulge_mass));
    _file_type.insert(datatype::native_float, "hot_gas",               offsetof(sage::galaxy, hot_gas));
    _file_type.insert(datatype::native_float, "ejected_mass",          offsetof(sage::galaxy, ejected_mass));
    _file_type.insert(datatype::native_float, "blackhole_mass",        offsetof(sage::galaxy, blackhole_mass));
    _file_type.insert(datatype::native_float, "ics",                   offsetof(sage::galaxy, ics));
    _file_type.insert(datatype::native_float, "metals_cold_gas",       offsetof(sage::galaxy, metals_cold_gas));
    _file_type.insert(datatype::native_float, "metals_stellar_mass",   offsetof(sage::galaxy, metals_stellar_mass));
    _file_type.insert(datatype::native_float, "metals_bulge_mass",     offsetof(sage::galaxy, metals_bulge_mass));
    _file_type.insert(datatype::native_float, "metals_hot_gas",        offsetof(sage::galaxy, metals_hot_gas));
    _file_type.insert(datatype::native_float, "metals_ejected_mass",   offsetof(sage::galaxy, metals_ejected_mass));
    _file_type.insert(datatype::native_float, "metals_ics",            offsetof(sage::galaxy, metals_ics));
    _file_type.insert(datatype::native_float, "sfr_disk",              offsetof(sage::galaxy, sfr_disk));
    _file_type.insert(datatype::native_float, "sfr_bulge",             offsetof(sage::galaxy, sfr_bulge));
    _file_type.insert(datatype::native_float, "sfr_disk_z",            offsetof(sage::galaxy, sfr_disk_z));
    _file_type.insert(datatype::native_float, "sfr_bulge_z",           offsetof(sage::galaxy, sfr_bulge_z));
    _file_type.insert(datatype::native_float, "disk_scale_radius",     offsetof(sage::galaxy, disk_scale_radius));
    _file_type.insert(datatype::native_float, "cooling",               offsetof(sage::galaxy, cooling));
    _file_type.insert(datatype::native_float, "heating",               offsetof(sage::galaxy, heating));
    _file_type.insert(datatype::native_float, "quasar_mode_bh_accretion_mass", offsetof(sage::galaxy, quasar_mode_bh_accretion_mass));
    _file_type.insert(datatype::native_float, "time_of_last_major_merger", offsetof(sage::galaxy, time_of_last_major_merger));
    _file_type.insert(datatype::native_float, "time_of_last_minor_merger", offsetof(sage::galaxy, time_of_last_minor_merger));
    _file_type.insert(datatype::native_float, "outflow_rate",          offsetof(sage::galaxy, outflow_rate));
    _file_type.insert(datatype::native_float, "infall_mvir",           offsetof(sage::galaxy, infall_mvir));
    _file_type.insert(datatype::native_float, "infall_vvir",           offsetof(sage::galaxy, infall_vvir));
    _file_type.insert(datatype::native_float, "infall_vmax",           offsetof(sage::galaxy, infall_vmax));
}

//==============================================================================
// Phase 1: SAGE HDF5 → Depth-First Ordered
// Based on sageh5toh5.cc lines 59-601
//==============================================================================

void sage2kdtree_application::phase1_sageh5_to_depthfirst() {
    _load_param(_param_fn);
    _load_redshifts(_alist_fn);

    // 1. Identify input files
    if (_verb && _comm->rank() == 0) {
        std::cout << "Searching for files in " << _sage_dir << std::endl;
    }
    std::vector<hpc::fs::path> input_files;
    if(hpc::fs::is_directory(_sage_dir)) {
        hpc::fs::directory_iterator end_iter;
        for(hpc::fs::directory_iterator dir_itr(_sage_dir); dir_itr != end_iter; ++dir_itr) {
            if(hpc::fs::is_regular_file(dir_itr->status())) {
                std::string ext = dir_itr->path().extension().string();
                std::string stem = dir_itr->path().stem().string();
                if(ext == ".hdf5" && (stem.find("model_hdf5_") == 0 || stem.find("model_") == 0)) {
                    input_files.push_back(dir_itr->path());
                }
            }
        }
    } else {
        input_files.push_back(_sage_dir);
    }
    if (_verb && _comm->rank() == 0) {
        std::cout << "Found " << input_files.size() << " input files." << std::endl;
    }
    std::sort(input_files.begin(), input_files.end()); // Ensure deterministic order

    // Distribute files among ranks (round-robin)
    std::vector<hpc::fs::path> my_files;
    for(size_t i = 0; i < input_files.size(); ++i) {
        if(i % _comm->size() == _comm->rank()) {
            my_files.push_back(input_files[i]);
        }
    }

    // 2. Scan Phase: Count trees and galaxies
    if (_verb && _comm->rank() == 0) {
        std::cout << "  → Scanning files to count trees and galaxies..." << std::endl;
    }

    unsigned long long my_tot_trees = 0;
    unsigned long long my_tot_gals = 0;
    std::vector<unsigned long long> file_tree_counts;
    std::vector<unsigned long long> file_gal_counts;

    for(size_t file_idx = 0; file_idx < my_files.size(); ++file_idx) {
        const auto& fn = my_files[file_idx];

        if (_verb >= 2 && _comm->rank() == 0) {
            std::cout << "    Scanning file " << (file_idx + 1) << "/" << my_files.size()
                      << ": " << fn.filename() << std::endl;
        }

        hpc::h5::file f;
        try {
            f.open(fn.native(), H5F_ACC_RDONLY);
        } catch (std::exception& e) {
            throw PipelineException(1, "Failed to open SAGE HDF5 file '" + fn.string() + "': " +
                                    std::string(e.what()));
        }

        // Determine number of trees from Snap_63/SAGETreeIndex
        size_t n_trees = 0;
        if(f.has_link("Snap_63/SAGETreeIndex")) {
            hpc::h5::dataset dset(f, "Snap_63/SAGETreeIndex");
            std::vector<int> tree_indices;
            tree_indices.resize(dset.extent());
            dset.read(tree_indices);
            if(!tree_indices.empty()) {
                int max_tree_idx = *std::max_element(tree_indices.begin(), tree_indices.end());
                n_trees = max_tree_idx + 1;
            }
        }

        // Fallback to TreeInfo if Snap_63 didn't give us trees
        if (n_trees == 0 && f.has_link("TreeInfo/Snap_0/NumGalsPerTreePerSnap")) {
            hpc::h5::dataset dset(f, "TreeInfo/Snap_0/NumGalsPerTreePerSnap");
            n_trees = dset.dataspace().size();
        }

        std::vector<int> tree_sizes(n_trees, 0);

        // Sum galaxies for each tree across all snapshots
        for(int snap = 0; snap < 64; ++snap) { // 64 snapshots
            std::string group_name = "Snap_" + std::to_string(snap);
            if(!f.has_link(group_name)) continue;

            hpc::h5::group g;
            f.open_group(group_name, g);

            if(g.has_link("SAGETreeIndex")) {
                std::vector<int> tree_indices;
                hpc::h5::dataset dset = g.dataset("SAGETreeIndex");
                tree_indices.resize(dset.extent());
                dset.read(tree_indices);

                for(int tree_idx : tree_indices) {
                    if(tree_idx >= 0 && (size_t)tree_idx < n_trees) {
                        tree_sizes[tree_idx]++;
                    }
                }
            }
        }

        unsigned long long n_gals = 0;
        for(int s : tree_sizes) n_gals += s;

        my_tot_trees += n_trees;
        my_tot_gals += n_gals;
        file_tree_counts.push_back(n_trees);
        file_gal_counts.push_back(n_gals);
    }

    // Global offsets via MPI scan
    unsigned long long global_tree_offset = _comm->scan(my_tot_trees);
    unsigned long long global_gal_offset = _comm->scan(my_tot_gals);
    unsigned long long total_trees = _comm->all_reduce(my_tot_trees);
    unsigned long long total_gals = _comm->all_reduce(my_tot_gals);

    if (_verb && _comm->rank() == 0) {
        std::cout << "  → Found " << total_trees << " trees with "
                  << total_gals << " galaxies" << std::endl;
    }
    LOGILN("Total Trees: ", total_trees);
    LOGILN("Total Galaxies: ", total_gals);

    // 3. Create Output File
    _make_output_hdf5_types();  // Use custom types with big-endian and correct precision
    sage::make_hdf5_sidecar("Millennium", _depthfirst_fn.native(),
                           std::set<double>(_redshifts.begin(), _redshifts.end()),
                           _hubble, _box_size);

    hpc::h5::file out_file(_depthfirst_fn.native(), H5F_ACC_TRUNC, *_comm);
    hpc::h5::dataset gals_dset(out_file, "galaxies", _file_type, hpc::h5::dataspace(total_gals));
    hpc::h5::dataset tree_displs_dset(out_file, "tree_displs", hpc::h5::datatype::native_ullong,
                                     hpc::h5::dataspace(total_trees + 1));
    hpc::h5::dataset tree_cnts_dset(out_file, "tree_counts", hpc::h5::datatype::native_uint,
                                   hpc::h5::dataspace(total_trees));

    // Buffers
    hpc::h5::buffer<sage::galaxy> gals_buf;
    gals_buf.create(gals_dset, _mem_type, hpc::h5::buffer_default_size, global_gal_offset);
    hpc::h5::buffer<unsigned long long> tree_displs_buf;
    tree_displs_buf.create(tree_displs_dset, hpc::h5::datatype::native_ullong,
                          hpc::h5::buffer_default_size, global_tree_offset);
    hpc::h5::buffer<unsigned> tree_cnts_buf;
    tree_cnts_buf.create(tree_cnts_dset, hpc::h5::datatype::native_uint,
                        hpc::h5::buffer_default_size, global_tree_offset);

    // 4. Process Files
    if (_verb && _comm->rank() == 0) {
        std::cout << "  → Processing " << my_files.size() << " files and applying depth-first ordering..." << std::endl;
    }

    unsigned long long current_gal_global_idx = global_gal_offset;
    unsigned long long current_tree_displ = global_gal_offset;

    for(size_t fidx = 0; fidx < my_files.size(); ++fidx) {
        hpc::fs::path fn = my_files[fidx];

        if (_verb >= 2 && _comm->rank() == 0) {
            std::cout << "    Processing file " << (fidx + 1) << "/" << my_files.size()
                      << ": " << fn.filename() << " (" << file_tree_counts[fidx]
                      << " trees, " << file_gal_counts[fidx] << " galaxies)" << std::endl;
        }
        LOGILN("Processing file: ", fn);

        hpc::h5::file f(fn.native(), H5F_ACC_RDONLY);
        size_t n_trees = file_tree_counts[fidx];

        // Read all data into memory, organized by tree
        std::vector<std::vector<sage::galaxy>> trees(n_trees);

        for(int snap = 0; snap < 64; ++snap) {
            std::string group_name = "Snap_" + std::to_string(snap);
            if(!f.has_link(group_name)) continue;

            hpc::h5::group g;
            f.open_group(group_name, g);

            // Read SAGETreeIndex to know where each galaxy goes
            std::vector<int> tree_indices;
            if(g.has_link("SAGETreeIndex")) {
                hpc::h5::dataset dset = g.dataset("SAGETreeIndex");
                tree_indices.resize(dset.extent());
                dset.read(tree_indices);
            }
            size_t n_gals_snap = tree_indices.size();
            if(n_gals_snap == 0) continue;

            // Helper to read column
            auto read_col = [&](const char* name, auto& vec) {
                vec.resize(n_gals_snap);
                if(g.has_link(name)) {
                    hpc::h5::dataset dset = g.dataset(name);
                    dset.read(vec);
                } else {
                    std::fill(vec.begin(), vec.end(), 0); // Handle missing columns with 0
                }
            };

            std::vector<float> posx, posy, posz, velx, vely, velz, spinx, spiny, spinz;
            std::vector<float> mvir, central_mvir, rvir, vvir, vmax, vel_disp;
            std::vector<float> cold_gas, stellar_mass, bulge_mass, hot_gas, ejected_mass, blackhole_mass, ics;
            std::vector<float> metals_cold_gas, metals_stellar_mass, metals_bulge_mass, metals_hot_gas, metals_ejected_mass, metals_ics;
            std::vector<float> sfr_disk, sfr_bulge, sfr_disk_z, sfr_bulge_z;
            std::vector<float> disk_radius, cooling, heating, quasar_acc, last_major, last_minor, outflow;
            std::vector<float> infall_mvir, infall_vmax, infall_vvir;
            std::vector<long long> gal_idx, central_gal_idx, sim_halo_idx;
            std::vector<int> sage_halo_idx, type, merge_type, merge_into_id, merge_into_snap;
            std::vector<float> dt;
            std::vector<int> len;

            read_col("Posx", posx); read_col("Posy", posy); read_col("Posz", posz);
            read_col("Velx", velx); read_col("Vely", vely); read_col("Velz", velz);
            read_col("Spinx", spinx); read_col("Spiny", spiny); read_col("Spinz", spinz);
            read_col("Mvir", mvir); read_col("CentralMvir", central_mvir);
            read_col("Rvir", rvir); read_col("Vvir", vvir); read_col("Vmax", vmax); read_col("VelDisp", vel_disp);
            read_col("ColdGas", cold_gas); read_col("StellarMass", stellar_mass); read_col("BulgeMass", bulge_mass);
            read_col("HotGas", hot_gas); read_col("EjectedMass", ejected_mass); read_col("BlackHoleMass", blackhole_mass); read_col("IntraClusterStars", ics);
            read_col("MetalsColdGas", metals_cold_gas); read_col("MetalsStellarMass", metals_stellar_mass); read_col("MetalsBulgeMass", metals_bulge_mass);
            read_col("MetalsHotGas", metals_hot_gas); read_col("MetalsEjectedMass", metals_ejected_mass); read_col("MetalsIntraClusterStars", metals_ics);
            read_col("SfrDisk", sfr_disk); read_col("SfrBulge", sfr_bulge); read_col("SfrDiskZ", sfr_disk_z); read_col("SfrBulgeZ", sfr_bulge_z);
            read_col("DiskRadius", disk_radius); read_col("Cooling", cooling); read_col("Heating", heating);
            read_col("QuasarModeBHaccretionMass", quasar_acc); read_col("TimeOfLastMajorMerger", last_major); read_col("TimeOfLastMinorMerger", last_minor);
            read_col("OutflowRate", outflow);
            read_col("infallMvir", infall_mvir); read_col("infallVmax", infall_vmax); read_col("infallVvir", infall_vvir);
            read_col("GalaxyIndex", gal_idx); read_col("CentralGalaxyIndex", central_gal_idx); read_col("SimulationHaloIndex", sim_halo_idx);
            read_col("SAGEHaloIndex", sage_halo_idx); read_col("Type", type);
            read_col("mergeType", merge_type); read_col("mergeIntoID", merge_into_id); read_col("mergeIntoSnapNum", merge_into_snap);
            read_col("dT", dt);
            read_col("Len", len);

            for(size_t i=0; i<n_gals_snap; ++i) {
                sage::galaxy g;
                g.snapshot = snap;
                g.type = type[i];
                g.galaxy_idx = gal_idx[i];
                g.central_galaxy_idx = central_gal_idx[i];
                g.sage_halo_idx = sage_halo_idx[i];
                g.sage_tree_idx = tree_indices[i];
                g.simulation_halo_idx = sim_halo_idx[i];

                g.pos[0] = posx[i]; g.pos[1] = posy[i]; g.pos[2] = posz[i];
                g.vel[0] = velx[i]; g.vel[1] = vely[i]; g.vel[2] = velz[i];
                g.spin[0] = spinx[i]; g.spin[1] = spiny[i]; g.spin[2] = spinz[i];

                g.mvir = mvir[i]; g.central_mvir = central_mvir[i];
                g.rvir = rvir[i]; g.vvir = vvir[i]; g.vmax = vmax[i]; g.vel_disp = vel_disp[i];

                g.cold_gas = cold_gas[i]; g.stellar_mass = stellar_mass[i]; g.bulge_mass = bulge_mass[i];
                g.hot_gas = hot_gas[i]; g.ejected_mass = ejected_mass[i]; g.blackhole_mass = blackhole_mass[i]; g.ics = ics[i];

                g.metals_cold_gas = metals_cold_gas[i]; g.metals_stellar_mass = metals_stellar_mass[i]; g.metals_bulge_mass = metals_bulge_mass[i];
                g.metals_hot_gas = metals_hot_gas[i]; g.metals_ejected_mass = metals_ejected_mass[i]; g.metals_ics = metals_ics[i];

                g.sfr_disk = sfr_disk[i]; g.sfr_bulge = sfr_bulge[i]; g.sfr_disk_z = sfr_disk_z[i]; g.sfr_bulge_z = sfr_bulge_z[i];

                g.disk_scale_radius = disk_radius[i]; g.cooling = cooling[i]; g.heating = heating[i];
                g.quasar_mode_bh_accretion_mass = quasar_acc[i];
                g.time_of_last_major_merger = last_major[i]; g.time_of_last_minor_merger = last_minor[i];
                g.outflow_rate = outflow[i];

                g.infall_mvir = infall_mvir[i]; g.infall_vmax = infall_vmax[i]; g.infall_vvir = infall_vvir[i];

                g.merge_type = merge_type[i];
                g.merge_into_id = merge_into_id[i];
                g.merge_into_snapshot = merge_into_snap[i];
                g.dt = dt[i] / 1000.0; // Convert to Gyrs

                g.num_particles = len[i];

                // Initialize others
                g.descendant = -1;
                g.global_descendant = -1;
                g.local_index = i; // Store snapshot index for merging logic
                g.global_index = -1; // Will be set later

                trees[tree_indices[i]].push_back(g);
            }
        }

        // Process trees
        for(size_t t=0; t<n_trees; ++t) {
            auto& tree_gals = trees[t];
            if(tree_gals.empty()) {
                tree_displs_buf.write(current_tree_displ);
                tree_cnts_buf.write(0);
                continue;
            }

            // Link descendants
            std::unordered_map<long long, int> desc_map;
            std::unordered_multimap<hpc::varray<unsigned, 2>, unsigned> merge_map;

            for(size_t i=0; i<tree_gals.size(); ++i) {
                auto& g = tree_gals[i];
                unsigned file_gal_idx = g.local_index;

                // Update global index for output
                g.global_index = current_gal_global_idx + i;

                // 1. If I merge into someone, register myself in merge_map
                if(g.merge_into_id != -1) {
                    hpc::varray<unsigned, 2> id{(unsigned)g.merge_into_id, (unsigned)g.merge_into_snapshot};
                    merge_map.emplace(id, i);
                }

                // 2. If I am a descendant of someone (via GalaxyIndex), link them
                if(desc_map.count(g.galaxy_idx)) {
                    int par = desc_map[g.galaxy_idx];
                    tree_gals[par].descendant = i;
                    tree_gals[par].global_descendant = g.global_index;
                }

                // 3. If someone merged into me (via merge_map), link them
                hpc::varray<unsigned, 2> my_id{file_gal_idx, (unsigned)g.snapshot};
                auto range = merge_map.equal_range(my_id);
                for(auto it = range.first; it != range.second; ++it) {
                    int par = it->second;
                    tree_gals[par].descendant = i;
                    tree_gals[par].global_descendant = g.global_index;
                }
                merge_map.erase(my_id);

                // 4. Register myself for future descendants
                if(g.merge_into_id != -1) {
                    desc_map[g.galaxy_idx] = -1; // Merged, so I'm not a main progenitor
                } else {
                    desc_map[g.galaxy_idx] = i;
                }
            }

            // Reset local_index for reordering
            for(size_t i=0; i<tree_gals.size(); ++i) tree_gals[i].local_index = i;

            // Reorder
            _depthfirst_ordering(hpc::view<std::vector<sage::galaxy>>(tree_gals));

            // Write
            gals_buf.write(tree_gals);
            tree_displs_buf.write(current_tree_displ);
            tree_cnts_buf.write((unsigned)tree_gals.size());

            current_gal_global_idx += tree_gals.size();
            current_tree_displ += tree_gals.size();
        }
    }

    // Close buffers
    gals_buf.close();
    tree_displs_buf.close();
    tree_cnts_buf.close();

    // Write final displacement
    if(_comm->rank() == _comm->size() - 1) {
        hpc::h5::dataspace mem_space(1);
        hpc::h5::dataspace file_space = tree_displs_dset.dataspace();
        hsize_t start = total_trees;
        hsize_t count = 1;
        file_space.select_hyperslab(H5S_SELECT_SET, count, start);
        tree_displs_dset.write(&total_gals, hpc::h5::datatype::native_ullong, mem_space, file_space);
    }

    // Write metadata
    if(_comm->rank() == 0) {
        std::vector<double> zs(_redshifts.size());
        std::copy(_redshifts.begin(), _redshifts.end(), zs.begin());
        out_file.write_serial("snapshot_redshifts", zs);
        out_file.write<double>("cosmology/box_size", _box_size);
        out_file.write<double>("cosmology/hubble", _hubble);
        out_file.write<double>("cosmology/omega_l", _omega_l);
        out_file.write<double>("cosmology/omega_m", _omega_m);
    }
}

void sage2kdtree_application::_load_param(hpc::fs::path const &fn) {
    if (_verb >= 2 && _comm->rank() == 0) {
        std::cout << "  → Loading parameters from " << fn.filename() << std::endl;
    }

    std::ifstream file(fn.native());
    if (!file.good()) {
        throw PipelineException(1, "Failed to open parameter file: " + fn.string() +
                                "\n      Check that the file exists and is readable");
    }

    std::string line;
    int line_num = 0;
    while(std::getline(file, line)) {
        line_num++;
        if(line.empty() || line[0] == '%')
            continue;

        std::stringstream ss(line);
        std::string       key;
        ss >> key;

        if(key == "BoxSize")
            ss >> _box_size;
        else if(key == "Hubble_h") {
            ss >> _hubble;
            _hubble *= 100.0;
        } else if(key == "Omega_Lambda" || key == "OmegaLambda")
            ss >> _omega_l;
        else if(key == "Omega_m" || key == "Omega")
            ss >> _omega_m;
    }

    // Validate that we got the critical parameters
    if (_box_size <= 0.0 || _hubble <= 0.0) {
        throw PipelineException(1, "Invalid or missing cosmology parameters in " + fn.filename().string() +
                                "\n      Required: BoxSize > 0, Hubble_h > 0" +
                                "\n      Got: BoxSize=" + std::to_string(_box_size) +
                                ", Hubble_h=" + std::to_string(_hubble/100.0));
    }

    if (_verb >= 2 && _comm->rank() == 0) {
        std::cout << "    Loaded: BoxSize=" << _box_size
                  << ", Hubble_h=" << _hubble/100.0
                  << ", Omega_m=" << _omega_m
                  << ", Omega_Lambda=" << _omega_l << std::endl;
    }
}

void sage2kdtree_application::_load_redshifts(hpc::fs::path const &fn) {
    if (_verb >= 2 && _comm->rank() == 0) {
        std::cout << "  → Loading redshifts from " << fn.filename() << std::endl;
    }

    std::ifstream file(fn.native());
    if (!file.good()) {
        throw PipelineException(1, "Failed to open expansion factor list: " + fn.string() +
                                "\n      Check that the file exists and is readable");
    }

    std::string line;
    int line_num = 0;
    while(std::getline(file, line)) {
        line_num++;
        if(line.empty() || line[0] == '%')
            continue;

        std::stringstream ss(line);
        double            a;
        ss >> a;

        if (a <= 0.0 || a > 1.0) {
            throw PipelineException(1, "Invalid scale factor at line " + std::to_string(line_num) +
                                    " in " + fn.filename().string() +
                                    ": a=" + std::to_string(a) +
                                    "\n      Scale factors must be in range (0, 1]");
        }

        _redshifts.push_back(1.0 / a - 1.0);
    }

    if (_redshifts.empty()) {
        throw PipelineException(1, "No valid scale factors found in " + fn.string() +
                                "\n      File should contain one scale factor per line");
    }

    if (_verb >= 2 && _comm->rank() == 0) {
        std::cout << "    Loaded " << _redshifts.size() << " redshifts "
                  << "(z_min=" << *std::min_element(_redshifts.begin(), _redshifts.end())
                  << ", z_max=" << *std::max_element(_redshifts.begin(), _redshifts.end())
                  << ")" << std::endl;
    }
}

void sage2kdtree_application::_depthfirst_ordering(hpc::view<std::vector<sage::galaxy>> gals) {
    // Cache the smallest and largest global index
    std::array<unsigned long long, 2> gidxs{(unsigned long long)gals.front().global_index,
                                            (unsigned long long)gals.back().global_index};

    // Begin by finding the roots and parents
    std::list<unsigned>               roots;
    std::multimap<unsigned, unsigned> parents;
    for(unsigned ii = 0; ii < gals.size(); ++ii) {
        if(gals[ii].descendant == -1)
            roots.push_back(ii);
        else
            parents.emplace(gals[ii].descendant, ii);
    }

    // Process roots one at a time to create a new ordering
    unsigned              dfi = 0;
    std::vector<unsigned> map(gals.size());
    for(auto const &root : roots)
        _depthfirst_recurse(root, gals, map, dfi, parents);
    ASSERT(dfi == gals.size());

    // Permute the galaxies
    hpc::permute(gals.begin(), gals.end(), map.begin());

    // Reset global indices
    for(unsigned ii = 0; ii < gals.size(); ++ii) {
        gals[ii].global_index = gidxs[0] + ii;
    }

    // Remap local indices and descendants
    for(unsigned ii = 0; ii < gals.size(); ++ii) {
        gals[ii].local_index = map[gals[ii].local_index];
        if(gals[ii].descendant != -1) {
            gals[ii].descendant        = map[gals[ii].descendant];
            gals[ii].global_descendant = gals[gals[ii].descendant].global_index;
        }
    }
}

unsigned sage2kdtree_application::_depthfirst_recurse(
    unsigned idx,
    hpc::view<std::vector<sage::galaxy>> &gals,
    std::vector<unsigned> &map,
    unsigned &dfi,
    std::multimap<unsigned, unsigned> const &parents) {
    map[idx]          = dfi++;
    gals[idx].subsize = 1;
    auto rng          = parents.equal_range(idx);
    while(rng.first != rng.second) {
        gals[idx].subsize += _depthfirst_recurse(rng.first->second, gals, map, dfi, parents);
        ++rng.first;
    }
    return gals[idx].subsize;
}

//==============================================================================
// Phase 2: Add Traversal Metadata
// Based on sageimport.cc lines 106-439
//==============================================================================

void sage2kdtree_application::phase2_add_traversal_metadata() {
    int my_rank = _comm->rank();
    int n_ranks = _comm->size();

    // Rank 0 prepares the output file and copies metadata
    if (my_rank == 0) {
        hpc::h5::file in_file(_depthfirst_fn.native(), H5F_ACC_RDONLY);
        hpc::h5::file out_file(_enhanced_fn.native(), H5F_ACC_TRUNC);

        LOGILN("Copying metadata groups...");
        try { hpc::h5::copy(in_file, "cosmology", out_file); } catch(...) {}
        try { hpc::h5::copy(in_file, "snapshot_redshifts", out_file); } catch(...) {}
        try { hpc::h5::copy(in_file, "tree_counts", out_file); } catch(...) {}
        try { hpc::h5::copy(in_file, "tree_displs", out_file); } catch(...) {}

        _create_import_settings_sidecar();
        in_file.close();
        out_file.close();
    }

    // Wait for rank 0
    _comm->barrier();

    // All ranks process trees in parallel
    hpc::h5::file in_file(_depthfirst_fn.native(), H5F_ACC_RDONLY, *_comm);
    hpc::h5::file out_file(_enhanced_fn.native(), H5F_ACC_RDWR, *_comm);

    hpc::h5::dataset in_ds(in_file, "galaxies");
    hpc::h5::dataset out_ds;

    // Read tree counts
    hpc::h5::dataset counts_ds(in_file, "tree_counts");
    hsize_t n_trees = counts_ds.dataspace().size();

    std::vector<int> tree_counts(n_trees);
    counts_ds.read(tree_counts.data(), hpc::h5::datatype::native_int);

    // Calculate tree offsets
    std::vector<long long> tree_offsets(n_trees + 1, 0);
    for(size_t i=0; i<n_trees; ++i) {
        tree_offsets[i+1] = tree_offsets[i] + tree_counts[i];
    }

    // Determine which trees this rank processes (round-robin)
    size_t trees_per_rank = n_trees / n_ranks;
    size_t start_tree = my_rank * trees_per_rank;
    size_t end_tree = (my_rank == n_ranks - 1) ? n_trees : (my_rank + 1) * trees_per_rank;

    // Prepare memory types
    hpc::h5::datatype in_file_type = in_ds.datatype();
    hid_t in_mem_type_id = H5Tget_native_type(in_file_type.id(), H5T_DIR_DEFAULT);
    hpc::h5::datatype in_mem_type(in_mem_type_id);
    size_t in_struct_size = in_mem_type.size();

    // Find descendant field offset
    size_t descendant_offset = 0;
    bool found_descendant = false;
    for(unsigned i=0; i<in_mem_type.n_members(); ++i) {
        if(in_mem_type.member_name(i) == "descendant") {
            descendant_offset = in_mem_type.member_offset(i);
            found_descendant = true;
            break;
        }
    }
    if(!found_descendant) {
        throw PipelineException(2, "Descendant field not found in Phase 1 output");
    }

    // Create output memory type with 5 additional fields
    size_t out_struct_size = in_struct_size + 4 * sizeof(long long) + sizeof(int);
    hid_t out_mem_type_id = H5Tcreate(H5T_COMPOUND, out_struct_size);

    // Copy members from input type
    for(unsigned i=0; i<in_mem_type.n_members(); ++i) {
        std::string name = in_mem_type.member_name(i);
        size_t offset = in_mem_type.member_offset(i);
        hid_t member_type = H5Tget_member_type(in_mem_type.id(), i);
        H5Tinsert(out_mem_type_id, name.c_str(), offset, member_type);
        H5Tclose(member_type);
    }

    // Add new fields
    H5Tinsert(out_mem_type_id, "globaltreeid", in_struct_size, H5T_NATIVE_LLONG);
    H5Tinsert(out_mem_type_id, "breadthfirst_traversalorder", in_struct_size + sizeof(long long), H5T_NATIVE_LLONG);
    H5Tinsert(out_mem_type_id, "depthfirst_traversalorder", in_struct_size + 2 * sizeof(long long), H5T_NATIVE_LLONG);
    H5Tinsert(out_mem_type_id, "subtree_count", in_struct_size + 3 * sizeof(long long), H5T_NATIVE_LLONG);
    H5Tinsert(out_mem_type_id, "localgalaxyid", in_struct_size + 4 * sizeof(long long), H5T_NATIVE_INT);

    hpc::h5::datatype out_mem_type(out_mem_type_id);

    // Create output dataset
    out_ds.create(out_file, "galaxies", out_mem_type, (hsize_t)tree_offsets[n_trees]);

    // Process trees
    std::vector<char> in_buffer;
    std::vector<char> out_buffer;

    for(size_t t = start_tree; t < end_tree; ++t) {
        long long count = tree_counts[t];
        if (count == 0) continue;

        long long start = tree_offsets[t];

        in_buffer.resize(count * in_struct_size);
        out_buffer.resize(count * out_struct_size);

        // Read galaxies for this tree
        hpc::h5::dataspace mem_space(count);
        hpc::h5::dataspace file_space = in_ds.dataspace();
        file_space.select_hyperslab(H5S_SELECT_SET, (hsize_t)count, (hsize_t)start);
        in_ds.read(in_buffer.data(), in_mem_type, mem_space, file_space);

        // Build tree structure
        std::vector<Node> nodes(count);

        for(int i=0; i<count; ++i) {
            nodes[i].id = i;
            nodes[i].original_index = start + i;

            // Read Descendant (local index relative to tree start)
            int descendant_local = *reinterpret_cast<int*>(in_buffer.data() + i * in_struct_size + descendant_offset);

            if (descendant_local != -1) {
                if (descendant_local >= 0 && descendant_local < count) {
                    nodes[i].descendant = &nodes[descendant_local];
                }
            }
        }

        // Link progenitors
        for(int i=0; i<count; ++i) {
            if (nodes[i].descendant) {
                nodes[i].descendant->progenitors.push_back(&nodes[i]);
            }
        }

        // Find roots (nodes with no descendant)
        std::vector<Node*> roots;
        for(int i=0; i<count; ++i) {
            if (!nodes[i].descendant) {
                roots.push_back(&nodes[i]);
            }
        }

        // BFS traversal
        long long bfs_counter = 0;
        std::deque<Node*> queue;
        for(auto* root : roots) queue.push_back(root);

        while(!queue.empty()) {
            Node* u = queue.front();
            queue.pop_front();
            u->bfs_idx = bfs_counter++;

            // Sort progenitors by original index for determinism
            std::sort(u->progenitors.begin(), u->progenitors.end(), [](Node* a, Node* b){
                return a->original_index < b->original_index;
            });

            for(auto* v : u->progenitors) {
                queue.push_back(v);
            }
        }

        // DFS traversal (reverse root order to match Python stack-based DFS)
        long long dfs_counter = 0;
        for(auto it = roots.rbegin(); it != roots.rend(); ++it) {
            _dfs_visit(*it, dfs_counter);
        }

        // Pack data into output buffer
        for(int i=0; i<count; ++i) {
            char* in_ptr = in_buffer.data() + i * in_struct_size;
            char* out_ptr = out_buffer.data() + i * out_struct_size;

            // Copy original data
            std::memcpy(out_ptr, in_ptr, in_struct_size);

            // Append new fields
            char* extra_ptr = out_ptr + in_struct_size;
            long long* extra_ll = reinterpret_cast<long long*>(extra_ptr);
            extra_ll[0] = start + i; // globaltreeid (global galaxy index)
            extra_ll[1] = nodes[i].bfs_idx;
            extra_ll[2] = nodes[i].dfs_idx;
            extra_ll[3] = nodes[i].subtree_count;

            int* extra_int = reinterpret_cast<int*>(extra_ptr + 4 * sizeof(long long));
            extra_int[0] = i; // localgalaxyid
        }

        // Write back
        out_ds.write(out_buffer.data(), out_mem_type, mem_space, file_space);
    }

    // Create XML sidecar for Phase 2 output (needed by Phase 3)
    if (_comm->rank() == 0) {
        // Read original fields from Phase 1 XML
        hpc::fs::path phase1_xml = _depthfirst_fn;
        phase1_xml.replace_extension(".xml");
        tao::xml_dict phase1_sidecar = tao::data_dict::getSidecar(phase1_xml.native(), "/sageinput");
        std::vector<tao::data_dict_field> expanded_fields = tao::data_dict::getFieldsANY(phase1_sidecar, "/sageinput/Field");

        // Add the 5 new fields that Phase 2 added
        int next_order = expanded_fields.size() + 1;

        tao::data_dict_field globaltreeid_field;
        globaltreeid_field._name = "globaltreeid";
        globaltreeid_field._label = "globaltreeid";
        globaltreeid_field._description = "globaltreeid";
        globaltreeid_field._order = std::to_string(next_order++);
        globaltreeid_field._units = "";
        globaltreeid_field._group = "general";
        globaltreeid_field._type = tao::batch<real_type>::LONG_LONG;
        expanded_fields.push_back(globaltreeid_field);

        tao::data_dict_field bfs_field;
        bfs_field._name = "breadthfirst_traversalorder";
        bfs_field._label = "breadthfirst_traversalorder";
        bfs_field._description = "breadthfirst_traversalorder";
        bfs_field._order = std::to_string(next_order++);
        bfs_field._units = "";
        bfs_field._group = "general";
        bfs_field._type = tao::batch<real_type>::LONG_LONG;
        expanded_fields.push_back(bfs_field);

        tao::data_dict_field dfs_field;
        dfs_field._name = "depthfirst_traversalorder";
        dfs_field._label = "depthfirst_traversalorder";
        dfs_field._description = "depthfirst_traversalorder";
        dfs_field._order = std::to_string(next_order++);
        dfs_field._units = "";
        dfs_field._group = "general";
        dfs_field._type = tao::batch<real_type>::LONG_LONG;
        expanded_fields.push_back(dfs_field);

        tao::data_dict_field subtree_field;
        subtree_field._name = "subtree_count";
        subtree_field._label = "subtree_count";
        subtree_field._description = "subtree_count";
        subtree_field._order = std::to_string(next_order++);
        subtree_field._units = "";
        subtree_field._group = "general";
        subtree_field._type = tao::batch<real_type>::LONG_LONG;
        expanded_fields.push_back(subtree_field);

        tao::data_dict_field localgalid_field;
        localgalid_field._name = "localgalaxyid";
        localgalid_field._label = "localgalaxyid";
        localgalid_field._description = "localgalaxyid";
        localgalid_field._order = std::to_string(next_order++);
        localgalid_field._units = "";
        localgalid_field._group = "general";
        localgalid_field._type = tao::batch<real_type>::INTEGER;
        expanded_fields.push_back(localgalid_field);

        // Save XML sidecar with expanded field list
        tao::data_dict::saveSidecar(_enhanced_fn.native(), expanded_fields, false);
    }

    _comm->barrier();
}

void sage2kdtree_application::_dfs_visit(Node* u, long long& counter) {
    u->dfs_idx = counter++;
    u->subtree_count = 1;

    // Sort progenitors descending to match Python's stack-based DFS (LIFO)
    std::sort(u->progenitors.begin(), u->progenitors.end(), [](Node* a, Node* b){
        return a->original_index > b->original_index;
    });

    for(auto* v : u->progenitors) {
        _dfs_visit(v, counter);
        u->subtree_count += v->subtree_count;
    }
}

void sage2kdtree_application::_create_import_settings_sidecar() {
    std::string xml_fn = _depthfirst_fn.string();
    xml_fn = xml_fn.substr(0, xml_fn.length() - 3) + "_import_settings.xml";

    std::ofstream xml(xml_fn);
    xml << "<settings>\n";
    xml << "  <sageinput>\n";

    // Write field definitions (would need to read from Phase 1 sidecar, but simplified here)
    xml << "  </sageinput>\n";
    xml << "  <RunningSettings>\n";
    xml << "    <InputFile>" << _depthfirst_fn.string() << "</InputFile>\n";
    xml << "    <OutputFile>" << _enhanced_fn.string() << "</OutputFile>\n";
    xml << "    <BoxSize>" << _box_size << "</BoxSize>\n";
    xml << "    <Hubble_h>" << (_hubble / 100.0) << "</Hubble_h>\n";
    xml << "  </RunningSettings>\n";
    xml << "  <TreeMappings>\n";
    xml << "    <GlobalIndex>global_index</GlobalIndex>\n";
    xml << "    <Descendant>descendant</Descendant>\n";
    xml << "    <SnapNum>snapshot</SnapNum>\n";
    xml << "  </TreeMappings>\n";
    xml << "</settings>\n";
}

//==============================================================================
// Phase 3: Tree → Snapshot
// Adapted from dstreeinit.cc tree2sageANY (lines 509-591)
//==============================================================================

void sage2kdtree_application::phase3_tree_to_snapshot() {
    if (_verb >= 1 && _comm->rank() == 0) {
        std::cout << "Phase 3: Converting tree order to snapshot order..." << std::endl;
    }

    // Open enhanced tree file and read metadata
    hpc::h5::file tree_file(_enhanced_fn.string(), H5F_ACC_RDONLY, *_comm);
    hpc::h5::dataset galaxies_ds(tree_file, "galaxies");

    // Get datatype info using tao::data_dict
    hpc::fs::path xml_fn = _enhanced_fn;
    xml_fn.replace_extension(".xml");
    tao::xml_dict sidecar_xml = tao::data_dict::getSidecar(xml_fn.native(), "/sageinput");
    std::vector<tao::data_dict_field> fields = tao::data_dict::getFieldsANY(sidecar_xml, "/sageinput/Field");

    int nfields = fields.size();
    int snapnum_idx = tao::data_dict::findField(fields, "snapnum");
    if (snapnum_idx < 0) {
        snapnum_idx = tao::data_dict::findField(fields, "snapshot");
    }
    if (snapnum_idx < 0) {
        snapnum_idx = tao::data_dict::findField(fields, "SnapNum");
    }
    if (snapnum_idx < 0) {
        throw PipelineException(3, "Could not find snapnum field in tree file");
    }

    // Count galaxies per snapshot
    std::vector<unsigned long long> snap_counts = _count_galaxies_by_snapshot(tree_file, nfields, snapnum_idx, fields);

    unsigned long long total_gals = 0;
    for (auto cnt : snap_counts) total_gals += cnt;

    if (_verb >= 2 && _comm->rank() == 0) {
        std::cout << "  Total galaxies: " << total_gals << std::endl;
        std::cout << "  Redistributing into 64 snapshots..." << std::endl;
    }

    // Create output file
    hpc::h5::file out_file(_bysnap_fn.string(), H5F_ACC_TRUNC, *_comm);

    // Get total galaxy count
    hpc::h5::dataspace file_space = galaxies_ds.dataspace();
    std::vector<hsize_t> dims_vec(1);
    file_space.simple_extent_dims<std::vector<hsize_t>>(dims_vec);
    unsigned long long n_gals = dims_vec[0];

    // Distribute galaxies across ranks
    unsigned long long n_loc = n_gals / _comm->size();
    unsigned long long offs = _comm->rank() * n_loc;
    if (_comm->rank() == _comm->size() - 1) {
        n_loc = n_gals - offs;
    }

    // Read chunk of galaxies
    const unsigned long long chunk_size = 10000;
    hpc::h5::datatype mem_type = galaxies_ds.datatype();
    size_t type_size = mem_type.size();

    // Prepare snapshot datasets (create all 64 snapshots)
    std::vector<hpc::h5::dataset> snap_datasets(64);
    for (int snap = 0; snap < 64; ++snap) {
        std::ostringstream ss;
        ss << "snapshot" << std::setfill('0') << std::setw(3) << snap;

        // Create dataset with proper dimensions
        // hsize_t snap_dims[1] = {snap_counts[snap]};
        hpc::h5::dataspace snap_space(snap_counts[snap]);
        snap_datasets[snap] = hpc::h5::dataset(out_file, ss.str(), mem_type, snap_space, hpc::h5::property_list());
    }

    // Buffers for each snapshot
    std::vector<std::vector<char>> snap_buffers(64);
    std::vector<unsigned long long> snap_offsets(64, 0);

    // Process in chunks
    for (unsigned long long chunk_offs = offs; chunk_offs < offs + n_loc; chunk_offs += chunk_size) {
        unsigned long long chunk_n = std::min(chunk_size, offs + n_loc - chunk_offs);

        // Read chunk
        std::vector<char> chunk_data(chunk_n * type_size);
        hpc::h5::dataspace mem_space(chunk_n);

        file_space.select_hyperslab(H5S_SELECT_SET, (hsize_t)chunk_n, (hsize_t)chunk_offs);

        galaxies_ds.read(chunk_data.data(), mem_type, mem_space, file_space);

        // Distribute to snapshot buffers based on snapnum
        for (unsigned long long ii = 0; ii < chunk_n; ++ii) {
            char *gal_ptr = chunk_data.data() + ii * type_size;

            // Extract snapnum field
            int snapnum;
            memcpy(&snapnum, gal_ptr + fields[snapnum_idx]._offset, sizeof(int));

            if (snapnum >= 0 && snapnum < 64) {
                snap_buffers[snapnum].insert(
                    snap_buffers[snapnum].end(),
                    gal_ptr,
                    gal_ptr + type_size
                );
            }
        }

        // Flush buffers when they get large
        for (int snap = 0; snap < 64; ++snap) {
            if (snap_buffers[snap].size() >= chunk_size * type_size) {
                unsigned long long n_write = snap_buffers[snap].size() / type_size;

                hpc::h5::dataspace write_mem_space(n_write);

                hpc::h5::dataspace snap_file_space = snap_datasets[snap].dataspace();
                snap_file_space.select_hyperslab(H5S_SELECT_SET, (hsize_t)n_write, (hsize_t)snap_offsets[snap]);

                snap_datasets[snap].write(snap_buffers[snap].data(), mem_type, write_mem_space, snap_file_space);

                snap_offsets[snap] += n_write;
                snap_buffers[snap].clear();
            }
        }
    }

    // Flush remaining buffers
    for (int snap = 0; snap < 64; ++snap) {
        if (!snap_buffers[snap].empty()) {
            unsigned long long n_write = snap_buffers[snap].size() / type_size;

            hpc::h5::dataspace write_mem_space(n_write);

            hpc::h5::dataspace snap_file_space = snap_datasets[snap].dataspace();
            snap_file_space.select_hyperslab(H5S_SELECT_SET, (hsize_t)n_write, (hsize_t)snap_offsets[snap]);

            snap_datasets[snap].write(snap_buffers[snap].data(), mem_type, write_mem_space, snap_file_space);
        }
    }

    _comm->barrier();

    // Copy cosmology and snapshot_redshifts metadata
    if (_comm->rank() == 0) {
        try { hpc::h5::copy(tree_file, "cosmology", out_file); } catch(...) {}
        try { hpc::h5::copy(tree_file, "snapshot_redshifts", out_file); } catch(...) {}
    }

    _comm->barrier();

    // Write XML sidecar for bysnap file
    if (_comm->rank() == 0) {
        hpc::fs::path xml_out = _bysnap_fn;
        xml_out.replace_extension(".xml");
        std::ofstream xml(xml_out.string());
        xml << "<?xml version=\"1.0\"?>\n";
        xml << "<settings>\n";
        xml << "  <snapshots>\n";
        for (int snap = 0; snap < 64; ++snap) {
            xml << "    <snapshot" << std::setfill('0') << std::setw(3) << snap << ">" << snap_counts[snap] << "</snapshot" << std::setfill('0') << std::setw(3) << snap << ">\n";
        }
        xml << "  </snapshots>\n";
        xml << "</settings>\n";
    }

    if (_verb >= 1 && _comm->rank() == 0) {
        std::cout << "  Phase 3 complete: " << _bysnap_fn << std::endl;
    }
}

std::vector<unsigned long long> sage2kdtree_application::_count_galaxies_by_snapshot(
    hpc::h5::file &file,
    int nfields,
    int snapnum_index,
    std::vector<tao::data_dict_field> const &fields) {

    std::vector<unsigned long long> counts(64, 0);

    hpc::h5::dataset galaxies_ds(file, "galaxies");
    hpc::h5::dataspace file_space = galaxies_ds.dataspace();
    std::vector<hsize_t> dims_vec(1);
    file_space.simple_extent_dims<std::vector<hsize_t>>(dims_vec);
    unsigned long long n_gals = dims_vec[0];

    // Distribute galaxies across ranks
    unsigned long long n_loc = n_gals / _comm->size();
    unsigned long long offs = _comm->rank() * n_loc;
    if (_comm->rank() == _comm->size() - 1) {
        n_loc = n_gals - offs;
    }

    const unsigned long long chunk_size = 100000;
    std::vector<unsigned long long> local_counts(64, 0);

    // Create a compound type with just the snapnum field
    std::string snapnum_field_name = fields[snapnum_index]._name;
    hpc::h5::datatype snapnum_type = hpc::h5::datatype::native_int;
    hpc::h5::datatype read_type(H5Tcreate(H5T_COMPOUND, sizeof(int)));
    read_type.insert(snapnum_type, snapnum_field_name, 0);

    for (unsigned long long chunk_offs = offs; chunk_offs < offs + n_loc; chunk_offs += chunk_size) {
        unsigned long long chunk_n = std::min(chunk_size, offs + n_loc - chunk_offs);

        std::vector<int> snapnums(chunk_n);
        hpc::h5::dataspace mem_space(chunk_n);

        file_space.select_hyperslab(H5S_SELECT_SET, (hsize_t)chunk_n, (hsize_t)chunk_offs);

        galaxies_ds.read(snapnums.data(), read_type, mem_space, file_space);

        for (unsigned long long ii = 0; ii < chunk_n; ++ii) {
            if (snapnums[ii] >= 0 && snapnums[ii] < 64) {
                local_counts[snapnums[ii]]++;
            }
        }
    }

    // Reduce across all ranks
    counts = local_counts;
    // _comm->all_reduce(counts, MPI_SUM);
    MPI_Allreduce(MPI_IN_PLACE, counts.data(), counts.size(), MPI_UNSIGNED_LONG_LONG, MPI_SUM, _comm->mpi_comm());

    return counts;
}

//==============================================================================
// Phase 4: Build KD-Tree Index
// Adapted from dstreeinit.cc initANY (lines 127-208)
//==============================================================================

void sage2kdtree_application::phase4_build_kdtree_index() {
    if (_verb >= 1 && _comm->rank() == 0) {
        std::cout << "Phase 4: Building KD-tree spatial index..." << std::endl;
    }

    // Open snapshot-organized file from Phase 3
    hpc::h5::file snap_file(_bysnap_fn.string(), H5F_ACC_RDONLY, *_comm);

    // Create output KD-tree file
    hpc::h5::file out_file(_output_fn.string(), H5F_ACC_TRUNC, *_comm);

    // Note: Don't create groups explicitly - HDF5 will auto-create them
    // when we create datasets with paths like "lightcone/data" or "data/posx"

    // Get number of snapshots from snapshot_redshifts dataset
    int n_snapshots = 64;  // Default value
    if (snap_file.has_link("snapshot_redshifts")) {
        hpc::h5::dataset redshift_ds = snap_file.dataset("snapshot_redshifts");
        hpc::h5::dataspace redshift_sp = redshift_ds.dataspace();
        std::vector<hsize_t> redshift_dims(1);
        redshift_sp.simple_extent_dims<std::vector<hsize_t>>(redshift_dims);
        n_snapshots = redshift_dims[0];
    }

    // Get total galaxy count across all snapshots
    unsigned long long total_gals = 0;
    std::vector<unsigned long long> snap_counts(n_snapshots);
    std::vector<unsigned long long> snap_displs(n_snapshots + 1);  // n_snapshots + 1 final total

    for (int snap = 0; snap < n_snapshots; ++snap) {
        std::ostringstream ss;
        ss << "snapshot" << std::setfill('0') << std::setw(3) << snap;

        if (snap_file.has_link(ss.str())) {
            hpc::h5::dataset ds = snap_file.dataset(ss.str());
            hpc::h5::dataspace sp = ds.dataspace();
            std::vector<hsize_t> dims_vec(1);
            sp.simple_extent_dims<std::vector<hsize_t>>(dims_vec);
            snap_counts[snap] = dims_vec[0];
        } else {
            snap_counts[snap] = 0;
        }

        snap_displs[snap] = total_gals;
        total_gals += snap_counts[snap];
    }

    // Add final total as the last element
    snap_displs[n_snapshots] = total_gals;

    if (_verb >= 2 && _comm->rank() == 0) {
        std::cout << "  Total galaxies across all snapshots: " << total_gals << std::endl;
    }

    // Create lightcone/data compound dataset
    hpc::h5::datatype lc_type(H5Tcreate(H5T_COMPOUND, sizeof(double) * 3 + sizeof(unsigned long long) + sizeof(unsigned)));
    hpc::h5::datatype double_type = hpc::h5::datatype::ieee_f64le;
    hpc::h5::datatype ulonglong_type = hpc::h5::datatype::native_ullong;
    hpc::h5::datatype uint_type = hpc::h5::datatype::native_uint;

    lc_type.insert(double_type, "x", 0);
    lc_type.insert(double_type, "y", sizeof(double));
    lc_type.insert(double_type, "z", sizeof(double) * 2);
    lc_type.insert(ulonglong_type, "global_index", sizeof(double) * 3);
    lc_type.insert(uint_type, "subsize", sizeof(double) * 3 + sizeof(unsigned long long));

    // hsize_t lc_dims[1] = {total_gals};
    hpc::h5::dataspace lc_space(total_gals);
    hpc::h5::dataset lc_data(out_file, "lightcone/data", lc_type, lc_space, hpc::h5::property_list());

    // Create all data/* datasets upfront (before processing snapshots)
    // Get field info from the first non-empty snapshot
    if (_comm->rank() == 0 && total_gals > 0) {
        for (int snap = 0; snap < n_snapshots; ++snap) {
            std::ostringstream ss;
            ss << "snapshot" << std::setfill('0') << std::setw(3) << snap;

            if (snap_counts[snap] > 0 && snap_file.has_link(ss.str())) {
                hpc::h5::dataset snap_ds = snap_file.dataset(ss.str());
                hpc::h5::datatype compound_type = snap_ds.datatype();
                int nfields = compound_type.n_members();

                for (int ii = 0; ii < nfields; ++ii) {
                    std::string field_name = compound_type.member_name(ii);
                    std::string out_name = field_name;
                    std::transform(out_name.begin(), out_name.end(), out_name.begin(), ::tolower);

                    // Create data/* dataset for all fields (matching dstreeinit.cc)
                    hpc::h5::datatype field_type = compound_type.member_type(ii);
                    hpc::h5::dataspace field_space(total_gals);
                    hpc::h5::dataset field_ds(out_file, "data/" + out_name, field_type, field_space, hpc::h5::property_list());
                }
                break;  // Only need to process one snapshot to get field info
            }
        }
    }
    _comm->barrier();

    // Process each snapshot (including empty ones)
    unsigned long long displ = 0;
    for (int snap = 0; snap < n_snapshots; ++snap) {
        std::ostringstream ss;
        ss << "snapshot" << std::setfill('0') << std::setw(3) << snap;

        if (_verb >= 2 && _comm->rank() == 0) {
            std::cout << "  Processing " << ss.str() << " (" << snap_counts[snap] << " galaxies)..." << std::endl;
        }

        // Process all snapshots, even if empty
        if (snap_file.has_link(ss.str())) {
            _process_snapshot(snap_file, ss.str(), out_file, lc_data, displ);
        } else {
            // Write empty KD-tree structure for missing snapshot
            _write_empty_kdtree(out_file, snap);
        }
    }

    _comm->barrier();

    // Write SED data from tree file
    if (_comm->rank() == 0) {
        if (_verb >= 2) {
            std::cout << "  Writing SED data from tree file..." << std::endl;
        }
        _write_sed_data(out_file, _enhanced_fn);
    }

    _comm->barrier();

    // Write snapshot metadata
    if (_comm->rank() == 0) {
        hpc::h5::dataspace counts_space(n_snapshots);
        hpc::h5::dataset counts_ds(out_file, "snapshot_counts", ulonglong_type, counts_space, hpc::h5::property_list());
        counts_ds.write(snap_counts);

        hpc::h5::dataspace displs_space(n_snapshots + 1);  // n_snapshots + 1 final total
        hpc::h5::dataset displs_ds(out_file, "snapshot_displs", ulonglong_type, displs_space, hpc::h5::property_list());
        displs_ds.write(snap_displs);

        // Copy cosmology and snapshot_redshifts
        try { hpc::h5::copy(snap_file, "cosmology", out_file); } catch(...) {}
        try { hpc::h5::copy(snap_file, "snapshot_redshifts", out_file); } catch(...) {}

        // Create XML sidecar for final output file
        hpc::fs::path enhanced_xml = _enhanced_fn;
        enhanced_xml.replace_extension(".xml");
        tao::xml_dict sidecar_xml = tao::data_dict::getSidecar(enhanced_xml.native(), "/sageinput");
        std::vector<tao::data_dict_field> fields = tao::data_dict::getFieldsANY(sidecar_xml, "/sageinput/Field");
        tao::data_dict::saveSidecar(_output_fn.native(), fields, false);
    }

    _comm->barrier();

    if (_verb >= 1 && _comm->rank() == 0) {
        std::cout << "  Phase 4 complete: " << _output_fn << std::endl;
    }
}

void sage2kdtree_application::_process_snapshot(
    hpc::h5::file &file,
    std::string const &name,
    hpc::h5::file &out_file,
    hpc::h5::dataset &data,
    unsigned long long &displ) {

    // Read coordinates and metadata
    hpc::h5::dataset in_ds(file, name);
    std::array<std::vector<double>, 3> crds;
    std::vector<unsigned long long> tree_idxs;
    std::vector<unsigned> subsize;

    _read_coords(in_ds, crds, tree_idxs, subsize);

    unsigned long long n_gals = crds[0].size();

    // Extract snapshot number from name (e.g., "snapshot063" -> 63)
    int snap_num = std::stoi(name.substr(8));

    if (n_gals == 0) {
        // Write empty KD-tree structure for this snapshot
        _write_empty_kdtree(out_file, snap_num);
        return;
    }

    std::vector<unsigned long long> idxs(n_gals);
    std::iota(idxs.begin(), idxs.end(), 0);

    // Build KD-tree
    data_permuter dp(crds, tree_idxs, subsize, idxs);
    hpc::binary_partitioner<std::array<std::vector<double>, 3>::iterator, data_permuter> bp(crds.begin(), crds.end(), dp, _ppc);
    hpc::kdtree<> kdt;
    kdt.construct(bp);

    // Write lightcone/data fields separately (matching dstreeinit.cc lines 241-264)
    // Write x, y, z coordinates
    static const char *const coord_names[3] = {"x", "y", "z"};
    for (unsigned ii = 0; ii < 3; ++ii) {
        hpc::h5::datatype dt;
        dt.compound(sizeof(double));
        dt.insert(hpc::h5::datatype::native_double, coord_names[ii], 0);
        hpc::h5::property_list props(H5P_DATASET_XFER);
        props.set_preserve();
        data.write(crds[ii].data(), dt, n_gals, displ, hpc::mpi::comm::self, props);
    }

    // Write global_index
    {
        hpc::h5::datatype dt;
        dt.compound(sizeof(unsigned long long));
        dt.insert(hpc::h5::datatype::native_ullong, "global_index", 0);
        hpc::h5::property_list props(H5P_DATASET_XFER);
        props.set_preserve();
        data.write(tree_idxs.data(), dt, n_gals, displ, hpc::mpi::comm::self, props);
    }

    // Write subsize
    {
        hpc::h5::datatype dt;
        dt.compound(sizeof(unsigned));
        dt.insert(hpc::h5::datatype::native_uint, "subsize", 0);
        hpc::h5::property_list props(H5P_DATASET_XFER);
        props.set_preserve();
        data.write(subsize.data(), dt, n_gals, displ, hpc::mpi::comm::self, props);
    }

    // Write all other attributes to data/* datasets
    _write_attributes(in_ds, out_file, idxs, displ);

    // Write KD-tree structure
    _write_kdtree(out_file, name, kdt, bp, crds, displ);

    displ += n_gals;
}

void sage2kdtree_application::_read_coords(
    hpc::h5::dataset const &dset,
    std::array<std::vector<double>, 3> &crds,
    std::vector<unsigned long long> &tree_idxs,
    std::vector<unsigned> &subsize) {

    // Get galaxy count
    hpc::h5::dataspace sp = dset.dataspace();
    std::vector<hsize_t> dims_vec(1);
    sp.simple_extent_dims<std::vector<hsize_t>>(dims_vec);
    unsigned long long n_gals = dims_vec[0];

    if (n_gals == 0) {
        crds[0].clear();
        crds[1].clear();
        crds[2].clear();
        tree_idxs.clear();
        subsize.clear();
        return;
    }

    crds[0].resize(n_gals);
    crds[1].resize(n_gals);
    crds[2].resize(n_gals);
    tree_idxs.resize(n_gals);
    subsize.resize(n_gals);

    // Get compound type info and find field indices (matching dstreeinit.cc approach)
    hpc::h5::datatype compound_type = dset.datatype();
    int n_members = compound_type.n_members();

    // Find field indices by name
    int posx_index = -1, posy_index = -1, posz_index = -1;
    int globaltreeid_index = -1, subtree_count_index = -1;

    for (int ii = 0; ii < n_members; ++ii) {
        std::string name = compound_type.member_name(ii);
        if (name == "posx") posx_index = ii;
        else if (name == "posy") posy_index = ii;
        else if (name == "posz") posz_index = ii;
        else if (name == "globaltreeid") globaltreeid_index = ii;
        else if (name == "subtree_count") subtree_count_index = ii;
    }

    // Read all data as raw bytes
    size_t type_size = compound_type.size();
    std::vector<char> buffer(n_gals * type_size);
    hpc::h5::dataspace mem_space(n_gals);
    dset.read(buffer.data(), compound_type, mem_space, sp);

    // Extract fields by index and offset (coordinates are float in file, double in memory)
    for (unsigned long long ii = 0; ii < n_gals; ++ii) {
        char* record = buffer.data() + ii * type_size;

        if (posx_index >= 0) {
            float val;
            size_t offset = compound_type.member_offset(posx_index);
            memcpy(&val, record + offset, sizeof(float));
            crds[0][ii] = static_cast<double>(val);
        }
        if (posy_index >= 0) {
            float val;
            size_t offset = compound_type.member_offset(posy_index);
            memcpy(&val, record + offset, sizeof(float));
            crds[1][ii] = static_cast<double>(val);
        }
        if (posz_index >= 0) {
            float val;
            size_t offset = compound_type.member_offset(posz_index);
            memcpy(&val, record + offset, sizeof(float));
            crds[2][ii] = static_cast<double>(val);
        }
        if (globaltreeid_index >= 0) {
            size_t offset = compound_type.member_offset(globaltreeid_index);
            memcpy(&tree_idxs[ii], record + offset, sizeof(unsigned long long));
        }
        if (subtree_count_index >= 0) {
            size_t offset = compound_type.member_offset(subtree_count_index);
            memcpy(&subsize[ii], record + offset, sizeof(unsigned));
        }
    }
}

void sage2kdtree_application::_write_attributes(
    hpc::h5::dataset &in_ds,
    hpc::h5::file &out_file,
    std::vector<unsigned long long> const &idxs,
    unsigned long long displ) {

    // Matching dstreeinit.cc lines 283-315
    // Read entire snapshot using dynamic memory type (to include all fields), then write each field

    // Derive memory type from input dataset to ensure we get ALL fields (including Phase 2 additions)
    hpc::h5::datatype file_type = in_ds.datatype();
    hid_t mem_type_id = H5Tget_native_type(file_type.id(), H5T_DIR_DEFAULT);
    hpc::h5::datatype mem_type(mem_type_id);

    unsigned n_mems = mem_type.n_members();
    unsigned long long n_gals = in_ds.extent();

    if (n_gals == 0) return;

    // Read entire snapshot into memory using mem_type
    std::vector<char> buf(n_gals * mem_type.size());
    hpc::h5::dataspace mem_space(n_gals);
    hpc::h5::dataspace file_space = in_ds.dataspace();
    in_ds.read(buf.data(), mem_type, mem_space, file_space);

    // Permute to KD-tree spatial order (matching dstreeinit.cc line 297)
    std::vector<char> permuted_buf(n_gals * mem_type.size());
    for (unsigned long long jj = 0; jj < n_gals; ++jj) {
        unsigned long long src_idx = idxs[jj];
        memcpy(permuted_buf.data() + jj * mem_type.size(),
               buf.data() + src_idx * mem_type.size(),
               mem_type.size());
    }

    // Write each field (matching dstreeinit.cc lines 298-314)
    for (unsigned ii = 0; ii < n_mems; ++ii) {
        std::string name = mem_type.member_name(ii);
        std::string name_lowercase = name;
        std::transform(name_lowercase.begin(), name_lowercase.end(), name_lowercase.begin(), ::tolower);

        unsigned dt_offs = mem_type.member_offset(ii);
        hpc::h5::datatype dt = mem_type.member_type(ii);  // Native (little-endian) type

        // Extract field from permuted buffer
        std::vector<char> tmp(n_gals * dt.size());
        for (unsigned long long jj = 0; jj < n_gals; ++jj) {
            memcpy(tmp.data() + jj * dt.size(),
                   permuted_buf.data() + jj * mem_type.size() + dt_offs,
                   dt.size());
        }

        // Write to /data/{field_name_lowercase}
        hpc::h5::dataset out_ds(out_file, "data/" + name_lowercase);
        hpc::h5::dataspace out_mem_space(n_gals);
        hpc::h5::dataspace out_file_space(out_ds);
        out_file_space.select_range(displ, displ + n_gals);

        hpc::h5::property_list props(H5P_DATASET_XFER);
        props.set_preserve();
        out_ds.write(tmp.data(), dt, out_mem_space, out_file_space, hpc::mpi::comm::self, props);
    }
}

void sage2kdtree_application::_write_kdtree(
    hpc::h5::file &file,
    std::string const &snap_name,
    hpc::kdtree<> const &kdt,
    hpc::binary_partitioner<std::array<std::vector<double>, 3>::iterator, data_permuter> const &bp,
    std::array<std::vector<double>, 3> const &crds,
    unsigned long long displ) {

    // Build path prefix: "lightcone/snapshot###"
    std::string path = std::string("lightcone/") + snap_name;

    // Write bounds (matching dstreeinit.cc lines 327-334)
    {
        hpc::h5::derive der(2 * sizeof(double));
        der.add(hpc::h5::datatype::native_double, 0, hpc::h5::datatype::ieee_f64be, "minimum");
        der.add(hpc::h5::datatype::native_double, sizeof(double), hpc::h5::datatype::ieee_f64be, "maximum");
        hpc::h5::datatype mdt, fdt;
        der.commit(mdt, fdt);

        hpc::h5::dataset dset(file, path + "/bounds", fdt, 3);
        dset.write(kdt.bounds().data(), mdt, 3, 0);
    }

    // Write splits.
    {
        hpc::h5::derive der(sizeof(hpc::kdtree<>::split_type));
        der.add(hpc::h5::datatype::native_double,
                HOFFSET(hpc::kdtree<>::split_type, pos),
                hpc::h5::datatype::ieee_f64be,
                "position");
        der.add(hpc::h5::datatype::native_uint,
                HOFFSET(hpc::kdtree<>::split_type, dim),
                hpc::h5::datatype::std_i32be,
                "dimension");
        hpc::h5::datatype mdt, fdt;
        der.commit(mdt, fdt);

        hpc::h5::dataset dset(file, path + "/splits", fdt, kdt.splits().size());
        dset.write(kdt.splits().data(), mdt, kdt.splits().size(), 0);
    }

    // Write counts and offsets (matching dstreeinit.cc lines 356-362)
    // Use file.write() which writes with native byte order (little-endian)
    file.write(path + "/cell_counts", bp.counts());  // std::vector<unsigned> = 32-bit

    {
        std::vector<unsigned long long> offs(bp.offsets().size());
        for(unsigned ii = 0; ii < offs.size(); ++ii)
            offs[ii] = displ + bp.offsets()[ii];
        file.write(path + "/cell_offs", offs);  // std::vector<unsigned long long> = 64-bit
    }
}

void sage2kdtree_application::_write_empty_kdtree(hpc::h5::file &file, int snap_num) {
    // Write empty KD-tree structure for snapshots with no galaxies
    if (_comm->rank() != 0) return;  // Only rank 0 writes

    std::ostringstream ss;
    ss << "lightcone/snapshot" << std::setfill('0') << std::setw(3) << snap_num;
    std::string path = ss.str();

    hpc::h5::datatype ulonglong_type = hpc::h5::datatype::native_ullong;
    hpc::h5::datatype uint_type = hpc::h5::datatype::native_uint;
    hpc::h5::datatype double_type = hpc::h5::datatype::native_double;

    // Write empty bounds (all zeros for empty snapshots)
    {
        hpc::h5::derive der(2 * sizeof(double));
        der.add(hpc::h5::datatype::native_double, 0, hpc::h5::datatype::ieee_f64be, "minimum");
        der.add(hpc::h5::datatype::native_double, sizeof(double), hpc::h5::datatype::ieee_f64be, "maximum");

        hpc::h5::datatype mdt, fdt;
        der.commit(mdt, fdt);

        struct bound_type {
            double minimum, maximum;
        };
        bound_type bounds[3] = {
            {0.0, 0.0},
            {0.0, 0.0},
            {0.0, 0.0}
        };

        hpc::h5::dataset dset(file, path + "/bounds", fdt, 3);
        dset.write(bounds, mdt, 3, 0);
    }

    // Write empty splits (size 0)
    {
        hpc::h5::derive der(sizeof(double) + sizeof(unsigned));
        der.add(hpc::h5::datatype::native_double, 0, hpc::h5::datatype::ieee_f64be, "position");
        der.add(hpc::h5::datatype::native_uint, sizeof(double), hpc::h5::datatype::std_i32be, "dimension");

        hpc::h5::datatype mdt, fdt;
        der.commit(mdt, fdt);

        hpc::h5::dataset dset(file, path + "/splits", fdt, 0);
    }

    // Write cell_counts (1 root cell with 0 galaxies) - use file.write() for native byte order
    {
        std::vector<unsigned> counts(1, 0);  // 32-bit unsigned
        file.write(path + "/cell_counts", counts);
    }

    // Write cell_offs (1 offset = 0)
    {
        std::vector<unsigned long long> offs(1, 0);  // 64-bit unsigned
        file.write(path + "/cell_offs", offs);
    }
}

void sage2kdtree_application::_write_sed_data(hpc::h5::file &out_file, hpc::fs::path const &tree_fn) {
    // Open tree file to extract SED data
    hpc::h5::file tree_file(tree_fn.string(), H5F_ACC_RDONLY);

    if (!tree_file.has_link("galaxies")) {
        return;
    }

    hpc::h5::dataset tree_ds(tree_file, "galaxies");
    hpc::h5::dataspace tree_space = tree_ds.dataspace();
    std::vector<hsize_t> dims(tree_space.simple_extent_num_dims());
    tree_space.simple_extent_dims<std::vector<hsize_t>>(dims);
    unsigned long long n_gals = dims[0];

    // Create sed/data group
    hpc::h5::group sed_grp;
    sed_grp.create(out_file, "sed");

    // Define SED fields (from dstreeinit.cc sed_data_t)
    std::vector<std::string> sed_fields = {
        "DiskMass", "BulgeMass", "HotGas", "ColdGas", "EjectedMass",
        "BlackHoleMass", "MetalsDiskMass", "MetalsBulgeMass", "MetalsHotGas",
        "MetalsColdGas", "MetalsEjectedMass", "SfrDisk", "SfrBulge"
    };

    // Get compound type info
    hpc::h5::datatype compound_type = tree_ds.datatype();
    int nfields = compound_type.n_members();

    // Create SED compound type
    size_t sed_size = 0;
    std::vector<size_t> sed_offsets;
    std::vector<hpc::h5::datatype> sed_types;

    for (const auto &field_name : sed_fields) {
        for (int ii = 0; ii < nfields; ++ii) {
            if (compound_type.member_name(ii) == field_name) {
                hpc::h5::datatype field_type = compound_type.member_type(ii);
                sed_offsets.push_back(sed_size);
                sed_types.push_back(field_type);
                sed_size += field_type.size();
                break;
            }
        }
    }

    if (sed_offsets.empty()) {
        return;  // No SED fields found
    }

    // Create compound type for SED data
    hpc::h5::datatype sed_type;
    sed_type.compound(sed_size);
    
    for (size_t ii = 0; ii < sed_fields.size() && ii < sed_offsets.size(); ++ii) {
        sed_type.insert(sed_types[ii], sed_fields[ii], sed_offsets[ii]);
    }

    // Create dataset
    hpc::h5::dataspace sed_space(n_gals);
    hpc::h5::dataset sed_ds(sed_grp, "data", sed_type, sed_space);

    // Read and write SED data
    std::vector<char> sed_data(n_gals * sed_size);
    tree_ds.read((void *)sed_data.data(), sed_type, sed_space, tree_space);
    sed_ds.write((void const *)sed_data.data(), sed_type, sed_space, sed_space);
}

//==============================================================================
// Main Entry Point
//==============================================================================

#define HPC_APP_CLASS sage2kdtree_application
#include <libhpc/mpi/main.hh>
