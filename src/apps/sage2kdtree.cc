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
#include <memory>
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
#include <libtao/base/mandatory_fields.hh>
#include <chrono>
#include <sys/resource.h>
#include <iomanip>

#ifdef __APPLE__
#include <mach/mach.h>
#else
#include <unistd.h>
#include <ios>
#endif

using namespace ::tao;
using namespace ::hpc;
using namespace ::sage;

// Helper for memory reporting
long get_peak_rss_kb() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
#ifdef __APPLE__
    return usage.ru_maxrss / 1024; // macOS reports bytes
#else
    return usage.ru_maxrss;        // Linux reports KB
#endif
}

long get_current_rss_kb() {
#ifdef __APPLE__
    struct mach_task_basic_info info;
    mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
    if ( task_info( mach_task_self(), MACH_TASK_BASIC_INFO,
                    (task_info_t)&info, &infoCount ) != KERN_SUCCESS )
        return 0L;
    return (long)info.resident_size / 1024;
#else
    long rss = 0L;
    std::ifstream statm("/proc/self/statm");
    if (statm.is_open()) {
        long size, resident, share, text, lib, data, dt;
        statm >> size >> resident >> share >> text >> lib >> data >> dt;
        rss = resident * (sysconf(_SC_PAGESIZE) / 1024);
    }
    return rss;
#endif
}

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

// Phase 1 Tree Node (replaces sage::galaxy usage)
struct TreeNode {
    int snapshot;
    long long galaxy_idx;
    int merge_into_id;
    int merge_into_snapshot;
    int original_index; // Index in the tree vector before reordering
    int local_index;    // Index within the snapshot (needed for merging logic)
    
    // Tree structure fields
    int descendant = -1;
    long long global_descendant = -1;
    long long global_index = -1;
    unsigned subsize = 0;
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
// HDF5 Attribute Helper Functions
//==============================================================================

// Read string attribute from HDF5 dataset (by path)
// Returns empty string if attribute doesn't exist
std::string read_dataset_string_attribute(hpc::h5::location& loc, const std::string& dataset_path,
                                          const std::string& attr_name) {
    hid_t dset_id = H5Dopen2(loc.id(), dataset_path.c_str(), H5P_DEFAULT);
    if (dset_id < 0) return "";

    std::string result;
    if (H5Aexists(dset_id, attr_name.c_str()) > 0) {
        hid_t attr_id = H5Aopen(dset_id, attr_name.c_str(), H5P_DEFAULT);
        hid_t type_id = H5Aget_type(attr_id);
        size_t size = H5Tget_size(type_id);

        std::vector<char> buf(size + 1, 0);
        H5Aread(attr_id, type_id, buf.data());

        H5Tclose(type_id);
        H5Aclose(attr_id);
        result = std::string(buf.data());
    }

    H5Dclose(dset_id);
    return result;
}

// Write string attribute to HDF5 dataset (by path)
void write_dataset_string_attribute(hpc::h5::location& loc, const std::string& dataset_path,
                                    const std::string& attr_name, const std::string& value) {
    hid_t dset_id = H5Dopen2(loc.id(), dataset_path.c_str(), H5P_DEFAULT);
    if (dset_id < 0) return;

    hid_t str_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(str_type, value.length() + 1);
    H5Tset_strpad(str_type, H5T_STR_NULLTERM);

    hid_t space_id = H5Screate(H5S_SCALAR);
    hid_t attr_id = H5Acreate2(dset_id, attr_name.c_str(), str_type,
                               space_id, H5P_DEFAULT, H5P_DEFAULT);

    H5Awrite(attr_id, str_type, value.c_str());

    H5Aclose(attr_id);
    H5Sclose(space_id);
    H5Tclose(str_type);
    H5Dclose(dset_id);
}

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
    void _depthfirst_ordering(hpc::view<std::vector<TreeNode>> gals);
    unsigned _depthfirst_recurse(unsigned idx,
                                 hpc::view<std::vector<TreeNode>> &gals,
                                 std::vector<unsigned> &map,
                                 unsigned &dfi,
                                 std::multimap<unsigned, unsigned> const &parents);

    //==========================================================================
    // Phase 2: Add Traversal Metadata (from sageimport.cc)
    //==========================================================================
    void phase2_add_traversal_metadata();
    void _dfs_visit(Node* u, long long& counter);

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
    void _read_coords(hpc::h5::file &file,
                     std::string const &snap_name,
                     std::array<std::vector<double>, 3> &crds,
                     std::vector<unsigned long long> &tree_idxs,
                     std::vector<unsigned> &subsize);
    void _write_attributes(hpc::h5::file &file,
                          std::string const &snap_name,
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
    struct PhaseMetrics {
        std::string name;
        double duration;
        long peak_rss;
    };
    std::vector<PhaseMetrics> metrics;

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
        auto p1_start = std::chrono::high_resolution_clock::now();
        phase1_sageh5_to_depthfirst();
        _comm->barrier();
        auto p1_end = std::chrono::high_resolution_clock::now();
        metrics.push_back({"Phase 1", std::chrono::duration<double>(p1_end - p1_start).count(), get_peak_rss_kb()});

        if (_verb && _comm->rank() == 0) {
            std::cout << "  → Written: " << _depthfirst_fn << std::endl;
        }

        // Phase 2: Add Traversal Metadata
        if (_verb && _comm->rank() == 0) {
            std::cout << "\n=== Phase 2: Adding traversal metadata ===" << std::endl;
        }
        auto p2_start = std::chrono::high_resolution_clock::now();
        phase2_add_traversal_metadata();
        _comm->barrier();
        auto p2_end = std::chrono::high_resolution_clock::now();
        metrics.push_back({"Phase 2", std::chrono::duration<double>(p2_end - p2_start).count(), get_peak_rss_kb()});

        if (_verb && _comm->rank() == 0) {
            std::cout << "  → Written: " << _enhanced_fn << std::endl;
        }

        // Phase 2 complete - continue to Phase 3

        // Phase 3: Tree Order → Snapshot Order
        if (_verb && _comm->rank() == 0) {
            std::cout << "\n=== Phase 3: Tree order → Snapshot order ===" << std::endl;
        }
        auto p3_start = std::chrono::high_resolution_clock::now();
        phase3_tree_to_snapshot();
        _comm->barrier();
        auto p3_end = std::chrono::high_resolution_clock::now();
        metrics.push_back({"Phase 3", std::chrono::duration<double>(p3_end - p3_start).count(), get_peak_rss_kb()});

        if (_verb && _comm->rank() == 0) {
            std::cout << "  → Written: " << _bysnap_fn << std::endl;
        }

        // Phase 4: Build KD-Tree Index
        if (_verb && _comm->rank() == 0) {
            std::cout << "\n=== Phase 4: Building KD-tree spatial index ===" << std::endl;
        }
        auto p4_start = std::chrono::high_resolution_clock::now();
        phase4_build_kdtree_index();
        _comm->barrier();
        auto p4_end = std::chrono::high_resolution_clock::now();
        metrics.push_back({"Phase 4", std::chrono::duration<double>(p4_end - p4_start).count(), get_peak_rss_kb()});

        // Success summary
        if (_verb && _comm->rank() == 0) {
            std::cout << "\n=== Conversion Complete ===" << std::endl;
            std::cout << "Final output: " << _output_fn << std::endl;
            std::cout << "\nIntermediate files (for debugging):" << std::endl;
            std::cout << "  - " << _depthfirst_fn << std::endl;
            std::cout << "  - " << _enhanced_fn << std::endl;
            std::cout << "  - " << _bysnap_fn << std::endl;

            // Performance Report
            std::cout << "\n=== Performance Breakdown ===" << std::endl;
            std::cout << std::left << std::setw(10) << "Phase" 
                      << std::right << std::setw(15) << "Time (s)" 
                      << std::right << std::setw(15) << "Peak RSS (KB)" << std::endl;
            std::cout << std::string(40, '-') << std::endl;
            
            double total_time = 0;
            for (const auto& m : metrics) {
                std::cout << std::left << std::setw(10) << m.name 
                          << std::right << std::setw(15) << std::fixed << std::setprecision(2) << m.duration 
                          << std::right << std::setw(15) << m.peak_rss << std::endl;
                total_time += m.duration;
            }
            std::cout << std::string(40, '-') << std::endl;
            std::cout << std::left << std::setw(10) << "Total" 
                      << std::right << std::setw(15) << std::fixed << std::setprecision(2) << total_time 
                      << std::right << std::setw(15) << get_peak_rss_kb() << std::endl;
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
// Phase 1: SAGE HDF5 → Depth-First Ordered
// Based on sageh5toh5.cc lines 59-601
//==============================================================================

void sage2kdtree_application::phase1_sageh5_to_depthfirst() {
    auto log_mem = [&](const std::string& stage) {
        if (_verb >= 1 && _comm->rank() == 0) {
            long curr = get_current_rss_kb();
            long peak = get_peak_rss_kb();
            std::cout << "  [MEM] " << stage << " - Current: " << (curr/1024.0) << " MB, Peak: " << (peak/1024.0) << " MB" << std::endl;
        }
    };

    log_mem("Start Phase 1");

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

    log_mem("After scanning files");

    // Read field metadata from SAGE HDF5 (all ranks read independently)
    std::map<std::string, std::pair<std::string, std::string>> field_metadata;

    if (!input_files.empty()) {
        if (_verb >= 2 && _comm->rank() == 0) {
            std::cout << "  → Reading field metadata from SAGE HDF5..." << std::endl;
        }

        // All ranks read metadata from first file
        hpc::h5::file first_file;
        first_file.open(input_files[0].native(), H5F_ACC_RDONLY);

        // Get temporary field list to know what fields to look for
        auto temp_field_list = sage::get_galaxy_field_list_full();

        // Find first non-empty snapshot
        for (int snap = 0; snap < 64; ++snap) {
            std::ostringstream ss;
            ss << "Snap_" << snap;
            if (!first_file.has_link(ss.str())) continue;

            hpc::h5::group snap_group;
            snap_group.open(first_file, ss.str());

            // Iterate through all datasets in the snapshot group to find SAGE fields
            hsize_t num_objs = 0;
            H5Gget_num_objs(snap_group.id(), &num_objs);

            for (hsize_t i = 0; i < num_objs; ++i) {
                char name_buf[256];
                H5Gget_objname_by_idx(snap_group.id(), i, name_buf, sizeof(name_buf));
                std::string sage_field_name(name_buf);

                // Read attributes from this SAGE field
                std::string desc = read_dataset_string_attribute(snap_group, sage_field_name, "Description");
                std::string units = read_dataset_string_attribute(snap_group, sage_field_name, "Units");

                // Use SAGE field name directly (preserve CamelCase)
                // Use SAGE field name as fallback for description
                if (desc.empty()) desc = sage_field_name;

                field_metadata[sage_field_name] = {desc, units};
            }
            break;  // Only need first snapshot
        }

        if (_verb >= 2 && _comm->rank() == 0) {
            std::cout << "  → Loaded metadata for " << field_metadata.size() << " fields" << std::endl;
        }

        // Validate that all mandatory SAGE fields are present
        std::vector<std::string> available_fields;
        for (const auto& entry : field_metadata) {
            available_fields.push_back(entry.first);
        }

        std::vector<std::string> missing_fields;
        if (!tao::validate_sage_fields(available_fields, missing_fields)) {
            if (_comm->rank() == 0) {
                std::cerr << "\n❌ ERROR: SAGE HDF5 input is missing mandatory fields:" << std::endl;
                for (const auto& field : missing_fields) {
                    std::cerr << "  - " << field << std::endl;
                }
                std::cerr << "\nPlease check that your SAGE HDF5 files contain all required fields." << std::endl;
                std::cerr << "See FIELD_NAMING.md for the complete list of mandatory fields." << std::endl;
            }
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        if (_verb && _comm->rank() == 0) {
            std::cout << "  ✓ All mandatory SAGE fields present" << std::endl;
        }
    }

    // 3. Create Output File - Columnar Storage
    // Build field list dynamically from SAGE HDF5 + computed fields
    std::vector<sage::FieldInfo> field_list;

    // Add SAGE fields (read directly from HDF5 metadata)
    for (const auto& entry : field_metadata) {
        const std::string& field_name = entry.first;
        const std::string& desc = entry.second.first;
        const std::string& units = entry.second.second;

        // Determine HDF5 datatype by reading from first file
        hpc::h5::datatype field_type = hpc::h5::datatype::native_float; // default
        if (!input_files.empty()) {
            hpc::h5::file first_file(input_files[0].native(), H5F_ACC_RDONLY);
            for (int snap = 0; snap < 64; ++snap) {
                std::ostringstream ss;
                ss << "Snap_" << snap;
                if (first_file.has_link(ss.str() + "/" + field_name)) {
                    hpc::h5::group snap_group;
                    snap_group.open(first_file, ss.str());
                    hpc::h5::dataset dset = snap_group.dataset(field_name);
                    H5T_class_t type_class = dset.type_class();

                    if (type_class == H5T_FLOAT) {
                        field_type = hpc::h5::datatype::native_float;
                    } else if (type_class == H5T_INTEGER) {
                        // Use heuristic: fields with "Index" are typically long long
                        if (field_name.find("Index") != std::string::npos ||
                            field_name.find("index") != std::string::npos) {
                            field_type = hpc::h5::datatype::native_llong;
                        } else {
                            field_type = hpc::h5::datatype::native_int;
                        }
                    }
                    break;
                }
            }
        }

        field_list.push_back({field_name, field_type, 0, desc, units});
    }

    // Add computed fields (created by pipeline, not in SAGE)
    // Note: SnapNum from SAGE will be accessible via lowercase "snapnum" key in the map
    field_list.push_back({"local_index", hpc::h5::datatype::native_int, 0, "Local index within snapshot", ""});
    field_list.push_back({"global_index", hpc::h5::datatype::native_llong, 0, "Global index across all galaxies", ""});
    field_list.push_back({"descendant", hpc::h5::datatype::native_int, 0, "Local index of descendant", ""});
    field_list.push_back({"global_descendant", hpc::h5::datatype::native_llong, 0, "Global index of descendant", ""});
    // NOTE: "subsize" removed - Phase 1 doesn't write it. Phase 2 creates "subtree_count" instead.

    if (_verb >= 2 && _comm->rank() == 0) {
        std::cout << "  → Using columnar storage with " << field_list.size() << " fields" << std::endl;
    }

    hpc::h5::file out_file(_depthfirst_fn.native(), H5F_ACC_TRUNC, *_comm);

    // Create datasets for each field - NO BUFFERS, we'll write directly
    // Use shared_ptr to avoid issues with dataset handle invalidation during map operations
    std::map<std::string, std::shared_ptr<hpc::h5::dataset>> field_datasets;

    for (const auto& field : field_list) {
        // Store in map with lowercase key for lookup, but HDF5 dataset uses actual name
        std::string lookup_key = field.name;
        std::transform(lookup_key.begin(), lookup_key.end(), lookup_key.begin(), ::tolower);

        auto dataset_ptr = std::make_shared<hpc::h5::dataset>(
            out_file, "fields/" + field.name, field.type,
            hpc::h5::dataspace(total_gals)
        );

        auto result = field_datasets.emplace(lookup_key, dataset_ptr);

        if (!result.second && _verb >= 2 && _comm->rank() == 0) {
            std::cerr << "  ⚠  Duplicate lowercase key '" << lookup_key
                      << "' for field '" << field.name << "' - skipping" << std::endl;
        }

        // Write attributes to dataset
        if (!field.description.empty()) {
            write_dataset_string_attribute(out_file, "fields/" + field.name, "Description", field.description);
        }
        if (!field.units.empty()) {
            write_dataset_string_attribute(out_file, "fields/" + field.name, "Units", field.units);
        }
    }

    if (_verb >= 2 && _comm->rank() == 0) {
        std::cout << "  → Created " << field_datasets.size() << " dataset entries in map" << std::endl;
    }

    log_mem("After creating output datasets");

    // Accumulate ALL field data in memory - will write in one batch per rank
    // Use map to store all fields dynamically (key = lowercase field name)
    // REMOVED: all_int_fields, all_llong_fields, all_float_fields (now using batching)


    // Tree metadata datasets (unchanged)
    hpc::h5::dataset tree_displs_dset(out_file, "tree_displs", hpc::h5::datatype::native_ullong,
                                     hpc::h5::dataspace(total_trees + 1));
    hpc::h5::dataset tree_cnts_dset(out_file, "tree_counts", hpc::h5::datatype::native_uint,
                                   hpc::h5::dataspace(total_trees));

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
    unsigned long long current_gal_write_offset = global_gal_offset; // For writing fields

    // Memory constraint: Process trees in batches
    const size_t MAX_GALAXIES = 50000; // Adjust as needed (e.g. 20M galaxies per batch)

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

        // 1. Determine tree sizes for batching
        // We need to re-scan SAGETreeIndex to know how many galaxies are in each tree
        // This adds a pass over the file but is necessary for memory-constrained batching
        std::vector<int> tree_sizes(n_trees, 0);
        
        if (_verb >= 2 && _comm->rank() == 0) {
            std::cout << "      Scanning tree sizes for batching..." << std::endl;
        }

        for(int snap = 0; snap < 64; ++snap) {
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

        // 2. Create batches
        struct Batch {
            size_t start_tree;
            size_t end_tree;
            size_t total_gals;
        };
        std::vector<Batch> batches;
        
        size_t current_batch_gals = 0;
        size_t batch_start = 0;
        
        for(size_t t = 0; t < n_trees; ++t) {
            current_batch_gals += tree_sizes[t];
            if (current_batch_gals >= MAX_GALAXIES) {
                batches.push_back({batch_start, t + 1, current_batch_gals});
                batch_start = t + 1;
                current_batch_gals = 0;
            }
        }
        if (batch_start < n_trees) {
            batches.push_back({batch_start, n_trees, current_batch_gals});
        }

        if (_verb >= 2 && _comm->rank() == 0) {
            std::cout << "      Split into " << batches.size() << " batches (max " 
                      << (MAX_GALAXIES/1000000.0) << "M galaxies/batch)" << std::endl;
        }

        // 3. Process batches
        for (size_t b_idx = 0; b_idx < batches.size(); ++b_idx) {
            const auto& batch = batches[b_idx];
            size_t batch_n_trees = batch.end_tree - batch.start_tree;
            
            if (_verb >= 2 && _comm->rank() == 0) {
                log_mem("Start Batch " + std::to_string(b_idx+1) + "/" + std::to_string(batches.size()));
            }

            // Read data for this batch
            std::vector<std::vector<TreeNode>> trees(batch_n_trees);

            // Store ALL SAGE field data for each galaxy
            struct GalaxyFields {
                std::map<std::string, float> float_fields;
                std::map<std::string, int> int_fields;
                std::map<std::string, long long> llong_fields;
            };
            std::vector<std::vector<GalaxyFields>> trees_all_fields(batch_n_trees);

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

                // Only read fields needed for tree construction
                std::vector<long long> gal_idx;
                std::vector<int> merge_into_id, merge_into_snap;
                
                read_col("GalaxyIndex", gal_idx);
                read_col("mergeIntoID", merge_into_id); 
                read_col("mergeIntoSnapNum", merge_into_snap);

                // Read ALL SAGE fields dynamically
                std::map<std::string, std::vector<float>> snap_float_fields;
                std::map<std::string, std::vector<int>> snap_int_fields;
                std::map<std::string, std::vector<long long>> snap_llong_fields;

                // Build a map of field names to types from field_list for lookup
                std::map<std::string, hpc::h5::datatype> field_types;
                for (const auto& field : field_list) {
                    field_types[field.name] = field.type;
                }

                // Read each field from metadata
                for (const auto& entry : field_metadata) {
                    const std::string& field_name = entry.first;

                    // Skip if not in snapshot or not a SAGE field
                    if (!g.has_link(field_name)) continue;
                    if (field_types.find(field_name) == field_types.end()) continue;

                    hpc::h5::dataset dset = g.dataset(field_name);
                    hpc::h5::datatype dtype = field_types[field_name];

                    if (H5Tequal(dtype.id(), hpc::h5::datatype::native_int.id()) > 0 ||
                        H5Tequal(dtype.id(), hpc::h5::datatype::native_uint.id()) > 0) {
                        snap_int_fields[field_name].resize(n_gals_snap);
                        dset.read(snap_int_fields[field_name]);
                    } else if (H5Tequal(dtype.id(), hpc::h5::datatype::native_llong.id()) > 0 ||
                            H5Tequal(dtype.id(), hpc::h5::datatype::native_ullong.id()) > 0) {
                        snap_llong_fields[field_name].resize(n_gals_snap);
                        dset.read(snap_llong_fields[field_name]);
                    } else {
                        // Assume float for everything else
                        snap_float_fields[field_name].resize(n_gals_snap);
                        dset.read(snap_float_fields[field_name]);
                    }
                }

                for(size_t i=0; i<n_gals_snap; ++i) {
                    int tree_idx = tree_indices[i];
                    
                    // Filter: Only process galaxies belonging to trees in this batch
                    if (tree_idx >= (int)batch.start_tree && tree_idx < (int)batch.end_tree) {
                        size_t local_tree_idx = tree_idx - batch.start_tree;

                        TreeNode node;
                        node.snapshot = snap;
                        node.galaxy_idx = gal_idx[i];
                        node.merge_into_id = merge_into_id[i];
                        node.merge_into_snapshot = merge_into_snap[i];
                        node.local_index = i; // Store snapshot index for merging logic
                        
                        // Initialize others
                        node.descendant = -1;
                        node.global_descendant = -1;
                        node.global_index = -1; // Will be set later
                        node.original_index = -1; // Will be set later

                        trees[local_tree_idx].push_back(node);

                        // Also store ALL field data for this galaxy
                        size_t gal_idx_in_tree = trees[local_tree_idx].size() - 1; // Just added

                        trees_all_fields[local_tree_idx].resize(gal_idx_in_tree + 1);
                        GalaxyFields& gf = trees_all_fields[local_tree_idx][gal_idx_in_tree];

                        // Store all float fields
                        for (const auto& entry : snap_float_fields) {
                            gf.float_fields[entry.first] = entry.second[i];
                        }

                        // Store all int fields
                        for (const auto& entry : snap_int_fields) {
                            gf.int_fields[entry.first] = entry.second[i];
                        }

                        // Store all long long fields
                        for (const auto& entry : snap_llong_fields) {
                            gf.llong_fields[entry.first] = entry.second[i];
                        }
                    }
                }
            }

            // Prepare batch output buffers
            std::map<std::string, std::vector<int>> batch_int_fields;
            std::map<std::string, std::vector<long long>> batch_llong_fields;
            std::map<std::string, std::vector<float>> batch_float_fields;

            // Initialize vectors for all fields
            for (const auto& field : field_list) {
                std::string key = field.name;
                std::transform(key.begin(), key.end(), key.begin(), ::tolower);

                if (H5Tequal(field.type.id(), hpc::h5::datatype::native_int.id()) > 0 ||
                    H5Tequal(field.type.id(), hpc::h5::datatype::native_uint.id()) > 0) {
                    batch_int_fields[key].reserve(batch.total_gals);
                } else if (H5Tequal(field.type.id(), hpc::h5::datatype::native_llong.id()) > 0 ||
                        H5Tequal(field.type.id(), hpc::h5::datatype::native_ullong.id()) > 0) {
                    batch_llong_fields[key].reserve(batch.total_gals);
                } else {
                    batch_float_fields[key].reserve(batch.total_gals);
                }
            }

            // Track actual number of galaxies processed (may differ from batch.total_gals if trees are empty)
            size_t actual_batch_gals = 0;

            // Process trees in batch
            for(size_t t=0; t<batch_n_trees; ++t) {
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

                // Set original_index for reordering
                for(size_t i=0; i<tree_gals.size(); ++i) tree_gals[i].original_index = i;

                // Reorder galaxies
                _depthfirst_ordering(hpc::view<std::vector<TreeNode>>(tree_gals));

                // Accumulate ALL field data for this tree into memory (in reordered sequence)
                size_t n_gals = tree_gals.size();

                for (size_t new_idx = 0; new_idx < n_gals; ++new_idx) {
                    const auto& g = tree_gals[new_idx];
                    size_t orig_idx = g.original_index;

                    // Add computed fields
                    batch_int_fields["snapnum"].push_back(g.snapshot);
                    batch_int_fields["descendant"].push_back(g.descendant);
                    batch_int_fields["local_index"].push_back((int)new_idx); // Match original behavior (reordered index 0..N-1)
                    batch_llong_fields["global_index"].push_back(g.global_index);
                    batch_llong_fields["global_descendant"].push_back(g.global_descendant);

                    // Add ALL SAGE fields from stored data (in reordered sequence)
                    // Skip computed field names to avoid overwriting them
                    const GalaxyFields& gf = trees_all_fields[t][orig_idx];

                    // List of computed field names that should NOT be overwritten
                    static const std::set<std::string> computed_fields = {
                        "snapnum", "descendant", "local_index", "global_index", "global_descendant"
                    };

                    for (const auto& entry : gf.float_fields) {
                        std::string key = entry.first;
                        std::transform(key.begin(), key.end(), key.begin(), ::tolower);
                        batch_float_fields[key].push_back(entry.second);
                    }

                    for (const auto& entry : gf.int_fields) {
                        std::string key = entry.first;
                        std::transform(key.begin(), key.end(), key.begin(), ::tolower);
                        if (computed_fields.count(key) == 0) {  // Skip if it's a computed field
                            batch_int_fields[key].push_back(entry.second);
                        }
                    }

                    for (const auto& entry : gf.llong_fields) {
                        std::string key = entry.first;
                        std::transform(key.begin(), key.end(), key.begin(), ::tolower);
                        if (computed_fields.count(key) == 0) {  // Skip if it's a computed field
                            batch_llong_fields[key].push_back(entry.second);
                        }
                    }
                }

                // Write tree metadata
                tree_displs_buf.write(current_tree_displ);
                tree_cnts_buf.write((unsigned)tree_gals.size());

                current_gal_global_idx += tree_gals.size();
                current_tree_displ += tree_gals.size();
                actual_batch_gals += tree_gals.size();  // Track actual galaxies written
            }

            // Write batch to HDF5
            if (_verb >= 2 && _comm->rank() == 0) {
                std::cout << "      Writing batch " << (b_idx+1) << " (" << actual_batch_gals << " galaxies)..." << std::endl;
            }

            // Diagnostic: Check if we're skipping empty trees
            if (actual_batch_gals != batch.total_gals && _comm->rank() == 0) {
                std::cout << "      ⚠ WARNING: actual_batch_gals (" << actual_batch_gals
                          << ") != batch.total_gals (" << batch.total_gals << ")" << std::endl;
            }

            // Use hyperslab selection to write this batch
            hpc::h5::dataspace mem_space(actual_batch_gals);

            // Write ALL int fields
            for (const auto& entry : batch_int_fields) {
                const std::string& field_key = entry.first;
                const std::vector<int>& data = entry.second;

                auto it = field_datasets.find(field_key);
                if (it != field_datasets.end() && it->second) {
                    hpc::h5::dataspace file_space = it->second->dataspace();
                    file_space.select_hyperslab(H5S_SELECT_SET, (hsize_t)actual_batch_gals, (hsize_t)current_gal_write_offset);
                    it->second->write(data.data(), hpc::h5::datatype::native_int, mem_space, file_space);
                }
            }

            // Write ALL long long fields
            for (const auto& entry : batch_llong_fields) {
                const std::string& field_key = entry.first;
                const std::vector<long long>& data = entry.second;

                auto it = field_datasets.find(field_key);
                if (it != field_datasets.end() && it->second) {
                    hpc::h5::dataspace file_space = it->second->dataspace();
                    file_space.select_hyperslab(H5S_SELECT_SET, (hsize_t)actual_batch_gals, (hsize_t)current_gal_write_offset);
                    it->second->write(data.data(), hpc::h5::datatype::native_llong, mem_space, file_space);
                }
            }

            // Write ALL float fields
            for (const auto& entry : batch_float_fields) {
                const std::string& field_key = entry.first;
                const std::vector<float>& data = entry.second;

                auto it = field_datasets.find(field_key);
                if (it != field_datasets.end() && it->second) {
                    hpc::h5::dataspace file_space = it->second->dataspace();
                    file_space.select_hyperslab(H5S_SELECT_SET, (hsize_t)actual_batch_gals, (hsize_t)current_gal_write_offset);
                    it->second->write(data.data(), hpc::h5::datatype::native_float, mem_space, file_space);
                }
            }

            current_gal_write_offset += actual_batch_gals;

            if (_verb >= 2 && _comm->rank() == 0) {
                log_mem("End Batch " + std::to_string(b_idx+1));
            }
        } // End batch loop

        if (_verb >= 2 && _comm->rank() == 0) {
             log_mem("After processing file " + std::to_string(fidx + 1) + "/" + std::to_string(my_files.size()));
        }
    }

    // Close tree metadata buffers
    tree_displs_buf.close();
    tree_cnts_buf.close();

    _comm->barrier();

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

void sage2kdtree_application::_depthfirst_ordering(hpc::view<std::vector<TreeNode>> gals) {
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

    // Remap descendants (original_index is preserved for field lookup)
    for(unsigned ii = 0; ii < gals.size(); ++ii) {
        if(gals[ii].descendant != -1) {
            gals[ii].descendant        = map[gals[ii].descendant];
            gals[ii].global_descendant = gals[gals[ii].descendant].global_index;
        }
    }
}

unsigned sage2kdtree_application::_depthfirst_recurse(
    unsigned idx,
    hpc::view<std::vector<TreeNode>> &gals,
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
    auto log_mem = [&](const std::string& stage) {
        if (_verb >= 1 && _comm->rank() == 0) {
            long curr = get_current_rss_kb();
            long peak = get_peak_rss_kb();
            std::cout << "  [MEM] " << stage << " - Current: " << (curr/1024.0) << " MB, Peak: " << (peak/1024.0) << " MB" << std::endl;
        }
    };

    log_mem("Start Phase 2");

    int my_rank = _comm->rank();
    int n_ranks = _comm->size();

    if (_verb >= 2 && my_rank == 0) {
        std::cout << "  → Phase 2 using columnar storage" << std::endl;
    }

    // Rank 0 prepares the output file and copies metadata
    if (my_rank == 0) {
        hpc::h5::file in_file(_depthfirst_fn.native(), H5F_ACC_RDONLY);
        hpc::h5::file out_file(_enhanced_fn.native(), H5F_ACC_TRUNC);

        LOGILN("Copying metadata groups...");
        try { hpc::h5::copy(in_file, "cosmology", out_file); } catch(...) {}
        try { hpc::h5::copy(in_file, "snapshot_redshifts", out_file); } catch(...) {}
        try { hpc::h5::copy(in_file, "tree_counts", out_file); } catch(...) {}
        try { hpc::h5::copy(in_file, "tree_displs", out_file); } catch(...) {}

        in_file.close();
        out_file.close();
    }

    // Wait for rank 0
    _comm->barrier();

    // All ranks process trees in parallel
    hpc::h5::file in_file(_depthfirst_fn.native(), H5F_ACC_RDONLY, *_comm);
    hpc::h5::file out_file(_enhanced_fn.native(), H5F_ACC_RDWR, *_comm);

    // Read tree counts
    hpc::h5::dataset counts_ds(in_file, "tree_counts");
    hsize_t n_trees = counts_ds.dataspace().size();

    // Enumerate actual fields from Phase 1 output (dynamic, not hardcoded)
    hpc::h5::group fields_group;
    in_file.open_group("fields", fields_group);

    H5G_info_t group_info;
    H5Gget_info(fields_group.id(), &group_info);
    hsize_t num_fields = group_info.nlinks;

    std::vector<sage::FieldInfo> phase1_fields;
    hsize_t total_gals = 0;

    for (hsize_t i = 0; i < num_fields; i++) {
        char name_buf[256];
        H5Gget_objname_by_idx(fields_group.id(), i, name_buf, sizeof(name_buf));
        std::string field_name(name_buf);

        hpc::h5::dataset dset = fields_group.dataset(field_name);

        // Determine datatype from dataset using provided methods
        H5T_class_t type_class = dset.type_class();
        hpc::h5::datatype dset_dtype = dset.datatype();
        size_t type_size = H5Tget_size(dset_dtype.id());

        hpc::h5::datatype dtype;
        if (type_class == H5T_INTEGER) {
            if (type_size == 8) {
                dtype = hpc::h5::datatype::native_llong;
            } else {
                dtype = hpc::h5::datatype::native_int;
            }
        } else {
            dtype = hpc::h5::datatype::native_float;
        }

        phase1_fields.push_back({field_name, dtype, 0, "", ""});

        if (total_gals == 0) {
            total_gals = dset.dataspace().size();
        }
    }

    if (_verb >= 2 && my_rank == 0) {
        std::cout << "  → Found " << phase1_fields.size() << " fields in Phase 1 output" << std::endl;
    }

    // Open input columnar datasets
    std::map<std::string, hpc::h5::dataset> in_datasets;
    std::map<std::string, std::string> field_name_map;  // lowercase -> actual name

    for (const auto& field : phase1_fields) {
        in_datasets[field.name] = in_file.dataset("fields/" + field.name);

        // Create lowercase mapping for case-insensitive lookup
        std::string lower_name = field.name;
        std::transform(lower_name.begin(), lower_name.end(), lower_name.begin(), ::tolower);
        field_name_map[lower_name] = field.name;
    }

    // Helper lambda to get dataset by case-insensitive name
    auto get_dataset = [&](const char* name) -> hpc::h5::dataset& {
        std::string lower_name(name);
        std::transform(lower_name.begin(), lower_name.end(), lower_name.begin(), ::tolower);
        auto it = field_name_map.find(lower_name);
        if (it == field_name_map.end()) {
            throw std::runtime_error(std::string("Field not found: ") + name);
        }
        return in_datasets.at(it->second);
    };

    std::vector<int> tree_counts(n_trees);
    counts_ds.read(tree_counts.data(), hpc::h5::datatype::native_int);

    // Create output columnar datasets (6 original + 5 new traversal fields = 11 total)
    std::map<std::string, std::shared_ptr<hpc::h5::dataset>> out_datasets;

    // Create property list with "no fill" to avoid Lustre zeroing overhead
    hpc::h5::property_list no_fill_props;
    no_fill_props.create(H5P_DATASET_CREATE);
    H5Pset_fill_time(no_fill_props.id(), H5D_FILL_TIME_NEVER);  // Critical for Lustre!

    if (_verb >= 2 && my_rank == 0) {
        std::cout << "  → Creating contiguous datasets (no fill, no chunking)" << std::endl;
    }

    // Copy original fields
    for (const auto& field : phase1_fields) {
        out_datasets[field.name] = std::make_shared<hpc::h5::dataset>(
            out_file, "fields/" + field.name, field.type, hpc::h5::dataspace(total_gals), no_fill_props
        );

        // Copy attributes from Phase 1 input
        std::string desc = read_dataset_string_attribute(in_file, "fields/" + field.name, "Description");
        std::string units = read_dataset_string_attribute(in_file, "fields/" + field.name, "Units");

        if (!desc.empty()) write_dataset_string_attribute(out_file, "fields/" + field.name, "Description", desc);
        if (!units.empty()) write_dataset_string_attribute(out_file, "fields/" + field.name, "Units", units);
    }

    // Add 5 new traversal fields
    out_datasets["globaltreeid"] = std::make_shared<hpc::h5::dataset>(
        out_file, "fields/globaltreeid", hpc::h5::datatype::native_llong, hpc::h5::dataspace(total_gals), no_fill_props
    );
    out_datasets["breadthfirst_traversalorder"] = std::make_shared<hpc::h5::dataset>(
        out_file, "fields/breadthfirst_traversalorder", hpc::h5::datatype::native_llong, hpc::h5::dataspace(total_gals), no_fill_props
    );
    out_datasets["depthfirst_traversalorder"] = std::make_shared<hpc::h5::dataset>(
        out_file, "fields/depthfirst_traversalorder", hpc::h5::datatype::native_llong, hpc::h5::dataspace(total_gals), no_fill_props
    );
    out_datasets["subtree_count"] = std::make_shared<hpc::h5::dataset>(
        out_file, "fields/subtree_count", hpc::h5::datatype::native_llong, hpc::h5::dataspace(total_gals), no_fill_props
    );
    out_datasets["localgalaxyid"] = std::make_shared<hpc::h5::dataset>(
        out_file, "fields/localgalaxyid", hpc::h5::datatype::native_int, hpc::h5::dataspace(total_gals), no_fill_props
    );

    // Define metadata for traversal fields and write attributes
    std::map<std::string, std::pair<std::string, std::string>> traversal_metadata = {
        {"globaltreeid", {"Global tree identifier", ""}},
        {"breadthfirst_traversalorder", {"Breadth-first traversal index", ""}},
        {"depthfirst_traversalorder", {"Depth-first traversal index", ""}},
        {"subtree_count", {"Number of descendants in subtree", ""}},
        {"localgalaxyid", {"Local galaxy identifier within tree", ""}}
    };

    for (const auto& kv : traversal_metadata) {
        write_dataset_string_attribute(out_file, "fields/" + kv.first, "Description", kv.second.first);
        if (!kv.second.second.empty()) {
            write_dataset_string_attribute(out_file, "fields/" + kv.first, "Units", kv.second.second);
        }
    }

    // Calculate tree offsets
    std::vector<long long> tree_offsets(n_trees + 1, 0);
    for(size_t i=0; i<n_trees; ++i) {
        tree_offsets[i+1] = tree_offsets[i] + tree_counts[i];
    }

    log_mem("After creating output datasets");

    // First: Copy ALL Phase 1 fields in bulk (no tree-by-tree processing)
    if (_verb >= 2 && _comm->rank() == 0) {
        std::cout << "  → Copying " << phase1_fields.size() << " fields from Phase 1 output..." << std::endl;
    }

    for (const auto& field : phase1_fields) {
        const std::string& fname = field.name;
        hpc::h5::datatype dtype = field.type;

        // Allocate buffer for entire dataset (each rank handles portion)
        unsigned long long n_loc = total_gals / n_ranks;
        unsigned long long offs = my_rank * n_loc;
        if (my_rank == n_ranks - 1) {
            n_loc = total_gals - offs;
        }

        // Get actual type info directly from dataset for diagnostic
        H5T_class_t type_class = in_datasets.at(fname).type_class();
        size_t type_size = H5Tget_size(in_datasets.at(fname).datatype().id());

        if (_verb >= 2 && _comm->rank() == 0) {
            const char* type_name = "unknown";
            if (type_class == H5T_INTEGER) {
                if (type_size == 8) type_name = "llong";
                else if (type_size == 4) type_name = "int";
                else type_name = "int?";
            } else if (type_class == H5T_FLOAT) {
                if (type_size == 4) type_name = "float";
                else if (type_size == 8) type_name = "double";
                else type_name = "float?";
            }
            std::cout << "    Copying field: " << fname << " [" << type_name
                      << ", " << type_size << " bytes, " << (n_loc * type_size / 1024.0 / 1024.0)
                      << " MB] ..." << std::flush;
        }

        hpc::h5::dataspace mem_space(n_loc);
        hpc::h5::dataspace file_space(total_gals);
        file_space.select_hyperslab(H5S_SELECT_SET, (hsize_t)n_loc, (hsize_t)offs);

        // Copy field data based on type
        if (type_class == H5T_INTEGER && type_size == 4) {
            std::vector<int> field_data(n_loc);
            if (_verb >= 2 && _comm->rank() == 0) {
                std::cout << " read..." << std::flush;
            }
            in_datasets.at(fname).read(field_data.data(), dtype, mem_space, file_space);
            if (_verb >= 2 && _comm->rank() == 0) {
                std::cout << " write..." << std::flush;
            }
            out_datasets[fname]->write(field_data.data(), dtype, mem_space, file_space);
        } else if (type_class == H5T_INTEGER && type_size == 8) {
            std::vector<long long> field_data(n_loc);
            if (_verb >= 2 && _comm->rank() == 0) {
                std::cout << " read..." << std::flush;
            }
            in_datasets.at(fname).read(field_data.data(), dtype, mem_space, file_space);
            if (_verb >= 2 && _comm->rank() == 0) {
                std::cout << " write..." << std::flush;
            }
            out_datasets[fname]->write(field_data.data(), dtype, mem_space, file_space);
        } else {
            std::vector<float> field_data(n_loc);
            if (_verb >= 2 && _comm->rank() == 0) {
                std::cout << " read..." << std::flush;
            }
            in_datasets.at(fname).read(field_data.data(), dtype, mem_space, file_space);
            if (_verb >= 2 && _comm->rank() == 0) {
                std::cout << " write..." << std::flush;
            }
            out_datasets[fname]->write(field_data.data(), dtype, mem_space, file_space);
        }

        if (_verb >= 2 && _comm->rank() == 0) {
            std::cout << " DONE" << std::endl;
        }
    }

    log_mem("After copying Phase 1 fields");

    if (_verb >= 2 && _comm->rank() == 0) {
        std::cout << "  → Computing traversal metadata for trees..." << std::endl;
    }

    // Determine which trees this rank processes (round-robin)
    size_t trees_per_rank = n_trees / n_ranks;
    size_t start_tree = my_rank * trees_per_rank;
    size_t end_tree = (my_rank == n_ranks - 1) ? n_trees : (my_rank + 1) * trees_per_rank;

    // Memory constraint: Process trees in batches (same strategy as Phase 1)
    const size_t MAX_GALAXIES = 50000; // Adjust as needed for memory constraints

    // Create batches based on total galaxy count
    struct Batch {
        size_t start_tree;
        size_t end_tree;
        size_t total_gals;
    };
    std::vector<Batch> batches;

    size_t current_batch_gals = 0;
    size_t batch_start = start_tree;

    for(size_t t = start_tree; t < end_tree; ++t) {
        current_batch_gals += tree_counts[t];
        if (current_batch_gals >= MAX_GALAXIES) {
            batches.push_back({batch_start, t + 1, current_batch_gals});
            batch_start = t + 1;
            current_batch_gals = 0;
        }
    }
    if (batch_start < end_tree) {
        batches.push_back({batch_start, end_tree, current_batch_gals});
    }

    if (_verb >= 2 && _comm->rank() == 0) {
        std::cout << "  → Split into " << batches.size() << " batches (max "
                  << (MAX_GALAXIES/1000000.0) << "M galaxies/batch)" << std::endl;
    }

    // Process each batch
    for (size_t b_idx = 0; b_idx < batches.size(); ++b_idx) {
        const auto& batch = batches[b_idx];
        size_t batch_n_trees = batch.end_tree - batch.start_tree;

        if (_verb >= 2 && _comm->rank() == 0) {
            if (batches.size() > 1) {
                std::cout << "    Processing batch " << (b_idx+1) << "/" << batches.size()
                          << " (" << batch_n_trees << " trees, "
                          << (batch.total_gals/1000.0) << "K galaxies)" << std::endl;
            }
            log_mem("Start Batch " + std::to_string(b_idx+1) + "/" + std::to_string(batches.size()));
        }

        // Calculate batch extent in global galaxy indices
        long long batch_start_offset = tree_offsets[batch.start_tree];
        long long batch_total_gals = tree_offsets[batch.end_tree] - batch_start_offset;

        // Read descendant data for entire batch in one operation
        std::vector<int> descendant_batch(batch_total_gals);
        hpc::h5::dataspace batch_mem_space(batch_total_gals);
        hpc::h5::dataspace batch_file_space(total_gals);
        batch_file_space.select_hyperslab(H5S_SELECT_SET, (hsize_t)batch_total_gals, (hsize_t)batch_start_offset);
        get_dataset("descendant").read(descendant_batch.data(), hpc::h5::datatype::native_int, batch_mem_space, batch_file_space);

        // Allocate output arrays for entire batch
        std::vector<long long> globaltreeid_batch(batch_total_gals);
        std::vector<long long> bfs_order_batch(batch_total_gals);
        std::vector<long long> dfs_order_batch(batch_total_gals);
        std::vector<long long> subtree_count_batch(batch_total_gals);
        std::vector<int> localgalaxyid_batch(batch_total_gals);

        // Process each tree within the batch
        for(size_t t = batch.start_tree; t < batch.end_tree; ++t) {
            long long count = tree_counts[t];
            if (count == 0) continue;

            long long tree_start = tree_offsets[t];
            long long tree_offset_in_batch = tree_start - batch_start_offset;

            // Build tree structure using descendant links (local to this tree within batch)
            std::vector<Node> nodes(count);

            for(int i=0; i<count; ++i) {
                nodes[i].id = i;
                nodes[i].original_index = tree_start + i;

                // Link descendant (local index relative to tree start)
                int descendant_local = descendant_batch[tree_offset_in_batch + i];

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

            // Store results in batch arrays
            for(int i=0; i<count; ++i) {
                long long batch_idx = tree_offset_in_batch + i;
                globaltreeid_batch[batch_idx] = tree_start + i; // global galaxy index
                bfs_order_batch[batch_idx] = nodes[i].bfs_idx;
                dfs_order_batch[batch_idx] = nodes[i].dfs_idx;
                subtree_count_batch[batch_idx] = nodes[i].subtree_count;
                localgalaxyid_batch[batch_idx] = i; // local index within tree
            }
        }

        // Write entire batch in 5 operations (instead of 5 * n_trees operations)
        out_datasets["globaltreeid"]->write(globaltreeid_batch.data(), hpc::h5::datatype::native_llong, batch_mem_space, batch_file_space);
        out_datasets["breadthfirst_traversalorder"]->write(bfs_order_batch.data(), hpc::h5::datatype::native_llong, batch_mem_space, batch_file_space);
        out_datasets["depthfirst_traversalorder"]->write(dfs_order_batch.data(), hpc::h5::datatype::native_llong, batch_mem_space, batch_file_space);
        out_datasets["subtree_count"]->write(subtree_count_batch.data(), hpc::h5::datatype::native_llong, batch_mem_space, batch_file_space);
        out_datasets["localgalaxyid"]->write(localgalaxyid_batch.data(), hpc::h5::datatype::native_int, batch_mem_space, batch_file_space);

        if (_verb >= 2 && _comm->rank() == 0) {
            log_mem("End Batch " + std::to_string(b_idx+1));
        }
    }

    log_mem("End Phase 2");

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

//==============================================================================
// Phase 3: Tree → Snapshot
// Adapted from dstreeinit.cc tree2sageANY (lines 509-591)
//==============================================================================

void sage2kdtree_application::phase3_tree_to_snapshot() {
    if (_verb >= 1 && _comm->rank() == 0) {
        std::cout << "Phase 3: Converting tree order to snapshot order (columnar)..." << std::endl;
    }

    // Open enhanced tree file (output of Phase 2)
    hpc::h5::file tree_file(_enhanced_fn.string(), H5F_ACC_RDONLY, *_comm);

    // Enumerate actual fields from Phase 2 output (dynamic, not hardcoded)
    hpc::h5::group fields_group;
    tree_file.open_group("fields", fields_group);

    H5G_info_t group_info;
    H5Gget_info(fields_group.id(), &group_info);
    hsize_t num_fields = group_info.nlinks;

    std::vector<sage::FieldInfo> all_fields;
    std::vector<std::string> all_field_names;
    hsize_t total_gals = 0;

    if (_verb >= 2 && _comm->rank() == 0) {
        std::cout << "  → Found " << num_fields << " fields in Phase 2 output" << std::endl;
    }

    for (hsize_t i = 0; i < num_fields; i++) {
        char name_buf[256];
        H5Gget_objname_by_idx(fields_group.id(), i, name_buf, sizeof(name_buf));
        std::string field_name(name_buf);

        hpc::h5::dataset dset = fields_group.dataset(field_name);

        // Determine datatype from dataset
        H5T_class_t type_class = dset.type_class();
        hpc::h5::datatype dset_dtype = dset.datatype();
        size_t type_size = H5Tget_size(dset_dtype.id());

        hpc::h5::datatype dtype;
        if (type_class == H5T_INTEGER) {
            if (type_size == 8) {
                dtype = hpc::h5::datatype::native_llong;
            } else {
                dtype = hpc::h5::datatype::native_int;
            }
        } else {
            dtype = hpc::h5::datatype::native_float;
        }

        all_fields.push_back({field_name, dtype, 0, "", ""});
        all_field_names.push_back(field_name);

        if (total_gals == 0) {
            total_gals = dset.dataspace().size();
        }
    }

    // Open all input field datasets
    std::map<std::string, hpc::h5::dataset> in_datasets;
    for (const auto& fname : all_field_names) {
        in_datasets[fname] = tree_file.dataset("fields/" + fname);
    }

    // Create case-insensitive field name map
    std::map<std::string, std::string> field_name_map;  // lowercase -> actual name
    for (const auto& fname : all_field_names) {
        std::string lower_name = fname;
        std::transform(lower_name.begin(), lower_name.end(), lower_name.begin(), ::tolower);
        field_name_map[lower_name] = fname;
    }

    // Helper lambda to get dataset by case-insensitive name
    auto get_dataset = [&](const char* name) -> hpc::h5::dataset& {
        std::string lower_name(name);
        std::transform(lower_name.begin(), lower_name.end(), lower_name.begin(), ::tolower);
        auto it = field_name_map.find(lower_name);
        if (it == field_name_map.end()) {
            throw std::runtime_error(std::string("Field not found: ") + name);
        }
        return in_datasets.at(it->second);
    };

    // Count galaxies per snapshot by reading snapnum field
    std::vector<unsigned long long> snap_counts(64, 0);
    {
        // Read snapnum field to count galaxies per snapshot
        unsigned long long n_loc = total_gals / _comm->size();
        unsigned long long offs = _comm->rank() * n_loc;
        if (_comm->rank() == _comm->size() - 1) {
            n_loc = total_gals - offs;
        }

        std::vector<int> snapnums(n_loc);
        hpc::h5::dataspace mem_space(n_loc);
        hpc::h5::dataspace file_space = get_dataset("snapnum").dataspace();
        file_space.select_hyperslab(H5S_SELECT_SET, (hsize_t)n_loc, (hsize_t)offs);
        get_dataset("snapnum").read(snapnums.data(), hpc::h5::datatype::native_int, mem_space, file_space);

        std::vector<unsigned long long> local_counts(64, 0);
        for (int sn : snapnums) {
            if (sn >= 0 && sn < 64) {
                local_counts[sn]++;
            }
        }

        // Reduce counts across all ranks (copy first, then reduce in place)
        snap_counts = local_counts;
        _comm->all_reduce(hpc::view<std::vector<unsigned long long>>(snap_counts), MPI_SUM);
    }

    unsigned long long total_from_counts = 0;
    for (auto cnt : snap_counts) total_from_counts += cnt;

    if (_verb >= 2 && _comm->rank() == 0) {
        std::cout << "  Total galaxies: " << total_from_counts << std::endl;
        std::cout << "  Redistributing into 64 snapshots..." << std::endl;
    }

    // Create output file with 64 snapshot groups, each with all columnar datasets
    hpc::h5::file out_file(_bysnap_fn.string(), H5F_ACC_TRUNC, *_comm);

    // Map field names to their HDF5 datatypes (from dynamic field list)
    std::map<std::string, hpc::h5::datatype> field_types;
    for (const auto& f : all_fields) {
        field_types[f.name] = f.type;
    }

    // Create snapshot groups and datasets
    std::vector<std::map<std::string, hpc::h5::dataset>> snap_datasets(64);
    for (int snap = 0; snap < 64; ++snap) {
        std::ostringstream ss;
        ss << "snapshot" << std::setfill('0') << std::setw(3) << snap;
        std::string snap_name = ss.str();

        // Create group for this snapshot
        hpc::h5::group snap_group(out_file, snap_name);

        // Create columnar dataset for each field
        for (const auto& fname : all_field_names) {
            hpc::h5::dataspace dspace(snap_counts[snap]);
            snap_datasets[snap][fname] = hpc::h5::dataset(
                snap_group, fname, field_types[fname], dspace
            );

            // Copy attributes from Phase 2 input
            std::string desc = read_dataset_string_attribute(tree_file, "fields/" + fname, "Description");
            std::string units = read_dataset_string_attribute(tree_file, "fields/" + fname, "Units");

            std::string snap_path = ss.str() + "/" + fname;
            if (!desc.empty()) write_dataset_string_attribute(out_file, snap_path, "Description", desc);
            if (!units.empty()) write_dataset_string_attribute(out_file, snap_path, "Units", units);
        }
    }

    // Read all fields from tree file and redistribute to snapshots
    // Process each field separately for memory efficiency
    for (const auto& fname : all_field_names) {
        if (_verb >= 2 && _comm->rank() == 0) {
            std::cout << "  Redistributing field: " << fname << std::endl;
        }

        // Determine field datatype and element size
        hpc::h5::datatype dtype = field_types[fname];
        size_t elem_size = dtype.size();

        // Read entire field from tree file (all ranks read their portion)
        unsigned long long n_loc = total_gals / _comm->size();
        unsigned long long offs = _comm->rank() * n_loc;
        if (_comm->rank() == _comm->size() - 1) {
            n_loc = total_gals - offs;
        }

        // Allocate buffer based on datatype
        std::vector<char> field_data(n_loc * elem_size);

        hpc::h5::dataspace mem_space(n_loc);
        hpc::h5::dataspace file_space = in_datasets[fname].dataspace();
        file_space.select_hyperslab(H5S_SELECT_SET, (hsize_t)n_loc, (hsize_t)offs);
        in_datasets[fname].read(field_data.data(), dtype, mem_space, file_space);

        // Also need to read snapnum to know where each galaxy goes
        std::vector<int> snapnums(n_loc);
        hpc::h5::dataspace snapnum_file_space = get_dataset("snapnum").dataspace();
        snapnum_file_space.select_hyperslab(H5S_SELECT_SET, (hsize_t)n_loc, (hsize_t)offs);
        get_dataset("snapnum").read(snapnums.data(), hpc::h5::datatype::native_int, mem_space, snapnum_file_space);

        // Build snapshot buffers for this field
        std::vector<std::vector<char>> snap_buffers(64);
        for (int snap = 0; snap < 64; ++snap) {
            snap_buffers[snap].reserve(snap_counts[snap] * elem_size / _comm->size());
        }

        // Distribute field data to snapshot buffers
        for (unsigned long long i = 0; i < n_loc; ++i) {
            int snap = snapnums[i];
            if (snap >= 0 && snap < 64) {
                char* elem_ptr = field_data.data() + i * elem_size;
                snap_buffers[snap].insert(
                    snap_buffers[snap].end(),
                    elem_ptr,
                    elem_ptr + elem_size
                );
            }
        }

        // Write each snapshot's portion of this field
        // Need to use collective write across ranks
        for (int snap = 0; snap < 64; ++snap) {
            unsigned long long local_count = snap_buffers[snap].size() / elem_size;

            // Gather counts from all ranks to determine offsets
            std::vector<unsigned long long> all_counts(_comm->size());
            _comm->all_gather(local_count, all_counts);

            unsigned long long my_offset = 0;
            for (int r = 0; r < _comm->rank(); ++r) {
                my_offset += all_counts[r];
            }

            if (local_count > 0) {
                hpc::h5::dataspace write_mem_space(local_count);
                hpc::h5::dataspace write_file_space = snap_datasets[snap][fname].dataspace();
                write_file_space.select_hyperslab(H5S_SELECT_SET, (hsize_t)local_count, (hsize_t)my_offset);
                snap_datasets[snap][fname].write(snap_buffers[snap].data(), dtype, write_mem_space, write_file_space);
            } else {
                // Still need to participate in collective write with empty selection
                hpc::h5::dataspace write_mem_space(0);
                hpc::h5::dataspace write_file_space = snap_datasets[snap][fname].dataspace();
                write_file_space.select_none();
                snap_datasets[snap][fname].write(nullptr, dtype, write_mem_space, write_file_space);
            }
        }
    }

    _comm->barrier();

    // Copy cosmology and snapshot_redshifts metadata
    if (_comm->rank() == 0) {
        try { hpc::h5::copy(tree_file, "cosmology", out_file); } catch(...) {}
        try { hpc::h5::copy(tree_file, "snapshot_redshifts", out_file); } catch(...) {}
    }

    _comm->barrier();

    if (_verb >= 1 && _comm->rank() == 0) {
        std::cout << "Phase 3 complete: " << _bysnap_fn << std::endl;
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
            // Phase 3 output: snapshot groups contain columnar field datasets
            // Try to find Posx field (case-sensitive first, then case-insensitive)
            std::string field_path;
            if (snap_file.has_link(ss.str() + "/Posx")) {
                field_path = ss.str() + "/Posx";
            } else if (snap_file.has_link(ss.str() + "/posx")) {
                field_path = ss.str() + "/posx";
            } else {
                // If Posx not found, open the group and use the first field
                hpc::h5::group snap_group;
                snap_group.open(snap_file, ss.str());
                hsize_t num_objs = 0;
                H5Gget_num_objs(snap_group.id(), &num_objs);
                if (num_objs > 0) {
                    char name_buf[256];
                    H5Gget_objname_by_idx(snap_group.id(), 0, name_buf, sizeof(name_buf));
                    field_path = ss.str() + "/" + std::string(name_buf);
                }
            }

            if (!field_path.empty() && snap_file.has_link(field_path)) {
                hpc::h5::dataset ds = snap_file.dataset(field_path);
                hpc::h5::dataspace sp = ds.dataspace();
                std::vector<hsize_t> dims_vec(1);
                sp.simple_extent_dims<std::vector<hsize_t>>(dims_vec);
                snap_counts[snap] = dims_vec[0];
            } else {
                snap_counts[snap] = 0;
            }
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

    // Create lightcone/data dataset - HDF5 will auto-create "lightcone" group
    if (_verb >= 2 && _comm->rank() == 0) {
        std::cout << "  Creating lightcone/data dataset..." << std::endl;
    }
    hpc::h5::dataspace lc_space(total_gals);
    hpc::h5::dataset lc_data(out_file, "lightcone/data", lc_type, lc_space, hpc::h5::property_list());

    // Add attributes describing the compound structure
    if (_comm->rank() == 0) {
        write_dataset_string_attribute(out_file, "lightcone/data", "Description",
            "Compound dataset with galaxy coordinates and metadata");
        write_dataset_string_attribute(out_file, "lightcone/data", "x_units", "Mpc/h");
        write_dataset_string_attribute(out_file, "lightcone/data", "y_units", "Mpc/h");
        write_dataset_string_attribute(out_file, "lightcone/data", "z_units", "Mpc/h");
    }

    if (_verb >= 2 && _comm->rank() == 0) {
        std::cout << "  lightcone/data created successfully" << std::endl;
    }

    // Discover fields and create data/* datasets - HDF5 will auto-create "data" group
    if (_comm->rank() == 0 && total_gals > 0) {
        if (_verb >= 2) {
            std::cout << "  Discovering fields and creating data/* datasets..." << std::endl;
        }
        for (int snap = 0; snap < n_snapshots; ++snap) {
            std::ostringstream ss;
            ss << "snapshot" << std::setfill('0') << std::setw(3) << snap;

            if (snap_counts[snap] > 0 && snap_file.has_link(ss.str())) {
                if (_verb >= 2) {
                    std::cout << "    Opening snapshot group: " << ss.str() << std::endl;
                }
                // Open snapshot group and iterate through field datasets
                hpc::h5::group snap_group;
                snap_group.open(snap_file, ss.str());
                if (_verb >= 2) {
                    std::cout << "    Snapshot group opened successfully" << std::endl;
                }
                hsize_t num_objs = 0;
                H5Gget_num_objs(snap_group.id(), &num_objs);

                for (hsize_t ii = 0; ii < num_objs; ++ii) {
                    char name_buf[256];
                    H5Gget_objname_by_idx(snap_group.id(), ii, name_buf, sizeof(name_buf));
                    std::string field_name(name_buf);
                    std::string field_path = ss.str() + "/" + field_name;

                    // Open the field dataset to get its type
                    hpc::h5::dataset field_ds = snap_file.dataset(field_path);
                    hpc::h5::datatype field_type = field_ds.datatype();

                    // Create data/* dataset - HDF5 auto-creates "data" group (preserve original case)
                    std::string out_name = field_name;
                    hpc::h5::dataspace field_space(total_gals);
                    hpc::h5::dataset out_ds(out_file, "data/" + out_name, field_type, field_space, hpc::h5::property_list());

                    // Copy attributes from Phase 3 input
                    std::string desc = read_dataset_string_attribute(snap_file, field_path, "Description");
                    std::string units = read_dataset_string_attribute(snap_file, field_path, "Units");

                    if (!desc.empty()) write_dataset_string_attribute(out_file, "data/" + out_name, "Description", desc);
                    if (!units.empty()) write_dataset_string_attribute(out_file, "data/" + out_name, "Units", units);
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

    // Read coordinates and metadata from columnar storage
    std::array<std::vector<double>, 3> crds;
    std::vector<unsigned long long> tree_idxs;
    std::vector<unsigned> subsize;

    _read_coords(file, name, crds, tree_idxs, subsize);

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

    // Write all other attributes to data/* datasets (from columnar storage)
    _write_attributes(file, name, out_file, idxs, displ);

    // Write KD-tree structure
    _write_kdtree(out_file, name, kdt, bp, crds, displ);

    displ += n_gals;
}

void sage2kdtree_application::_read_coords(
    hpc::h5::file &file,
    std::string const &snap_name,
    std::array<std::vector<double>, 3> &crds,
    std::vector<unsigned long long> &tree_idxs,
    std::vector<unsigned> &subsize) {

    // Helper to find field with case-insensitive fallback
    auto find_field = [&](const std::string& base_name, const std::string& alt_name) -> std::string {
        std::string path1 = snap_name + "/" + base_name;
        std::string path2 = snap_name + "/" + alt_name;
        if (file.has_link(path1)) return path1;
        if (file.has_link(path2)) return path2;
        return "";
    };

    // Try CamelCase first, then lowercase
    std::string posx_path = find_field("Posx", "posx");
    std::string posy_path = find_field("Posy", "posy");
    std::string posz_path = find_field("Posz", "posz");
    std::string globaltreeid_path = find_field("globaltreeid", "globaltreeid");  // already lowercase
    std::string subtree_count_path = find_field("subtree_count", "subtree_count");  // already lowercase

    // Get galaxy count from posx dataset
    if (posx_path.empty()) {
        crds[0].clear();
        crds[1].clear();
        crds[2].clear();
        tree_idxs.clear();
        subsize.clear();
        return;
    }

    hpc::h5::dataset posx_ds = file.dataset(posx_path);
    hpc::h5::dataspace sp = posx_ds.dataspace();
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

    // Resize output vectors
    crds[0].resize(n_gals);
    crds[1].resize(n_gals);
    crds[2].resize(n_gals);
    tree_idxs.resize(n_gals);
    subsize.resize(n_gals);

    // Read columnar fields (coordinates are float in file, double in memory)
    std::vector<float> posx_float(n_gals);
    std::vector<float> posy_float(n_gals);
    std::vector<float> posz_float(n_gals);

    hpc::h5::dataspace mem_space(n_gals);

    posx_ds.read(posx_float.data(), hpc::h5::datatype::native_float, mem_space, sp);
    file.dataset(posy_path).read(posy_float.data(), hpc::h5::datatype::native_float, mem_space, sp);
    file.dataset(posz_path).read(posz_float.data(), hpc::h5::datatype::native_float, mem_space, sp);

    // Convert float to double
    for (unsigned long long ii = 0; ii < n_gals; ++ii) {
        crds[0][ii] = static_cast<double>(posx_float[ii]);
        crds[1][ii] = static_cast<double>(posy_float[ii]);
        crds[2][ii] = static_cast<double>(posz_float[ii]);
    }

    // Read tree indices and subsize (both stored as long long in Phase 3)
    std::vector<unsigned long long> globaltreeid_data(n_gals);
    std::vector<unsigned long long> subtree_count_data(n_gals);

    file.dataset(globaltreeid_path).read(globaltreeid_data.data(), hpc::h5::datatype::native_llong, mem_space, sp);
    file.dataset(subtree_count_path).read(subtree_count_data.data(), hpc::h5::datatype::native_llong, mem_space, sp);

    // Copy to output arrays
    for (unsigned long long ii = 0; ii < n_gals; ++ii) {
        tree_idxs[ii] = globaltreeid_data[ii];
        subsize[ii] = static_cast<unsigned>(subtree_count_data[ii]);
    }
}

void sage2kdtree_application::_write_attributes(
    hpc::h5::file &file,
    std::string const &snap_name,
    hpc::h5::file &out_file,
    std::vector<unsigned long long> const &idxs,
    unsigned long long displ) {

    // Read from columnar storage: snapshot group contains separate field datasets
    unsigned long long n_gals = idxs.size();
    if (n_gals == 0) return;

    // Open snapshot group and iterate through all field datasets
    hpc::h5::group snap_group;
    snap_group.open(file, snap_name);
    hsize_t num_objs = 0;
    H5Gget_num_objs(snap_group.id(), &num_objs);

    for (hsize_t ii = 0; ii < num_objs; ++ii) {
        char name_buf[256];
        H5Gget_objname_by_idx(snap_group.id(), ii, name_buf, sizeof(name_buf));
        std::string field_name(name_buf);
        std::string field_path = snap_name + "/" + field_name;

        // Open the field dataset
        hpc::h5::dataset field_ds = file.dataset(field_path);
        hpc::h5::datatype field_type = field_ds.datatype();
        size_t elem_size = field_type.size();

        // Read entire field
        std::vector<char> field_data(n_gals * elem_size);
        hpc::h5::dataspace mem_space(n_gals);
        hpc::h5::dataspace file_space = field_ds.dataspace();
        field_ds.read(field_data.data(), field_type, mem_space, file_space);

        // Permute to KD-tree spatial order
        std::vector<char> permuted_data(n_gals * elem_size);
        for (unsigned long long jj = 0; jj < n_gals; ++jj) {
            unsigned long long src_idx = idxs[jj];
            memcpy(permuted_data.data() + jj * elem_size,
                   field_data.data() + src_idx * elem_size,
                   elem_size);
        }

        // Write to /data/{field_name} (preserve original case)
        std::string out_name = field_name;
        hpc::h5::dataset out_ds(out_file, "data/" + out_name);
        hpc::h5::dataspace out_mem_space(n_gals);
        hpc::h5::dataspace out_file_space(out_ds);
        out_file_space.select_range(displ, displ + n_gals);

        hpc::h5::property_list props(H5P_DATASET_XFER);
        props.set_preserve();
        out_ds.write(permuted_data.data(), field_type, out_mem_space, out_file_space, hpc::mpi::comm::self, props);
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
    // Open tree file to extract SED data (columnar storage)
    hpc::h5::file tree_file(tree_fn.string(), H5F_ACC_RDONLY);

    // For columnar storage, SED fields are not available in the Phase 2 output
    // Skip SED data creation for now
    // TODO: In future, could read SED fields directly from SAGE HDF5 input
    return;

    /* Original compound-type code - disabled for columnar storage
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
    */
}

//==============================================================================
// Main Entry Point
//==============================================================================

#define HPC_APP_CLASS sage2kdtree_application
#include <libhpc/mpi/main.hh>
