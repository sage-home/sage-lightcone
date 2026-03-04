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

#include <algorithm>
#include <boost/range/algorithm/fill.hpp>
#include <cctype>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <libhpc/libhpc.hh>
#include <libhpc/system/type_traits.hh>
#include <libtao/base/batch.hh>
#include <libtao/base/mandatory_fields.hh>
#include <libtao/base/sage.hh>
#include <libtao/base/utils.hh>
#include <list>
#include <map>
#include <memory>
// #include <pugixml.hpp>
#include <set>
#include <string>
#include <sys/resource.h>
#include <unordered_map>
#include <vector>

#ifdef __APPLE__
#include <mach/mach.h>
#else
#include <ios>
#include <unistd.h>
#endif

using namespace ::tao;
using namespace ::hpc;
using namespace ::sage;

// Helper for memory reporting
long get_peak_rss_kb()
{
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
#ifdef __APPLE__
    return usage.ru_maxrss / 1024; // macOS reports bytes
#else
    return usage.ru_maxrss; // Linux reports KB
#endif
}

long get_current_rss_kb()
{
#ifdef __APPLE__
    struct mach_task_basic_info info;
    mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&info, &infoCount) !=
        KERN_SUCCESS)
        return 0L;
    return (long)info.resident_size / 1024;
#else
    long rss = 0L;
    std::ifstream statm("/proc/self/statm");
    if (statm.is_open())
    {
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
struct Node
{
    size_t id;
    long long original_index;
    Node* descendant = nullptr;
    std::vector<Node*> progenitors;

    long long bfs_idx = 0;
    long long dfs_idx = 0;
    long long subtree_count = 0;
};

// Phase 1 Tree Node (replaces sage::galaxy usage)
struct TreeNode
{
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
};

// From dstreeinit.cc - Functor for KD-tree spatial reordering
struct data_permuter
{
    data_permuter()
        : crds(nullptr)
        , tree_idxs(nullptr)
        , idxs(nullptr)
    {
    }

    data_permuter(std::array<std::vector<double>, 3>& crds,
                  std::vector<unsigned long long>& tree_idxs, std::vector<unsigned long long>& idxs)
        : crds(&crds)
        , tree_idxs(&tree_idxs)
        , idxs(&idxs)
    {
    }

    void operator()(hpc::mpi::balanced_partition const& part)
    {
        for (unsigned ii = 0; ii < 3; ++ii)
            part.transfer((*crds)[ii]);
        part.transfer(*tree_idxs);
        part.transfer(*idxs);
    }

    std::array<std::vector<double>, 3>* crds;
    std::vector<unsigned long long>* tree_idxs;
    std::vector<unsigned long long>* idxs;
};

// From dstreeinit.cc - SED data structure
struct sed_data_t
{
    int descendant;
    int snapshot;
    int local_index;
    int merge_type;
    double dt;
    double disk_sfr;
    double bulge_sfr;
    double disk_sfr_z;
    double bulge_sfr_z;
};

// From sageimport.cc - Field metadata
struct SageField
{
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

class PipelineException : public std::runtime_error
{
public:
    PipelineException(int phase, const std::string& msg)
        : std::runtime_error("Phase " + std::to_string(phase) + ": " + msg)
        , _phase(phase)
    {
    }
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
                                          const std::string& attr_name)
{
    hid_t dset_id = H5Dopen2(loc.id(), dataset_path.c_str(), H5P_DEFAULT);
    if (dset_id < 0)
        return "";

    std::string result;
    if (H5Aexists(dset_id, attr_name.c_str()) > 0)
    {
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
                                    const std::string& attr_name, const std::string& value)
{
    hid_t dset_id = H5Dopen2(loc.id(), dataset_path.c_str(), H5P_DEFAULT);
    if (dset_id < 0)
        return;

    hid_t str_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(str_type, value.length() + 1);
    H5Tset_strpad(str_type, H5T_STR_NULLTERM);

    hid_t space_id = H5Screate(H5S_SCALAR);
    hid_t attr_id =
        H5Acreate2(dset_id, attr_name.c_str(), str_type, space_id, H5P_DEFAULT, H5P_DEFAULT);

    H5Awrite(attr_id, str_type, value.c_str());

    H5Aclose(attr_id);
    H5Sclose(space_id);
    H5Tclose(str_type);
    H5Dclose(dset_id);
}

//==============================================================================
// Main Application Class
//==============================================================================

class sage2kdtree_application : public hpc::mpi::application
{
public:
    sage2kdtree_application(int argc, char* argv[]);
    void operator()();

private:
    //==========================================================================
    // Phase 1 (Direct): SAGE HDF5 → Snapshot Organized (Bypassing Tree
    // Reordering)
    //==========================================================================
    void phase1_direct_sage_to_snapshot();
    void _load_param(hpc::fs::path const& fn);
    void _load_redshifts(hpc::fs::path const& fn);
    void _depthfirst_ordering(hpc::view<std::vector<TreeNode>> gals);
    void _depthfirst_recurse(unsigned idx, hpc::view<std::vector<TreeNode>>& gals,
                             std::vector<unsigned>& map, unsigned& dfi,
                             std::multimap<unsigned, unsigned> const& parents);

    //==========================================================================
    // Phase 3: Tree → Snapshot (from dstreeinit.cc tree2sage)
    //==========================================================================
    void phase3_tree_to_snapshot();

    //==========================================================================
    // Phase 4: Build KD-Tree Index (from dstreeinit.cc init)
    //==========================================================================
    void phase4_build_kdtree_index();
    void _process_snapshot(hpc::h5::file& file, std::string const& name, hpc::h5::file& out_file,
                           hpc::h5::dataset& data, unsigned long long& displ);
    void _copy_header_to_output();
    void _read_coords(hpc::h5::file& file, std::string const& snap_name,
                      std::array<std::vector<double>, 3>& crds,
                      std::vector<unsigned long long>& tree_idxs);
    void _write_attributes(hpc::h5::file& file, std::string const& snap_name,
                           hpc::h5::file& out_file, std::vector<unsigned long long> const& idxs,
                           unsigned long long displ);
    void _write_kdtree(hpc::h5::file& file, std::string const& snap_name, hpc::kdtree<> const& kdt,
                       hpc::binary_partitioner<std::array<std::vector<double>, 3>::iterator,
                                               data_permuter> const& bp,
                       std::array<std::vector<double>, 3> const& crds, unsigned long long displ);
    void _write_empty_kdtree(hpc::h5::file& file, int snap_num);

    //==========================================================================
    // Utility Methods
    //==========================================================================
    void _validate_inputs(); // Validate all inputs before processing (Priority 1.1)
    void _setup_intermediate_paths();

    //==========================================================================
    // Member Variables
    //==========================================================================

    // Command-line arguments
    hpc::fs::path _sage_dir;  // Input: SAGE HDF5 directory
    hpc::fs::path _param_fn;  // Input: SAGE parameter file
    hpc::fs::path _alist_fn;  // Input: expansion factor list
    hpc::fs::path _output_fn; // Output: final KD-tree HDF5

    // Intermediate file paths (for debugging)
    hpc::fs::path _depthfirst_fn; // Phase 1 output
    hpc::fs::path _enhanced_fn;   // Phase 2 output
    hpc::fs::path _bysnap_fn;     // Phase 3 output

    unsigned _ppc; // Particles per cell for KD-tree
    int _verb;     // Verbosity level

    // Cosmology parameters (loaded in Phase 1)
    double _box_size;
    double _hubble;
    double _omega_l;
    double _omega_m;
    std::vector<double> _redshifts; // Redshift per snapshot

    // Field metadata (used in Phases 2-4)
    std::vector<SageField> _fields;

    // MPI communicator
    hpc::mpi::comm* _comm;
};

//==============================================================================
// Constructor: Parse command-line arguments
//==============================================================================

sage2kdtree_application::sage2kdtree_application(int argc, char* argv[])
    : hpc::mpi::application(argc, argv)
    , _comm(const_cast<hpc::mpi::comm*>(&hpc::mpi::comm::world))
    , _box_size(0)
    , _hubble(0)
    , _omega_l(0)
    , _omega_m(0)
    , _ppc(1000)
    , _verb(1)
{

    // Setup command-line options
    options().add_options()("sage,s", hpc::po::value<hpc::fs::path>(&_sage_dir)->required(),
                            "SAGE HDF5 output directory")(
        "param,p", hpc::po::value<hpc::fs::path>(&_param_fn)->required(),
        "SAGE parameter file")("alist,a", hpc::po::value<hpc::fs::path>(&_alist_fn)->required(),
                               "SAGE expansion factor list file")(
        "output,o", hpc::po::value<hpc::fs::path>(&_output_fn)->required(),
        "Output KD-tree HDF5 file")("ppc", hpc::po::value<unsigned>(&_ppc)->default_value(1000),
                                    "Particles per cell for KD-tree")(
        "verbose,v", hpc::po::value<int>(&_verb)->default_value(1),
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
    if (_verb > 0)
    {
        hpc::log::levels_type lvl;
        if (_verb == 1)
            lvl = hpc::log::info;
        else if (_verb == 2)
            lvl = hpc::log::debug;
        else if (_verb == 3)
            lvl = hpc::log::trivial;

        if (_comm->size() > 1)
            LOG_PUSH(new hpc::mpi::logger("sage2kdtree.log.", lvl));
        else
            LOG_PUSH(new hpc::log::stdout(lvl));
    }
}

//==============================================================================
// Main Operator: Execute four-phase pipeline
//==============================================================================

void sage2kdtree_application::phase1_direct_sage_to_snapshot()
{
    if (_verb && _comm->rank() == 0)
    {
        std::cout << "=== Phase 1 (Direct): Aggregating SAGE snapshots ===" << std::endl;
    }

    _load_param(_param_fn);
    _load_redshifts(_alist_fn);

    std::vector<hpc::fs::path> sage_files;
    if (hpc::fs::is_directory(_sage_dir))
    {
        // Collect files matching SAGE output pattern: model_N.hdf5 where N >= 0
        // This whitelist approach prevents processing of:
        // - Master link files (model.hdf5)
        // - Intermediate output files (*-bysnap.h5, *-kdtree.h5, etc.)
        // - Any other HDF5 files in the directory
        for (const auto& entry : hpc::fs::directory_iterator(_sage_dir))
        {
            std::string filename = entry.path().filename().string();
            std::string extension = entry.path().extension().string();

            // Must have .hdf5 or .h5 extension
            if (extension != ".hdf5" && extension != ".h5")
                continue;

            // Must match pattern: model_N.hdf5 where N is a non-negative integer
            // Examples: model_0.hdf5, model_1.hdf5, model_42.hdf5
            std::string stem = entry.path().stem().string();
            if (stem.length() >= 7 && stem.substr(0, 6) == "model_")
            {
                // Check if remaining characters are all digits
                std::string suffix = stem.substr(6);
                bool all_digits =
                    !suffix.empty() && std::all_of(suffix.begin(), suffix.end(), ::isdigit);

                if (all_digits)
                {
                    sage_files.push_back(entry.path());
                }
                else if (_verb >= 2 && _comm->rank() == 0)
                {
                    std::cout << "  → Skipping non-SAGE file: \"" << filename
                              << "\" (expected model_N.hdf5 format)" << std::endl;
                }
            }
            else if (_verb >= 2 && _comm->rank() == 0)
            {
                // Verbose logging for filtered files
                if (stem == "model")
                {
                    std::cout << "  → Skipping master file: \"" << filename << "\"" << std::endl;
                }
                else
                {
                    std::cout << "  → Skipping non-SAGE file: \"" << filename
                              << "\" (expected model_N.hdf5 format)" << std::endl;
                }
            }
        }
    }
    else
    {
        sage_files.push_back(_sage_dir);
    }

    std::sort(sage_files.begin(), sage_files.end());

    int n_files = sage_files.size();
    int size = _comm->size();
    int rank = _comm->rank();

    int files_per_rank = n_files / size;
    int remainder = n_files % size;
    int my_start = rank * files_per_rank + std::min(rank, remainder);
    int my_count = files_per_rank + (rank < remainder ? 1 : 0);
    int my_end = my_start + my_count;

    if (_verb && rank == 0)
    {
        std::cout << "  → Found " << n_files << " input files" << std::endl;
        std::cout << "  → Writing to: " << _bysnap_fn << std::endl;
    }

    hpc::h5::file out_file(_bysnap_fn.string(), H5F_ACC_TRUNC, *_comm);

    // Write snapshot_redshifts (Rank 0 only)
    if (rank == 0)
    {
        hsize_t n_snaps = _redshifts.size();
        hpc::h5::dataspace z_space(n_snaps);
        hpc::h5::dataset z_dset(out_file, "snapshot_redshifts", hpc::h5::datatype::native_double,
                                z_space, hpc::h5::property_list());
        z_dset.write(_redshifts);

        // Create cosmology group and write datasets (required by kdtree_backend)
        // kdtree_backend expects: cosmology/box_size, cosmology/hubble,
        // cosmology/omega_m, cosmology/omega_l as scalar datasets
        hpc::h5::group cosmo_group;
        cosmo_group.create(out_file, "cosmology");

        out_file.write<double>("cosmology/box_size", _box_size);
        out_file.write<double>("cosmology/hubble", _hubble);
        out_file.write<double>("cosmology/omega_m", _omega_m);
        out_file.write<double>("cosmology/omega_l", _omega_l);

        // Copy other attributes from Header/Simulation if available
        // These are copied as attributes on the 'cosmology' group itself
        if (!sage_files.empty())
        {
            try
            {
                hpc::h5::file sage_input(sage_files[0].string(), H5F_ACC_RDONLY);
                if (sage_input.has_link("Header/Simulation"))
                {
                    hpc::h5::group sage_sim = sage_input.group("Header/Simulation");
                    hid_t gid = sage_sim.id();

                    H5O_info_t oinfo;
#if H5Oget_info_vers < 3
                    H5Oget_info(gid, &oinfo);
#else
                    H5Oget_info(gid, &oinfo, H5O_INFO_NUM_ATTRS);
#endif

                    for (hsize_t i = 0; i < oinfo.num_attrs; i++)
                    {
                        const size_t MAX_NAME_LEN = 256;
                        char attr_name[MAX_NAME_LEN];
                        H5Aget_name_by_idx(gid, ".", H5_INDEX_NAME, H5_ITER_NATIVE, i, attr_name,
                                           MAX_NAME_LEN, H5P_DEFAULT);
                        std::string name(attr_name);

                        // transform to lower for comparison
                        std::string name_lower = name;
                        std::transform(name_lower.begin(), name_lower.end(), name_lower.begin(),
                                       ::tolower);

                        // Skip already handled attributes (these are stored as datasets in
                        // the cosmology group)
                        if (name_lower == "boxsize" || name_lower == "box_size" ||
                            name_lower == "hubble_h" || name_lower == "hubbleparam" ||
                            name_lower == "hubble" || name_lower == "omega_m" ||
                            name_lower == "omega" || name_lower == "omega_matter" ||
                            name_lower == "omega_l" || name_lower == "omegalambda" ||
                            name_lower == "omega_lambda")
                        {
                            continue;
                        }

                        // Copy attribute
                        hid_t attr_id = H5Aopen_by_idx(gid, ".", H5_INDEX_NAME, H5_ITER_NATIVE, i,
                                                       H5P_DEFAULT, H5P_DEFAULT);
                        if (attr_id < 0)
                            continue;

                        hid_t type_id = H5Aget_type(attr_id);
                        hid_t space_id = H5Aget_space(attr_id);

                        // Create attribute on output group 'cosmology'
                        // Use the original name
                        hid_t out_attr_id = H5Acreate2(cosmo_group.id(), name.c_str(), type_id,
                                                       space_id, H5P_DEFAULT, H5P_DEFAULT);

                        if (out_attr_id >= 0)
                        {
                            // Read and Write
                            // For reading, we need a buffer large enough.
                            // H5Tget_size gives the size in bytes of the type.
                            // H5Sget_simple_extent_npoints gives number of elements.
                            size_t type_size = H5Tget_size(type_id);
                            hssize_t n_points = H5Sget_simple_extent_npoints(space_id);
                            if (n_points < 1)
                                n_points = 1;

                            std::vector<char> data(type_size * n_points);

                            H5Aread(attr_id, type_id, data.data());
                            H5Awrite(out_attr_id, type_id, data.data());
                            H5Aclose(out_attr_id);
                        }

                        H5Sclose(space_id);
                        H5Tclose(type_id);
                        H5Aclose(attr_id);

                        if (_verb >= 2 && rank == 0)
                            std::cout << "    Copied cosmology attribute: " << name << std::endl;
                    }
                }
            }
            catch (std::exception& e)
            {
                if (_verb >= 1 && rank == 0)
                    std::cout << "Warning: Failed to copy extra cosmology attributes: " << e.what()
                              << std::endl;
            }
        }
    }

    for (size_t snap = 0; snap < _redshifts.size(); ++snap)
    {
        std::stringstream ss_in;
        ss_in << "Snap_" << snap;
        std::string group_in = ss_in.str();

        std::stringstream ss_out;
        ss_out << "snapshot" << std::setfill('0') << std::setw(3) << snap;
        std::string group_out = ss_out.str();

        unsigned long long local_gals = 0;
        std::vector<std::shared_ptr<hpc::h5::file>> open_files;
        std::vector<int> valid_file_indices;
        std::string representative_field; // Track which field was used for counting

        for (int i = my_start; i < my_end; ++i)
        {
            try
            {
                auto f = std::make_shared<hpc::h5::file>(sage_files[i].string(), H5F_ACC_RDONLY);
                if (f->has_link(group_in))
                {
                    hpc::h5::group g;
                    g.open(*f, group_in);
                    std::string count_ds_name;
                    if (g.has_link("GalaxyIndex"))
                        count_ds_name = "GalaxyIndex";
                    else if (g.has_link("Posx"))
                        count_ds_name = "Posx";
                    else
                    {
                        hsize_t n = 0;
                        H5Gget_num_objs(g.id(), &n);
                        if (n > 0)
                        {
                            char name[256];
                            H5Gget_objname_by_idx(g.id(), 0, name, 256);
                            count_ds_name = name;
                        }
                    }

                    if (!count_ds_name.empty())
                    {
                        hpc::h5::dataset d = g.dataset(count_ds_name);
                        hsize_t this_count = d.dataspace().size();
                        local_gals += this_count;
                        open_files.push_back(f);
                        valid_file_indices.push_back(i);

                        if (representative_field.empty())
                        {
                            representative_field = count_ds_name;
                        }

                        if (_verb >= 3 && rank == 0)
                        {
                            std::cout << "      File " << i << ": counted " << this_count
                                      << " galaxies from field '" << count_ds_name << "'"
                                      << std::endl;
                        }
                    }
                }
            }
            catch (...)
            {
            }
        }

        unsigned long long total_gals = 0;
        unsigned long long my_offset = 0;

        MPI_Allreduce(&local_gals, &total_gals, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                      _comm->mpi_comm());
        MPI_Exscan(&local_gals, &my_offset, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, _comm->mpi_comm());
        if (rank == 0)
            my_offset = 0;

        if (total_gals == 0)
            continue;

        if (_verb >= 2 && rank == 0)
        {
            std::cout << "    Processing " << group_out << " (" << total_gals << " gals)";
            if (!representative_field.empty())
            {
                std::cout << " [counted from '" << representative_field << "']";
            }
            std::cout << std::endl;
        }

        hpc::h5::group out_group;
        out_group.create(out_file, group_out);

        std::vector<std::string> fields;
        if (!open_files.empty())
        {
            hpc::h5::group g;
            g.open(*open_files[0], group_in);
            hsize_t num_objs = 0;
            H5Gget_num_objs(g.id(), &num_objs);
            for (hsize_t o = 0; o < num_objs; ++o)
            {
                char name[256];
                H5Gget_objname_by_idx(g.id(), o, name, 256);
                fields.push_back(name);
            }
        }

        int info_rank = -1;
        if (!fields.empty())
            info_rank = rank;

        int global_info_rank = -1;
        MPI_Allreduce(&info_rank, &global_info_rank, 1, MPI_INT, MPI_MAX, _comm->mpi_comm());

        if (global_info_rank == -1)
            continue;

        int n_fields = fields.size();
        MPI_Bcast(&n_fields, 1, MPI_INT, global_info_rank, _comm->mpi_comm());

        fields.resize(n_fields);
        for (int k = 0; k < n_fields; ++k)
        {
            int len = fields[k].size();
            MPI_Bcast(&len, 1, MPI_INT, global_info_rank, _comm->mpi_comm());
            fields[k].resize(len);
            MPI_Bcast(const_cast<char*>(fields[k].data()), len, MPI_CHAR, global_info_rank,
                      _comm->mpi_comm());
        }

        for (const auto& fname : fields)
        {
            // Skip array fields (only process 1D scalar fields)
            // Array fields like SfrBulgeSTEPS have shape (n_galaxies, n_timesteps)
            // We only want fields with shape (n_galaxies,)
            bool is_array_field = false;
            if (rank == global_info_rank)
            {
                hpc::h5::group g;
                g.open(*open_files[0], group_in);
                hpc::h5::dataset ds = g.dataset(fname);
                hpc::h5::dataspace sp = ds.dataspace();
                int ndims = sp.simple_extent_num_dims();

                // Only process 1D fields (ndims == 1)
                if (ndims != 1)
                {
                    is_array_field = true;
                    if (_verb >= 1)
                    {
                        std::cout << "    → Skipping array field: " << fname << " (ndims=" << ndims
                                  << ")" << std::endl;
                    }
                }
            }

            int skip_field = is_array_field ? 1 : 0;
            MPI_Bcast(&skip_field, 1, MPI_INT, global_info_rank, _comm->mpi_comm());

            if (skip_field)
                continue;

            hpc::h5::datatype type = hpc::h5::datatype::native_float;

            int type_code = 0;
            if (rank == global_info_rank)
            {
                hpc::h5::group g;
                g.open(*open_files[0], group_in);
                H5T_class_t cls = g.dataset(fname).type_class();
                size_t sz = H5Tget_size(g.dataset(fname).datatype().id());
                if (cls == H5T_INTEGER)
                {
                    if (sz == 8)
                        type_code = 2;
                    else
                        type_code = 1;
                }
            }
            MPI_Bcast(&type_code, 1, MPI_INT, global_info_rank, _comm->mpi_comm());

            if (type_code == 1)
                type = hpc::h5::datatype::native_int;
            else if (type_code == 2)
                type = hpc::h5::datatype::native_llong;
            else
                type = hpc::h5::datatype::native_float;

            hpc::h5::dataspace dspace(total_gals);
            hpc::h5::dataset dset(out_group, fname, type, dspace, hpc::h5::property_list());

            if (my_count > 0 && !open_files.empty())
            {
                unsigned long long current_offset = my_offset;
                for (auto& fptr : open_files)
                {
                    hpc::h5::group g;
                    g.open(*fptr, group_in);
                    if (g.has_link(fname))
                    {
                        hpc::h5::dataset src = g.dataset(fname);
                        hsize_t count = src.dataspace().size();

                        // DEBUG: Verify dataset size matches expectation
                        hpc::h5::dataspace check_space = dset.dataspace();
                        hsize_t actual_dataset_size = check_space.size();

                        if (_verb >= 2 && rank == 0)
                        {
                            std::cout << "    DEBUG: " << group_out << "/" << fname << std::endl;
                            std::cout << "      total_gals=" << total_gals
                                      << ", actual_dataset_size=" << actual_dataset_size
                                      << ", count=" << count
                                      << ", current_offset=" << current_offset
                                      << ", offset+count=" << (current_offset + count) << std::endl;
                        }

                        if (current_offset + count > actual_dataset_size)
                        {
                            std::cerr << "ERROR: Write would exceed dataset bounds!" << std::endl;
                            std::cerr << "  Snapshot: " << group_out << std::endl;
                            std::cerr << "  Field: " << fname << std::endl;
                            std::cerr << "  Dataset size: " << actual_dataset_size << std::endl;
                            std::cerr << "  Requested total_gals: " << total_gals << std::endl;
                            std::cerr << "  Write offset: " << current_offset << std::endl;
                            std::cerr << "  Write count: " << count << std::endl;
                            std::cerr << "  offset + count: " << (current_offset + count)
                                      << std::endl;
                            throw std::runtime_error("Dataset size mismatch detected");
                        }

                        size_t el_size = H5Tget_size(type.id());
                        std::vector<char> buf(count * el_size);
                        src.read(buf.data(), type);

                        std::vector<hsize_t> start = {current_offset};
                        std::vector<hsize_t> count_vec = {count};
                        std::vector<hsize_t> empty_v;
                        hpc::h5::dataspace file_space = dset.dataspace();
                        file_space.select_hyperslab<std::vector<hsize_t>>(H5S_SELECT_SET, count_vec,
                                                                          start, empty_v, empty_v);

                        hpc::h5::dataspace mem(count);
                        dset.write(buf.data(), type, mem, file_space);

                        current_offset += count;
                    }
                }
            }
        }
    }
}

void sage2kdtree_application::operator()()
{
    struct PhaseMetrics
    {
        std::string name;
        double duration;
        long peak_rss;
    };
    std::vector<PhaseMetrics> metrics;

    try
    {
        // Validate inputs before processing
        if (_verb && _comm->rank() == 0)
        {
            std::cout << "=== Validating Inputs ===" << std::endl;
        }
        _validate_inputs();
        if (_verb && _comm->rank() == 0)
        {
            std::cout << "  ✓ All inputs validated successfully\n" << std::endl;
        }

        // Setup intermediate file paths
        _setup_intermediate_paths();

        // Phase 1 (Direct): SAGE HDF5 → Snapshot Organized
        // Replaces Phase 1 & 3
        if (_verb && _comm->rank() == 0)
        {
            std::cout << "=== Phase 1 (Direct): SAGE HDF5 → Snapshot Organized ===" << std::endl;
        }
        auto p1_start = std::chrono::high_resolution_clock::now();
        phase1_direct_sage_to_snapshot();
        _comm->barrier();
        auto p1_end = std::chrono::high_resolution_clock::now();
        metrics.push_back({"Phase 1 (Direct)",
                           std::chrono::duration<double>(p1_end - p1_start).count(),
                           get_peak_rss_kb()});

        if (_verb && _comm->rank() == 0)
        {
            std::cout << "  → Written: " << _bysnap_fn << std::endl;
        }

        // Phase 4: Build KD-Tree Index
        if (_verb && _comm->rank() == 0)
        {
            std::cout << "\n=== Phase 4: Building KD-tree spatial index ===" << std::endl;
        }
        auto p4_start = std::chrono::high_resolution_clock::now();
        phase4_build_kdtree_index();
        _comm->barrier();
        auto p4_end = std::chrono::high_resolution_clock::now();
        metrics.push_back({"Phase 4", std::chrono::duration<double>(p4_end - p4_start).count(),
                           get_peak_rss_kb()});

        // Success summary
        if (_verb && _comm->rank() == 0)
        {
            std::cout << "\n=== Conversion Complete ===" << std::endl;
            std::cout << "Final output: " << _output_fn << std::endl;
            std::cout << "\nIntermediate files (for debugging):" << std::endl;
            std::cout << "  - " << _depthfirst_fn << std::endl;
            std::cout << "  - " << _enhanced_fn << std::endl;
            std::cout << "  - " << _bysnap_fn << std::endl;

            // Performance Report
            std::cout << "\n=== Performance Breakdown ===" << std::endl;
            std::cout << std::left << std::setw(10) << "Phase" << std::right << std::setw(15)
                      << "Time (s)" << std::right << std::setw(15) << "Peak RSS (KB)" << std::endl;
            std::cout << std::string(40, '-') << std::endl;

            double total_time = 0;
            for (const auto& m : metrics)
            {
                std::cout << std::left << std::setw(10) << m.name << std::right << std::setw(15)
                          << std::fixed << std::setprecision(2) << m.duration << std::right
                          << std::setw(15) << m.peak_rss << std::endl;
                total_time += m.duration;
            }
            std::cout << std::string(40, '-') << std::endl;
            std::cout << std::left << std::setw(10) << "Total" << std::right << std::setw(15)
                      << std::fixed << std::setprecision(2) << total_time << std::right
                      << std::setw(15) << get_peak_rss_kb() << std::endl;
        }
    }
    catch (PipelineException& e)
    {
        if (_comm->rank() == 0)
        {
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
    }
    catch (std::exception& e)
    {
        if (_comm->rank() == 0)
        {
            std::cerr << "\n=== UNEXPECTED ERROR ===" << std::endl;
            std::cerr << "Error: " << e.what() << std::endl;
        }
        _comm->abort(1);
    }
}

//==============================================================================
// Utility: Setup intermediate file paths
//==============================================================================

void sage2kdtree_application::_setup_intermediate_paths()
{
    hpc::fs::path output_dir = _output_fn.parent_path();
    std::string base_stem = _output_fn.stem().string();

    // Remove -kdtree suffix if present
    if (base_stem.size() > 7 && base_stem.substr(base_stem.size() - 7) == "-kdtree")
    {
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

void sage2kdtree_application::_validate_inputs()
{
    // Only rank 0 performs validation to avoid redundant checks
    if (_comm->rank() != 0)
        return;

    std::vector<std::string> errors;
    std::vector<std::string> warnings;

    // 1. Check file existence
    if (_verb >= 2)
        std::cout << "  → Checking input file existence..." << std::endl;

    if (!hpc::fs::exists(_param_fn))
    {
        errors.push_back("Parameter file not found: " + _param_fn.string() +
                         "\n      Suggestion: Check the path to your SAGE "
                         "parameter file (.par)");
    }

    if (!hpc::fs::exists(_alist_fn))
    {
        errors.push_back("Expansion factor list not found: " + _alist_fn.string() +
                         "\n      Suggestion: This file should contain scale "
                         "factors (a_list) from your simulation");
    }

    if (!hpc::fs::exists(_sage_dir))
    {
        errors.push_back("SAGE output directory not found: " + _sage_dir.string() +
                         "\n      Suggestion: Check the path to your SAGE HDF5 "
                         "output directory");
    }
    else if (!hpc::fs::is_directory(_sage_dir))
    {
        errors.push_back("SAGE output path is not a directory: " + _sage_dir.string());
    }

    // Early exit if files don't exist
    if (!errors.empty())
    {
        std::string msg = "\n=== INPUT VALIDATION FAILED ===\n";
        for (const auto& err : errors)
        {
            msg += "  ✗ " + err + "\n";
        }
        throw PipelineException(0, msg);
    }

    // 2. Validate SAGE HDF5 directory has .hdf5 or .h5 files
    if (_verb >= 2)
        std::cout << "  → Checking for SAGE HDF5 files..." << std::endl;

    bool found_hdf5_files = false;
    hpc::fs::directory_iterator end_iter;
    for (hpc::fs::directory_iterator dir_itr(_sage_dir); dir_itr != end_iter; ++dir_itr)
    {
        if (hpc::fs::is_regular_file(dir_itr->status()))
        {
            std::string ext = dir_itr->path().extension().string();
            if (ext == ".hdf5" || ext == ".h5")
            {
                found_hdf5_files = true;
                break;
            }
        }
    }

    if (!found_hdf5_files)
    {
        errors.push_back("No HDF5 files found in SAGE directory: " + _sage_dir.string() +
                         "\n      Suggestion: SAGE HDF5 output should have .hdf5 or .h5 "
                         "extension");
    }

    // 3. Validate first SAGE HDF5 file structure
    if (_verb >= 2)
        std::cout << "  → Validating SAGE HDF5 structure..." << std::endl;

    // Find first substantial HDF5 file (skip small master files)
    // SAGE typically creates model_0.hdf5, model_1.hdf5, etc. (data files)
    // and model.hdf5 (small master file with just header)
    hpc::fs::path first_hdf5;
    const size_t min_file_size = 100000; // 100KB - data files are typically much larger

    for (hpc::fs::directory_iterator dir_itr(_sage_dir); dir_itr != end_iter; ++dir_itr)
    {
        if (hpc::fs::is_regular_file(dir_itr->status()))
        {
            std::string ext = dir_itr->path().extension().string();
            std::string stem = dir_itr->path().stem().string();

            // Look for .hdf5 or .h5 files
            if (ext == ".hdf5" || ext == ".h5")
            {
                // Skip small files (likely master files without data)
                size_t file_size = hpc::fs::file_size(dir_itr->path());
                if (file_size < min_file_size)
                {
                    if (_verb >= 2)
                    {
                        std::cout << "    Skipping small file (likely master): "
                                  << dir_itr->path().filename() << " (" << file_size << " bytes)"
                                  << std::endl;
                    }
                    continue;
                }

                // Prefer files matching model_N.hdf5 pattern (actual data files)
                if (stem.find("model_") == 0 && stem.length() > 6)
                {
                    first_hdf5 = dir_itr->path();
                    break;
                }

                // Otherwise accept any substantial file
                if (first_hdf5.empty())
                {
                    first_hdf5 = dir_itr->path();
                }
            }
        }
    }

    if (!first_hdf5.empty())
    {
        if (_verb >= 2)
        {
            std::cout << "    Validating file: " << first_hdf5.filename() << std::endl;
        }

        try
        {
            hpc::h5::file test_file(first_hdf5.string(), H5F_ACC_RDONLY);

            // Check for Header/Simulation group
            bool has_header_simulation = test_file.has_link("Header/Simulation");
            if (!has_header_simulation)
            {
                warnings.push_back("Header/Simulation group not found in " +
                                   first_hdf5.filename().string() +
                                   "\n      Cosmology parameters may not be available");
            }
            else
            {
                // Validate cosmology attributes
                auto sim_group = test_file.group("Header/Simulation");
                std::vector<std::string> required_attrs = {"BoxSize", "Hubble_h", "Omega_m",
                                                           "Omega_lambda"};
                std::vector<std::string> missing_attrs;

                for (const auto& attr : required_attrs)
                {
                    bool found = false;
                    // Check for PascalCase/CamelCase (standard)
                    if (H5Aexists(sim_group.id(), attr.c_str()) > 0)
                    {
                        found = true;
                    }
                    // Check for lowercase (observed in some SAGE outputs)
                    if (!found)
                    {
                        std::string lower_attr = attr;
                        if (attr == "Omega_m")
                            lower_attr = "omega_matter";
                        else if (attr == "Omega_lambda")
                            lower_attr = "omega_lambda";
                        else if (attr == "BoxSize")
                            lower_attr = "box_size";
                        else if (attr == "Hubble_h")
                            lower_attr = "hubble_h";
                        else
                            std::transform(lower_attr.begin(), lower_attr.end(), lower_attr.begin(),
                                           ::tolower);

                        if (H5Aexists(sim_group.id(), lower_attr.c_str()) > 0)
                        {
                            found = true;
                        }
                    }

                    if (!found)
                    {
                        missing_attrs.push_back(attr);
                    }
                }

                if (!missing_attrs.empty())
                {
                    std::string missing_list;
                    for (const auto& attr : missing_attrs)
                    {
                        missing_list += attr + " ";
                    }
                    warnings.push_back(
                        "Missing cosmology attributes in Header/Simulation: " + missing_list +
                        "\n      These will be read from parameter file instead");
                }
            }

            // Check for at least one Snap_N group
            int snap_count = 0;
            while (test_file.has_link("Snap_" + std::to_string(snap_count)))
            {
                snap_count++;
            }

            if (snap_count == 0)
            {
                errors.push_back("No Snap_N groups found in " + first_hdf5.filename().string() +
                                 "\n      Suggestion: SAGE HDF5 output should contain "
                                 "Snap_0, Snap_1, etc.");
            }
            else
            {
                if (_verb >= 2)
                    std::cout << "    Found " << snap_count << " snapshots" << std::endl;

                // Validate required fields in Snap_0
                auto snap0 = test_file.group("Snap_0");
                std::vector<std::string> required_fields = {"Posx", "Posy", "Posz",
                                                            "SAGETreeIndex"};
                std::vector<std::string> missing_fields;

                for (const auto& field : required_fields)
                {
                    if (!snap0.has_link(field))
                    {
                        missing_fields.push_back(field);
                    }
                }

                if (!missing_fields.empty())
                {
                    std::string missing_list;
                    for (const auto& field : missing_fields)
                    {
                        missing_list += field + " ";
                    }
                    errors.push_back("Missing required fields in Snap_0: " + missing_list +
                                     "\n      Suggestion: SAGE HDF5 output must contain position "
                                     "(Posx/y/z) and tree index fields");
                }
            }
        }
        catch (std::exception& e)
        {
            errors.push_back("Failed to open/validate SAGE HDF5 file " +
                             first_hdf5.filename().string() + ": " + std::string(e.what()) +
                             "\n      Suggestion: Ensure SAGE output is in HDF5 "
                             "format (OutputFormat sage_hdf5)");
        }
    }

    // 4. Validate parameter file has required keys
    if (_verb >= 2)
        std::cout << "  → Validating parameter file..." << std::endl;

    std::ifstream param_file(_param_fn.native());
    if (param_file.good())
    {
        std::set<std::string> found_params;
        std::string line;
        while (std::getline(param_file, line))
        {
            if (line.empty() || line[0] == '%')
                continue;
            std::stringstream ss(line);
            std::string key;
            ss >> key;
            found_params.insert(key);
        }

        std::vector<std::string> required_params = {"BoxSize", "Hubble_h"};
        std::vector<std::string> recommended_params = {"Omega_Lambda", "OmegaLambda", "Omega_m",
                                                       "Omega"};

        for (const auto& param : required_params)
        {
            if (found_params.find(param) == found_params.end())
            {
                warnings.push_back("Parameter '" + param + "' not found in " +
                                   _param_fn.filename().string());
            }
        }

        bool has_omega_lambda = (found_params.find("Omega_Lambda") != found_params.end() ||
                                 found_params.find("OmegaLambda") != found_params.end());
        bool has_omega_m = (found_params.find("Omega_m") != found_params.end() ||
                            found_params.find("Omega") != found_params.end());

        if (!has_omega_lambda || !has_omega_m)
        {
            warnings.push_back("Cosmology parameters (Omega_Lambda, Omega_m) "
                               "incomplete in parameter file");
        }
    }

    // 5. Validate expansion factor list
    if (_verb >= 2)
        std::cout << "  → Validating expansion factor list..." << std::endl;

    std::ifstream alist_file(_alist_fn.native());
    if (alist_file.good())
    {
        int line_count = 0;
        std::string line;
        while (std::getline(alist_file, line))
        {
            if (line.empty() || line[0] == '%')
                continue;
            line_count++;
        }

        if (line_count == 0)
        {
            errors.push_back("Expansion factor list is empty: " + _alist_fn.filename().string() +
                             "\n      Suggestion: File should contain one scale factor per "
                             "snapshot");
        }
        else
        {
            if (_verb >= 2)
                std::cout << "    Found " << line_count << " scale factors" << std::endl;
        }
    }

    // Report results
    if (!warnings.empty() && _verb >= 1)
    {
        std::cout << "\n  Warnings:" << std::endl;
        for (const auto& warn : warnings)
        {
            std::cout << "  ⚠  " << warn << std::endl;
        }
        std::cout << std::endl;
    }

    if (!errors.empty())
    {
        std::string msg = "\n=== INPUT VALIDATION FAILED ===\n";
        for (const auto& err : errors)
        {
            msg += "  ✗ " + err + "\n";
        }
        throw PipelineException(0, msg);
    }
}

//==============================================================================
// Phase 1: SAGE HDF5 → Depth-First Ordered
// Based on sageh5toh5.cc lines 59-601
//==============================================================================

void sage2kdtree_application::_load_param(hpc::fs::path const& fn)
{
    if (_verb >= 2 && _comm->rank() == 0)
    {
        std::cout << "  → Loading parameters from " << fn.filename() << std::endl;
    }

    std::ifstream file(fn.native());
    if (!file.good())
    {
        throw PipelineException(1, "Failed to open parameter file: " + fn.string() +
                                       "\n      Check that the file exists and is readable");
    }

    std::string line;
    int line_num = 0;
    while (std::getline(file, line))
    {
        line_num++;
        if (line.empty() || line[0] == '%')
            continue;

        std::stringstream ss(line);
        std::string key;
        ss >> key;

        if (key == "BoxSize")
            ss >> _box_size;
        else if (key == "Hubble_h")
        {
            ss >> _hubble;
            _hubble *= 100.0;
        }
        else if (key == "Omega_Lambda" || key == "OmegaLambda")
            ss >> _omega_l;
        else if (key == "Omega_m" || key == "Omega")
            ss >> _omega_m;
    }

    // Validate that we got the critical parameters
    if (_box_size <= 0.0 || _hubble <= 0.0)
    {
        throw PipelineException(1, "Invalid or missing cosmology parameters in " +
                                       fn.filename().string() +
                                       "\n      Required: BoxSize > 0, Hubble_h > 0" +
                                       "\n      Got: BoxSize=" + std::to_string(_box_size) +
                                       ", Hubble_h=" + std::to_string(_hubble / 100.0));
    }

    if (_verb >= 2 && _comm->rank() == 0)
    {
        std::cout << "    Loaded: BoxSize=" << _box_size << ", Hubble_h=" << _hubble / 100.0
                  << ", Omega_m=" << _omega_m << ", Omega_Lambda=" << _omega_l << std::endl;
    }
}

void sage2kdtree_application::_load_redshifts(hpc::fs::path const& fn)
{
    if (_verb >= 2 && _comm->rank() == 0)
    {
        std::cout << "  → Loading redshifts from " << fn.filename() << std::endl;
    }

    std::ifstream file(fn.native());
    if (!file.good())
    {
        throw PipelineException(1, "Failed to open expansion factor list: " + fn.string() +
                                       "\n      Check that the file exists and is readable");
    }

    std::string line;
    int line_num = 0;
    while (std::getline(file, line))
    {
        line_num++;
        if (line.empty() || line[0] == '%')
            continue;

        std::stringstream ss(line);
        double a;
        ss >> a;

        if (a <= 0.0 || a > 1.0)
        {
            throw PipelineException(1, "Invalid scale factor at line " + std::to_string(line_num) +
                                           " in " + fn.filename().string() +
                                           ": a=" + std::to_string(a) +
                                           "\n      Scale factors must be in range (0, 1]");
        }

        _redshifts.push_back(1.0 / a - 1.0);
    }

    if (_redshifts.empty())
    {
        throw PipelineException(1, "No valid scale factors found in " + fn.string() +
                                       "\n      File should contain one scale factor per line");
    }

    if (_verb >= 2 && _comm->rank() == 0)
    {
        std::cout << "    Loaded " << _redshifts.size() << " redshifts "
                  << "(z_min=" << *std::min_element(_redshifts.begin(), _redshifts.end())
                  << ", z_max=" << *std::max_element(_redshifts.begin(), _redshifts.end()) << ")"
                  << std::endl;
    }
}

void sage2kdtree_application::_copy_header_to_output()
{
    if (_verb >= 1)
    {
        std::cout << "  → Initializing output file and copying metadata..." << std::endl;
    }

    // 1. Create output file (truncate mode) - Sequential access (Rank 0 only)
    // This ensures file exists for subsequent collective open
    try
    {
        hpc::h5::file out_file(_output_fn.string(), H5F_ACC_TRUNC);

        // 2. Find a valid SAGE HDF5 input file that contains the Header group
        hpc::fs::path sage_file_path;
        hpc::fs::directory_iterator end_iter;
        for (hpc::fs::directory_iterator dir_itr(_sage_dir); dir_itr != end_iter; ++dir_itr)
        {
            if (hpc::fs::is_regular_file(dir_itr->status()))
            {
                std::string ext = dir_itr->path().extension().string();
                std::string stem = dir_itr->path().stem().string();
                if (ext == ".hdf5" || ext == ".h5")
                {
                    try
                    {
                        // Check if file has "Header" group
                        hpc::h5::file check_file(dir_itr->path().string(), H5F_ACC_RDONLY);
                        if (check_file.has_link("Header"))
                        {
                            // Prefer model.hdf5 (master) or model_0.hdf5
                            if (stem == "model" || stem == "model_0")
                            {
                                sage_file_path = dir_itr->path();
                                break;
                            }
                            if (sage_file_path.empty())
                            {
                                sage_file_path = dir_itr->path();
                            }
                        }
                    }
                    catch (...)
                    {
                        // Ignore files we can't open
                    }
                }
            }
        }

        if (sage_file_path.empty())
        {
            if (_verb >= 1)
                std::cout << "    Warning: No SAGE HDF5 files found with Header group."
                          << std::endl;
            return;
        }

        // 3. Open input file and copy Header group
        if (_verb >= 2)
            std::cout << "    Copying Header from " << sage_file_path.filename()
                      << " to SageOutputHeader" << std::endl;

        hpc::h5::file in_file(sage_file_path.string(), H5F_ACC_RDONLY);

        // H5Ocopy(src_loc_id, src_name, dst_loc_id, dst_name, ocpypl_id,
        // lcpl_id)
        hid_t status = H5Ocopy(in_file.id(), "Header", out_file.id(), "SageOutputHeader",
                               H5P_DEFAULT, H5P_DEFAULT);
        if (status < 0)
        {
            if (_verb >= 1)
                std::cout << "    Warning: Failed to copy Header group." << std::endl;
        }
        else
        {
            if (_verb >= 2)
                std::cout << "    Successfully copied Header group." << std::endl;
        }
    }
    catch (std::exception& e)
    {
        std::cerr << "Error initializing output file: " << e.what() << std::endl;
        throw;
    }
}

void sage2kdtree_application::phase4_build_kdtree_index()
{
    if (_verb >= 1 && _comm->rank() == 0)
    {
        std::cout << "Phase 4: Building KD-tree spatial index..." << std::endl;
    }

    // Open snapshot-organized file from Phase 3
    hpc::h5::file snap_file(_bysnap_fn.string(), H5F_ACC_RDONLY, *_comm);

    // Create output KD-tree file
    // Rank 0 creates file sequentially and copies Header first
    if (_comm->rank() == 0)
    {
        _copy_header_to_output();
    }
    _comm->barrier();

    // Open collectively with RDWR (file must exist now)
    // Use H5F_ACC_RDWR to preserve the Header we just copied
    hpc::h5::file out_file(_output_fn.string(), H5F_ACC_RDWR, *_comm);

    // Note: Don't create groups explicitly - HDF5 will auto-create them
    // when we create datasets with paths like "lightcone/data" or "data/posx"

    // Get number of snapshots from snapshot_redshifts dataset
    int n_snapshots = 64; // Default value
    if (snap_file.has_link("snapshot_redshifts"))
    {
        hpc::h5::dataset redshift_ds = snap_file.dataset("snapshot_redshifts");
        hpc::h5::dataspace redshift_sp = redshift_ds.dataspace();
        std::vector<hsize_t> redshift_dims(1);
        redshift_sp.simple_extent_dims<std::vector<hsize_t>>(redshift_dims);
        n_snapshots = redshift_dims[0];
    }

    // Get total galaxy count across all snapshots
    unsigned long long total_gals = 0;
    std::vector<unsigned long long> snap_counts(n_snapshots);
    std::vector<unsigned long long> snap_displs(n_snapshots + 1); // n_snapshots + 1 final total

    for (int snap = 0; snap < n_snapshots; ++snap)
    {
        std::ostringstream ss;
        ss << "snapshot" << std::setfill('0') << std::setw(3) << snap;

        if (snap_file.has_link(ss.str()))
        {
            // Phase 3 output: snapshot groups contain columnar field datasets
            // Try to find Posx field (case-sensitive first, then case-insensitive)
            std::string field_path;
            if (snap_file.has_link(ss.str() + "/Posx"))
            {
                field_path = ss.str() + "/Posx";
            }
            else if (snap_file.has_link(ss.str() + "/posx"))
            {
                field_path = ss.str() + "/posx";
            }
            else
            {
                // If Posx not found, open the group and use the first field
                hpc::h5::group snap_group;
                snap_group.open(snap_file, ss.str());
                hsize_t num_objs = 0;
                H5Gget_num_objs(snap_group.id(), &num_objs);
                if (num_objs > 0)
                {
                    char name_buf[256];
                    H5Gget_objname_by_idx(snap_group.id(), 0, name_buf, sizeof(name_buf));
                    field_path = ss.str() + "/" + std::string(name_buf);
                }
            }

            if (!field_path.empty() && snap_file.has_link(field_path))
            {
                hpc::h5::dataset ds = snap_file.dataset(field_path);
                hpc::h5::dataspace sp = ds.dataspace();
                std::vector<hsize_t> dims_vec(1);
                sp.simple_extent_dims<std::vector<hsize_t>>(dims_vec);
                snap_counts[snap] = dims_vec[0];
            }
            else
            {
                snap_counts[snap] = 0;
            }
        }
        else
        {
            snap_counts[snap] = 0;
        }

        snap_displs[snap] = total_gals;
        total_gals += snap_counts[snap];
    }

    // Add final total as the last element
    snap_displs[n_snapshots] = total_gals;

    if (_verb >= 2 && _comm->rank() == 0)
    {
        std::cout << "  Total galaxies across all snapshots: " << total_gals << std::endl;
    }

    // Create lightcone/data compound dataset
    hpc::h5::datatype lc_type(H5Tcreate(H5T_COMPOUND, sizeof(double) * 3));
    hpc::h5::datatype double_type = hpc::h5::datatype::ieee_f64le;
    hpc::h5::datatype ulonglong_type = hpc::h5::datatype::native_ullong;

    lc_type.insert(double_type, "x", 0);
    lc_type.insert(double_type, "y", sizeof(double));
    lc_type.insert(double_type, "z", sizeof(double) * 2);

    // Create lightcone/data dataset - HDF5 will auto-create "lightcone" group
    if (_verb >= 2 && _comm->rank() == 0)
    {
        std::cout << "  Creating lightcone/data dataset..." << std::endl;
    }
    hpc::h5::dataspace lc_space(total_gals);
    hpc::h5::dataset lc_data(out_file, "lightcone/data", lc_type, lc_space,
                             hpc::h5::property_list());

    // Add attributes describing the compound structure
    if (_comm->rank() == 0)
    {
        write_dataset_string_attribute(out_file, "lightcone/data", "Description",
                                       "Compound dataset with galaxy coordinates and metadata");
        write_dataset_string_attribute(out_file, "lightcone/data", "x_units", "Mpc/h");
        write_dataset_string_attribute(out_file, "lightcone/data", "y_units", "Mpc/h");
        write_dataset_string_attribute(out_file, "lightcone/data", "z_units", "Mpc/h");
    }

    if (_verb >= 2 && _comm->rank() == 0)
    {
        std::cout << "  lightcone/data created successfully" << std::endl;
    }

    // Discover fields and create data/* datasets - HDF5 will auto-create "data"
    // group
    if (_comm->rank() == 0 && total_gals > 0)
    {
        if (_verb >= 2)
        {
            std::cout << "  Discovering fields and creating data/* datasets..." << std::endl;
        }
        for (int snap = 0; snap < n_snapshots; ++snap)
        {
            std::ostringstream ss;
            ss << "snapshot" << std::setfill('0') << std::setw(3) << snap;

            if (snap_counts[snap] > 0 && snap_file.has_link(ss.str()))
            {
                if (_verb >= 2)
                {
                    std::cout << "    Opening snapshot group: " << ss.str() << std::endl;
                }
                // Open snapshot group and iterate through field datasets
                hpc::h5::group snap_group;
                snap_group.open(snap_file, ss.str());
                if (_verb >= 2)
                {
                    std::cout << "    Snapshot group opened successfully" << std::endl;
                }
                hsize_t num_objs = 0;
                H5Gget_num_objs(snap_group.id(), &num_objs);

                // Computed fields to exclude from data/* output
                static const std::set<std::string> computed_fields = {
                    "descendant",   "global_descendant", "global_index",
                    "globaltreeid", "local_index",       "localgalaxyid"};

                for (hsize_t ii = 0; ii < num_objs; ++ii)
                {
                    char name_buf[256];
                    H5Gget_objname_by_idx(snap_group.id(), ii, name_buf, sizeof(name_buf));
                    std::string field_name(name_buf);

                    // Skip computed fields - they should not appear in final output
                    std::string field_name_lower = field_name;
                    std::transform(field_name_lower.begin(), field_name_lower.end(),
                                   field_name_lower.begin(), ::tolower);
                    if (computed_fields.count(field_name_lower) > 0)
                    {
                        if (_verb >= 2)
                        {
                            std::cout << "      Skipping computed field: " << field_name
                                      << std::endl;
                        }
                        continue;
                    }

                    std::string field_path = ss.str() + "/" + field_name;

                    // Open the field dataset to get its type
                    hpc::h5::dataset field_ds = snap_file.dataset(field_path);
                    hpc::h5::datatype field_type = field_ds.datatype();

                    // Create data/* dataset - HDF5 auto-creates "data" group (preserve
                    // original case)
                    std::string out_name = field_name;
                    hpc::h5::dataspace field_space(total_gals);
                    hpc::h5::dataset out_ds(out_file, "data/" + out_name, field_type, field_space,
                                            hpc::h5::property_list());

                    // Copy attributes from Phase 3 input
                    std::string desc =
                        read_dataset_string_attribute(snap_file, field_path, "Description");
                    std::string units =
                        read_dataset_string_attribute(snap_file, field_path, "Units");

                    if (!desc.empty())
                        write_dataset_string_attribute(out_file, "data/" + out_name, "Description",
                                                       desc);
                    if (!units.empty())
                        write_dataset_string_attribute(out_file, "data/" + out_name, "Units",
                                                       units);
                }
                break; // Only need to process one snapshot to get field info
            }
        }
    }
    _comm->barrier();

    // Process each snapshot (including empty ones)
    unsigned long long displ = 0;
    for (int snap = 0; snap < n_snapshots; ++snap)
    {
        std::ostringstream ss;
        ss << "snapshot" << std::setfill('0') << std::setw(3) << snap;

        if (_verb >= 2 && _comm->rank() == 0)
        {
            std::cout << "  Processing " << ss.str() << " (" << snap_counts[snap] << " galaxies)..."
                      << std::endl;
        }

        // Process all snapshots, even if empty
        if (snap_file.has_link(ss.str()))
        {
            _process_snapshot(snap_file, ss.str(), out_file, lc_data, displ);
        }
        else
        {
            // Write empty KD-tree structure for missing snapshot
            _write_empty_kdtree(out_file, snap);
        }
    }

    _comm->barrier();

    _comm->barrier();

    // Write snapshot metadata
    if (_comm->rank() == 0)
    {
        hpc::h5::dataspace counts_space(n_snapshots);
        hpc::h5::dataset counts_ds(out_file, "snapshot_counts", ulonglong_type, counts_space,
                                   hpc::h5::property_list());
        counts_ds.write(snap_counts);

        hpc::h5::dataspace displs_space(n_snapshots + 1); // n_snapshots + 1 final total
        hpc::h5::dataset displs_ds(out_file, "snapshot_displs", ulonglong_type, displs_space,
                                   hpc::h5::property_list());
        displs_ds.write(snap_displs);

        // Copy cosmology and snapshot_redshifts
        try
        {
            hpc::h5::copy(snap_file, "cosmology", out_file);
        }
        catch (...)
        {
        }
        try
        {
            hpc::h5::copy(snap_file, "snapshot_redshifts", out_file);
        }
        catch (...)
        {
        }
    }

    _comm->barrier();

    if (_verb >= 1 && _comm->rank() == 0)
    {
        std::cout << "  Phase 4 complete: " << _output_fn << std::endl;
    }
}

void sage2kdtree_application::_process_snapshot(hpc::h5::file& file, std::string const& name,
                                                hpc::h5::file& out_file, hpc::h5::dataset& data,
                                                unsigned long long& displ)
{

    // Read coordinates and metadata from columnar storage
    std::array<std::vector<double>, 3> crds;
    std::vector<unsigned long long> tree_idxs;

    _read_coords(file, name, crds, tree_idxs);

    unsigned long long n_gals = crds[0].size();

    // Extract snapshot number from name (e.g., "snapshot063" -> 63)
    int snap_num = std::stoi(name.substr(8));

    if (n_gals == 0)
    {
        // Write empty KD-tree structure for this snapshot
        _write_empty_kdtree(out_file, snap_num);
        return;
    }

    std::vector<unsigned long long> idxs(n_gals);
    std::iota(idxs.begin(), idxs.end(), 0);

    // Build KD-tree
    data_permuter dp(crds, tree_idxs, idxs);
    hpc::binary_partitioner<std::array<std::vector<double>, 3>::iterator, data_permuter> bp(
        crds.begin(), crds.end(), dp, _ppc);
    hpc::kdtree<> kdt;
    kdt.construct(bp);

    // Write lightcone/data fields separately (matching dstreeinit.cc lines
    // 241-264) Write x, y, z coordinates
    static const char* const coord_names[3] = {"x", "y", "z"};
    for (unsigned ii = 0; ii < 3; ++ii)
    {
        hpc::h5::datatype dt;
        dt.compound(sizeof(double));
        dt.insert(hpc::h5::datatype::native_double, coord_names[ii], 0);
        hpc::h5::property_list props(H5P_DATASET_XFER);
        props.set_preserve();
        data.write(crds[ii].data(), dt, n_gals, displ, hpc::mpi::comm::self, props);
    }

    // Write all other attributes to data/* datasets (from columnar storage)
    _write_attributes(file, name, out_file, idxs, displ);

    // Write KD-tree structure
    _write_kdtree(out_file, name, kdt, bp, crds, displ);

    displ += n_gals;
}

void sage2kdtree_application::_read_coords(hpc::h5::file& file, std::string const& snap_name,
                                           std::array<std::vector<double>, 3>& crds,
                                           std::vector<unsigned long long>& tree_idxs)
{

    // Helper to find field with case-insensitive fallback
    auto find_field = [&](const std::string& base_name,
                          const std::string& alt_name) -> std::string {
        std::string path1 = snap_name + "/" + base_name;
        std::string path2 = snap_name + "/" + alt_name;
        if (file.has_link(path1))
            return path1;
        if (file.has_link(path2))
            return path2;
        return "";
    };

    // Try CamelCase first, then lowercase
    std::string posx_path = find_field("Posx", "posx");
    std::string posy_path = find_field("Posy", "posy");
    std::string posz_path = find_field("Posz", "posz");
    std::string globaltreeid_path = find_field("globaltreeid", "globaltreeid"); // already lowercase

    // Get galaxy count from posx dataset
    if (posx_path.empty())
    {
        crds[0].clear();
        crds[1].clear();
        crds[2].clear();
        tree_idxs.clear();
        return;
    }

    hpc::h5::dataset posx_ds = file.dataset(posx_path);
    hpc::h5::dataspace sp = posx_ds.dataspace();
    std::vector<hsize_t> dims_vec(1);
    sp.simple_extent_dims<std::vector<hsize_t>>(dims_vec);
    unsigned long long n_gals = dims_vec[0];

    if (n_gals == 0)
    {
        crds[0].clear();
        crds[1].clear();
        crds[2].clear();
        tree_idxs.clear();
        return;
    }

    // Resize output vectors
    crds[0].resize(n_gals);
    crds[1].resize(n_gals);
    crds[2].resize(n_gals);
    tree_idxs.resize(n_gals);

    // Read columnar fields (coordinates are float in file, double in memory)
    std::vector<float> posx_float(n_gals);
    std::vector<float> posy_float(n_gals);
    std::vector<float> posz_float(n_gals);

    hpc::h5::dataspace mem_space(n_gals);

    posx_ds.read(posx_float.data(), hpc::h5::datatype::native_float, mem_space, sp);
    file.dataset(posy_path).read(posy_float.data(), hpc::h5::datatype::native_float, mem_space, sp);
    file.dataset(posz_path).read(posz_float.data(), hpc::h5::datatype::native_float, mem_space, sp);

    // Convert float to double
    for (unsigned long long ii = 0; ii < n_gals; ++ii)
    {
        crds[0][ii] = static_cast<double>(posx_float[ii]);
        crds[1][ii] = static_cast<double>(posy_float[ii]);
        crds[2][ii] = static_cast<double>(posz_float[ii]);
    }

    // Read tree indices (stored as long long in Phase 3)
    // If globaltreeid is missing (e.g. Phase 2 skipped or using direct Mode),
    // fill with zeros This prevents crash when dataset does not exist
    if (!globaltreeid_path.empty())
    {
        std::vector<unsigned long long> globaltreeid_data(n_gals);

        file.dataset(globaltreeid_path)
            .read(globaltreeid_data.data(), hpc::h5::datatype::native_llong, mem_space, sp);

        // Copy to output arrays
        for (unsigned long long ii = 0; ii < n_gals; ++ii)
        {
            tree_idxs[ii] = globaltreeid_data[ii];
        }
    }
    else
    {
        std::fill(tree_idxs.begin(), tree_idxs.end(), 0);
    }
}

void sage2kdtree_application::_write_attributes(hpc::h5::file& file, std::string const& snap_name,
                                                hpc::h5::file& out_file,
                                                std::vector<unsigned long long> const& idxs,
                                                unsigned long long displ)
{

    // Read from columnar storage: snapshot group contains separate field datasets
    unsigned long long n_gals = idxs.size();
    if (n_gals == 0)
        return;

    // Open snapshot group and iterate through all field datasets
    hpc::h5::group snap_group;
    snap_group.open(file, snap_name);
    hsize_t num_objs = 0;
    H5Gget_num_objs(snap_group.id(), &num_objs);

    // Computed fields to exclude from data/* output
    static const std::set<std::string> computed_fields = {"descendant",   "global_descendant",
                                                          "global_index", "globaltreeid",
                                                          "local_index",  "localgalaxyid"};

    for (hsize_t ii = 0; ii < num_objs; ++ii)
    {
        char name_buf[256];
        H5Gget_objname_by_idx(snap_group.id(), ii, name_buf, sizeof(name_buf));
        std::string field_name(name_buf);

        // Skip computed fields - they should not appear in final output
        std::string field_name_lower = field_name;
        std::transform(field_name_lower.begin(), field_name_lower.end(), field_name_lower.begin(),
                       ::tolower);
        if (computed_fields.count(field_name_lower) > 0)
        {
            continue;
        }

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
        for (unsigned long long jj = 0; jj < n_gals; ++jj)
        {
            unsigned long long src_idx = idxs[jj];
            memcpy(permuted_data.data() + jj * elem_size, field_data.data() + src_idx * elem_size,
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
        out_ds.write(permuted_data.data(), field_type, out_mem_space, out_file_space,
                     hpc::mpi::comm::self, props);
    }
}

void sage2kdtree_application::_write_kdtree(
    hpc::h5::file& file, std::string const& snap_name, hpc::kdtree<> const& kdt,
    hpc::binary_partitioner<std::array<std::vector<double>, 3>::iterator, data_permuter> const& bp,
    std::array<std::vector<double>, 3> const& crds, unsigned long long displ)
{

    // Build path prefix: "lightcone/snapshot###"
    std::string path = std::string("lightcone/") + snap_name;

    // Write bounds (matching dstreeinit.cc lines 327-334)
    {
        hpc::h5::derive der(2 * sizeof(double));
        der.add(hpc::h5::datatype::native_double, 0, hpc::h5::datatype::ieee_f64be, "minimum");
        der.add(hpc::h5::datatype::native_double, sizeof(double), hpc::h5::datatype::ieee_f64be,
                "maximum");
        hpc::h5::datatype mdt, fdt;
        der.commit(mdt, fdt);

        hpc::h5::dataset dset(file, path + "/bounds", fdt, 3);
        dset.write(kdt.bounds().data(), mdt, 3, 0);
    }

    // Write splits.
    {
        hpc::h5::derive der(sizeof(hpc::kdtree<>::split_type));
        der.add(hpc::h5::datatype::native_double, HOFFSET(hpc::kdtree<>::split_type, pos),
                hpc::h5::datatype::ieee_f64be, "position");
        der.add(hpc::h5::datatype::native_uint, HOFFSET(hpc::kdtree<>::split_type, dim),
                hpc::h5::datatype::std_i32be, "dimension");
        hpc::h5::datatype mdt, fdt;
        der.commit(mdt, fdt);

        hpc::h5::dataset dset(file, path + "/splits", fdt, kdt.splits().size());
        dset.write(kdt.splits().data(), mdt, kdt.splits().size(), 0);
    }

    // Write counts and offsets (matching dstreeinit.cc lines 356-362)
    // Use file.write() which writes with native byte order (little-endian)
    file.write(path + "/cell_counts",
               bp.counts()); // std::vector<unsigned> = 32-bit

    {
        std::vector<unsigned long long> offs(bp.offsets().size());
        for (unsigned ii = 0; ii < offs.size(); ++ii)
            offs[ii] = displ + bp.offsets()[ii];
        file.write(path + "/cell_offs",
                   offs); // std::vector<unsigned long long> = 64-bit
    }
}

void sage2kdtree_application::_write_empty_kdtree(hpc::h5::file& file, int snap_num)
{
    // Write empty KD-tree structure for snapshots with no galaxies
    if (_comm->rank() != 0)
        return; // Only rank 0 writes

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
        der.add(hpc::h5::datatype::native_double, sizeof(double), hpc::h5::datatype::ieee_f64be,
                "maximum");

        hpc::h5::datatype mdt, fdt;
        der.commit(mdt, fdt);

        struct bound_type
        {
            double minimum, maximum;
        };
        bound_type bounds[3] = {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};

        hpc::h5::dataset dset(file, path + "/bounds", fdt, 3);
        dset.write(bounds, mdt, 3, 0);
    }

    // Write empty splits (size 0)
    {
        hpc::h5::derive der(sizeof(double) + sizeof(unsigned));
        der.add(hpc::h5::datatype::native_double, 0, hpc::h5::datatype::ieee_f64be, "position");
        der.add(hpc::h5::datatype::native_uint, sizeof(double), hpc::h5::datatype::std_i32be,
                "dimension");

        hpc::h5::datatype mdt, fdt;
        der.commit(mdt, fdt);

        hpc::h5::dataset dset(file, path + "/splits", fdt, 0);
    }

    // Write cell_counts (1 root cell with 0 galaxies) - use file.write() for
    // native byte order
    {
        std::vector<unsigned> counts(1, 0); // 32-bit unsigned
        file.write(path + "/cell_counts", counts);
    }

    // Write cell_offs (1 offset = 0)
    {
        std::vector<unsigned long long> offs(1, 0); // 64-bit unsigned
        file.write(path + "/cell_offs", offs);
    }
}

//==============================================================================
// Main Entry Point
//==============================================================================

#define HPC_APP_CLASS sage2kdtree_application
#include <libhpc/mpi/main.hh>
