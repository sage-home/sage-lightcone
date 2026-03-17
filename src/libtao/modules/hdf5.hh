#ifndef tao_modules_hdf5_hh
#define tao_modules_hdf5_hh

#include "libtao/base/batch.hh"
#include "libtao/base/filter.hh"
#include "libtao/base/module.hh"
#include <libhpc/libhpc.hh>

constexpr hsize_t DEFAULT_CHUNK_SIZE = 100000;

// Helper function to write string attribute to HDF5 dataset
static void write_dataset_string_attribute(hid_t dset_id, const std::string& attr_name,
                                           const std::string& value)
{
    if (value.empty())
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
}

// Structure to hold metadata for calculated fields
struct FieldMetadata
{
    std::string description;
    std::string units;
};

// Map of calculated field names to their metadata
static const std::map<std::string, FieldMetadata> calculated_fields_metadata = {
    {"redshift_cosmological", {"Cosmological redshift calculated from distance", "dimensionless"}},
    {"redshift_observed", {"Observed redshift including peculiar velocity", "dimensionless"}},
    {"ra", {"Right Ascension", "degrees"}},
    {"dec", {"Declination", "degrees"}},
    {"distance", {"Comoving distance from observer", "Mpc/h"}},
    {"sfr", {"Total Star Formation Rate (disk + bulge)", "Msun/yr"}},
    {"central_spatial_index",
     {"Absolute snapshot index of host central galaxy; -1 for centrals", "dimensionless"}}};

// Helper to copy SageOutputHeader from source to dest file
static void copy_sage_header(const std::string& src_path, hid_t dst_file_id)
{
    try
    {
        hpc::h5::file src_file(src_path, H5F_ACC_RDONLY);
        if (src_file.has_link("SageOutputHeader"))
        {
            // H5Ocopy(src_loc_id, src_name, dst_loc_id, dst_name, ocpypl_id,
            // lcpl_id)
            H5Ocopy(src_file.id(), "SageOutputHeader", dst_file_id, "SageOutputHeader", H5P_DEFAULT,
                    H5P_DEFAULT);
        }
    }
    catch (...)
    {
        // Ignore errors
    }
}

// Helpers for writing scalar attributes
static void write_string_attribute(hid_t loc_id, const std::string& attr_name,
                                   const std::string& value)
{
    if (value.empty())
        return;
    hid_t str_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(str_type, value.length() + 1);
    H5Tset_strpad(str_type, H5T_STR_NULLTERM);
    hid_t space_id = H5Screate(H5S_SCALAR);
    hid_t attr_id =
        H5Acreate2(loc_id, attr_name.c_str(), str_type, space_id, H5P_DEFAULT, H5P_DEFAULT);
    if (attr_id >= 0)
    {
        H5Awrite(attr_id, str_type, value.c_str());
        H5Aclose(attr_id);
    }
    H5Sclose(space_id);
    H5Tclose(str_type);
}

static void write_double_attribute(hid_t loc_id, const std::string& attr_name, double value)
{
    hid_t space_id = H5Screate(H5S_SCALAR);
    hid_t attr_id = H5Acreate2(loc_id, attr_name.c_str(), H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT,
                               H5P_DEFAULT);
    if (attr_id >= 0)
    {
        H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &value);
        H5Aclose(attr_id);
    }
    H5Sclose(space_id);
}

static void write_int_attribute(hid_t loc_id, const std::string& attr_name, int value)
{
    hid_t space_id = H5Screate(H5S_SCALAR);
    hid_t attr_id =
        H5Acreate2(loc_id, attr_name.c_str(), H5T_NATIVE_INT, space_id, H5P_DEFAULT, H5P_DEFAULT);
    if (attr_id >= 0)
    {
        H5Awrite(attr_id, H5T_NATIVE_INT, &value);
        H5Aclose(attr_id);
    }
    H5Sclose(space_id);
}

static void write_bool_attribute(hid_t loc_id, const std::string& attr_name, bool value)
{
    // HDF5 doesn't have a native bool type, use int (0/1) or string
    // ("true"/"false") Using int for standard compliance
    int int_val = value ? 1 : 0;
    write_int_attribute(loc_id, attr_name, int_val);
}

// Function to write LightconeOutputHeader group with attributes from cli_dict
static void write_lightcone_header(hid_t file_id, const tao::cli_dict& dict)
{
    // Create group "LightconeOutputHeader"
    hid_t group_id =
        H5Gcreate2(file_id, "LightconeOutputHeader", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (group_id < 0)
    {
        // If group already exists, open it? Assuming new file or truncated, so
        // create should work unless permissions/corruption.
        // If it exists (e.g. appended), maybe open it.
        group_id = H5Gopen2(file_id, "LightconeOutputHeader", H5P_DEFAULT);
    }

    if (group_id >= 0)
    {
        write_string_attribute(group_id, "dataset", dict._dataset);
        write_double_attribute(group_id, "decmin", dict._decmin);
        write_double_attribute(group_id, "decmax", dict._decmax);
        write_double_attribute(group_id, "ramin", dict._ramin);
        write_double_attribute(group_id, "ramax", dict._ramax);
        write_double_attribute(group_id, "zmin", dict._zmin);
        write_double_attribute(group_id, "zmax", dict._zmax);
        write_string_attribute(group_id, "outfile", dict._outfile);
        write_string_attribute(group_id, "outdir", dict._outdir);
        write_string_attribute(group_id, "filter_field", dict._filter_field);
        write_string_attribute(group_id, "filter_min", dict._filter_min);
        write_string_attribute(group_id, "filter_max", dict._filter_max);
        write_bool_attribute(group_id, "unique", dict._unique);
        write_int_attribute(group_id, "rng_seed", dict._rng_seed);

        // Join output fields into a single string
        std::string outfields_str;
        for (size_t i = 0; i < dict._output_fields.size(); ++i)
        {
            outfields_str += dict._output_fields[i];
            if (i < dict._output_fields.size() - 1)
            {
                outfields_str += ",";
            }
        }
        if (!outfields_str.empty())
        {
            write_string_attribute(group_id, "output_fields", outfields_str);
        }

        H5Gclose(group_id);
    }
}

namespace tao
{
namespace modules
{

template <class Backend>
class hdf5 : public module<Backend>
{
public:
    typedef Backend backend_type;
    typedef module<backend_type> module_type;

    static module_type* factory(const std::string& name) { return new hdf5(name); }

    static void process_cli_options(const cli_dict& global_cli_dict)
    {
        // Access global CLI options
        std::string outfile = global_cli_dict._outfile;
        std::string outdir = global_cli_dict._outdir;

        // Create output directory if it doesn't exist
        if (!boost::filesystem::exists(outdir))
        {
            boost::filesystem::create_directories(outdir);
        }
        // gathering results from all ranks
        if (mpi::comm::world.rank() == 0)
        {
            // If we are the root rank, we need to gather results from all ranks.
            LOGILN("Gathering results from all ", mpi::comm::world.size(), " ranks.");
            std::string file0 = global_cli_dict._outdir + "/" +
                                _construct_filename_static(global_cli_dict._outfile, -1);
            // Create a new file to store the gathered results.
            // This will overwrite the existing file0 if it exists.
            // Create a new file to store the gathered results.
            // This will overwrite the existing file0 if it exists.
            // If the file does not exist, it will be created.
            // If the file exists, it will be truncated.
            hpc::h5::file _file0(file0, H5F_ACC_TRUNC);
            if (!_file0.is_open())
            {
                std::cerr << "Error: Could not open file for gathering results." << std::endl;
                throw silent_terminate();
            }

            // Copy SageOutputHeader from input dataset
            copy_sage_header(global_cli_dict._dataset, _file0.id());

            // Write LightconeOutputHeader with CLI arguments
            write_lightcone_header(_file0.id(), global_cli_dict);

            int rank = 0;
            std::vector<std::string> result_files;
            std::string file1 = global_cli_dict._outdir + "/" +
                                _construct_filename_static(global_cli_dict._outfile, rank);
            while (boost::filesystem::exists(file1))
            {
                result_files.push_back(file1);
                rank++;
                file1 = global_cli_dict._outdir + "/" +
                        _construct_filename_static(global_cli_dict._outfile, rank);
            }
            if (!result_files.empty())
            {
                // Get field list from first rank file by reading dataset names
                std::vector<std::string> field_names;
                {
                    hpc::h5::file first_h5file(result_files[0], H5F_ACC_RDONLY);
                    if (!first_h5file.is_open())
                    {
                        std::cerr << "Error: Could not open file for gathering results."
                                  << std::endl;
                        throw silent_terminate();
                    }

                    // Iterate through all datasets in the file
                    H5G_info_t group_info;
                    H5Gget_info(first_h5file.id(), &group_info);

                    for (hsize_t i = 0; i < group_info.nlinks; i++)
                    {
                        char name_buf[256];
                        H5Lget_name_by_idx(first_h5file.id(), ".", H5_INDEX_NAME, H5_ITER_NATIVE, i,
                                           name_buf, sizeof(name_buf), H5P_DEFAULT);
                        std::string name(name_buf);

                        // Assume all links are datasets (we're iterating at the file root
                        // level)
                        field_names.push_back(name);
                    }
                }

                int file_index = 0;
                for (const auto& field_name : field_names)
                {
                    bool first_file = true;
                    for (const auto& file : result_files)
                    {
                        hpc::h5::file h5file(file, H5F_ACC_RDONLY);
                        if (!h5file.is_open())
                        {
                            std::cerr << "Error: Could not open file for gathering results."
                                      << std::endl;
                            throw silent_terminate();
                        }
                        std::string dataset_name = field_name;
                        // to_lower(dataset_name);
                        //  Open source dataset
                        hpc::h5::dataset src_dset(h5file, dataset_name);

                        // Get datatype and dataspace from source
                        hpc::h5::datatype dtype = src_dset.datatype();

                        // Read data from source
                        hsize_t size = src_dset.extent();
                        std::vector<char> buffer(size * dtype.size());
                        src_dset.read(buffer.data(), dtype);

                        // Create destination dataset with same properties

                        // Create destination dataset with same properties
                        if (first_file)
                        {
                            first_file = false;

                            // Create property list for chunking
                            hpc::h5::property_list props(H5P_DATASET_CREATE);
                            props.set_chunk_size(DEFAULT_CHUNK_SIZE);

                            // Create dataspace with initial size and maxsize
                            hpc::h5::dataspace dspace;
                            dspace.create(1, true);

                            // Create dataset
                            hpc::h5::dataset dst_dset(_file0, dataset_name, dtype, dspace, props);
                            dst_dset.set_extent(size);

                            // Copy HDF5 attributes from source to destination
                            hid_t src_id = H5Dopen2(h5file.id(), dataset_name.c_str(), H5P_DEFAULT);
                            hid_t dst_id = H5Dopen2(_file0.id(), dataset_name.c_str(), H5P_DEFAULT);
                            if (src_id >= 0 && dst_id >= 0)
                            {
                                // Copy Description attribute
                                if (H5Aexists(src_id, "Description") > 0)
                                {
                                    hid_t attr_id = H5Aopen(src_id, "Description", H5P_DEFAULT);
                                    hid_t type_id = H5Aget_type(attr_id);
                                    size_t attr_size = H5Tget_size(type_id);
                                    std::vector<char> attr_buf(attr_size + 1, 0);
                                    H5Aread(attr_id, type_id, attr_buf.data());
                                    H5Aclose(attr_id);

                                    write_dataset_string_attribute(dst_id, "Description",
                                                                   std::string(attr_buf.data()));
                                    H5Tclose(type_id);
                                }

                                // Copy Units attribute
                                if (H5Aexists(src_id, "Units") > 0)
                                {
                                    hid_t attr_id = H5Aopen(src_id, "Units", H5P_DEFAULT);
                                    hid_t type_id = H5Aget_type(attr_id);
                                    size_t attr_size = H5Tget_size(type_id);
                                    std::vector<char> attr_buf(attr_size + 1, 0);
                                    H5Aread(attr_id, type_id, attr_buf.data());
                                    H5Aclose(attr_id);

                                    write_dataset_string_attribute(dst_id, "Units",
                                                                   std::string(attr_buf.data()));
                                    H5Tclose(type_id);
                                }

                                H5Dclose(src_id);
                                H5Dclose(dst_id);
                            }

                            // Write data using chunking
                            write_chunked_dataset(dst_dset, buffer.data(), size, dtype);

                            dst_dset.close();
                            src_dset.close();
                        }
                        else
                        {
                            hpc::h5::dataset dst_dset(_file0, dataset_name);

                            // Get current size and extend
                            hsize_t current_size = dst_dset.extent();
                            dst_dset.set_extent(current_size + size);

                            // Write data using chunking
                            write_chunked_dataset(dst_dset, buffer.data(), size, dtype,
                                                  current_size);

                            dst_dset.close();
                            src_dset.close();
                        }
                        h5file.close();
                        file_index++;
                    }
                }
                for (const auto& file : result_files)
                {
                    boost::filesystem::remove(file);
                }
            }
            _file0.close();
            LOGILN("Output file: ", outfile);
            LOGILN("Output directory: ", outdir);
        }
    }

public:
    hdf5(const std::string& name = std::string())
        : module_type(name)
        , _chunk_size(10000)
        , _ready(false)
        , _records(0)
    {
    }

    virtual ~hdf5() {}

    ///
    ///
    ///
    virtual void initialise(const cli_dict& global_cli_dict,
                            boost::optional<boost::property_tree::ptree> checkpoint =
                                boost::optional<boost::property_tree::ptree>())
    {
        // Don't initialise if we're already doing so.
        if (this->_init)
            return;
        module_type::initialise(global_cli_dict, checkpoint);

        LOGILN(".Initialising HDF5 module.", setindent(2));

        // Get our information.

        std::list<boost::optional<std::string>> field_labels;
        std::string filename = global_cli_dict._outfile;
        _fn = global_cli_dict._outdir + "/" + _construct_filename(filename);

        // Store input kdtree file path for attribute reading
        _input_kdtree_file = global_cli_dict._dataset;

        auto const& input_hdf5_file = global_cli_dict._dataset;

        // Get field list directly from HDF5 file /data/* datasets
        std::vector<std::string> available_fields;
        {
            hpc::h5::file input_file(input_hdf5_file, H5F_ACC_RDONLY);
            hpc::h5::group data_group;
            input_file.open_group("data", data_group);

            // Iterate through datasets in /data group
            H5G_info_t group_info;
            H5Gget_info(data_group.id(), &group_info);

            for (hsize_t i = 0; i < group_info.nlinks; i++)
            {
                char name_buf[256];
                H5Lget_name_by_idx(data_group.id(), ".", H5_INDEX_NAME, H5_ITER_NATIVE, i, name_buf,
                                   sizeof(name_buf), H5P_DEFAULT);
                std::string actual_name = std::string(name_buf);
                available_fields.push_back(actual_name);

                std::string lower_name = actual_name;
                std::transform(lower_name.begin(), lower_name.end(), lower_name.begin(), ::tolower);
                _dataset_name_map[lower_name] = actual_name;
            }
        }

        size_t nfields = available_fields.size();
        LOGILN("Initialising HDF5 module. nfields from input HDF5 file: ", nfields, setindent(2));

        if (global_cli_dict._output_fields.size() == 0)
        {
            LOGILN("No output fields specified via CLI, using all fields from the "
                   "input hdf5 file.");
            for (auto const& field_name : available_fields)
            {
                _fields.push_back(field_name);
                _labels.push_back(field_name); // Use field name as label
                field_labels.push_back(field_name);
            }

            // Also add calculated fields by default
            LOGILN("Adding calculated fields by default.");
            for (const auto& pair : calculated_fields_metadata)
            {
                std::string field_name = pair.first;
                // Check if already added to avoid duplicates (though calculated fields
                // shouldn't be in input)
                bool exists = false;
                for (const auto& existing : _fields)
                {
                    if (existing == field_name)
                    {
                        exists = true;
                        break;
                    }
                }

                // Skip central_spatial_index unless central galaxies mode is active
                if (field_name == "central_spatial_index" && !global_cli_dict._central_galaxies)
                    continue;

                if (!exists)
                {
                    _fields.push_back(field_name);
                    _labels.push_back(field_name);
                    field_labels.push_back(field_name);
                }
            }
        }
        else
        {
            LOGILN("Using output fields specified in the command line.");
            for (auto const& field : global_cli_dict._output_fields)
            {
                _fields.push_back(field);
                _labels.push_back(field);
            }
        }

        // In central galaxies mode, always add central_spatial_index field
        if (global_cli_dict._central_galaxies)
        {
            bool exists = false;
            for (const auto& f : _fields)
            {
                if (f == "central_spatial_index")
                {
                    exists = true;
                    break;
                }
            }
            if (!exists)
            {
                _fields.push_back("central_spatial_index");
                _labels.push_back("central_spatial_index");
            }
        }

        // Note: Alias handling (X, Y, Z, COLOUR) has been removed with XML sidecar
        // removal Note: XML sidecar metadata creation has been removed - using HDF5
        // attributes instead

        // Ensure the full parent directory of the output file exists,
        // in case --outfile contains a path component (e.g. "subdir/out.h5").
        {
            boost::filesystem::path parent = boost::filesystem::path(_fn).parent_path();
            if (!parent.empty() && !boost::filesystem::exists(parent))
                boost::filesystem::create_directories(parent);
        }

        // Open the file. Truncate if we are not reloading.
        _file.open(_fn, H5F_ACC_TRUNC);

        // Copy SageOutputHeader if running in single process mode
        // (In multi-process mode, this is handled during gather in
        // process_cli_options)
        if (mpi::comm::world.size() == 1)
        {
            copy_sage_header(global_cli_dict._dataset, _file.id());
            write_lightcone_header(_file.id(), global_cli_dict);
        }

        _records = 0;
        _ready = false;

        // Get the filter from the lightcone module
        _filt = this->template attribute<tao::filter const*>("filter");

        LOGILN("Done.", setindent(-2));
        // LOGILN(" Exit Now
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        // exit(0);
    }

    static std::string _construct_filename_static(const std::string& base_name, int rank)
    {
        // Construct the filename.
        std::string filename = base_name;
        size_t last_dot = filename.find_last_of('.');
        if (last_dot != std::string::npos)
        {
            std::string ext = filename.substr(last_dot);
            if (ext == ".h5" || ext == ".hdf5")
            {
                filename = filename.substr(0, last_dot);
            }
        }
        if (rank >= 0)
        {
            char chbuff[10];
            snprintf(chbuff, sizeof(chbuff), "%05d", rank);
            std::string rank_str(chbuff);
            filename += "." + rank_str;
        }
        filename += ".h5";
        return filename;
    }

    std::string _construct_filename(const std::string& base_name)
    {
        // Construct the filename.
        std::string filename = base_name;
        size_t last_dot = filename.find_last_of('.');
        if (last_dot != std::string::npos)
        {
            std::string ext = filename.substr(last_dot);
            if (ext == ".h5" || ext == ".hdf5")
            {
                filename = filename.substr(0, last_dot);
            }
        }
        if (mpi::comm::world.size() > 1)
        {
            filename += "." + mpi::rank_string();
        }
        filename += ".h5";
        return filename;
    }

    std::string _encode(std::string _toencode_string)
    {
        std::map<char, std::string> transformations;
        transformations['&'] = std::string("_");
        transformations[' '] = std::string("_");
        transformations['\''] = std::string("_");
        transformations['"'] = std::string("_");
        transformations['>'] = std::string("_");
        transformations['<'] = std::string("_");
        transformations['/'] = std::string("_");
        transformations['('] = std::string("");
        transformations[')'] = std::string("");

        std::string reserved_chars;
        for (auto ti = transformations.begin(); ti != transformations.end(); ti++)
        {
            reserved_chars += ti->first;
        }

        size_t pos = 0;
        while (std::string::npos != (pos = _toencode_string.find_first_of(reserved_chars, pos)))
        {
            _toencode_string.replace(pos, 1, transformations[_toencode_string[pos]]);
            pos++;
        }

        return _toencode_string;
    }

    ///
    ///
    ///
    virtual void execute()
    {
        ASSERT(this->parents().size() == 1);

        // Grab the batch from the parent object.
        const tao::batch<real_type>& bat = this->parents().front()->batch();
        // count the number of entries just once - it will be the same for all the
        // fields
        uint32_t thistime = 0;
        for (auto it = _filt->begin(bat, true); it != _filt->end(bat, true); ++it)
        {
            thistime++;
        }

        // Create datasets for the field names. Note that I can only
        // do this when I have a galaxy object, hence doing it once
        // here.
        bool alt = true;
        h5::dataspace mem_space;
        mem_space.create(1);
        mem_space.select_all();
        if (!_ready)
        {

            auto lblit = _labels.cbegin();

            for (const auto& field : _fields)
            {
                hpc::to_lower((std::string&)field); // input field name of hdf5 must be lowercase
                h5::datatype dtype = _field_type(bat, field);
                std::string field_name = _encode(*lblit); // use label as field name
                if (lblit->empty())
                    field_name = field;

                // Check if this is a 2D array field (VECTOR rank in batch)
                auto bfield = bat.field(field);
                bool is_vector = (std::get<1>(bfield) == tao::batch<real_type>::VECTOR);

                h5::dataset* dset = nullptr;
                hsize_t field_ncols = 0;
                if (is_vector)
                {
                    // Create 2D chunked extendable dataset via raw HDF5 API
                    auto fvtype = std::get<2>(bfield);
                    hsize_t n_cols = 0;
                    hid_t h5_elem_type = H5T_NATIVE_DOUBLE;
                    if (fvtype == tao::batch<real_type>::DOUBLE)
                    {
                        n_cols = bat.template vector<double>(field).n_cols();
                        h5_elem_type = H5T_NATIVE_DOUBLE;
                    }
                    else if (fvtype == tao::batch<real_type>::LONG_LONG)
                    {
                        n_cols = bat.template vector<long long>(field).n_cols();
                        h5_elem_type = H5T_NATIVE_LLONG;
                    }
                    field_ncols = n_cols;
                    hsize_t dims2[2]    = {static_cast<hsize_t>(thistime), n_cols};
                    hsize_t maxdims2[2] = {H5S_UNLIMITED, n_cols};
                    hsize_t chunk2[2]   = {static_cast<hsize_t>(_chunk_size), n_cols};
                    hid_t fspace_id = H5Screate_simple(2, dims2, maxdims2);
                    hid_t plist_id  = H5Pcreate(H5P_DATASET_CREATE);
                    H5Pset_chunk(plist_id, 2, chunk2);
                    hid_t dset_id = H5Dcreate2(_file.id(), field_name.c_str(), h5_elem_type,
                                               fspace_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);
                    H5Pclose(plist_id);
                    H5Sclose(fspace_id);
                    dset = new h5::dataset(dset_id, true);
                }
                else
                {
                    h5::dataspace dspace;
                    dspace.create(1, true);
                    h5::property_list props(H5P_DATASET_CREATE);
                    props.set_chunk_size(_chunk_size);
                    dset = new hpc::h5::dataset(_file, field_name, dtype, dspace, props);
                    dset->set_extent(thistime);
                }
                _dsets.push_back(std::unique_ptr<hpc::h5::dataset>(dset));
                _dset_ncols.push_back(field_ncols);
                _dset_names.push_back(field_name);

                // Copy HDF5 attributes from input kdtree file
                // Skip calculated fields - they don't exist as datasets in the input file
                std::string field_lower_check = field;
                std::transform(field_lower_check.begin(), field_lower_check.end(),
                               field_lower_check.begin(), ::tolower);
                bool is_calculated = calculated_fields_metadata.count(field_lower_check) > 0;

                if (!_input_kdtree_file.empty() && !is_calculated)
                {
                    hid_t input_file_id =
                        H5Fopen(_input_kdtree_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
                    if (input_file_id >= 0)
                    {
                        // Resolve the actual dataset name using our map (handles case
                        // insensitivity)
                        std::string field_lower = field;
                        std::transform(field_lower.begin(), field_lower.end(), field_lower.begin(),
                                       ::tolower);

                        std::string input_dataset_name;
                        if (_dataset_name_map.count(field_lower))
                        {
                            input_dataset_name = _dataset_name_map[field_lower];
                        }
                        else
                        {
                            // Fallback to original name if not found in map (e.g. calculated
                            // fields)
                            input_dataset_name = field;
                        }

                        std::string input_dataset_path = "data/" + input_dataset_name;
                        hid_t input_dset_id =
                            H5Dopen2(input_file_id, input_dataset_path.c_str(), H5P_DEFAULT);
                        hid_t output_dset_id =
                            H5Dopen2(_file.id(), field_name.c_str(), H5P_DEFAULT);

                        if (input_dset_id >= 0 && output_dset_id >= 0)
                        {
                            // Copy Description attribute
                            if (H5Aexists(input_dset_id, "Description") > 0)
                            {
                                hid_t attr_id = H5Aopen(input_dset_id, "Description", H5P_DEFAULT);
                                hid_t type_id = H5Aget_type(attr_id);
                                size_t attr_size = H5Tget_size(type_id);
                                std::vector<char> attr_buf(attr_size + 1, 0);
                                H5Aread(attr_id, type_id, attr_buf.data());
                                H5Aclose(attr_id);
                                write_dataset_string_attribute(output_dset_id, "Description",
                                                               std::string(attr_buf.data()));
                                H5Tclose(type_id);
                            }

                            // Copy Units attribute
                            if (H5Aexists(input_dset_id, "Units") > 0)
                            {
                                hid_t attr_id = H5Aopen(input_dset_id, "Units", H5P_DEFAULT);
                                hid_t type_id = H5Aget_type(attr_id);
                                size_t attr_size = H5Tget_size(type_id);
                                std::vector<char> attr_buf(attr_size + 1, 0);
                                H5Aread(attr_id, type_id, attr_buf.data());
                                H5Aclose(attr_id);
                                write_dataset_string_attribute(output_dset_id, "Units",
                                                               std::string(attr_buf.data()));
                                H5Tclose(type_id);
                            }

                            if (input_dset_id >= 0)
                                H5Dclose(input_dset_id);
                            if (output_dset_id >= 0)
                                H5Dclose(output_dset_id);
                        }

                        H5Fclose(input_file_id);
                    }
                }

                // Add metadata for calculated fields if they exist
                // If attributes were already copied from input file, they will be
                // overwritten if we write them again, but usually calculated fields
                // won't exist in the input file.
                std::string field_lower = field_name;
                std::transform(field_lower.begin(), field_lower.end(), field_lower.begin(),
                               ::tolower);

                if (calculated_fields_metadata.count(field_lower))
                {
                    const auto& meta = calculated_fields_metadata.at(field_lower);
                    hid_t output_dset_id = H5Dopen2(_file.id(), field_name.c_str(), H5P_DEFAULT);
                    if (output_dset_id >= 0)
                    {
                        // Only write if not already present (prefer input file metadata)
                        if (H5Aexists(output_dset_id, "Description") <= 0)
                        {
                            write_dataset_string_attribute(output_dset_id, "Description",
                                                           meta.description);
                        }
                        if (H5Aexists(output_dset_id, "Units") <= 0)
                        {
                            write_dataset_string_attribute(output_dset_id, "Units", meta.units);
                        }
                        H5Dclose(output_dset_id);
                    }
                }

                // Dump first chunk.
                auto val = bat.field(field);

                if (is_vector)
                {
                    // 2D array field: write rows via raw HDF5
                    auto fvtype = std::get<2>(val);
                    _write_2d_batch(dset, bat, field, 0, fvtype);
                }
                else
                {
                    switch (std::get<2>(val))
                    {
                    case tao::batch<real_type>::STRING: {
                        const hpc::view<std::vector<std::string>> thedata =
                            bat.scalar<std::string>(field);
                        unsigned si = 0;
                        for (auto it = _filt->begin(bat, true); it != _filt->end(bat, true);
                             ++it, ++si)
                        {
                            h5::dataspace dspace = dset->dataspace();
                            dspace.select_one(si);
                            _write_field_STRING(bat, *it, field, *dset, dspace, mem_space, thedata);
                        }
                    }
                    break;
                    case tao::batch<real_type>::DOUBLE: {
                        const hpc::view<std::vector<double>> thedata = bat.scalar<double>(field);
                        std::vector<double> buf;
                        buf.reserve(thistime);
                        for (auto it = _filt->begin(bat, true); it != _filt->end(bat, true); ++it)
                            buf.push_back(thedata[*it]);
                        if (!buf.empty())
                            dset->write(buf, (hsize_t)0);
                    }
                    break;
                    case tao::batch<real_type>::FLOAT: {
                        const hpc::view<std::vector<float>> thedata = bat.scalar<float>(field);
                        std::vector<float> buf;
                        buf.reserve(thistime);
                        for (auto it = _filt->begin(bat, true); it != _filt->end(bat, true); ++it)
                            buf.push_back(thedata[*it]);
                        if (!buf.empty())
                            dset->write(buf, (hsize_t)0);
                    }
                    break;
                    case tao::batch<real_type>::INTEGER: {
                        const hpc::view<std::vector<int>> thedata = bat.scalar<int>(field);
                        std::vector<int> buf;
                        buf.reserve(thistime);
                        for (auto it = _filt->begin(bat, true); it != _filt->end(bat, true); ++it)
                            buf.push_back(thedata[*it]);
                        if (!buf.empty())
                            dset->write(buf, (hsize_t)0);
                    }
                    break;
                    case tao::batch<real_type>::UNSIGNED_LONG_LONG: {
                        const hpc::view<std::vector<unsigned long long>> thedata =
                            bat.scalar<unsigned long long>(field);
                        std::vector<long long> buf; // preserve existing signed reinterpretation
                        buf.reserve(thistime);
                        for (auto it = _filt->begin(bat, true); it != _filt->end(bat, true); ++it)
                            buf.push_back(static_cast<long long>(thedata[*it]));
                        if (!buf.empty())
                            dset->write(buf, (hsize_t)0);
                    }
                    break;
                    case tao::batch<real_type>::LONG_LONG: {
                        const hpc::view<std::vector<long long>> thedata =
                            bat.scalar<long long>(field);
                        std::vector<long long> buf;
                        buf.reserve(thistime);
                        for (auto it = _filt->begin(bat, true); it != _filt->end(bat, true); ++it)
                            buf.push_back(thedata[*it]);
                        if (!buf.empty())
                            dset->write(buf, (hsize_t)0);
                    }
                    break;
                    default:
                        ASSERT(0);
                    }
                }
                lblit++;

                // Set number written for first batch.
                _records = thistime;
            }

            // Flag as complete.
            _ready = true;
        }
        else if (alt)
        {
            auto dset_it   = _dsets.begin();
            auto ncols_it  = _dset_ncols.begin();
            auto name_it   = _dset_names.begin();
            uint32_t offset = _records;
            for (const auto& field : _fields)
            {
                h5::dataset* dset      = (*dset_it++).get();
                hsize_t      n_cols    = *ncols_it++;
                const std::string& dset_name = *name_it++;

                if (n_cols > 0)
                {
                    // 2D array field: extend the row dimension then write.
                    hid_t raw_id = H5Dopen2(_file.id(), dset_name.c_str(), H5P_DEFAULT);
                    hsize_t cur_dims[2] = {0, 0};
                    hid_t fspace_id = H5Dget_space(raw_id);
                    H5Sget_simple_extent_dims(fspace_id, cur_dims, nullptr);
                    H5Sclose(fspace_id);
                    hsize_t new_dims[2] = {cur_dims[0] + (hsize_t)thistime, n_cols};
                    H5Dset_extent(raw_id, new_dims);
                    H5Dclose(raw_id);

                    auto val = bat.field(field);
                    _write_2d_batch(dset, bat, field, (hsize_t)offset, std::get<2>(val));
                }
                else
                {
                    {
                        h5::dataspace dspace(*dset);
                        dset->set_extent(dspace.size() + thistime);
                    }
                    h5::dataspace dspace(*dset);

                    // Dump current chunk.
                    auto val = bat.field(field);

                    switch (std::get<2>(val))
                    {
                    case tao::batch<real_type>::STRING: {
                        const hpc::view<std::vector<std::string>> thedata =
                            bat.scalar<std::string>(field);
                        unsigned si = offset;
                        for (auto it = _filt->begin(bat, true); it != _filt->end(bat, true);
                             ++it, ++si)
                        {
                            dspace.select_one(si);
                            _write_field_STRING(bat, *it, field, *dset, dspace, mem_space,
                                                thedata);
                        }
                    }
                    break;
                    case tao::batch<real_type>::DOUBLE: {
                        const hpc::view<std::vector<double>> thedata = bat.scalar<double>(field);
                        std::vector<double> buf;
                        buf.reserve(thistime);
                        for (auto it = _filt->begin(bat, true); it != _filt->end(bat, true); ++it)
                            buf.push_back(thedata[*it]);
                        if (!buf.empty())
                            dset->write(buf, (hsize_t)offset);
                    }
                    break;
                    case tao::batch<real_type>::FLOAT: {
                        const hpc::view<std::vector<float>> thedata = bat.scalar<float>(field);
                        std::vector<float> buf;
                        buf.reserve(thistime);
                        for (auto it = _filt->begin(bat, true); it != _filt->end(bat, true); ++it)
                            buf.push_back(thedata[*it]);
                        if (!buf.empty())
                            dset->write(buf, (hsize_t)offset);
                    }
                    break;
                    case tao::batch<real_type>::INTEGER: {
                        const hpc::view<std::vector<int>> thedata = bat.scalar<int>(field);
                        std::vector<int> buf;
                        buf.reserve(thistime);
                        for (auto it = _filt->begin(bat, true); it != _filt->end(bat, true); ++it)
                            buf.push_back(thedata[*it]);
                        if (!buf.empty())
                            dset->write(buf, (hsize_t)offset);
                    }
                    break;
                    case tao::batch<real_type>::UNSIGNED_LONG_LONG: {
                        const hpc::view<std::vector<unsigned long long>> thedata =
                            bat.scalar<unsigned long long>(field);
                        std::vector<long long> buf; // preserve existing signed reinterpretation
                        buf.reserve(thistime);
                        for (auto it = _filt->begin(bat, true); it != _filt->end(bat, true); ++it)
                            buf.push_back(static_cast<long long>(thedata[*it]));
                        if (!buf.empty())
                            dset->write(buf, (hsize_t)offset);
                    }
                    break;
                    case tao::batch<real_type>::LONG_LONG: {
                        const hpc::view<std::vector<long long>> thedata =
                            bat.scalar<long long>(field);
                        std::vector<long long> buf;
                        buf.reserve(thistime);
                        for (auto it = _filt->begin(bat, true); it != _filt->end(bat, true); ++it)
                            buf.push_back(thedata[*it]);
                        if (!buf.empty())
                            dset->write(buf, (hsize_t)offset);
                    }
                    break;
                    default:
                        ASSERT(0);
                    }
                }
            }
            _records += thistime;
        }
        else
        {
            // Process each element.
            for (auto it = _filt->begin(bat, true); it != _filt->end(bat, true); ++it)
            {
                // Write the fields.
                auto dset_it = _dsets.begin();
                for (const auto& field : _fields)
                {
                    h5::dataset* dset = (*dset_it++).get();
                    hsize_t old_size;
                    {
                        h5::dataspace dspace(*dset);
                        old_size = dspace.size();
                        if (_records == old_size)
                        {
                            dset->set_extent(old_size + 1);
                        }
                    }
                    h5::dataspace dspace(*dset);
                    dspace.select_one(_records);
                    _write_field(bat, *it, field, *dset, dspace);
                }

                // Increment records.
                ++_records;
            }
        }

        // Flush the file here. I want to be sure that all records have
        // been written before potentially checkpointing.
        _file.flush();
    }

    virtual void log_metrics()
    {
        module_type::log_metrics();
        LOGILN(this->_name, " number of records written: ", _records);
    }

    virtual void do_checkpoint(boost::property_tree::ptree& pt)
    {
        std::string fn(_fn);
        std::replace(fn.begin(), fn.end(), '.', '-');
        std::replace(fn.begin(), fn.end(), '/', '-');
        pt.put(std::string("hdf5.") + fn, std::to_string(_records));
    }

protected:
    // Write a 2D VECTOR batch field to an already-sized HDF5 dataset.
    // row_offset: first row index in the file to write to.
    // fvtype: the element type (DOUBLE or LONG_LONG).
    void _write_2d_batch(h5::dataset* dset, const tao::batch<real_type>& bat,
                         const std::string& field, hsize_t row_offset,
                         typename tao::batch<real_type>::field_value_type fvtype)
    {
        if (fvtype == tao::batch<real_type>::DOUBLE)
        {
            const auto& mat = bat.template vector<double>(field);
            hsize_t n_cols = mat.n_cols();
            std::vector<double> buf;
            for (auto it = _filt->begin(bat, true); it != _filt->end(bat, true); ++it)
                for (hsize_t c = 0; c < n_cols; ++c)
                    buf.push_back(mat(*it, c));
            if (buf.empty())
                return;
            hsize_t nrows = buf.size() / n_cols;
            std::vector<hsize_t> file_start = {row_offset, 0};
            std::vector<hsize_t> file_count = {nrows, n_cols};
            h5::dataspace file_space = dset->dataspace();
            file_space.select_hyperslab<std::vector<hsize_t>>(H5S_SELECT_SET, file_count,
                                                              file_start);
            std::vector<hsize_t> mem_dims = {nrows, n_cols};
            h5::dataspace mem_space(mem_dims);
            dset->write(buf.data(), h5::datatype::native_double, mem_space, file_space);
        }
        else if (fvtype == tao::batch<real_type>::LONG_LONG)
        {
            const auto& mat = bat.template vector<long long>(field);
            hsize_t n_cols = mat.n_cols();
            std::vector<long long> buf;
            for (auto it = _filt->begin(bat, true); it != _filt->end(bat, true); ++it)
                for (hsize_t c = 0; c < n_cols; ++c)
                    buf.push_back(mat(*it, c));
            if (buf.empty())
                return;
            hsize_t nrows = buf.size() / n_cols;
            std::vector<hsize_t> file_start = {row_offset, 0};
            std::vector<hsize_t> file_count = {nrows, n_cols};
            h5::dataspace file_space = dset->dataspace();
            file_space.select_hyperslab<std::vector<hsize_t>>(H5S_SELECT_SET, file_count,
                                                              file_start);
            std::vector<hsize_t> mem_dims = {nrows, n_cols};
            h5::dataspace mem_space(mem_dims);
            dset->write(buf.data(), h5::datatype::native_llong, mem_space, file_space);
        }
    }

    // Inside the hdf5 class
    static bool _has_dataset(hpc::h5::file& file, const std::string& dataset_name)
    {
        // Use H5Lexists to check if the link exists
        htri_t link_exists = H5Lexists(file.id(), dataset_name.c_str(), H5P_DEFAULT);
        return (link_exists > 0);
    }

    static void write_chunked_dataset(hpc::h5::dataset& dst_dset, const char* buffer, hsize_t size,
                                      hpc::h5::datatype& dtype, hsize_t current_size = 0)
    {
        const hsize_t chunk_size = DEFAULT_CHUNK_SIZE;
        const hsize_t num_chunks = (size + chunk_size - 1) / chunk_size;

        for (hsize_t chunk = 0; chunk < num_chunks; chunk++)
        {
            // Calculate chunk boundaries
            hsize_t chunk_start = chunk * chunk_size;
            hsize_t chunk_end = std::min(chunk_start + chunk_size, size);
            hsize_t chunk_length = chunk_end - chunk_start;

            // Create memory space for chunk
            hpc::h5::dataspace mem_space;
            mem_space.create(chunk_length, false); // Create with exact size

            // Select hyperslab in file
            hpc::h5::dataspace file_space = dst_dset.dataspace();
            file_space.select_hyperslab(H5S_SELECT_SET,
                                        chunk_length,               // count
                                        current_size + chunk_start, // start
                                        1,                          // stride
                                        1);                         // block

            // Write chunk
            dst_dset.write(buffer + (chunk_start * dtype.size()), dtype, mem_space, file_space);
        }
    }

    h5::datatype _field_type(const tao::batch<real_type>& bat, const std::string& field)
    {
        auto val = bat.field(field);
        switch (std::get<2>(val))
        {
        case tao::batch<real_type>::STRING:
            return h5::datatype::string;

        case tao::batch<real_type>::DOUBLE:
            return h5::datatype::native_double;

        case tao::batch<real_type>::FLOAT:
            return h5::datatype::native_float;

        case tao::batch<real_type>::INTEGER:
            return h5::datatype::native_int;

        case tao::batch<real_type>::UNSIGNED_LONG_LONG:
            return h5::datatype::native_ullong;

        case tao::batch<real_type>::LONG_LONG:
            return h5::datatype::native_llong;

        default:
            ASSERT(0);
        }
    }

    void _write_field(const tao::batch<real_type>& bat, unsigned idx, const std::string& field,
                      h5::dataset& dset, h5::dataspace& dspace)
    {
        auto val = bat.field(field);
        h5::dataspace mem_space;
        mem_space.create(1);
        mem_space.select_all();
        switch (std::get<2>(val))
        {
        case tao::batch<real_type>::STRING: {
            std::string data = bat.scalar<std::string>(field)[idx];
            dset.write(data.c_str(), h5::datatype::string, mem_space, dspace);
            break;
        }

        case tao::batch<real_type>::DOUBLE: {
            double data = bat.scalar<double>(field)[idx];
            dset.write(&data, h5::datatype::native_double, mem_space, dspace);
            break;
        }

        case tao::batch<real_type>::INTEGER: {
            int data = bat.scalar<int>(field)[idx];
            dset.write(&data, h5::datatype::native_int, mem_space, dspace);
            break;
        }

        case tao::batch<real_type>::UNSIGNED_LONG_LONG: {
            unsigned long long data = bat.scalar<unsigned long long>(field)[idx];
            dset.write(&data, h5::datatype::native_ullong, mem_space, dspace);
            break;
        }

        case tao::batch<real_type>::LONG_LONG: {
            long long data = bat.scalar<long long>(field)[idx];
            dset.write(&data, h5::datatype::native_llong, mem_space, dspace);
            break;
        }

        default:
            ASSERT(0);
        }
    }

    void _write_field_STRING(const tao::batch<real_type>& bat, unsigned idx,
                             const std::string& field, h5::dataset& dset, h5::dataspace& dspace,
                             h5::dataspace& mem_space,
                             const hpc::view<std::vector<std::string>>& thedata)
    {
        {
            std::string data = thedata[idx];
            dset.write(data.c_str(), h5::datatype::string, mem_space, dspace);
        }
    }

    void _write_field_DOUBLE(const tao::batch<real_type>& bat, unsigned idx,
                             const std::string& field, h5::dataset& dset, h5::dataspace& dspace,
                             h5::dataspace& mem_space,
                             const hpc::view<std::vector<double>>& thedata)
    {
        {
            double data = thedata[idx];
            dset.write(&data, h5::datatype::native_double, mem_space, dspace);
        }
    }

    void _write_field_FLOAT(const tao::batch<real_type>& bat, unsigned idx,
                            const std::string& field, h5::dataset& dset, h5::dataspace& dspace,
                            h5::dataspace& mem_space, const hpc::view<std::vector<float>>& thedata)
    {
        {
            float data = thedata[idx];
            dset.write(&data, h5::datatype::native_float, mem_space, dspace);
        }
    }

    void _write_field_INTEGER(const tao::batch<real_type>& bat, unsigned idx,
                              const std::string& field, h5::dataset& dset, h5::dataspace& dspace,
                              h5::dataspace& mem_space, const hpc::view<std::vector<int>>& thedata)
    {
        {
            int data = thedata[idx];
            dset.write(&data, h5::datatype::native_int, mem_space, dspace);
        }
    }

    void _write_field_UNSIGNED_LONG_LONG(const tao::batch<real_type>& bat, unsigned idx,
                                         const std::string& field, h5::dataset& dset,
                                         h5::dataspace& dspace, h5::dataspace& mem_space,
                                         const hpc::view<std::vector<unsigned long long>>& thedata)
    {
        {
            long long data = thedata[idx];
            dset.write(&data, h5::datatype::native_llong, mem_space, dspace);
        }
    }

    void _write_field_LONG_LONG(const tao::batch<real_type>& bat, unsigned idx,
                                const std::string& field, h5::dataset& dset, h5::dataspace& dspace,
                                h5::dataspace& mem_space,
                                const hpc::view<std::vector<long long>>& thedata)
    {
        {
            signed long long data = thedata[idx];
            dset.write(&data, h5::datatype::native_llong, mem_space, dspace);
        }
    }

    typename tao::batch<tao::real_type>::field_value_type _datatypeAsTAO(hpc::h5::datatype datatype)
    {
        H5T_class_t dt1 = H5Tget_class(datatype.id());
        int dt1size = H5Tget_size(datatype.id());
        switch (dt1)
        {
        case H5T_class_t::H5T_INTEGER:
            if (dt1size == 4)
            {
                return tao::batch<tao::real_type>::INTEGER;
            }
            else if (dt1size == 8)
            {
                return tao::batch<tao::real_type>::LONG_LONG;
            }
            break;
        case H5T_class_t::H5T_FLOAT:
            if (dt1size == 4)
            {
                return tao::batch<tao::real_type>::DOUBLE;
            }
            break;
        default:
            ASSERT(0);
            break;
        }
        return tao::batch<tao::real_type>::DOUBLE;
    }

protected:
    h5::file _file;
    std::string _fn;
    std::string _input_kdtree_file; // Path to input kdtree file for attribute reading
    std::map<std::string, std::string>
        _dataset_name_map; // Map lowercase field name to actual HDF5 dataset name
    std::list<std::string> _fields;
    unsigned long long _records;
    std::list<std::string> _labels;
    std::list<std::unique_ptr<h5::dataset>> _dsets;
    std::list<hsize_t> _dset_ncols;    // 0 for scalar; n_cols for 2D array fields
    std::list<std::string> _dset_names; // encoded dataset name for each field
    hsize_t _chunk_size;
    bool _ready;
    tao::filter const* _filt;
};

} // namespace modules
} // namespace tao

#endif
