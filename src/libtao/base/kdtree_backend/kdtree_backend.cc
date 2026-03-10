#include "kdtree_backend.hh"
#include "../mandatory_fields.hh"
#include <sstream>

namespace tao
{
namespace backends
{

// Helper function to read string attribute from HDF5 dataset
static std::string read_dataset_attribute(hid_t loc_id, const std::string& dataset_path,
                                          const std::string& attr_name)
{
    hid_t dset_id = H5Dopen2(loc_id, dataset_path.c_str(), H5P_DEFAULT);
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

kdtree_backend::kdtree_backend(hpc::mpi::comm const& comm)
    : _kdt(hpc::mpi::comm::self)
    , _snap(std::numeric_limits<unsigned>::max())
    , _comm(&comm)
{
}

kdtree_backend::kdtree_backend(hpc::fs::path const& fn, hpc::mpi::comm const& comm)
    : _kdt(hpc::mpi::comm::self)
    , _snap(std::numeric_limits<unsigned>::max())
    , _comm(&comm)
{
    open(fn);
}

// For some backends there are different fields expected.
// Note: subsize is read directly from lightcone/data compound dataset,
// so subtree_count is not needed as a separate field in data/
void kdtree_backend::add_conditional_fields(query<real_type>& qry)
{
    if (_central_galaxies_mode)
    {
        qry.add_output_field("central_spatial_index");
    }
}

hpc::mpi::comm const& kdtree_backend::comm() const { return *_comm; }

herr_t static getFields(hid_t g_id, const char* name, const H5L_info_t* info, void* op_data)
{
    std::vector<std::string>* fields = (std::vector<std::string>*)op_data;
    std::string name_lowercase = name;
    std::transform(name_lowercase.begin(), name_lowercase.end(), name_lowercase.begin(), ::tolower);
    fields->push_back(name);
    return 0;
}

void kdtree_backend::connect(const cli_dict& global_cli_dict)
{
    std::cerr << "DEBUG: kdtree_backend::connect() called" << std::endl;
    std::string fn = global_cli_dict._dataset;
    _field_types.clear();
    _central_galaxies_mode = global_cli_dict._central_galaxies;
    _include_orphan_satellites_mode = global_cli_dict._central_galaxies;

    // Open HDF5 file and read field information
    this->open(fn);
    std::cerr << "DEBUG: After open(), _field_types has " << _field_types.size() << " entries"
              << std::endl;
    hpc::h5::group data = kdtree_file().group("data");
    std::vector<std::string> Hfields;
    H5Lvisit(data.id(), H5_INDEX_CRT_ORDER, H5_ITER_NATIVE, getFields, &Hfields);

    //  for (std::string field : Hfields)
    //  {
    //      auto ds = data.dataset(field);

    //      std::cout << field << "
    //      dt="<<ds.datatype()<<","<<ds.type_class()<<std::endl;
    //  }

    // Read field metadata from HDF5 attributes (override XML if present)
    hid_t file_id = kdtree_file().id();
    for (const std::string& field : Hfields)
    {
        std::string field_path = "data/" + field;

        // Read Description attribute
        std::string desc = read_dataset_attribute(file_id, field_path, "Description");
        if (!desc.empty())
        {
            std::string field_lower = field;
            to_lower(field_lower);
            _field_description[field_lower] = desc;
        }

        // Read Units attribute
        std::string units = read_dataset_attribute(file_id, field_path, "Units");
        if (!units.empty())
        {
            std::string field_lower = field;
            to_lower(field_lower);
            _field_units[field_lower] = units;
        }
    }

    // Field name aliases - kdtree files use SAGE CamelCase for SAGE fields
    // Computed fields remain lowercase_with_underscores

    // Aliases for backward compatibility with different naming conventions
    this->_field_map["diskscaleradius"] = "DiskRadius";
    this->_field_map["disk_scale_radius"] = "DiskRadius";
    this->_field_map["stellar_mass"] = "StellarMass";
    this->_field_map["snapnum"] = "SnapNum";

    // Add calculated types.
    this->_field_types.emplace("redshift_cosmological", batch<real_type>::DOUBLE);
    this->_field_types.emplace("redshift_observed", batch<real_type>::DOUBLE);
    this->_field_types.emplace("ra", batch<real_type>::DOUBLE);
    this->_field_types.emplace("dec", batch<real_type>::DOUBLE);
    this->_field_types.emplace("distance", batch<real_type>::DOUBLE);
    this->_field_types.emplace("sfr", batch<real_type>::DOUBLE);

    if (_central_galaxies_mode)
    {
        this->_field_types.emplace("central_spatial_index", batch<real_type>::LONG_LONG);
        if (!_has_central_index)
        {
            std::cerr << "Warning: --centralgalaxies specified but input file has no "
                         "central galaxy index. Run sage2kdtree --centralgalaxies to build it."
                      << std::endl;
        }
    }

#ifndef NDEBUG
    // In debug mode we add some custom types (using SAGE CamelCase).
    this->_field_map["original_x"] = "Posx";
    this->_field_map["original_y"] = "Posy";
    this->_field_map["original_z"] = "Posz";
    this->_field_types.emplace("original_x", batch<real_type>::DOUBLE);
    this->_field_types.emplace("original_y", batch<real_type>::DOUBLE);
    this->_field_types.emplace("original_z", batch<real_type>::DOUBLE);
#endif

    // Make sure we have all the essential fields available. Do this by
    // checking that all the mapped fields exist in the field types.
    // if (this->_init_tbls) {
    //    for (const auto &item : this->_field_map)
    //        EXCEPT(hpc::has(this->_field_types, item.second), "Database is
    //        missing essential field: ",
    //               item.second);
    //}

    std::cerr << "DEBUG: At end of connect(), _field_types has " << _field_types.size()
              << " entries" << std::endl;
    std::cerr << "DEBUG: _field_map has " << _field_map.size() << " entries" << std::endl;
}

void kdtree_backend::close() {}

void kdtree_backend::open(hpc::fs::path const& fn)
{
    LOGBLOCKD("Opening kdtree backend: ", fn.native());
    close();

    // Suppress HDF5 diagnostics while checking accessibility, so that
    // missing/invalid files produce a clean error message instead of an
    // assertion crash.
    H5E_auto2_t old_func;
    void* old_client_data;
    H5Eget_auto2(H5E_DEFAULT, &old_func, &old_client_data);
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
    htri_t accessible = H5Fis_accessible(fn.native().c_str(), H5P_DEFAULT);
    H5Eset_auto2(H5E_DEFAULT, old_func, old_client_data);

    if (accessible <= 0)
    {
        std::cerr << "Error: Dataset file '" << fn.native()
                  << "' does not exist or is not a valid HDF5 file." << std::endl;
        exit(1);
    }

    _file.open(fn.native(), H5F_ACC_RDONLY);
    _snap = std::numeric_limits<unsigned>::max();

    // Check for central galaxy index attribute
    _has_central_index = (H5Aexists(_file.id(), "has_central_galaxy_index") > 0);

    // Load snapshot_displs for central galaxies satellite lookup
    if (_file.has_link("snapshot_displs"))
    {
        _snap_displs = _file.read<std::vector<unsigned long long>>("snapshot_displs");
    }

    _lc_mem_type.compound(sizeof(lightcone_data));
    _lc_mem_type.insert(hpc::h5::datatype::native_double, "x", 0);
    _lc_mem_type.insert(hpc::h5::datatype::native_double, "y", sizeof(double));
    _lc_mem_type.insert(hpc::h5::datatype::native_double, "z", 2 * sizeof(double));

    // Check for required groups before trying to open them, so we get a clean
    // error if the file is valid HDF5 but not a kdtree-indexed dataset.
    {
        H5E_auto2_t old_func2;
        void* old_client_data2;
        H5Eget_auto2(H5E_DEFAULT, &old_func2, &old_client_data2);
        H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
        htri_t data_exists = H5Lexists(_file.id(), "data", H5P_DEFAULT);
        htri_t lc_exists = H5Lexists(_file.id(), "lightcone", H5P_DEFAULT);
        H5Eset_auto2(H5E_DEFAULT, old_func2, old_client_data2);

        if (data_exists <= 0 || lc_exists <= 0)
        {
            std::cerr << "Error: '" << fn.native()
                      << "' is a valid HDF5 file but is not a KD-tree indexed dataset (missing";
            if (data_exists <= 0)
                std::cerr << " '/data'";
            if (lc_exists <= 0)
                std::cerr << " '/lightcone'";
            std::cerr << " group).\n"
                      << "Please pass a file generated by sage2kdtree." << std::endl;
            exit(1);
        }
    }

    // Populate field types from HDF5 /data group datasets
    hpc::h5::group data_group;
    _file.open_group("data", data_group);

    H5G_info_t group_info;
    H5Gget_info(data_group.id(), &group_info);

    std::vector<std::string> available_fields;
    for (hsize_t i = 0; i < group_info.nlinks; i++)
    {
        char name_buf[256];
        H5Lget_name_by_idx(data_group.id(), ".", H5_INDEX_NAME, H5_ITER_NATIVE, i, name_buf,
                           sizeof(name_buf), H5P_DEFAULT);
        std::string field_name(name_buf);
        available_fields.push_back(field_name);

        // Convert to lowercase for backward compatibility alias
        std::string field_name_lower = field_name;
        std::transform(field_name_lower.begin(), field_name_lower.end(), field_name_lower.begin(),
                       ::tolower);

        // Open dataset and get its type class
        hpc::h5::dataset ds(data_group, field_name);
        H5T_class_t type_class = ds.type_class();

        // Map HDF5 type class to batch field_value_type
        batch<real_type>::field_value_type ftype;
        if (type_class == H5T_FLOAT)
        {
            ftype = batch<real_type>::DOUBLE;
        }
        else if (type_class == H5T_INTEGER)
        {
            // For integers, just use LONG_LONG by default
            // (could check size/signedness but this is simpler)
            ftype = batch<real_type>::LONG_LONG;
        }
        else
        {
            // Default to DOUBLE for unknown types
            ftype = batch<real_type>::DOUBLE;
        }

        // Store field type using actual field name (CamelCase from HDF5)
        this->_field_types.emplace(field_name, ftype);

        // Create lowercase alias for backward compatibility (if different from
        // actual name)
        if (field_name_lower != field_name)
        {
            this->_field_map[field_name_lower] = field_name;
        }
    }

    // Validate that KD-tree file has required spatial fields
    auto required_fields = tao::get_kdtree_required_fields();
    std::vector<std::string> missing_fields;

    for (const auto& required : required_fields)
    {
        bool found = false;
        for (const auto& available : available_fields)
        {
            if (tao::iequals(required, available))
            {
                found = true;
                break;
            }
        }
        if (!found)
        {
            missing_fields.push_back(required);
        }
    }

    if (!missing_fields.empty())
    {
        std::ostringstream error_msg;
        error_msg << "KD-tree file is missing required fields:\n";
        for (const auto& field : missing_fields)
        {
            error_msg << "  - " << field << "\n";
        }
        error_msg << "This file may have been created with an older version of "
                     "sage2kdtree.\n";
        error_msg << "Please regenerate the KD-tree file from SAGE HDF5 input.";
        throw std::runtime_error(error_msg.str());
    }
}

void chround(char* a, int ndigits)
{
    int i;
    int len = 0;
    for (i = 0; i < 1000; i++)
    {
        if (a[i] == 0)
            break;
        len++;
    }
    int count_ndigits = 0;
    for (i = 0; i < len; i++)
    {
        if (a[i] == ' ')
            continue;
        if (a[i] == '.')
            continue;
        if (a[i] == '0' && count_ndigits == 0)
            continue;
        count_ndigits++;
        if (count_ndigits == ndigits)
            break;
    }
    char last = a[i];
    bool up = false;
    bool done = false;
    for (int j = i + 1; j < len; j++)
    {
        switch (a[j])
        {
        case '0':
            done = true;
            break;
        case '1':
            done = true;
            break;
        case '2':
            done = true;
            break;
        case '3':
            done = true;
            break;
        case '4':
            done = true;
            break;
        case '5':
            up = true;
            done = true;
            break;
        case '6':
            up = true;
            done = true;
            break;
        case '7':
            up = true;
            done = true;
            break;
        case '8':
            up = true;
            done = true;
            break;
        case '9':
            up = true;
            done = true;
            break;
        }
        if (done)
            break;
    }
    int endofstring_index = i + 1;
    if (up)
    {
        done = false;
        while (!done)
        {
            switch (a[i])
            {
            case '0':
                a[i] = '1';
                done = true;
                break;
            case '1':
                a[i] = '2';
                done = true;
                break;
            case '2':
                a[i] = '3';
                done = true;
                break;
            case '3':
                a[i] = '4';
                done = true;
                break;
            case '4':
                a[i] = '5';
                done = true;
                break;
            case '5':
                a[i] = '6';
                done = true;
                break;
            case '6':
                a[i] = '7';
                done = true;
                break;
            case '7':
                a[i] = '8';
                done = true;
                break;
            case '8':
                a[i] = '9';
                done = true;
                break;
            case '9':
                a[i] = '0';
                break;
            }
            i--;
            if (i < 0 && !done)
            {
                /*
                 * will need to shift all the characters one to the right and put a '1'
                 * at the start
                 */

                for (int j = endofstring_index; j > 0; j--)
                {
                    a[j] = a[j - 1];
                }
                a[0] = '1';
                endofstring_index++;
                done = true;
            }
        }
    }
    a[endofstring_index] = 0;
}

tao::simulation const* kdtree_backend::load_simulation()
{
    LOGBLOCKD("Loading simulation.");

    std::vector<real_type> redshifts = _file.read<std::vector<real_type>>("snapshot_redshifts");
    std::map<int, real_type> zs;
    std::cout.precision(17);
    for (size_t ii = 0; ii < redshifts.size(); ++ii)
    {
        zs[ii] = redshifts[ii];
        char chredshift[1000];
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(12) << zs[ii];
        std::string temp = oss.str();
        strncpy(chredshift, temp.c_str(), 999);
        chredshift[999] = '\0';
        chround(chredshift, 6);
        // std::cout << "replace " << zs[ii] << " with " << chredshift << std::endl;
        zs[ii] = atof(chredshift);
    }

    auto* sim = new tao::simulation(
        _file.read<real_type>("cosmology/box_size"), _file.read<real_type>("cosmology/hubble"),
        _file.read<real_type>("cosmology/omega_m"), _file.read<real_type>("cosmology/omega_l"), zs);
    set_simulation(sim);
    return sim;
}

void kdtree_backend::load_snapshot(unsigned snap)
{
    if (snap != _snap)
    {

        std::string name = make_snapshot_name(snap);
        // LOGBLOCKD("Loading kdtree snapshot: ", name);
        ASSERT(_file.is_open(), "Have not opened an HDF5 kdtree file.");

        hpc::h5::group grp = _file.group(name);
        grp >> _kdt;
        grp.dataset("cell_counts") >> _cell_cnts;
        grp.dataset("cell_offs") >> _cell_offs;
        // LOGILN("Loading kdtree snapshot: ", name);

        // Debug: Print kdtree structure information
        // LOGILN("Kdtree loaded - Dimensions: ", _kdt.n_dims());
        // LOGILN("Kdtree loaded - Branches: ", _kdt.n_branches());
        // LOGILN("Kdtree loaded - Leafs: ", _kdt.n_leafs());
        // LOGILN("Kdtree loaded - Total cells: ", _kdt.n_cells());
        // LOGILN("Kdtree loaded - Bounds size: ", _kdt.bounds().size());
        // LOGILN("Kdtree loaded - Splits size: ", _kdt.splits().size());

        // Load central galaxies satellite index if available
        if (_has_central_index)
        {
            // The centralgalaxies group uses "snapshot###" (not "lightcone/snapshot###")
            std::ostringstream cg_ss;
            cg_ss << "centralgalaxies/snapshot" << hpc::index_string(snap, 3);
            std::string cg_path = cg_ss.str();

            _sat_offs.clear();
            _sat_list.clear();

            std::string offs_path = cg_path + "/satellite_offsets";
            std::string list_path = cg_path + "/satellite_list";

            if (_file.has_link(offs_path))
            {
                hpc::h5::dataset offs_ds = _file.dataset(offs_path);
                hpc::h5::dataspace offs_sp = offs_ds.dataspace();
                std::vector<hsize_t> dims(1);
                offs_sp.simple_extent_dims<std::vector<hsize_t>>(dims);
                _sat_offs.resize(dims[0]);
                offs_ds.read(_sat_offs.data(), hpc::h5::datatype::native_ullong, dims[0], 0);
            }

            if (_file.has_link(list_path))
            {
                hpc::h5::dataset list_ds = _file.dataset(list_path);
                hpc::h5::dataspace list_sp = list_ds.dataspace();
                std::vector<hsize_t> dims(1);
                list_sp.simple_extent_dims<std::vector<hsize_t>>(dims);
                if (dims[0] > 0)
                {
                    _sat_list.resize(dims[0]);
                    list_ds.read(_sat_list.data(), hpc::h5::datatype::native_ullong, dims[0], 0);
                }
            }

            // Build reverse map (satellite abs_idx -> central abs_idx)
            if (!_sat_offs.empty())
            {
                _sat_to_central.clear();
                unsigned long long snap_displ =
                    (snap < _snap_displs.size()) ? _snap_displs[snap] : 0ULL;
                for (unsigned long long ci = 0; ci + 1 < _sat_offs.size(); ++ci)
                {
                    unsigned long long sat_start = _sat_offs[ci];
                    unsigned long long sat_end = _sat_offs[ci + 1];
                    for (unsigned long long k = sat_start; k < sat_end; ++k)
                    {
                        _sat_to_central[snap_displ + _sat_list[k]] = snap_displ + ci;
                    }
                }
            }
        }

        _snap = snap;
    }
}

real_type kdtree_backend::get_max_smoothing() { ASSERT(0); }

std::string kdtree_backend::make_snapshot_name(unsigned snap) const
{
    return std::string("lightcone/snapshot") + hpc::index_string(snap, 3);
}

kdtree_backend::box_galaxy_iterator kdtree_backend::galaxy_begin(query<real_type>& qry,
                                                                 box<real_type> const& box,
                                                                 batch<real_type>* bat,
                                                                 filter const* flt)
{
    load_snapshot(box.snapshot());
    return box_galaxy_iterator(*this, box, bat, box.snapshot());
}

kdtree_backend::box_galaxy_iterator kdtree_backend::galaxy_end(query<real_type>& qry,
                                                               box<real_type> const& box) const
{
    return box_galaxy_iterator(*this);
}

kdtree_backend::lightcone_galaxy_iterator kdtree_backend::galaxy_begin(query<real_type>& qry,
                                                                       lightcone const& lc,
                                                                       batch<real_type>* bat,
                                                                       filter const* flt)
{
    kdtree_backend::lightcone_galaxy_iterator result =
        lightcone_galaxy_iterator(lc, *this, qry, bat, flt);
    return result;
}

kdtree_backend::lightcone_galaxy_iterator kdtree_backend::galaxy_end(query<real_type> const& qry,
                                                                     lightcone const& lc) const
{
    return lightcone_galaxy_iterator();
}

kdtree_backend::tile_galaxy_iterator kdtree_backend::galaxy_begin(query<real_type>& qry,
                                                                  tile<real_type> const& tile,
                                                                  batch<real_type>* bat,
                                                                  filter const* flt, bool first)
{
    return tile_galaxy_iterator(*this, tile, tile.lightcone(), bat, _snap);
}

kdtree_backend::tile_galaxy_iterator kdtree_backend::galaxy_end(query<real_type> const& qry,
                                                                tile<real_type> const& tile) const
{
    return tile_galaxy_iterator(*this);
}

void kdtree_backend::load_lightcone_data(unsigned cell, std::vector<lightcone_data>& data) const
{
    unsigned n_elems = _cell_cnts[cell];
    unsigned long long offs = _cell_offs[cell];
    data.resize(n_elems);
    _file.dataset("lightcone/data").read(data.data(), _lc_mem_type, n_elems, offs);
}

hpc::h5::file const& kdtree_backend::kdtree_file() const { return _file; }

hpc::kdtree<real_type> const& kdtree_backend::kdtree() const
{
    // Debug: Print kdtree status when accessed
    // std::cout << "DEBUG: Accessing kdtree - bounds size: " <<
    // _kdt.bounds().size()
    //           << ", splits size: " << _kdt.splits().size() << std::endl;
    // if (_kdt.bounds().size() > 0) {
    //    std::cout << "DEBUG: kdtree dimensions: " << _kdt.n_dims()
    //              << ", branches: " << _kdt.n_branches()
    //              << ", leafs: " << _kdt.n_leafs() << std::endl;
    // } else {
    //    std::cout << "DEBUG: kdtree appears to be empty/uninitialized" <<
    //    std::endl;
    // }
    return _kdt;
}

std::vector<unsigned> const& kdtree_backend::cell_counts() const { return _cell_cnts; }

std::vector<unsigned long long> const& kdtree_backend::cell_offs() const { return _cell_offs; }

// Debug methods for inspecting kdtree structure
bool kdtree_backend::is_kdtree_loaded() const
{
    return _kdt.bounds().size() > 0 || _kdt.splits().size() > 0;
}

unsigned kdtree_backend::get_kdtree_dims() const { return _kdt.n_dims(); }

unsigned kdtree_backend::get_kdtree_branches() const { return _kdt.n_branches(); }

unsigned kdtree_backend::get_kdtree_leafs() const { return _kdt.n_leafs(); }

std::string kdtree_backend::debug_kdtree_info() const
{
    std::ostringstream oss;
    oss << "Kdtree info: bounds=" << _kdt.bounds().size() << ", splits=" << _kdt.splits().size();
    if (_kdt.bounds().size() > 0)
    {
        oss << ", dims=" << _kdt.n_dims() << ", branches=" << _kdt.n_branches()
            << ", leafs=" << _kdt.n_leafs() << ", cells=" << _kdt.n_cells();
    }
    else
    {
        oss << " [EMPTY/UNINITIALIZED]";
    }
    return oss.str();
}

} // namespace backends
} // namespace tao
