#ifndef tao_modules_hdf5_hh
#define tao_modules_hdf5_hh

#include <libhpc/libhpc.hh>
#include "libtao/base/module.hh"
#include "libtao/base/batch.hh"
#include "libtao/base/filter.hh"

constexpr hsize_t DEFAULT_CHUNK_SIZE = 100000;

// Helper function to write string attribute to HDF5 dataset
static void write_dataset_string_attribute(hid_t dset_id, const std::string& attr_name, const std::string& value) {
    if (value.empty()) return;

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
}

namespace tao {
   namespace modules {

      template< class Backend >
      class hdf5
         : public module<Backend>
      {
      public:

         typedef Backend backend_type;
         typedef module<backend_type> module_type;

         static
         module_type*
         factory( const std::string& name,
                  pugi::xml_node base )
         {
            return new hdf5( name, base );
         }

    static void process_cli_options(const cli_dict& global_cli_dict) 
    {
        // Access global CLI options
        std::string outfile = global_cli_dict._outfile;
        std::string outdir = global_cli_dict._outdir;
        
        // Create output directory if it doesn't exist
        if (!boost::filesystem::exists(outdir)) {
            boost::filesystem::create_directories(outdir);
        }
        // gathering results from all ranks
        if (mpi::comm::world.rank() == 0) {
            // If we are the root rank, we need to gather results from all ranks.
            LOGILN("Gathering results from all ",mpi::comm::world.size()," ranks.");
            std::string file0=global_cli_dict._outdir + "/"+_construct_filename_static(global_cli_dict._outfile, -1);
            // Create a new file to store the gathered results.
            // This will overwrite the existing file0 if it exists.
            // If the file does not exist, it will be created.
            // If the file exists, it will be truncated.
            hpc::h5::file _file0(file0, H5F_ACC_TRUNC);
            if (!_file0.is_open()) {
                std::cerr << "Error: Could not open file for gathering results." << std::endl;
                throw silent_terminate();
            }

            int rank = 0;
            std::vector<std::string> result_files;
            std::string file1=global_cli_dict._outdir + "/"+_construct_filename_static(global_cli_dict._outfile, rank);
            while (boost::filesystem::exists(file1)) {
                result_files.push_back(file1);
                rank++;
                file1 = global_cli_dict._outdir + "/" + _construct_filename_static(global_cli_dict._outfile, rank);
            }
            if (!result_files.empty()) {
                // Get field list from first rank file by reading dataset names
                std::vector<std::string> field_names;
                {
                    hpc::h5::file first_h5file(result_files[0], H5F_ACC_RDONLY);
                    if (!first_h5file.is_open()) {
                        std::cerr << "Error: Could not open file for gathering results." << std::endl;
                        throw silent_terminate();
                    }

                    // Iterate through all datasets in the file
                    H5G_info_t group_info;
                    H5Gget_info(first_h5file.id(), &group_info);

                    for (hsize_t i = 0; i < group_info.nlinks; i++) {
                        char name_buf[256];
                        H5Lget_name_by_idx(first_h5file.id(), ".", H5_INDEX_NAME, H5_ITER_NATIVE, i, name_buf, sizeof(name_buf), H5P_DEFAULT);
                        std::string name(name_buf);

                        // Assume all links are datasets (we're iterating at the file root level)
                        field_names.push_back(name);
                    }
                }

                int file_index=0;
                for (const auto& field_name : field_names) {
                    bool first_file = true;
                    for (const auto& file : result_files) {
                        hpc::h5::file h5file(file, H5F_ACC_RDONLY);
                        if (!h5file.is_open()) {
                            std::cerr << "Error: Could not open file for gathering results." << std::endl;
                            throw silent_terminate();
                        }
                        std::string dataset_name = field_name;
                        //to_lower(dataset_name);
                        // Open source dataset
                        hpc::h5::dataset src_dset(h5file, dataset_name);
                        
                        // Get datatype and dataspace from source
                        hpc::h5::datatype dtype = src_dset.datatype();
                        
                        // Read data from source
                        hsize_t size = src_dset.extent();
                        std::vector<char> buffer(size * dtype.size());
                        src_dset.read(buffer.data(), dtype);

                        // Create destination dataset with same properties
                        
                        // Create destination dataset with same properties
                        if (first_file) {
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
                            if (src_id >= 0 && dst_id >= 0) {
                                // Copy Description attribute
                                if (H5Aexists(src_id, "Description") > 0) {
                                    hid_t attr_id = H5Aopen(src_id, "Description", H5P_DEFAULT);
                                    hid_t type_id = H5Aget_type(attr_id);
                                    size_t attr_size = H5Tget_size(type_id);
                                    std::vector<char> attr_buf(attr_size + 1, 0);
                                    H5Aread(attr_id, type_id, attr_buf.data());
                                    H5Aclose(attr_id);

                                    write_dataset_string_attribute(dst_id, "Description", std::string(attr_buf.data()));
                                    H5Tclose(type_id);
                                }

                                // Copy Units attribute
                                if (H5Aexists(src_id, "Units") > 0) {
                                    hid_t attr_id = H5Aopen(src_id, "Units", H5P_DEFAULT);
                                    hid_t type_id = H5Aget_type(attr_id);
                                    size_t attr_size = H5Tget_size(type_id);
                                    std::vector<char> attr_buf(attr_size + 1, 0);
                                    H5Aread(attr_id, type_id, attr_buf.data());
                                    H5Aclose(attr_id);

                                    write_dataset_string_attribute(dst_id, "Units", std::string(attr_buf.data()));
                                    H5Tclose(type_id);
                                }

                                H5Dclose(src_id);
                                H5Dclose(dst_id);
                            }

                            // Write data using chunking
                            write_chunked_dataset(dst_dset, buffer.data(), size, dtype);
                            
                            dst_dset.close();
                            src_dset.close();
                        } else {
                            hpc::h5::dataset dst_dset(_file0, dataset_name);
                            
                            // Get current size and extend
                            hsize_t current_size = dst_dset.extent();
                            dst_dset.set_extent(current_size + size);
                            
                            // Write data using chunking
                            write_chunked_dataset(dst_dset, buffer.data(), size, dtype, current_size);
                            
                            dst_dset.close();
                            src_dset.close();
                        }
                        h5file.close();
                        file_index++;
                    }
                }
                for (const auto& file : result_files) {
                    boost::filesystem::remove(file);
                }
            }
            _file0.close();
            LOGILN("Output file: ", outfile);
            LOGILN("Output directory: ", outdir);
        }
    }

      public:

         hdf5( const std::string& name = std::string(),
               pugi::xml_node base = pugi::xml_node() )
            : module_type( name, base ),
              _chunk_size( 10000 ),
              _ready( false ),
              _records( 0 )
         {
         }

         virtual
         ~hdf5()
         {
         }

         ///
         ///
         ///
         virtual
         void
         initialise( const cli_dict& global_cli_dict,
                     boost::optional<boost::property_tree::ptree> checkpoint = boost::optional<boost::property_tree::ptree>() )
         {
            // Don't initialise if we're already doing so.
            if( this->_init )
               return;
            module_type::initialise( global_cli_dict, checkpoint );

            LOGILN( ".Initialising HDF5 module.", setindent( 2 ) );


            // Get our information.

            std::list<boost::optional<std::string>> field_labels;
            std::string filename=global_cli_dict._outfile;
            _fn = global_cli_dict._outdir + "/" + _construct_filename( filename );

            // Store input kdtree file path for attribute reading
            _input_kdtree_file = global_cli_dict._dataset;

            auto const &input_hdf5_file = global_cli_dict._dataset;

            // Get field list directly from HDF5 file /data/* datasets
            std::vector<std::string> available_fields;
            {
                hpc::h5::file input_file(input_hdf5_file, H5F_ACC_RDONLY);
                hpc::h5::group data_group;
                input_file.open_group("data", data_group);

                // Iterate through datasets in /data group
                H5G_info_t group_info;
                H5Gget_info(data_group.id(), &group_info);

                for (hsize_t i = 0; i < group_info.nlinks; i++) {
                    char name_buf[256];
                    H5Lget_name_by_idx(data_group.id(), ".", H5_INDEX_NAME, H5_ITER_NATIVE, i, name_buf, sizeof(name_buf), H5P_DEFAULT);
                    available_fields.push_back(std::string(name_buf));
                }
            }

            size_t nfields = available_fields.size();
            LOGILN( "Initialising HDF5 module. nfields from input HDF5 file: ", nfields, setindent( 2 ) );

            if (global_cli_dict._output_fields.size() == 0)
            {
                LOGILN( "No output fields specified via CLI, using all fields from the input hdf5 file." );
                for( auto const& field_name : available_fields )
                {
                    _fields.push_back( field_name );
                    _labels.push_back( field_name );  // Use field name as label
                    field_labels.push_back( field_name );
                    LOGILN( "CLI.Adding field ", field_name, " with label ", field_name );
                }
            }
            else
            {
                LOGILN( "Using output fields specified in the command line." );
                for( auto const& field : global_cli_dict._output_fields )
                {
                    LOGILN( "CLI.Adding field ", field );
                    _fields.push_back( field );
                    _labels.push_back( field );
                }
            }

             // Note: Alias handling (X, Y, Z, COLOUR) has been removed with XML sidecar removal
             // Note: XML sidecar metadata creation has been removed - using HDF5 attributes instead

            // Open the file. Truncate if we are not reloading.
            _file.open( _fn, H5F_ACC_TRUNC );
            _records = 0;
            _ready = false;

            // Get the filter from the lightcone module
            _filt = this->template attribute<tao::filter const*>( "filter" );

            LOGILN( "Done.", setindent( -2 ) );
            // LOGILN(" Exit Now !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
            // exit(0);
         }

         static std::string
         _construct_filename_static( const std::string& base_name, int rank )
         {
            // Construct the filename.
            std::string filename = base_name;
            size_t last_dot = filename.find_last_of('.');
            if (last_dot != std::string::npos) {
                std::string ext = filename.substr(last_dot);
                if (ext == ".h5" || ext == ".hdf5") {
                    filename = filename.substr(0, last_dot);
                }
            }
            if( rank >= 0 )
            {
                char chbuff[10];
                snprintf(chbuff, sizeof(chbuff), "%05d", rank);
                std::string rank_str(chbuff);
               filename += "." + rank_str;
            }
            filename += ".h5";
            return filename;
         }

         std::string
         _construct_filename( const std::string& base_name )
         {
            // Construct the filename.
            std::string filename = base_name;
            size_t last_dot = filename.find_last_of('.');
            if (last_dot != std::string::npos) {
                std::string ext = filename.substr(last_dot);
                if (ext == ".h5" || ext == ".hdf5") {
                    filename = filename.substr(0, last_dot);
                }
            }
            if( mpi::comm::world.size() > 1 )
            {
               filename += "." + mpi::rank_string();
            }
            filename += ".h5";
            return filename;
         }

         std::string
         _encode( std::string _toencode_string )
         {
            std::map<char, std::string> transformations;
            transformations['&']  = std::string("_");
            transformations[' ']  = std::string("_");
            transformations['\''] = std::string("_");
            transformations['"']  = std::string("_");
            transformations['>']  = std::string("_");
            transformations['<']  = std::string("_");
            transformations['/']  = std::string("_");
            transformations['(']  = std::string("");
            transformations[')']  = std::string("");


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
         virtual
         void
         execute()
         {
            ASSERT( this->parents().size() == 1 );

            // Grab the batch from the parent object.
            const tao::batch<real_type>& bat = this->parents().front()->batch();
             // count the number of entries just once - it will be the same for all the fields
             uint32_t thistime=0;
             for( auto it = _filt->begin( bat, true ); it != _filt->end( bat, true ); ++it)
             {
                 thistime++;
             }

            // Create datasets for the field names. Note that I can only
            // do this when I have a galaxy object, hence doing it once
            // here.
            bool alt=true;
             h5::dataspace mem_space;
             mem_space.create( 1 );
             mem_space.select_all();
            if( !_ready )
            {

            	auto lblit = _labels.cbegin();

                for( const auto& field : _fields )
                {
                    hpc::to_lower( (std::string&)field );   // input field name of hdf5 must be lowercase
                   h5::datatype dtype = _field_type( bat, field );
                   h5::dataspace dspace;
                   dspace.create( 1, true );
                   h5::property_list props( H5P_DATASET_CREATE );
                   props.set_chunk_size( _chunk_size );
                   //props.set_deflate();
                    std::string field_name = _encode( *lblit ); // use label as field name
                   if (lblit->empty())
                       field_name = field;
                    //h5::dataset* dset = new h5::dataset( _file, field, dtype, dspace, none, false, props );
                   //std::cout << "CREATING HDF5 field "<<field_name<<" with dtype="<<_datatypeAsTAO(dtype)<<std::endl;
                   h5::dataset* dset = new hpc::h5::dataset( _file, field_name, dtype, dspace, props );
                   dset->set_extent(thistime);
                   _dsets.push_back( std::unique_ptr<hpc::h5::dataset>( dset ) );

                   // Copy HDF5 attributes from input kdtree file
                   if (!_input_kdtree_file.empty()) {
                       hid_t input_file_id = H5Fopen(_input_kdtree_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
                       if (input_file_id >= 0) {
                           std::string input_dataset_path = "data/" + field;
                           hid_t input_dset_id = H5Dopen2(input_file_id, input_dataset_path.c_str(), H5P_DEFAULT);
                           hid_t output_dset_id = H5Dopen2(_file.id(), field_name.c_str(), H5P_DEFAULT);

                           if (input_dset_id >= 0 && output_dset_id >= 0) {
                               // Copy Description attribute
                               if (H5Aexists(input_dset_id, "Description") > 0) {
                                   hid_t attr_id = H5Aopen(input_dset_id, "Description", H5P_DEFAULT);
                                   hid_t type_id = H5Aget_type(attr_id);
                                   size_t attr_size = H5Tget_size(type_id);
                                   std::vector<char> attr_buf(attr_size + 1, 0);
                                   H5Aread(attr_id, type_id, attr_buf.data());
                                   H5Aclose(attr_id);
                                   write_dataset_string_attribute(output_dset_id, "Description", std::string(attr_buf.data()));
                                   H5Tclose(type_id);
                               }

                               // Copy Units attribute
                               if (H5Aexists(input_dset_id, "Units") > 0) {
                                   hid_t attr_id = H5Aopen(input_dset_id, "Units", H5P_DEFAULT);
                                   hid_t type_id = H5Aget_type(attr_id);
                                   size_t attr_size = H5Tget_size(type_id);
                                   std::vector<char> attr_buf(attr_size + 1, 0);
                                   H5Aread(attr_id, type_id, attr_buf.data());
                                   H5Aclose(attr_id);
                                   write_dataset_string_attribute(output_dset_id, "Units", std::string(attr_buf.data()));
                                   H5Tclose(type_id);
                               }

                               if (input_dset_id >= 0) H5Dclose(input_dset_id);
                               if (output_dset_id >= 0) H5Dclose(output_dset_id);
                           }

                           H5Fclose(input_file_id);
                       }
                   }

                   // Dump first chunk.
                    auto val = bat.field( field );

                   unsigned ii = 0;
                    switch( std::get<2>( val ) )
                    {
                        case tao::batch<real_type>::STRING: {
                            const hpc::view<std::vector<std::string> > thedata = bat.scalar<std::string>(field);
                            for (auto it = _filt->begin(bat, true); it != _filt->end(bat, true); ++it, ++ii) {
                                dspace = dset->dataspace();
                                dspace.select_one(ii);
                                _write_field_STRING(bat, *it, field, *dset, dspace, mem_space, thedata);
                            }
                        }
                            break;
                        case tao::batch<real_type>::DOUBLE: {
                            const hpc::view<std::vector<double> > thedata = bat.scalar<double>(field);
                            for (auto it = _filt->begin(bat, true); it != _filt->end(bat, true); ++it, ++ii) {
                                dspace = dset->dataspace();
                                dspace.select_one(ii);
                                _write_field_DOUBLE(bat, *it, field, *dset, dspace, mem_space, thedata);
                            }
                        }
                            break;
                        case tao::batch<real_type>::FLOAT: {
                            const hpc::view<std::vector<float> > thedata = bat.scalar<float>(field);
                            for (auto it = _filt->begin(bat, true); it != _filt->end(bat, true); ++it, ++ii) {
                                dspace = dset->dataspace();
                                dspace.select_one(ii);
                                _write_field_FLOAT(bat, *it, field, *dset, dspace, mem_space, thedata);
                            }
                        }
                        break;
                        case tao::batch<real_type>::INTEGER: {
                            const hpc::view<std::vector<int> > thedata = bat.scalar<int>(field);
                            for (auto it = _filt->begin(bat, true); it != _filt->end(bat, true); ++it, ++ii) {
                                dspace = dset->dataspace();
                                dspace.select_one(ii);
                                _write_field_INTEGER(bat, *it, field, *dset, dspace, mem_space, thedata);
                            }
                        }
                            break;
                        case tao::batch<real_type>::UNSIGNED_LONG_LONG: {
                            const hpc::view<std::vector<unsigned long long> > thedata = bat.scalar<unsigned long long>(
                                    field);
                            for (auto it = _filt->begin(bat, true); it != _filt->end(bat, true); ++it, ++ii) {
                                dspace = dset->dataspace();
                                dspace.select_one(ii);
                                _write_field_UNSIGNED_LONG_LONG(bat, *it, field, *dset, dspace, mem_space, thedata);
                            }
                        }
                            break;
                        case tao::batch<real_type>::LONG_LONG: {
                            const hpc::view<std::vector<long long> > thedata = bat.scalar<long long>(field);
                            for (auto it = _filt->begin(bat, true); it != _filt->end(bat, true); ++it, ++ii) {
                                dspace = dset->dataspace();
                                dspace.select_one(ii);
                                _write_field_LONG_LONG(bat, *it, field, *dset, dspace, mem_space, thedata);
                            }
                        }
                            break;
                        default:
                            ASSERT( 0 );
                    }
                   lblit++;

                   // Set number written for first batch.
                   _records = ii;

                }

                // Flag as complete.
                _ready = true;
            }
            else if (alt) {
                auto dset_it = _dsets.begin();
                uint32_t offset=_records;
                for( const auto& field : _fields )
                {
                    h5::dataset* dset = (*dset_it++).get();
                    {
                        h5::dataspace dspace(*dset);
                        dset->set_extent(dspace.size() + thistime);
                    }
                    h5::dataspace dspace(*dset);

                    // Dump current chunk.
                    auto val = bat.field( field );

                    unsigned ii = offset;

                    switch( std::get<2>( val ) )
                    {
                        case tao::batch<real_type>::STRING: {
                            const hpc::view<std::vector<std::string> > thedata = bat.scalar<std::string>(field);
                            for (auto it = _filt->begin(bat, true); it != _filt->end(bat, true); ++it, ++ii) {
                                dspace.select_one(ii);
                                _write_field_STRING(bat, *it, field, *dset, dspace, mem_space, thedata);
                            }
                        }
                            break;
                        case tao::batch<real_type>::DOUBLE: {
                            const hpc::view<std::vector<double> > thedata = bat.scalar<double>(field);
                            for (auto it = _filt->begin(bat, true); it != _filt->end(bat, true); ++it, ++ii) {
                                //dspace = dset->dataspace();
                                dspace.select_one(ii);
                                _write_field_DOUBLE(bat, *it, field, *dset, dspace, mem_space, thedata);
                            }
                        }
                            break;
                        case tao::batch<real_type>::FLOAT: {
                            const hpc::view<std::vector<float> > thedata = bat.scalar<float>(field);
                            for (auto it = _filt->begin(bat, true); it != _filt->end(bat, true); ++it, ++ii) {
                                //dspace = dset->dataspace();
                                dspace.select_one(ii);
                                _write_field_FLOAT(bat, *it, field, *dset, dspace, mem_space, thedata);
                            }
                        }
                        break;
                        case tao::batch<real_type>::INTEGER: {
                            const hpc::view<std::vector<int> > thedata = bat.scalar<int>(field);
                            for (auto it = _filt->begin(bat, true); it != _filt->end(bat, true); ++it, ++ii) {
                                dspace.select_one(ii);
                                _write_field_INTEGER(bat, *it, field, *dset, dspace, mem_space, thedata);
                            }
                        }
                            break;
                        case tao::batch<real_type>::UNSIGNED_LONG_LONG: {
                            const hpc::view<std::vector<unsigned long long> > thedata = bat.scalar<unsigned long long>(field);
                            for (auto it = _filt->begin(bat, true); it != _filt->end(bat, true); ++it, ++ii) {
                                dspace.select_one(ii);
                                _write_field_UNSIGNED_LONG_LONG(bat, *it, field, *dset, dspace, mem_space, thedata);
                            }
                        }
                            break;
                        case tao::batch<real_type>::LONG_LONG: {
                            const hpc::view<std::vector<long long> > thedata = bat.scalar<long long>(field);
                            for (auto it = _filt->begin(bat, true); it != _filt->end(bat, true); ++it, ++ii) {
                                dspace.select_one(ii);
                                _write_field_LONG_LONG(bat, *it, field, *dset, dspace, mem_space, thedata);
                            }
                        }
                            break;
                        default:
                            ASSERT( 0 );
                    }

                }
                _records += thistime;
            }
            else {
               // Process each element.
               for( auto it = _filt->begin( bat, true ); it != _filt->end( bat, true ); ++it )
               {
                  // Write the fields.
                  auto dset_it = _dsets.begin();
                  for( const auto& field : _fields )
                  {
                     h5::dataset* dset = (*dset_it++).get();
                     hsize_t old_size;
                     {
                        h5::dataspace dspace( *dset );
                        old_size = dspace.size();
                        if( _records == old_size )
                        {
                           dset->set_extent( old_size + 1 );
                        }
                     }
                     h5::dataspace dspace( *dset );
                     dspace.select_one( _records );
                     _write_field( bat, *it, field, *dset, dspace );
                  }

		  // Increment records.
		  ++_records;
               }
            }

            // Flush the file here. I want to be sure that all records have
            // been written before potentially checkpointing.
            _file.flush();
         }

         virtual
         void
         log_metrics()
         {
            module_type::log_metrics();
            LOGILN( this->_name, " number of records written: ", _records );
         }

         virtual
         void
         do_checkpoint( boost::property_tree::ptree& pt )
         {
            std::string fn( _fn );
            std::replace( fn.begin(), fn.end(), '.', '-' );
            std::replace( fn.begin(), fn.end(), '/', '-' );
            pt.put( std::string( "hdf5." ) + fn, std::to_string( _records ) );
         }

      protected:

      // Inside the hdf5 class
        static bool _has_dataset(hpc::h5::file& file, const std::string& dataset_name) {
            // Use H5Lexists to check if the link exists
            htri_t link_exists = H5Lexists(file.id(), dataset_name.c_str(), H5P_DEFAULT);
            return (link_exists > 0);
        }

        static void write_chunked_dataset(hpc::h5::dataset& dst_dset, 
                                const char* buffer,
                                hsize_t size,
                                hpc::h5::datatype& dtype,
                                hsize_t current_size = 0) {
            const hsize_t chunk_size = DEFAULT_CHUNK_SIZE;
            const hsize_t num_chunks = (size + chunk_size - 1) / chunk_size;

            for (hsize_t chunk = 0; chunk < num_chunks; chunk++) {
                // Calculate chunk boundaries
                hsize_t chunk_start = chunk * chunk_size;
                hsize_t chunk_end = std::min(chunk_start + chunk_size, size);
                hsize_t chunk_length = chunk_end - chunk_start;

                // Create memory space for chunk
                hpc::h5::dataspace mem_space;
                mem_space.create(chunk_length, false);  // Create with exact size

                // Select hyperslab in file
                hpc::h5::dataspace file_space = dst_dset.dataspace();
                file_space.select_hyperslab(H5S_SELECT_SET,
                                  chunk_length,    // count
                                  current_size + chunk_start,  // start
                                  1,    // stride
                                  1);   // block

                // Write chunk
                dst_dset.write(buffer + (chunk_start * dtype.size()),
                            dtype,
                            mem_space,
                            file_space);
            }
        }

         h5::datatype
         _field_type( const tao::batch<real_type>& bat,
                      const std::string& field )
         {
            auto val = bat.field( field );
            switch( std::get<2>( val ) )
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
                  ASSERT( 0 );
            }
         }

         void
         _write_field( const tao::batch<real_type>& bat,
                       unsigned idx,
                       const std::string& field,
                       h5::dataset& dset,
                       h5::dataspace& dspace )
         {
            auto val = bat.field( field );
            h5::dataspace mem_space;
            mem_space.create( 1 );
            mem_space.select_all();
            switch( std::get<2>( val ) )
            {
               case tao::batch<real_type>::STRING:
               {
                  std::string data = bat.scalar<std::string>( field )[idx];
                  dset.write( data.c_str(), h5::datatype::string, mem_space, dspace );
                  break;
               }

               case tao::batch<real_type>::DOUBLE:
               {
                  double data = bat.scalar<double>( field )[idx];
                  dset.write( &data, h5::datatype::native_double, mem_space, dspace );
                  break;
               }

               case tao::batch<real_type>::INTEGER:
               {
                  int data = bat.scalar<int>( field )[idx];
                  dset.write( &data, h5::datatype::native_int, mem_space, dspace );
                  break;
               }

               case tao::batch<real_type>::UNSIGNED_LONG_LONG:
               {
                  unsigned long long data = bat.scalar<unsigned long long>( field )[idx];
                  dset.write( &data, h5::datatype::native_ullong, mem_space, dspace );
                  break;
               }

               case tao::batch<real_type>::LONG_LONG:
               {
                  long long data = bat.scalar<long long>( field )[idx];
                  dset.write( &data, h5::datatype::native_llong, mem_space, dspace );
                  break;
               }

               default:
                  ASSERT( 0 );
            }
         }

          void
          _write_field_STRING( const tao::batch<real_type>& bat,
                               unsigned idx,
                               const std::string& field,
                               h5::dataset& dset,
                               h5::dataspace& dspace,
                               h5::dataspace& mem_space,
                               const hpc::view<std::vector<std::string> >& thedata)
          {
              {
                  std::string data = thedata[idx];
                  dset.write( data.c_str(), h5::datatype::string, mem_space, dspace );
              }
          }

          void
          _write_field_DOUBLE( const tao::batch<real_type>& bat,
                        unsigned idx,
                        const std::string& field,
                        h5::dataset& dset,
                        h5::dataspace& dspace,
                        h5::dataspace& mem_space,
                        const hpc::view<std::vector<double> >& thedata)
          {
              {
                      double data = thedata[idx];
                      dset.write( &data, h5::datatype::native_double, mem_space, dspace );
              }
          }

          void
          _write_field_FLOAT( const tao::batch<real_type>& bat,
                               unsigned idx,
                               const std::string& field,
                               h5::dataset& dset,
                               h5::dataspace& dspace,
                               h5::dataspace& mem_space,
                               const hpc::view<std::vector<float> >& thedata)
                               {
             {
                 float data = thedata[idx];
                 dset.write( &data, h5::datatype::native_float, mem_space, dspace );
             }
                               }

          void
          _write_field_INTEGER( const tao::batch<real_type>& bat,
                               unsigned idx,
                               const std::string& field,
                               h5::dataset& dset,
                               h5::dataspace& dspace,
                               h5::dataspace& mem_space,
                                const hpc::view<std::vector<int> >& thedata)
          {
              {
                  int data = thedata[idx];
                  dset.write( &data, h5::datatype::native_int, mem_space, dspace );
              }
          }

          void
          _write_field_UNSIGNED_LONG_LONG( const tao::batch<real_type>& bat,
                                unsigned idx,
                                const std::string& field,
                                h5::dataset& dset,
                                h5::dataspace& dspace,
                                h5::dataspace& mem_space,
                                           const hpc::view<std::vector<unsigned long long> >& thedata)
          {
              {
                  long long data = thedata[idx];
                  dset.write( &data, h5::datatype::native_llong, mem_space, dspace );
              }
          }

          void
          _write_field_LONG_LONG( const tao::batch<real_type>& bat,
                                           unsigned idx,
                                           const std::string& field,
                                           h5::dataset& dset,
                                           h5::dataspace& dspace,
                                            h5::dataspace& mem_space,
                                  const hpc::view<std::vector<long long> >& thedata)
          {
              {
                  signed long long data = thedata[idx];
                  dset.write( &data, h5::datatype::native_llong, mem_space, dspace );
              }
          }

         typename tao::batch<tao::real_type>::field_value_type _datatypeAsTAO(hpc::h5::datatype datatype) {
            H5T_class_t dt1 = H5Tget_class(datatype.id());
            int dt1size = H5Tget_size(datatype.id());
            switch (dt1) {
            case H5T_class_t::H5T_INTEGER:
               if (dt1size == 4) {
                  return tao::batch<tao::real_type>::INTEGER;
               } else if (dt1size == 8) {
                  return tao::batch<tao::real_type>::LONG_LONG;
               }
               break;
            case H5T_class_t::H5T_FLOAT:
               if (dt1size == 4) {
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
         std::string _input_kdtree_file;  // Path to input kdtree file for attribute reading
         std::list<std::string> _fields;
         unsigned long long _records;
         std::list<std::string> _labels;
         std::list<std::unique_ptr<h5::dataset>> _dsets;
         hsize_t _chunk_size;
         bool _ready;
         tao::filter const* _filt;
      };

   }
}

#endif
