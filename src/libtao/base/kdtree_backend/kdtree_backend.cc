#include "kdtree_backend.hh"
#include <sstream>

namespace tao {
   namespace backends {

      kdtree_backend::kdtree_backend( hpc::mpi::comm const& comm )
         : _kdt( hpc::mpi::comm::self ),
           _snap( std::numeric_limits<unsigned>::max() ),
           _comm( &comm )
      {
      }

      kdtree_backend::kdtree_backend( hpc::fs::path const& fn,
                                      hpc::mpi::comm const& comm )
         : _kdt( hpc::mpi::comm::self ),
           _snap( std::numeric_limits<unsigned>::max() ),
           _comm( &comm )
      {
         open( fn );
      }

      // For some backends there are different fields expected.
      // For example, the kdtree requires a field subtree_count
      void
      kdtree_backend::add_conditional_fields(query<real_type> &qry) {
         qry.add_output_field("subtree_count");
      }

      hpc::mpi::comm const&
      kdtree_backend::comm() const
      {
         return *_comm;
      }

       herr_t
       static getFields( hid_t g_id, const char *name, const H5L_info_t *info, void *op_data) {
           std::vector<std::string> *fields = (std::vector<std::string>*)op_data;
           std::string name_lowercase=name;
           std::transform(name_lowercase.begin(), name_lowercase.end(), name_lowercase.begin(), ::tolower);
           fields->push_back(name);
           return 0;
       }

      void
      kdtree_backend::connect(const cli_dict&  global_cli_dict) {
         std::string fn = global_cli_dict._dataset;
         std::string xmlfile=fn;
         xmlfile.erase(xmlfile.length()-3,3);
         xmlfile.append(".xml");
         xml_dict h5xml;
         h5xml.read(xmlfile ,"/sageinput");
         // Iterate over the module nodes.
         xpath_node_set fields = h5xml.get_nodes("/sageinput/Field");
         std::map<std::string, int> _field_order;

         _field_types.clear();

         for (const xpath_node *it = fields.begin(); it != fields.end(); ++it) {
               xml_node cur = it->node();

               std::string element = cur.name();
               std::string field_str = cur.text().as_string();
               to_lower(field_str);
               std::string label = cur.attribute("label").value();
               std::string description = cur.attribute("description").value();
               int order = cur.attribute("order").as_int();
               std::string units = cur.attribute("units").value();
               std::string group = cur.attribute("group").value();
               std::string type_str = cur.attribute("Type").value();
               //std::cout << "Field[" << field_str << "]=" << label << "," << type_str << std::endl;
               typename batch<real_type>::field_value_type type;
               if (type_str == "int" || type_str == "short")
                  type = batch<real_type>::INTEGER;
               else if (type_str == "long long")
                  type = batch<real_type>::LONG_LONG;
               else if (type_str == "float")
                  type = batch<real_type>::DOUBLE;
               else {
                  EXCEPT(0, "Unknown field type for field '", field_str, "': ", type_str);
               }
               _field_types.emplace(field_str, type);
               _field_map[field_str] = field_str;
               _field_description[field_str] = description;
               _field_units[field_str] = units;
               _field_group[field_str] = group;
               _field_order[field_str] = order;
               //_field_order[field_str] = order;
               //std::cout << "Add ["<<field_str<<"]="<<_field_description[field_str]<<std::endl;

         }

         // Get it via the meta-data

          this->open(fn);
          hpc::h5::group data = kdtree_file().group("data");
          std::vector<std::string> Hfields;
          H5Lvisit(data.id(),H5_INDEX_CRT_ORDER,H5_ITER_NATIVE,getFields,&Hfields);

         //  for (std::string field : Hfields)
         //  {
         //      auto ds = data.dataset(field);

         //      std::cout << field << " dt="<<ds.datatype()<<","<<ds.type_class()<<std::endl;
         //  }

                  // Before finishing, insert the known mapping conversions.
                  // TODO: Generalise.
                  this->_field_map["pos_x"] = "posx";
                  this->_field_map["pos_y"] = "posy";
                  this->_field_map["pos_z"] = "posz";
                  this->_field_map["vel_x"] = "velx";
                  this->_field_map["vel_y"] = "vely";
                  this->_field_map["vel_z"] = "velz";
                  this->_field_map["snapshot"] = "snapnum";
                  this->_field_map["global_index"] = "global_index";
                  this->_field_map["global_tree_id"] = "globaltreeid";
                  this->_field_map["localgalaxyid"] = "localgalaxyid";

                  // RS Maybe for SED?

                  this->_field_map["subtree_count"] = "subtree_count";

                  // Add calculated types.
                  this->_field_types.emplace("redshift_cosmological", batch<real_type>::DOUBLE);
                  this->_field_types.emplace("column_density", batch<real_type>::DOUBLE);
                  this->_field_types.emplace("redshift_observed", batch<real_type>::DOUBLE);
                  this->_field_types.emplace("ra", batch<real_type>::DOUBLE);
                  this->_field_types.emplace("dec", batch<real_type>::DOUBLE);
                  this->_field_types.emplace("distance", batch<real_type>::DOUBLE);
                  this->_field_types.emplace("distance_from_beam", batch<real_type>::DOUBLE);
                  this->_field_types.emplace("sfr", batch<real_type>::DOUBLE);

#ifndef NDEBUG
                  // In debug mode we add some custom types.
                  this->_field_map["original_x"] = "posx";
                  this->_field_map["original_y"] = "posy";
                  this->_field_map["original_z"] = "posz";
                  this->_field_types.emplace("original_x", batch<real_type>::DOUBLE);
                  this->_field_types.emplace("original_y", batch<real_type>::DOUBLE);
                  this->_field_types.emplace("original_z", batch<real_type>::DOUBLE);
#endif

                  // Make sure we have all the essential fields available. Do this by
                  // checking that all the mapped fields exist in the field types.
                  //if (this->_init_tbls) {
                  //    for (const auto &item : this->_field_map)
                  //        EXCEPT(hpc::has(this->_field_types, item.second), "Database is missing essential field: ",
                  //               item.second);
                  //}
      }

      void
      kdtree_backend::close()
      {
      }

      void
      kdtree_backend::open( hpc::fs::path const& fn )
      {
         LOGBLOCKD( "Opening kdtree backend: ", fn.native() );
         close();
         _file.open( fn.native(), H5F_ACC_RDONLY );
         _snap = std::numeric_limits<unsigned>::max();

         _lc_mem_type.compound( sizeof(lightcone_data) );
         _lc_mem_type.insert( hpc::h5::datatype::native_double, "x", 0 );
         _lc_mem_type.insert( hpc::h5::datatype::native_double, "y", sizeof(double) );
         _lc_mem_type.insert( hpc::h5::datatype::native_double, "z", 2*sizeof(double) );
         _lc_mem_type.insert( hpc::h5::datatype::native_ullong, "global_index", HOFFSET( lightcone_data, gidx ) );
         _lc_mem_type.insert( hpc::h5::datatype::native_uint, "subsize", HOFFSET( lightcone_data, subsize ) );
      }

      void chround(char* a,int ndigits) {
        	int i;
        	int len=0;
        	for (i=0;i<1000;i++) {
                	if (a[i]==0)
                        	break;
                	len++;
        	}
        	int count_ndigits=0;
        	for (i=0;i<len;i++) {
                	if (a[i]==' ') continue;
                	if (a[i]=='.') continue;
                	if (a[i]=='0' && count_ndigits==0) continue;
                	count_ndigits++;
                	if (count_ndigits==ndigits) break;
        	}
        	char last = a[i];
        	bool up = false;
        	bool done = false;
        	for (int j=i+1;j<len;j++) {
                	switch(a[j]) {
                	case '0': done = true; break;
                	case '1': done = true; break;
                	case '2': done = true; break;
                	case '3': done = true; break;
                	case '4': done = true; break;
                	case '5': up = true; done = true; break;
                	case '6': up = true; done = true; break;
                	case '7': up = true; done = true; break;
                	case '8': up = true; done = true; break;
                	case '9': up = true; done = true; break;
                	}
                	if (done)
                        	break;
        	}
        	int endofstring_index = i+1;
        	if (up) {
                	done = false;
                	while (!done) {
                        	switch(a[i]) {
                        	case '0': a[i] = '1'; done = true; break;
                        	case '1': a[i] = '2'; done = true;  break;
                        	case '2': a[i] = '3'; done = true;  break;
                        	case '3': a[i] = '4'; done = true;  break;
                        	case '4': a[i] = '5'; done = true;  break;
                        	case '5': a[i] = '6'; done = true;  break;
                        	case '6': a[i] = '7'; done = true;  break;
                        	case '7': a[i] = '8'; done = true;  break;
                        	case '8': a[i] = '9'; done = true;  break;
                        	case '9': a[i] = '0'; break;
                        	}
                        	i--;
                        	if (i<0 && !done) {
					/*
 					* will need to shift all the characters one to the right and put a '1' at the start
 					*/

                                	for (int j=endofstring_index;j>0;j--) {
                                        	a[j] = a[j-1];
                                	}
                                	a[0] = '1';
                                	endofstring_index++;
                                	done = true;
                        	}
                	}
        	}
        	a[endofstring_index] = 0;
     }

      tao::simulation const*
      kdtree_backend::load_simulation()
      {
         LOGBLOCKD( "Loading simulation." );

         std::vector<real_type> redshifts = _file.read<std::vector<real_type> >("snapshot_redshifts");
         std::map<int, real_type> zs;
         std::cout.precision(17);
         for( size_t ii = 0; ii < redshifts.size(); ++ii ) {
            zs[ii] = redshifts[ii];
            char chredshift[1000];
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(12) << zs[ii];
            std::string temp = oss.str();
            strncpy(chredshift, temp.c_str(), 999);
            chredshift[999] = '\0';
            chround(chredshift, 6);
            //std::cout << "replace " << zs[ii] << " with " << chredshift << std::endl;
	         zs[ii] = atof(chredshift);
         }

         auto* sim = new tao::simulation(
            _file.read<real_type>( "cosmology/box_size" ),
            _file.read<real_type>( "cosmology/hubble" ),
            _file.read<real_type>( "cosmology/omega_m" ),
            _file.read<real_type>( "cosmology/omega_l" ),
            zs
            );
         set_simulation( sim );
         return sim;
      }

      void
      kdtree_backend::load_snapshot( unsigned snap )
      {
         if( snap != _snap )
         {


            std::string name = make_snapshot_name( snap );
            LOGBLOCKD( "Loading kdtree snapshot: ", name );
            ASSERT( _file.is_open(), "Have not opened an HDF5 kdtree file." );

            hpc::h5::group grp = _file.group( name );
            grp >> _kdt;
            grp.dataset( "cell_counts" ) >> _cell_cnts;
            grp.dataset( "cell_offs" ) >> _cell_offs;
            LOGILN( "Loading kdtree snapshot: ", name );
            
            // Debug: Print kdtree structure information
            LOGILN( "Kdtree loaded - Dimensions: ", _kdt.n_dims() );
            LOGILN( "Kdtree loaded - Branches: ", _kdt.n_branches() );
            LOGILN( "Kdtree loaded - Leafs: ", _kdt.n_leafs() );
            LOGILN( "Kdtree loaded - Total cells: ", _kdt.n_cells() );
            LOGILN( "Kdtree loaded - Bounds size: ", _kdt.bounds().size() );
            LOGILN( "Kdtree loaded - Splits size: ", _kdt.splits().size() );
            
            _snap = snap;
         }
      }

      real_type
      kdtree_backend::get_max_smoothing()
      {
         ASSERT( 0 );
      }

      std::string
      kdtree_backend::make_snapshot_name( unsigned snap ) const
      {
         return std::string( "lightcone/snapshot") + hpc::index_string( snap, 3 );
      }

      kdtree_backend::box_galaxy_iterator
      kdtree_backend::galaxy_begin( query<real_type>& qry,
                                    box<real_type> const& box,
                                    batch<real_type>* bat,
                                    filter const* flt )
      {
         load_snapshot( box.snapshot() );
         return box_galaxy_iterator( *this, box, bat, box.snapshot() );
      }

      kdtree_backend::box_galaxy_iterator
      kdtree_backend::galaxy_end( query<real_type>& qry,
                                  box<real_type> const& box ) const
      {
         return box_galaxy_iterator( *this );
      }

      kdtree_backend::lightcone_galaxy_iterator
      kdtree_backend::galaxy_begin( query<real_type>& qry,
                                    lightcone const& lc,
                                    batch<real_type>* bat,
                                    filter const* flt )
      {
          kdtree_backend::lightcone_galaxy_iterator result = lightcone_galaxy_iterator( lc, *this, qry, bat, flt );
         return result;
      }

      kdtree_backend::lightcone_galaxy_iterator
      kdtree_backend::galaxy_end( query<real_type> const& qry,
                                  lightcone const& lc ) const
      {
         return lightcone_galaxy_iterator();
      }

      kdtree_backend::tile_galaxy_iterator
      kdtree_backend::galaxy_begin( query<real_type>& qry,
                                    tile<real_type> const& tile,
                                    batch<real_type>* bat,
                                    filter const* flt,
                                    bool first )
      {
         return tile_galaxy_iterator( *this, tile, tile.lightcone(), bat, _snap );
      }

      kdtree_backend::tile_galaxy_iterator
      kdtree_backend::galaxy_end( query<real_type> const& qry,
                                  tile<real_type> const& tile ) const
      {
         return tile_galaxy_iterator( *this );
      }

      void
      kdtree_backend::load_lightcone_data( unsigned cell,
                                           std::vector<lightcone_data>& data ) const
      {
         unsigned n_elems = _cell_cnts[cell];
         unsigned long long offs = _cell_offs[cell];
         data.resize( n_elems );
         _file.dataset( "lightcone/data" ).read( data.data(), _lc_mem_type, n_elems, offs );
      }

      hpc::h5::file const&
      kdtree_backend::kdtree_file() const
      {
         return _file;
      }

      hpc::kdtree<real_type> const&
      kdtree_backend::kdtree() const
      {
         // Debug: Print kdtree status when accessed
         // std::cout << "DEBUG: Accessing kdtree - bounds size: " << _kdt.bounds().size() 
         //           << ", splits size: " << _kdt.splits().size() << std::endl;
         // if (_kdt.bounds().size() > 0) {
         //    std::cout << "DEBUG: kdtree dimensions: " << _kdt.n_dims() 
         //              << ", branches: " << _kdt.n_branches() 
         //              << ", leafs: " << _kdt.n_leafs() << std::endl;
         // } else {
         //    std::cout << "DEBUG: kdtree appears to be empty/uninitialized" << std::endl;
         // }
         return _kdt;
      }

      std::vector<unsigned> const&
      kdtree_backend::cell_counts() const
      {
         return _cell_cnts;
      }

      std::vector<unsigned long long> const&
      kdtree_backend::cell_offs() const
      {
         return _cell_offs;
      }

      // Debug methods for inspecting kdtree structure
      bool
      kdtree_backend::is_kdtree_loaded() const
      {
         return _kdt.bounds().size() > 0 || _kdt.splits().size() > 0;
      }

      unsigned
      kdtree_backend::get_kdtree_dims() const
      {
         return _kdt.n_dims();
      }

      unsigned
      kdtree_backend::get_kdtree_branches() const
      {
         return _kdt.n_branches();
      }

      unsigned
      kdtree_backend::get_kdtree_leafs() const
      {
         return _kdt.n_leafs();
      }

      std::string
      kdtree_backend::debug_kdtree_info() const
      {
         std::ostringstream oss;
         oss << "Kdtree info: bounds=" << _kdt.bounds().size() 
             << ", splits=" << _kdt.splits().size();
         if (_kdt.bounds().size() > 0) {
            oss << ", dims=" << _kdt.n_dims() 
                << ", branches=" << _kdt.n_branches() 
                << ", leafs=" << _kdt.n_leafs()
                << ", cells=" << _kdt.n_cells();
         } else {
            oss << " [EMPTY/UNINITIALIZED]";
         }
         return oss.str();
      }

   }
}
