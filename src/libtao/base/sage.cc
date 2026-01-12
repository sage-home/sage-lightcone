#include <iostream>
#include <libhpc/h5/derive.hh>
#include "sage.hh"

namespace sage {

   using namespace hpc;

   void
   make_hdf5_types( h5::datatype& mem_type,
		    h5::datatype& file_type )
   {
      // Use native endian for both memory and file types (simpler, more efficient)
      // Changed from big-endian (std_i32be, ieee_f32be) to native endian
      h5::derive der( sizeof(galaxy) );
      der.add( h5::datatype::native_int,   HOFFSET( galaxy, snapshot ),            h5::datatype::native_int,   "snapnum" );
      der.add( h5::datatype::native_int,   HOFFSET( galaxy, type ),                h5::datatype::native_int,   "type" );
      der.add( h5::datatype::native_llong, HOFFSET( galaxy, galaxy_idx ),          h5::datatype::native_llong, "galaxy_index" );
      der.add( h5::datatype::native_llong, HOFFSET( galaxy, central_galaxy_idx ),  h5::datatype::native_llong, "central_galaxy_index" );
      der.add( h5::datatype::native_int,   HOFFSET( galaxy, sage_halo_idx ),       h5::datatype::native_int,   "sage_halo_index" );
      der.add( h5::datatype::native_int,   HOFFSET( galaxy, sage_tree_idx ),       h5::datatype::native_int,   "sage_tree_index" );
      der.add( h5::datatype::native_llong, HOFFSET( galaxy, simulation_halo_idx ), h5::datatype::native_llong, "simulation_halo_index" );
      der.add( h5::datatype::native_int,   HOFFSET( galaxy, local_index ),         h5::datatype::native_int,   "local_index" );
      der.add( h5::datatype::native_llong, HOFFSET( galaxy, global_index ),        h5::datatype::native_llong, "global_index" );
      der.add( h5::datatype::native_int,   HOFFSET( galaxy, descendant ),          h5::datatype::native_int,   "descendant" );
      der.add( h5::datatype::native_llong, HOFFSET( galaxy, global_descendant ),   h5::datatype::native_llong, "global_descendant" );
      der.add( h5::datatype::native_int,   HOFFSET( galaxy, subsize ),             h5::datatype::native_int,   "subsize" );
      der.add( h5::datatype::native_int,   HOFFSET( galaxy, merge_type ),          h5::datatype::native_int,   "merge_type" );
      der.add( h5::datatype::native_int,   HOFFSET( galaxy, merge_into_id ),       h5::datatype::native_int,   "merge_into_id" );
      der.add( h5::datatype::native_int,   HOFFSET( galaxy, merge_into_snapshot ), h5::datatype::native_int,   "merge_into_snapshot" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, dt ),                  h5::datatype::native_float, "dt" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, pos[0] ),              h5::datatype::native_float, "posx" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, pos[1] ),              h5::datatype::native_float, "posy" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, pos[2] ),              h5::datatype::native_float, "posz" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, vel[0] ),              h5::datatype::native_float, "velx" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, vel[1] ),              h5::datatype::native_float, "vely" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, vel[2] ),              h5::datatype::native_float, "velz" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, spin[0] ),             h5::datatype::native_float, "spin_x" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, spin[1] ),             h5::datatype::native_float, "spin_y" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, spin[2] ),             h5::datatype::native_float, "spin_z" );
      der.add( h5::datatype::native_int,   HOFFSET( galaxy, num_particles ),       h5::datatype::native_int,   "n_darkmatter_particles" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, mvir ),                h5::datatype::native_float, "virial_mass" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, central_mvir ),        h5::datatype::native_float, "central_galaxy_mvir" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, rvir ),                h5::datatype::native_float, "virial_radius" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, vvir ),                h5::datatype::native_float, "virial_velocity" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, vmax ),                h5::datatype::native_float, "max_velocity" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, vel_disp ),            h5::datatype::native_float, "velocity_dispersion" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, cold_gas ),            h5::datatype::native_float, "cold_gas" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, stellar_mass ),        h5::datatype::native_float, "stellar_mass" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, bulge_mass ),          h5::datatype::native_float, "bulge_mass" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, hot_gas ),             h5::datatype::native_float, "hot_gas" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, ejected_mass ),        h5::datatype::native_float, "ejected_mass" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, blackhole_mass ),      h5::datatype::native_float, "blackhole_mass" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, ics ),                 h5::datatype::native_float, "ics" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, metals_cold_gas ),     h5::datatype::native_float, "metals_cold_gas" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, metals_stellar_mass ), h5::datatype::native_float, "metals_stellar_mass" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, metals_bulge_mass ),   h5::datatype::native_float, "metals_bulge_mass" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, metals_hot_gas ),      h5::datatype::native_float, "metals_hot_gas" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, metals_ejected_mass ), h5::datatype::native_float, "metals_ejected_mass" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, metals_ics ),          h5::datatype::native_float, "metals_ics" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, sfr_disk ),            h5::datatype::native_float, "sfr_disk" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, sfr_bulge ),           h5::datatype::native_float, "sfr_bulge" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, sfr_disk_z ),          h5::datatype::native_float, "sfr_disk_z" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, sfr_bulge_z ),         h5::datatype::native_float, "sfr_bulge_z" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, disk_scale_radius ),   h5::datatype::native_float, "disk_scale_radius" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, cooling ),             h5::datatype::native_float, "cooling" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, heating ),             h5::datatype::native_float, "heating" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, quasar_mode_bh_accretion_mass ), h5::datatype::native_float, "quasar_mode_bh_accretion_mass" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, time_of_last_major_merger ), h5::datatype::native_float, "time_of_last_major_merger" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, time_of_last_minor_merger ), h5::datatype::native_float, "time_of_last_minor_merger" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, outflow_rate ),        h5::datatype::native_float, "outflow_rate" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, infall_mvir ),         h5::datatype::native_float, "infall_mvir" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, infall_vvir ),         h5::datatype::native_float, "infall_vvir" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, infall_vmax ),         h5::datatype::native_float, "infall_vmax" );
      der.commit( mem_type, file_type );
   }

   // Get subset of galaxy fields for proof of concept
   // Minimal set required for full pipeline: tree traversal + spatial indexing + validation
   std::vector<FieldInfo> get_galaxy_field_list_subset() {
      using namespace hpc::h5;

      std::vector<FieldInfo> fields;

      // Snapshot number (required for Phase 3 snapshot redistribution)
      fields.push_back({"snapnum", datatype::native_int, offsetof(galaxy, snapshot)});

      // Tree structure (required for Phase 2 traversal metadata)
      fields.push_back({"descendant", datatype::native_int, offsetof(galaxy, descendant)});
      fields.push_back({"local_index", datatype::native_int, offsetof(galaxy, local_index)});

      // Spatial coordinates (required for KD-tree in Phase 4)
      fields.push_back({"posx", datatype::native_float, offsetof(galaxy, pos[0])});
      fields.push_back({"posy", datatype::native_float, offsetof(galaxy, pos[1])});
      fields.push_back({"posz", datatype::native_float, offsetof(galaxy, pos[2])});

      // Stellar mass (useful for validation)
      fields.push_back({"stellar_mass", datatype::native_float, offsetof(galaxy, stellar_mass)});

      return fields;
   }

   // Get full list of galaxy fields for columnar storage
   std::vector<FieldInfo> get_galaxy_field_list_full() {
      using namespace hpc::h5;

      std::vector<FieldInfo> fields;

      // SAGE fields use CamelCase (matching SAGE HDF5 output)
      // Computed/pipeline fields use lowercase

      // Basic identifiers - SAGE CamelCase
      fields.push_back({"SnapNum", datatype::native_int, offsetof(galaxy, snapshot)});
      fields.push_back({"Type", datatype::native_int, offsetof(galaxy, type)});
      fields.push_back({"GalaxyIndex", datatype::native_llong, offsetof(galaxy, galaxy_idx)});
      fields.push_back({"CentralGalaxyIndex", datatype::native_llong, offsetof(galaxy, central_galaxy_idx)});
      fields.push_back({"SAGEHaloIndex", datatype::native_int, offsetof(galaxy, sage_halo_idx)});
      fields.push_back({"SAGETreeIndex", datatype::native_int, offsetof(galaxy, sage_tree_idx)});
      fields.push_back({"SimulationHaloIndex", datatype::native_llong, offsetof(galaxy, simulation_halo_idx)});

      // Traversal metadata - computed fields (lowercase)
      fields.push_back({"local_index", datatype::native_int, offsetof(galaxy, local_index)});
      fields.push_back({"global_index", datatype::native_llong, offsetof(galaxy, global_index)});
      fields.push_back({"descendant", datatype::native_int, offsetof(galaxy, descendant)});
      fields.push_back({"global_descendant", datatype::native_llong, offsetof(galaxy, global_descendant)});
      fields.push_back({"subsize", datatype::native_int, offsetof(galaxy, subsize)});

      // Merger information - SAGE mixed case
      fields.push_back({"mergeType", datatype::native_int, offsetof(galaxy, merge_type)});
      fields.push_back({"mergeIntoID", datatype::native_int, offsetof(galaxy, merge_into_id)});
      fields.push_back({"mergeIntoSnapNum", datatype::native_int, offsetof(galaxy, merge_into_snapshot)});
      fields.push_back({"dT", datatype::native_float, offsetof(galaxy, dt)});

      // Spatial coordinates - SAGE CamelCase
      fields.push_back({"Posx", datatype::native_float, offsetof(galaxy, pos[0])});
      fields.push_back({"Posy", datatype::native_float, offsetof(galaxy, pos[1])});
      fields.push_back({"Posz", datatype::native_float, offsetof(galaxy, pos[2])});

      // Velocities - SAGE CamelCase
      fields.push_back({"Velx", datatype::native_float, offsetof(galaxy, vel[0])});
      fields.push_back({"Vely", datatype::native_float, offsetof(galaxy, vel[1])});
      fields.push_back({"Velz", datatype::native_float, offsetof(galaxy, vel[2])});

      // Spin - SAGE CamelCase
      fields.push_back({"Spinx", datatype::native_float, offsetof(galaxy, spin[0])});
      fields.push_back({"Spiny", datatype::native_float, offsetof(galaxy, spin[1])});
      fields.push_back({"Spinz", datatype::native_float, offsetof(galaxy, spin[2])});

      // Halo properties - SAGE CamelCase
      fields.push_back({"Len", datatype::native_int, offsetof(galaxy, num_particles)});
      fields.push_back({"Mvir", datatype::native_float, offsetof(galaxy, mvir)});
      fields.push_back({"CentralMvir", datatype::native_float, offsetof(galaxy, central_mvir)});
      fields.push_back({"Rvir", datatype::native_float, offsetof(galaxy, rvir)});
      fields.push_back({"Vvir", datatype::native_float, offsetof(galaxy, vvir)});
      fields.push_back({"Vmax", datatype::native_float, offsetof(galaxy, vmax)});
      fields.push_back({"VelDisp", datatype::native_float, offsetof(galaxy, vel_disp)});

      // Baryonic reservoirs - SAGE CamelCase
      fields.push_back({"ColdGas", datatype::native_float, offsetof(galaxy, cold_gas)});
      fields.push_back({"StellarMass", datatype::native_float, offsetof(galaxy, stellar_mass)});
      fields.push_back({"BulgeMass", datatype::native_float, offsetof(galaxy, bulge_mass)});
      fields.push_back({"HotGas", datatype::native_float, offsetof(galaxy, hot_gas)});
      fields.push_back({"EjectedMass", datatype::native_float, offsetof(galaxy, ejected_mass)});
      fields.push_back({"BlackHoleMass", datatype::native_float, offsetof(galaxy, blackhole_mass)});
      fields.push_back({"IntraClusterStars", datatype::native_float, offsetof(galaxy, ics)});

      // Metals - SAGE CamelCase
      fields.push_back({"MetalsColdGas", datatype::native_float, offsetof(galaxy, metals_cold_gas)});
      fields.push_back({"MetalsStellarMass", datatype::native_float, offsetof(galaxy, metals_stellar_mass)});
      fields.push_back({"MetalsBulgeMass", datatype::native_float, offsetof(galaxy, metals_bulge_mass)});
      fields.push_back({"MetalsHotGas", datatype::native_float, offsetof(galaxy, metals_hot_gas)});
      fields.push_back({"MetalsEjectedMass", datatype::native_float, offsetof(galaxy, metals_ejected_mass)});
      fields.push_back({"MetalsIntraClusterStars", datatype::native_float, offsetof(galaxy, metals_ics)});

      // Star formation rates - SAGE CamelCase
      fields.push_back({"SfrDisk", datatype::native_float, offsetof(galaxy, sfr_disk)});
      fields.push_back({"SfrBulge", datatype::native_float, offsetof(galaxy, sfr_bulge)});
      fields.push_back({"SfrDiskZ", datatype::native_float, offsetof(galaxy, sfr_disk_z)});
      fields.push_back({"SfrBulgeZ", datatype::native_float, offsetof(galaxy, sfr_bulge_z)});

      // Miscellaneous - SAGE CamelCase
      fields.push_back({"DiskRadius", datatype::native_float, offsetof(galaxy, disk_scale_radius)});
      fields.push_back({"Cooling", datatype::native_float, offsetof(galaxy, cooling)});
      fields.push_back({"Heating", datatype::native_float, offsetof(galaxy, heating)});
      fields.push_back({"QuasarModeBHaccretionMass", datatype::native_float, offsetof(galaxy, quasar_mode_bh_accretion_mass)});
      fields.push_back({"TimeOfLastMajorMerger", datatype::native_float, offsetof(galaxy, time_of_last_major_merger)});
      fields.push_back({"TimeOfLastMinorMerger", datatype::native_float, offsetof(galaxy, time_of_last_minor_merger)});
      fields.push_back({"OutflowRate", datatype::native_float, offsetof(galaxy, outflow_rate)});

      // Infall properties - SAGE mixed case
      fields.push_back({"infallMvir", datatype::native_float, offsetof(galaxy, infall_mvir)});
      fields.push_back({"infallVvir", datatype::native_float, offsetof(galaxy, infall_vvir)});
      fields.push_back({"infallVmax", datatype::native_float, offsetof(galaxy, infall_vmax)});

      return fields;
   }

   void addField(std::ofstream &outfile, h5::datatype &datatype, std::string name, int32_t order) {
    std::string type;
    
    if (datatype == h5::datatype::native_int) {
        type = "int";
    } else if (datatype == h5::datatype::native_float) {
        type = "float";
    } else if (datatype == h5::datatype::native_llong) {
        type = "long long";
    } else {
        type = "Unknown";
    }

    outfile << "    <Field Type=\"" << type << "\"" << std::endl
            << "     label=\"" << name << "\"" << std::endl
            << "     description=\"" << name << "\"" << std::endl
            << "     units=\"" << "" << "\"" << std::endl
            << "     order=\"" << order << "\"" << std::endl
            << "     group=\"" << "general" << "\">" << name << "</Field>" << std::endl;
    }

    void writeSageInputFields_tree(std::ofstream& outfile) {
        outfile << "  <sageinput>" << std::endl;

        int32_t order = 1;
        addField( outfile, h5::datatype::native_int, "objecttype", order );
        addField( outfile, h5::datatype::native_float, "posx", order );
        addField( outfile, h5::datatype::native_float, "posy", order );
        addField( outfile, h5::datatype::native_float, "posz", order );
        addField( outfile, h5::datatype::native_float, "velx", order );
        addField( outfile, h5::datatype::native_float, "vely", order );
        addField( outfile, h5::datatype::native_float, "velz", order );
        addField( outfile, h5::datatype::native_int, "snapnum", order );
        addField( outfile, h5::datatype::native_llong, "globalindex", order );
        addField( outfile, h5::datatype::native_int, "descendant", order );
        addField( outfile, h5::datatype::native_int, "mergetype", order );
        addField( outfile, h5::datatype::native_float, "dt", order );
        addField( outfile, h5::datatype::native_float, "sfrdisk", order );
        addField( outfile, h5::datatype::native_float, "sfrbulge", order );
        addField( outfile, h5::datatype::native_float, "sfrdiskz", order );
        addField( outfile, h5::datatype::native_float, "sfrbulgez", order );
        addField( outfile, h5::datatype::native_int, "treeindex", order );
        addField( outfile, h5::datatype::native_int, "localindex", order );
        addField( outfile, h5::datatype::native_llong, "globaldescendant", order );
        addField( outfile, h5::datatype::native_int, "subsize", order );
        addField( outfile, h5::datatype::native_float, "coldgas", order );
        addField( outfile, h5::datatype::native_float, "metalscoldgas", order );
        addField( outfile, h5::datatype::native_float, "diskscaleradius", order );
        addField( outfile, h5::datatype::native_llong, "galaxyindex", order );
        addField( outfile, h5::datatype::native_llong, "centralgalaxyindex", order );
        addField( outfile, h5::datatype::native_llong, "simulationhaloindex", order );
        addField( outfile, h5::datatype::native_int, "mergeintoid", order );
        addField( outfile, h5::datatype::native_int, "mergeintosnapnum", order );
        addField( outfile, h5::datatype::native_float, "spin_x", order );
        addField( outfile, h5::datatype::native_float, "spin_y", order );
        addField( outfile, h5::datatype::native_float, "spin_z", order );
        addField( outfile, h5::datatype::native_int, "len", order );
        addField( outfile, h5::datatype::native_float, "mvir", order );
        addField( outfile, h5::datatype::native_float, "centralmvir", order );
        addField( outfile, h5::datatype::native_float, "rvir", order );
        addField( outfile, h5::datatype::native_float, "vvir", order );
        addField( outfile, h5::datatype::native_float, "vmax", order );
        addField( outfile, h5::datatype::native_float, "veldisp", order );
        addField( outfile, h5::datatype::native_float, "stellarmass", order );
        addField( outfile, h5::datatype::native_float, "bulgemass", order );
        addField( outfile, h5::datatype::native_float, "hotgas", order );
        addField( outfile, h5::datatype::native_float, "ejectedmass", order );
        addField( outfile, h5::datatype::native_float, "blackholemass", order );
        addField( outfile, h5::datatype::native_float, "ics", order );
        addField( outfile, h5::datatype::native_float, "metalsstellarmass", order );
        addField( outfile, h5::datatype::native_float, "metalsbulgemass", order );
        addField( outfile, h5::datatype::native_float, "metalshotgas", order );
        addField( outfile, h5::datatype::native_float, "metalsejectedmass", order );
        addField( outfile, h5::datatype::native_float, "metalsics", order );
        addField( outfile, h5::datatype::native_float, "cooling", order );
        addField( outfile, h5::datatype::native_float, "heating", order );
        addField( outfile, h5::datatype::native_float, "quasarmodebhaccretionmass", order );
        addField( outfile, h5::datatype::native_float, "timeoflastmajormerger", order );
        addField( outfile, h5::datatype::native_float, "timeoflastminormerger", order );
        addField( outfile, h5::datatype::native_float, "outflowrate", order );
        addField( outfile, h5::datatype::native_float, "infallmvir", order );
        addField( outfile, h5::datatype::native_float, "infallvvir", order );
        addField( outfile, h5::datatype::native_float, "infallvmax", order );
        addField( outfile, h5::datatype::native_float, "totsfr", order );
        addField( outfile, h5::datatype::native_float, "vpeak", order );
        addField( outfile, h5::datatype::native_int,   "isflyby", order );

        outfile << "  </sageinput>" << std::endl;
    }   
    void writeSageInputFields(std::ofstream& outfile) {
        outfile << "  <sageinput>" << std::endl;

        int32_t order = 1;
        addField( outfile, h5::datatype::native_int,   "snapnum" , order++);
        addField( outfile, h5::datatype::native_int,   "type" , order++);
        addField( outfile, h5::datatype::native_llong, "galaxy_index" , order++);
        addField( outfile, h5::datatype::native_llong, "central_galaxy_index" , order++);
        addField( outfile, h5::datatype::native_int,   "sage_halo_index" , order++);
        addField( outfile, h5::datatype::native_int,   "sage_tree_index" , order++);
        addField( outfile, h5::datatype::native_llong, "simulation_halo_index" , order++);
        addField( outfile, h5::datatype::native_int,   "local_index" , order++);
        addField( outfile, h5::datatype::native_llong, "global_index" , order++);
        addField( outfile, h5::datatype::native_int,   "descendant" , order++);
        addField( outfile, h5::datatype::native_llong, "global_descendant" , order++);
        addField( outfile, h5::datatype::native_int,   "subsize" , order++);
        addField( outfile, h5::datatype::native_int,   "merge_type" , order++);
        addField( outfile, h5::datatype::native_int,   "merge_into_id" , order++);
        addField( outfile, h5::datatype::native_int,   "merge_into_snapshot" , order++);
        addField( outfile, h5::datatype::native_float, "dt" , order++);
        addField( outfile, h5::datatype::native_float, "posx" , order++);
        addField( outfile, h5::datatype::native_float, "posy" , order++);
        addField( outfile, h5::datatype::native_float, "posz" , order++);
        addField( outfile, h5::datatype::native_float, "velx" , order++);
        addField( outfile, h5::datatype::native_float, "vely" , order++);
        addField( outfile, h5::datatype::native_float, "velz" , order++);
        addField( outfile, h5::datatype::native_float, "spin_x" , order++);
        addField( outfile, h5::datatype::native_float, "spin_y" , order++);
        addField( outfile, h5::datatype::native_float, "spin_z" , order++);
        addField( outfile, h5::datatype::native_int,   "n_darkmatter_particles" , order++);
        addField( outfile, h5::datatype::native_float, "virial_mass" , order++);
        addField( outfile, h5::datatype::native_float, "central_galaxy_mvir" , order++);
        addField( outfile, h5::datatype::native_float, "virial_radius" , order++);
        addField( outfile, h5::datatype::native_float, "virial_velocity" , order++);
        addField( outfile, h5::datatype::native_float, "max_velocity" , order++);
        addField( outfile, h5::datatype::native_float, "velocity_dispersion" , order++);
        addField( outfile, h5::datatype::native_float, "cold_gas" , order++);
        addField( outfile, h5::datatype::native_float, "stellar_mass" , order++);
        addField( outfile, h5::datatype::native_float, "bulge_mass" , order++);
        addField( outfile, h5::datatype::native_float, "hot_gas" , order++);
        addField( outfile, h5::datatype::native_float, "ejected_mass" , order++);
        addField( outfile, h5::datatype::native_float, "blackhole_mass" , order++);
        addField( outfile, h5::datatype::native_float, "ics" , order++);
        addField( outfile, h5::datatype::native_float, "metals_cold_gas" , order++);
        addField( outfile, h5::datatype::native_float, "metals_stellar_mass" , order++);
        addField( outfile, h5::datatype::native_float, "metals_bulge_mass" , order++);
        addField( outfile, h5::datatype::native_float, "metals_hot_gas" , order++);
        addField( outfile, h5::datatype::native_float, "metals_ejected_mass" , order++);
        addField( outfile, h5::datatype::native_float, "metals_ics" , order++);
        addField( outfile, h5::datatype::native_float, "sfr_disk" , order++);
        addField( outfile, h5::datatype::native_float, "sfr_bulge" , order++);
        addField( outfile, h5::datatype::native_float, "sfr_disk_z" , order++);
        addField( outfile, h5::datatype::native_float, "sfr_bulge_z" , order++);
        addField( outfile, h5::datatype::native_float, "disk_scale_radius" , order++);
        addField( outfile, h5::datatype::native_float, "cooling" , order++);
        addField( outfile, h5::datatype::native_float, "heating" , order++);
        addField( outfile, h5::datatype::native_float, "quasar_mode_bh_accretion_mass" , order++);
        addField( outfile, h5::datatype::native_float, "time_of_last_major_merger" , order++);
        addField( outfile, h5::datatype::native_float, "time_of_last_minor_merger" , order++);
        addField( outfile, h5::datatype::native_float, "outflow_rate" , order++);
        addField( outfile, h5::datatype::native_float, "infall_mvir" , order++);
        addField( outfile, h5::datatype::native_float, "infall_vvir" , order++);
        addField( outfile, h5::datatype::native_float, "infall_vmax" , order++);
        outfile << "  </sageinput>" << std::endl;
    }


  void
  make_hdf5_sidecar( std::string simulation_name, std::string inname, const std::set<double>& redshifts, const double& hubble, const double& box_size )
  {
    std::string xml_filename = inname.substr(0, inname.find_last_of('.')) + ".xml";
    std::ofstream outfile(xml_filename.c_str());
    if (!outfile.is_open()) {
      std::cerr << "Error: Could not open file " << xml_filename << std::endl;
      return;
    }
    std::cout << "Creating sidecar XML file: " << xml_filename << std::endl;

    

    writeSageInputFields(outfile);
    outfile.close();

    std::string xml_settings_filename = inname.substr(0, inname.find_last_of('.')) + "_import_settings.xml";
    std::ofstream outfile_settings(xml_settings_filename.c_str());
    if (!outfile_settings.is_open()) {
      std::cerr << "Error: Could not open file " << xml_settings_filename << std::endl;
      return;
    }
    std::cout << "Creating settings XML file: " << xml_settings_filename << std::endl;

    outfile_settings << "<settings>" << std::endl;
    writeSageInputFields(outfile_settings);
    outfile_settings << "  <PGDB>" << std::endl;
    outfile_settings << "    <ScienceModuleDBUserName>" << "taoadmin" << "</ScienceModuleDBUserName>" << std::endl;
    outfile_settings << "    <TreeTablePrefix>" << "tree_" << "</TreeTablePrefix>" << std::endl;
    outfile_settings << "    <NewDBName>" << "tmpray1" << "</NewDBName>" << std::endl;
    outfile_settings << "    <NewDBAlias>" << "tmpray1" << "</NewDBAlias>" << std::endl;
    outfile_settings << "    <ServersCount>" << "1" << "</ServersCount>" << std::endl;
    outfile_settings << "    <snapshots>" << std::endl;
    for (double z : redshifts) {
        outfile_settings << "      <snapshot>" << std::fixed << std::setprecision(12-std::to_string((int)z).length()) << z << "</snapshot>" << std::endl;
    }
    outfile_settings << "    </snapshots>" << std::endl;
    outfile_settings << "  </PGDB>" << std::endl;
    outfile_settings << "  <RunningSettings>" << std::endl;
    outfile_settings << "    <Simulation>" << simulation_name << "</Simulation>" << std::endl;
    outfile_settings << "    <GalaxyModel>" << "SAGE" << "</GalaxyModel>" << std::endl;
    outfile_settings << "    <hubble>" << hubble << "</hubble>" << std::endl;
    outfile_settings << "    <InputFile>" << inname << "</InputFile>" << std::endl;
    // Split inname into directory and filename
    size_t last_slash = inname.find_last_of("/\\");
    std::string directory = "";
    std::string filename = inname;
    
    if (last_slash != std::string::npos) {
      directory = inname.substr(0, last_slash + 1);
      filename = inname.substr(last_slash + 1);
    }
    
    std::string output_file = directory + "out_" + filename;
    outfile_settings << "    <OutputFile>" << output_file << "</OutputFile>" << std::endl;
    outfile_settings << "    <SimulationBoxX>" << box_size << "</SimulationBoxX>" << std::endl;
    outfile_settings << "    <SimulationBoxY>" << box_size << "</SimulationBoxY>" << std::endl;
    outfile_settings << "    <BSPCellSize>" << "10" << "</BSPCellSize>" << std::endl;
    outfile_settings << "    <GalaxiesPerTable>" << "50000" << "</GalaxiesPerTable>" << std::endl;
    outfile_settings << "    <OverWriteDB>" << "yes" << "</OverWriteDB>" << std::endl;
    outfile_settings << "    <RegenerateTables>" << "yes" << "</RegenerateTables>" << std::endl;
    outfile_settings << "  </RunningSettings>" << std::endl;
    outfile_settings << "  <TreeTraversal>" << std::endl;
    outfile_settings << "    <item>" << "global_index" << "</item>" << std::endl;
    outfile_settings << "    <item>" << "descendant" << "</item>" << std::endl;
    outfile_settings << "    <item>" << "snapnum" << "</item>" << std::endl;
    outfile_settings << "  </TreeTraversal>" << std::endl;
    outfile_settings << "</settings>" << std::endl;
    outfile_settings.close();
  }

   void
  make_hdf5_sidecar_tree( std::string simulation_name, std::string inname, const std::set<double>& redshifts, const double& hubble, const double& box_size )
  {
    std::string xml_filename = inname.substr(0, inname.find_last_of('.')) + ".xml";
    std::ofstream outfile(xml_filename.c_str());
    if (!outfile.is_open()) {
      std::cerr << "Error: Could not open file " << xml_filename << std::endl;
      return;
    }
    std::cout << "Creating sidecar XML file: " << xml_filename << std::endl;

    

    writeSageInputFields_tree(outfile);
    outfile.close();

    std::string xml_settings_filename = inname.substr(0, inname.find_last_of('.')) + "_import_settings.xml";
    std::ofstream outfile_settings(xml_settings_filename.c_str());
    if (!outfile_settings.is_open()) {
      std::cerr << "Error: Could not open file " << xml_settings_filename << std::endl;
      return;
    }
    std::cout << "Creating settings XML file: " << xml_settings_filename << std::endl;

    outfile_settings << "<settings>" << std::endl;
    writeSageInputFields(outfile_settings);
    outfile_settings << "  <PGDB>" << std::endl;
    outfile_settings << "    <ScienceModuleDBUserName>" << "taoadmin" << "</ScienceModuleDBUserName>" << std::endl;
    outfile_settings << "    <TreeTablePrefix>" << "tree_" << "</TreeTablePrefix>" << std::endl;
    outfile_settings << "    <NewDBName>" << "tmpray1" << "</NewDBName>" << std::endl;
    outfile_settings << "    <NewDBAlias>" << "tmpray1" << "</NewDBAlias>" << std::endl;
    outfile_settings << "    <ServersCount>" << "1" << "</ServersCount>" << std::endl;
    outfile_settings << "    <snapshots>" << std::endl;
    for (double z : redshifts) {
        outfile_settings << "      <snapshot>" << std::fixed << std::setprecision(12-std::to_string((int)z).length()) << z << "</snapshot>" << std::endl;
    }
    outfile_settings << "    </snapshots>" << std::endl;
    outfile_settings << "  </PGDB>" << std::endl;
    outfile_settings << "  <RunningSettings>" << std::endl;
    outfile_settings << "    <Simulation>" << simulation_name << "</Simulation>" << std::endl;
    outfile_settings << "    <GalaxyModel>" << "SAGE" << "</GalaxyModel>" << std::endl;
    outfile_settings << "    <hubble>" << hubble << "</hubble>" << std::endl;
    outfile_settings << "    <InputFile>" << inname << "</InputFile>" << std::endl;
    // Split inname into directory and filename
    size_t last_slash = inname.find_last_of("/\\");
    std::string directory = "";
    std::string filename = inname;
    
    if (last_slash != std::string::npos) {
      directory = inname.substr(0, last_slash + 1);
      filename = inname.substr(last_slash + 1);
    }
    
    std::string output_file = directory + "out_" + filename;
    outfile_settings << "    <OutputFile>" << output_file << "</OutputFile>" << std::endl;
    outfile_settings << "    <SimulationBoxX>" << box_size << "</SimulationBoxX>" << std::endl;
    outfile_settings << "    <SimulationBoxY>" << box_size << "</SimulationBoxY>" << std::endl;
    outfile_settings << "    <BSPCellSize>" << "10" << "</BSPCellSize>" << std::endl;
    outfile_settings << "    <GalaxiesPerTable>" << "50000" << "</GalaxiesPerTable>" << std::endl;
    outfile_settings << "    <OverWriteDB>" << "yes" << "</OverWriteDB>" << std::endl;
    outfile_settings << "    <RegenerateTables>" << "yes" << "</RegenerateTables>" << std::endl;
    outfile_settings << "  </RunningSettings>" << std::endl;
    outfile_settings << "  <TreeTraversal>" << std::endl;
    outfile_settings << "    <item>" << "global_index" << "</item>" << std::endl;
    outfile_settings << "    <item>" << "descendant" << "</item>" << std::endl;
    outfile_settings << "    <item>" << "snapshot" << "</item>" << std::endl;
    outfile_settings << "  </TreeTraversal>" << std::endl;
    outfile_settings << "</settings>" << std::endl;
    outfile_settings.close();
  }

}
