#include <iostream>
#include <libhpc/h5/derive.hh>
#include "sage.hh"

namespace sage {

   using namespace hpc;

   void
   make_hdf5_types( h5::datatype& mem_type,
		    h5::datatype& file_type )
   {
      h5::derive der( sizeof(galaxy) );
      der.add( h5::datatype::native_int,   HOFFSET( galaxy, snapshot ),            h5::datatype::std_i32be,  "snapnum" );
      der.add( h5::datatype::native_int,   HOFFSET( galaxy, type ),                h5::datatype::std_i32be,  "type" );
      der.add( h5::datatype::native_llong, HOFFSET( galaxy, galaxy_idx ),          h5::datatype::std_i64be,  "galaxy_index" );
      der.add( h5::datatype::native_llong, HOFFSET( galaxy, central_galaxy_idx ),  h5::datatype::std_i64be,  "central_galaxy_index" );
      der.add( h5::datatype::native_int,   HOFFSET( galaxy, sage_halo_idx ),       h5::datatype::std_i32be,  "sage_halo_index" );
      der.add( h5::datatype::native_int,   HOFFSET( galaxy, sage_tree_idx ),       h5::datatype::std_i32be,  "sage_tree_index" );
      der.add( h5::datatype::native_llong, HOFFSET( galaxy, simulation_halo_idx ), h5::datatype::std_i32be,  "simulation_halo_index" );
      der.add( h5::datatype::native_int,   HOFFSET( galaxy, local_index ),         h5::datatype::std_i32be,  "local_index" );
      der.add( h5::datatype::native_llong, HOFFSET( galaxy, global_index ),        h5::datatype::std_i64be,  "global_index" );
      der.add( h5::datatype::native_int,   HOFFSET( galaxy, descendant ),          h5::datatype::std_i32be,  "descendant" );
      der.add( h5::datatype::native_llong, HOFFSET( galaxy, global_descendant ),   h5::datatype::std_i64be,  "global_descendant" );
      der.add( h5::datatype::native_int,   HOFFSET( galaxy, subsize ),             h5::datatype::std_i32be,  "subsize" );
      der.add( h5::datatype::native_int,   HOFFSET( galaxy, merge_type ),          h5::datatype::std_i32be,  "merge_type" );
      der.add( h5::datatype::native_int,   HOFFSET( galaxy, merge_into_id ),       h5::datatype::std_i32be,  "merge_into_id" );
      der.add( h5::datatype::native_int,   HOFFSET( galaxy, merge_into_snapshot ), h5::datatype::std_i32be,  "merge_into_snapshot" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, dt ),                  h5::datatype::ieee_f32be, "dt" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, pos[0] ),              h5::datatype::ieee_f32be, "posx" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, pos[1] ),              h5::datatype::ieee_f32be, "posy" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, pos[2] ),              h5::datatype::ieee_f32be, "posz" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, vel[0] ),              h5::datatype::ieee_f32be, "velx" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, vel[1] ),              h5::datatype::ieee_f32be, "vely" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, vel[2] ),              h5::datatype::ieee_f32be, "velz" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, spin[0] ),             h5::datatype::ieee_f32be, "spin_x" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, spin[1] ),             h5::datatype::ieee_f32be, "spin_y" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, spin[2] ),             h5::datatype::ieee_f32be, "spin_z" );
      der.add( h5::datatype::native_int,   HOFFSET( galaxy, num_particles ),       h5::datatype::std_i32be,  "n_darkmatter_particles" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, mvir ),                h5::datatype::ieee_f32be, "virial_mass" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, central_mvir ),        h5::datatype::ieee_f32be, "central_galaxy_mvir" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, rvir ),                h5::datatype::ieee_f32be, "virial_radius" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, vvir ),                h5::datatype::ieee_f32be, "virial_velocity" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, vmax ),                h5::datatype::ieee_f32be, "max_velocity" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, vel_disp ),            h5::datatype::ieee_f32be, "velocity_dispersion" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, cold_gas ),            h5::datatype::ieee_f32be, "cold_gas" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, stellar_mass ),        h5::datatype::ieee_f32be, "stellar_mass" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, bulge_mass ),          h5::datatype::ieee_f32be, "bulge_mass" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, hot_gas ),             h5::datatype::ieee_f32be, "hot_gas" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, ejected_mass ),        h5::datatype::ieee_f32be, "ejected_mass" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, blackhole_mass ),      h5::datatype::ieee_f32be, "blackhole_mass" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, ics ),                 h5::datatype::ieee_f32be, "ics" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, metals_cold_gas ),     h5::datatype::ieee_f32be, "metals_cold_gas" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, metals_stellar_mass ), h5::datatype::ieee_f32be, "metals_stellar_mass" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, metals_bulge_mass ),   h5::datatype::ieee_f32be, "metals_bulge_mass" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, metals_hot_gas ),      h5::datatype::ieee_f32be, "metals_hot_gas" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, metals_ejected_mass ), h5::datatype::ieee_f32be, "metals_ejected_mass" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, metals_ics ),          h5::datatype::ieee_f32be, "metals_ics" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, sfr_disk ),            h5::datatype::ieee_f32be, "sfr_disk" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, sfr_bulge ),           h5::datatype::ieee_f32be, "sfr_bulge" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, sfr_disk_z ),          h5::datatype::ieee_f32be, "sfr_disk_z" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, sfr_bulge_z ),         h5::datatype::ieee_f32be, "sfr_bulge_z" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, disk_scale_radius ),   h5::datatype::ieee_f32be, "disk_scale_radius" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, cooling ),             h5::datatype::ieee_f32be, "cooling" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, heating ),             h5::datatype::ieee_f32be, "heating" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, quasar_mode_bh_accretion_mass ), h5::datatype::ieee_f32be, "quasar_mode_bh_accretion_mass" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, time_of_last_major_merger ), h5::datatype::ieee_f32be, "time_of_last_major_merger" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, time_of_last_minor_merger ), h5::datatype::ieee_f32be, "time_of_last_minor_merger" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, outflow_rate ),        h5::datatype::ieee_f32be, "outflow_rate" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, infall_mvir ),         h5::datatype::ieee_f32be, "infall_mvir" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, infall_vvir ),         h5::datatype::ieee_f32be, "infall_vvir" );
      der.add( h5::datatype::native_float, HOFFSET( galaxy, infall_vmax ),         h5::datatype::ieee_f32be, "infall_vmax" );
      der.commit( mem_type, file_type );
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
    outfile_settings << "    <OutputFile>" << "out_" << inname << "</OutputFile>" << std::endl;
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

   void
   make_hdf5_types_tree( h5::datatype& mem_type,
                         h5::datatype& file_type )
   {
       h5::derive der( sizeof(galaxyTree) );
       der.add(h5::datatype::native_int,HOFFSET(galaxyTree,objecttype),h5::datatype::std_i32be,"objecttype");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,pos[0]),h5::datatype::ieee_f32be,"posx");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,pos[1]),h5::datatype::ieee_f32be,"posy");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,pos[2]),h5::datatype::ieee_f32be,"posz");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,velx),h5::datatype::ieee_f32be,"velx");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,vely),h5::datatype::ieee_f32be,"vely");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,velz),h5::datatype::ieee_f32be,"velz");
       der.add(h5::datatype::native_int,HOFFSET(galaxyTree,snapnum),h5::datatype::std_i32be,"snapnum");
       der.add(h5::datatype::native_llong,HOFFSET(galaxyTree,globalindex),h5::datatype::std_i64be,"globalindex");
       der.add(h5::datatype::native_int,HOFFSET(galaxyTree,descendant),h5::datatype::std_i32be,"descendant");
       der.add(h5::datatype::native_int,HOFFSET(galaxyTree,mergetype),h5::datatype::std_i32be,"mergetype");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,dt),h5::datatype::ieee_f32be,"dt");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,sfrdisk),h5::datatype::ieee_f32be,"sfrdisk");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,sfrbulge),h5::datatype::ieee_f32be,"sfrbulge");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,sfrdiskz),h5::datatype::ieee_f32be,"sfrdiskz");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,sfrbulgez),h5::datatype::ieee_f32be,"sfrbulgez");
       der.add(h5::datatype::native_int,HOFFSET(galaxyTree,treeindex),h5::datatype::std_i32be,"treeindex");
       der.add(h5::datatype::native_int,HOFFSET(galaxyTree,localindex),h5::datatype::std_i32be,"localindex");
       der.add(h5::datatype::native_llong,HOFFSET(galaxyTree,globaldescendant),h5::datatype::std_i64be,"globaldescendant");
       der.add(h5::datatype::native_int,HOFFSET(galaxyTree,subsize),h5::datatype::std_i32be,"subsize");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,coldgas),h5::datatype::ieee_f32be,"coldgas");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,metalscoldgas),h5::datatype::ieee_f32be,"metalscoldgas");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,diskscaleradius),h5::datatype::ieee_f32be,"diskscaleradius");
       der.add(h5::datatype::native_llong,HOFFSET(galaxyTree,galaxyindex),h5::datatype::std_i64be,"GalaxyIndex");
       der.add(h5::datatype::native_llong,HOFFSET(galaxyTree,centralgalaxyindex),h5::datatype::std_i64be,"CentralGalaxyIndex");
       der.add(h5::datatype::native_llong,HOFFSET(galaxyTree,simulationhaloindex),h5::datatype::std_i64be,"SimulationHaloIndex");
       der.add(h5::datatype::native_int,HOFFSET(galaxyTree,mergeintoid),h5::datatype::std_i32be,"mergeIntoID");
       der.add(h5::datatype::native_int,HOFFSET(galaxyTree,mergeintosnapnum),h5::datatype::std_i32be,"mergeIntoSnapNum");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,spin_x),h5::datatype::ieee_f32be,"Spin_x");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,spin_y),h5::datatype::ieee_f32be,"Spin_y");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,spin_z),h5::datatype::ieee_f32be,"Spin_z");
       der.add(h5::datatype::native_int,HOFFSET(galaxyTree,len),h5::datatype::std_i32be,"Len");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,mvir),h5::datatype::ieee_f32be,"Mvir");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,centralmvir),h5::datatype::ieee_f32be,"CentralMvir");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,rvir),h5::datatype::ieee_f32be,"Rvir");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,vvir),h5::datatype::ieee_f32be,"Vvir");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,vmax),h5::datatype::ieee_f32be,"Vmax");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,veldisp),h5::datatype::ieee_f32be,"VelDisp");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,stellarmass),h5::datatype::ieee_f32be,"StellarMass");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,bulgemass),h5::datatype::ieee_f32be,"BulgeMass");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,hotgas),h5::datatype::ieee_f32be,"HotGas");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,ejectedmass),h5::datatype::ieee_f32be,"EjectedMass");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,blackholemass),h5::datatype::ieee_f32be,"BlackHoleMass");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,ics),h5::datatype::ieee_f32be,"ICS");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,metalsstellarmass),h5::datatype::ieee_f32be,"MetalsStellarMass");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,metalsbulgemass),h5::datatype::ieee_f32be,"MetalsBulgeMass");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,metalshotgas),h5::datatype::ieee_f32be,"MetalsHotGas");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,metalsejectedmass),h5::datatype::ieee_f32be,"MetalsEjectedMass");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,metalsics),h5::datatype::ieee_f32be,"MetalsICS");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,cooling),h5::datatype::ieee_f32be,"Cooling");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,heating),h5::datatype::ieee_f32be,"Heating");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,quasarmodebhaccretionmass),h5::datatype::ieee_f32be,"QuasarModeBHaccretionMass");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,timeoflastmajormerger),h5::datatype::ieee_f32be,"TimeofLastMajorMerger");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,timeoflastminormerger),h5::datatype::ieee_f32be,"TimeofLastMinorMerger");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,outflowrate),h5::datatype::ieee_f32be,"OutflowRate");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,infallmvir),h5::datatype::ieee_f32be,"infallMvir");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,infallvvir),h5::datatype::ieee_f32be,"infallVvir");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,infallvmax),h5::datatype::ieee_f32be,"infallVmax");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,totsfr),h5::datatype::ieee_f32be,"TotSfr");
       der.add(h5::datatype::native_float,HOFFSET(galaxyTree,vpeak),h5::datatype::ieee_f32be,"Vpeak");
       der.add(h5::datatype::native_int,HOFFSET(galaxyTree,isflyby),h5::datatype::std_i32be,"isFlyby");
       der.commit( mem_type, file_type );

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
    outfile_settings << "    <OutputFile>" << "out_" << inname << "</OutputFile>" << std::endl;
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

    void
    make_hdf5_types_treesed( h5::datatype& mem_type,
                          h5::datatype& file_type )
    {
        h5::derive der( sizeof(galaxyTreeSed) );
        der.add(h5::datatype::native_int,HOFFSET(galaxyTreeSed,objecttype),h5::datatype::std_i32be,"objecttype");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,pos[0]),h5::datatype::ieee_f32be,"posx");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,pos[1]),h5::datatype::ieee_f32be,"posy");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,pos[2]),h5::datatype::ieee_f32be,"posz");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,velx),h5::datatype::ieee_f32be,"velx");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,vely),h5::datatype::ieee_f32be,"vely");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,velz),h5::datatype::ieee_f32be,"velz");
        der.add(h5::datatype::native_int,HOFFSET(galaxyTreeSed,snapnum),h5::datatype::std_i32be,"snapnum");
        der.add(h5::datatype::native_llong,HOFFSET(galaxyTreeSed,globalindex),h5::datatype::std_i64be,"globalindex");
        der.add(h5::datatype::native_int,HOFFSET(galaxyTreeSed,descendant),h5::datatype::std_i32be,"descendant");
        der.add(h5::datatype::native_int,HOFFSET(galaxyTreeSed,mergetype),h5::datatype::std_i32be,"mergetype");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,dt),h5::datatype::ieee_f32be,"dt");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,sfrdisk),h5::datatype::ieee_f32be,"sfrdisk");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,sfrbulge),h5::datatype::ieee_f32be,"sfrbulge");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,sfrdiskz),h5::datatype::ieee_f32be,"sfrdiskz");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,sfrbulgez),h5::datatype::ieee_f32be,"sfrbulgez");
        der.add(h5::datatype::native_int,HOFFSET(galaxyTreeSed,treeindex),h5::datatype::std_i32be,"treeindex");
        der.add(h5::datatype::native_int,HOFFSET(galaxyTreeSed,localindex),h5::datatype::std_i32be,"localindex");
        der.add(h5::datatype::native_llong,HOFFSET(galaxyTreeSed,globaldescendant),h5::datatype::std_i64be,"globaldescendant");
        der.add(h5::datatype::native_int,HOFFSET(galaxyTreeSed,subsize),h5::datatype::std_i32be,"subsize");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,coldgas),h5::datatype::ieee_f32be,"coldgas");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,metalscoldgas),h5::datatype::ieee_f32be,"metalscoldgas");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,diskscaleradius),h5::datatype::ieee_f32be,"diskscaleradius");
        der.add(h5::datatype::native_llong,HOFFSET(galaxyTreeSed,galaxyindex),h5::datatype::std_i64be,"GalaxyIndex");
        der.add(h5::datatype::native_llong,HOFFSET(galaxyTreeSed,centralgalaxyindex),h5::datatype::std_i64be,"CentralGalaxyIndex");
        der.add(h5::datatype::native_llong,HOFFSET(galaxyTreeSed,simulationhaloindex),h5::datatype::std_i64be,"SimulationHaloIndex");
        der.add(h5::datatype::native_int,HOFFSET(galaxyTreeSed,mergeintoid),h5::datatype::std_i32be,"mergeIntoID");
        der.add(h5::datatype::native_int,HOFFSET(galaxyTreeSed,mergeintosnapnum),h5::datatype::std_i32be,"mergeIntoSnapNum");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,spin_x),h5::datatype::ieee_f32be,"Spin_x");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,spin_y),h5::datatype::ieee_f32be,"Spin_y");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,spin_z),h5::datatype::ieee_f32be,"Spin_z");
        der.add(h5::datatype::native_int,HOFFSET(galaxyTreeSed,len),h5::datatype::std_i32be,"Len");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,mvir),h5::datatype::ieee_f32be,"Mvir");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,centralmvir),h5::datatype::ieee_f32be,"CentralMvir");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,rvir),h5::datatype::ieee_f32be,"Rvir");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,vvir),h5::datatype::ieee_f32be,"Vvir");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,vmax),h5::datatype::ieee_f32be,"Vmax");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,veldisp),h5::datatype::ieee_f32be,"VelDisp");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,stellarmass),h5::datatype::ieee_f32be,"StellarMass");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,bulgemass),h5::datatype::ieee_f32be,"BulgeMass");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,hotgas),h5::datatype::ieee_f32be,"HotGas");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,ejectedmass),h5::datatype::ieee_f32be,"EjectedMass");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,blackholemass),h5::datatype::ieee_f32be,"BlackHoleMass");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,ics),h5::datatype::ieee_f32be,"ICS");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,metalsstellarmass),h5::datatype::ieee_f32be,"MetalsStellarMass");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,metalsbulgemass),h5::datatype::ieee_f32be,"MetalsBulgeMass");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,metalshotgas),h5::datatype::ieee_f32be,"MetalsHotGas");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,metalsejectedmass),h5::datatype::ieee_f32be,"MetalsEjectedMass");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,metalsics),h5::datatype::ieee_f32be,"MetalsICS");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,cooling),h5::datatype::ieee_f32be,"Cooling");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,heating),h5::datatype::ieee_f32be,"Heating");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,quasarmodebhaccretionmass),h5::datatype::ieee_f32be,"QuasarModeBHaccretionMass");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,timeoflastmajormerger),h5::datatype::ieee_f32be,"TimeofLastMajorMerger");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,timeoflastminormerger),h5::datatype::ieee_f32be,"TimeofLastMinorMerger");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,outflowrate),h5::datatype::ieee_f32be,"OutflowRate");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,infallmvir),h5::datatype::ieee_f32be,"infallMvir");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,infallvvir),h5::datatype::ieee_f32be,"infallVvir");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,infallvmax),h5::datatype::ieee_f32be,"infallVmax");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,totsfr),h5::datatype::ieee_f32be,"TotSfr");
        der.add(h5::datatype::native_float,HOFFSET(galaxyTreeSed,vpeak),h5::datatype::ieee_f32be,"Vpeak");
        der.add(h5::datatype::native_int,HOFFSET(galaxyTreeSed,isflyby),h5::datatype::std_i32be,"isFlyby");

        der.add(h5::datatype::native_llong,HOFFSET(galaxyTreeSed,globaltreeid),h5::datatype::std_i64be,"globaltreeid");
        der.add(h5::datatype::native_llong,HOFFSET(galaxyTreeSed,breadthfirst_traversalorder),h5::datatype::std_i64be,"breadthfirst_traversalorder");
        der.add(h5::datatype::native_llong,HOFFSET(galaxyTreeSed,depthfirst_traversalorder),h5::datatype::std_i64be,"depthfirst_traversalorder");
        der.add(h5::datatype::native_llong,HOFFSET(galaxyTreeSed,subtree_count),h5::datatype::std_i64be,"subtree_count");
        der.add(h5::datatype::native_int,HOFFSET(galaxyTreeSed,localgalaxyid),h5::datatype::std_i32be,"localgalaxyid");
        der.commit( mem_type, file_type );

    }

}
