#ifndef tao_sage_hh
#define tao_sage_hh

#include <libhpc/h5/datatype.hh>
#include <string>
#include <vector>

struct GALAXY_OUTPUT
{
    int SnapNum;
    int Type;

    long long GalaxyIndex;
    long long CentralGalaxyIndex;
    int SAGEHaloIndex;
    int SAGETreeIndex;
    long long SimulationHaloIndex;

    int mergeType; // 0=none; 1=minor merger; 2=major merger; 3=disk instability;
                   // 4=disrupt to ICS
    int mergeIntoID;
    int mergeIntoSnapNum;
    float dT;

    // (sub)halo properties
    float Pos[3];
    float Vel[3];
    float Spin[3];
    int Len;
    float Mvir;
    float CentralMvir;
    float Rvir;
    float Vvir;
    float Vmax;
    float VelDisp;

    // baryonic reservoirs
    float ColdGas;
    float StellarMass;
    float BulgeMass;
    float HotGas;
    float EjectedMass;
    float BlackHoleMass;
    float ICS;

    // metals
    float MetalsColdGas;
    float MetalsStellarMass;
    float MetalsBulgeMass;
    float MetalsHotGas;
    float MetalsEjectedMass;
    float MetalsICS;

    // to calculate magnitudes
    float SfrDisk;
    float SfrBulge;
    float SfrDiskZ;
    float SfrBulgeZ;

    // misc
    float DiskScaleRadius;
    float Cooling;
    float Heating;
    float QuasarModeBHaccretionMass;
    float TimeOfLastMajorMerger;
    float TimeOfLastMinorMerger;
    float OutflowRate;

    // infall properties
    float infallMvir;
    float infallVvir;
    float infallVmax;
};

namespace sage
{

// Field information for columnar storage
struct FieldInfo
{
    std::string name;        // HDF5 dataset name
    hpc::h5::datatype type;  // HDF5 datatype
    size_t offset;           // Offset in sage::galaxy struct
    std::string description; // Field description (from HDF5 attributes)
    std::string units;       // Field units (from HDF5 attributes)
};

struct galaxy
{
    int snapshot;
    int type;

    long long galaxy_idx;
    long long central_galaxy_idx;
    int sage_halo_idx;
    int sage_tree_idx;
    long long simulation_halo_idx;

    // LUKE: See struct GALAXY.
    int local_index;
    long long global_index;
    int descendant;
    long long global_descendant;
    int subsize;

    int merge_type; // 0=none; 1=minor merger; 2=major merger; 3=disk instability;
                    // 4=disrupt to ICS
    int merge_into_id;
    int merge_into_snapshot;
    float dt;

    // properties of subhalo at the last time this galaxy was a central galaaxy
    float pos[3];
    float vel[3];
    float spin[3];
    int num_particles;
    float mvir;
    float central_mvir;
    float rvir;
    float vvir;
    float vmax;
    float vel_disp;

    // baryonic reservoirs
    float cold_gas;
    float stellar_mass;
    float bulge_mass;
    float hot_gas;
    float ejected_mass;
    float blackhole_mass;
    float ics;

    // metals
    float metals_cold_gas;
    float metals_stellar_mass;
    float metals_bulge_mass;
    float metals_hot_gas;
    float metals_ejected_mass;
    float metals_ics;

    // to calculate magnitudes
    float sfr_disk;
    float sfr_bulge;
    float sfr_disk_z;
    float sfr_bulge_z;

    // misc
    float disk_scale_radius;
    float cooling;
    float heating;
    float quasar_mode_bh_accretion_mass;
    float time_of_last_major_merger;
    float time_of_last_minor_merger;
    float outflow_rate;

    float infall_mvir;
    float infall_vvir;
    float infall_vmax;
};

struct galaxyTree
{
    int objecttype;
    float pos[3];
    float velx;
    float vely;
    float velz;
    int snapnum;
    long long globalindex;
    int descendant;
    int mergetype;
    float dt;
    float sfrdisk;
    float sfrbulge;
    float sfrdiskz;
    float sfrbulgez;
    int treeindex;
    int localindex;
    long long globaldescendant;
    int subsize;
    float coldgas;
    float metalscoldgas;
    float diskscaleradius;
    long long galaxyindex;
    long long centralgalaxyindex;
    long long simulationhaloindex;
    int mergeintoid;
    int mergeintosnapnum;
    float spin_x;
    float spin_y;
    float spin_z;
    int len;
    float mvir;
    float centralmvir;
    float rvir;
    float vvir;
    float vmax;
    float veldisp;
    float stellarmass;
    float bulgemass;
    float hotgas;
    float ejectedmass;
    float blackholemass;
    float ics;
    float metalsstellarmass;
    float metalsbulgemass;
    float metalshotgas;
    float metalsejectedmass;
    float metalsics;
    float cooling;
    float heating;
    float quasarmodebhaccretionmass;
    float timeoflastmajormerger;
    float timeoflastminormerger;
    float outflowrate;
    float infallmvir;
    float infallvvir;
    float infallvmax;
    float totsfr;
    float vpeak;
    int isflyby;
};

/*
 * Derived from existing RDB case (That is, these are the columns in the RDB
 * Tree tables)
 */
struct galaxyTreeSed
{
    int objecttype;
    float pos[3];
    float velx;
    float vely;
    float velz;
    int snapnum;
    long long globalindex;
    int descendant;
    int mergetype;
    float dt;
    float sfrdisk;
    float sfrbulge;
    float sfrdiskz;
    float sfrbulgez;
    int treeindex;
    int localindex;
    long long globaldescendant;
    int subsize;
    float coldgas;
    float metalscoldgas;
    float diskscaleradius;
    long long galaxyindex;
    long long centralgalaxyindex;
    long long simulationhaloindex;
    int mergeintoid;
    int mergeintosnapnum;
    float spin_x;
    float spin_y;
    float spin_z;
    int len;
    float mvir;
    float centralmvir;
    float rvir;
    float vvir;
    float vmax;
    float veldisp;
    float stellarmass;
    float bulgemass;
    float hotgas;
    float ejectedmass;
    float blackholemass;
    float ics;
    float metalsstellarmass;
    float metalsbulgemass;
    float metalshotgas;
    float metalsejectedmass;
    float metalsics;
    float cooling;
    float heating;
    float quasarmodebhaccretionmass;
    float timeoflastmajormerger;
    float timeoflastminormerger;
    float outflowrate;
    float infallmvir;
    float infallvvir;
    float infallvmax;
    float totsfr;
    float vpeak;
    int isflyby;

    long long globaltreeid;
    long long breadthfirst_traversalorder;
    long long depthfirst_traversalorder;
    long long subtree_count;
    int localgalaxyid;
};

} // namespace sage

#endif
