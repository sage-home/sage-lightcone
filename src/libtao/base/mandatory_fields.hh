#ifndef tao_mandatory_fields_hh
#define tao_mandatory_fields_hh

#include <set>
#include <string>
#include <vector>

namespace tao
{

/// SAGE field names (exact CamelCase names as they appear in SAGE HDF5 output)
/// These fields are required for the sage2kdtree -> cli_lightcone pipeline to
/// function.
namespace sage_fields
{
// Position (critical for KD-tree spatial indexing)
constexpr const char* POSX = "Posx";
constexpr const char* POSY = "Posy";
constexpr const char* POSZ = "Posz";

// Velocity (critical for redshift calculations in lightcone)
constexpr const char* VELX = "Velx";
constexpr const char* VELY = "Vely";
constexpr const char* VELZ = "Velz";

// Snapshot number
constexpr const char* SNAPNUM = "SnapNum";

// Galaxy identifiers (needed for tree structure)
constexpr const char* GALAXY_INDEX = "GalaxyIndex";
constexpr const char* CENTRAL_GALAXY_INDEX = "CentralGalaxyIndex";

} // namespace sage_fields

/// Fields computed by the pipeline (not in SAGE output)
/// These use lowercase names by convention.
namespace computed_fields
{
// Tree traversal indices (added in Phase 2)
constexpr const char* GLOBAL_TREE_ID = "globaltreeid";
constexpr const char* BREADTH_FIRST_ORDER = "breadthfirst_traversalorder";
constexpr const char* DEPTH_FIRST_ORDER = "depthfirst_traversalorder";
constexpr const char* LOCAL_GALAXY_ID = "localgalaxyid";
constexpr const char* SUBTREE_COUNT = "subtree_count";

// Index fields (added during traversal)
constexpr const char* LOCAL_INDEX = "local_index";
constexpr const char* GLOBAL_INDEX = "global_index";
constexpr const char* DESCENDANT = "descendant";
constexpr const char* GLOBAL_DESCENDANT = "global_descendant";
constexpr const char* SUBSIZE = "subsize";

// Lowercase alias for snapnum (for backward compatibility)
constexpr const char* SNAPNUM_LOWER = "snapnum";

// Derived spatial fields (computed during lightcone extraction)
constexpr const char* DISTANCE = "distance";
constexpr const char* RA = "ra";
constexpr const char* DEC = "dec";

// Derived redshift fields
constexpr const char* REDSHIFT_COSMO = "redshift_cosmological";
constexpr const char* REDSHIFT_OBS = "redshift_observed";

// Derived star formation rate (SfrDisk + SfrBulge)
constexpr const char* SFR = "sfr";
} // namespace computed_fields

/// Validate that all mandatory SAGE fields are present in the available fields
/// list. Uses case-insensitive matching to be flexible with field name
/// variations.
///
/// @param available_fields List of field names available in the HDF5 file
/// @param missing_fields Output vector of missing field names (will be cleared
/// first)
/// @return true if all mandatory fields are present, false otherwise
bool validate_sage_fields(const std::vector<std::string>& available_fields,
                          std::vector<std::string>& missing_fields);

/// Get all mandatory SAGE field names as a set
std::set<std::string> get_mandatory_sage_fields();

/// Get KD-tree-critical field names (subset of mandatory fields)
std::set<std::string> get_kdtree_required_fields();

/// Case-insensitive string comparison utility
bool iequals(const std::string& a, const std::string& b);

} // namespace tao

#endif
