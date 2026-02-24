#ifndef tao_base_query_hh
#define tao_base_query_hh

#include <libhpc/system/string.hh>
#include <libhpc/system/view.hh>
#include <set>
#include <string>
#include <vector>

namespace tao
{

/// Query objects are used to indicate which fields should be extracted
/// from data sources.
///
/// @tparam T The type of real values (to be removed soon).
///
template <class T>
class query
{
public:
    typedef T real_type;

public:
    /// Initialise an instance. By default, a suite of base fields are given
    /// to the object that are known to be required:
    ///
    /// @sa add_base_output_fields
    ///
    query()
    {
        ///         add_base_output_fields();
    }

    /// Clear all fields from the query.
    ///
    void clear()
    {
        _of_set.clear();
        _out_fields.clear();
    }

    /// Add the set of base output fields known to be required. The fields
    /// that are added are:
    ///
    ///   - posx
    ///   - posy
    ///   - posz
    ///   - velx
    ///   - vely
    ///   - velz
    ///   - snapnum
    ///   - sfrdisk
    ///   - sfrbulge
    ///   - sfrdiskz
    ///   - sfrbulgez
    ///   - diskradius
    ///   - coldgas
    ///   - metalscoldgas
    ///
    /// Note: globaltreeid, localgalaxyid, globalindex are computed fields
    /// that are not stored in data/ group - they come from lightcone/data
    /// compound.
    ///
    void add_base_output_fields()
    {
        add_output_field("posx");
        add_output_field("posy");
        add_output_field("posz");
        add_output_field("velx");
        add_output_field("vely");
        add_output_field("velz");
        add_output_field("snapnum");
        // Note: globaltreeid, localgalaxyid, globalindex removed - computed fields
        // not in data/
        add_output_field("sfrdisk");
        add_output_field("sfrbulge");
        add_output_field("sfrdiskz");
        add_output_field("sfrbulgez");
        add_output_field("diskradius"); // Note: SAGE produces DiskRadius, not DiskScaleRadius
        add_output_field("coldgas");
        add_output_field("metalscoldgas");
        add_output_field("distance");
    }

    /// Add the set of base output fields for particle-based data. The fields
    /// that are added are:
    ///
    ///   - posx
    ///   - posy
    ///   - posz
    ///   - velx
    ///   - vely
    ///   - velz
    ///   - snapshot
    ///   - distance
    ///
    void add_base_output_fields_pb()
    {
        add_output_field("posx");
        add_output_field("posy");
        add_output_field("posz");
        add_output_field("velx");
        add_output_field("vely");
        add_output_field("velz");
        add_output_field("snapshot");
        add_output_field("distance");
    }

    /// Add an output field to the query.
    ///
    /// @param[in] field The name of the field to add.
    ///
    void add_output_field(std::string const& field)
    {
        _out_fields.clear();
        _of_set.insert(hpc::to_lower_copy(field));
    }

    const hpc::view<std::vector<std::string>> output_fields()
    {
        if (_out_fields.empty())
        {
            _out_fields.resize(_of_set.size());
            std::copy(_of_set.begin(), _of_set.end(), _out_fields.begin());
        }
        return _out_fields;
    }

protected:
    std::set<std::string> _of_set;
    std::vector<std::string> _out_fields;
};

} // namespace tao

#endif
