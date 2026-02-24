#ifndef tao_base_backend_hh
#define tao_base_backend_hh

#include "batch.hh"
#include "query.hh"
#include "simulation.hh"
#include "types.hh"
#include <boost/property_tree/ptree.hpp>

namespace tao
{

///
/// Base class for TAO backends. A backend defines access to
/// data sources for astronomical datasets. For example, the
/// datasets could be stored in a database (relational or
/// otherwise), in HDF5 files or even simple text files. The
/// backends abstract the storage from the processing.
///
class backend
{
public:
    ///
    /// Construct with a simulation.
    ///
    /// @param[in] sim The simulation this backend refers to.
    ///
    backend(tao::simulation const* sim = nullptr);

    ///
    /// Set the simulation. Each backend is constructed to
    /// represent a particular dataset, which itself must refer
    /// to a simulation.
    ///
    /// @param[in] sim The simulation this backend refers to.
    ///
    virtual void set_simulation(tao::simulation const* sim);

    ///
    /// Get the simulation.
    ///
    /// @returns The simulation this backend refers to.
    ///
    tao::simulation const* simulation() const;

    ///
    /// Load the simulation details from the dataset. In most
    /// circumstances the simulation metadata will be stored
    /// along with the simulation data.
    ///
    /// @returns The loaded simulation.
    ///
    virtual tao::simulation const* load_simulation() = 0;

    virtual void load_checkpoint(boost::property_tree::ptree const& pt);

    // For some backends there are different fields expected.
    // Override in derived classes if additional fields are needed.
    virtual void add_conditional_fields(query<real_type>& qry);

    void init_batch(batch<real_type>& bat, query<real_type>& qry) const;

    std::string getFieldDescription(std::string& field)
    {
        std::string ask(field);
        std::cout << "ASKING " << ask << std::endl;
        std::cout << "ASKINGs " << _field_description.size() << std::endl;
        return _field_description[ask];
    }
    std::string getFieldUnits(std::string& field) { return _field_units[field]; }
    std::string getFieldGroup(std::string& field) { return _field_group[field]; }
    int getFieldOrder(std::string& field) { return _field_order[field]; }

protected:
    tao::simulation const* _sim;
    std::unordered_map<std::string, std::string> _field_map;
    std::map<std::string, typename batch<real_type>::field_value_type> _field_types;
    std::unordered_map<std::string, std::string> _field_description;
    std::unordered_map<std::string, std::string> _field_units;
    std::unordered_map<std::string, std::string> _field_group;
    std::unordered_map<std::string, int> _field_order;
};

} // namespace tao

#endif
