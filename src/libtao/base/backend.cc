#include "backend.hh"

namespace tao {

    backend::backend(tao::simulation const *sim)
            : _sim(sim) {
    }

    void
    backend::set_simulation(tao::simulation const *sim) {
        _sim = sim;
    }

    tao::simulation const *
    backend::simulation() const {
        return _sim;
    }

    void
    backend::load_checkpoint(boost::property_tree::ptree const &pt) {
    }

    // For some backends there are different fields expected.
    // For example, the kdtree requires a field subtree_count
    void
    backend::add_conditional_fields(query<real_type> &qry) {
    }


    /// Initialise a batch object from a query object. The batch object
    /// will be populated with entries for all known fields.
    ///
    /// @param[in] bat The batch object to be populated.
    /// @param[in] qry A query object to take fields from.
    ///
    void
    backend::init_batch(batch<real_type> &bat,
                        query<real_type> &qry) const {
        // Add fields from the query object.
        for (auto const &field : qry.output_fields()) {
            // Check that the field actually exists. Due to calculated fields
            // it may not actually be on the database.
            if (this->_field_map.find(field) != this->_field_map.end()) {
                // Field is an alias - verify the target field exists before using it
                const std::string& mapped_field = this->_field_map.at(field);
                if (this->_field_types.find(mapped_field) != this->_field_types.end()) {
                    bat.set_scalar(field, _field_types.at(mapped_field));
                } else {
                    std::cerr << "WARNING: Alias '" << field << "' maps to '" << mapped_field
                              << "' but that field doesn't exist in _field_types" << std::endl;
                }
            } else if (this->_field_types.find(field) != this->_field_types.end()) {
                // Field is a direct field name - use it directly to look up type
                bat.set_scalar(field, _field_types.at(field));
            }
        }

        // Add fields from the density object.
        for (auto const &field : qry.density_fields()) {
            // All column density base fields are already added as output fields
            // So here, only add the column density itself, which is calculated
            auto density_field = field + "_column_density";
            bat.template set_scalar<real_type>(density_field);
        }

        // Add fields that will need to be calculated.
        bat.template set_scalar<real_type>("redshift_cosmological");
        bat.template set_scalar<real_type>("redshift_observed");
        bat.template set_scalar<real_type>("column_density");
        bat.template set_scalar<real_type>("ra");
        bat.template set_scalar<real_type>("dec");
        // we need to keep it here for boxes / lightcones for now.
        bat.template set_scalar<real_type>("distance");
        bat.template set_scalar<real_type>("sfr");
        bat.template set_scalar<real_type>("distance_from_beam");
        // Both required for column density calculations
        bat.template set_scalar<real_type>("mass");
        bat.template set_scalar<real_type>("density");
    }

}
