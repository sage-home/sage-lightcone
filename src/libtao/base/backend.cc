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
    // Override in derived classes if additional fields are needed.
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

        // Add fields that will need to be calculated.
        bat.template set_scalar<real_type>("redshift_cosmological");
        bat.template set_scalar<real_type>("cosmological_redshift");
        bat.template set_scalar<real_type>("redshift_observed");
        bat.template set_scalar<real_type>("observed_redshift");
        bat.template set_scalar<real_type>("ra");
        bat.template set_scalar<real_type>("dec");
        // we need to keep it here for boxes / lightcones for now.
        bat.template set_scalar<real_type>("distance");
        bat.template set_scalar<real_type>("sfr");

    }

}
