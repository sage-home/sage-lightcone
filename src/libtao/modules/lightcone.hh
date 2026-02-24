#ifndef tao_modules_lightcone_hh
#define tao_modules_lightcone_hh

#include "libtao/base/base.hh"
#include <boost/optional.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <iostream>
#include <libhpc/system/math.hh>
#include <set>

namespace tao
{
namespace modules
{

///
/// Lightcone science module.
///
template <class Backend>
class lightcone : public module<Backend>
{
public:
    typedef Backend backend_type;
    typedef module<backend_type> module_type;

    // Type of geometry to use.
    enum geometry_type
    {
        CONE,
        BOX,
        PENCIL
    };

    // Factory function used to create a new module.
    static module_type* factory(std::string const& name) { return new lightcone(name); }

public:
    lightcone(std::string const& name = std::string())
        : module_type(name)
        , _my_be(false)
        , _be(0)
    {
    }

    virtual ~lightcone()
    {
        if (_my_be)
            delete _be;
    }

    void set_backend(backend_type* be)
    {
        _my_be = (be == 0);
        _be = be;
    }

    virtual backend_type* backend() { return _be; }

    tao::simulation const* simulation() const { return _sim; }

    ///
    ///
    ///
    virtual void initialise(const cli_dict& global_cli_dict,
                            boost::optional<boost::property_tree::ptree> checkpoint =
                                boost::optional<boost::property_tree::ptree>())
    {
        // Don't initialise if we're already doing so.
        if (this->_init)
            return;
        module_type::initialise(global_cli_dict, checkpoint);

        LOGILN("Initialising lightcone module.", setindent(2));

        // If we have been given an existing backend then use that.
        if (!_be)
        {
            _my_be = true;
            _be = new backend_type;
        }
        // Always call connect() to ensure backend is properly initialized
        _be->connect(global_cli_dict);

        // Read all my options.
        _read_options(global_cli_dict);

        // Make sure the local batch object is prepared. I need to
        // do this here because other modules will likely need
        // to access and cache values.

        _qry.add_base_output_fields();
        _be->add_conditional_fields(_qry);
        _be->init_batch(_bat, _qry);

        // Set the random seed now, so that we actually start from
        // the seed given.
        _eng.seed(_rng_seed);
        LOGILN("SETTING SEED: ", _rng_seed);

        // Show in the logs what we're querying.
        LOGILN("Querying the following fields: ", _qry.output_fields());

        // We need to calculate the total number of tiles in the
        // cone that we will be querying so we can give a meaningful
        // progress indicator (I say meaningful...).
        _num_gals = 0;
        _num_gals = mpi::comm::world.all_reduce(_num_gals);

        // Initialise the progress tracker.
#ifndef NO_PROGRESS
        _prog.start(_num_gals);
#endif
        // Set the seed again after iterating over the tables once.
        _eng.seed(_rng_seed);

        // Show in the logs what we're querying.
        LOGILN("Querying the following fields: ", _qry.output_fields());

        // Fast-forward the checkpoint on our backend.
        if ((bool)checkpoint)
            _be->load_checkpoint(*checkpoint);

        LOGILN("Done.", setindent(-2));
    }

    ///
    /// Run the module.
    ///
    virtual void execute()
    {
        // Is this my first time through? If so begin iterating.
        if (this->_it == 0)
        {
            // First of all, fix up my filter. I need to postpone it until
            // here because the batch object is only full of all the possible
            // fields now.
            if (!_filt.field_name().empty() && _bat.has_field(_filt.field_name()))
                _filt.set_type(_bat.get_field_type(_filt.field_name()));

            // Restart the iterators.
            _c_it = _lc.galaxy_begin(_qry, *_be, &_bat, &_filt);

            // Don't forget to include any tables initially skipped in the
            // progress count.
            _n_processed_gals = _c_it.n_processed_gals();

            // Send through the first set of galaxies to the progress tracker.
#ifndef NO_PROGRESS
            if (_n_processed_gals)
                _prog.append(_n_processed_gals);
#endif
        }
        else
        {
            ++_c_it; // Increment the iterator to the next galaxy doing the work as we
                     // go.

            // Track progress.
            double diff = (double)(_c_it.n_processed_gals() - _n_processed_gals);
            _n_processed_gals += diff;
#ifndef NO_PROGRESS
            if (diff)
                _prog.append(diff);
#endif
            // Update checkpointing info.
            if (_c_it.should_checkpoint())
            {
                this->save_checkpoint();
            }
        }

        // Check for completion.
        if (_c_it == _lc.galaxy_end(_qry, *_be))
        {
            this->_complete = true;
        }

        // Once complete, terminate the progress tracker.
#ifndef NO_PROGRESS
        if (this->_complete)
            _prog.stop();
#endif
    }

    ///
    ///
    ///
    virtual tao::batch<real_type>& batch() { return _bat; }

    virtual boost::optional<boost::any> find_attribute(std::string const& name)
    {
        if (name == "simulation")
            return boost::any((tao::simulation const*)_sim);
        else if (name == "filter")
        {
            return boost::any(&((filter const&)_filt));
        }
        else
            return module_type::find_attribute(name);
    }

    geometry_type geometry() const { return _geom; }

    int random_seed() const { return _rng_seed; }

    bool tile_repetition_random() const { return !_unique; }

    tao::lightcone const& base_lightcone() const { return (_lc); }

    real_type box_size() const { return _box_size; }

    real_type box_redshift() const { return _box_z; }

    virtual void log_metrics()
    {
        module_type::log_metrics();
        LOGILN(this->_name, " number of tiles: ", _num_tiles);
    }

    void save_checkpoint()
    {
#ifndef NO_CHECKPOINT
        boost::property_tree::ptree pt;
        this->checkpoint(pt);
        std::string fn(global_cli_dict._outdir + "/log/checkpoint.");
        fn += std::to_string(mpi::comm::world.rank());
        fn += ".json";
        std::ofstream out(fn, std::ios::trunc);
        boost::property_tree::write_json(out, pt);
#endif
    }

    virtual void do_checkpoint(boost::property_tree::ptree& pt)
    {
        if (_geom == CONE)
            _c_it.save_checkpoint(pt);
        // else if( _geom == PENCIL )
        //    _p_it.save_checkpoint( pt );
        else
            _b_it.save_checkpoint(pt);
    }

protected:
    void _read_options(const cli_dict& global_cli_dict)
    {
        bool use_cli_exclusively = true;

        // Cache the log directory for checkpointing.
        _logdir = global_cli_dict._outdir + "/log";

        // Load simulation information from database.
        _sim = _be->load_simulation();

        // Extract the random number generator seed and set it.
        // Look for single cone rng-seed style
        boost::optional<int> rng_seed = global_cli_dict._rng_seed;
        if (rng_seed != 0)
            _rng_seed = *rng_seed;
        else
            _rng_seed = rand();

        mpi::comm::world.bcast<int>(_rng_seed, 0);
        LOGILN("Random seed: ", _rng_seed);

        // Get tile repetition type.
        _unique = global_cli_dict._unique;

        _geom = CONE;

        // Redshift ranges.
        auto const& redshifts = _sim->redshifts();
        real_type z1 = redshifts.begin()->second;
        real_type z2 = redshifts.rbegin()->second;
        real_type snap_z_max = std::max(z1, z2);
        real_type snap_z_min = std::min(z1, z2);
        real_type max_z;
        real_type min_z;
        max_z = std::min(global_cli_dict._zmax, snap_z_max);
        min_z = std::max(global_cli_dict._zmin, snap_z_min);

        LOGILN("Redshift range: [", min_z, ", ", max_z, ")");

        // Right ascension.
        real_type min_ra;
        real_type max_ra;
        min_ra = global_cli_dict._ramin;
        max_ra = global_cli_dict._ramax;

        LOGILN("Right ascension range: [", min_ra, ", ", max_ra, ")");

        // Declination.
        real_type min_dec;
        real_type max_dec;
        min_dec = global_cli_dict._decmin;
        max_dec = global_cli_dict._decmax;

        LOGILN("Declination range: [", min_dec, ", ", max_dec, ")");

        _lc.set_simulation(_sim);
        _lc.set_geometry(min_ra, max_ra, min_dec, max_dec, max_z, min_z);

        _lc.set_random(!_unique, _rng_seed, &_eng);

        // Prepare the origin if we're running a unique cone.
        if (_unique)
        {
            if (!_calc_origin())
            {
                LOGILN("No unique cones can be computed for ", global_cli_dict._dataset);
                LOGILN("with the given ranges:");
                LOGILN("Dec: [", min_dec, ", ", max_dec);
                LOGILN("RA: [", min_ra, ", ", max_ra);
                LOGILN("Redshift: [", min_z, ", ", max_z, ")");
                LOGILN("Please check the ranges and/or the input hdf5 file.");
                LOGILN("Try smaller ranges for Dec, RA, and/or redshift");
                std::cerr << "Error: No unique cones can be computed for "
                          << global_cli_dict._dataset << std::endl;
                std::cerr << "  with the given ranges:" << std::endl;
                std::cerr << "  Dec: [" << min_dec << ", " << max_dec << ")" << std::endl;
                std::cerr << "  RA: [" << min_ra << ", " << max_ra << ")" << std::endl;
                std::cerr << "  Redshift: [" << min_z << ", " << max_z << ")" << std::endl;
                std::cerr << "  Please check the ranges and/or the input hdf5 file." << std::endl;
                std::cerr << "  Try smaller ranges for Dec, RA, and/or redshift." << std::endl;
                exit(1);
            }
        }

        // Output field information.
        // Get all field names directly from HDF5 file /data/* datasets
        std::string input_hdf5_file = global_cli_dict._dataset;
        std::vector<std::string> all_field_names;
        {
            hpc::h5::file input_file(input_hdf5_file, H5F_ACC_RDONLY);
            hpc::h5::group data_group;
            input_file.open_group("data", data_group);

            H5G_info_t group_info;
            H5Gget_info(data_group.id(), &group_info);

            for (hsize_t i = 0; i < group_info.nlinks; i++)
            {
                char name_buf[256];
                H5Lget_name_by_idx(data_group.id(), ".", H5_INDEX_NAME, H5_ITER_NATIVE, i, name_buf,
                                   sizeof(name_buf), H5P_DEFAULT);
                std::string field_name(name_buf);
                // Convert to lowercase for case-insensitive comparison
                std::transform(field_name.begin(), field_name.end(), field_name.begin(), ::tolower);
                all_field_names.push_back(field_name);
            }
        }

        LOGILN("Using output fields specified in the command line.");

        // Known calculated fields that are computed during lightcone extraction
        // (not present in input HDF5 file)
        static const std::set<std::string> calculated_fields = {
            "ra", "dec", "distance", "redshift_cosmological", "redshift_observed", "sfr"};

        // Use the output fields specified in the command line.
        if (global_cli_dict._output_fields.size() == 0)
        {
            LOGILN("No output fields specified via CLI, using all fields from the "
                   "input hdf5 file.");
            // If no output fields specified, use all fields from the input hdf5 file.
            for (auto const& field_name : all_field_names)
            {
                // Note: Without XML, we can't filter by group. All fields from /data
                // are used.
                _qry.add_output_field(field_name);
                LOGILN("CLI.Adding field ", field_name);
            }

            // Also add calculated fields by default
            LOGILN("Adding calculated fields by default.");
            for (const auto& field : calculated_fields)
            {
                _qry.add_output_field(field);
                LOGILN("CLI.Adding calculated field ", field);
            }
        }
        else
        {
            std::vector<std::string> invalid_fields;

            // Helper to process a single field
            auto process_field = [&](const std::string& field) {
                // Convert requested field to lowercase for case-insensitive comparison
                std::string field_lower = field;
                std::transform(field_lower.begin(), field_lower.end(), field_lower.begin(),
                               ::tolower);

                // Allow known calculated fields even if not in input file
                if (calculated_fields.count(field_lower) > 0)
                {
                    _qry.add_output_field(field);
                    LOGILN("CLI.Adding calculated field ", field);
                    return;
                }

                auto it = std::find(all_field_names.begin(), all_field_names.end(), field_lower);
                if (it == all_field_names.end())
                {
                    LOGILN("Field '", field, "' not found in input hdf5 file.");
                    invalid_fields.push_back(field);
                    return;
                }
                _qry.add_output_field(field);
                LOGILN("CLI.Adding field ", field);
            };

            for (auto const& field_arg : global_cli_dict._output_fields)
            {
                // Handle space-separated fields in a single argument
                // (e.g., --outfields "ra dec sfr")
                std::istringstream iss(field_arg);
                std::string field;
                while (iss >> field)
                {
                    process_field(field);
                }
            }

            if (!invalid_fields.empty())
            {
                std::cerr << "Error: The following requested output fields were not "
                             "found:"
                          << std::endl;
                for (const auto& f : invalid_fields)
                {
                    std::cerr << "  - " << f << std::endl;
                }
                std::cerr << std::endl;
                std::cerr << "Available fields in input HDF5 file:" << std::endl;
                for (const auto& f : all_field_names)
                {
                    std::cerr << "  - " << f << std::endl;
                }
                std::cerr << std::endl;
                std::cerr << "Available calculated fields:" << std::endl;
                for (const auto& f : calculated_fields)
                {
                    std::cerr << "  - " << f << std::endl;
                }
                exit(1);
            }
        }

        // TODO: CLI option hardcoded for now
        // Filter information.
        _filt.set_field(global_cli_dict._filter_field, global_cli_dict._filter_min,
                        global_cli_dict._filter_max);
        LOGILN("CLI.Filter name: ", global_cli_dict._filter_field);
        LOGILN("CLI.Filter range: [", global_cli_dict._filter_min, ", ",
               global_cli_dict._filter_max, ")");
        LOGILN("Extracting fields: ", _qry.output_fields());

        LOGILN("Summary of Extracting fields: ", _qry.output_fields().size(), " ",
               use_cli_exclusively);
    }

    bool _calc_origin()
    {
        // EXCEPT( _lc.min_ra() == 0.0, "Cannot compute multiple unique cones when "
        //         "minimum RA is not 0." );
        // EXCEPT( _lc.min_dec() == 0.0, "Cannot compute multiple unique cones when
        // "
        //         "minimum DEC is not 0." );

        // Get the subjobindex.
        unsigned sub_idx = 0;
        LOGILN("Subcone index: ", sub_idx);

        // Check this index is okay.
        unsigned max_subcones = tao::calc_max_subcones(_lc);
        if (max_subcones == 0)
            return false;
        EXCEPT(sub_idx < max_subcones, "Subcone with index ", sub_idx,
               " is outside the maximum range of ", max_subcones);

        // Calculate subcone viewing angle and origin.
        real_type view_angle = *tao::calc_subcone_angle(_lc);
        std::array<real_type, 3> ori = tao::calc_subcone_origin<real_type>(_lc, sub_idx);

        // Set angle and origin.
        _lc.set_ra_offset(view_angle);
        _lc.set_dec_offset(-_lc.min_dec());
        _lc.set_origin(ori);
        LOGILN("RA offset set to: ", to_degrees(_lc.ra_offset()));
        LOGILN("Dec offset set to: ", to_degrees(_lc.dec_offset()));
        LOGILN("Origin set to: (", ori[0], ", ", ori[1], ", ", ori[2], ")");
        return true;
    }

protected:
    geometry_type _geom;
    real_type _box_size;
    real_type _box_z;
    int _rng_seed;
    engine_type _eng;
    bool _unique;

    tao::simulation const* _sim;
    query<real_type> _qry;
    filter _filt;
    tao::lightcone _lc;
    box<real_type> _box;
    bool _my_be;
    backend_type* _be;
    typename backend_type::lightcone_galaxy_iterator _c_it;
    typename backend_type::box_galaxy_iterator _b_it;
    tao::batch<real_type> _bat;

    unsigned _num_tiles;
    unsigned _num_tbls;
    double _num_gals;
    unsigned long long _n_processed_gals;
    progress _prog;

    std::string _logdir;
};

} // namespace modules
} // namespace tao

#endif
