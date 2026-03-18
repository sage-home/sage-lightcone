#include "kdapplication.hh"
#include "libtao/modules/modules.hh"
#include <boost/algorithm/string/trim.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <stdexcept>
#include <string>
// #include <pugixml.hpp>

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

using namespace hpc;
// using namespace pugi;
using namespace boost;
using namespace std;
using namespace boost::filesystem;

namespace tao
{

KdApplication::KdApplication(int argc, char* argv[])
    : mpi::application(argc, argv)
{
    // Setup some options.
    hpc::application::options().add_options()(
        "dataset,d", hpc::po::value<std::string>(&_global_cli_dict._dataset), "dataset file")(
        "decmin", hpc::po::value<double>(&_global_cli_dict._decmin)->default_value(0.0),
        "minimum declination")(
        "decmax", hpc::po::value<double>(&_global_cli_dict._decmax)->default_value(10.0),
        "maximum declination")("ramin",
                               hpc::po::value<double>(&_global_cli_dict._ramin)->default_value(0.0),
                               "minimum right ascension")(
        "ramax", hpc::po::value<double>(&_global_cli_dict._ramax)->default_value(10.0),
        "maximum right ascension")(
        "zmin", hpc::po::value<double>(&_global_cli_dict._zmin)->default_value(0.0),
        "minimum redshift")("zmax",
                            hpc::po::value<double>(&_global_cli_dict._zmax)->default_value(1.0),
                            "maximum redshift")(
        "unique", hpc::po::value<bool>(&_global_cli_dict._unique)->default_value(false),
        "unique lightcone if true, random lightcone if false")(
        "seed", hpc::po::value<int>(&_global_cli_dict._rng_seed)->default_value(0),
        "random number generator seed")(
        "filterfield,f",
        hpc::po::value<std::string>(&_global_cli_dict._filter_field)->default_value("stellar_mass"),
        "field to filter on")(
        "filtermin", hpc::po::value<std::string>(&_global_cli_dict._filter_min)->default_value("0"),
        "minimum value for filter field")(
        "filtermax", hpc::po::value<std::string>(&_global_cli_dict._filter_max)->default_value(""),
        "maximum value for filter field")(
        "outfile,o",
        hpc::po::value<std::string>(&_global_cli_dict._outfile)->default_value("output.hdf5"),
        "output file name")(
        "outfields",
        hpc::po::value<std::vector<std::string>>(&_global_cli_dict._output_fields)->multitoken(),
        "fields to output (unquoted and space separated)")(
        "outdir", hpc::po::value<std::string>(&_global_cli_dict._outdir)->default_value("output"),
        "directory for output files")

        ("centralgalaxies",
         hpc::po::bool_switch(&_global_cli_dict._central_galaxies)->default_value(false),
         "Central galaxies mode: include all satellites when their central is in lightcone")(
            "verbose,v",
            "verbose output")("debug", "debug output")("version", "print version number and exit");

    // Parse options.
    parse_options(argc, argv);

    if (_vm.count("version"))
    {
        std::cout << "1.0.0" << std::endl;
        exit(0);
    }
    // Validation of the command line arguments.
    if (_vm.count("dataset") == 0)
    {
        std::cerr << "Error: --dataset is required." << std::endl;
        throw silent_terminate();
    }

    if (mpi::comm::world.rank() == 0)
    {
        // validation of the declination and right ascension ranges
        if (validate(_global_cli_dict))
        {
            if (!exists(_global_cli_dict._outdir))
            {
                std::cout << "Creating output directory: " << _global_cli_dict._outdir << std::endl;
                create_directories(_global_cli_dict._outdir);
            }
        }
        else
        {
            std::cerr << "Error: Invalid command line arguments." << std::endl;
            throw silent_terminate();
        }
        std::cout << "dataset=" << _global_cli_dict._dataset << std::endl;
        std::cout << "decmin=" << _global_cli_dict._decmin << std::endl;
        std::cout << "decmax=" << _global_cli_dict._decmax << std::endl;
        std::cout << "ramin=" << _global_cli_dict._ramin << std::endl;
        std::cout << "ramax=" << _global_cli_dict._ramax << std::endl;
        std::cout << "zmin=" << _global_cli_dict._zmin << std::endl;
        std::cout << "zmax=" << _global_cli_dict._zmax << std::endl;
        std::cout << "unique=" << _global_cli_dict._unique << std::endl;
        std::cout << "seed=" << _global_cli_dict._rng_seed << std::endl;
        std::cout << "filterfield=" << _global_cli_dict._filter_field << std::endl;
        std::cout << "filtermin=" << _global_cli_dict._filter_min << std::endl;
        std::cout << "filtermax=" << _global_cli_dict._filter_max << std::endl;
        std::cout << "outfile=" << _global_cli_dict._outfile << std::endl;
        std::cout << "outdir=" << _global_cli_dict._outdir << std::endl;
        std::cout << "centralgalaxies=" << _global_cli_dict._central_galaxies << std::endl;
    }
    if (_global_cli_dict._rng_seed == 0)
    {
        _global_cli_dict._rng_seed = rand();
    }

    if (mpi::comm::world.rank() == 0)
    {
        if (_global_cli_dict._unique)
            std::cout << "unique lightcone" << std::endl;
        else
            std::cout << "random lightcone "
                      << " using seed=" << _global_cli_dict._rng_seed << std::endl;
    }

    // Check for missing arguments.
    EXCEPT(argc >= 2, "Insufficient arguments. "
                      "see --help for usage.");

    // Dump information to the console if there is just
    // one rank.
#ifndef NLOG
    if (mpi::comm::world.size() == 1)
        LOG_PUSH(new hpc::log::stdout(hpc::log::info));
#endif

    // Cache the filenames.

    {
        char buf[501];
        readlink("/proc/self/exe", buf, 500);
    }

    // First thing, dump the version number. (Keep it simple)
    // LOGILN("TAO version: ", TOSTRING(VERSION));
    LOGILN("TAO version: 5.0");
}

bool KdApplication::validate(const cli_dict& _global_cli_dict)
{
    // Validate the command line arguments.
    // decimation validation
    if (_global_cli_dict._decmin < -90.0 || _global_cli_dict._decmin > 90.0)
    {
        std::cerr << "Minimum DEC cannot be less than -90 or greater than 90 degrees" << std::endl;
        return false;
    }
    if (_global_cli_dict._decmax < -90.0 || _global_cli_dict._decmax > 90.0)
    {
        std::cerr << "Maximum DEC cannot be less than -90 or greater than 90 degrees" << std::endl;
        return false;
    }
    if (_global_cli_dict._decmin >= _global_cli_dict._decmax)
    {
        std::cerr << "Minimum DEC must be less than Maximum DEC" << std::endl;
        return false;
    }
    // right ascension validation
    if (_global_cli_dict._ramin < 0.0 || _global_cli_dict._ramin > 360.0)
    {
        std::cerr << "Minimum RA cannot be less than zero or greater than 360 degrees" << std::endl;
        return false;
    }
    if (_global_cli_dict._ramax < 0.0 || _global_cli_dict._ramax > 360.0)
    {
        std::cerr << "Maximum RA cannot be less than zero or greater than 360 degrees" << std::endl;
        return false;
    }
    if (_global_cli_dict._ramin >= _global_cli_dict._ramax)
    {
        std::cerr << "Minimum RA must be less than Maximum RA" << std::endl;
        return false;
    }

    // redshift validation
    if (_global_cli_dict._zmin < 0.0)
    {
        std::cerr << "Minimum redshift must be greater than or equal to zero." << std::endl;
        return false;
    }
    if (_global_cli_dict._zmax <= _global_cli_dict._zmin)
    {
        std::cerr << "Minimum redshift must be less than Maximum redshift." << std::endl;
        return false;
    }

    if (_global_cli_dict._outdir.empty())
    {
        std::cerr << "Error: --outdir cannot be empty." << std::endl;
        return false;
    }

    if (_global_cli_dict._outfile.empty())
    {
        std::cerr << "Error: --outfile cannot be empty." << std::endl;
        return false;
    }

    if (_global_cli_dict._filter_field.empty())
    {
        std::cerr << "Error: --filterfield cannot be empty." << std::endl;
        return false;
    }

    if (_global_cli_dict._filter_min.empty() && _global_cli_dict._filter_max.empty())
    {
        std::cerr << "Error: --filtermin and --filtermax cannot both be empty." << std::endl;
        return false;
    }
    if (_global_cli_dict._filter_min == _global_cli_dict._filter_max)
    {
        std::cerr << "Error: --filtermin and --filtermax cannot be the same." << std::endl;
        return false;
    }

    // When filtermax is set, both values must parse as numbers and filtermax must exceed filtermin
    if (!_global_cli_dict._filter_max.empty())
    {
        try
        {
            double fmin = std::stod(_global_cli_dict._filter_min);
            double fmax = std::stod(_global_cli_dict._filter_max);
            if (fmax <= fmin)
            {
                std::cerr << "Error: --filtermax must be greater than --filtermin." << std::endl;
                return false;
            }
        }
        catch (const std::invalid_argument&)
        {
            std::cerr << "Error: --filtermin and --filtermax must be numeric values." << std::endl;
            return false;
        }
        catch (const std::out_of_range&)
        {
            std::cerr << "Error: --filtermin or --filtermax value is out of range." << std::endl;
            return false;
        }
    }

    return true;
}

void KdApplication::operator()()
{
    // Preprocess the incoming XML file, only if we're
    // the root process, as we don't want any conflicts
    // in writing the processed file.
    mpi::comm::world.bcast(_currentxml_version);
    LOGDLN("Current XMl Schema Version: ", _currentxml_version);

    // Load all the modules first up.
    load_and_connect_modules();

    // Prepare the logging.
    string subjobindex = "0";
    setup_log(_global_cli_dict._outdir + "/log/tao.log." + subjobindex);
    LOG_PUSH(
        new mpi::logger(_global_cli_dict._outdir + "/log/debug.log." + subjobindex, log::debug));

    // Check if there are any checkpoint files available for use.
    boost::optional<boost::property_tree::ptree> checkpoint;

    // Initialise all the modules.
    try
    {
        for (auto module : _fact)
            module->initialise(_global_cli_dict, checkpoint);
    }
    catch (const std::runtime_error& e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        exit(1);
    }

    // Mark the beginning of the run.
    LOG_PUSH_TAG("progress");
    LOGILN(runtime(), ",start");
    LOG_POP_TAG("progress");

    // Run.
    _execute();

    // Finalise all the modules.
    for (auto module : _fact)
        module->finalise();

    // Mark the conclusion of the run.
    mpi::comm::world.barrier();
    LOG_PUSH_TAG("progress");
    LOGILN(runtime(), ",end,successful");
    LOG_POP_TAG("progress");

    // Dump timing information to the end of the info file.
    LOGILN("Module metrics:", setindent(2));
    for (auto module : _fact)
        module->log_metrics();

    if (mpi::comm::world.rank() == 0)
    {
        if (mpi::comm::world.size() > 1)
        {
            // If we are the root rank, we need to gather results from all ranks.
            tao::modules::hdf5<backends::kdtree_backend>::process_cli_options(_global_cli_dict);
        }
        // If we are the root rank, we need to gather results from all ranks.
        // Close the log file.
        LOG_POP();
    }
}

///
/// Load modules.
///
void KdApplication::load_and_connect_modules()
{
    LOGILN("Loading modules from file: ", setindent(2));

    // Register all the available science modules.
    tao::register_modules(_fact);
    // Register the modules from the command line.
    LOGILN("Registering modules from command line. Hardcoded to be lightcone and "
           "hdf5.");
    _fact.create_module("light-cone", "lightcone");
    _fact.create_module("hdf5", "hdf5");
    _fact["hdf5"]->add_parent(*_fact["lightcone"]);
    LOGILN("Done.", setindent(-2));
}

///
/// Prepare log file.
///
void KdApplication::setup_log(const string& filename)
{
    if (mpi::comm::world.rank() == 0)
    {
#ifndef NLOG
        // Create the log directory if it does not exist.
        if (!exists(_global_cli_dict._outdir + "/log"))
        {
            create_directories(_global_cli_dict._outdir + "/log");
        }
        LOGDLN("Setting logging file to: ", filename);
        log::file* Logf = new log::file(filename, log::info);
        Logf->add_tag("progress");
        LOG_PUSH(Logf);
#endif
    }
}

///
/// Execute the application.
///
void KdApplication::_execute()
{
    // Keep looping over modules until all report being complete.
    bool complete;
    unsigned long long it = 1;

    do
    {
        LOGDLN("Beginning iteration: ", it, setindent(2));

        // Reset the complete flag.
        complete = true;

        // Loop over the modules.
        for (auto module : _fact)
        {
            module->process(it);
            if (!module->complete())
                complete = false;
        }

        // Advance the counter.
        ++it;

        LOGDLN("Done.", setindent(-2));
    } while (!complete);
}

} // namespace tao
