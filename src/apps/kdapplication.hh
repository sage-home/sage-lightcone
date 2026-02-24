#ifndef tao_apps_application_hh
#define tao_apps_application_hh

#include <libhpc/libhpc.hh>
#include <libhpc/mpi/application.hh>
#include <libtao/base/base.hh>

namespace tao
{
using namespace hpc;

///
///
///
class KdApplication : public mpi::application
{
public:
    KdApplication(int argc, char* argv[]);

    bool validate(const cli_dict& _global_cli_dict);

    void arguments(int argc, char* argv[]);

    void operator()();

protected:
    ///
    /// Load modules.
    ///
    void load_and_connect_modules();

    ///
    /// Prepare log file.
    ///
    void setup_log(const std::string& filename);

    ///
    /// Execute the application.
    ///
    void _execute();

protected:
    std::string _currentxml_version;

    cli_dict _global_cli_dict;

    double _decmin = -90.0;
    double _decmax = 90.0;
    double _ramin = 0.0;
    double _ramax = 360.0;
    double _zmin = 0.0;
    double _zmax = 10.0;
    factory<backends::kdtree_backend> _fact;
};
} // namespace tao

#endif
