#ifndef tao_base_types_hh
#define tao_base_types_hh

#include <string>
#include <vector>

namespace tao
{

typedef double real_type;

struct cli_dict
{
    std::string _dataset;
    double _decmin;
    double _decmax;
    double _ramin;
    double _ramax;
    double _zmin;
    double _zmax;
    bool _unique;
    int _rng_seed;
    std::string _filter_field;
    std::string _filter_min;
    std::string _filter_max;
    std::string _outfile;
    std::vector<std::string> _output_fields;
    std::string _outdir;
    bool _central_galaxies = false;
};

} // namespace tao

#endif
