#ifndef tao_base_globals_hh
#define tao_base_globals_hh

#include <libhpc/system/timer.hh>
#include "simulation.hh"

namespace tao {

   typedef double real_type;

   hpc::real_timer::time_type
   runtime();

   // Simulations.
   extern simulation millennium;
   extern simulation mini_millennium;

}

#endif
