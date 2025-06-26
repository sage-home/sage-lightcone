#ifndef tao_base_integration_hh
#define tao_base_integration_hh

#include <libhpc/numerics/spline_integrator.hh>
#include "types.hh"

#define M_C 2.99792458e18 // angstrom/s

namespace tao {

   template< class Spline >
   real_type
   integrate( hpc::num::spline_integrator<real_type> const& integ,
              Spline const& sp )
   {
      // Do it in reverse to sum smaller values first.
      return integ(
         sp,
         []( double x, double val )
         {
            val = std::max( val, 0.0 );
            return val;
         },
	 []( double x )
	 {
	    return 1.0/x;
	 }
         );
   }

   template< class Spline >
   real_type
   integrate( const Spline& sp )
   {
      hpc::num::spline_integrator<real_type> integ;
      return -M_C*integrate( integ, sp );
   }

   template< class Spline0,
             class Spline1 >
   real_type
   integrate( const hpc::num::spline_spline_integrator<real_type>& integ,
              const Spline0& sp0,
              const Spline1& sp1 )
   {
      return integ(
         sp0, sp1,
         []( double x, double val0, double val1 )
         {
            val0 = std::max( val0, 0.0 );
            val1 = std::max( val1, 0.0 );
            return val0*val1;
         }
         );
   }

   template< class Spline0,
             class Spline1 >
   real_type
   integrate( const Spline0& sp0,
              const Spline1& sp1 )
   {
      hpc::num::spline_spline_integrator<real_type> integ;
      return integrate( integ, sp0, sp1 );
   }

}

#endif
