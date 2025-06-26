#ifndef tao_base_utils_hh
#define tao_base_utils_hh

#include <boost/filesystem.hpp>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <libhpc/system/varray.hh>
#include <libhpc/numerics/constants.hh>
#include "types.hh"

namespace tao {

   class lightcone;

   template< class T >
   struct cosmology
   {
      T hubble;
      T omega_m;
      T omega_l;
      T omega_k;
      T omega_r;
   };

   template< class T >
   auto
   expansion_to_redshift( T ef ) -> T
   {
      return 1.0/ef - 1.0;
   }

   template< class T >
   T
   redshift_to_age_func( double x,
                         void* param )
   {
      cosmology<T>* cosmo = (cosmology<T>*)param;
      return 1.0/sqrt( cosmo->omega_k +
                       cosmo->omega_m/x +
                       cosmo->omega_r/(x*x) +
                       cosmo->omega_l*x*x );
   }

   template< class T >
   T
   redshift_to_age( T redshift,
                    T hubble = 73.0,
                    T omega_m = 0.25,
                    T omega_l = 0.75 )
   {
      // Prepare the GSL function.
      cosmology<T> cosmo;
      cosmo.hubble = hubble;
      cosmo.omega_m = omega_m;
      cosmo.omega_l = omega_l;
      cosmo.omega_r = 4.165e-5/((hubble/100.0)*(hubble/100.0));
      cosmo.omega_k = 1.0 - cosmo.omega_m - cosmo.omega_l - cosmo.omega_r;
      gsl_function func;
      func.function = &redshift_to_age_func<T>;
      func.params = &cosmo;

      // Create some work space.
      gsl_integration_workspace* work = gsl_integration_workspace_alloc( 1000 );

      // Perform the integration.
      T res, abserr;
      T upp = 1.0/(1.0 + redshift);
      gsl_integration_qag( &func, 0.0, upp, 1e-5, 1e-8, 1000, GSL_INTEG_GAUSS21, work, &res, &abserr );
      res *= (977.8/hubble); // convert to Gyrs

      // Don't forget to free work space.
      gsl_integration_workspace_free( work );

      return res;
   }

   real_type
   observed_redshift( real_type z,
                      hpc::varray<real_type,3> const& pos,
                      hpc::varray<real_type,3> const& vel,
                      real_type h0,
                      real_type c = hpc::constant::c_km_s );

   real_type
   approx_observed_redshift( lightcone const& lc,
                             hpc::varray<real_type,3> const& pos,
                             hpc::varray<real_type,3> const& vel );

   boost::filesystem::path
   data_prefix();

   template< class BoxSeqT,
             class T>
   int GetIntersection( T fDst1, T fDst2,
		                BoxSeqT const&  pb1,
		                BoxSeqT const&  pb2,
		                BoxSeqT &  Hit) {

	  if ( (fDst1 * fDst2) >= 0.0f) return 0;
      if ( fDst1 == fDst2) return 0;
      Hit[0] = pb1[0] + (pb2[0]-pb1[0]) * ( -fDst1/(fDst2-fDst1) );
      Hit[1] = pb1[1] + (pb2[1]-pb1[1]) * ( -fDst1/(fDst2-fDst1) );
      Hit[2] = pb1[2] + (pb2[2]-pb1[2]) * ( -fDst1/(fDst2-fDst1) );
      return 1;
   }

   template< class BoxSeqT>
   int InBox( BoxSeqT &  Hit,
		             BoxSeqT const&  Box1,
					 BoxSeqT const&  Box2,
					 const int Axis) {

      if ( Axis==1 && Hit[2] >= Box1[2] && Hit[2] <= Box2[2] && Hit[1] >= Box1[1] && Hit[1] <= Box2[1])
    	  return 1;
      if ( Axis==2 && Hit[2] >= Box1[2] && Hit[2] <= Box2[2] && Hit[0] >= Box1[0] && Hit[0] <= Box2[0])
    	  return 1;
      if ( Axis==3 && Hit[0] >= Box1[0] && Hit[0] <= Box2[0] && Hit[1] >= Box1[1] && Hit[1] <= Box2[1])
    	  return 1;
      return 0;
   }

   // returns true if line (L1, L2) intersects with the box (B1, B2)
   // returns intersection point in Hit
   template< class BoxSeqT,
	         class LineSeqT >
   int checklinebox( BoxSeqT const&  B1,
		             BoxSeqT const&  B2,
					 LineSeqT const & L1,
					 LineSeqT const & L2,
					 LineSeqT &Hit )
   {
      typedef typename BoxSeqT::value_type real_type;

      if (L2[0] < B1[0] && L1[0] < B1[0]) return false;
      if (L2[0] > B2[0] && L1[0] > B2[0]) return false;
      if (L2[1] < B1[1] && L1[1] < B1[1]) return false;
      if (L2[1] > B2[1] && L1[1] > B2[1]) return false;
      if (L2[2] < B1[2] && L1[2] < B1[2]) return false;
      if (L2[2] > B2[2] && L1[2] > B2[2]) return false;
      if (L1[0] >= B1[0] && L1[0] <= B2[0] &&
         L1[1] >= B1[1] && L1[1] <= B2[1] &&
         L1[2] >= B1[2] && L1[2] <= B2[2])
         {
    	    Hit = L1;
            return true;
         }
      if (
         (GetIntersection<BoxSeqT, real_type>( L1[0]-B1[0], L2[0]-B1[0], L1, L2, Hit) && InBox( Hit, B1, B2, 1 ))
      || (GetIntersection<BoxSeqT, real_type>( L1[1]-B1[1], L2[1]-B1[1], L1, L2, Hit) && InBox( Hit, B1, B2, 2 ))
      || (GetIntersection<BoxSeqT, real_type>( L1[2]-B1[2], L2[2]-B1[2], L1, L2, Hit) && InBox( Hit, B1, B2, 3 ))
      || (GetIntersection<BoxSeqT, real_type>( L1[0]-B2[0], L2[0]-B2[0], L1, L2, Hit) && InBox( Hit, B1, B2, 1 ))
      || (GetIntersection<BoxSeqT, real_type>( L1[1]-B2[1], L2[1]-B2[1], L1, L2, Hit) && InBox( Hit, B1, B2, 2 ))
      || (GetIntersection<BoxSeqT, real_type>( L1[2]-B2[2], L2[2]-B2[2], L1, L2, Hit) && InBox( Hit, B1, B2, 3 )))
   	     return true;

      return false;
   }
}

#endif
