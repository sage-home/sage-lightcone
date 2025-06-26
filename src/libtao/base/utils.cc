#include <libhpc/system/filesystem.hh>
#include "utils.hh"
#include "lightcone.hh"

namespace tao {

   real_type
   observed_redshift( real_type z,
                      hpc::varray<real_type,3> const& pos,
                      hpc::varray<real_type,3> const& vel,
                      real_type h0,
                      real_type c )
   {
      real_type dist = pos.magnitude();
      if( dist > 0.0 )
         return (1.0 + z)*(1.0 + (pos/dist).dot( vel )/c) - 1.0;
      else
         return 0.0;
   }

   real_type
   approx_observed_redshift( lightcone const& lc,
                             hpc::varray<real_type,3> const& pos,
                             hpc::varray<real_type,3> const& vel )
   {
      real_type dist = pos.magnitude();
      if( dist > 0.0 )
         return lc.distance_to_redshift( dist + (pos/dist).dot( vel )/lc.simulation()->hubble() );
      else
         return 0.0;
   }

   hpc::fs::path
   data_prefix()
   {
      return hpc::executable_path().parent_path().parent_path()/"data";
   }

}
