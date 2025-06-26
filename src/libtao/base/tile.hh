#ifndef tao_base_tile_hh
#define tao_base_tile_hh

#include <array>
#include <libhpc/numerics/coords.hh>
#include "box.hh"
#include "query.hh"
#include "lightcone.hh"

namespace tao {
   using namespace hpc;

   ///
   ///
   ///
   template< class T >
   class tile
      : public box<T>
   {
   public:

      typedef T real_type;

   public:

      tile( const tao::lightcone* lc = 0,
            std::array<real_type,3> offs = { { 0, 0, 0 } } )
         : box<T>( 0 ),
           _lc( 0 )
      {
	      set_lightcone( lc, offs );
      }

      void
      set_lightcone( const lightcone* lc,
                     std::array<real_type,3> offs = { { 0, 0, 0 } } )
      {
         if( lc )
         {
            _lc = lc;
            this->set_simulation( lc->simulation() );
            set_offset( offs );
            this->set_random( _lc->random(), _lc->rng_engine() );
            this->set_origin( _lc->origin() );
         }
         else
            this->set_simulation( 0 );
      }

      void
      set_offset( const std::array<real_type,3>& offs )
      {
         this->_min = offs;
         for( unsigned ii = 0; ii < 3; ++ii )
            this->_max[ii] = this->_min[ii] + _lc->simulation()->box_size();
      }

      const tao::lightcone*
      lightcone() const
      {
         return _lc;
      }

      template< class Backend >
      typename Backend::tile_galaxy_iterator
      galaxy_begin( tao::query<real_type>& query,
                    Backend& be ) const
      {
         return be.galaxy_begin( query, *this );
      }

      template< class Backend >
      typename Backend::tile_galaxy_iterator
      galaxy_end( tao::query<real_type>& query,
                  Backend& be ) const
      {
         return be.galaxy_end( query, *this );
      }

      bool
      contains( std::array<real_type,3> const& crd,
                unsigned snap ) const
      {
         // TODO: Include rotation and shifting.
         std::array<real_type,3> ecs;
         hpc::num::cartesian_to_ecs<real_type>(
            crd[0] + this->_min[0], // include tile offset
            crd[1] + this->_min[1], // include tile offset
            crd[2] + this->_min[2], // include tile offset
            ecs[0], ecs[1], ecs[2] );
         return ecs[0] >= _lc->min_ra()         && ecs[0] < _lc->max_ra() &&
                ecs[1] >= _lc->min_dec()        && ecs[1] < _lc->max_dec() &&
                ecs[2] >= _lc->min_dist( snap ) && ecs[2] < _lc->max_dist( snap );
      }

      template< class Iter >
      bool
      box_collision( Iter bnds_begin,
                     unsigned snap ) const
      {
          real_type pmin[3],pmax[3];
          real_type simsize=this->_sim->box_size();
         // TODO: Check rotation and shifting. (RS cannot ignore the wrapping)
         std::array<real_type,3> min, max, ecs_min, ecs_max;

         // The cone
          ecs_min[0] = _lc->min_ra();         ecs_max[0] = _lc->max_ra();
          ecs_min[1] = _lc->min_dec();        ecs_max[1] = _lc->max_dec();
          ecs_min[2] = _lc->min_dist( snap ); ecs_max[2] = _lc->max_dist( snap );

         // Apply the rotation and translation
         for( unsigned ii = 0; ii < 3; ++ii, ++bnds_begin )
         {
            pmin[ii] = (*bnds_begin)[0];
            pmax[ii] = (*bnds_begin)[1];
         }
         for( unsigned ii = 0; ii < 3; ++ii)
         {
             int irotation=this->_rot[ii];

             min[ii] = pmin[irotation]+this->translation()[irotation];
             max[ii] = pmax[irotation]+this->translation()[irotation];
         }
         // Because of wrapping need to check up to 8 rectangular boxes
         // Apply wrap
         for (unsigned ii = 0; ii < 3; ++ii)
         {
             if (min[ii]>simsize)
                 min[ii] -= simsize;
             if (max[ii]>simsize)
                 max[ii] -= simsize;
         }

         // Test if we need to split
          int ix=0;
          int iy=1;
          int iz=2;
          real_type xorigin=(this->_min[0] - this->_orig[0]);
          real_type yorigin=(this->_min[1] - this->_orig[1]);
          real_type zorigin=(this->_min[2] - this->_orig[2]);
         bool split_in_x=max[ix]<=min[ix];
         bool split_in_y=max[iy]<=min[iy];
         bool split_in_z=max[iz]<=min[iz];
         min[ix] += xorigin;
         min[iy] += yorigin;
         min[iz] += zorigin;

         max[ix] += xorigin;
         max[iy] += yorigin;
         max[iz] += zorigin;
         if (split_in_x)
         {
             std::array<real_type,3> xminL, xmaxL, xminH, xmaxH;
             xminL[ix] = min[ix];
             xmaxL[ix] = xorigin+simsize;
             xminH[ix] = xorigin;
             xmaxH[ix] = max[ix];
             if (split_in_y)
             {
                 std::array<real_type,3> yminL, ymaxL, yminH, ymaxH;
                 yminL[iy] = min[iy];
                 ymaxL[iy] = yorigin+simsize;
                 yminH[iy] = yorigin;
                 ymaxH[iy] = max[iy];
                 if (split_in_z)
                 {
                     // Split in x,y,z: 8 rectangle boxes to check
                     std::array<real_type,3> zminL, zmaxL, zminH, zmaxH;
                     zminL[iz] = min[iz];
                     zmaxL[iz] = zorigin+simsize;
                     zminH[iz] = zorigin;
                     zmaxH[iz] = max[iz];

                     //1.
                     min[ix] = xminL[ix];
                     min[iy] = yminL[iy];
                     min[iz] = zminL[iz];
                     max[ix] = xmaxL[ix];
                     max[iy] = ymaxL[iy];
                     max[iz] = zmaxL[iz];
                     if (ecs_box_collision( ecs_min, ecs_max, min, max ))
                         return true;
                     //2.
                     min[ix] = xminH[ix];
                     min[iy] = yminL[iy];
                     min[iz] = zminL[iz];
                     max[ix] = xmaxH[ix];
                     max[iy] = ymaxL[iy];
                     max[iz] = zmaxL[iz];
                     if (ecs_box_collision( ecs_min, ecs_max, min, max ))
                         return true;
                     //3.
                     min[ix] = xminL[ix];
                     min[iy] = yminH[iy];
                     min[iz] = zminL[iz];
                     max[ix] = xmaxL[ix];
                     max[iy] = ymaxH[iy];
                     max[iz] = zmaxL[iz];
                     if (ecs_box_collision( ecs_min, ecs_max, min, max ))
                         return true;
                     //4.
                     min[ix] = xminH[ix];
                     min[iy] = yminH[iy];
                     min[iz] = zminL[iz];
                     max[ix] = xmaxH[ix];
                     max[iy] = ymaxH[iy];
                     max[iz] = zmaxL[iz];
                     if (ecs_box_collision( ecs_min, ecs_max, min, max ))
                         return true;
                     //5.
                     min[ix] = xminL[ix];
                     min[iy] = yminL[iy];
                     min[iz] = zminH[iz];
                     max[ix] = xmaxL[ix];
                     max[iy] = ymaxL[iy];
                     max[iz] = zmaxH[iz];
                     if (ecs_box_collision( ecs_min, ecs_max, min, max ))
                         return true;
                     //6.
                     min[ix] = xminH[ix];
                     min[iy] = yminL[iy];
                     min[iz] = zminH[iz];
                     max[ix] = xmaxH[ix];
                     max[iy] = ymaxL[iy];
                     max[iz] = zmaxH[iz];
                     if (ecs_box_collision( ecs_min, ecs_max, min, max ))
                         return true;
                     //7.
                     min[ix] = xminL[ix];
                     min[iy] = yminH[iy];
                     min[iz] = zminH[iz];
                     max[ix] = xmaxL[ix];
                     max[iy] = ymaxH[iy];
                     max[iz] = zmaxH[iz];
                     if (ecs_box_collision( ecs_min, ecs_max, min, max ))
                         return true;
                     //8.
                     min[ix] = xminH[ix];
                     min[iy] = yminH[iy];
                     min[iz] = zminH[iz];
                     max[ix] = xmaxH[ix];
                     max[iy] = ymaxH[iy];
                     max[iz] = zmaxH[iz];
                     if (ecs_box_collision( ecs_min, ecs_max, min, max ))
                         return true;
                 } else {
                     // Split in x,y: Not split in z: 4 rectangle boxes to check

                     //1.
                     min[ix] = xminL[ix];
                     min[iy] = yminL[iy];
                     max[ix] = xmaxL[ix];
                     max[iy] = ymaxL[iy];
                     if (ecs_box_collision( ecs_min, ecs_max, min, max ))
                         return true;
                     //2.
                     min[ix] = xminH[ix];
                     min[iy] = yminL[iy];
                     max[ix] = xmaxH[ix];
                     max[iy] = ymaxL[iy];
                     if (ecs_box_collision( ecs_min, ecs_max, min, max ))
                         return true;
                     //3.
                     min[ix] = xminL[ix];
                     min[iy] = yminH[iy];
                     max[ix] = xmaxL[ix];
                     max[iy] = ymaxH[iy];
                     if (ecs_box_collision( ecs_min, ecs_max, min, max ))
                         return true;
                     //4.
                     min[ix] = xminH[ix];
                     min[iy] = yminH[iy];
                     max[ix] = xmaxH[ix];
                     max[iy] = ymaxH[iy];
                     if (ecs_box_collision( ecs_min, ecs_max, min, max ))
                         return true;
                 }
             } else {
                 if (split_in_z) {
                     // Split in x,z: Not splitin y: 4 rectangle boxes to check
                     std::array<real_type, 3> zminL, zmaxL, zminH, zmaxH;
                     zminL[iz] = min[iz];
                     zmaxL[iz] = zorigin+simsize;
                     zminH[iz] = zorigin;
                     zmaxH[iz] = max[iz];

                     //1.
                     min[ix] = xminL[ix];
                     min[iz] = zminL[iz];
                     max[ix] = xmaxL[ix];
                     max[iz] = zmaxL[iz];
                     if (ecs_box_collision(ecs_min, ecs_max, min, max))
                         return true;
                     //2.
                     min[ix] = xminH[ix];
                     min[iz] = zminL[iz];
                     max[ix] = xmaxH[ix];
                     max[iz] = zmaxL[iz];
                     if (ecs_box_collision(ecs_min, ecs_max, min, max))
                         return true;
                     //3.
                     min[ix] = xminL[ix];
                     min[iz] = zminH[iz];
                     max[ix] = xmaxL[ix];
                     max[iz] = zmaxH[iz];
                     if (ecs_box_collision(ecs_min, ecs_max, min, max))
                         return true;
                     //4.
                     min[ix] = xminH[ix];
                     min[iz] = zminH[iz];
                     max[ix] = xmaxH[ix];
                     max[iz] = zmaxH[iz];
                     if (ecs_box_collision(ecs_min, ecs_max, min, max))
                         return true;
                 } else {
                     // Split in x: Not split in y,z: 2 rectangle boxes to check
                     //1.
                     min[ix] = xminL[ix];
                     max[ix] = xmaxL[ix];
                     if (ecs_box_collision(ecs_min, ecs_max, min, max))
                         return true;
                     //2.
                     min[ix] = xminH[ix];
                     max[ix] = xmaxH[ix];
                     if (ecs_box_collision(ecs_min, ecs_max, min, max))
                         return true;
                 }
             }
         } else if (split_in_y) {
             std::array<real_type,3> yminL, ymaxL, yminH, ymaxH;
             yminL[iy] = min[iy];
             ymaxL[iy] = yorigin+simsize;
             yminH[iy] = yorigin;
             ymaxH[iy] = max[iy];
             if (split_in_z)
             {
                 // Split in y,z: Not split in x: 4 rectangle boxes to check
                 std::array<real_type,3> zminL, zmaxL, zminH, zmaxH;
                 zminL[iz] = min[iz];
                 zmaxL[iz] = zorigin+simsize;
                 zminH[iz] = zorigin;
                 zmaxH[iz] = max[iz];

                 //1.
                 min[iy] = yminL[iy];
                 min[iz] = zminL[iz];
                 max[iy] = ymaxL[iy];
                 max[iz] = zmaxL[iz];
                 if (ecs_box_collision( ecs_min, ecs_max, min, max ))
                     return true;
                 //2.
                 min[iy] = yminH[iy];
                 min[iz] = zminL[iz];
                 max[iy] = ymaxH[iy];
                 max[iz] = zmaxL[iz];
                 if (ecs_box_collision( ecs_min, ecs_max, min, max ))
                     return true;
                 //3.
                 min[iy] = yminL[iy];
                 min[iz] = zminH[iz];
                 max[iy] = ymaxL[iy];
                 max[iz] = zmaxH[iz];
                 if (ecs_box_collision( ecs_min, ecs_max, min, max ))
                     return true;
                 //4.
                 min[iy] = yminH[iy];
                 min[iz] = zminH[iz];
                 max[iy] = ymaxH[iy];
                 max[iz] = zmaxH[iz];
                 if (ecs_box_collision( ecs_min, ecs_max, min, max ))
                     return true;
             } else {
                 // Split in y: Not split in x,z: 2 rectangle boxes to check
                 //1.
                 min[iy] = yminL[iy];
                 max[iy] = ymaxL[iy];
                 if (ecs_box_collision( ecs_min, ecs_max, min, max ))
                     return true;
                 //2.
                 min[iy] = yminH[iy];
                 max[iy] = ymaxH[iy];
                 if (ecs_box_collision( ecs_min, ecs_max, min, max ))
                     return true;
             }
         } else if (split_in_z) {
            // Split in z: Not split in x,y: 2 rectangle boxes to check
             std::array<real_type,3> zminL, zmaxL, zminH, zmaxH;
             zminL[iz] = min[iz];
             zmaxL[iz] = zorigin+simsize;
             zminH[iz] = zorigin;
             zmaxH[iz] = max[iz];
             //1.
             min[iz] = zminL[iz];
             max[iz] = zmaxL[iz];
             if (ecs_box_collision( ecs_min, ecs_max, min, max ))
                 return true;
             //2.
             min[iz] = zminH[iz];
             max[iz] = zmaxH[iz];
             if (ecs_box_collision( ecs_min, ecs_max, min, max ))
                 return true;
         } else {
             // only one rectangle box to check
             return ecs_box_collision( ecs_min, ecs_max, min, max );
         }
         return false;
      }

   protected:

      tao::lightcone const* _lc;
   };

}

#endif
