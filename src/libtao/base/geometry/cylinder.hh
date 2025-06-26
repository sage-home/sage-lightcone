#ifndef tao_base_geometry_generators_hh
#define tao_base_geometry_generators_hh

#include <array>
#include <boost/iterator/iterator_facade.hpp>
#include "../types.hh"

namespace tao {

   /// Generate a cylinder based on ECS coordinates. The output is
   /// in cartesian coordinates.
   ///
   template< class T = tao::real_type >
   class cylinder_iterator
      : public boost::iterator_facade< cylinder_iterator<T>,
                                       std::array<T,3>,
                                       std::forward_iterator_tag,
                                       std::array<T,3> >
   {
      friend class boost::iterator_core_access;

   public:

      typedef T real_type;
      typedef std::array<T,3> value_type;
      typedef value_type reference_type;

   public:

      cylinder_iterator()
         : _done( true )
      {
      }

      cylinder_iterator( real_type ra,
                         real_type dec,
                         real_type offset,
                         real_type length,
                         real_type radius,
                         unsigned n_segs,
                         unsigned n_rings,
                         unsigned n_slices )
         : _ra( ra ),
           _dec( -dec ),
           _offs( offset ),
           _len( length ),
           _r( radius ),
           _n_segs( n_segs ),
           _n_rings( n_rings ),
           _n_slices( n_slices ),
           _crd{ 0, 0, 0 },
           _done( false )
      {
         _seps[0] = 2.0*M_PI/(real_type)n_segs;
         _seps[1] = radius/(real_type)n_rings;
         _seps[2] = length/(real_type)n_slices;
      }

      // lightcone_galaxy_iterator&
      // operator=( const lightcone_galaxy_iterator& op )
      // {
      // }

      // reference_type
      // operator*()
      // {
      // }

   protected:

      void
      increment()
      {
         if( ++_crd[0] >= _n_segs ) {
            _crd[0] = 0;
            if( ++_crd[1] >= _n_rings ) {
               _crd[1] = 0;
               if( ++_crd[2] >= _n_slices + 1 )
                  _done = true;
            }
         }
      }

      bool
      equal( cylinder_iterator const& op ) const
      {
         return _done == op._done;
      }

      reference_type
      dereference() const
      {
         // Calculate base point.
         real_type r = _seps[1]*(real_type)(_crd[1] + 1);
         real_type x = _offs + _seps[2]*(real_type)_crd[2];
         real_type y = r*cos( _seps[0]*(real_type)_crd[0] );
         real_type z = r*sin( _seps[0]*(real_type)_crd[0] );

         // Rotation around the y-axis.
         real_type xp = x*cos( _dec ) + z*sin( _dec );
         real_type yp = y;
         real_type zp = -x*sin( _dec ) + z*cos( _dec );

         // Rotation around the z-axis.
         x = xp*cos( _ra ) - yp*sin( _ra );
         y = xp*sin( _ra ) + yp*cos( _ra );
         z = zp;

         return value_type{ x, y, z };
      }

   protected:

      real_type _ra, _dec;
      real_type _offs, _len;
      real_type _r;
      unsigned _n_segs, _n_rings, _n_slices;
      std::array<real_type,3> _seps;
      std::array<unsigned,3> _crd;
      bool _done;
   };

}

#endif
