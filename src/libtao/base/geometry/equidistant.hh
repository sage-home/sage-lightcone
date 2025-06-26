#ifndef tao_base_geometry_equidistant_hh
#define tao_base_geometry_equidistant_hh

#include <array>
#include <boost/iterator/iterator_facade.hpp>
#include "../types.hh"

namespace tao {

   /// Generate equidistance points in a box.
   ///
   template< class T = tao::real_type >
   class equidistant_iterator
      : public boost::iterator_facade< equidistant_iterator<T>,
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

      equidistant_iterator()
         : _done( true )
      {
      }

      equidistant_iterator( real_type size,
                            unsigned density )
         : _size( size ),
           _dens( density ),
           _sep( size / (real_type)density ),
           _crd{ 0, 0, 0 },
           _done( false )
      {
      }

   protected:

      void
      increment()
      {
         if( ++_crd[0] >= _dens )
         {
            _crd[0] = 0;
            if( ++_crd[1] >= _dens )
            {
               _crd[1] = 0;
               if( ++_crd[2] >= _dens )
                  _done = true;
            }
         }
      }

      bool
      equal( equidistant_iterator const& op ) const
      {
         return _done == op._done;
      }

      reference_type
      dereference() const
      {
         return value_type{
            0.5 * _sep + (_crd[0] * _sep),
            0.5 * _sep + (_crd[1] * _sep),
            0.5 * _sep + (_crd[2] * _sep)
         };
      }

   protected:

      real_type _size;
      unsigned _dens;
      real_type _sep;
      std::array<unsigned,3> _crd;
      bool _done;
   };

}

#endif
