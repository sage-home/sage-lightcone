#ifndef tao_base_galaxy_iterator
#define tao_base_galaxy_iterator

#include "galaxy.hh"

namespace tao {

   template< class Source,
             class Backend = backends::postgresql<typename Source::real_type> >
   class galaxy_iterator
      : public boost::iterator_facade< galaxy_iterator<Source,Backend>,
                                       galaxy<typename Source::real_type>&,
				       std::forward_iterator_tag,
                                       galaxy<typename Source::real_type>& >
   {
      friend class boost::iterator_core_access;

   public:

      typedef Source source_type;
      typedef typename source_type::real_type real_type;
      typedef Backend backend_type;
      typedef galaxy<real_type> value_type;
      typedef value_type reference_type;

   public:

      galaxy_iterator( source_type& src,
                       backend_type& be )
         : _src( src ),
           _be( be )
      {
      }

   protected:

      void
      increment()
      {
         _be.fetch( _src, _gal );
      }

      bool
      equal( const galaxy_iterator& op ) const
      {
         return ;
      }

      reference_type
      dereference() const
      {
         return _gal;
      }

   protected:

      source_type& _src;
      backend_type& _be;
      galaxy<real_type> _gal;
   };

}

#endif
