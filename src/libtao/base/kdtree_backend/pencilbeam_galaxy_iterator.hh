#ifndef tao_base_kdtree_backend_pencilbeam_galaxy_iterator_hh
#define tao_base_kdtree_backend_pencilbeam_galaxy_iterator_hh

#include "../pencilbeam_galaxy_iterator.hh"
#include "kdtree_backend.hh"

namespace tao {

   template< class BackendT >
   class kdtree_pencilbeam_galaxy_iterator
   {
   public:

      typedef          BackendT                           backend_type;
      typedef typename pencilbeam::pb_tile_iterator       tile_iterator;
      typedef typename backend_type::tile_galaxy_iterator tile_galaxy_iterator;
      typedef          batch<real_type>&                  value_type;
      typedef          value_type                         reference_type;

   public:

      kdtree_pencilbeam_galaxy_iterator()
      {
      }

      kdtree_pencilbeam_galaxy_iterator( pencilbeam const& pb,
                                         backend_type& be,
                                         query<real_type>& qry,
                                         batch<real_type>* bat = 0,
                                         filter const* flt = 0 )
         : _pb( &pb ),
           _be( &be ),
           _qry( &qry ),
           _bat( bat ),
           _flt( flt ),
           _snap_pb(pb)
      {
         _snap_ii = 0;
         do
         {
            _snap_pb.set_geometry( hpc::to_degrees( pb.ra_center() ), hpc::to_degrees( pb.dec_center() ),
                                   hpc::to_degrees( pb.pb_radius() ),
                                   pb.redshift_bins()[_snap_ii], pb.redshift_bins()[_snap_ii + 1] );
            be.load_snapshot( pb.snapshot_bins()[_snap_ii] );
            _gal_it = pencilbeam_galaxy_iterator<backend_type>( _snap_pb, be, qry, bat, flt );
            _gal_end = pencilbeam_galaxy_iterator<backend_type>();
            if( _gal_it == _gal_end )
               ++_snap_ii;
         }
         while( _snap_ii < pb.snapshot_bins().size() && _gal_it == _gal_end );
      }

      reference_type
      operator*()
      {
         return *_gal_it;
      }

      void
      operator++()
      {
         if( ++_gal_it == _gal_end )
         {
            if( ++_snap_ii < _pb->snapshot_bins().size() )
            {
               _snap_pb.set_geometry( hpc::to_degrees( _pb->ra_center() ), hpc::to_degrees( _pb->dec_center() ),
                                      hpc::to_degrees( _pb->pb_radius() ),
                                      _pb->redshift_bins()[_snap_ii], _pb->redshift_bins()[_snap_ii + 1] );
               _be->load_snapshot( _pb->snapshot_bins()[_snap_ii] );
               _gal_it = pencilbeam_galaxy_iterator<backend_type>( _snap_pb, *_be, *_qry, _bat, _flt );
               _gal_end = pencilbeam_galaxy_iterator<backend_type>();
            }
         }
      }

      bool
      operator==( kdtree_pencilbeam_galaxy_iterator const& op ) const
      {
         return _snap_ii >= _pb->snapshot_bins().size();
      }

      bool
      operator!=( kdtree_pencilbeam_galaxy_iterator const& op ) const
      {
         return _snap_ii < _pb->snapshot_bins().size();
      }

      unsigned long long
      n_processed_gals() const
      {
         return _gal_it.n_processed_gals();
      }

      unsigned
      tile_index() const
      {
         return _gal_it.tile_index();
      }

      bool
      should_checkpoint()
      {
         return _gal_it.should_checkpoint();
      }

      void
      save_checkpoint( boost::property_tree::ptree& pt )
      {
         _gal_it.save_checkpoint( pt );
      }

   protected:

      tao::pencilbeam const* _pb;
      backend_type* _be;
      query<real_type>* _qry;
      tao::batch<real_type>* _bat;
      filter const* _flt;
      unsigned _snap_ii;
      tao::pencilbeam _snap_pb;
      pencilbeam_galaxy_iterator<backend_type> _gal_it, _gal_end;
   };

}

#endif
