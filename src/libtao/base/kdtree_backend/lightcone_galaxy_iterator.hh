#ifndef tao_base_kdtree_backend_lightcone_galaxy_iterator_hh
#define tao_base_kdtree_backend_lightcone_galaxy_iterator_hh

#include "../lightcone_galaxy_iterator.hh"
#include "kdtree_backend.hh"

namespace tao {

   template< class BackendT >
   class kdtree_lightcone_galaxy_iterator
   {
   public:

      typedef          BackendT                 backend_type;
      typedef typename lightcone::tile_iterator tile_iterator;
      typedef typename backend_type::tile_galaxy_iterator tile_galaxy_iterator;
      typedef          batch<real_type>& value_type;
      typedef          value_type        reference_type;

   public:

      kdtree_lightcone_galaxy_iterator()
      {
      }

      kdtree_lightcone_galaxy_iterator( lightcone const& lc,
                                        backend_type& be,
                                        query<real_type>& qry,
                                        batch<real_type>* bat = 0,
                                        filter const* flt = 0 )
         : _lc( &lc ),
           _be( &be ),
           _qry( &qry ),
           _bat( bat ),
           _flt( flt ),
           _n_processed_gals(0)
      {
         _snap_ii = 0;
         _snap_lc = new lightcone();
         _snap_lc->copy(lc);
         _tile_it_common_to_all_snapshots = _snap_lc->tile_begin();
         // The rdb_vs_kd needs the do loop gone
         //do
         {
             lightcone* mylc = new lightcone();
             // mylc = _snap_lc;
             mylc->copy(*_snap_lc);
             mylc->set_geometry( hpc::to_degrees( lc.min_ra() ), hpc::to_degrees( lc.max_ra() ),
                                    hpc::to_degrees( lc.min_dec() ), hpc::to_degrees( lc.max_dec() ),
                                    lc.redshift_bins()[_snap_ii], lc.redshift_bins()[_snap_ii + 1] );
             //_snap_lc = mylc;
            // Keep it broader than it needs to be to allow comparison with RDB (Also will need a separate instance of _snap_lc for each snapshot eventually)
            //_snap_lc->set_geometry( hpc::to_degrees( lc.min_ra() ), hpc::to_degrees( lc.max_ra() ),
            //                       hpc::to_degrees( lc.min_dec() ), hpc::to_degrees( lc.max_dec() ),
            //                       lc.redshift_bins()[_snap_ii], lc.redshift_bins()[_snap_ii + 1] );
            be.load_snapshot( lc.snapshot_bins()[_snap_ii] );
            // TODO: leaking code
            _gal_it = new lightcone_galaxy_iterator<backend_type>( *mylc, be, qry, _tile_it_common_to_all_snapshots, bat, flt );
            _gal_end = new lightcone_galaxy_iterator<backend_type>();
            //if( *_gal_it == *_gal_end )
            //   ++_snap_ii;
         }
         //while( _snap_ii < lc.snapshot_bins().size() && *_gal_it == *_gal_end );

      }

      reference_type
      operator*()
      {
         return *_gal_it;
      }

      void
      operator++()
      {
         if( ++(*_gal_it) == *_gal_end )
         {
            if( ++_snap_ii < _lc->snapshot_bins().size() )
            {
                _lc->rng_reset();
                _tile_it_common_to_all_snapshots = _snap_lc->tile_begin();

                lightcone* mylc = new lightcone();
                mylc->copy(*_snap_lc);
                mylc->set_geometry( hpc::to_degrees( _lc->min_ra() ), hpc::to_degrees( _lc->max_ra() ),
                                    hpc::to_degrees( _lc->min_dec() ), hpc::to_degrees( _lc->max_dec() ),
                                    _lc->redshift_bins()[_snap_ii], _lc->redshift_bins()[_snap_ii + 1] );
                // Keep it broader than it needs to be to allow comparison with RDB (Also will need a separate instance of _snap_lc for each snapshot eventually)
               //_snap_lc->set_geometry( hpc::to_degrees( _lc->min_ra() ), hpc::to_degrees( _lc->max_ra() ),
               //                       hpc::to_degrees( _lc->min_dec() ), hpc::to_degrees( _lc->max_dec() ),
                //                      _lc->redshift_bins()[_snap_ii], _lc->redshift_bins()[_snap_ii + 1] );
               _be->load_snapshot( _lc->snapshot_bins()[_snap_ii] );
               // TODO: leaking code
               _gal_it = new lightcone_galaxy_iterator<backend_type>( *mylc, *_be, *_qry, _tile_it_common_to_all_snapshots, _bat, _flt );
               _gal_end = new lightcone_galaxy_iterator<backend_type>();
            }
         }
      }

      bool
      operator==( kdtree_lightcone_galaxy_iterator const& op ) const
      {
         return _snap_ii >= _lc->snapshot_bins().size();
      }

      bool
      operator!=( kdtree_lightcone_galaxy_iterator const& op ) const
      {
         return _snap_ii < _lc->snapshot_bins().size();
      }

      unsigned long long
      n_processed_gals() const
      {
         return _gal_it->n_processed_gals();
      }

      unsigned
      tile_index() const
      {
         return _gal_it->tile_index();
      }

      bool
      should_checkpoint()
      {
         return _gal_it->should_checkpoint();
      }

      void
      save_checkpoint( boost::property_tree::ptree& pt )
      {
         _gal_it->save_checkpoint( pt );
      }
      

   protected:

      tao::lightcone const* _lc;
      backend_type* _be;
      query<real_type>* _qry;
      tao::batch<real_type>* _bat;
      filter const* _flt;
      unsigned _snap_ii;
      lightcone* _snap_lc;
      tile_iterator _tile_it_common_to_all_snapshots;
      lightcone_galaxy_iterator<backend_type> *_gal_it, *_gal_end;
      unsigned long long _n_processed_gals;
   };

}

#endif
