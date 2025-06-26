#ifndef tao_base_lightcone_hh
#define tao_base_lightcone_hh

#include <vector>
#include <array>
#include <libhpc/system/random.hh>
#include <libhpc/system/view.hh>
#include <libhpc/numerics/interp.hh>
#include "simulation.hh"
#include "filter.hh"
#include "query.hh"

namespace tao {

   class lightcone_tile_iterator;

   /// Lightcone representation. Describes the geometry of a lightcone,
   /// including redshift distance mappings.
   ///
   class lightcone
   {
   public:

      typedef lightcone_tile_iterator tile_iterator;

   public:

      /// Constructor.
      /// @param[in] sim  Simulation for the lightcone.
      ///
      lightcone( tao::simulation const* sim = nullptr );
      /// Copier
      void copy( const lightcone& from);

      /// Set lightcone simulation.
      /// @param[in] sim  Simulation for the lightcone.
      ///
      void
      set_simulation( tao::simulation const* sim );

      /// Get lightcone simulation.
      /// @returns Simulation.
      ///
      tao::simulation const*
      simulation() const;

      /// Set lightcone geometry.
      /// @param[in] ra_min  Minimum right-ascension.
      /// @param[in] ra_max  Maximum right-ascension.
      /// @param[in] dec_min  Minimum declination.
      /// @param[in] dec_max  Maximum declination.
      /// @param[in] z_min  Minimum redshift.
      /// @param[in] z_max  Maximum redshift.
      ///
      void
      set_geometry( real_type ra_min,
                    real_type ra_max,
                    real_type dec_min,
                    real_type dec_max,
                    real_type z_max,
                    real_type z_min = 0 );

      void
      set_min_ra( real_type ra_min );

      void
      set_max_ra( real_type ra_max );

      void
      set_min_dec( real_type dec_min );

      void
      set_max_dec( real_type dec_max );

      void
      set_min_redshift( real_type z_min );

      void
      set_max_redshift( real_type z_max );

      real_type
      min_ra() const;

      real_type
      max_ra() const;

      real_type
      min_dec() const;

      real_type
      max_dec() const;

      real_type
      min_redshift() const;

      real_type
      max_redshift() const;

      real_type
      min_dist() const;

      real_type
      max_dist() const;

      real_type
      min_dist( unsigned snap ) const;

      real_type
      max_dist( unsigned snap ) const;

      /// Set randomised lightcone. Lightcones will pass through a
      /// number of tiles, typically. These tiles can be optionally
      /// randomly rotated and shifted.
      /// @param[in] rand  Flag indicating randomness on or off.
      /// @param[in] engine  The random engine to use.
      ///
      void
      set_random( bool rand, int rng_seed=-1,
		  hpc::engine_type* engine = &hpc::engine );

      /// Get lightcone random state.
      /// @returns Lightcone random state.
      ///
      bool
      random() const;

      /// Set viewing inclination angle. This angle decides an right-
      /// ascension offset used in calculating unique cones.
      /// @param[in] angle  Viewing offset angle.
      ///
      void
      set_ra_offset( real_type angle );

      void
      set_dec_offset( real_type angle );

      /// Get viewing angle.
      /// @returns Viewing inclination angle.
      ///
      real_type
      ra_offset() const;

      real_type
      dec_offset() const;

      /// Set lightcone origin.
      /// @param orig  Lightcone origin.
      ///
      void
      set_origin( std::array<real_type,3> const& orig );

      /// Get lightcone origin.
      /// @returns Lightcone origin.
      ///
      std::array<real_type,3> const&
      origin() const;

      void
      set_single_snapshot( bool state );

      bool
      single_snapshot() const;

      void
      set_snapshot( unsigned snap );

      unsigned
      snapshot() const;

      tile_iterator
      tile_begin() const;

      tile_iterator
      tile_end() const;

      template< class Backend >
      typename Backend::lightcone_galaxy_iterator
      galaxy_begin( query<real_type>& qry,
                    Backend& be,
                    tao::batch<real_type>* bat = 0,
                    filter const* filt = 0)
      {
          typename Backend::lightcone_galaxy_iterator result=be.galaxy_begin( qry, *this, bat, filt );
    	  return result;
      }

      template< class Backend >
      typename Backend::lightcone_galaxy_iterator
      galaxy_end( query<real_type>& qry,
                  Backend& be )
      {
         return be.galaxy_end( qry, *this );
      }

      std::vector<unsigned> const&
      snapshot_bins() const;

      std::vector<real_type> const&
      distance_bins() const;

      std::vector<real_type> const&
      redshift_bins() const;

      real_type
      distance_to_redshift( real_type dist ) const;

      hpc::engine_type*
      rng_engine() const;

      void
      rng_reset() const {
          if (_rng_seed!=-1)
          {
              rng_engine()->seed(_rng_seed);
          }
      }

      bool
      overlap( std::array<real_type,3> const& min,
               std::array<real_type,3> const& max ) const;

      lightcone
      dereference() const
      {
          return *this;
      }

   protected:

      void
	  _recalc();

   protected:
      tao::simulation const* _sim;
      std::array<real_type,2> _ra;
      std::array<real_type,2> _dec;
      std::array<real_type,2> _z;
      std::array<real_type,2> _dist;
      std::vector<real_type> _dist_bins;
      std::vector<unsigned> _snap_bins;
      std::vector<real_type> _z_bins;
      hpc::num::interp<real_type> _dist_to_z;
      bool _rand;
      real_type _ra_offs, _dec_offs;
      std::array<real_type,3> _orig;
      bool _sng_snap;
      unsigned _snap;
      hpc::engine_type* _eng;
      int _rng_seed;    // for resetting the sequence in kdtree case
   };

}

#endif
