#ifndef tao_base_age_line_hh
#define tao_base_age_line_hh

#include <vector>
#include <fstream>

// Ignore soci warnings
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
//#include <soci.h>
#pragma GCC diagnostic pop

#include <libhpc/debug/except.hh>
#include <libhpc/logging.hh>
#include <libhpc/system/filesystem.hh>
#include <libhpc/system/view.hh>
#include <libhpc/system/assign.hh>
#include <libhpc/system/deallocate.hh>
#include <libhpc/system/reallocate.hh>
#include "utils.hh"
#include "simulation.hh"

namespace tao {

   /// The `age_line` utility class represents a sequence of points along a timeline.
   /// It differs from a simple array of values in that it handle calculation of the
   /// "dual", and provides support for locating which age bucket an arbitrary time-point
   /// falls.
   ///
   template< class T >
   class age_line
   {
   public:

      typedef T real_type;

   public:

      /// Default constructor.
      ///
      age_line()
      {
      }

      /// Construct an age-line object from a sequence of ages. Ages must be specified in giga-years
      /// in order to be consistent with other parts of the TAO calculations. The template `Seq` may
      /// be any forward-iterable sequence.
      ///
      /// @tparam Seq  A forward-iterator sequence.
      /// @param[in] ages  A sequence of ages.
      ///
      template< class Seq >
      age_line( Seq&& ages )
      {
         set_ages( std::forward<Seq>( ages ) );
      }

      /// Construct an age-line object from the contents an SQL database. This is a helper function
      /// and was never really intended for use outside of a special application in TAO.
      /// The database must contain a table named `snap_redshift`, which itself is composed of two
      /// columns. The first column is `redshifts`, listing the redshifts of each snapshot in the
      /// simulation. The second column is `snapnum`, which refers to the ordinal of the snapshot.
      ///
      /// @param[in] sql  SOCI session connected to a database.
      /// @param[in] hubble  Hubble constant.
      ///
      //age_line(soci::session& sql,
      //          real_type hubble = 73.0,
      //          real_type omega_m = 0.25,
      //          real_type omega_l = 0.75 )
      //{
      //   load_ages( sql, hubble, omega_m, omega_l );
      //}

      /// Construct an age-line object from the contents of a file. The file should begin with
      /// a single line containing the number of ages. Each subsequent line should contain an
      /// age value in gigayears.
      ///
      age_line( hpc::fs::path const& path )
      {
         load( path );
      }

      //age_line( soci::session& sql,
      //          const simulation& sim )
      //{
      //   load_ages( sql, sim.hubble(), sim.omega_m(), sim.omega_l() );
      //}

      void
      clear()
      {
         hpc::deallocate( _ages );
         hpc::deallocate( _dual );
      }

      unsigned
      size() const
      {
         return _ages.size();
      }

      real_type
      dual( unsigned idx ) const
      {
         return _dual[idx];
      }

      typename hpc::view<std::vector<real_type>> const
      dual() const
      {
         return _dual;
      }

      typename hpc::view<std::vector<real_type>> const
      ages() const
      {
         return _ages;
      }

      /// @sa age_line::age_line( Seq&& ages )
      ///
      template< class Seq >
      void
      set_ages( Seq&& ages )
      {
         LOGBLOCKT( "Setting ages on age line." );

         // Clear existing values.
         clear();

         // Don't do anything if we have no ages.
         if( ages.size() )
         {
            // Assign age values.
            hpc::assign( _ages, std::forward<Seq>( ages ) );

            // Must be sure ages are sorted correctly.
#ifndef NDEBUG
            for( unsigned ii = 1; ii < _ages.size(); ++ii )
               ASSERT( _ages[ii] >= _ages[ii - 1], "Ages must be sorted in ascending order." );
#endif

            // Calculate dual.
            _calc_dual();
         }
      }

      /// @sa age_line::age_line( hpc::fs::path const& path )
      ///
      void
      load( hpc::fs::path const& path )
      {
         LOGBLOCKI( "Loading ages from: ", path );

         // Open the file.
         std::ifstream file( path.c_str() );
         EXCEPT( file.is_open(), "Couldn't find ages file: ", path );

         // Read the number of ages.
         unsigned num_ages;
         file >> num_ages;
         EXCEPT( file.good(), "Error reading ages file." );
         LOGILN( "Number of ages: ", num_ages );

         // Read the ages.
         std::vector<real_type> bin_ages( num_ages );
         for( unsigned ii = 0; ii < num_ages; ++ii )
            file >> bin_ages[ii];
         EXCEPT( file.good(), "Error reading ages file." );

         // Setup the bin ages.
         set_ages( bin_ages );
      }

      /// @sa age_line::age_line( soci::session&, real_type, real_type, real_type )
      ///
      //void
      //load_ages( soci::session& sql,
      //           real_type hubble = 73.0,
      //           real_type omega_m = 0.25,
      //           real_type omega_l = 0.75 )
      //{
      //   LOGBLOCKD( "Loading ages from database." );

      //   // Clear existing values.
      //   clear();

      //   // Find number of snapshots and resize the containers.
      //   unsigned num_snaps;
      //   sql << "SELECT COUNT(*) FROM snap_redshift", soci::into( num_snaps );
      //   LOGDLN( "Number of snapshots: ", num_snaps );

      //   // Need space to store the snapshots.
      //   hpc::reallocate( _ages, num_snaps );

      //   // Read meta data.
      //   sql << "SELECT redshift FROM snap_redshift ORDER BY snapnum",
      //      soci::into( _ages );
      //   LOGTLN( "Redshifts: ", _ages );

         // Convert to ages.
      //   for( unsigned ii = 0; ii < _ages.size(); ++ii )
      //      _ages[ii] = redshift_to_age<real_type>( _ages[ii], hubble, omega_m, omega_l );
      //   LOGDLN( "Snapshot ages: ", _ages );

      //   // Calculate the dual.
      //   _calc_dual();
      //}

      /// Scan the age-bins currently set for which bin contains the time-point `age`. If the point
      /// exists before the first age-bin, the first bin is returned. Similarly, if the point
      /// exists beyond the final age-bin, the last bin is returned.
      ///
      unsigned
      find_bin( real_type age ) const
      {
         LOGBLOCKT( "Searching for bin using age: ", age );
         unsigned bin;
         {
            // Use binary search to find first element greater.
            auto it = std::lower_bound( _dual.begin(), _dual.end(), age );
            if( it == _dual.end() )
               bin = _dual.size();
            else
               bin = it - _dual.begin();
         }
         LOGTLN( "Found bin ", bin, " with age of ", _ages[bin], "." );
         return bin;
      }

      real_type
      operator[]( unsigned idx ) const
      {
         return _ages[idx];
      }

   protected:

      void
      _calc_dual()
      {
         if( _ages.size() )
         {
            _dual.resize( _ages.size() - 1 );
            for( unsigned ii = 1; ii < _ages.size(); ++ii )
               _dual[ii - 1] = 0.5*(_ages[ii] + _ages[ii - 1]);
         }
         else
            hpc::deallocate( _dual );
         LOGTLN( "Dual: ", _dual );
      }

   protected:

      std::vector<real_type> _ages, _dual;
   };
}

#endif
