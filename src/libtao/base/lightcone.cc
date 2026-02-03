#include "lightcone.hh"
#include "lightcone_tile_iterator.hh"
#include <libhpc/numerics/clip.hh>
#include <libhpc/numerics/coords.hh>
#include <libhpc/system/deallocate.hh>
#include <libhpc/system/reallocate.hh>

namespace tao {

lightcone::lightcone(const tao::simulation *sim)
    : _sim(NULL), _rand(false), _ra_offs(0.0), _dec_offs(0.0),
      _orig{{0.0, 0.0, 0.0}}, _sng_snap(false), _snap(0), _eng(&hpc::engine) {
  set_geometry(0, 10, 0, 10, 0.06);
  set_simulation(sim);
}

void lightcone::set_simulation(const tao::simulation *sim) {
  _sim = sim;
  _recalc();
}

void lightcone::copy(const lightcone &from) {
  _sim = from._sim;
  _ra = from._ra;
  _dec = from._dec;
  _z = from._z;
  _dist = from._dist;
  _dist_bins = from._dist_bins;
  _snap_bins = from._snap_bins;
  _z_bins = from._z_bins;
  _dist_to_z = from._dist_to_z;
  _rand = from._rand;
  _ra_offs = from._ra_offs;
  _dec_offs = from._dec_offs;
  _orig = from._orig;
  _sng_snap = from._sng_snap;
  _snap = from._snap;
  _eng = from._eng;
}

void lightcone::set_geometry(real_type ra_min, real_type ra_max,
                             real_type dec_min, real_type dec_max,
                             real_type z_max, real_type z_min) {
  // Check values make sense.
  ASSERT(ra_min >= 0.0 && ra_min <= 360.0,
         "Minimum RA cannot be less than zero or greater than 360 degrees.");
  ASSERT(ra_max >= 0.0 && ra_max <= 360.0,
         "Maximum RA cannot be less than zero or greater than 360 degrees.");
  ASSERT(ra_min < ra_max, "Minimum RA must be less than maximum RA.");
  ASSERT(dec_min >= -90 && dec_min <= 90.0,
         "Minimum DEC cannot be less than -90 or greater than 90 degrees.");
  ASSERT(dec_max >= -90.0 && dec_max <= 90.0,
         "Maximum DEC cannot be less than -90 or greater than 90 degrees.");
  ASSERT(dec_min < dec_max, "Minimum DEC must be less than maximum DEC.");
  ASSERT(z_min >= 0.0,
         "Minimum redshift must be greater than or equal to zero.");
  ASSERT(z_min < z_max, "Minimum redshift must be less than maximum redshift.");

  // Set values.
  _ra[0] = to_radians(ra_min);
  _ra[1] = to_radians(ra_max);
  _dec[0] = to_radians(dec_min);
  _dec[1] = to_radians(dec_max);
  _z[0] = z_min;
  _z[1] = z_max;
  _recalc();
}

void lightcone::set_random(bool rand, int rng_seed, engine_type *engine) {
  _rand = rand;
  _rng_seed = rng_seed;
  _eng = engine;
}

void lightcone::set_ra_offset(real_type angle) { _ra_offs = angle; }

void lightcone::set_dec_offset(real_type angle) { _dec_offs = angle; }

real_type lightcone::ra_offset() const { return _ra_offs; }

real_type lightcone::dec_offset() const { return _dec_offs; }

void lightcone::set_origin(std::array<real_type, 3> const &orig) {
  _orig = orig;
}

std::array<real_type, 3> const &lightcone::origin() const { return _orig; }

void lightcone::set_single_snapshot(bool state) { _sng_snap = state; }

bool lightcone::single_snapshot() const { return _sng_snap; }

void lightcone::set_snapshot(unsigned snap) { _snap = snap; }

unsigned lightcone::snapshot() const { return _snap; }

const tao::simulation *lightcone::simulation() const { return _sim; }

void lightcone::set_min_ra(real_type ra_min) {
  ASSERT(ra_min >= 0.0 && ra_min <= 360.0,
         "Minimum RA cannot be less than zero or greater than 360 degrees.");
  ASSERT(to_radians(ra_min) < _ra[1],
         "Minimum RA must be less than maximum RA.");
  _ra[0] = to_radians(ra_min);
  _recalc();
}

void lightcone::set_max_ra(real_type ra_max) {
  ASSERT(ra_max >= 0.0 && ra_max <= 360.0,
         "Maximum RA cannot be less than zero or greater than 360 degrees.");
  ASSERT(to_radians(ra_max) > _ra[0],
         "Minimum RA must be less than maximum RA.");
  _ra[1] = to_radians(ra_max);
  _recalc();
}

void lightcone::set_min_dec(real_type dec_min) {
  ASSERT(dec_min >= -90.0 && dec_min <= 90.0,
         "Minimum DEC cannot be less than -90 or greater than 90 degrees.");
  ASSERT(to_radians(dec_min) < _dec[1],
         "Minimum DEC must be less than maximum DEC.");
  _dec[0] = to_radians(dec_min);
  _recalc();
}

void lightcone::set_max_dec(real_type dec_max) {
  ASSERT(dec_max >= -90.0 && dec_max <= 90.0,
         "Maximum DEC cannot be less than -90 or greater than 90 degrees.");
  ASSERT(to_radians(dec_max) > _dec[0],
         "Minimum DEC must be less than maximum DEC.");
  _dec[1] = to_radians(dec_max);
  _recalc();
}

void lightcone::set_min_redshift(real_type z_min) {
  ASSERT(z_min >= 0.0,
         "Minimum redshift must be greater than or equal to zero.");
  ASSERT(z_min < _z[1], "Minimum redshift must be less than maximum redshift.");
  _z[0] = z_min;
  _recalc();
}

void lightcone::set_max_redshift(real_type z_max) {
  ASSERT(_z[0] < z_max, "Minimum redshift must be less than maximum redshift.");
  _z[1] = z_max;
  _recalc();
}

lightcone::tile_iterator lightcone::tile_begin() const {
  return tile_iterator(*this);
}

lightcone::tile_iterator lightcone::tile_end() const { return tile_iterator(); }

real_type lightcone::min_ra() const { return _ra[0]; }

real_type lightcone::max_ra() const { return _ra[1]; }

real_type lightcone::min_dec() const { return _dec[0]; }

real_type lightcone::max_dec() const { return _dec[1]; }

real_type lightcone::min_redshift() const { return _z[0]; }

real_type lightcone::max_redshift() const { return _z[1]; }

real_type lightcone::min_dist() const { return _dist[0]; }

real_type lightcone::max_dist() const { return _dist[1]; }

real_type lightcone::min_dist(unsigned snap) const {
  ASSERT(snap >= _snap_bins.front() && snap <= _snap_bins.back(),
         "Invalid snapshot.");
  return _dist_bins[snap - _snap_bins.front() + 1];
}

real_type lightcone::max_dist(unsigned snap) const {
  ASSERT(snap >= _snap_bins.front() && snap <= _snap_bins.back(),
         "Invalid snapshot.");
  return _dist_bins[snap - _snap_bins.front()];
}

std::vector<unsigned> const &lightcone::snapshot_bins() const {
  return _snap_bins;
}

std::vector<real_type> const &lightcone::distance_bins() const {
  return _dist_bins;
}

std::vector<real_type> const &lightcone::redshift_bins() const {
  return _z_bins;
}

real_type lightcone::distance_to_redshift(real_type dist) const {
  return _dist_to_z[dist];
}

bool lightcone::random() const { return _rand; }

engine_type *lightcone::rng_engine() const { return _eng; }

bool lightcone::overlap(std::array<real_type, 3> const &min,
                        std::array<real_type, 3> const &max) const {
  LOGBLOCKT("Checking overlap in box and lightcone: ", min, " - ", max);

  // Precalculate some things for distance evaluations.
  std::array<real_type, 3> l, u, ecs_min, ecs_max;
  ecs_min[0] = _ra[0] + _ra_offs;
  ecs_max[0] = _ra[1] + _ra_offs;
  ecs_min[1] = _dec[0] + _dec_offs;
  ecs_max[1] = _dec[1] + _dec_offs;
  ecs_min[2] = _dist[0];
  ecs_max[2] = _dist[1];
  for (unsigned ii = 0; ii < 3; ++ii) {
    l[ii] = min[ii] - _orig[ii];
    u[ii] = max[ii] - _orig[ii];
  }

  // Call out to the collision routine.
  return hpc::ecs_box_collision(ecs_min, ecs_max, l, u);
}

void lightcone::_recalc() {
  if (_sim) {
    LOGBLOCKD("Recalculating lightcone information.");

    _dist[0] =
        hpc::num::redshift_to_comoving_distance(
            _z[0], 1000, _sim->hubble(), _sim->omega_l(), _sim->omega_m()) *
        _sim->h();
    _dist[1] =
        hpc::num::redshift_to_comoving_distance(
            _z[1], 1000, _sim->hubble(), _sim->omega_l(), _sim->omega_m()) *
        _sim->h();
    // LOGDLN( "Distance range: [", _dist[0], ", ", _dist[1], ")" );
    // std::cout <<"Distance range: ["<<_dist[0]<<", "<<_dist[1]<<
    // ")"<<std::endl;

    // Prepare the redshift distance bins. Note that I will incorporate the
    // minimum and maximum redshift here. First calculate the
    // number of bins in the redshift range.
    hpc::deallocate(_dist_bins);
    hpc::deallocate(_snap_bins);
    hpc::deallocate(_z_bins);
    unsigned first, last;
    {
      unsigned ii = 0;
      unsigned ns = _sim->num_snapshots();
      while (ii < ns && _sim->redshift(ii) > _z[1])
        ++ii;
      first = ii;
      while (ii < ns && _sim->redshift(ii) > _z[0])
        ++ii;
      last = ii + 1;
    }

    // Store the distances.
    hpc::reallocate(_dist_bins, last - first + 1);
    hpc::reallocate(_snap_bins, last - first);
    hpc::reallocate(_z_bins, last - first + 1);
    _z_bins.back() = std::min(_sim->redshift(0), _z[0]);
    _dist_bins.back() = hpc::num::redshift_to_comoving_distance(
        _z_bins.back(), 1000, _sim->hubble(), _sim->omega_l(), _sim->omega_m());
    _z_bins.front() =
        std::max(_sim->redshift(_sim->num_snapshots() - 1), _z[1]);
    _dist_bins.front() = hpc::num::redshift_to_comoving_distance(
        _z_bins.front(), 1000, _sim->hubble(), _sim->omega_l(),
        _sim->omega_m());
    for (unsigned ii = 1; ii < _dist_bins.size() - 1; ++ii) {
      _z_bins[ii] = _sim->redshift(first + ii - 1);
      _dist_bins[ii] = hpc::num::redshift_to_comoving_distance(
          _z_bins[ii], 1000, _sim->hubble(), _sim->omega_l(), _sim->omega_m());
    }
    for (unsigned ii = 0; ii < _snap_bins.size(); ++ii)
      _snap_bins[ii] = first + ii;

    // The distances I've calculated are all in Mpc. I really want them in
    // Mpc/h.
    std::transform(_dist_bins.begin(), _dist_bins.end(), _dist_bins.begin(),
                   [this](real_type d) { return d * this->_sim->h(); });

    // Log to debugging stream.
    LOGDLN("Distance bins: ", _dist_bins);
    LOGDLN("Snapshot bins: ", _snap_bins);

    // Build the distance to redshift interpolator.
    {
      // How many points do I need?
      unsigned num_points =
          std::max<unsigned>((unsigned)((_z[1] - _z[0]) / 0.001), 2);

      // Setup arrays.
      std::vector<real_type> dists(num_points);
      std::vector<real_type> zs(num_points);
      for (unsigned ii = 0; ii < num_points; ++ii) {
        zs[ii] = _z[0] + (_z[1] - _z[0]) *
                             ((real_type)ii / (real_type)(num_points - 1));
        dists[ii] = hpc::num::redshift_to_comoving_distance(
                        zs[ii], 1000, _sim->hubble(), _sim->omega_l(),
                        _sim->omega_m()) *
                    _sim->h(); // don't forget to put in Mpc/h
      }

      // Transfer to interpolator.
      _dist_to_z.set_abscissa(dists);
      _dist_to_z.set_values(zs);
    }
  }
}
} // namespace tao
