#include "simulation.hh"
#include <libhpc/debug/assert.hh>
#include <libhpc/logging.hh>

namespace tao {

simulation::simulation()
    : _box_size(0), _hubble(0.0), _h(0.0), _omega_m(0.0), _omega_l(0.0),
      _omega_r(0.0), _omega_k(0.0) {}

simulation::simulation(real_type box_size, real_type hubble, real_type omega_m,
                       real_type omega_l,
                       std::map<int, real_type> const &snap_zs)
    : _box_size(box_size), _zs(snap_zs) {
  set_cosmology(hubble, omega_m, omega_l);
}

simulation::simulation(real_type box_size, real_type hubble, real_type omega_m,
                       real_type omega_l, unsigned num_snaps, ...)
    : _box_size(box_size) {
  set_cosmology(hubble, omega_m, omega_l);

  // Extract the expansion list and convert to redshifts.
  va_list vl;
  va_start(vl, num_snaps);
  for (unsigned ii = 0; ii < num_snaps; ++ii) {
    real_type v = va_arg(vl, real_type);
    _zs[ii] = expansion_to_redshift(v);
  }

  // We really want these to be ordered by snapshots, and that
  // usually means oldest (largest redshifts) first.
#ifndef NDEBUG
  for (unsigned ii = 1; ii < _zs.size(); ++ii) {
    ASSERT(_zs[ii] < _zs[ii - 1],
           "Expansion factor list must be supplied in snapshot order, "
           "which must be oldest first.");
  }
#endif
}

void simulation::set_box_size(real_type box_size) { _box_size = box_size; }

real_type simulation::box_size() const { return _box_size; }

void simulation::set_cosmology(real_type hubble, real_type omega_m,
                               real_type omega_l) {
  // LOGILN( "Setting simulation cosmology:", setindent( 2 ) );
  _hubble = hubble;
  // LOGILN( "Hubble: ", hubble );
  _h = _hubble / 100.0;
  _omega_m = omega_m;
  // LOGILN( "Omega M: ", omega_m );
  _omega_l = omega_l;
  // LOGILN( "Omega L: ", omega_l );
  _omega_r = 4.165e-5 / (_h * _h);
  // LOGILN( "Omega R: ", _omega_r );
  _omega_k = 1 - _omega_m - _omega_l - _omega_r;
  // LOGILN( "Omega K: ", _omega_k );
  // LOGILN( "Done.", setindent( -2 ) );
}

real_type simulation::hubble() const { return _hubble; }

real_type simulation::h() const { return _h; }

real_type simulation::omega_m() const { return _omega_m; }

real_type simulation::omega_l() const { return _omega_l; }

real_type simulation::omega_r() const { return _omega_r; }

real_type simulation::omega_k() const { return _omega_k; }

unsigned simulation::num_snapshots() const { return _zs.size(); }

real_type simulation::redshift(unsigned snap) const { return _zs.at(snap); }

std::map<int, real_type> const &simulation::redshifts() const { return _zs; }

} // namespace tao
