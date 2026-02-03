#ifndef tao_base_box_hh
#define tao_base_box_hh

#include "batch.hh"
#include "filter.hh"
#include "query.hh"
#include <array>
#include <libhpc/logging.hh>
#include <libhpc/mpi/init.hh>
#include <libhpc/numerics/clip.hh>
#include <libhpc/system/random.hh>

namespace tao {
using namespace hpc;

class simulation;

template <class T> class box {
public:
  typedef T real_type;

public:
  box(const tao::simulation *sim = 0)
      : _sim(sim), _snap(0), _rand(false), _eng(0), _orig{
                                                        {0.0, 0.0, 0.0},
                                                    } {
    _update();
  }

  void set_simulation(const tao::simulation *sim) {
    _sim = sim;
    _update();
  }

  void set_size(real_type size) {
    std::fill(_min.begin(), _min.end(), 0.0);
    std::fill(_max.begin(), _max.end(), size);
    for (int ii = 0; ii < 3; ++ii) {
      _bnds[ii][0] = _min[ii];
      _bnds[ii][1] = _max[ii];
    }
  }

  void set_snapshot(unsigned snap) { _snap = snap; }

  unsigned snapshot() const { return _snap; }

  real_type redshift() const { return _sim->redshift(_snap); }

  void set_random(bool rand, engine_type *engine = &hpc::engine) {
    _rand = rand;
    _eng = engine;
    _update_random();
  }

  void set_origin(std::array<real_type, 3> const &orig) { _orig = orig; }

  std::array<real_type, 3> const &origin() const { return _orig; }

  const tao::simulation *simulation() const { return _sim; }

  const std::array<real_type, 3> &min() const { return _min; }

  const std::array<real_type, 3> &max() const { return _max; }

  std::array<std::array<real_type, 2>, 3> const &bounds() const {
    return _bnds;
  }

  bool random() const { return _rand; }

  const std::array<unsigned, 3> &rotation() const { return _rot; }

  const std::array<real_type, 3> &translation() const { return _trans; }

  template <class Backend>
  typename Backend::box_galaxy_iterator
  galaxy_begin(tao::query<real_type> &query, Backend &be,
               tao::batch<real_type> *bat = 0, tao::filter const *filt = 0) {
    return be.galaxy_begin(query, *this, bat, filt);
  }

  template <class Backend>
  typename Backend::box_galaxy_iterator galaxy_end(tao::query<real_type> &query,
                                                   Backend &be) const {
    return be.galaxy_end(query, *this);
  }

  bool contains(std::array<real_type, 3> const &crd, unsigned snap = 0) const {
    return crd[0] >= _min[0] && crd[0] < _max[0] && crd[1] >= _min[1] &&
           crd[1] < _max[1] && crd[2] >= _min[2] && crd[2] < _max[2];
  }

  bool contains(double const &x, double const &y, double const &z) const {
    return x >= _min[0] && x < _max[0] && y >= _min[1] && y < _max[1] &&
           z >= _min[2] && z < _max[2];
  }

  template <class Iter>
  bool box_collision(Iter const &bnds_begin, unsigned snap = 0) const {
    return box_box_collision(_bnds.begin(), _bnds.end(), bnds_begin);
  }

  void set_calls(uint64_t calls) {
    for (uint64_t i = 0; i < calls; i++) {
      randomise();
    }
  }

  void randomise() {
    // Only do this if we have a simulation set.
    if (_sim) {
      hpc::mpi::set_rng(hpc::mpi::get_rng() + 1); // need a better place

      LOGBLOCKD("Randomising box.");

      // Generate translation.
      for (unsigned ii = 0; ii < 3; ++ii)
        _trans[ii] = generate_uniform<real_type>(0.0, _sim->box_size(), *_eng);
      // RS: TODO TAKE THESE NEXT TWO LINE OUT!!!!!! JUST FOR TESTING
      // for( unsigned ii = 0; ii < 3; ++ii )
      //    _trans[ii] = 0;
      LOGDLN("Translation: ", _trans);

      // Generate rotations.
      int rnd = generate_uniform<int>(0, 5, *_eng);
      // RS: TODO TAKE THIS NEXT LINE OUT!!!!!! JUST FOR TESTING
      //        rnd = 1;
      switch (rnd) {
      case 0:
        _rot[0] = 0;
        _rot[1] = 1;
        _rot[2] = 2;
        break;
      case 1:
        _rot[0] = 2;
        _rot[1] = 0;
        _rot[2] = 1;
        break;
      case 2:
        _rot[0] = 1;
        _rot[1] = 2;
        _rot[2] = 0;
        break;
      case 3:
        _rot[0] = 0;
        _rot[1] = 2;
        _rot[2] = 1;
        break;
      case 4:
        _rot[0] = 1;
        _rot[1] = 0;
        _rot[2] = 2;
        break;
      case 5:
        _rot[0] = 2;
        _rot[1] = 1;
        _rot[2] = 0;
        break;
      };
      LOGDLN("Rotation: ", _rot);
    }
  }

protected:
  void _update() {
    _update_random();

    // If simulation is set, update min and max.
    if (_sim) {
      std::fill(_min.begin(), _min.end(), 0);
      std::fill(_max.begin(), _max.end(), _sim->box_size());
      for (int ii = 0; ii < 3; ++ii) {
        _bnds[ii][0] = _min[ii];
        _bnds[ii][1] = _max[ii];
      }
    }
  }

  void _update_random() {
    // If randomised, then randomise!
    if (_rand) {
      randomise();
    } else {
      std::iota(_rot.begin(), _rot.end(), 0);
      std::fill(_trans.begin(), _trans.end(), 0);
    }
  }

protected:
  const tao::simulation *_sim;
  std::array<real_type, 3> _min, _max;
  std::array<std::array<real_type, 2>, 3> _bnds;
  bool _rand;
  engine_type *_eng;
  std::array<unsigned, 3> _rot;
  std::array<real_type, 3> _trans;
  std::array<real_type, 3> _orig;
  unsigned _snap;
};

} // namespace tao

#endif
