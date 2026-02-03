#ifndef tao_base_simulation_hh
#define tao_base_simulation_hh

#include "types.hh"
#include "utils.hh"
#include <libhpc/system/assign.hh>
#include <stdarg.h>
#include <vector>

namespace tao {

///
/// Represents a dataset simulation. Simulations define various
/// characterstics of each dataset, those begin:
///
///   * simulation box size,
///   * hubble constant,
///   * OmegaM,
///   * OmegaL,
///   * OmegaK,
///   * OmegaR,
///   * simulation snapshot redshifts.
///
class simulation {
public:
  ///
  /// Default constructor.
  ///
  simulation();

  ///
  /// Construct with cosmology and redshifts.
  ///
  /// @param[in] box_size Simulation box size.
  /// @param[in] hubble Hubble constant.
  /// @param[in] omega_m Cosmology constant.
  /// @param[in] omega_l Cosmology constant.
  /// @param[in] snap_zs Snapshot redshifts.
  ///
  simulation(real_type box_size, real_type hubble, real_type omega_m,
             real_type omega_l, std::map<int, real_type> const &snap_zs);

  ///
  /// Construct with cosmology and redshifts.
  ///
  /// @param[in] box_size Simulation box size.
  /// @param[in] hubble Hubble constant.
  /// @param[in] omega_m Cosmology constant.
  /// @param[in] omega_l Cosmology constant.
  /// @param[in] num_snaps Number of snapshots.
  /// @param[in] ... Expansion factor of each snapshot.
  ///
  simulation(real_type box_size, real_type hubble, real_type omega_m,
             real_type omega_l, unsigned num_snaps, ...);

  ///
  /// Set simulation box size.
  ///
  /// @param[in] box_size Simulation box size.
  ///
  void set_box_size(real_type box_size);

  real_type box_size() const;

  ///
  /// Set simulation cosmology.
  ///
  /// @param[in] hubble Hubble constant.
  /// @param[in] omega_m Cosmology constant.
  /// @param[in] omega_l Cosmology constant.
  ///
  void set_cosmology(real_type hubble, real_type omega_m, real_type omega_l);

  real_type hubble() const;

  real_type h() const;

  real_type omega_m() const;

  real_type omega_l() const;

  real_type omega_r() const;

  real_type omega_k() const;

  template <class Seq> void set_snapshot_redshifts(Seq &&redshifts) {
    _zs.clear();
    _zs.insert(redshifts.begin(), redshifts.end());
  }

  unsigned num_snapshots() const;

  real_type redshift(unsigned snap) const;

  std::map<int, real_type> const &redshifts() const;

protected:
  real_type _box_size;
  real_type _hubble;
  real_type _h;
  real_type _omega_m;
  real_type _omega_l;
  real_type _omega_r;
  real_type _omega_k;
  std::map<int, real_type> _zs;
};

} // namespace tao

#endif
