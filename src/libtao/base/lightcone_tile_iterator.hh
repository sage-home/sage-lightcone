#ifndef tao_base_lightcone_tile_iterator_hh
#define tao_base_lightcone_tile_iterator_hh

#include "simulation.hh"
#include "tile.hh"
#include "types.hh"
#include <array>
#include <boost/property_tree/ptree.hpp>
#include <libhpc/numerics/coords.hh>

namespace tao {

class lightcone_tile_iterator
    : public boost::iterator_facade<lightcone_tile_iterator, tile<real_type>,
                                    std::forward_iterator_tag,
                                    tile<real_type>> {
  friend class boost::iterator_core_access;

public:
  typedef tile<real_type> value_type;
  typedef value_type reference_type;

public:
  lightcone_tile_iterator();

  lightcone_tile_iterator(tao::lightcone const &lc);

  tao::lightcone const *lightcone() const;

  bool done() const;

  unsigned index() const;

  std::list<std::array<int, 3>> const &remaining_tiles() const;

  std::map<std::array<int, 3>, int> const &done_tiles() const;

  void save_checkpoint(boost::property_tree::ptree &pt) const;

protected:
  void increment();

  bool equal(lightcone_tile_iterator const &op) const;

  reference_type dereference() const;

protected:
  tao::lightcone const *_lc;
  std::array<real_type, 3> _cur;
  std::list<std::array<int, 3>> _rem_tiles;
  std::map<std::array<int, 3>, int> _done_tiles;
  bool _done;
  unsigned _idx;
  value_type _tile;
};

} // namespace tao

#endif
