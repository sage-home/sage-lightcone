#ifndef tao_base_lightcone_galaxy_iterator_hh
#define tao_base_lightcone_galaxy_iterator_hh

#include "batch.hh"
#include "lightcone.hh"
#include "lightcone_tile_iterator.hh"
#include <libhpc/system/view.hh>

namespace tao {
using namespace hpc;

template <class BackendT>
class lightcone_galaxy_iterator
    : public boost::iterator_facade<lightcone_galaxy_iterator<BackendT>,
                                    batch<typename BackendT::real_type> &,
                                    std::forward_iterator_tag,
                                    batch<typename BackendT::real_type> &> {
  friend class boost::iterator_core_access;

public:
  typedef BackendT backend_type;
  typedef typename backend_type::real_type real_type;
  typedef typename tao::lightcone::tile_iterator tile_iterator;
  typedef typename backend_type::tile_galaxy_iterator tile_galaxy_iterator;
  typedef batch<real_type> &value_type;
  typedef value_type reference_type;

public:
  lightcone_galaxy_iterator() : _done(true) {}

  lightcone_galaxy_iterator(tao::lightcone const &lc, backend_type &be,
                            query<real_type> &qry, tile_iterator tile_it,
                            tao::batch<real_type> *bat = 0,
                            filter const *filt = 0, unsigned first_tile = 0)
      : _lc(&lc), _be(&be), _qry(&qry), _tile_it(tile_it), _bat(bat),
        _tbl_idx(0), _n_processed_gals(0), _filt(filt), _done(false),
        _first_tile(first_tile), _last_tile(first_tile) {

    // Skip over first tiles to get us back to where we
    // checkpointed.
    if (_first_tile) {
      LOGILN("Skipping first ", _first_tile, " tiles.");
    }
    for (unsigned ii = 0; ii < _first_tile && !_tile_it.done();
         ++_tile_it, ++ii)
      ;

    if (!_tile_it.done()) {
      // Pepare an iterator.
      // LOGILN("Processing tile at: [",
      //        _tile_it->min()[0] - _tile_it->origin()[0], ", ",
      //        _tile_it->min()[1] - _tile_it->origin()[1], ", ",
      //        _tile_it->min()[2] - _tile_it->origin()[2], "]");
      /*
      std::cout << "QUERY=["<<_tile_it->translation()[0]<<"
      "<<_tile_it->translation()[1]<<" "<<_tile_it->translation()[2]<<"] ["<<
                _tile_it->min()[0]<<" "<<_tile_it->min()[1]<<"
      "<<_tile_it->min()[2]<<"] ["   << _tile_it->origin()[0]<<"
      "<<_tile_it->origin()[1]<<" "<<_tile_it->origin()[2]<<"] ROT["   <<
                _tile_it->rotation()[0]<<" "<<_tile_it->rotation()[1]<<"
      "<<_tile_it->rotation()[2]<<"]"   << std::endl;
      */
      _gal_it = _be->galaxy_begin(*_qry, *_tile_it, _bat, _filt, true);
      _settle();
    } else
      _done = true;
  }

  lightcone_galaxy_iterator &operator=(const lightcone_galaxy_iterator &op) {
    _lc = op._lc;
    _be = op._be;
    _qry = op._qry;
    _bat = op._bat;
    _filt = op._filt;
    _done = op._done;
    _tile_it = op._tile_it;
    _gal_it = op._gal_it;
    _tbl_idx = op._tbl_idx;
    _n_processed_gals = op._n_processed_gals;
    _done = op._done;
    return *this;
  }

  reference_type operator*() { return _gal_it.deref(_lc); }

  bool done() const { return _done; }

  tao::lightcone const *lightcone() const { return _lc; }

  unsigned tile_index() const { return _tile_it.index(); }

  // unsigned
  // table_index() const
  // {
  //    return _tbl_idx + _gal_it.table_index();
  // }

  unsigned long long n_processed_gals() const {
    return _n_processed_gals + _gal_it.n_processed_gals();
  }

  bool should_checkpoint() {
    bool res = false;
    if (_last_tile != tile_index()) {
      _last_tile = tile_index();
      res = true;
    }
    if (!res) {
      res = _gal_it.should_checkpoint();
    }
    return res;
  }

  void save_checkpoint(boost::property_tree::ptree &pt) {
    pt.put("tile", std::to_string(_tile_it.index()));
    _gal_it.save_checkpoint(pt);
  }

protected:
  void increment() {
    ++_gal_it;
    _settle();
  }

  bool equal(const lightcone_galaxy_iterator &op) const {
    return _done == op._done;
  }

  reference_type dereference() { return *_gal_it; }

  void _settle() {
    if (_gal_it.done()) {
      do {
        // LOGILN( "Done.", setindent( -2 ) );
        ++_tile_it;
        if (_tile_it.done()) {
          _done = true;
          break;
        }

        // Update our progress through both tables and galaxies.
        // _tbl_idx += _gal_it.table_index();
        _n_processed_gals += _gal_it.n_processed_gals();

        // LOGILN("Processing tile at: [",
        //        _tile_it->min()[0] - _tile_it->origin()[0], ", ",
        //        _tile_it->min()[1] - _tile_it->origin()[1], ", ",
        //        _tile_it->min()[2] - _tile_it->origin()[2], "]");
        // std::cout << "QUERY=["<<_tile_it->translation()[0]<<"
        // "<<_tile_it->translation()[1]<<" "<<_tile_it->translation()[2]<<"]
        // ["<<
        //           _tile_it->min()[0]<<" "<<_tile_it->min()[1]<<"
        //           "<<_tile_it->min()[2]<<"] ["   << _tile_it->origin()[0]<<"
        //           "<<_tile_it->origin()[1]<<" "<<_tile_it->origin()[2]<<"]
        //           ROT["   << _tile_it->rotation()[0]<<"
        //           "<<_tile_it->rotation()[1]<<"
        //           "<<_tile_it->rotation()[2]<<"]"   << std::endl;
        _gal_it = _be->galaxy_begin(*_qry, *_tile_it, _bat, _filt);
      } while (_gal_it.done());
    }
  }

protected:
  tao::lightcone const *_lc;
  backend_type *_be;
  query<real_type> *_qry;
  tao::batch<real_type> *_bat;
  filter const *_filt;
  unsigned _first_tile, _last_tile;
  tile_iterator _tile_it;
  tile_galaxy_iterator _gal_it;
  unsigned _tbl_idx;
  unsigned long long _n_processed_gals;
  bool _done;
};

} // namespace tao

#endif
