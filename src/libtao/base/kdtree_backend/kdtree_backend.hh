#ifndef tao_base_kdtree_backend_kdtree_backend_hh
#define tao_base_kdtree_backend_kdtree_backend_hh

#include "../backend.hh"
#include "../batch.hh"
#include "../box.hh"
#include "../filter.hh"
#include "../query.hh"
#include "../tile.hh"
#include "../xml_dict.hh"
#include "lightcone_galaxy_iterator.hh"
#include <libhpc/algorithm/kdtree.hh>
#include <libhpc/h5/file.hh>
#include <libhpc/system/filesystem.hh>

// extern int iiii;
namespace tao {
namespace backends {

template <class> class kdtree_galaxy_iterator;

// class kdtree_box_galaxy_iterator;
// class kdtree_tile_galaxy_iterator;

///
/// Backend for use with kd-tree structured data.
///
class kdtree_backend : public backend {
public:
  typedef tao::real_type real_type;
  typedef kdtree_galaxy_iterator<box<real_type>> box_galaxy_iterator;
  typedef kdtree_galaxy_iterator<tile<real_type>> tile_galaxy_iterator;
  typedef kdtree_lightcone_galaxy_iterator<kdtree_backend>
      lightcone_galaxy_iterator;

  struct lightcone_data {
    std::array<real_type, 3> crd;
    unsigned long long gidx;
    unsigned subsize;
  };

  struct sed_data {
    int descendant;
    int snapshot;
    int local_index;
    int merge_type;
    real_type dt;
    real_type disk_sfr;
    real_type bulge_sfr;
    real_type disk_sfr_z;
    real_type bulge_sfr_z;
  };

  struct mydata {
    real_type stellarmass;
    real_type bulgemass;
    real_type blackholemass;
    real_type coldgas;
    real_type hotgas;
    real_type ejectedmass;
    real_type ics;
    real_type metalsstellarmass;
    real_type metalsbulgemass;
    real_type metalscoldgas;
    real_type metalshotgas;
    real_type metalsejectedmass;
    real_type metalsics;
    int objecttype;
    real_type diskscaleradius;
    real_type totsfr;
    real_type cooling;
    real_type heating;
    real_type quasarmodebhaccretionmass;
    real_type timeoflastmajormerger;
    real_type timeoflastminormerger;
    real_type outflowrate;
    real_type mvir;
    real_type rvir;
    real_type vvir;
    real_type vmax;
    real_type veldisp;
    real_type spin_x;
    real_type spin_y;
    real_type spin_z;
    int len;
    real_type centralmvir;
    real_type infallmvir;
    real_type infallvvir;
    real_type infallvmax;
    real_type posx;
    real_type posy;
    real_type posz;
    real_type velx;
    real_type vely;
    real_type velz;
    int snapnum;
    unsigned long long galaxyindex;
    unsigned long long centralgalaxyindex;
    unsigned long long simulationhaloindex;
    int mergeintoid;
    int mergeintosnapnum;
    int mergetype;
    real_type dt;
    real_type vpeak;
    int isflyby;
    unsigned long long globalindex;
    int descendant;
    real_type sfrdisk;
    real_type sfrbulge;
    real_type sfrdiskz;
    real_type sfrbulgez;
    int treeindex;
    int localindex;
    unsigned long long globaldescendant;
    int subsize;
  };

public:
  hpc::h5::datatype mem_type;

  kdtree_backend(hpc::mpi::comm const &comm = hpc::mpi::comm::world);

  kdtree_backend(hpc::fs::path const &fn,
                 hpc::mpi::comm const &comm = hpc::mpi::comm::world);

  hpc::mpi::comm const &comm() const;

  void connect(const cli_dict &global_cli_dict);

  void close();

  void open(hpc::fs::path const &fn);

  virtual tao::simulation const *load_simulation();

  virtual void add_conditional_fields(query<real_type> &qry);

  void load_snapshot(unsigned snap);

  // TODO: This should be removed in favour of a metadata
  //  source. Perhaps something like "get_metadata", which
  //  takes a name and returns the details?
  real_type get_max_smoothing();

  std::string make_snapshot_name(unsigned snap) const;

  box_galaxy_iterator galaxy_begin(query<real_type> &qry,
                                   box<real_type> const &box,
                                   batch<real_type> *bat = nullptr,
                                   filter const *flt = nullptr);

  box_galaxy_iterator galaxy_end(query<real_type> &qry,
                                 box<real_type> const &box) const;

  lightcone_galaxy_iterator galaxy_begin(query<real_type> &qry,
                                         lightcone const &lc,
                                         batch<real_type> *bat = nullptr,
                                         filter const *flt = nullptr);

  lightcone_galaxy_iterator galaxy_end(query<real_type> const &qry,
                                       lightcone const &lc) const;

  tile_galaxy_iterator galaxy_begin(query<real_type> &qry,
                                    tile<real_type> const &tile,
                                    batch<real_type> *bat = nullptr,
                                    filter const *flt = nullptr,
                                    bool first = false);

  tile_galaxy_iterator galaxy_end(query<real_type> const &qry,
                                  tile<real_type> const &tile) const;

  void load_lightcone_data(unsigned cell,
                           std::vector<lightcone_data> &data) const;

  hpc::h5::file const &kdtree_file() const;

  hpc::kdtree<real_type> const &kdtree() const;

  std::vector<unsigned> const &cell_counts() const;

  std::vector<unsigned long long> const &cell_offs() const;

  // Debug methods for inspecting kdtree structure
  bool is_kdtree_loaded() const;
  unsigned get_kdtree_dims() const;
  unsigned get_kdtree_branches() const;
  unsigned get_kdtree_leafs() const;
  std::string debug_kdtree_info() const;

protected:
  std::string _fn;
  hpc::h5::file _file;
  hpc::h5::dataset _sed_ds;
  unsigned _snap;
  hpc::kdtree<real_type> _kdt;
  std::vector<unsigned> _cell_cnts;
  std::vector<unsigned long long> _cell_offs;
  hpc::h5::datatype _lc_mem_type;
  hpc::mpi::comm const *_comm;
};

template <class DomainT> class kdtree_galaxy_iterator {
public:
  typedef DomainT domain_type;

public:
  kdtree_galaxy_iterator() : _be(nullptr) {}

  kdtree_galaxy_iterator(kdtree_backend const &be)
      : _be(&be), _kdt(&be.kdtree()), _it(be.kdtree().end()) {}

  kdtree_galaxy_iterator(kdtree_backend const &be, domain_type const &dom,
                         batch<real_type> *bat, unsigned snap)
      : _be(&be), _kdt(&be.kdtree()), _dom(dom), _bat(bat), _lc(NULL), _ii(-1),
        _last_cp(-1), _inside(false), _it(be.kdtree().begin()),
        _n_ranks(be.comm().size()), _rank(be.comm().rank()), _snap(snap),
        _todo(0), _done_so_far(0) {
    _count = 0;
    _find();
    _fetch();
  }

  kdtree_galaxy_iterator(kdtree_backend const &be, domain_type const &dom,
                         const lightcone *lc, batch<real_type> *bat,
                         unsigned snap)
      : _be(&be), _kdt(&be.kdtree()), _dom(dom), _bat(bat), _lc(lc), _ii(-1),
        _last_cp(-1), _inside(false), _it(be.kdtree().begin()),
        _n_ranks(be.comm().size()), _rank(be.comm().rank()), _snap(snap),
        _todo(0), _done_so_far(0) {
    _count = 0;
    _find();
    _fetch();
  }

  herr_t static getFields(hid_t g_id, const char *name, const H5L_info_t *info,
                          void *op_data) {
    std::vector<std::string> *fields = (std::vector<std::string> *)op_data;
    std::string name_lowercase = name;
    std::transform(name_lowercase.begin(), name_lowercase.end(),
                   name_lowercase.begin(), ::tolower);
    fields->push_back(name);
    return 0;
  }

  bool operator==(kdtree_galaxy_iterator const &op) const {
    return _it == op._it;
  }

  bool operator!=(kdtree_galaxy_iterator const &op) const {
    return _it != op._it;
  }

  batch<real_type> &operator*() { return deref(); }

  batch<real_type> &deref(lightcone const *lc = nullptr) {
    ASSERT(_it != _kdt->end());
    if (!_ready) {
      _be->load_lightcone_data(_it.cell(), _lc_data);
      _copy_lightcone_data(lc);
      _ready = true;
    }
    return *_bat;
  }

  void operator++() {
    if (_todo == 0) {
      if (!done()) {
        ++_it;
        _find();
      } /*else if (_ii%_n_ranks != _rank) {
              LOGILN( "extend find");
              _find();
      }*/
    }
    _fetch();
  }

  bool done() const {
    // return !(_it != _be->kdtree().end() && _ii%_n_ranks != _rank);
    return _it == _kdt->end();
  }

  bool should_checkpoint() {
    bool res = false;
    if (_last_cp != _ii) {
      _last_cp = _ii;
      res = true;
    }
    return res;
  }

  void save_checkpoint(boost::property_tree::ptree &pt) {
    // TODO: Hmm, what goes here?
  }

  unsigned long long n_processed_gals() const {
    // TODO
    return 0;
  }

protected:
  void _copy_lightcone_data(lightcone const *lc = nullptr) {
    _bat->set_size(_lc_data.size());

    auto x = _bat->set_scalar<real_type>("Posx"); // SAGE CamelCase
    auto y = _bat->set_scalar<real_type>("Posy"); // SAGE CamelCase
    auto z = _bat->set_scalar<real_type>("Posz"); // SAGE CamelCase
    auto subsize = _bat->set_scalar<unsigned>("subsize");
    auto gidx = _bat->set_scalar<unsigned long long>("global_index");

    hpc::view<std::vector<real_type>> dist, redshift;
    if (lc) {
      dist = _bat->set_scalar<real_type>("distance");
      redshift = _bat->set_scalar<real_type>("redshift_cosmological");
    }

    for (unsigned ii = 0; ii < _lc_data.size(); ++ii) {
      if (_dom.contains(_lc_data[ii].crd, _snap)) {
        x[ii] = _lc_data[ii].crd[0] + _dom.min()[0];
        y[ii] = _lc_data[ii].crd[1] + _dom.min()[1];
        z[ii] = _lc_data[ii].crd[2] + _dom.min()[2];
        subsize[ii] = _lc_data[ii].subsize;
        gidx[ii] = _lc_data[ii].gidx;

        if (lc) {
          dist[ii] = sqrt(x[ii] * x[ii] + y[ii] * y[ii] + z[ii] * z[ii]);
          redshift[ii] = lc->distance_to_redshift(dist[ii]);
        }
      } else
        _bat->mask(ii);
    }
  }

  void _find() {
    _done_so_far = 0;
    unsigned first = _ii;
    _ready = false;
    // if (_ii==-1)
    //	LOGILN( "rank ",_rank," inits at ");
    // else
    //	LOGILN( "rank ",_rank," starts at ", _ii%_n_ranks);
    do {
      _todo = 0;
      // Each time we skip a cell because of parallel distribution
      // we need to increment the cell iterator.
      if (first != _ii) {
        ++_it;
        // LOGILN( "b.block: ", _it.cell(), " ",_ii," ", ((_ii+1)%_n_ranks),"
        // more=",(_it != _be->kdtree().end()));
      }

      while (_it != _be->kdtree().end()) {
        auto res = _dom.box_collision(_it.bounds().begin(), _snap);
        /*
        if (_it.is_leaf())
        {
            double xmin=_it.bounds()[0][0]+_dom.min()[0];
            double xmax=_it.bounds()[0][1]+_dom.min()[0];
            double ymin=_it.bounds()[1][0]+_dom.min()[1];
            double ymax=_it.bounds()[1][1]+_dom.min()[1];
            double zmin=_it.bounds()[2][0]+_dom.min()[2];
            double zmax=_it.bounds()[2][1]+_dom.min()[2];
            //std::cout<<"@"<<_snap<<"
        "<<res<<"="<<_it.bounds()[0][0]+_dom.min()[0]<<","<<_it.bounds()[0][1]+_dom.min()[0]<<
            //         ","
        <<_it.bounds()[1][0]+_dom.min()[1]<<","<<_it.bounds()[1][1]+_dom.min()[1]<<
            //
        ","<<_it.bounds()[2][0]+_dom.min()[2]<<","<<_it.bounds()[2][1]+_dom.min()[2]<<std::endl;
            std::cout<<"@
        "<<res<<","<<","<<_snap<<","<<xmin<<","<<ymin<<","<<zmin<<std::endl;
            std::cout<<"@
        "<<res<<","<<","<<_snap<<","<<xmax<<","<<ymin<<","<<zmin<<std::endl;
            std::cout<<"@
        "<<res<<","<<","<<_snap<<","<<xmax<<","<<ymax<<","<<zmin<<std::endl;
            std::cout<<"@
        "<<res<<","<<","<<_snap<<","<<xmin<<","<<ymax<<","<<zmin<<std::endl;
            std::cout<<"@
        "<<res<<","<<","<<_snap<<","<<xmin<<","<<ymax<<","<<zmax<<std::endl;
            std::cout<<"@
        "<<res<<","<<","<<_snap<<","<<xmin<<","<<ymin<<","<<zmax<<std::endl;
            std::cout<<"@
        "<<res<<","<<","<<_snap<<","<<xmax<<","<<ymin<<","<<zmax<<std::endl;
            std::cout<<"@
        "<<res<<","<<","<<_snap<<","<<xmax<<","<<ymax<<","<<zmax<<std::endl;
            std::cout<<"@
        "<<res<<","<<","<<_snap<<","<<xmin<<","<<ymax<<","<<zmax<<std::endl;
            _count++;
        }
         */
        if (!res) {
          _it.skip();
          // LOGILN( "skip: ", _it.cell() );
          continue;
        }

        /* RS: WTF?
        if( res == 2 )
           _inside = true;
        else
           */
        if (_it.is_leaf()) {
          _todo = _be->cell_counts()[_it.cell()];
          // LOGILN( "break: ", _it.cell() ," ",_todo," more=",(_it !=
          // _be->kdtree().end()));
          break;
        }
        ++_it;
        // LOGILN( "c.block: ", _it.cell(), " ",_ii," ", ((_ii+1)%_n_ranks !=
        // _rank)," more=",(_it != _be->kdtree().end()));
      }
      // LOGILN("more",((_ii+1)%_n_ranks != _rank)," more=",(_it !=
      // _be->kdtree().end()));
    } while (_it != _be->kdtree().end() && ++_ii % _n_ranks != _rank);
    // if (_todo>0)
    //     LOGILN( "found: ", _it.cell() );
    // if (_it == _be->kdtree().end() && _ii!=-1)
    //     LOGILN( "rank ",_rank," ends at ", _ii%_n_ranks);
  }
  void _fetch() {
    // Fetch actual data
    if (_todo == 0)
      return;
    auto count = _todo;
    if (_todo > _bat->max_size()) {
      count = _bat->max_size();
    }
    auto off = _be->cell_offs()[_it.cell()] + _done_so_far;

    _bat->set_size(count);

    hpc::h5::group data = _be->kdtree_file().group("data");
    std::vector<std::string> fields;
    H5Lvisit(data.id(), H5_INDEX_CRT_ORDER, H5_ITER_NATIVE, getFields, &fields);

    for (std::string field : fields) {
      // Convert HDF5 field name (CamelCase) to lowercase for batch lookup
      std::string field_lower = field;
      std::transform(field_lower.begin(), field_lower.end(),
                     field_lower.begin(), ::tolower);

      if (_bat->has_field(field_lower)) {
        _bat->field(field_lower);
        auto ds = _be->kdtree_file().dataset("data/" + field);
        switch (static_cast<tao::batch<real_type>::field_value_type>(
            _bat->get_field_type(field_lower))) {
        case tao::batch<real_type>::DOUBLE: {
          ds.read(_bat->scalar<double>(field_lower).data(),
                  hpc::h5::datatype::native_double, count, off);
          view<std::vector<double>> abc = _bat->scalar<double>(field_lower);
          double abcv = abc[0];
          // std::cout << "["<<field<<"][0]="<<abcv<<std::endl;
        } break;
        case tao::batch<real_type>::INTEGER: {
          ds.read(_bat->scalar<int>(field_lower).data(),
                  hpc::h5::datatype::native_int, count, off);
          view<std::vector<int>> abc = _bat->scalar<int>(field_lower);
          int abcv = abc[0];
          // std::cout << "[" << field << "][0]=" << abcv << std::endl;
        } break;
        case tao::batch<real_type>::LONG_LONG: {
          ds.read(_bat->scalar<long long>(field_lower).data(),
                  hpc::h5::datatype::native_llong, count, off);
          // view<std::vector<long long>> abc = _bat->scalar<long long>(field);
          // long long abcv = abc[0];
          // std::cout << "[" << field << "][0]=" << abcv << std::endl;
        } break;
        default:
          std::cout << "datatype not handled" << std::endl;
          break;
        }
      }
    }

    auto posx_ds = _be->kdtree_file().dataset("data/Posx"); // SAGE CamelCase
    auto posy_ds = _be->kdtree_file().dataset("data/Posy"); // SAGE CamelCase
    auto posz_ds = _be->kdtree_file().dataset("data/Posz"); // SAGE CamelCase
    std::vector<double> posx(count);
    std::vector<double> posy(count);
    std::vector<double> posz(count);
    posx_ds.read(posx.data(), hpc::h5::datatype::native_double, count, off);
    posy_ds.read(posy.data(), hpc::h5::datatype::native_double, count, off);
    posz_ds.read(posz.data(), hpc::h5::datatype::native_double, count, off);
    //_bat->update_size();
    std::array<real_type, 3> crd;
    for (unsigned int i = 0; i < count; i++) {
      crd[0] = posx[i];
      crd[1] = posy[i];
      crd[2] = posz[i];
      if (!_dom.contains(crd, _snap)) {
        _bat->mask(i);
        // std::cout << "mask["<<i<<"]=" << posx[i] << "," << posy[i] << "," <<
        // posz[i] << std::endl;
      }
    }

    _calc_fields();

    _todo -= count;
    _done_so_far += count;
  }

  void _calc_fields() {
    if (_lc) { // TODO Better way to decide the context (this is, avoid when
               // using box rather than light cone)

      real_type pos[3];
      real_type vel[3];

      std::array<real_type, 3> ecs_min, ecs_max, min;
      ecs_min[0] = _lc->min_ra();
      ecs_max[0] = _lc->max_ra();
      ecs_min[1] = _lc->min_dec();
      ecs_max[1] = _lc->max_dec();
      ecs_min[2] = _lc->min_dist(_snap);
      ecs_max[2] = _lc->max_dist(_snap);

      auto pos_x = _bat->template scalar<real_type>("posx");
      auto pos_y = _bat->template scalar<real_type>("posy");
      auto pos_z = _bat->template scalar<real_type>("posz");
      auto vel_x = _bat->template scalar<real_type>("velx");
      auto vel_y = _bat->template scalar<real_type>("vely");
      auto vel_z = _bat->template scalar<real_type>("velz");
      auto z_cos = _bat->template scalar<real_type>("redshift_cosmological");
      auto z_obs = _bat->template scalar<real_type>("redshift_observed");
      auto ra = _bat->template scalar<real_type>("ra");
      auto dec = _bat->template scalar<real_type>("dec");
      auto dist = _bat->template scalar<real_type>("distance");
      auto disk_sfr = _bat->template scalar<real_type>("sfrdisk");
      auto bulge_sfr = _bat->template scalar<real_type>("sfrbulge");
      auto sfr = _bat->template scalar<real_type>("sfr");
      // Note: global_index is only used in commented-out debug code
      // It may not be in data/ group if computed fields are filtered out
      // auto gid = _bat->template scalar<long long>("global_index");
      // auto galaxyindex = _bat->template scalar<long long>("galaxyindex");

      int ix = _dom.rotation()[0];
      int iy = _dom.rotation()[1];
      int iz = _dom.rotation()[2];

      // std::cout<<"+++COUNT="<<_bat->size()<<std::endl;
      auto off = _be->cell_offs()[_it.cell()] + _done_so_far;
      for (unsigned ii = 0; ii < _bat->size(); ++ii) {
        _bat->unmask(ii);
        // LOGILN( "galaxyindex: ", galaxyindex[ii],"=off=",_it.cell(),"
        // ",off+ii );

        if (_lc) {
          pos[0] = pos_x[ii];
          pos[1] = pos_y[ii];
          pos[2] = pos_z[ii];
          vel[0] = vel_x[ii];
          vel[1] = vel_y[ii];
          vel[2] = vel_z[ii];
          pos_x[ii] = pos[ix];
          pos_y[ii] = pos[iy];
          pos_z[ii] = pos[iz];
          vel_x[ii] = vel[ix];
          vel_y[ii] = vel[iy];
          vel_z[ii] = vel[iz];
          real_type h0 = _lc->simulation()->hubble();
          // apply the translation and rotation magic now

          if (pos_x[ii] + _dom.translation()[ix] <
              _lc->simulation()->box_size()) {
            pos_x[ii] += _dom.translation()[ix];
          } else {
            pos_x[ii] -= _lc->simulation()->box_size();
            pos_x[ii] += _dom.translation()[ix];
          }
          pos_x[ii] += (_dom.min()[0] - _dom.origin()[0]);

          if (pos_y[ii] + _dom.translation()[iy] <
              _lc->simulation()->box_size()) {
            pos_y[ii] += _dom.translation()[iy];
          } else {
            pos_y[ii] -= _lc->simulation()->box_size();
            pos_y[ii] += _dom.translation()[iy];
          }
          pos_y[ii] += (_dom.min()[1] - _dom.origin()[1]);

          if (pos_z[ii] + _dom.translation()[iz] <
              _lc->simulation()->box_size()) {
            pos_z[ii] += _dom.translation()[iz];
          } else {
            pos_z[ii] -= _lc->simulation()->box_size();
            pos_z[ii] += _dom.translation()[iz];
          }
          pos_z[ii] += (_dom.min()[2] - _dom.origin()[2]);

          // std::cout<<"x,y,z="<<pos_x[ii]<<","<<pos_y[ii]<<","<<pos_z[ii]<<"="<<gid[ii]<<"
          // "<<disk_sfr[ii]<<" "<<vel_x[ii]<<" "<<vel_y[ii]<<"
          // "<<vel_z[ii]<<std::endl;

          {
            min[0] = pos_x[ii];
            min[1] = pos_y[ii];
            min[2] = pos_z[ii];
            bool btest = ecs_box_collision(ecs_min, ecs_max, min, min);
            // if (galaxyindex[ii]==27044853000000652L)
            //{
            // LOGILN( "CHECK IT: ", galaxyindex[ii],"=off=",off );
            // LOGILN( "pos ", min[0]," ",min[1]," ",min[2]);
            // LOGILN( "in? ", btest);
            //}
            if (!btest) {
              _bat->mask(ii);
              continue;
            }
          }

          // std::cout<<"x,y,z="<<pos_x[ii]<<","<<pos_y[ii]<<","<<pos_z[ii]<<"="<<gid[ii]<<"
          // "<<disk_sfr[ii]<<" "<<vel_x[ii]<<" "<<vel_y[ii]<<"
          // "<<vel_z[ii]<<std::endl;

          // Compute distance.
          dist[ii] = sqrt(pos_x[ii] * pos_x[ii] + pos_y[ii] * pos_y[ii] +
                          pos_z[ii] * pos_z[ii]);
          // ASSERT(dist[ii] <= _lc->max_dist(), "Calculated distance exceeds
          // maximum: ", dist[ii]); ASSERT(dist[ii] >= _lc->min_dist(),
          // "Calculated distance below minimum: ", dist[ii]);

          // Compute cosmological redshift.
          z_cos[ii] = _lc->distance_to_redshift(dist[ii]);

          // Compute RA and DEC.
          hpc::num::cartesian_to_ecs(pos_x[ii], pos_y[ii], pos_z[ii], ra[ii],
                                     dec[ii]);

          // If the lightcone is being generated with unique cones, we may need
          // to offset the RA and DEC, then recalculate the positions.
          bool recalc = false;
          if (_lc->dec_offset() != 0.0) {
            real_type eff_ra_min = _lc->min_ra() + _lc->ra_offset();
            real_type eff_ra_max = _lc->max_ra() + _lc->ra_offset();
            real_type eff_ra_mid = 0.5 * (eff_ra_max + eff_ra_min);
            real_type eff_dec = dec[ii] - _lc->dec_offset();
            ra[ii] = (ra[ii] - eff_ra_mid) * cos(dec[ii]) / cos(eff_dec) +
                     eff_ra_mid;
            dec[ii] = eff_dec;
            recalc = true;
          }
          if (_lc->ra_offset() != 0.0) {
            ra[ii] -= _lc->ra_offset();
            recalc = true;
          }
          if (recalc)
            hpc::num::ecs_to_cartesian<real_type>(
                ra[ii], dec[ii], pos_x[ii], pos_y[ii], pos_z[ii], dist[ii]);

          // Check angles.
          // ASSERT(ra[ii] >= _dom._lc->min_ra() - EPSILON && ra[ii] <=
          // _dom._lc->max_ra() + EPSILON,
          //       "Calculated RA exceeds limits: ", to_degrees(ra[ii]));
          // ASSERT(dec[ii] >= _dom._lc->min_dec() - EPSILON && dec[ii] <=
          // _dom._lc->max_dec() + EPSILON,
          //       "Calculated RA exceeds limits: ", to_degrees(dec[ii]));

          // Return the angles in degrees.
          ra[ii] = to_degrees(ra[ii]);
          dec[ii] = to_degrees(dec[ii]);

          // Calculate observed redshift.
          z_obs[ii] = observed_redshift(
              z_cos[ii], {pos_x[ii], pos_y[ii], pos_z[ii]},
              {vel_x[ii], vel_y[ii], vel_z[ii]}, _lc->simulation()->hubble());

          // z_obs[ii] = approx_observed_redshift(
          //    *_lc,
          //    { pos_x[ii], pos_y[ii], pos_z[ii] },
          //    { vel_x[ii], vel_y[ii], vel_z[ii] }
          //    );
        } /* else if (_dom._box) {
        ASSERT(_dom._box, "Must have either lightcone or box selected.");

        z_cos[ii] = _dom._box->redshift();
        z_obs[ii] = _dom._box->redshift();

     }*/

        // Combine disk and bulge SFRs.
        sfr[ii] = disk_sfr[ii] + bulge_sfr[ii];
      }
    }
  }

protected:
  kdtree_backend const *_be;
  hpc::kdtree<real_type> const *_kdt;
  domain_type _dom;
  unsigned _snap;
  unsigned _ii;
  unsigned _last_cp;
  bool _inside;
  hpc::kdtree<real_type>::iterator _it;
  std::vector<kdtree_backend::lightcone_data> _lc_data;
  batch<real_type> *_bat;
  bool _ready;
  unsigned _n_ranks, _rank;
  unsigned _todo;
  unsigned _done_so_far;
  const lightcone *_lc;
  int _count;
};

} // namespace backends
} // namespace tao

#endif
