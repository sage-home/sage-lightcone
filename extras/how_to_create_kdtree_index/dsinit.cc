#include <iostream>
#include <vector>
#include <boost/range/algorithm/fill.hpp>
#include <libhpc/libhpc.hh>
#include <libtao/base/utils.hh>
#include <libtao/base/sage.hh>
#include "libhpc/algorithm/xdmf_writer.hh"

static unsigned long long const chunk_size = 10000;

struct sed_data_t {
    int    descendant;
    int    snapshot;
    int    local_index;
    int    merge_type;
    double dt;
    double disk_sfr;
    double bulge_sfr;
    double disk_sfr_z;
    double bulge_sfr_z;
};

struct data_permuter {
    data_permuter() : crds(nullptr), tree_idxs(nullptr), subsize(nullptr), idxs(nullptr) {
    }

    data_permuter(std::array<std::vector<double>, 3> &crds,
                  std::vector<unsigned long long> &   tree_idxs,
                  std::vector<unsigned> &             subsize,
                  std::vector<unsigned long long> &   idxs) :
        crds(&crds),
        tree_idxs(&tree_idxs), subsize(&subsize), idxs(&idxs) {
    }

    void operator()(hpc::mpi::balanced_partition const &part) {
        for(unsigned ii = 0; ii < 3; ++ii)
            part.transfer((*crds)[ii]);
        part.transfer(*tree_idxs);
        part.transfer(*subsize);
        part.transfer(*idxs);
    }

    std::array<std::vector<double>, 3> *crds;
    std::vector<unsigned long long> *   tree_idxs;
    std::vector<unsigned> *             subsize;
    std::vector<unsigned long long> *   idxs;
};

class application : public hpc::mpi::application {
  public:
    application(int argc, char *argv[]) : hpc::mpi::application(argc, argv) {
        // LOG_CONSOLE();
        LOG_PUSH(new hpc::mpi::logger("log.", hpc::log::info));

        // Setup some options.
        options().add_options()(
            "mode,m", hpc::po::value<std::string>(&_mode)->default_value("init"), "mode of operation")(
            "sage,s", hpc::po::value<hpc::fs::path>(&_sage_fn), "SAGE HDF5 file")(
            "tree,t", hpc::po::value<hpc::fs::path>(&_tree_fn), "SAGE trees file")(
            "ppc,p", hpc::po::value<unsigned>(&_ppc)->default_value(1000), "particles per cell")(
            "output,o", hpc::po::value<hpc::fs::path>(&_out_fn), "output HDF5 file");

        // Parse options.
        parse_options(argc, argv);

        // Check validity.
        hpc::to_lower(_mode);
        if(_mode == "init") {
            EXCEPT(!_sage_fn.empty(), "No SAGE filename given.");
            EXCEPT(!_tree_fn.empty(), "No SAGE trees file given.");
            EXCEPT(!_out_fn.empty(), "No output filename given.");
        } else if(_mode == "count")
            EXCEPT(!_sage_fn.empty(), "No SAGE filename given.");
    }

    void operator()() {
        if(_mode == "init")
            init();
        else if(_mode == "flatten")
            flatten();
        else if(_mode == "count")
            count();
    }

    void init() {
        LOGBLOCKI("Initialising new file.");

        sage::make_hdf5_types(_mem_type, _file_type);
        hpc::h5::file         file(_sage_fn.native(), H5F_ACC_RDONLY, hpc::mpi::comm::world);
        std::set<std::string> snap_names;
        {
            std::vector<std::string> lnks = file.links();
            for(auto const &name : lnks) {
                if(name.length() == 11 && name.compare(0, 8, "snapshot") == 0 && isdigit(name[8]) && isdigit(name[9]) &&
                   isdigit(name[10])) {
                    snap_names.insert(name);
                }
            }
            LOGILN("Snapshots to process: ", snap_names);
        }

        // Prepare the output file.
        hpc::h5::file                   out_file(_out_fn.native(), H5F_ACC_TRUNC, hpc::mpi::comm::world);
        std::vector<unsigned long long> n_gals_per_snap;
        n_gals_per_snap.reserve(snap_names.size());
        unsigned long long n_gals = 0;
        for(auto const &name : snap_names) {
            unsigned long long size = file.dataset(name).extent();
            n_gals_per_snap.push_back(size);
            n_gals += size;
        }

        // Write out snapshot counts and displacements.
        out_file.write_serial("snapshot_counts", n_gals_per_snap);
        out_file.write_serial("snapshot_displs", hpc::counts_to_displs(n_gals_per_snap));

        // Prepare output datatypes.
        {
            hsize_t         offs = 0;
            hpc::h5::derive der(3 * sizeof(double) + sizeof(unsigned long long) + sizeof(unsigned));
            der.add(hpc::h5::datatype::native_double, offs, hpc::h5::datatype::ieee_f64be, "x");
            offs += sizeof(double);
            der.add(hpc::h5::datatype::native_double, offs, hpc::h5::datatype::ieee_f64be, "y");
            offs += sizeof(double);
            der.add(hpc::h5::datatype::native_double, offs, hpc::h5::datatype::ieee_f64be, "z");
            offs += sizeof(double);
            der.add(hpc::h5::datatype::native_ullong, offs, hpc::h5::datatype::std_u64be, "global_index");
            offs += sizeof(unsigned long long);
            der.add(hpc::h5::datatype::native_uint, offs, hpc::h5::datatype::std_u32be, "subsize");
            der.commit(_lc_mem_type, _lc_file_type);
        }

        // Create coordinate dataset.
        hpc::h5::dataset lc_data_ds(out_file, "lightcone/data", _lc_file_type, n_gals);

        // Create attributes datasets.
        unsigned n_mems = _mem_type.n_members();
        for(unsigned ii = 0; ii < n_mems; ++ii) {
            std::string       name = _mem_type.member_name(ii);
            hpc::h5::datatype dt   = _mem_type.member_type(ii);
            hpc::h5::dataset  ds(out_file, "data/" + name, dt, n_gals);
        }

        // Process each snapshot.
        unsigned long long displ = 0;
        for(auto const &name : snap_names)
            process_snapshot(file, name, out_file, lc_data_ds, displ);

        // Write the SED data.
        write_sed_data(out_file);

        // Add in cosmology stuff.
        hpc::h5::copy(file, "cosmology", out_file);
        hpc::h5::copy(file, "snapshot_redshifts", out_file);
    }

    void process_snapshot(hpc::h5::file &     file,
                          std::string const & name,
                          hpc::h5::file &     out_file,
                          hpc::h5::dataset &  data,
                          unsigned long long &displ) {
        LOGBLOCKI("Processing snapshot ", name, ".");

        std::array<std::vector<double>, 3> crds;
        std::vector<unsigned long long>    tree_idxs;
        std::vector<unsigned>              subsize;
        read_coords(file.dataset(name), crds, tree_idxs, subsize);
        // read_coords( file.dataset( name ), file.dataset( name + "_map" ), crds, tree_idxs, subsize );
        // if( crds[0].size() )
        // {
        std::vector<unsigned long long> idxs(crds[0].size());
        boost::iota(idxs, 0);
        data_permuter perm(crds, tree_idxs, subsize, idxs);
        auto          part = hpc::make_binary_partitioner(crds.begin(), crds.end(), perm, _ppc);
        hpc::kdtree<> tree;
        {
            LOGBLOCKI("Creating kdtree.");
            tree.construct(part);
            LOGILN("Depth: ", tree.max_depth());
        }

        // Dump coordinates.
        {
            LOGBLOCKI("Dumping data.");
            static const char *const names[3] = {"x", "y", "z"};
            for(unsigned ii = 0; ii < 3; ++ii) {
                hpc::h5::datatype dt;
                dt.compound(sizeof(double));
                dt.insert(hpc::h5::datatype::native_double, names[ii], 0);
                hpc::h5::property_list props(H5P_DATASET_XFER);
                props.set_preserve();
                data.write(crds[ii].data(), dt, crds[0].size(), displ, hpc::mpi::comm::self, props);
            }
            {
                hpc::h5::datatype dt;
                dt.compound(sizeof(unsigned long long));
                dt.insert(hpc::h5::datatype::native_ullong, "global_index", 0);
                hpc::h5::property_list props(H5P_DATASET_XFER);
                props.set_preserve();
                data.write(tree_idxs.data(), dt, tree_idxs.size(), displ, hpc::mpi::comm::self, props);
            }
            {
                hpc::h5::datatype dt;
                dt.compound(sizeof(unsigned));
                dt.insert(hpc::h5::datatype::native_uint, "subsize", 0);
                hpc::h5::property_list props(H5P_DATASET_XFER);
                props.set_preserve();
                data.write(subsize.data(), dt, subsize.size(), displ, hpc::mpi::comm::self, props);
            }

            // Rearrange the indices and dump remaining attributes.
            std::vector<unsigned long long> inv_idxs(idxs.size());
            boost::iota(inv_idxs, 0);
            hpc::sort_by_key(idxs, inv_idxs);
            hpc::deallocate(idxs);
            hpc::h5::dataset in_ds(file, name);
            write_attributes(in_ds, out_file, inv_idxs, displ);

            // Dump kdtree.
            write_kdtree(out_file, name, tree, part, crds, displ);
        }

        displ += crds[0].size();
        // }
    }

    void write_attributes(hpc::h5::dataset &                     in_ds,
                          hpc::h5::file &                        out_file,
                          std::vector<unsigned long long> const &idxs,
                          unsigned long long                     displ) {
        unsigned           n_mems = _mem_type.n_members();
        unsigned long long n_gals = hpc::mpi::balanced_local_sizeg(in_ds.extent());
        unsigned long long offs   = hpc::mpi::comm::world.scan(n_gals);

        // Read entire snapshot into memory.
        std::vector<sage::galaxy> buf(n_gals);
        in_ds.read(buf.data(), _mem_type, buf.size(), offs);

        // Permute to snapshot order and store.
        hpc::permute(buf.begin(), buf.end(), idxs.begin());
        for(unsigned ii = 0; ii < n_mems; ++ii) {
            std::string       name    = _mem_type.member_name(ii);
            unsigned          dt_offs = _mem_type.member_offset(ii);
            hpc::h5::datatype dt      = _mem_type.member_type(ii);
            std::vector<char> tmp(buf.size() * dt.size());
            for(unsigned jj = 0; jj < buf.size(); ++jj)
                memcpy((char *)tmp.data() + dt.size() * jj,
                       (char *)buf.data() + _mem_type.size() * jj + dt_offs,
                       dt.size());
            hpc::h5::dataset   out_ds(out_file, "data/" + name);
            hpc::h5::dataspace mem_space(buf.size());
            hpc::h5::dataspace file_space(out_ds);
            file_space.select_range(displ + offs, displ + offs + buf.size());
            out_ds.write(tmp.data(), dt, mem_space, file_space, hpc::mpi::comm::self);
        }

        // std::vector<sage::galaxy> buf( chunk_size );
        // while( n_gals )
        // {
        //    buf.resize( std::min( n_gals, chunk_size ) );
        //    in_ds.read( buf.data(), _mem_type, buf.size(), offs );
        //    for( unsigned ii = 0; ii < n_mems; ++ii )
        //    {
        //       std::string name = _mem_type.member_name( ii );
        //       unsigned dt_offs = _mem_type.member_offset( ii );
        //       hpc::h5::datatype dt = _mem_type.member_type( ii );
        //       std::vector<char> tmp( buf.size()*dt.size() );
        //       for( unsigned jj = 0; jj < buf.size(); ++jj )
        //          memcpy( (char*)tmp.data() + dt.size()*jj, (char*)buf.data() + _mem_type.size()*jj + dt_offs,
        //          dt.size() );
        //       hpc::h5::dataset out_ds( out_file, "data/" + name );
        //       hpc::h5::dataspace mem_space( buf.size() );
        //       hpc::h5::dataspace file_space( out_ds );
        //       file_space.select_elements(
        //          hpc::view<std::vector<unsigned long long> >( idxs, buf.size(), offs )
        //          );
        //       // out_ds.write( tmp.data(), dt, buf.size(), offs, hpc::mpi::comm::self );
        //       out_ds.write( tmp.data(), dt, mem_space, file_space, hpc::mpi::comm::self );
        //    }
        //    offs += buf.size();
        //    n_gals -= buf.size();
        // }
    }

    void write_kdtree(hpc::h5::file &      file,
                      std::string const &  snap_name,
                      hpc::kdtree<> const &kdt,
                      hpc::binary_partitioner<std::array<std::vector<double>, 3>::iterator, data_permuter> const &bp,
                      std::array<std::vector<double>, 3> const &                                                  crds,
                      unsigned long long displ) {
        std::string name = std::string("lightcone/") + snap_name;

        // Write bounds.
        {
            hpc::h5::derive der(2 * sizeof(double));
            der.add(hpc::h5::datatype::native_double, 0, hpc::h5::datatype::ieee_f64be, "minimum");
            der.add(hpc::h5::datatype::native_double, sizeof(double), hpc::h5::datatype::ieee_f64be, "maximum");
            hpc::h5::datatype mdt, fdt;
            der.commit(mdt, fdt);

            hpc::h5::dataset dset(file, name + "/bounds", fdt, 3);
            dset.write(kdt.bounds().data(), mdt, 3, 0);
        }

        // Write splits.
        {
            hpc::h5::derive der(sizeof(hpc::kdtree<>::split_type));
            der.add(hpc::h5::datatype::native_double,
                    HOFFSET(hpc::kdtree<>::split_type, pos),
                    hpc::h5::datatype::ieee_f64be,
                    "position");
            der.add(hpc::h5::datatype::native_uint,
                    HOFFSET(hpc::kdtree<>::split_type, dim),
                    hpc::h5::datatype::std_i32be,
                    "dimension");
            hpc::h5::datatype mdt, fdt;
            der.commit(mdt, fdt);

            hpc::h5::dataset dset(file, name + "/splits", fdt, kdt.splits().size());
            dset.write(kdt.splits().data(), mdt, kdt.splits().size(), 0);
        }

        // Write counts and offsets.
        file.write(name + "/cell_counts", bp.counts());
        {
            std::vector<unsigned long long> offs(bp.offsets().size());
            for(unsigned ii = 0; ii < offs.size(); ++ii)
                offs[ii] = displ + bp.offsets()[ii];
            file.write(name + "/cell_offs", offs);
        }

        // Dump to XDMF, but only the last snapshot.
        if(snap_name == "snapshot063") {
            hpc::xdmf_writer xdmf(snap_name);

            // Dump tree.
            xdmf.open_grid(snap_name + "_tree", kdt);
            std::vector<std::array<unsigned, 1>> cells(kdt.n_leafs());
            for(unsigned ii = 0; ii < cells.size(); ++ii)
                cells[ii][0] = ii;
            xdmf.write_field("cells", cells.begin(), cells.end());
            xdmf.close_grid();

            // Dump particles.
            std::vector<std::array<double, 3>> real_crds(crds[0].size());
            for(unsigned ii = 0; ii < crds[0].size(); ++ii) {
                std::get<0>(real_crds[ii]) = crds[0][ii];
                std::get<1>(real_crds[ii]) = crds[1][ii];
                std::get<2>(real_crds[ii]) = crds[2][ii];
            }
            xdmf.open_grid(snap_name + "_particles", real_crds.begin(), real_crds.end());
            cells.resize(real_crds.size());
            for(unsigned ii = 0, jj = 0; ii < kdt.n_leafs(); ++ii) {
                unsigned cell = kdt.leaf_to_cell(ii);
                for(unsigned kk = 0; kk < bp.counts()[cell]; ++kk, ++jj)
                    cells[jj][0] = ii;
            }
            xdmf.write_field("cells", cells.begin(), cells.end());
            xdmf.close_grid();
        }
    }

    void write_sed_data(hpc::h5::file &out_file) {
        hpc::h5::datatype sed_mem_type, sed_file_type;
        hpc::h5::derive   der(sizeof(sed_data_t));
        der.add(
            hpc::h5::datatype::native_int, HOFFSET(sed_data_t, descendant), hpc::h5::datatype::std_i32be, "descendant");
        der.add(hpc::h5::datatype::native_int, HOFFSET(sed_data_t, snapshot), hpc::h5::datatype::std_i32be, "snapshot");
        der.add(hpc::h5::datatype::native_int,
                HOFFSET(sed_data_t, local_index),
                hpc::h5::datatype::std_i32be,
                "local_index");
        der.add(
            hpc::h5::datatype::native_int, HOFFSET(sed_data_t, merge_type), hpc::h5::datatype::std_i32be, "merge_type");
        der.add(hpc::h5::datatype::native_double, HOFFSET(sed_data_t, dt), hpc::h5::datatype::ieee_f64be, "dt");
        der.add(
            hpc::h5::datatype::native_double, HOFFSET(sed_data_t, disk_sfr), hpc::h5::datatype::ieee_f64be, "disk_sfr");
        der.add(hpc::h5::datatype::native_double,
                HOFFSET(sed_data_t, bulge_sfr),
                hpc::h5::datatype::ieee_f64be,
                "bulge_sfr");
        der.add(hpc::h5::datatype::native_double,
                HOFFSET(sed_data_t, disk_sfr_z),
                hpc::h5::datatype::ieee_f64be,
                "disk_sfr_z");
        der.add(hpc::h5::datatype::native_double,
                HOFFSET(sed_data_t, bulge_sfr_z),
                hpc::h5::datatype::ieee_f64be,
                "bulge_sfr_z");
        der.commit(sed_mem_type, sed_file_type);

        hpc::h5::file             in_file(_tree_fn.native(), H5F_ACC_RDONLY);
        hpc::h5::dataset          in_ds(in_file, "galaxies");
        unsigned long long        n_gals = in_ds.extent();
        hpc::h5::dataset          out_ds(out_file, "sed/data", sed_file_type, n_gals);
        std::vector<sage::galaxy> buf(chunk_size);
        std::vector<sed_data_t>   out_buf(chunk_size);
        unsigned long long        my_n_gals = hpc::mpi::balanced_local_sizeg(n_gals);
        unsigned long long        my_displ  = hpc::mpi::comm::world.scan(my_n_gals);

        while(my_n_gals) {
            buf.resize(std::min(my_n_gals, chunk_size));
            in_ds.read(buf.data(), _mem_type, buf.size(), my_displ);
            for(unsigned ii = 0; ii < buf.size(); ++ii) {
                out_buf[ii].descendant  = buf[ii].descendant;
                out_buf[ii].snapshot    = buf[ii].snapshot;
                out_buf[ii].local_index = buf[ii].local_index;
                out_buf[ii].merge_type  = buf[ii].merge_type;
                out_buf[ii].dt          = buf[ii].dt;
                out_buf[ii].disk_sfr    = buf[ii].sfr_disk;
                out_buf[ii].bulge_sfr   = buf[ii].sfr_bulge;
                out_buf[ii].disk_sfr_z  = buf[ii].sfr_disk_z;
                out_buf[ii].bulge_sfr_z = buf[ii].sfr_bulge_z;
            }

            out_ds.write(out_buf.data(), sed_mem_type, buf.size(), my_displ);

            my_displ += buf.size();
            my_n_gals -= buf.size();
        }
    }

    void flatten() {
        LOGBLOCKI("Flattening.");

        sage::make_hdf5_types(_mem_type, _file_type);
        hpc::h5::file                   file(_sage_fn.native(), H5F_ACC_RDONLY, hpc::mpi::comm::world);
        std::vector<unsigned long long> gal_cnts   = count_galaxies(file);
        std::vector<unsigned long long> gal_displs = hpc::mpi::comm::world.scan(gal_cnts);
        LOGDLN("Counts: ", gal_cnts);
        LOGDLN("Displs: ", gal_displs);

        hpc::h5::file                 out_file(_out_fn.native(), H5F_ACC_TRUNC, hpc::mpi::comm::world);
        unsigned                      n_snaps = gal_cnts.size();
        std::vector<hpc::h5::dataset> dsets(n_snaps);
        // std::vector<hpc::h5::dataset> ord_dsets( n_snaps );
        for(unsigned ii = 0; ii < n_snaps; ++ii) {
            unsigned long long gsize = hpc::mpi::comm::world.all_reduce(gal_cnts[ii]);
            dsets[ii].create(out_file, std::string("snapshot") + hpc::index_string(ii, 3), _file_type, gsize);
            // ord_dsets[ii].create( out_file, std::string( "snapshot" ) + hpc::index_string( ii, 3 ) + "_map",
            //                       hpc::h5::datatype::native_ullong, gsize );
        }

        LOGILN("Flattening ", n_snaps, " snapshots.");
        hpc::h5::dataset          gals_ds(file, "galaxies");
        unsigned long long        n_gals = gals_ds.extent();
        std::vector<sage::galaxy> buf(chunk_size);
        hpc::matrix<sage::galaxy> out_buf(n_snaps, chunk_size);
        // hpc::matrix<unsigned long long> ord_out_buf( n_snaps, chunk_size );
        std::vector<unsigned> out_buf_sizes(n_snaps);
        unsigned long long    my_n_gals = hpc::mpi::balanced_local_sizeg(n_gals);
        unsigned long long    my_displ  = hpc::mpi::comm::world.scan(my_n_gals);
        LOGDLN("Processing ", my_n_gals, " galaxies.");
        LOGDLN("Starting at ", my_displ, ".");

        while(my_n_gals) {
            buf.resize(std::min(my_n_gals, chunk_size));
            LOGILN("Writing chunk from ", my_displ, " to ", my_displ + buf.size(), " (", my_n_gals, " remaining)");
            gals_ds.read(buf.data(), _mem_type, buf.size(), my_displ);
            boost::fill(out_buf_sizes, 0);
            for(unsigned ii = 0; ii < buf.size(); ++ii) {
                out_buf(buf[ii].snapshot, out_buf_sizes[buf[ii].snapshot]) = buf[ii];
                // ord_out_buf( buf[ii].snapshot, out_buf_sizes[buf[ii].snapshot] ) = my_displ + ii;
                ++out_buf_sizes[buf[ii].snapshot];
            }
            my_displ += buf.size();
            my_n_gals -= buf.size();

            for(unsigned ii = 0; ii < n_snaps; ++ii) {
                if(out_buf_sizes[ii]) {
                    dsets[ii].write(out_buf[ii].data(), _mem_type, out_buf_sizes[ii], gal_displs[ii]);
                    // ord_dsets[ii].write( ord_out_buf[ii].data(), hpc::h5::datatype::native_ullong, out_buf_sizes[ii],
                    // gal_displs[ii] );
                    gal_displs[ii] += out_buf_sizes[ii];
                }
            }
        }

        LOGDLN("Writing auxilliary data.");
        hpc::h5::copy(file, "cosmology", out_file);
        hpc::h5::copy(file, "snapshot_redshifts", out_file);

        hpc::mpi::comm::world.barrier();
    }

    void read_coords(hpc::h5::dataset const &dset,
                     // hpc::h5::dataset const& tree_idxs_dset,
                     std::array<std::vector<double>, 3> &crds,
                     std::vector<unsigned long long> &   tree_idxs,
                     std::vector<unsigned> &             subsize) {
        LOGBLOCKI("Reading coordinates.");

        unsigned long long n_crds = hpc::mpi::balanced_local_sizeg(dset.extent());
        unsigned long long offs   = hpc::mpi::comm::world.scan(n_crds);
        for(unsigned ii = 0; ii < 3; ++ii) {
            crds[ii].reserve(n_crds);
            crds[ii].resize(0);
        }
        tree_idxs.resize(n_crds);
        subsize.resize(n_crds);
        std::vector<sage::galaxy> buf(chunk_size);
        while(n_crds) {
            buf.resize(std::min(n_crds, chunk_size));
            dset.read(buf.data(), _mem_type, buf.size(), offs);
            for(unsigned ii = 0; ii < buf.size(); ++ii) {
                for(unsigned jj = 0; jj < 3; ++jj)
                    crds[jj].push_back(buf[ii].pos[jj]);
                subsize[offs + ii]   = buf[ii].subsize;
                tree_idxs[offs + ii] = buf[ii].global_index;
                ASSERT(subsize[ii] > 0, "Subsize can never be zero.");
            }
            offs += buf.size();
            n_crds -= buf.size();
        }
        // tree_idxs_dset.read( tree_idxs );
        LOGILN("Read ", crds[0].size(), " coordinates.");
    }

    std::vector<unsigned long long> count_galaxies(hpc::h5::file &file) {
        LOGBLOCKI("Counting galaxies.");
        std::vector<unsigned long long> cnts(file.dataset("snapshot_redshifts").extent());
        boost::fill(cnts, 0);
        hpc::h5::dataset          gals_ds(file, "galaxies");
        unsigned long long        n_gals = gals_ds.extent();
        std::vector<sage::galaxy> buf(chunk_size);
        unsigned long long        my_n_gals = hpc::mpi::balanced_local_sizeg(n_gals);
        unsigned long long        my_displ  = hpc::mpi::comm::world.scan(my_n_gals);
        while(my_n_gals) {
            buf.resize(std::min(my_n_gals, chunk_size));
            LOGILN("Counting chunk from ", my_displ, " to ", my_displ + buf.size(), " (", my_n_gals, " remaining)");
            gals_ds.read(buf.data(), _mem_type, buf.size(), my_displ);
            for(unsigned ii = 0; ii < buf.size(); ++ii)
                ++cnts[buf[ii].snapshot];
            my_displ += buf.size();
            my_n_gals -= buf.size();
        }
        return cnts;
    }

    void count() {
        hpc::h5::file    file(_sage_fn.native(), H5F_ACC_RDONLY, hpc::mpi::comm::world);
        hpc::h5::dataset gals_dset(file, "galaxies");
        std::cout << gals_dset.extent() << " galaxies\n";
        std::cout << gals_dset.extent() * sizeof(double) * 3 << " bytes\n";
    }

  protected:
    std::string       _mode;
    hpc::fs::path     _sage_fn;
    hpc::fs::path     _out_fn;
    hpc::fs::path     _tree_fn;
    hpc::h5::datatype _mem_type, _file_type;
    hpc::h5::datatype _lc_mem_type, _lc_file_type;
    unsigned          _ppc;
};

#define HPC_APP_CLASS application
#include <libhpc/mpi/main.hh>
