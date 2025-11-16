/**
 *  @file    dstreeinit.cc
 *  @author  Ray Seikel (rseikel@bigpond.com)
 *  @date    30/7/2018
 *
 *  @brief General utility and testing program used for KDTREE implementation
 *
 *  @section DESCRIPTION
 *
 *
 *  Command line modes:
 *  init - actually creates the KDTREE indexed HDF5 file (requires the original galaxyTree HDF5
 *  so it can pick out the columns needed for SED and put them in a "sed/data" group in such a way that
 *  globaltreeid can be used as an index to access)
 *  flatten - legacy
 *  tree2sage - convert TAOImport output hdf5 file with tree structure to KDTree indexed hdf5
 *  dump - legacy
 *  count - legacy
 *
 *
 */


#include <iostream>
#include <vector>
#include <boost/range/algorithm/fill.hpp>
#include <pugixml.hpp>
#include <libhpc/libhpc.hh>
#include <libtao/base/utils.hh>
#include <libtao/base/sage.hh>
#include <libtao/base/xml_dict.hh>
#include <libtao/base/data_dict.hh>
#include <libtao/base/batch.hh>
#include "libhpc/algorithm/xdmf_writer.hh"
#include "H5LTpublic.h"

using namespace ::tao;
using namespace ::hpc;
using namespace ::sage;

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

class myapplication : public hpc::mpi::application {
  public:
    myapplication(int argc, char *argv[]) : hpc::mpi::application(argc, argv) {
        // LOG_CONSOLE();
        LOG_PUSH(new hpc::mpi::logger("log.", hpc::log::info));

        // Setup some options.
        options().add_options()(
            "mode,m", hpc::po::value<std::string>(&_mode)->default_value("init"), "mode of operation")(
            "sage,s", hpc::po::value<hpc::fs::path>(&_sage_fn), "SAGE HDF5 file")(
            "tree,t", hpc::po::value<hpc::fs::path>(&_tree_fn), "SAGE trees file")(
            "txt,i", hpc::po::value<hpc::fs::path>(&_txt_fn), "txt file")(
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
        } else if (_mode == "tree2sage") {
            EXCEPT(!_tree_fn.empty(), "No SAGE trees file given.");
            EXCEPT(!_out_fn.empty(), "No output filename given.");
        } else if (_mode == "dump") {
            EXCEPT(!_tree_fn.empty(), "No SAGE trees file given.");
        } else if(_mode == "count")
            EXCEPT(!_sage_fn.empty(), "No SAGE filename given.");
    }

    void operator()() {
        if(_mode == "init")
            initANY();
        else if(_mode == "flatten")
            flatten();
        else if(_mode == "tree2sage")
            tree2sageANY();
        else if(_mode == "dump")
            dump();
        else if(_mode == "count")
            count();
    }

    void init() {
        LOGBLOCKI("Initialising new file.");

        sage::make_hdf5_types_treesed(_mem_type, _file_type);
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
            std::string       name_lowercase = name;
            std::transform(name_lowercase.begin(), name_lowercase.end(), name_lowercase.begin(), ::tolower);
            hpc::h5::datatype dt   = _mem_type.member_type(ii);
            hpc::h5::dataset  ds(out_file, "data/" + name_lowercase, dt, n_gals);
            std::cout << name<< " bytes="<<dt.size()<<std::endl;
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
    void initANY() {
            LOGBLOCKI("Initialising new file.");

            tao::xml_dict sidecar_xml = tao::data_dict::getSidecar(_sage_fn.string(),"/sageinput");
            _fields = tao::data_dict::getFieldsANY(sidecar_xml,"/sageinput/Field");
            size_t nfields=_fields.size();
            hpc::h5::derive der = tao::data_dict::getCompositeANY(_fields, _mem_type, _file_type);

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
                std::string       name_lowercase = name;
                std::transform(name_lowercase.begin(), name_lowercase.end(), name_lowercase.begin(), ::tolower);
                hpc::h5::datatype dt   = _mem_type.member_type(ii);
                hpc::h5::dataset  ds(out_file, "data/" + name_lowercase, dt, n_gals);
            }

            // Process each snapshot.
            unsigned long long displ = 0;
            for(auto const &name : snap_names) {
                std::cout << "Snapshop:"<< name<<std::endl;
                process_snapshotANY(file, name, out_file, lc_data_ds, displ);
            }

            // Write the SED data.
            write_sed_dataANY(out_file);

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

    void process_snapshotANY(hpc::h5::file &     file,
                              std::string const & name,
                              hpc::h5::file &     out_file,
                              hpc::h5::dataset &  data,
                              unsigned long long &displ) {
            LOGBLOCKI("Processing snapshot ", name, ".");
            std::cout << "Processing snapshot "<<name<<std::endl;

            std::array<std::vector<double>, 3> crds;
            std::vector<unsigned long long>    tree_idxs;
            std::vector<unsigned>              subsize;
            read_coordsANY(file.dataset(name), crds, tree_idxs, subsize);
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
                write_attributesANY(in_ds, out_file, inv_idxs, displ);

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
        std::vector<sage::galaxyTreeSed> buf(n_gals);
        in_ds.read(buf.data(), _mem_type, buf.size(), offs);

        // Permute to snapshot order and store.
        hpc::permute(buf.begin(), buf.end(), idxs.begin());
        for(unsigned ii = 0; ii < n_mems; ++ii) {
            std::string       name    = _mem_type.member_name(ii);
            std::string       name_lowercase = name;
            std::transform(name_lowercase.begin(), name_lowercase.end(), name_lowercase.begin(), ::tolower);
            unsigned          dt_offs = _mem_type.member_offset(ii);
            hpc::h5::datatype dt      = _mem_type.member_type(ii);
            std::vector<char> tmp(buf.size() * dt.size());
            for(unsigned jj = 0; jj < buf.size(); ++jj)
                memcpy((char *)tmp.data() + dt.size() * jj,
                       (char *)buf.data() + _mem_type.size() * jj + dt_offs,
                       dt.size());
            hpc::h5::dataset   out_ds(out_file, "data/" + name_lowercase);
            hpc::h5::dataspace mem_space(buf.size());
            hpc::h5::dataspace file_space(out_ds);
            file_space.select_range(displ + offs, displ + offs + buf.size());
            out_ds.write(tmp.data(), dt, mem_space, file_space, hpc::mpi::comm::self);
        }

    }

    void write_attributesANY(hpc::h5::dataset &                     in_ds,
                              hpc::h5::file &                        out_file,
                              std::vector<unsigned long long> const &idxs,
                              unsigned long long                     displ) {
            unsigned           n_mems = _mem_type.n_members();
            unsigned long long n_gals = hpc::mpi::balanced_local_sizeg(in_ds.extent());
            unsigned long long offs   = hpc::mpi::comm::world.scan(n_gals);

            // Read entire snapshot into memory.
            size_t todo=n_gals;
            std::vector<tao::ANY> buf(n_gals*n_mems);
            in_ds.read(buf.data(), _mem_type, todo, offs);

            // Permute to snapshot order and store.
            hpc::permuteANY(buf.begin(), buf.end(), idxs.begin(), n_mems); // TODO: RS: I don't understand this but looks like it needs to change
            for(unsigned ii = 0; ii < n_mems; ++ii) {
                std::string       name    = _mem_type.member_name(ii);
                std::string       name_lowercase = name;
                std::transform(name_lowercase.begin(), name_lowercase.end(), name_lowercase.begin(), ::tolower);
                unsigned          dt_offs = _mem_type.member_offset(ii);
                hpc::h5::datatype dt      = _mem_type.member_type(ii);
                std::vector<char> tmp(todo * dt.size());
                for(unsigned jj = 0; jj < todo; ++jj)
                    memcpy((char *)tmp.data() + dt.size() * jj,
                           (char *)buf.data() + _mem_type.size() * jj + dt_offs,
                           dt.size());
                hpc::h5::dataset   out_ds(out_file, "data/" + name_lowercase);
                hpc::h5::dataspace mem_space(todo);
                hpc::h5::dataspace file_space(out_ds);
                file_space.select_range(displ + offs, displ + offs + todo);
                out_ds.write(tmp.data(), dt, mem_space, file_space, hpc::mpi::comm::self);
            }
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
        std::vector<sage::galaxyTreeSed> buf(chunk_size);
        std::vector<sed_data_t>   out_buf(chunk_size);
        unsigned long long        my_n_gals = hpc::mpi::balanced_local_sizeg(n_gals);
        unsigned long long        my_displ  = hpc::mpi::comm::world.scan(my_n_gals);

        while(my_n_gals) {
            buf.resize(std::min(my_n_gals, chunk_size));
            in_ds.read(buf.data(), _mem_type, buf.size(), my_displ);
            for(unsigned ii = 0; ii < buf.size(); ++ii) {
                out_buf[ii].descendant  = buf[ii].descendant;
                out_buf[ii].snapshot    = buf[ii].snapnum;
                out_buf[ii].local_index = buf[ii].localindex;
                out_buf[ii].merge_type  = buf[ii].mergetype;
                out_buf[ii].dt          = buf[ii].dt;
                out_buf[ii].disk_sfr    = buf[ii].sfrdisk;
                out_buf[ii].bulge_sfr   = buf[ii].sfrbulge;
                out_buf[ii].disk_sfr_z  = buf[ii].sfrdiskz;
                out_buf[ii].bulge_sfr_z = buf[ii].sfrbulgez;
            }

            out_ds.write(out_buf.data(), sed_mem_type, buf.size(), my_displ);

            my_displ += buf.size();
            my_n_gals -= buf.size();
        }
    }

    void write_sed_dataANY(hpc::h5::file &out_file) {
            unsigned n_mems = _mem_type.n_members();
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
            std::vector<tao::ANY> buf(chunk_size*n_mems);
            std::vector<sed_data_t>   out_buf(chunk_size);
            unsigned long long        my_n_gals = hpc::mpi::balanced_local_sizeg(n_gals);
            unsigned long long        my_displ  = hpc::mpi::comm::world.scan(my_n_gals);

            uint16_t descendant_index=tao::data_dict::findField(_fields,"descendant");
            uint16_t snapnum_index=tao::data_dict::findField(_fields,"snapnum");
            uint16_t localindex_index=tao::data_dict::findField(_fields,"localindex");
            uint16_t mergetype_index=tao::data_dict::findField(_fields,"mergetype");
            uint16_t dt_index=tao::data_dict::findField(_fields,"dt");
            uint16_t sfrdisk_index=tao::data_dict::findField(_fields,"sfrdisk");
            uint16_t sfrbulge_index=tao::data_dict::findField(_fields,"sfrbulge");
            uint16_t sfrdiskz_index=tao::data_dict::findField(_fields,"sfrdiskz");
            uint16_t sfrbulgez_index=tao::data_dict::findField(_fields,"sfrbulgez");

            while(my_n_gals) {
                size_t thistime = std::min(my_n_gals, chunk_size);
                buf.resize(thistime*n_mems);
                in_ds.read(buf.data(), _mem_type, thistime, my_displ);
                for(unsigned ii = 0; ii < thistime; ++ii) {
                    out_buf[ii].descendant  = buf[ii*n_mems + descendant_index].i4;
                    out_buf[ii].snapshot    = buf[ii*n_mems + snapnum_index].i4;
                    out_buf[ii].local_index = buf[ii*n_mems + localindex_index].i4;
                    out_buf[ii].merge_type  = buf[ii*n_mems + mergetype_index].i4;
                    out_buf[ii].dt          = buf[ii*n_mems + dt_index].f4;
                    out_buf[ii].disk_sfr    = buf[ii*n_mems + sfrdisk_index].f4;
                    out_buf[ii].bulge_sfr   = buf[ii*n_mems + sfrbulge_index].f4;
                    out_buf[ii].disk_sfr_z  = buf[ii*n_mems + sfrdiskz_index].f4;
                    out_buf[ii].bulge_sfr_z = buf[ii*n_mems + sfrbulgez_index].f4;
                }

                out_ds.write(out_buf.data(), sed_mem_type, thistime, my_displ);

                my_displ += thistime;
                my_n_gals -= thistime;
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
        std::vector<sage::galaxyTreeSed> buf(chunk_size);
        hpc::matrix<sage::galaxyTreeSed> out_buf(n_snaps, chunk_size);
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
                out_buf(buf[ii].snapnum, out_buf_sizes[buf[ii].snapnum]) = buf[ii];
                // ord_out_buf( buf[ii].snapshot, out_buf_sizes[buf[ii].snapshot] ) = my_displ + ii;
                ++out_buf_sizes[buf[ii].snapnum];
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
        std::vector<sage::galaxyTreeSed> buf(chunk_size);
        while(n_crds) {
            buf.resize(std::min(n_crds, chunk_size));
            dset.read(buf.data(), _mem_type, buf.size(), offs);
            for(unsigned ii = 0; ii < buf.size(); ++ii) {
                for(unsigned jj = 0; jj < 3; ++jj)
                    crds[jj].push_back(buf[ii].pos[jj]);
                subsize[offs + ii]   = buf[ii].subsize;
                tree_idxs[offs + ii] = buf[ii].globalindex;
                ASSERT(subsize[offs + ii] > 0, "Subsize can never be zero.");
            }
            offs += buf.size();
            n_crds -= buf.size();
        }
        // tree_idxs_dset.read( tree_idxs );
        LOGILN("Read ", crds[0].size(), " coordinates.");
    }
    void read_coordsANY(hpc::h5::dataset const &dset,
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
            unsigned           n_mems = _mem_type.n_members();

            uint16_t posx_index = tao::data_dict::findField(_fields, "posx");
            uint16_t posy_index = tao::data_dict::findField(_fields, "posy");
            uint16_t posz_index = tao::data_dict::findField(_fields, "posz");
            uint16_t subsize_index = tao::data_dict::findField(_fields, "subsize");
            uint16_t globalindex_index = tao::data_dict::findField(_fields, "global_index");
            std::vector<tao::ANY> buf(chunk_size*n_mems);
            while(n_crds) {
                size_t thistime = std::min(n_crds, chunk_size);
                buf.resize(thistime*n_mems);
                dset.read(buf.data(), _mem_type, thistime, offs);
                for(unsigned ii = 0; ii < thistime; ++ii) {
                    crds[0].push_back(buf[ii*n_mems + posx_index].f4);
                    crds[1].push_back(buf[ii*n_mems + posy_index].f4);
                    crds[2].push_back(buf[ii*n_mems + posz_index].f4);
                    subsize[offs + ii]   = buf[ii*n_mems + subsize_index].i8;
                    tree_idxs[offs + ii] = buf[ii*n_mems + globalindex_index].i8;
                    //ASSERT(subsize[offs + ii] > 0, "Subsize can never be zero.");
                }
                offs += thistime;
                n_crds -= thistime;
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
        std::vector<sage::galaxyTreeSed> buf(chunk_size);
        unsigned long long        my_n_gals = hpc::mpi::balanced_local_sizeg(n_gals);
        LOGILN("my_n_gals=",my_n_gals);
        unsigned long long        my_displ  = hpc::mpi::comm::world.scan(my_n_gals);
        LOGILN("my_displ=",my_displ);
        while(my_n_gals) {
            buf.resize(std::min(my_n_gals, chunk_size));
            LOGILN("Counting chunk from ", my_displ, " to ", my_displ + buf.size(), " (", my_n_gals, " remaining)");
            gals_ds.read(buf.data(), _mem_type, buf.size(), my_displ);
            for(unsigned ii = 0; ii < buf.size(); ++ii)
                ++cnts[buf[ii].snapnum];
            my_displ += buf.size();
            my_n_gals -= buf.size();
        }
        return cnts;
    }
    std::vector<unsigned long long> count_galaxiesANY(hpc::h5::file &file, int nfields, int snapnum_index) {
            LOGBLOCKI("Counting galaxies.");
            std::vector<unsigned long long> cnts(file.dataset("snapshot_redshifts").extent());
            boost::fill(cnts, 0);
            hpc::h5::dataset          gals_ds(file, "galaxies");
            unsigned long long        n_gals = gals_ds.extent();
            std::vector<sage::galaxyTreeSed> buf(chunk_size);
            std::vector<tao::ANY> ubuf(chunk_size * nfields);
            unsigned long long        my_n_gals = hpc::mpi::balanced_local_sizeg(n_gals);
            LOGILN("my_n_gals=",my_n_gals);
            unsigned long long        my_displ  = hpc::mpi::comm::world.scan(my_n_gals);
            LOGILN("my_displ=",my_displ);
            while(my_n_gals) {
                int todo=std::min(my_n_gals, chunk_size);
                buf.resize(todo);
                ubuf.resize(todo * nfields);
                LOGILN("Counting chunk from ", my_displ, " to ", my_displ + buf.size(), " (", my_n_gals, " remaining)");
                gals_ds.read(ubuf.data(), _mem_type, todo, my_displ);
                for(unsigned ii = 0; ii < buf.size(); ++ii) {
                    int32_t snapnum = ubuf[nfields*ii+snapnum_index].i4;
                    ++cnts[snapnum];
                }
                my_displ += buf.size();
                my_n_gals -= buf.size();
            }
            return cnts;
        }

    void tree2sage() {
        LOGBLOCKI("Convert tree h5 file to sage hd5 file suitable for input to init");
        // Work out the columns
        tao::xml_dict sidecar_xml = tao::data_dict::getSidecar(_tree_fn.native(), "/sageinput");
        std::vector<tao::data_dict_field> fields = tao::data_dict::getFieldsANY(sidecar_xml, "/sageinput/Field");
        std::unordered_map<std::string,std::string> _field_map;
        std::map<std::string,typename tao::batch<tao::real_type>::field_value_type> _field_types;
        {

            tao::xml_dict xml;
            xml.read("settings.xml", "/settings");
            // Iterate over the module nodes.
            pugi::xpath_node_set fields = xml.get_nodes("/settings/sageinput/Field");
            std::map<std::string, int> _field_order;

            _field_types.clear();

            int order=0;
            for (const pugi::xpath_node *it = fields.begin(); it != fields.end(); ++it) {
                pugi::xml_node cur = it->node();

                std::string element = cur.name();
                std::string field_str = cur.text().as_string();
                std::transform(field_str.begin(), field_str.end(), field_str.begin(), ::tolower);
                std::string label = cur.attribute("label").value();
                std::string description = cur.attribute("description").value();
                //int order = cur.attribute("order").as_int();
                std::string units = cur.attribute("units").value();
                std::string group = cur.attribute("group").value();
                std::string type_str = cur.attribute("Type").value();
                std::cout << "Field[" << field_str << "]=" << label << "," << description << std::endl;
                typename tao::batch<tao::real_type>::field_value_type type;
                if (type_str == "int" || type_str == "short")
                    type = tao::batch<tao::real_type>::INTEGER;
                else if (type_str == "long long")
                    type = tao::batch<tao::real_type>::LONG_LONG;
                else if (type_str == "float")
                    type = tao::batch<tao::real_type>::DOUBLE;
                else {
                    EXCEPT(0, "Unknown field type for field '", field_str, "': ", type_str);
                }
                _field_types.emplace(field_str, type);
                _field_map[field_str] = field_str;
                _field_order[field_str] = order;
                order++;

            }

            //=========================
            // Setup the dictionary of column names versus the galaxy structure
            class Entry {
            public:
                Entry()
                {

                }
                Entry(std::string nt,int ioffset,std::string mt, std::string ilabel)
                {
                    native_type = nt;
                    offset = ioffset;
                    mem_type = mt;
                    label = ilabel;
                }
                Entry(const Entry& entry)
                {
                    native_type = entry.native_type;
                    offset = entry.offset;
                    mem_type = entry.mem_type;
                    label = entry.label;
                }
                std::string native_type;
                int offset;
                std::string mem_type;
                std::string label;
            };
            std::unordered_map<std::string, Entry> field_dict;

            field_dict["snapshot"] = Entry("h5::datatype::native_int",   HOFFSET( galaxy, snapshot ),            "h5::datatype::std_i32be",  "snapshot" );
            field_dict["type"] = Entry("h5::datatype::native_int",   HOFFSET( galaxy, type ),                "h5::datatype::std_i32be",  "type" );
            field_dict["galaxyindex"] = Entry("h5::datatype::native_llong", HOFFSET( galaxy, galaxy_idx ),          "h5::datatype::std_i64be",  "galaxy_index" );
            field_dict["centralgalaxyindex"] = {"h5::datatype::native_llong", HOFFSET( galaxy, central_galaxy_idx ),  "h5::datatype::std_i64be",  "central_galaxy_index"};
            field_dict["sagehaloindex"] = {"h5::datatype::native_int",   HOFFSET( galaxy, sage_halo_idx ),       "h5::datatype::std_i32be",  "sage_halo_index"};
            field_dict["sagetreeindex"] = {"h5::datatype::native_int",   HOFFSET( galaxy, sage_tree_idx ),       "h5::datatype::std_i32be",  "sage_tree_index" };
            field_dict["simulationhaloindex"] = {"h5::datatype::native_llong", HOFFSET( galaxy, simulation_halo_idx ), "h5::datatype::std_i32be",  "simulation_halo_index" };
            field_dict["localindex"] = {"h5::datatype::native_int",   HOFFSET( galaxy, local_index ),         "h5::datatype::std_i32be",  "local_index" };
            field_dict["globalindex"] = {"h5::datatype::native_llong", HOFFSET( galaxy, global_index ),        "h5::datatype::std_i64be",  "global_index" };
            field_dict["descendant"] = {"h5::datatype::native_int",   HOFFSET( galaxy, descendant ),          "h5::datatype::std_i32be",  "descendant" };

            field_dict["globaldescendant"] = {"h5::datatype::native_llong", HOFFSET( galaxy, global_descendant ),   "h5::datatype::std_i64be",  "global_descendant" };
            field_dict["subsize"] = {"h5::datatype::native_int",   HOFFSET( galaxy, subsize ),             "h5::datatype::std_i32be",  "subsize" };
            field_dict["mergetype"] = {"h5::datatype::native_int",   HOFFSET( galaxy, merge_type ),          "h5::datatype::std_i32be",  "merge_type"};
            field_dict["mergeintoid"] = {"h5::datatype::native_int",   HOFFSET( galaxy, merge_into_id ),       "h5::datatype::std_i32be",  "merge_into_id"};
            field_dict["mergeintosnapshot"] = {"h5::datatype::native_int",   HOFFSET( galaxy, merge_into_snapshot ), "h5::datatype::std_i32be",  "merge_into_snapshot" };
            field_dict["dt"] = {"h5::datatype::native_float", HOFFSET( galaxy, dt ),                  "h5::datatype::ieee_f32be", "dt" };
            field_dict["positionx"] = { "h5::datatype::native_float", HOFFSET( galaxy, pos[0] ),              "h5::datatype::ieee_f32be", "position_x"};
            field_dict["positiony"] = { "h5::datatype::native_float", HOFFSET( galaxy, pos[1] ),              "h5::datatype::ieee_f32be", "position_y" };
            field_dict["positionz"] = { "h5::datatype::native_float", HOFFSET( galaxy, pos[2] ),              "h5::datatype::ieee_f32be", "position_z" };
            field_dict["velocityx"] = { "h5::datatype::native_float", HOFFSET( galaxy, vel[0] ),              "h5::datatype::ieee_f32be", "velocity_x" };
            field_dict["velocityy"] = { "h5::datatype::native_float", HOFFSET( galaxy, vel[1] ),              "h5::datatype::ieee_f32be", "velocity_y" };
            field_dict["velocityz"] = { "h5::datatype::native_float", HOFFSET( galaxy, vel[2] ),              "h5::datatype::ieee_f32be", "velocity_z" };
            field_dict["spinx"] = { "h5::datatype::native_float", HOFFSET( galaxy, spin[0] ),             "h5::datatype::ieee_f32be", "spin_x" };
            field_dict["spiny"] = { "h5::datatype::native_float", HOFFSET( galaxy, spin[1] ),             "h5::datatype::ieee_f32be", "spin_y" };
            field_dict["spinz"] = { "h5::datatype::native_float", HOFFSET( galaxy, spin[2] ),             "h5::datatype::ieee_f32be", "spin_z" };
            field_dict["ndarkmatterparticles"] = {  "h5::datatype::native_int",   HOFFSET( galaxy, num_particles ),       "h5::datatype::std_i32be",  "n_darkmatter_particles" };
            field_dict["virialmass"] = {  "h5::datatype::native_float", HOFFSET( galaxy, mvir ),                "h5::datatype::ieee_f32be", "virial_mass" };
            field_dict["centralgalaxymvir"] = { "h5::datatype::native_float", HOFFSET( galaxy, central_mvir ),        "h5::datatype::ieee_f32be", "central_galaxy_mvir" };
            field_dict["virialradius"] = {  "h5::datatype::native_float", HOFFSET( galaxy, rvir ),                "h5::datatype::ieee_f32be", "virial_radius" };
            field_dict["virialvelocity"] = {  "h5::datatype::native_float", HOFFSET( galaxy, vvir ),                "h5::datatype::ieee_f32be", "virial_velocity" };
            field_dict["maxvelocity"] = {  "h5::datatype::native_float", HOFFSET( galaxy, vmax ),                "h5::datatype::ieee_f32be", "max_velocity" };
            field_dict["velocitydispersion"] = {  "h5::datatype::native_float", HOFFSET( galaxy, vel_disp ),            "h5::datatype::ieee_f32be", "velocity_dispersion" };
            field_dict["coldgas"] = {  "h5::datatype::native_float", HOFFSET( galaxy, cold_gas ),            "h5::datatype::ieee_f32be", "cold_gas" };
            field_dict["stellarmass"] = {  "h5::datatype::native_float", HOFFSET( galaxy, stellar_mass ),        "h5::datatype::ieee_f32be", "stellar_mass" };
            field_dict["bulgemass"] = {  "h5::datatype::native_float", HOFFSET( galaxy, bulge_mass ),          "h5::datatype::ieee_f32be", "bulge_mass" };
            field_dict["hotgas"] = { "h5::datatype::native_float", HOFFSET( galaxy, hot_gas ),             "h5::datatype::ieee_f32be", "hot_gas" };
            field_dict["ejectedmass"] = { "h5::datatype::native_float", HOFFSET( galaxy, ejected_mass ),        "h5::datatype::ieee_f32be", "ejected_mass" };
            field_dict["blackholemass"] = { "h5::datatype::native_float", HOFFSET( galaxy, blackhole_mass ),      "h5::datatype::ieee_f32be", "blackhole_mass" };
            field_dict["ics"] = {"h5::datatype::native_float", HOFFSET( galaxy, ics ),                 "h5::datatype::ieee_f32be", "ics" };
            field_dict["metalscoldgas"] = { "h5::datatype::native_float", HOFFSET( galaxy, metals_cold_gas ),     "h5::datatype::ieee_f32be", "metals_cold_gas" };
            field_dict["metalsstellarmass"] = { "h5::datatype::native_float", HOFFSET( galaxy, metals_stellar_mass ), "h5::datatype::ieee_f32be", "metals_stellar_mass"};
            field_dict["metalsbulgemass"] = { "h5::datatype::native_float", HOFFSET( galaxy, metals_bulge_mass ),   "h5::datatype::ieee_f32be", "metals_bulge_mass"};
            field_dict["metalshotgas"] = { "h5::datatype::native_float", HOFFSET( galaxy, metals_hot_gas ),      "h5::datatype::ieee_f32be", "metals_hot_gas"};
            field_dict["metalsejectedmass"] = { "h5::datatype::native_float", HOFFSET( galaxy, metals_ejected_mass ), "h5::datatype::ieee_f32be", "metals_ejected_mass"};
            field_dict["metalsics"] = { "h5::datatype::native_float", HOFFSET( galaxy, metals_ics ),          "h5::datatype::ieee_f32be", "metals_ics" };
            field_dict["sfrdisk"] = {  "h5::datatype::native_float", HOFFSET( galaxy, sfr_disk ),            "h5::datatype::ieee_f32be", "sfr_disk" };
            field_dict["sfrbulge"] = {  "h5::datatype::native_float", HOFFSET( galaxy, sfr_bulge ),           "h5::datatype::ieee_f32be", "sfr_bulge" };
            field_dict["sfrdisk"] = {  "h5::datatype::native_float", HOFFSET( galaxy, sfr_disk_z ),          "h5::datatype::ieee_f32be", "sfr_disk_z"};
            field_dict["sfrbulgez"] = {  "h5::datatype::native_float", HOFFSET( galaxy, sfr_bulge_z ),         "h5::datatype::ieee_f32be", "sfr_bulge_z" };
            field_dict["diskscaleradius"] = {  "h5::datatype::native_float", HOFFSET( galaxy, disk_scale_radius ),   "h5::datatype::ieee_f32be", "disk_scale_radius" };
            field_dict["cooling"] = {  "h5::datatype::native_float", HOFFSET( galaxy, cooling ),             "h5::datatype::ieee_f32be", "cooling" };
            field_dict["heating"] = {  "h5::datatype::native_float", HOFFSET( galaxy, heating ),             "h5::datatype::ieee_f32be", "heating" };
            field_dict["quasarmodebhaccretionmass"] = {  "h5::datatype::native_float", HOFFSET( galaxy, quasar_mode_bh_accretion_mass ), "h5::datatype::ieee_f32be", "quasar_mode_bh_accretion_mass"};
            field_dict["timeoflastmajormerger"] = {  "h5::datatype::native_float", HOFFSET( galaxy, time_of_last_major_merger ), "h5::datatype::ieee_f32be", "time_of_last_major_merger"};
            field_dict["timeoflastminormerger"] = {  "h5::datatype::native_float", HOFFSET( galaxy, time_of_last_minor_merger ), "h5::datatype::ieee_f32be", "time_of_last_minor_merger"};
            field_dict["outflowrate"] = {  "h5::datatype::native_float", HOFFSET( galaxy, outflow_rate ),        "h5::datatype::ieee_f32be", "outflow_rate"};
            field_dict["infallmvir"] = {  "h5::datatype::native_float", HOFFSET( galaxy, infall_mvir ),         "h5::datatype::ieee_f32be", "infall_mvir"};
            field_dict["infallvvir"] = {  "h5::datatype::native_float", HOFFSET( galaxy, infall_vvir ),         "h5::datatype::ieee_f32be", "infall_vvir"};
            field_dict["infallvmax"] = { "h5::datatype::native_float", HOFFSET( galaxy, infall_vmax ),         "h5::datatype::ieee_f32be", "infall_vmax"};

            // aliases and unknowns (dodgy)

            field_dict["snapnum"] = Entry("h5::datatype::native_int",   HOFFSET( galaxy, snapshot ),            "h5::datatype::std_i32be",  "snapshot" );
            field_dict["objecttype"] = {"h5::datatype::native_int",   HOFFSET( galaxy, type ),          "h5::datatype::std_i32be",  "type"};
            field_dict["totsfr"] = {"h5::datatype::native_float",   HOFFSET( galaxy, infall_vmax ),          "h5::datatype::ieee_f32be",  "infall_vmax"};
            field_dict["mvir"] = {"h5::datatype::native_float",   HOFFSET( galaxy, infall_mvir ),          "h5::datatype::ieee_f32be",  "infall_mvir"};
            field_dict["rvir"] = {"h5::datatype::native_float",   HOFFSET( galaxy, infall_vvir ),          "h5::datatype::ieee_f32be",  "infall_vvir"};
            field_dict["vvir"] = {"h5::datatype::native_float",   HOFFSET( galaxy, infall_vvir ),          "h5::datatype::ieee_f32be",  "infall_vvir"};
            field_dict["vmax"] = {"h5::datatype::native_float",   HOFFSET( galaxy, infall_vmax ),          "h5::datatype::ieee_f32be",  "infall_vmax"};
            field_dict["veldisp"] = {"h5::datatype::native_float",   HOFFSET( galaxy, infall_vmax ),          "h5::datatype::ieee_f32be",  "infall_vmax"};
            field_dict["veldisp"] = {  "h5::datatype::native_float", HOFFSET( galaxy, vel_disp ),            "h5::datatype::ieee_f32be", "velocity_dispersion" };
            field_dict["centralmvir"] = { "h5::datatype::native_float", HOFFSET( galaxy, central_mvir ),            "h5::datatype::ieee_f32be", "central_galaxy_mvir" };

            field_dict["posx"] = { "h5::datatype::native_float", HOFFSET( galaxy, pos[0] ),              "h5::datatype::ieee_f32be", "position_x"};
            field_dict["posy"] = { "h5::datatype::native_float", HOFFSET( galaxy, pos[1] ),              "h5::datatype::ieee_f32be", "position_y" };
            field_dict["posz"] = { "h5::datatype::native_float", HOFFSET( galaxy, pos[2] ),              "h5::datatype::ieee_f32be", "position_z" };
            field_dict["velx"] = { "h5::datatype::native_float", HOFFSET( galaxy, vel[0] ),              "h5::datatype::ieee_f32be", "velocity_x" };
            field_dict["vely"] = { "h5::datatype::native_float", HOFFSET( galaxy, vel[1] ),              "h5::datatype::ieee_f32be", "velocity_y" };
            field_dict["velz"] = { "h5::datatype::native_float", HOFFSET( galaxy, vel[2] ),              "h5::datatype::ieee_f32be", "velocity_z" };
            field_dict["mergeintosnapnum"] = {"h5::datatype::native_int",   HOFFSET( galaxy, merge_into_snapshot ), "h5::datatype::std_i32be",  "merge_into_snapshot" };
            field_dict["vpeak"] = {"h5::datatype::native_float",   HOFFSET( galaxy, infall_vmax ),          "h5::datatype::ieee_f32be",  "infall_vmax"};
            field_dict["isflyby"] = {"h5::datatype::native_int",   HOFFSET( galaxy, type ),          "h5::datatype::std_i32be",  "type"};
            field_dict["sfrdiskz"] = {  "h5::datatype::native_float", HOFFSET( galaxy, sfr_disk_z ),          "h5::datatype::ieee_f32be", "sfr_disk_z"};
            field_dict["len"] = {"h5::datatype::native_float",   HOFFSET( galaxy, type ),          "h5::datatype::std_i32be",  "type"};
            field_dict["treeindex"] = {"h5::datatype::native_float",   HOFFSET( galaxy, type ),          "h5::datatype::std_i32be",  "type"};

            // now go through in order
            // Declaring the type of Predicate that accepts 2 pairs and return a bool
            typedef std::function<bool(std::pair<std::string, int>, std::pair<std::string, int>)> Comparator;

            // Defining a lambda function to compare two pairs. It will compare two pairs using second field
            Comparator compFunctor =
                    [](std::pair<std::string, int> elem1, std::pair<std::string, int> elem2) {
                        return elem1.second < elem2.second;
                    };

            // Declaring a set that will store the pairs using above comparision logic
            std::set<std::pair<std::string, int>, Comparator> orderedFields(
                    _field_order.begin(), _field_order.end(), compFunctor);

            // write out a structure definition
            for (std::pair<std::string, int> element : orderedFields) {
                std::string name=element.first;
                switch (_field_types[name])
                {
                    case tao::batch<tao::real_type>::INTEGER:
                        std::cout <<"\tint "<<name<<";"<<std::endl;
                        break;
                    case tao::batch<tao::real_type>::LONG_LONG:
                        std::cout <<"\tlong long "<<name<<";"<<std::endl;
                        break;
                    case tao::batch<tao::real_type>::DOUBLE:
                        std::cout <<"\tfloat "<<name<<";"<<std::endl;
                        break;
                    default:
                        std::cout <<"\tunknown "<<name<<";"<<std::endl;
                        break;
                }

            }
            // write out a mem definition
            for (std::pair<std::string, int> element : orderedFields) {
                std::string name=element.first;
                switch (_field_types[name])
                {
                    case tao::batch<tao::real_type>::INTEGER:
                        std::cout <<"\tder.add(h5::datatype::native_int,HOFFSET(sage::galaxyTreeSed,"<<name<<"),h5::datatype::std_i32be,\""<<name<<"\");"<<std::endl;
                        break;
                    case tao::batch<tao::real_type>::LONG_LONG:
                        std::cout <<"\tder.add(h5::datatype::native_llong,HOFFSET(sage::galaxyTreeSed,"<<name<<"),h5::datatype::std_i64be,\""<<name<<"\");"<<std::endl;
                        break;
                    case tao::batch<tao::real_type>::DOUBLE:
                        std::cout <<"\tder.add(h5::datatype::native_float,HOFFSET(sage::galaxyTreeSed,"<<name<<"),h5::datatype::ieee_f32be,\""<<name<<"\");"<<std::endl;
                        break;
                    default:
                        std::cout <<"\tunknown "<<name<<";"<<std::endl;
                        break;
                }

            }

            // Iterate over a set of fields in order building the record definition as we go
            //hpc::h5::datatype mem_type;
            h5::derive der( sizeof(sage::galaxyTreeSed) );
            for (std::pair<std::string, int> element : orderedFields) {
                std::string name=element.first;
                boost::erase_all(name, "_");
                if (!field_dict[name].label.empty()) {
                    Entry entry = field_dict[name];
                    std::cout <<"\tder.add("<<entry.native_type<<","<<entry.offset<<","<<entry.mem_type<<",\""<<entry.label<<"\");"<<std::endl;
                } else {

                    //std::cout << "Unmapped "<<element.first << std::endl;
                }

            }

            sage::make_hdf5_types_treesed(_mem_type, _file_type);
            //=========================
        }

        //sage::make_hdf5_types(_mem_type, _file_type);
	std::cout << "OPEN MPI HDF5 INPUT"<<std::endl;
        hpc::h5::file                   file(_tree_fn.native(), H5F_ACC_RDONLY, hpc::mpi::comm::world);
        hpc::mpi::comm::world.barrier();
        std::vector<unsigned long long> gal_cnts   = count_galaxies(file);
        std::vector<unsigned long long> gal_displs = hpc::mpi::comm::world.scan(gal_cnts);
        LOGDLN("Counts: ", gal_cnts);
        LOGDLN("Displs: ", gal_displs);

	std::cout << "OPEN MPI HDF5"<<std::endl;
        hpc::h5::file                 out_file(_out_fn.native(), H5F_ACC_TRUNC, hpc::mpi::comm::world);
	    tao::data_dict::saveSidecar(_out_fn.native(), fields);
        unsigned                      n_snaps = gal_cnts.size();
	std::cout << "OPEN MPI HDF5 OUTPUT done."<<std::endl;
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
        std::vector<sage::galaxyTreeSed> buf(chunk_size);
        hpc::matrix<sage::galaxyTreeSed> out_buf(n_snaps, chunk_size);
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
                out_buf(buf[ii].snapnum, out_buf_sizes[buf[ii].snapnum]) = buf[ii];
                // ord_out_buf( buf[ii].snapshot, out_buf_sizes[buf[ii].snapshot] ) = my_displ + ii;
                ++out_buf_sizes[buf[ii].snapnum];
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
    void tree2sageANY() {
        LOGBLOCKI("Convert tree h5 file to sage hd5 file suitable for input to init");
        tao::xml_dict sidecar_xml = tao::data_dict::getSidecar(_tree_fn.native(), "/sageinput");
        std::vector<tao::data_dict_field> fields = tao::data_dict::getFieldsANY(sidecar_xml, "/sageinput/Field");
        size_t nfields=fields.size();
        int32_t snapnum_index = tao::data_dict::findField(fields,tao::SNAPNUM);

        tao::data_dict::getCompositeANY(fields, _mem_type, _file_type);
	std::cout << "OPEN MPI HDF5 INPUT"<<std::endl;
        hpc::h5::file                   file(_tree_fn.native(), H5F_ACC_RDONLY, hpc::mpi::comm::world);
        hpc::mpi::comm::world.barrier();
        std::vector<unsigned long long> gal_cnts   = count_galaxiesANY(file, nfields, snapnum_index);
        std::vector<unsigned long long> gal_displs = hpc::mpi::comm::world.scan(gal_cnts);
        LOGDLN("Counts: ", gal_cnts);
        LOGDLN("Displs: ", gal_displs);

	std::cout << "OPEN MPI HDF5"<<std::endl;
	tao::data_dict::saveSidecar(_out_fn.native(), fields, false);
        hpc::h5::file                 out_file(_out_fn.native(), H5F_ACC_TRUNC, hpc::mpi::comm::world);
        unsigned                      n_snaps = gal_cnts.size();
	    std::cout << "OPEN MPI HDF5 OUTPUT done."<<std::endl;
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
        std::vector<tao::ANY> ubuf(chunk_size*nfields);
        hpc::matrix<tao::ANY> uout_buf(n_snaps,chunk_size*nfields);
        std::vector<unsigned> out_buf_sizes(n_snaps);
        unsigned long long    my_n_gals = hpc::mpi::balanced_local_sizeg(n_gals);
        unsigned long long    my_displ  = hpc::mpi::comm::world.scan(my_n_gals);
        LOGDLN("Processing ", my_n_gals, " galaxies.");
        LOGDLN("Starting at ", my_displ, ".");

        int percentage=0;
        uint64_t todo_rank0=my_n_gals;
        while(my_n_gals) {
            
            int todo=std::min(my_n_gals, chunk_size);
            ubuf.resize(todo*nfields);
            LOGILN("Writing chunk from ", my_displ, " to ", my_displ + todo, " (", my_n_gals, " remaining)");
            gals_ds.read(ubuf.data(), _mem_type, todo, my_displ);
            boost::fill(out_buf_sizes, 0);
            for(unsigned ii = 0; ii < todo; ++ii) {
                int32_t snapnum = ubuf[nfields*ii+snapnum_index].i4;
                for (size_t iii=0; iii < nfields; ++iii)
                {
                    uout_buf(snapnum, out_buf_sizes[snapnum]*nfields+iii) = ubuf[ii*nfields+iii];
                }
                ++out_buf_sizes[snapnum];
            }
            my_displ += todo;
            my_n_gals -= todo;

            for(unsigned ii = 0; ii < n_snaps; ++ii) {
                if(out_buf_sizes[ii]) {
                    dsets[ii].write(uout_buf[ii].data(), _mem_type, out_buf_sizes[ii], gal_displs[ii]);
                    gal_displs[ii] += out_buf_sizes[ii];
                }
            }

            double sofar=todo_rank0-my_n_gals;
            sofar /= todo_rank0;
            sofar *= 100;
            if (std::rint(sofar)!=percentage) {
                percentage = std::rint(sofar);
                std::cout<< "Progress:"<<percentage<<std::endl;
            }
        }

        LOGDLN("Writing auxilliary data.");
        hpc::h5::copy(file, "cosmology", out_file);
        hpc::h5::copy(file, "snapshot_redshifts", out_file);

        hpc::mpi::comm::world.barrier();
    }
    void dump() {
        LOGBLOCKI("Convert tree h5 file to sage hd5 file suitable for input to init");

        // Work out the columns
        std::unordered_map<std::string,std::string> _field_map;
        std::map<std::string,typename tao::batch<tao::real_type>::field_value_type> _field_types;
        {

            tao::xml_dict xml;
            xml.read("settings.xml", "/settings");
            // Iterate over the module nodes.
            pugi::xpath_node_set fields = xml.get_nodes("/settings/sageinput/Field");
            std::map<std::string, int> _field_order;

            _field_types.clear();

            int order=0;
            for (const pugi::xpath_node *it = fields.begin(); it != fields.end(); ++it) {
                pugi::xml_node cur = it->node();

                std::string element = cur.name();
                std::string field_str = cur.text().as_string();
                std::string label = cur.attribute("label").value();
                std::string description = cur.attribute("description").value();
                std::string units = cur.attribute("units").value();
                std::string group = cur.attribute("group").value();
                std::string type_str = cur.attribute("Type").value();
                std::cout << "Field[" << field_str << "]=" << label << "," << description << std::endl;
                typename tao::batch<tao::real_type>::field_value_type type;
                if (type_str == "short")
                    type = tao::batch<tao::real_type>::INTEGER;
                else if (type_str == "int")
                    type = tao::batch<tao::real_type>::INTEGER;
                else if (type_str == "long long")
                    type = tao::batch<tao::real_type>::LONG_LONG;
                else if (type_str == "float")
                    type = tao::batch<tao::real_type>::DOUBLE;
                else {
                    EXCEPT(0, "Unknown field type for field '", field_str, "': ", type_str);
                }
                _field_types.emplace(field_str, type);
                _field_map[field_str] = field_str;
                _field_order[field_str] = order;
                order++;

            }

            // now go through in order
            // Declaring the type of Predicate that accepts 2 pairs and return a bool
            typedef std::function<bool(std::pair<std::string, int>, std::pair<std::string, int>)> Comparator;

            // Defining a lambda function to compare two pairs. It will compare two pairs using second field
            Comparator compFunctor =
                    [](std::pair<std::string, int> elem1, std::pair<std::string, int> elem2) {
                        return elem1.second < elem2.second;
                    };

            // Declaring a set that will store the pairs using above comparision logic
            std::set<std::pair<std::string, int>, Comparator> orderedFields(
                    _field_order.begin(), _field_order.end(), compFunctor);

            // write out a structure definition
            for (std::pair<std::string, int> element : orderedFields) {
                std::string name=element.first;
                std::string name_lowercase=name;
                std::transform(name_lowercase.begin(), name_lowercase.end(), name_lowercase.begin(), ::tolower);

                switch (_field_types[name])
                {
                    case tao::batch<tao::real_type>::INTEGER:
                        std::cout <<"\tint "<<name_lowercase<<";"<<std::endl;
                        break;
                    case tao::batch<tao::real_type>::LONG_LONG:
                        std::cout <<"\tlong long "<<name_lowercase<<";"<<std::endl;
                        break;
                    case tao::batch<tao::real_type>::DOUBLE:
                        std::cout <<"\tfloat "<<name_lowercase<<";"<<std::endl;
                        break;
                    default:
                        break;
                }

            }
            // write out a mem definition
            for (std::pair<std::string, int> element : orderedFields) {
                std::string name=element.first;
                std::string name_lowercase=name;
                std::transform(name_lowercase.begin(), name_lowercase.end(), name_lowercase.begin(), ::tolower);
                switch (_field_types[name])
                {
                    case tao::batch<tao::real_type>::INTEGER:
                        std::cout <<"\tder.add(h5::datatype::native_int,HOFFSET(sage::galaxyTreeSed,"<<name_lowercase<<"),h5::datatype::std_i32be,\""<<name<<"\");"<<std::endl;
                        break;
                    case tao::batch<tao::real_type>::LONG_LONG:
                        std::cout <<"\tder.add(h5::datatype::native_llong,HOFFSET(sage::galaxyTreeSed,"<<name_lowercase<<"),h5::datatype::std_i64be,\""<<name<<"\");"<<std::endl;
                        break;
                    case tao::batch<tao::real_type>::DOUBLE:
                        std::cout <<"\tder.add(h5::datatype::native_float,HOFFSET(sage::galaxyTreeSed,"<<name_lowercase<<"),h5::datatype::ieee_f32be,\""<<name<<"\");"<<std::endl;
                        break;
                    default:
                        break;
                }

            }
            //=========================
        }

        sage::make_hdf5_types_treesed(_mem_type, _file_type);
        hpc::h5::file                   file(_tree_fn.native(), H5F_ACC_RDONLY, hpc::mpi::comm::world);

        {
            std::vector<unsigned long long> cnts(file.dataset("snapshot_redshifts").extent());
            boost::fill(cnts, 0);
            hpc::h5::dataset          gals_ds(file, "galaxies");
            unsigned long long        n_gals = gals_ds.extent();
            std::vector<sage::galaxyTreeSed> buf(chunk_size);
            unsigned long long        my_n_gals = hpc::mpi::balanced_local_sizeg(n_gals);
            unsigned long long        my_displ  = hpc::mpi::comm::world.scan(my_n_gals);
            while(my_n_gals) {
                buf.resize(std::min(my_n_gals, chunk_size));
                LOGILN("Counting chunk from ", my_displ, " to ", my_displ + buf.size(), " (", my_n_gals, " remaining)");
                gals_ds.read(buf.data(), _mem_type, buf.size(), my_displ);
                for(unsigned ii = 0; ii < buf.size(); ++ii)
                {
                    if (buf[ii].globalindex==4460)
                        std::cout << "galaxy["<<(my_displ+ii)<<"]="<<buf[ii].stellarmass<<","<<buf[ii].globalindex<<","<<buf[ii].diskscaleradius<<std::endl;
                }
                my_displ += buf.size();
                my_n_gals -= buf.size();
            }
        }
        std::vector<unsigned long long> gal_cnts   = count_galaxies(file);
        std::vector<unsigned long long> gal_displs = hpc::mpi::comm::world.scan(gal_cnts);
        LOGDLN("Counts: ", gal_cnts);
        LOGDLN("Displs: ", gal_displs);

        unsigned                      n_snaps = gal_cnts.size();
        std::vector<hpc::h5::dataset> dsets(n_snaps);
        // std::vector<hpc::h5::dataset> ord_dsets( n_snaps );
        for(unsigned ii = 0; ii < n_snaps; ++ii) {
            unsigned long long gsize = hpc::mpi::comm::world.all_reduce(gal_cnts[ii]);
            //dsets[ii].create(out_file, std::string("snapshot") + hpc::index_string(ii, 3), _file_type, gsize);
            // ord_dsets[ii].create( out_file, std::string( "snapshot" ) + hpc::index_string( ii, 3 ) + "_map",
            //                       hpc::h5::datatype::native_ullong, gsize );
        }

        LOGILN("Flattening ", n_snaps, " snapshots.");
        hpc::h5::dataset          gals_ds(file, "galaxies");
        unsigned long long        n_gals = gals_ds.extent();
        std::vector<sage::galaxyTreeSed> buf(chunk_size);
        hpc::matrix<sage::galaxyTreeSed> out_buf(n_snaps, chunk_size);
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
                out_buf(buf[ii].snapnum, out_buf_sizes[buf[ii].snapnum]) = buf[ii];
                // ord_out_buf( buf[ii].snapshot, out_buf_sizes[buf[ii].snapshot] ) = my_displ + ii;
                ++out_buf_sizes[buf[ii].snapnum];
            }
            my_displ += buf.size();
            my_n_gals -= buf.size();
        }

        LOGDLN("Writing auxilliary data.");
        //hpc::h5::copy(file, "cosmology", out_file);
        //hpc::h5::copy(file, "snapshot_redshifts", out_file);

        hpc::mpi::comm::world.barrier();
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
    hpc::fs::path     _txt_fn;
    hpc::h5::datatype _mem_type, _file_type;
    hpc::h5::datatype _lc_mem_type, _lc_file_type;
    std::vector<tao::data_dict_field> _fields;
    unsigned          _ppc;
};

#define HPC_APP_CLASS myapplication
#include <libhpc/mpi/main.hh>


