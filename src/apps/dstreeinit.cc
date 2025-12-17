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
#include <set>
#include <algorithm>
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
        } else if (_mode == "kdtree") {
            EXCEPT(!_sage_fn.empty(), "No SAGE filename given.");
            EXCEPT(!_tree_fn.empty(), "No SAGE trees file given.");
            EXCEPT(!_out_fn.empty(), "No output filename given.");
        }
    }

    void operator()() {
        if(_mode == "init")
            initANY(_sage_fn, _tree_fn, _out_fn);
        else if(_mode == "tree2sage")
            tree2sageANY(_tree_fn, _out_fn);
        else if(_mode == "kdtree") {
            tree2sageANY(_tree_fn, _sage_fn);
            initANY(_sage_fn, _tree_fn, _out_fn);
        }
    }

    void initANY(hpc::fs::path const &sage_fn, hpc::fs::path const &tree_fn, hpc::fs::path const &out_fn) {
            LOGBLOCKI("Initialising new file.");

            tao::xml_dict sidecar_xml = tao::data_dict::getSidecar(sage_fn.string(),"/sageinput");
            _fields = tao::data_dict::getFieldsANY(sidecar_xml,"/sageinput/Field");
            size_t nfields=_fields.size();
            hpc::h5::derive der = tao::data_dict::getCompositeANY(_fields, _mem_type, _file_type);

            hpc::h5::file         file(sage_fn.native(), H5F_ACC_RDONLY, hpc::mpi::comm::world);
            
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

            // Prepare the output file and its sideCar.
            hpc::h5::file                   out_file(out_fn.native(), H5F_ACC_TRUNC, hpc::mpi::comm::world);
            tao::data_dict::saveSidecar(out_fn.native(), _fields, false);

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
            write_sed_dataANY(out_file, tree_fn);

            // Add in cosmology stuff.
            hpc::h5::copy(file, "cosmology", out_file);
            hpc::h5::copy(file, "snapshot_redshifts", out_file);
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


    void write_sed_dataANY(hpc::h5::file &out_file, hpc::fs::path const &tree_fn) {
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

            hpc::h5::file             in_file(tree_fn.native(), H5F_ACC_RDONLY);
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

    void tree2sageANY(hpc::fs::path const &tree_fn, hpc::fs::path const &out_fn) {
        LOGBLOCKI("Convert tree h5 file to sage hd5 file suitable for input to init");
        tao::xml_dict sidecar_xml = tao::data_dict::getSidecar(tree_fn.native(), "/sageinput");
        std::vector<tao::data_dict_field> fields = tao::data_dict::getFieldsANY(sidecar_xml, "/sageinput/Field");
        size_t nfields=fields.size();
        int32_t snapnum_index = tao::data_dict::findField(fields,tao::SNAPNUM);

        tao::data_dict::getCompositeANY(fields, _mem_type, _file_type);
	std::cout << "OPEN MPI HDF5 INPUT"<<std::endl;
        hpc::h5::file                   file(tree_fn.native(), H5F_ACC_RDONLY, hpc::mpi::comm::world);
        hpc::mpi::comm::world.barrier();
        std::vector<unsigned long long> gal_cnts   = count_galaxiesANY(file, nfields, snapnum_index);
        std::vector<unsigned long long> gal_displs = hpc::mpi::comm::world.scan(gal_cnts);
        LOGDLN("Counts: ", gal_cnts);
        LOGDLN("Displs: ", gal_displs);

	std::cout << "OPEN MPI HDF5"<<std::endl;
	tao::data_dict::saveSidecar(out_fn.native(), fields, false);
        hpc::h5::file                 out_file(out_fn.native(), H5F_ACC_TRUNC, hpc::mpi::comm::world);
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


  protected:
    std::string       _mode;
    hpc::fs::path     _sage_fn;
    hpc::fs::path     _out_fn;
    hpc::fs::path     _tree_fn;
    hpc::h5::datatype _mem_type, _file_type;
    hpc::h5::datatype _lc_mem_type, _lc_file_type;
    std::vector<tao::data_dict_field> _fields;
    unsigned          _ppc;
};

#define HPC_APP_CLASS myapplication
#include <libhpc/mpi/main.hh>


