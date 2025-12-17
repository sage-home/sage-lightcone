/**
 *  @file    sageh5tokdtree.cc
 *  @author  GitHub Copilot
 *  @date    6 December 2025
 *
 *  @brief   Converts SAGE HDF5 output to TAO KDTree HDF5 format.
 */

#include <iostream>
#include <vector>
#include <string>
#include <numeric>
#include <cstring>
#include <unordered_map>
#include <map>
#include <algorithm>
#include <boost/range/algorithm/fill.hpp>
#include <libhpc/libhpc.hh>
#include <libtao/base/utils.hh>
#include <libtao/base/sage.hh>
#include "libhpc/algorithm/xdmf_writer.hh"

using namespace ::tao;
using namespace ::hpc;

// Structure to handle data permutation for KDTree
struct data_permuter {
    data_permuter() : crds(nullptr), tree_idxs(nullptr), subsize(nullptr), idxs(nullptr), galaxy_idxs(nullptr) {
    }

    data_permuter(std::array<std::vector<double>, 3> &crds,
                  std::vector<unsigned long long> &   tree_idxs,
                  std::vector<unsigned> &             subsize,
                  std::vector<unsigned long long> &   idxs,
                  std::vector<long long> &            galaxy_idxs) :
        crds(&crds),
        tree_idxs(&tree_idxs), subsize(&subsize), idxs(&idxs), galaxy_idxs(&galaxy_idxs) {
    }

    void operator()(hpc::mpi::balanced_partition const &part) {
        for(unsigned ii = 0; ii < 3; ++ii)
            part.transfer((*crds)[ii]);
        part.transfer(*tree_idxs);
        part.transfer(*subsize);
        part.transfer(*idxs);
        part.transfer(*galaxy_idxs);
    }

    std::array<std::vector<double>, 3> *crds;
    std::vector<unsigned long long> *   tree_idxs;
    std::vector<unsigned> *             subsize;
    std::vector<unsigned long long> *   idxs;
    std::vector<long long> *            galaxy_idxs;
};

struct GalaxyInfo {
    int32_t tree_id;
    int32_t snap_num;
    int32_t file_index; // Original index in the snapshot file
    int64_t galaxy_index;
    int32_t merge_id;
    int32_t merge_snap;
    int32_t descendant_file_index = -1; 
    int32_t descendant_snap = -1;
    uint32_t subsize = 0;
};

class sageh5tokdtree_app : public hpc::mpi::application {
  public:
    sageh5tokdtree_app(int argc, char *argv[]) : hpc::mpi::application(argc, argv) {
        LOG_PUSH(new hpc::mpi::logger("sageh5tokdtree.log.", hpc::log::info));

        options().add_options()
            ("sage,s", hpc::po::value<hpc::fs::path>(&_sage_fn)->required(), "Input SAGE HDF5 file")
            ("output,o", hpc::po::value<hpc::fs::path>(&_out_fn)->required(), "Output KDTree HDF5 file")
            ("ppc,p", hpc::po::value<unsigned>(&_ppc)->default_value(1000), "Particles per cell");

        parse_options(argc, argv);
    }

    void operator()() {
        LOGBLOCKI("Starting SAGE HDF5 to KDTree conversion.");
        LOGILN("Input file: ", _sage_fn);
        LOGILN("Output file: ", _out_fn);

        hpc::h5::file in_file(_sage_fn.string(), H5F_ACC_RDONLY, hpc::mpi::comm::world);
        hpc::h5::file out_file(_out_fn.string(), H5F_ACC_TRUNC, hpc::mpi::comm::world);

        // Copy cosmology info
        copy_cosmology(in_file, out_file);

        // Create groups
        out_file.create_group("data");
        out_file.create_group("lightcone");

        // Determine number of snapshots
        int n_snaps = 0;
        while(in_file.has_link("Snap_" + std::to_string(n_snaps))) {
            n_snaps++;
        }
        LOGILN("Found ", n_snaps, " snapshots.");

        // Pass 1: Build Tree Structure and Calculate Subsize
        build_tree_structure(in_file, n_snaps);

        // Pass 2: Process each snapshot and write data
        _total_count = 0;
        
        for(int i = 0; i < n_snaps; ++i) {
            process_snapshot(in_file, i, out_file);
        }

        // Add final displacement
        _snapshot_displs.push_back(_total_count);

        // Write global metadata
        LOGILN("Writing global metadata...");
        out_file.write("snapshot_counts", _snapshot_counts);
        out_file.write("snapshot_displs", _snapshot_displs);
        out_file.write("snapshot_redshifts", _snapshot_redshifts);
        
        LOGILN("Conversion complete.");
    }

  private:
    void build_tree_structure(hpc::h5::file &in, int n_snaps) {
        LOGBLOCKI("Pass 1: Building Tree Structure...");
        
        std::vector<GalaxyInfo> all_gals;
        
        // Read structure from all snapshots
        for(int i = 0; i < n_snaps; ++i) {
            std::string snap_name = "Snap_" + std::to_string(i);
            auto snap_group = in.group(snap_name);
            
            std::vector<int> sage_tree_idxs;
            std::vector<long long> galaxy_idxs;
            std::vector<int> merge_snap;
            std::vector<int> merge_id;
            
            read_dataset(snap_group, "SAGETreeIndex", sage_tree_idxs);
            
            if(snap_group.has_link("GalaxyIndex")) {
                read_dataset(snap_group, "GalaxyIndex", galaxy_idxs);
            } else {
                galaxy_idxs.resize(sage_tree_idxs.size(), -1);
            }
            
            bool has_mergers = snap_group.has_link("mergeIntoSnapNum") && snap_group.has_link("mergeIntoID");
            if(has_mergers) {
                read_dataset(snap_group, "mergeIntoSnapNum", merge_snap);
                read_dataset(snap_group, "mergeIntoID", merge_id);
            } else {
                merge_snap.resize(sage_tree_idxs.size(), -1);
                merge_id.resize(sage_tree_idxs.size(), -1);
            }
            
            size_t n_gals = sage_tree_idxs.size();
            for(size_t k = 0; k < n_gals; ++k) {
                GalaxyInfo g;
                g.tree_id = sage_tree_idxs[k];
                g.snap_num = i;
                g.file_index = (int32_t)k;
                g.galaxy_index = galaxy_idxs[k];
                g.merge_id = merge_id[k];
                g.merge_snap = merge_snap[k];
                all_gals.push_back(g);
            }
        }
        
        LOGILN("Total galaxies loaded: ", all_gals.size());
        
        // Sort by TreeID, then SnapNum
        std::sort(all_gals.begin(), all_gals.end(), [](const GalaxyInfo &a, const GalaxyInfo &b) {
            if(a.tree_id != b.tree_id) return a.tree_id < b.tree_id;
            return a.snap_num < b.snap_num;
        });
        
        // Process trees
        _subsize_lookup.resize(n_snaps);
        
        size_t start_idx = 0;
        while(start_idx < all_gals.size()) {
            size_t end_idx = start_idx + 1;
            while(end_idx < all_gals.size() && all_gals[end_idx].tree_id == all_gals[start_idx].tree_id) {
                end_idx++;
            }
            
            process_single_tree(all_gals, start_idx, end_idx);
            start_idx = end_idx;
        }
        
        LOGILN("Tree structure built.");
    }
    
    void process_single_tree(std::vector<GalaxyInfo> &gals, size_t start, size_t end) {
        // Map GalaxyIndex -> Index in 'gals' (for main branch)
        // Map {Snap, FileIndex} -> Index in 'gals' (for mergers)
        
        std::unordered_map<long long, size_t> galaxy_map;
        std::map<std::pair<int, int>, size_t> file_map;
        
        for(size_t i = start; i < end; ++i) {
            if(gals[i].galaxy_index != -1) {
                galaxy_map[gals[i].galaxy_index] = i;
            }
            file_map[{gals[i].snap_num, gals[i].file_index}] = i;
        }
        
        // Initialize subsize to 1
        for(size_t i = start; i < end; ++i) {
            gals[i].subsize = 1;
        }
        
        // Pre-link main branches
        // Map GalaxyIndex -> Vector of indices (sorted by snap)
        std::unordered_map<long long, std::vector<size_t>> main_branches;
        for(size_t i = start; i < end; ++i) {
            if(gals[i].galaxy_index != -1) {
                main_branches[gals[i].galaxy_index].push_back(i);
            }
        }
        
        // Propagate subsize forward (Progenitor -> Descendant)
        for(size_t i = start; i < end; ++i) {
            // If merger, propagate to target
            if(gals[i].merge_id != -1) {
                auto key = std::make_pair(gals[i].merge_snap, gals[i].merge_id);
                if(file_map.count(key)) {
                    size_t target = file_map[key];
                    gals[target].subsize += gals[i].subsize;
                }
            } else {
                // Main branch: propagate to next in chain
                if(gals[i].galaxy_index != -1) {
                    auto &chain = main_branches[gals[i].galaxy_index];
                    // Find current 'i' in chain
                    auto it = std::lower_bound(chain.begin(), chain.end(), i);
                    if(it != chain.end() && *it == i) {
                        auto next = it + 1;
                        if(next != chain.end()) {
                            gals[*next].subsize += gals[i].subsize;
                        }
                    }
                }
            }
        }
        
        // Store in lookup
        for(size_t i = start; i < end; ++i) {
            if(_subsize_lookup[gals[i].snap_num].size() <= gals[i].file_index) {
                _subsize_lookup[gals[i].snap_num].resize(gals[i].file_index + 1000); // Resize with buffer
            }
            // Ensure size
            if(_subsize_lookup[gals[i].snap_num].size() <= gals[i].file_index) {
                 _subsize_lookup[gals[i].snap_num].resize(gals[i].file_index + 1);
            }
            _subsize_lookup[gals[i].snap_num][gals[i].file_index] = gals[i].subsize;
        }
    }

    void write_kdtree(hpc::h5::group &      lc_snap,
                      hpc::kdtree<> const &kdt,
                      hpc::binary_partitioner<std::array<std::vector<double>, 3>::iterator, data_permuter> const &bp,
                      unsigned long long displ) {
        std::string name = "."; // Relative to lc_snap

        // Write bounds.
        {
            hpc::h5::derive der(2 * sizeof(double));
            der.add(hpc::h5::datatype::native_double, 0, hpc::h5::datatype::ieee_f64be, "minimum");
            der.add(hpc::h5::datatype::native_double, sizeof(double), hpc::h5::datatype::ieee_f64be, "maximum");
            hpc::h5::datatype mdt, fdt;
            der.commit(mdt, fdt);

            hpc::h5::dataset dset(lc_snap, "bounds", fdt, 3);
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

            hpc::h5::dataset dset(lc_snap, "splits", fdt, kdt.splits().size());
            dset.write(kdt.splits().data(), mdt, kdt.splits().size(), 0);
        }

        // Write counts and offsets.
        {
            auto counts = bp.counts();
            hpc::h5::dataset dset(lc_snap, "cell_counts", hpc::h5::datatype::native_uint, counts.size());
            dset.write(counts);
        }
        {
            std::vector<unsigned long long> offs(bp.offsets().size());
            for(unsigned ii = 0; ii < offs.size(); ++ii)
                offs[ii] = displ + bp.offsets()[ii];
            
            hpc::h5::dataset dset(lc_snap, "cell_offs", hpc::h5::datatype::native_ullong, offs.size());
            dset.write(offs);
        }
    }

    void append_dataset_compound(hpc::h5::group &data_group, const std::string &name, 
                        hpc::h5::datatype &mem_type, hpc::h5::datatype &file_type, 
                        const void *data, size_t count, unsigned long long offset) {
        
        hpc::h5::dataset dset;
        if(data_group.has_link(name)) {
            dset = data_group.dataset(name);
            dset.set_extent(offset + count);
        } else {
            hsize_t dims[1] = {offset + count};
            hsize_t max_dims[1] = {H5S_UNLIMITED};
            hsize_t chunk_dims[1] = {1024 * 1024 / mem_type.size()}; 
            if(chunk_dims[0] == 0) chunk_dims[0] = 1;
            
            hid_t space_id = H5Screate_simple(1, dims, max_dims);
            hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
            H5Pset_chunk(plist_id, 1, chunk_dims);
            
            hid_t dset_id = H5Dcreate2(data_group.id(), name.c_str(), file_type.id(), space_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);
            dset = hpc::h5::dataset(dset_id, false); 
            
            H5Pclose(plist_id);
            H5Sclose(space_id);
        }
        
        dset.write(data, mem_type, count, offset);
    }

    void append_lightcone_data(hpc::h5::group &lc_group, 
                               const std::array<std::vector<double>, 3> &crds,
                               const std::vector<long long> &global_idxs,
                               const std::vector<long long> &galaxy_idxs,
                               const std::vector<unsigned> &subsize,
                               unsigned long long offset) {
        size_t count = crds[0].size();
        if(count == 0) return;

        struct particle_t {
            double x, y, z;
            long long global_index;
            long long galaxy_index;
            unsigned subsize;
        };

        std::vector<particle_t> data(count);
        for(size_t i=0; i<count; ++i) {
            data[i].x = crds[0][i];
            data[i].y = crds[1][i];
            data[i].z = crds[2][i];
            data[i].global_index = global_idxs[i];
            data[i].galaxy_index = galaxy_idxs[i];
            data[i].subsize = subsize[i];
        }

        hpc::h5::derive der(sizeof(particle_t));
        der.add(hpc::h5::datatype::native_double, HOFFSET(particle_t, x), hpc::h5::datatype::ieee_f64be, "x");
        der.add(hpc::h5::datatype::native_double, HOFFSET(particle_t, y), hpc::h5::datatype::ieee_f64be, "y");
        der.add(hpc::h5::datatype::native_double, HOFFSET(particle_t, z), hpc::h5::datatype::ieee_f64be, "z");
        der.add(hpc::h5::datatype::native_llong, HOFFSET(particle_t, global_index), hpc::h5::datatype::std_u64be, "global_index");
        der.add(hpc::h5::datatype::native_llong, HOFFSET(particle_t, galaxy_index), hpc::h5::datatype::std_u64be, "central_galaxy_index");
        der.add(hpc::h5::datatype::native_uint, HOFFSET(particle_t, subsize), hpc::h5::datatype::std_u32be, "subsize");
        
        hpc::h5::datatype mdt, fdt;
        der.commit(mdt, fdt);

        append_dataset_compound(lc_group, "data", mdt, fdt, data.data(), count, offset);
    }

    void append_dataset(hpc::h5::group &data_group, const std::string &name, 
                        hid_t type_id, size_t element_size, 
                        const void *data, size_t count, unsigned long long offset) {
        
        hpc::h5::dataset dset;
        if(data_group.has_link(name)) {
            dset = data_group.dataset(name);
            dset.set_extent(offset + count);
        } else {
            // Create chunked dataset
            hsize_t dims[1] = {offset + count};
            hsize_t max_dims[1] = {H5S_UNLIMITED};
            hsize_t chunk_dims[1] = {1024 * 1024 / element_size}; // ~1MB chunks
            if(chunk_dims[0] == 0) chunk_dims[0] = 1;
            
            hid_t space_id = H5Screate_simple(1, dims, max_dims);
            hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
            H5Pset_chunk(plist_id, 1, chunk_dims);
            
            hid_t dset_id = H5Dcreate2(data_group.id(), name.c_str(), type_id, space_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);
            dset = hpc::h5::dataset(dset_id, false); 
            
            H5Pclose(plist_id);
            H5Sclose(space_id);
        }
        
        // Write to hyperslab using libhpc wrapper
        // write(buf, type, size, offset)
        dset.write(data, hpc::h5::datatype(H5Tcopy(type_id)), count, offset);
    }

    void permute_and_write_generic(hpc::h5::group &in_group, const std::string &name, 
                                   hpc::h5::group &out_data_group, 
                                   const std::vector<unsigned long long> &idxs, 
                                   unsigned long long offset) {
        
        hid_t dset_id = H5Dopen2(in_group.id(), name.c_str(), H5P_DEFAULT);
        if(dset_id < 0) return;
        
        hid_t type_id = H5Dget_type(dset_id);
        hid_t space_id = H5Dget_space(dset_id);
        size_t size = H5Tget_size(type_id);
        hssize_t npoints = H5Sget_simple_extent_npoints(space_id);
        
        std::vector<char> buf(size * npoints);
        H5Dread(dset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
        
        // Permute
        std::vector<char> permuted(size * npoints);
        for(int i=0; i<npoints; ++i) {
            size_t old_idx = idxs[i];
            std::memcpy(&permuted[i*size], &buf[old_idx*size], size);
        }
        
        append_dataset(out_data_group, name, type_id, size, permuted.data(), npoints, offset);
        
        H5Sclose(space_id);
        H5Tclose(type_id);
        H5Dclose(dset_id);
    }

    void copy_cosmology(hpc::h5::file &in, hpc::h5::file &out) {
        hpc::h5::group out_cosmo(out, "cosmology");
        
        // Try to read from Header/Simulation
        if(in.has_link("Header/Simulation")) {
            auto sim = in.group("Header/Simulation");
            copy_attr(sim, out_cosmo, "Hubble_h");
            copy_attr(sim, out_cosmo, "Omega_b");
            copy_attr(sim, out_cosmo, "Omega_lambda");
            copy_attr(sim, out_cosmo, "Omega_m");
            copy_attr(sim, out_cosmo, "PartMass");
            copy_attr(sim, out_cosmo, "BoxSize");
        } else {
            LOGILN("Header/Simulation group not found. Cosmology info might be missing.");
        }
    }

    void copy_attr(hpc::h5::group &in, hpc::h5::group &out, const std::string &name, const std::string &out_name = "") {
        hid_t in_id = in.id();
        if(H5Aexists(in_id, name.c_str()) > 0) {
            hid_t attr_id = H5Aopen(in_id, name.c_str(), H5P_DEFAULT);
            if(attr_id >= 0) {
                hid_t type_id = H5Aget_type(attr_id);
                hid_t space_id = H5Aget_space(attr_id);
                
                size_t size = H5Tget_size(type_id);
                std::vector<char> buf(size); 
                hssize_t npoints = H5Sget_simple_extent_npoints(space_id);
                if(npoints > 0) buf.resize(size * npoints);
                
                H5Aread(attr_id, type_id, buf.data());
                
                hid_t out_id = out.id();
                std::string oname = out_name.empty() ? name : out_name;
                hid_t out_attr_id = H5Acreate2(out_id, oname.c_str(), type_id, space_id, H5P_DEFAULT, H5P_DEFAULT);
                if(out_attr_id >= 0) {
                    H5Awrite(out_attr_id, type_id, buf.data());
                    H5Aclose(out_attr_id);
                }
                
                H5Sclose(space_id);
                H5Tclose(type_id);
                H5Aclose(attr_id);
            }
        }
    }

    void process_snapshot(hpc::h5::file &in, int snap_idx, hpc::h5::file &out) {
        std::string snap_name = "Snap_" + std::to_string(snap_idx);
        LOGBLOCKI("Processing ", snap_name);
        
        auto snap_group = in.group(snap_name);
        
        // Read coordinates for KDTree
        std::array<std::vector<double>, 3> crds;
        std::vector<unsigned long long>    tree_idxs;
        std::vector<unsigned>              subsize;
        
        // Read datasets
        read_dataset(snap_group, "Posx", crds[0]);
        read_dataset(snap_group, "Posy", crds[1]);
        read_dataset(snap_group, "Posz", crds[2]);
        
        // Read SAGETreeIndex for grouping
        std::vector<int> sage_tree_idxs;
        if(snap_group.has_link("SAGETreeIndex")) {
            read_dataset(snap_group, "SAGETreeIndex", sage_tree_idxs);
        } else {
            sage_tree_idxs.resize(crds[0].size(), 0);
        }

        // Read GalaxyIndex
        std::vector<long long> galaxy_idxs;
        if(snap_group.has_link("GalaxyIndex")) {
            read_dataset(snap_group, "GalaxyIndex", galaxy_idxs);
        } else {
            galaxy_idxs.resize(crds[0].size(), -1);
        }
        
        size_t n_gals = crds[0].size();
        
        // Accumulate metadata
        _snapshot_counts.push_back(n_gals);
        _snapshot_displs.push_back(_total_count);
        
        double z = 0.0;
        if(H5Aexists(snap_group.id(), "redshift") > 0) {
            hid_t attr_id = H5Aopen(snap_group.id(), "redshift", H5P_DEFAULT);
            if(attr_id >= 0) {
                H5Aread(attr_id, H5T_NATIVE_DOUBLE, &z);
                H5Aclose(attr_id);
            }
        } else if(H5Aexists(snap_group.id(), "Redshift") > 0) {
            hid_t attr_id = H5Aopen(snap_group.id(), "Redshift", H5P_DEFAULT);
            if(attr_id >= 0) {
                H5Aread(attr_id, H5T_NATIVE_DOUBLE, &z);
                H5Aclose(attr_id);
            }
        }
        _snapshot_redshifts.push_back(z);

        // Create output group for this snapshot
        char buf[16];
        snprintf(buf, sizeof(buf), "%03d", snap_idx);
        std::string lc_snap_name = "lightcone/snapshot" + std::string(buf);
        hpc::h5::group lc_snap(out, lc_snap_name);

        // Check if we have data
        if(n_gals == 0) {
            LOGILN("Snapshot ", snap_idx, " is empty.");
            return;
        }

        // Sort indices by SAGETreeIndex to group trees
        std::vector<unsigned long long> idxs(n_gals);
        std::iota(idxs.begin(), idxs.end(), 0);
        std::stable_sort(idxs.begin(), idxs.end(), [&](unsigned long long a, unsigned long long b) {
            return sage_tree_idxs[a] < sage_tree_idxs[b];
        });
        
        // Retrieve pre-calculated subsize
        subsize.resize(n_gals);
        if(snap_idx < _subsize_lookup.size()) {
            const auto &lookup = _subsize_lookup[snap_idx];
            for(size_t i = 0; i < n_gals; ++i) {
                // idxs[i] is the original index (file index)
                size_t k = idxs[i];
                if(k < lookup.size()) {
                    subsize[k] = lookup[k];
                } else {
                    subsize[k] = 1; // Fallback
                }
            }
        } else {
            std::fill(subsize.begin(), subsize.end(), 1);
        }

        // Permute data arrays to match sorted order
        std::vector<double> px(n_gals), py(n_gals), pz(n_gals);
        std::vector<unsigned long long> p_tree_idxs(n_gals);
        std::vector<unsigned> p_subsize(n_gals);
        std::vector<long long> p_galaxy_idxs(n_gals);
        
        for(size_t i=0; i<n_gals; ++i) {
            size_t k = idxs[i];
            px[i] = crds[0][k];
            py[i] = crds[1][k];
            pz[i] = crds[2][k];
            // Generate sequential global index
            p_tree_idxs[i] = _total_count + i; 
            p_subsize[i] = subsize[k];
            p_galaxy_idxs[i] = galaxy_idxs[k];
        }
        
        // Update crds to permuted values
        crds[0] = px;
        crds[1] = py;
        crds[2] = pz;
        tree_idxs = p_tree_idxs;
        subsize = p_subsize;
        galaxy_idxs = p_galaxy_idxs;
        
        // Reset idxs for KDTree partitioner (it expects 0..N-1)
        std::iota(idxs.begin(), idxs.end(), 0);

        data_permuter perm(crds, tree_idxs, subsize, idxs, galaxy_idxs);
        auto          part = hpc::make_binary_partitioner(crds.begin(), crds.end(), perm, _ppc);
        hpc::kdtree<> tree;
        {
            tree.construct(part);
        }
        
        // Write KDTree info to lightcone/snapshotNNN
        write_kdtree(lc_snap, tree, part, _snapshot_displs.back());
        
        // Permute global_index (which is tree_idxs)
        std::vector<long long> permuted_global_idxs(n_gals);
        for(size_t i=0; i<n_gals; ++i) {
            permuted_global_idxs[i] = (long long)tree_idxs[i];
        }

        // Write to lightcone/data
        auto lc_group = out.group("lightcone");
        append_lightcone_data(lc_group, crds, permuted_global_idxs, galaxy_idxs, subsize, _snapshot_displs.back());
        
        // Write reordered data to data group
        auto data_group = out.group("data");
        unsigned long long offset = _snapshot_displs.back();
        
        // We need to write the data in the SORTED order.
        // The generic writer reads from input (unsorted) and writes to output.
        // We need to pass the permutation vector that maps OutputIndex -> InputIndex.
        // 'idxs' currently holds the KDTree permutation of the ALREADY SORTED data.
        // We need a vector that maps OutputIndex -> OriginalInputIndex.
        // Let's reconstruct it.
        // Step 1: Sort by SAGETreeIndex -> 'sorted_idxs' (maps Sorted -> Original)
        // Step 2: KDTree partition -> 'idxs' (maps KDTree -> Sorted)
        // Combined: 'final_idxs' (maps KDTree -> Original) = sorted_idxs[idxs[i]]
        
        // Re-create the sort permutation
        std::vector<unsigned long long> sorted_idxs(n_gals);
        std::iota(sorted_idxs.begin(), sorted_idxs.end(), 0);
        std::stable_sort(sorted_idxs.begin(), sorted_idxs.end(), [&](unsigned long long a, unsigned long long b) {
            return sage_tree_idxs[a] < sage_tree_idxs[b];
        });
        
        std::vector<unsigned long long> final_idxs(n_gals);
        for(size_t i=0; i<n_gals; ++i) {
            final_idxs[i] = sorted_idxs[idxs[i]];
        }
        
        struct callback_data {
            sageh5tokdtree_app *app;
            hpc::h5::group *snap_group;
            hpc::h5::group *data_group;
            const std::vector<unsigned long long> *idxs;
            unsigned long long offset;
        };
        
        callback_data cb_data = {this, &snap_group, &data_group, &final_idxs, offset};
        
        H5Literate(snap_group.id(), H5_INDEX_NAME, H5_ITER_NATIVE, NULL, 
            [](hid_t g_id, const char *name, const H5L_info_t *info, void *op_data) -> herr_t {
                auto *data = static_cast<callback_data*>(op_data);
                data->app->permute_and_write_generic(*data->snap_group, name, *data->data_group, *data->idxs, data->offset);
                return 0;
            }, &cb_data);
            
        // Write subsize (calculated, already permuted by KDTree)
        append_dataset(data_group, "subsize", H5T_NATIVE_UINT, sizeof(unsigned), subsize.data(), n_gals, offset);
        
        // Write global_index to data group (already permuted)
        append_dataset(data_group, "global_index", H5T_NATIVE_LLONG, sizeof(long long), permuted_global_idxs.data(), n_gals, offset);

        LOGILN("KDTree built with depth ", tree.max_depth());
        
        _total_count += n_gals;
    }
    
    template<typename T>
    void read_dataset(hpc::h5::group &g, const std::string &name, std::vector<T> &vec) {
        if(g.has_link(name)) {
            auto dset = g.dataset(name);
            vec.resize(g.extent(name));
            dset.read(vec);
        } else {
            LOGILN("Dataset ", name, " not found.");
        }
    }

    hpc::fs::path _sage_fn;
    hpc::fs::path _out_fn;
    unsigned      _ppc;
    
    std::vector<unsigned long long> _snapshot_counts;
    std::vector<unsigned long long> _snapshot_displs;
    std::vector<double>             _snapshot_redshifts;
    unsigned long long              _total_count;
    
    // SnapNum -> FileIndex -> Subsize
    std::vector<std::vector<uint32_t>> _subsize_lookup;
};

#define HPC_APP_CLASS sageh5tokdtree_app
#include <libhpc/mpi/main.hh>
