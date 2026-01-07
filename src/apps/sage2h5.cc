#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <boost/regex.hpp>
#include <libhpc/libhpc.hh>
#include <libtao/base/utils.hh>
#include <libtao/base/sage.hh>

class application : public hpc::mpi::application {
  public:
    application(int argc, char *argv[]) : hpc::mpi::application(argc, argv), _comm(&hpc::mpi::comm::world) {
        // Setup some options.
        options().add_options()(
            "mode,m", hpc::po::value<std::string>(&_mode)->default_value("convert"), "mode of operation")(
            "sage,s", hpc::po::value<hpc::fs::path>(&_sage_dir)->required(), "SAGE output directory")(
            "param,p", hpc::po::value<hpc::fs::path>(&_param_fn)->required(), "SAGE parameter file")(
            "alist,a", hpc::po::value<hpc::fs::path>(&_alist_fn)->required(), "SAGE expansion list file")(
            "output,o", hpc::po::value<hpc::fs::path>(&_out_fn), "output file")(
            "verbose,v", hpc::po::value<int>(&_verb)->default_value(0), "verbosity");
        positional_options().add("sage", 1);
        positional_options().add("param", 2);
        positional_options().add("alist", 3);
        positional_options().add("output", 4);

        // Parse options.
        parse_options(argc, argv);
        EXCEPT(!_sage_dir.empty(), "No SAGE output directory given.");
        EXCEPT(!_param_fn.empty(), "No SAGE parameter file given.");
        EXCEPT(!_alist_fn.empty(), "No SAGE expansion file given.");

        // Setup logging.
        hpc::log::levels_type lvl;
        if(_verb == 1)
            lvl = hpc::log::info;
        else if(_verb == 2)
            lvl = hpc::log::debug;
        else if(_verb == 3)
            lvl = hpc::log::trivial;
        if(_verb) {
            if(_comm->size() > 1)
                LOG_PUSH(new hpc::mpi::logger("sage2h5.log.", lvl));
            else
                LOG_PUSH(new hpc::log::stdout(lvl));
        }
    }

    void operator()() {
        if(_mode == "convert")
            convert();
    }

    void convert() {
        EXCEPT(!_out_fn.empty(), "No output file given.");

        _load_param(_param_fn);
        _load_redshifts(_alist_fn);
        hpc::array<size_t, 2> idx_rng = _sage_idx_rng();
        LOGILN("Processing range: ", idx_rng);

        // Sum total trees and galaxies.
        unsigned long long tot_trees = 0, tot_gals = 0;
        {
            LOGBLOCKD("Summing total trees/galaxies.");
            for(size_t ii = idx_rng[0]; ii < idx_rng[1]; ++ii) {
                int                        n_trees;
                std::vector<std::ifstream> files = _open_files(ii);
#ifndef NDEBUG
                int n_c_trees = -1;
#endif
                for(unsigned jj = 0; jj < files.size(); ++jj) {
                    int n_gals;
                    files[jj].read((char *)&n_trees, sizeof(n_trees));
                    files[jj].read((char *)&n_gals, sizeof(n_gals));
                    LOGDLN("File ", jj, " has ", n_gals, " galaxies and ", n_trees, " trees.");
#ifndef NDEBUG
                    if(n_c_trees == -1)
                        n_c_trees = n_trees;
#endif
                    ASSERT(n_trees == n_c_trees, "Number of trees don't match.");

                    // Sum galaxies here.
                    tot_gals += n_gals;
                }

                // Trees go here.
                tot_trees += n_trees;
            }

            // Scan to get my global offsets.
            _goffs[0] = _comm->scan(tot_trees);
            _goffs[1] = _comm->scan(tot_gals);

            // Reduce to get global sizes.
            tot_trees = _comm->all_reduce(tot_trees);
            tot_gals  = _comm->all_reduce(tot_gals);

            LOGILN("Trees:         ", tot_trees);
            LOGILN("Tree offset:   ", _goffs[0]);
            LOGILN("Galaxies:      ", tot_gals);
            LOGILN("Galaxy offset: ", _goffs[1]);
        }

        // Open output file for writing.
        sage::make_hdf5_types(_mem_type, _file_type);
        sage::make_hdf5_sidecar("Millennium", _out_fn.native(), std::set<double>(_redshifts.begin(), _redshifts.end()), _hubble, _box_size);
        _out_file.open(_out_fn.native(), H5F_ACC_TRUNC, *_comm);
        _gals_dset.create(_out_file, "galaxies", _file_type, tot_gals);
        _tree_displs_dset.create(_out_file, "tree_displs", hpc::h5::datatype::native_ullong, tot_trees + 1);
        _tree_cnts_dset.create(_out_file, "tree_counts", hpc::h5::datatype::native_uint, tot_trees);

        // Prepare my buffered output objects.
        _tree_displs_buf.create(
            _tree_displs_dset, hpc::h5::datatype::native_ullong, hpc::h5::buffer_default_size, _goffs[0]);
        _tree_cnts_buf.create(
            _tree_cnts_dset, hpc::h5::datatype::native_ullong, hpc::h5::buffer_default_size, _goffs[0]);
        _gals_buf.create(_gals_dset, _mem_type, hpc::h5::buffer_default_size, _goffs[1]);

        // Process my index range.
        for(size_t ii = idx_rng[0]; ii < idx_rng[1]; ++ii) {
            _fileidx = ii;
            process_index(ii);
        }

        // Be sure to close the buffers to flush them.
        _tree_displs_buf.close();
        _tree_cnts_buf.close();
        _gals_buf.close();

        // Write the final displacement last.
        if(_comm->rank() == 0)
            _tree_displs_dset.write<unsigned long long>(tot_gals, tot_trees);

        // Write out the final information.
        LOGILN("Writing parameter information.");
        std::vector<double> zs(_redshifts.size());
        std::copy(_redshifts.begin(), _redshifts.end(), zs.begin());
        std::reverse(zs.begin(), zs.end());
        _out_file.write_serial("snapshot_redshifts", zs);
        _out_file.write<double>("cosmology/box_size", _box_size);
        _out_file.write<double>("cosmology/hubble", _hubble);
        _out_file.write<double>("cosmology/omega_l", _omega_l);
        _out_file.write<double>("cosmology/omega_m", _omega_m);

        // Close everything down.
        _gals_dset.close();
        _tree_displs_dset.close();
        _tree_cnts_dset.close();
        _out_file.close();

        // Do some final checks.
        check();
    }

    void process_index(size_t idx) {
        LOGBLOCKI("Processing index: ", idx);

        std::vector<std::ifstream> files = _open_files(idx);

        // Skip the number of trees and galaxies, and calculate the
        // counts of galaxies in each tree in this index.
        std::vector<unsigned long long> tree_sizes;
        hpc::matrix<int>                file_tree_sizes;
        int                             n_trees;
        {
            for(unsigned ii = 0; ii < files.size(); ++ii) {
                int n_gals;
                files[ii].read((char *)&n_trees, sizeof(n_trees));
                files[ii].read((char *)&n_gals, sizeof(n_gals));
                LOGDLN("File ", ii, " has ", n_gals, " galaxies and ", n_trees, " trees.");

                // Allocate tree displacements on the first pass.
                if(tree_sizes.empty()) {
                    tree_sizes.resize(n_trees);
                    file_tree_sizes.resize(_redshifts.size(), n_trees);
                    std::fill(tree_sizes.begin(), tree_sizes.end(), 0);
                }

                // Load the sizes and add to displacements.
                files[ii].read((char *)file_tree_sizes[ii].data(), n_trees * sizeof(int));
                for(unsigned jj = 0; jj < n_trees; ++jj)
                    tree_sizes[jj] += file_tree_sizes(ii, jj);
            }
        }

        // Finish calculating displacements.
        std::vector<unsigned long long> tree_displs = hpc::counts_to_displs(tree_sizes);
        std::transform(tree_displs.begin(), tree_displs.end(), tree_displs.begin(), [this](unsigned long long x) {
            return x + _goffs[1];
        });

        // Process in tree major order.
        for(unsigned ii = 0; ii < n_trees; ++ii) {
            LOGBLOCKT("Processing tree: ", ii);
            _treeidx = ii;

            // Keep a mapping to build descendants. This maps from
            // galaxy index to snapshot.
            std::unordered_map<long long, int>                          desc_map;
            std::unordered_multimap<hpc::varray<unsigned, 2>, unsigned> merge_map;

            std::vector<sage::galaxy> h5_gals(tree_sizes[ii]);
            unsigned                  displ = 0;
            for(unsigned jj = 0; jj < files.size(); ++jj) {
                LOGBLOCKT("Processing file: ", jj);
                _filez = jj;

                process_tree_redshift(files[jj], file_tree_sizes(jj, ii), h5_gals, displ, desc_map, merge_map);
            }

            // Reorder the galaxies into depth first.
            depthfirst_ordering(hpc::view<std::vector<sage::galaxy>>(h5_gals));

            // Write galaxy data.
            _gals_buf.write(h5_gals);

            // Update galaxy global index base to the beginning of the
            // next tree.
            _goffs[1] += h5_gals.size();
        }

            // Check the counts and displacements.
#ifndef NDEBUG
        for(unsigned ii = 0; ii < tree_sizes.size(); ++ii) {
            ASSERT(tree_displs[ii + 1] - tree_displs[ii] == tree_sizes[ii], "Sizes and displacements don't match.");
        }
#endif

        // Write out tree details here, after galaxies, because I need
        // to add the global offset to the displacements.
        _tree_displs_buf.write<hpc::view<std::vector<unsigned long long>>>(
            hpc::view<std::vector<unsigned long long>>(tree_displs, n_trees));
        _tree_cnts_buf.write(tree_sizes);

        // Update the global offset values with this index.
        _goffs[0] += n_trees;
    }

    void process_tree_redshift(std::ifstream &                                              file,
                               unsigned                                                     n_file_gals,
                               std::vector<sage::galaxy> &                                  h5_gals,
                               unsigned &                                                   displ,
                               std::unordered_map<long long, int> &                         desc_map,
                               std::unordered_multimap<hpc::varray<unsigned, 2>, unsigned> &merge_map) {
        LOGTLN("Reading ", n_file_gals, " galaxies.");

        // Prepare memory.
        std::vector<GALAXY_OUTPUT> sage_gals(n_file_gals);

        // Load data and convert.
        file.read((char *)sage_gals.data(), sage_gals.size() * sizeof(GALAXY_OUTPUT));
        ASSERT(file.good(), "Error reading SAGE file.");
        for(unsigned ii = 0; ii < sage_gals.size(); ++ii)
            _convert(sage_gals[ii], h5_gals, displ + ii, ii, _goffs[1] + displ + ii, desc_map, merge_map);

        // Update the displacement.
        displ += n_file_gals;
    }

    void depthfirst_ordering(hpc::view<std::vector<sage::galaxy>> gals) {
        // Cache the smallest and largest global index.
        std::array<unsigned long long, 2> gidxs{(unsigned long long)gals.front().global_index,
                                                (unsigned long long)gals.back().global_index};

        // Begin by finding the roots and parents.
        std::list<unsigned>               roots;
        std::multimap<unsigned, unsigned> parents;
        for(unsigned ii = 0; ii < gals.size(); ++ii) {
            if(gals[ii].descendant == -1)
                roots.push_back(ii);
            else
                parents.emplace(gals[ii].descendant, ii);
        }

        // Process roots one at a time to create a new ordering.
        unsigned              dfi = 0;
        std::vector<unsigned> map(gals.size());
        for(auto const &root : roots)
            _depthfirst_recurse(root, gals, map, dfi, parents);
        ASSERT(dfi == gals.size());

        // Permute the galaxies.
        hpc::permute(gals.begin(), gals.end(), map.begin());

        // Reset global indices.
        for(unsigned ii = 0; ii < gals.size(); ++ii) {
            ASSERT(gidxs[0] + ii <= gidxs[1]);
            gals[ii].global_index = gidxs[0] + ii;
        }

        // Remap local indices and descendants.
        for(unsigned ii = 0; ii < gals.size(); ++ii) {
            gals[ii].local_index = map[gals[ii].local_index];
            if(gals[ii].descendant != -1) {
                gals[ii].descendant        = map[gals[ii].descendant];
                gals[ii].global_descendant = gals[gals[ii].descendant].global_index;
            }
        }
    }

    unsigned _depthfirst_recurse(unsigned                                 idx,
                                 hpc::view<std::vector<sage::galaxy>> &   gals,
                                 std::vector<unsigned> &                  map,
                                 unsigned &                               dfi,
                                 std::multimap<unsigned, unsigned> const &parents) {
        map[idx]          = dfi++;
        gals[idx].subsize = 1;
        auto rng          = parents.equal_range(idx);
        while(rng.first != rng.second) {
            gals[idx].subsize += _depthfirst_recurse(rng.first->second, gals, map, dfi, parents);
            ++rng.first;
        }
        return gals[idx].subsize;
    }

    void check() {
        LOGBLOCKI("Checking file.");

        sage::make_hdf5_types(_mem_type, _file_type);
        _out_file.open(_out_fn.native(), H5F_ACC_RDONLY, *_comm);
        _gals_dset.open(_out_file, "galaxies");
        _tree_displs_dset.open(_out_file, "tree_displs");
        _tree_cnts_dset.open(_out_file, "tree_counts");

        unsigned                n_trees  = _tree_cnts_dset.extent();
        hpc::array<unsigned, 2> tree_rng = hpc::mpi::modulo(n_trees);

        LOGILN("Tree range: ", tree_rng);

        for(unsigned ii = tree_rng[0]; ii < tree_rng[1]; ++ii) {
            LOGDLN("Checking tree: ", ii);

            unsigned           size   = _tree_cnts_dset.read<unsigned>(ii);
            unsigned long long displ  = _tree_displs_dset.read<unsigned long long>(ii);
            unsigned long long displ2 = _tree_displs_dset.read<unsigned long long>(ii + 1);
            ASSERT(displ2 - displ == size, "Sizes and displacements don't match.");

            std::vector<sage::galaxy> gals(size);
            _gals_dset.read(gals.data(), _mem_type, size, displ);
            check_tree(gals, displ, ii);
        }
    }



    void check_tree(std::vector<sage::galaxy> const &gals, unsigned long long displ, unsigned tree_idx) {
        std::set<unsigned long long> gids;
        for(unsigned jj = 0; jj < gals.size(); ++jj) {
            gids.insert(gals[jj].global_index);
            ASSERT(gals[jj].local_index == jj, "Invalid local index.");
            ASSERT(gals[jj].global_index >= displ && gals[jj].global_index < displ + gals.size(),
                   "Invalid global index: ",
                   gals[jj].global_index);
            ASSERT(gals[jj].global_descendant == -1 ||
                       (gals[jj].global_index >= displ && gals[jj].global_index < displ + gals.size()),
                   "Invalid global descendant: ",
                   gals[jj].global_descendant);
            ASSERT(gals[jj].descendant == -1 || gals[jj].descendant < gals.size(),
                   "Invalid descendant: ",
                   gals[jj].descendant);
            ASSERT((gals[jj].descendant == -1 && gals[jj].global_descendant == -1) ||
                       (gals[jj].descendant != -1 && gals[jj].global_descendant != -1),
                   "Inconsistent descendant and global descendant indices.");
            // ASSERT( gals[jj].merge_into_id == -1 ||
            // 	 (gals[jj].descendant != -1 &&
            // 	  gals[gals[jj].descendant].snapshot == gals[jj].merge_into_snapshot),
            // 	 "Incorrect merger: merge into ID=", gals[jj].merge_into_id,
            //         ", merge into snapshot=", gals[jj].merge_into_snapshot,
            //         ", global galaxy index=", gals[jj].global_index,
            // 	 ", tree index=", tree_idx, ", local galaxy index=", jj );
            if(!(gals[jj].merge_into_id == -1 ||
                 (gals[jj].descendant != -1 && gals[gals[jj].descendant].snapshot == gals[jj].merge_into_snapshot))) {
                std::cerr << "Warning: incorrect merger detected.\n";
                std::cerr << "  Claimed snapshot = " << gals[jj].merge_into_snapshot << "\n";
                std::cerr << "  Actual snapshot  = " << gals[gals[jj].descendant].snapshot << "\n";
            }
        }
        ASSERT(gids.size() == gals.size(), "Duplicate global indices: ", gids.size(), ", ", gals.size());
        ASSERT(*gids.begin() == displ, "Global indices begin at wrong index.");
        ASSERT(*gids.rbegin() == displ + gals.size() - 1, "Global indices end at wrong index.");
        for(unsigned jj = 0; jj < gals.size(); ++jj) {
            if(gals[jj].global_descendant != -1) {
                ASSERT(hpc::has(gids, gals[jj].global_descendant),
                       "Invalid global descendant: ",
                       gals[jj].global_descendant,
                       " from global index: ",
                       gals[jj].global_index);
            }
        }
    }

  protected:
    void _convert(GALAXY_OUTPUT const &                                        sage_gal,
                  std::vector<sage::galaxy> &                                  h5_gals,
                  unsigned                                                     gal_idx,
                  unsigned                                                     file_gal_idx,
                  unsigned long long                                           gidx,
                  std::unordered_map<long long, int> &                         desc_map,
                  std::unordered_multimap<hpc::varray<unsigned, 2>, unsigned> &merge_map) {
        LOGBLOCKT("Converting galaxy: ", gal_idx);

        ASSERT(sage_gal.Type < 2, "Found an invalid SAGE galaxy type: ", sage_gal.Type);

        h5_gals[gal_idx].snapshot = sage_gal.SnapNum;
        h5_gals[gal_idx].type     = sage_gal.Type;

        h5_gals[gal_idx].galaxy_idx              = sage_gal.GalaxyIndex;
        h5_gals[gal_idx].central_galaxy_idx      = sage_gal.CentralGalaxyIndex;
        h5_gals[gal_idx].sage_halo_idx           = sage_gal.SAGEHaloIndex;
        h5_gals[gal_idx].sage_tree_idx           = sage_gal.SAGETreeIndex;
        h5_gals[gal_idx].simulation_halo_idx = sage_gal.SimulationHaloIndex;

        h5_gals[gal_idx].local_index       = gal_idx;
        h5_gals[gal_idx].global_index      = gidx;
        h5_gals[gal_idx].descendant        = -1;
        h5_gals[gal_idx].global_descendant = -1;

        // Add to the merge map, if needed.
        if(sage_gal.mergeIntoID != -1) {
            hpc::varray<unsigned, 2> id{(unsigned)sage_gal.mergeIntoID, (unsigned)sage_gal.mergeIntoSnapNum};
            LOGDLN("Adding merge to map: ", id);
            merge_map.emplace(id, gal_idx);
        }

        ASSERT(h5_gals[gal_idx].global_index >= _goffs[1] && h5_gals[gal_idx].global_index < _goffs[1] + h5_gals.size(),
               "Invalid global index.");

        LOGTLN("Galaxy index: ", h5_gals[gal_idx].galaxy_idx);

        // Check for a parent.
        if(hpc::has(desc_map, h5_gals[gal_idx].galaxy_idx)) {
            int par = desc_map[h5_gals[gal_idx].galaxy_idx];

            ASSERT(par != -1 && h5_gals[par].descendant == -1, "Found a galaxy with a parent that has merged.");

            h5_gals[par].descendant        = gal_idx;
            h5_gals[par].global_descendant = h5_gals[gal_idx].global_index;

            ASSERT(h5_gals[par].global_descendant == -1 ||
                       (h5_gals[par].global_descendant >= _goffs[1] &&
                        h5_gals[par].global_descendant < _goffs[1] + h5_gals.size()),
                   "Invalid global descendant.");

            ASSERT(h5_gals[par].snapshot < sage_gal.SnapNum, "Out of order snapshot indices.");

            LOGTLN("Parent: ", par);
        }

        // Check for a merge.
        hpc::varray<unsigned, 2> merge_id{(unsigned)file_gal_idx, (unsigned)sage_gal.SnapNum};
        if(hpc::has(merge_map, merge_id)) {
            auto rng = merge_map.equal_range(merge_id);
            while(rng.first != rng.second) {
                int par = rng.first->second;
                ++rng.first;

                ASSERT(par != -1 && h5_gals[par].descendant == -1, "Found a galaxy with a parent that has merged.");

                h5_gals[par].descendant        = gal_idx;
                h5_gals[par].global_descendant = h5_gals[gal_idx].global_index;

                ASSERT(h5_gals[par].global_descendant == -1 ||
                           (h5_gals[par].global_descendant >= _goffs[1] &&
                            h5_gals[par].global_descendant < _goffs[1] + h5_gals.size()),
                       "Invalid global descendant.");

                ASSERT(h5_gals[par].snapshot < sage_gal.SnapNum, "Out of order snapshot indices.");

                LOGTLN("Parent: ", par);
            }

            LOGDLN("Erasing merge from map: ", merge_id);
            merge_map.erase(merge_id);
        }

        h5_gals[gal_idx].merge_type          = sage_gal.mergeType;
        h5_gals[gal_idx].merge_into_id       = sage_gal.mergeIntoID;
        h5_gals[gal_idx].merge_into_snapshot = sage_gal.mergeIntoSnapNum;
        h5_gals[gal_idx].dt                  = sage_gal.dT; // Preserve SAGE units (Myr)

        ASSERT(h5_gals[gal_idx].merge_into_id == -1 || h5_gals[gal_idx].merge_into_id < h5_gals.size(),
               "Invalid merged-into ID: ",
               h5_gals[gal_idx].merge_into_id);

        // If we merged indicate so in the map.
        if(h5_gals[gal_idx].merge_into_id != -1) {
            desc_map[h5_gals[gal_idx].galaxy_idx] = -1;

            LOGTLN("Merged into: ", h5_gals[gal_idx].merge_into_id);
        } else
            desc_map[h5_gals[gal_idx].galaxy_idx] = gal_idx;

        // printf("sage2h5[gal_idx=%u]: pos=(%f,%f,%f), vel=(%f,%f,%f), spin=(%f,%f,%f)\n",
        //        gal_idx,
        //        sage_gal.Pos[0],
        //        sage_gal.Pos[1],
        //        sage_gal.Pos[2],
        //        sage_gal.Vel[0],
        //        sage_gal.Vel[1],
        //        sage_gal.Vel[2],
        //        sage_gal.Spin[0],
        //        sage_gal.Spin[1],
        //        sage_gal.Spin[2]);
            for(unsigned kk = 0; kk < 3; ++kk) {
            h5_gals[gal_idx].pos[kk]  = sage_gal.Pos[kk];
            h5_gals[gal_idx].vel[kk]  = sage_gal.Vel[kk];
            h5_gals[gal_idx].spin[kk] = sage_gal.Spin[kk];

            ASSERT(h5_gals[gal_idx].pos[kk] == h5_gals[gal_idx].pos[kk],
                   "Bad spin: ",
                   "file index=",
                   _fileidx,
                   ", file z=",
                   _filez,
                   ", tree index=",
                   _treeidx,
                   ", galaxy index=",
                   gal_idx);
            ASSERT(h5_gals[gal_idx].vel[kk] == h5_gals[gal_idx].vel[kk],
                   "Bad spin: ",
                   "file index=",
                   _fileidx,
                   ", file z=",
                   _filez,
                   ", tree index=",
                   _treeidx,
                   ", galaxy index=",
                   gal_idx);
            // TODO: Bolshoi has NaN spins for some reason...
            // ASSERT( h5_gals[gal_idx].spin[kk] == h5_gals[gal_idx].spin[kk], "Bad spin: ",
            // 	 "file index=", _fileidx, ", file z=", _filez, ", tree index=", _treeidx,
            // 	 ", galaxy index=", file_gal_idx );
        }
        h5_gals[gal_idx].num_particles = sage_gal.Len;
        h5_gals[gal_idx].mvir          = sage_gal.Mvir;
        h5_gals[gal_idx].central_mvir  = sage_gal.CentralMvir;
        h5_gals[gal_idx].rvir          = sage_gal.Rvir;
        h5_gals[gal_idx].vvir          = sage_gal.Vvir;
        h5_gals[gal_idx].vmax          = sage_gal.Vmax;
        h5_gals[gal_idx].vel_disp      = sage_gal.VelDisp;

        h5_gals[gal_idx].cold_gas       = sage_gal.ColdGas;
        h5_gals[gal_idx].stellar_mass   = sage_gal.StellarMass;
        h5_gals[gal_idx].bulge_mass     = sage_gal.BulgeMass;
        h5_gals[gal_idx].hot_gas        = sage_gal.HotGas;
        h5_gals[gal_idx].ejected_mass   = sage_gal.EjectedMass;
        h5_gals[gal_idx].blackhole_mass = sage_gal.BlackHoleMass;
        h5_gals[gal_idx].ics            = sage_gal.ICS;

        h5_gals[gal_idx].metals_cold_gas     = sage_gal.MetalsColdGas;
        h5_gals[gal_idx].metals_stellar_mass = sage_gal.MetalsStellarMass;
        h5_gals[gal_idx].metals_bulge_mass   = sage_gal.MetalsBulgeMass;
        h5_gals[gal_idx].metals_hot_gas      = sage_gal.MetalsHotGas;
        h5_gals[gal_idx].metals_ejected_mass = sage_gal.MetalsEjectedMass;
        h5_gals[gal_idx].metals_ics          = sage_gal.MetalsICS;

        h5_gals[gal_idx].sfr_disk    = sage_gal.SfrDisk;
        h5_gals[gal_idx].sfr_bulge   = sage_gal.SfrBulge;
        h5_gals[gal_idx].sfr_disk_z  = sage_gal.SfrDiskZ;
        h5_gals[gal_idx].sfr_bulge_z = sage_gal.SfrBulgeZ;

        h5_gals[gal_idx].disk_scale_radius             = sage_gal.DiskScaleRadius;
        h5_gals[gal_idx].cooling                       = sage_gal.Cooling;
        h5_gals[gal_idx].heating                       = sage_gal.Heating;
        h5_gals[gal_idx].quasar_mode_bh_accretion_mass = sage_gal.QuasarModeBHaccretionMass;
        h5_gals[gal_idx].time_of_last_major_merger       = sage_gal.TimeOfLastMajorMerger;
        h5_gals[gal_idx].time_of_last_minor_merger       = sage_gal.TimeOfLastMinorMerger;
        h5_gals[gal_idx].outflow_rate                  = sage_gal.OutflowRate;

        h5_gals[gal_idx].infall_mvir = sage_gal.infallMvir;
        h5_gals[gal_idx].infall_vvir = sage_gal.infallVvir;
        h5_gals[gal_idx].infall_vmax = sage_gal.infallVmax;
    }

    std::vector<std::ifstream> _open_files(size_t idx) {
        std::vector<std::ifstream> files(_redshifts.size());
        size_t                     ii = files.size() - 1;
        for(auto z : _redshifts) {
            hpc::fs::path fn = _make_filename(idx, z);
            LOGDLN("Opening file: ", fn);
            printf("Opening file: %s\n", fn.c_str());
            files[ii].open(fn.native(), std::ios::binary);
            EXCEPT(files[ii].good(), "Failed to open file: ", fn);
            --ii;
        }
        return files;
    }

    void _load_param(hpc::fs::path const &fn) {
        LOGBLOCKD("Reading parameter file: ", fn);

        std::ifstream file(fn.native());
        EXCEPT(file.good(), "Failed to find parameter file: ", fn);

        boost::regex  first_prog("\\s*FirstFile\\s+(\\d+)\\s*");
        boost::regex  last_prog("\\s*LastFile\\s+(\\d+)\\s*");
        //boost::regex  box_size_prog("\\s*%\\s*BoxSize\\s+(\\d+[.]?\\d*)\\s*");
        boost::regex  box_size_prog("\\s*BoxSize\\s+(\\d+[.]?\\d*).*");
        boost::regex  hubble_prog("\\s*Hubble_h\\s+(\\d+[.]?\\d*).*");
        boost::regex  omega_m_prog("\\s*Omega\\s+(\\d+[.]?\\d*)\\s*");
        boost::regex  omega_l_prog("\\s*OmegaLambda\\s+(\\d+[.]?\\d*)\\s*");
        boost::cmatch match;

        _box_size = std::numeric_limits<double>::max();

        std::string line;
        int         both = 0;
        while(file.good()) {
            std::getline(file, line);
            printf("Processing line: %s\n", line.c_str());
            if(boost::regex_match(line.c_str(), match, first_prog)) {
                printf("matches first_prog: %s\n", match[1].first);
                _idx_rng[0] = boost::lexical_cast<size_t>(match[1].first, match[1].second - match[1].first);
                ++both;
            }
            if(boost::regex_match(line.c_str(), match, last_prog)) {
                printf("matches last_prog: %s\n", match[1].first);
                EXCEPT(_idx_rng[0] != std::numeric_limits<size_t>::max(), "Found LastFile before FirstFile.");
                _idx_rng[1] = boost::lexical_cast<size_t>(match[1].first, match[1].second - match[1].first) + 1;
                ++both;
            }
            if(boost::regex_match(line.c_str(), match, box_size_prog)) {
                printf("matches box_size_prog: %s\n", match[1].first);
                EXCEPT(_box_size == std::numeric_limits<double>::max(), "Found BoxSize multiple times.");
                _box_size = boost::lexical_cast<double>(match[1].first, match[1].second - match[1].first);
            }
            if(boost::regex_match(line.c_str(), match, hubble_prog)) {
                printf("matches hubble_prog: %s\n", match[1].first);
                EXCEPT(_hubble == 0.0, "Found Hubble_h multiple times.");
                _hubble = boost::lexical_cast<double>(match[1].first, match[1].second - match[1].first);
                _hubble *= 100.0;
            }
            if(boost::regex_match(line.c_str(), match, omega_m_prog)) {
                printf("matches omega_m_prog: %s\n", match[1].first);
                EXCEPT(_omega_m == 0.0, "Found Omega multiple times.");
                _omega_m = boost::lexical_cast<double>(match[1].first, match[1].second - match[1].first);
            }
            if(boost::regex_match(line.c_str(), match, omega_l_prog)) {
                printf("matches omega_l_prog: %s\n", match[1].first);
                _omega_l = boost::lexical_cast<double>(match[1].first, match[1].second - match[1].first);
            }
        }

        EXCEPT(both == 2, "Failed to find FirstFile and LastFile in parameters.");
        EXCEPT(_idx_rng[1] > _idx_rng[0], "Invalid index range.");
        EXCEPT(_box_size != std::numeric_limits<double>::max(), "Failed to find box size.");

        LOGDLN("Have index range: ", _idx_rng);
    }

    void _load_redshifts(hpc::fs::path const &fn) {
        LOGBLOCKD("Reading alist file: ", fn);

        std::ifstream file(fn.native());
        EXCEPT(file.good(), "Failed to open expansion file: ", fn);

        double exp;
        while(file.good()) {
            file >> exp;
            _redshifts.insert(tao::expansion_to_redshift(exp));
        }

        LOGDLN("Have redshifts: ", _redshifts);
    }

    hpc::fs::path _make_filename(size_t idx, double z) const {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(3) << "model_z" << z << "_" << idx;
        return _sage_dir / ss.str();
    }

    hpc::array<size_t, 2> _sage_idx_rng() {
        // Decide on my index range.
        hpc::array<size_t, 2> idx_rng = hpc::mpi::modulo(_idx_rng[1] - _idx_rng[0]);
        idx_rng[0] += _idx_rng[0];
        idx_rng[1] += _idx_rng[0];
        return idx_rng;
    }

  protected:
    std::set<double>                  _redshifts;
    double                            _box_size;
    double                            _hubble;
    double                            _omega_l;
    double                            _omega_m;
    hpc::array<size_t, 2>             _idx_rng;
    hpc::array<unsigned long long, 2> _goffs;

    hpc::h5::file                       _out_file;
    hpc::h5::datatype                   _mem_type;
    hpc::h5::datatype                   _file_type;
    hpc::h5::dataset                    _gals_dset;
    hpc::h5::dataset                    _tree_displs_dset;
    hpc::h5::dataset                    _tree_cnts_dset;
    hpc::h5::buffer<sage::galaxy>       _gals_buf;
    hpc::h5::buffer<unsigned long long> _tree_displs_buf;
    hpc::h5::buffer<unsigned long long> _tree_cnts_buf;

    std::string   _mode;
    hpc::fs::path _sage_dir;
    hpc::fs::path _param_fn;
    hpc::fs::path _alist_fn;
    hpc::fs::path _out_fn;
    unsigned      _treeidx;
    unsigned      _lidx;
    unsigned      _fileidx;
    unsigned      _filez;
    int           _verb;

    hpc::mpi::comm const *_comm;
};

#define HPC_APP_CLASS application
#include <libhpc/mpi/main.hh>
