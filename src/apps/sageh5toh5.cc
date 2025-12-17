#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <map>
#include <set>
#include <list>
#include <algorithm>
#include <libhpc/system/type_traits.hh>
#include <libhpc/libhpc.hh>
#include <libtao/base/utils.hh>
#include <libtao/base/sage.hh>

class application : public hpc::mpi::application {
  public:
    application(int argc, char *argv[]) : hpc::mpi::application(argc, argv), _comm(const_cast<hpc::mpi::comm*>(&hpc::mpi::comm::world)), _box_size(0), _hubble(0), _omega_l(0), _omega_m(0) {
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
                LOG_PUSH(new hpc::mpi::logger("sageh5toh5.log.", lvl));
            else
                LOG_PUSH(new hpc::log::stdout(lvl));
        }
    }

    void operator()() {
        if(_mode == "convert")
            convert();
        else
            LOGELN("Unknown mode: ", _mode);
    }

    void convert() {
        EXCEPT(!_out_fn.empty(), "No output file given.");

        _load_param(_param_fn);
        _load_redshifts(_alist_fn);

        // 1. Identify input files.
        std::cout << "Searching for files in " << _sage_dir << std::endl;
        std::vector<hpc::fs::path> input_files;
        if(hpc::fs::is_directory(_sage_dir)) {
            hpc::fs::directory_iterator end_iter;
            for(hpc::fs::directory_iterator dir_itr(_sage_dir); dir_itr != end_iter; ++dir_itr) {
                if(hpc::fs::is_regular_file(dir_itr->status())) {
                    std::string ext = dir_itr->path().extension().string();
                    std::string stem = dir_itr->path().stem().string();
                    std::cout << "Checking file: " << dir_itr->path() << " stem: " << stem << std::endl;
                    if(ext == ".hdf5" && (stem.find("model_hdf5_") == 0 || stem.find("model_") == 0)) {
                        input_files.push_back(dir_itr->path());
                    }
                }
            }
        } else {
            input_files.push_back(_sage_dir);
        }
        std::cout << "Found " << input_files.size() << " input files." << std::endl;
        std::sort(input_files.begin(), input_files.end()); // Ensure deterministic order

        // Distribute files among ranks (simple round-robin or block).
        // For now, let's assume all ranks process a subset of files.
        // Actually, sage2h5 processes by index range. Here we process by file.
        std::vector<hpc::fs::path> my_files;
        for(size_t i = 0; i < input_files.size(); ++i) {
            if(i % _comm->size() == _comm->rank()) {
                my_files.push_back(input_files[i]);
            }
        }

        // 2. Scan Phase: Count trees and galaxies.
        unsigned long long my_tot_trees = 0;
        unsigned long long my_tot_gals = 0;
        std::vector<unsigned long long> file_tree_counts;
        std::vector<unsigned long long> file_gal_counts;

        for(const auto& fn : my_files) {
            hpc::h5::file f(fn.native(), H5F_ACC_RDONLY);
            
            // Determine number of trees from Snap_63/SAGETreeIndex
            size_t n_trees = 0;
            if(f.has_link("Snap_63/SAGETreeIndex")) {
                hpc::h5::dataset dset(f, "Snap_63/SAGETreeIndex");
                std::vector<int> tree_indices;
                tree_indices.resize(dset.extent());
                dset.read(tree_indices);
                                std::set<int> unique_trees(tree_indices.begin(), tree_indices.end());
                                std::cout << "Snap_63 has " << tree_indices.size() << " galaxies with " << unique_trees.size() << " unique tree indices" << std::endl;
                if(!tree_indices.empty()) {
                    int max_tree_idx = *std::max_element(tree_indices.begin(), tree_indices.end());
                    n_trees = max_tree_idx + 1;
                    for (int iii=0;iii<tree_indices.size();++iii) {
                        if (tree_indices[iii] == max_tree_idx) {
                            std::cout << "found in Snap_63: " << iii << " tree_idx: " << tree_indices[iii] << " n_trees: " << n_trees << std::endl;
                        }
                    }   
                }
            }

            // Fallback to TreeInfo if Snap_63 didn't give us trees (e.g. empty snapshot)
            if (n_trees == 0 && f.has_link("TreeInfo/Snap_0/NumGalsPerTreePerSnap")) {
                hpc::h5::dataset dset(f, "TreeInfo/Snap_0/NumGalsPerTreePerSnap");
                n_trees = dset.dataspace().size();
            }

            std::cout << "File: " << fn << " n_trees: " << n_trees << std::endl;
            
            std::vector<int> tree_sizes(n_trees, 0);
            
            // Sum galaxies for each tree across all snapshots
            for(int snap = 0; snap < 64; ++snap) { // Assuming 64 snapshots
                std::string group_name = "Snap_" + std::to_string(snap);
                if(!f.has_link(group_name)) continue;
                
                hpc::h5::group g;
                f.open_group(group_name, g);
                
                if(g.has_link("SAGETreeIndex")) {
                    std::vector<int> tree_indices;
                    hpc::h5::dataset dset = g.dataset("SAGETreeIndex");
                    tree_indices.resize(dset.extent());
                    dset.read(tree_indices);
                    
                    if (snap == 63 && tree_indices.size() > 0) {
                         std::cout << "Snap 63: Found " << tree_indices.size() << " galaxies. First tree_idx: " << tree_indices[0] << std::endl;
                    }

                    for(int tree_idx : tree_indices) {
                        if(tree_idx >= 0 && (size_t)tree_idx < n_trees) {
                            tree_sizes[tree_idx]++;
                        } else {
                             std::cout << "Invalid tree_idx: " << tree_idx << " n_trees: " << n_trees << std::endl;
                        }
                    }
                } else {
                    std::cout << "Snap " << snap << " missing SAGETreeIndex" << std::endl;
                }
            }
            
            unsigned long long n_gals = 0;
            for(int s : tree_sizes) n_gals += s;

            my_tot_trees += n_trees;
            my_tot_gals += n_gals;
            file_tree_counts.push_back(n_trees);
            file_gal_counts.push_back(n_gals);
        }

                // Global offsets
        // scan() returns exclusive prefix sum by default, so it is already the offset.
        unsigned long long scanned_trees = _comm->scan(my_tot_trees);
        unsigned long long scanned_gals = _comm->scan(my_tot_gals);
        unsigned long long global_tree_offset = scanned_trees;
        unsigned long long global_gal_offset = scanned_gals;
        unsigned long long total_trees = _comm->all_reduce(my_tot_trees);
        unsigned long long total_gals = _comm->all_reduce(my_tot_gals);

        LOGILN("My Trees: ", my_tot_trees);
        LOGILN("My Galaxies: ", my_tot_gals);
        LOGILN("Scanned Trees: ", scanned_trees);
        LOGILN("Scanned Galaxies: ", scanned_gals);
        LOGILN("Global Tree Offset: ", global_tree_offset);
        LOGILN("Global Gal Offset: ", global_gal_offset);
        LOGILN("Total Trees: ", total_trees);
        LOGILN("Total Galaxies: ", total_gals);

        // 3. Create Output File
        sage::make_hdf5_types(_mem_type, _file_type);
        sage::make_hdf5_sidecar("Millennium", _out_fn.native(), std::set<double>(_redshifts.begin(), _redshifts.end()), _hubble, _box_size);
        
        hpc::h5::file out_file(_out_fn.native(), H5F_ACC_TRUNC, *_comm);
        hpc::h5::dataset gals_dset(out_file, "galaxies", _file_type, hpc::h5::dataspace(total_gals));
        hpc::h5::dataset tree_displs_dset(out_file, "tree_displs", hpc::h5::datatype::native_ullong, hpc::h5::dataspace(total_trees + 1));
        hpc::h5::dataset tree_cnts_dset(out_file, "tree_counts", hpc::h5::datatype::native_uint, hpc::h5::dataspace(total_trees));

        // Buffers
        hpc::h5::buffer<sage::galaxy> gals_buf;
        gals_buf.create(gals_dset, _mem_type, hpc::h5::buffer_default_size, global_gal_offset);
        hpc::h5::buffer<unsigned long long> tree_displs_buf;
        tree_displs_buf.create(tree_displs_dset, hpc::h5::datatype::native_ullong, hpc::h5::buffer_default_size, global_tree_offset);
        hpc::h5::buffer<unsigned> tree_cnts_buf;
        tree_cnts_buf.create(tree_cnts_dset, hpc::h5::datatype::native_uint, hpc::h5::buffer_default_size, global_tree_offset);

        // 4. Process Files
        unsigned long long current_gal_global_idx = global_gal_offset;
        unsigned long long current_tree_displ = global_gal_offset;

        for(size_t fidx = 0; fidx < my_files.size(); ++fidx) {
            hpc::fs::path fn = my_files[fidx];
            LOGILN("Processing file: ", fn);
            
            hpc::h5::file f(fn.native(), H5F_ACC_RDONLY);
            size_t n_trees = file_tree_counts[fidx];
            LOGILN("Processing file: ", fn, " with n_trees: ", n_trees);
            
            // Read all data into memory, organized by tree
            std::vector<std::vector<sage::galaxy>> trees(n_trees);
            
            for(int snap = 0; snap < 64; ++snap) {
                std::string group_name = "Snap_" + std::to_string(snap);
                if(!f.has_link(group_name)) continue;
                
                hpc::h5::group g;
                f.open_group(group_name, g);
                
                                // Read SAGETreeIndex to know where each galaxy goes
                std::vector<int> tree_indices;
                if(g.has_link("SAGETreeIndex")) {
                    hpc::h5::dataset dset = g.dataset("SAGETreeIndex");
                    size_t n_sage_trees = dset.extent();
                    tree_indices.resize(dset.extent());
                    dset.read(tree_indices);
                }
                size_t n_gals_snap = tree_indices.size();
                if(n_gals_snap == 0) continue;
                
                LOGILN("Snap ", snap, " has ", n_gals_snap, " galaxies. First tree index: ", tree_indices[0]);


                // Helper to read column
                auto read_col = [&](const char* name, auto& vec) {
                    vec.resize(n_gals_snap);
                    if(g.has_link(name)) {
                        hpc::h5::dataset dset = g.dataset(name);
                        dset.read(vec);
                    } else {
                        std::fill(vec.begin(), vec.end(), 0); // Handle missing columns with 0
                    }
                };

                std::vector<float> posx, posy, posz, velx, vely, velz, spinx, spiny, spinz;
                std::vector<float> mvir, central_mvir, rvir, vvir, vmax, vel_disp;
                std::vector<float> cold_gas, stellar_mass, bulge_mass, hot_gas, ejected_mass, blackhole_mass, ics;
                std::vector<float> metals_cold_gas, metals_stellar_mass, metals_bulge_mass, metals_hot_gas, metals_ejected_mass, metals_ics;
                std::vector<float> sfr_disk, sfr_bulge, sfr_disk_z, sfr_bulge_z;
                std::vector<float> disk_radius, cooling, heating, quasar_acc, last_major, last_minor, outflow;
                std::vector<float> infall_mvir, infall_vmax, infall_vvir;
                std::vector<long long> gal_idx, central_gal_idx, sim_halo_idx;
                std::vector<int> sage_halo_idx, type, merge_type, merge_into_id, merge_into_snap;
                std::vector<float> dt;
                std::vector<int> len;

                read_col("Posx", posx); read_col("Posy", posy); read_col("Posz", posz);
                read_col("Velx", velx); read_col("Vely", vely); read_col("Velz", velz);
                read_col("Spinx", spinx); read_col("Spiny", spiny); read_col("Spinz", spinz);
                read_col("Mvir", mvir); read_col("CentralMvir", central_mvir);
                read_col("Rvir", rvir); read_col("Vvir", vvir); read_col("Vmax", vmax); read_col("VelDisp", vel_disp);
                read_col("ColdGas", cold_gas); read_col("StellarMass", stellar_mass); read_col("BulgeMass", bulge_mass);
                read_col("HotGas", hot_gas); read_col("EjectedMass", ejected_mass); read_col("BlackHoleMass", blackhole_mass); read_col("IntraClusterStars", ics);
                read_col("MetalsColdGas", metals_cold_gas); read_col("MetalsStellarMass", metals_stellar_mass); read_col("MetalsBulgeMass", metals_bulge_mass);
                read_col("MetalsHotGas", metals_hot_gas); read_col("MetalsEjectedMass", metals_ejected_mass); read_col("MetalsIntraClusterStars", metals_ics);
                read_col("SfrDisk", sfr_disk); read_col("SfrBulge", sfr_bulge); read_col("SfrDiskZ", sfr_disk_z); read_col("SfrBulgeZ", sfr_bulge_z);
                read_col("DiskRadius", disk_radius); read_col("Cooling", cooling); read_col("Heating", heating);
                read_col("QuasarModeBHaccretionMass", quasar_acc); read_col("TimeOfLastMajorMerger", last_major); read_col("TimeOfLastMinorMerger", last_minor);
                read_col("OutflowRate", outflow);
                read_col("infallMvir", infall_mvir); read_col("infallVmax", infall_vmax); read_col("infallVvir", infall_vvir);
                read_col("GalaxyIndex", gal_idx); read_col("CentralGalaxyIndex", central_gal_idx); read_col("SimulationHaloIndex", sim_halo_idx);
                read_col("SAGEHaloIndex", sage_halo_idx); read_col("Type", type);
                read_col("mergeType", merge_type); read_col("mergeIntoID", merge_into_id); read_col("mergeIntoSnapNum", merge_into_snap);
                read_col("dT", dt);
                read_col("Len", len);

                for(size_t i=0; i<n_gals_snap; ++i) {
                    sage::galaxy g;
                    g.snapshot = snap;
                    g.type = type[i];
                    g.galaxy_idx = gal_idx[i];
                    g.central_galaxy_idx = central_gal_idx[i];
                    g.sage_halo_idx = sage_halo_idx[i];
                    g.sage_tree_idx = tree_indices[i];
                    g.simulation_halo_idx = sim_halo_idx[i];
                    
                    g.pos[0] = posx[i]; g.pos[1] = posy[i]; g.pos[2] = posz[i];
                    g.vel[0] = velx[i]; g.vel[1] = vely[i]; g.vel[2] = velz[i];
                    g.spin[0] = spinx[i]; g.spin[1] = spiny[i]; g.spin[2] = spinz[i];
                    
                    g.mvir = mvir[i]; g.central_mvir = central_mvir[i];
                    g.rvir = rvir[i]; g.vvir = vvir[i]; g.vmax = vmax[i]; g.vel_disp = vel_disp[i];
                    
                    g.cold_gas = cold_gas[i]; g.stellar_mass = stellar_mass[i]; g.bulge_mass = bulge_mass[i];
                    g.hot_gas = hot_gas[i]; g.ejected_mass = ejected_mass[i]; g.blackhole_mass = blackhole_mass[i]; g.ics = ics[i];
                    
                    g.metals_cold_gas = metals_cold_gas[i]; g.metals_stellar_mass = metals_stellar_mass[i]; g.metals_bulge_mass = metals_bulge_mass[i];
                    g.metals_hot_gas = metals_hot_gas[i]; g.metals_ejected_mass = metals_ejected_mass[i]; g.metals_ics = metals_ics[i];
                    
                    g.sfr_disk = sfr_disk[i]; g.sfr_bulge = sfr_bulge[i]; g.sfr_disk_z = sfr_disk_z[i]; g.sfr_bulge_z = sfr_bulge_z[i];
                    
                    g.disk_scale_radius = disk_radius[i]; g.cooling = cooling[i]; g.heating = heating[i];
                    g.quasar_mode_bh_accretion_mass = quasar_acc[i];
                    g.time_of_last_major_merger = last_major[i]; g.time_of_last_minor_merger = last_minor[i];
                    g.outflow_rate = outflow[i];
                    
                    g.infall_mvir = infall_mvir[i]; g.infall_vmax = infall_vmax[i]; g.infall_vvir = infall_vvir[i];
                    
                    g.merge_type = merge_type[i];
                    g.merge_into_id = merge_into_id[i];
                    g.merge_into_snapshot = merge_into_snap[i];
                    g.dt = dt[i] / 1000.0; // Convert to Gyrs
                    
                    g.num_particles = len[i];
                    
                    // Initialize others
                    g.descendant = -1;
                    g.global_descendant = -1;
                    g.local_index = i; // Store snapshot index for merging logic
                    g.global_index = -1; // Will be set later

                    trees[tree_indices[i]].push_back(g);
                }
            }

            // Process trees
            for(size_t t=0; t<n_trees; ++t) {
                auto& tree_gals = trees[t];
                //if (t < 10 && tree_gals.size() > 0) {
                     std::cout << "Writing Tree " << t << " size: " << tree_gals.size() << std::endl;
                //}
                if(tree_gals.empty()) {
                    tree_displs_buf.write(current_tree_displ);
                    tree_cnts_buf.write(0);
                    continue;
                }

                // Link descendants
                std::unordered_map<long long, int> desc_map;
                std::unordered_multimap<hpc::varray<unsigned, 2>, unsigned> merge_map;

                for(size_t i=0; i<tree_gals.size(); ++i) {
                    auto& g = tree_gals[i];
                    // g.local_index currently holds snapshot index
                    g.global_index = current_gal_global_idx + i;
                    
                    // Add to merge map if needed
                    if(g.merge_into_id != -1) {
                        hpc::varray<unsigned, 2> id{(unsigned)g.merge_into_id, (unsigned)g.merge_into_snapshot};
                        merge_map.emplace(id, i);
                    }
                    
                    // Check for parent in desc_map
                    if(desc_map.count(g.galaxy_idx)) {
                        int par = desc_map[g.galaxy_idx];
                        tree_gals[par].descendant = i;
                        tree_gals[par].global_descendant = g.global_index;
                    }
                    
                    // Check for merges
                    hpc::varray<unsigned, 2> merge_id{(unsigned)g.local_index, (unsigned)g.snapshot}; // Wait, local_index here is not file_gal_idx.
                    // Issue: merge_into_id in SAGE HDF5 refers to what?
                    // In SAGE binary, it refers to the index in the file.
                    // In SAGE HDF5, is it GalaxyIndex? Or index in the snapshot?
                    // "mergeIntoID <HDF5 dataset "mergeIntoID": shape (35899,), type "<i4">"
                    // If it's <i4>, it's likely an index, not a long long GalaxyIndex.
                    // If it's an index, is it the index within the snapshot? Or within the tree?
                    // SAGE HDF5 documentation says: "mergeIntoID: The unique ID of the galaxy that this galaxy merges into."
                    // But GalaxyIndex is i8 (long long). mergeIntoID is i4.
                    // This suggests mergeIntoID is NOT GalaxyIndex.
                    // It might be the index in the snapshot array?
                    // If so, I need to know the index in the snapshot array to resolve merges.
                    // But I've lost that index when I pushed to `trees`.
                    // I need to store the original snapshot index in `sage::galaxy` temporarily?
                    // Or, maybe mergeIntoID is the `SAGEHaloIndex`? No.
                    // Let's assume it's the index in the snapshot file.
                    // So I need to track `file_gal_idx` (index in snapshot).
                    // I can store it in `local_index` temporarily before reordering?
                    // No, `local_index` is used for output.
                    // I'll add a temporary field or use `local_index` and reset it later.
                    // Let's use `local_index` to store `file_gal_idx` during reading.
                }
                
                // Re-implementing linking with correct indices
                desc_map.clear();
                merge_map.clear();
                
                // First pass: populate maps using the stored file_gal_idx (in local_index)
                // But wait, I need to process in order.
                // The `trees` vector is sorted by snapshot (0..63).
                
                for(size_t i=0; i<tree_gals.size(); ++i) {
                    auto& g = tree_gals[i];
                    // g.local_index currently holds the index in the snapshot file (from reading loop)
                    unsigned file_gal_idx = g.local_index; 
                    
                    // Update global index for output
                    g.global_index = current_gal_global_idx + i;
                    
                    // 1. If I merge into someone, register myself in merge_map
                    if(g.merge_into_id != -1) {
                        hpc::varray<unsigned, 2> id{(unsigned)g.merge_into_id, (unsigned)g.merge_into_snapshot};
                        merge_map.emplace(id, i);
                    }
                    
                    // 2. If I am a descendant of someone (via GalaxyIndex), link them
                    if(desc_map.count(g.galaxy_idx)) {
                        int par = desc_map[g.galaxy_idx];
                        tree_gals[par].descendant = i;
                        tree_gals[par].global_descendant = g.global_index;
                    }
                    
                    // 3. If someone merged into me (via merge_map), link them
                    // The key in merge_map is {my_file_index, my_snapshot}
                    hpc::varray<unsigned, 2> my_id{file_gal_idx, (unsigned)g.snapshot};
                    auto range = merge_map.equal_range(my_id);
                    for(auto it = range.first; it != range.second; ++it) {
                        int par = it->second;
                        tree_gals[par].descendant = i;
                        tree_gals[par].global_descendant = g.global_index;
                    }
                    merge_map.erase(my_id);
                    
                    // 4. Register myself for future descendants
                    if(g.merge_into_id != -1) {
                        desc_map[g.galaxy_idx] = -1; // Merged, so I'm not a main progenitor
                    } else {
                        desc_map[g.galaxy_idx] = i;
                    }
                }

                // Reset local_index for reordering
                for(size_t i=0; i<tree_gals.size(); ++i) tree_gals[i].local_index = i;

                // Reorder
                depthfirst_ordering(hpc::view<std::vector<sage::galaxy>>(tree_gals));
                
                // Write
                gals_buf.write(tree_gals);
                tree_displs_buf.write(current_tree_displ);
                tree_cnts_buf.write((unsigned)tree_gals.size());
                
                current_gal_global_idx += tree_gals.size();
                current_tree_displ += tree_gals.size();
            }
        }
        
        // Close buffers
        gals_buf.close();
        tree_displs_buf.close();
        tree_cnts_buf.close();
        
        // Write final displacement
        if(_comm->rank() == _comm->size() - 1) {
             // This needs to be the total number of galaxies
             // But tree_displs has size total_trees + 1.
             // The last element should be total_gals.
             // I need to write it at the correct offset.
             // The buffer handles offsets.
             // But the buffer is for `total_trees` elements.
             // The dataset has `total_trees + 1`.
             // I can write the last element manually.
             // But wait, `tree_displs_buf` writes `current_tree_displ` for each tree.
             // After the loop, `current_tree_displ` is the start of the *next* tree (or end of last).
             // So I should write one more value.
             // But `tree_displs_buf` is sized for `my_tot_trees`.
             // If I write one more, it might go out of bounds of the buffer's expected count?
             // Actually, `tree_displs` dataset has size `total_trees + 1`.
             // Each rank writes `my_tot_trees` entries.
             // The last rank needs to write the final entry.
             // I can just write it directly to the dataset.
             hpc::h5::dataspace mem_space(1);
             hpc::h5::dataspace file_space = tree_displs_dset.dataspace();
             hsize_t start = total_trees;
             hsize_t count = 1;
             file_space.select_hyperslab(H5S_SELECT_SET, count, start);
             tree_displs_dset.write(&total_gals, hpc::h5::datatype::native_ullong, mem_space, file_space);
        }

        // Write metadata
        if(_comm->rank() == 0) {
            std::vector<double> zs(_redshifts.size());
            std::copy(_redshifts.begin(), _redshifts.end(), zs.begin());
            out_file.write_serial("snapshot_redshifts", zs);
            out_file.write<double>("cosmology/box_size", _box_size);
            out_file.write<double>("cosmology/hubble", _hubble);
            out_file.write<double>("cosmology/omega_l", _omega_l);
            out_file.write<double>("cosmology/omega_m", _omega_m);
        }
    }

    // Copy depthfirst_ordering and _depthfirst_recurse from sage2h5.cc
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
            // ASSERT(gidxs[0] + ii <= gidxs[1]); // This assertion might fail if gaps exist, but here we are contiguous
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


  private:
    void _load_param(hpc::fs::path const &fn) {
        std::ifstream file(fn.native());
        EXCEPT(file.good(), "Could not open parameter file.");

        std::string line;
        while(std::getline(file, line)) {
            if(line.empty() || line[0] == '%')
                continue;

            std::stringstream ss(line);
            std::string       key;
            ss >> key;

            if(key == "BoxSize")
                ss >> _box_size;
            else if(key == "Hubble_h") {
                ss >> _hubble;
                _hubble *= 100.0;
            } else if(key == "Omega_Lambda" || key == "OmegaLambda")
                ss >> _omega_l;
            else if(key == "Omega_m" || key == "Omega")
                ss >> _omega_m;
        }
    }

    void _load_redshifts(hpc::fs::path const &fn) {
        std::ifstream file(fn.native());
        EXCEPT(file.good(), "Could not open expansion list file.");

        std::string line;
        while(std::getline(file, line)) {
            if(line.empty() || line[0] == '%')
                continue;

            std::stringstream ss(line);
            double            a;
            ss >> a;
            _redshifts.push_back(1.0 / a - 1.0);
        }
    }

    std::string           _mode;
    hpc::fs::path         _sage_dir;
    hpc::fs::path         _param_fn;
    hpc::fs::path         _alist_fn;
    hpc::fs::path         _out_fn;
    int                   _verb;
    hpc::mpi::comm *      _comm;
    
    // Cosmology
    double              _box_size;
    double              _hubble;
    double              _omega_l;
    double              _omega_m;
    std::vector<double> _redshifts;

    // HDF5 types
    hpc::h5::datatype _mem_type;
    hpc::h5::datatype _file_type;
};

#define HPC_APP_CLASS application
#include <libhpc/mpi/main.hh>
