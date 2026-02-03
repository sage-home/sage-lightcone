// Copyright 2012 Luke Hodkinson

// This file is part of libhpc.
//
// libhpc is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// libhpc is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with libhpc.  If not, see <http://www.gnu.org/licenses/>.

#ifndef hpc_mpi_main_hh
#define hpc_mpi_main_hh

#ifndef HPC_APP_CLASS
#define HPC_APP_CLASS application
#endif

#include <libhpc/debug/assertions.hh>
#include <libhpc/mpi/comm.hh>
#include <libhpc/mpi/init.hh>

namespace hpc {

extern hpc::application *global_app;
}

int main(int argc, char *argv[]) {
  typedef HPC_APP_CLASS application_type;
  {
    int nchunks = 1;
    std::string option_chunks("--chunks=");
    std::vector<std::string> pass_argv;
    int i;
    for (i = 0; i < argc; i++) {
      std::string argv_string = argv[i];
      if (argv_string.compare(0, option_chunks.length(), option_chunks) == 0) {
        // std::cout << "found chunks option"<<std::endl;
        nchunks = std::atoi(argv_string.c_str() + option_chunks.length());
      } else {
        pass_argv.push_back(argv_string);
      }
    }
    char **argv_pass = new char *[pass_argv.size()];
    for (i = 0; i < (int)pass_argv.size(); i++) {
      argv_pass[i] = (char *)pass_argv[i].c_str();
    }
    // std::cout << ".Number of chunks "<<nchunks<<std::endl;
    hpc::mpi::initialise(argc, argv);
    // std::cout << "after initialise"<<std::endl;
    for (int chunk = 0; chunk < nchunks; chunk++) {
      // std::cout << " start_chunk"<<std::endl;
      // bool needs_to_be_done = hpc::mpi::start_chunk(chunk, nchunks);
      bool needs_to_be_done = true;
      if (needs_to_be_done) {
        try {
          // std::cout << " start_chunk create app "<<
          // pass_argv.size()<<std::endl;
          application_type app(pass_argv.size(), argv_pass);
          hpc::global_app = &app;

          // std::cout << " do start_chunk"<<std::endl;
          app();

          // std::cout << " done start_chunk"<<std::endl;
        } catch (hpc::silent_terminate &ex) {
          continue;
          // std::cout << " ?hpc::silent_terminate? "<<std::endl;
        }
#ifdef NDEBUG
        catch (std::exception &ex) {
          std::cerr << "\nError: " << ex.what() << "\n\n";
          hpc::mpi::comm::world.abort();
        }
#endif
        // std::cout << " finishing "<<std::endl;
        // hpc::mpi::finish_chunk(chunk, nchunks);
      } else {
        std::cout << "SKIP rank="
                  << hpc::mpi::comm::world.rank() * nchunks + chunk
                  << std::endl;
      }
    }
    hpc::mpi::comm::world.barrier();
    hpc::mpi::finalise();
  }
  return EXIT_SUCCESS;
}

#endif
