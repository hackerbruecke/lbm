#include "model.h"

#define D3Q 19

#if D3Q == 19
    using model = lbm::model::d3q19;
#elif D3Q == 15
    using model = lbm::model::d3q15;
#elif D3Q == 27
    using model = lbm::model::d3q27;
#else
    #error "Set lattice model via -DD3Q=15, -DD3Q=19 or -DD3Q=27"
#endif

#include <mpi.h>
#include <algorithm>
#include <iostream>
#include <memory>
#include <vector>

#include "parallel.h"
#include "lbmdefinitions.h"
#include "helper.h"
#include "collision.h"
#include "boundary.h"
#include "cell.h"
#include "domain.h"
#include "io/configuration.h"
#include "io/vtk.h"
#include "io/scenario.h"
#include "mpihelper.h"

int main(int argc, char** argv)
{
	int rank;
	int number_of_ranks;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_ranks);
    MPI_Comm comm;
    // TODO: Make this more beautiful :-)
	std::array<int, 3> cart_rank;
    std::array<size_t, 3> length { 0, 0, 0 };

    lbm::io::Config cfg;
    lbm::Domain_ptr<model> domain;
    lbm::Domain_ptr<model> subdomain;
    lbm::BGKCollision<model>* collision;
    std::array<int, 3> periods { 0, 0, 0 };
    lbm::double_array<3> origin{0,0,0}, spacing{0,0,0};
//    std::vector<int32_t> buffer;
//    int* buffer;

    if (rank == 0) {
        // Read configuration file
        cfg = lbm::io::Config(argc, argv);

        // Allocate collision operator
        collision = new lbm::BGKCollision<model>(cfg.tau());
        domain = lbm::io::parse_scenario_file<model>(cfg.scenario_xml(), cfg, *collision);
        length[0] = domain->xlength();
        length[1] = domain->ylength();
        length[2] = domain->zlength();
        // Print configuration file values and domain lengths
        std::cout << cfg << std::endl;
        std::cout << "> Domain lengths (x,y,z): " << domain->xlength() << ", " << domain->ylength()
                  << ", " << domain->zlength() << std::endl;
        origin = { domain->xorigin(), domain->yorigin(), domain->zorigin() };
        spacing = { domain->xspacing(), domain->yspacing(), domain->zspacing() };
    }
    // Broadcast config values values and create cartesian grid
    lbm::mpi::broadcast(cfg, length, origin, spacing, MPI_COMM_WORLD);
    if (rank != 0)
        collision = new lbm::BGKCollision<model>(cfg.tau());
    //Create cartesian grid
    std::array<int, 3> dimensions { (int)cfg.iproc(), (int)cfg.jproc(), (int)cfg.kproc() };
//    std::cout << rank << " Dim: " << dimensions[0] << " " << dimensions[1]  << " " << dimensions[2] << std::endl;
    MPI_Cart_create(MPI_COMM_WORLD, 3, dimensions.data(), periods.data(), 0, &comm);
    MPI_Cart_coords(comm, rank, 3, cart_rank.data());

    // Compute sublengths to create domains
    auto subl = lbm::mpi::compute_sublengths(length, cart_rank, cfg);
//    size_t bufsize = (subl[0]+2)*(subl[1]+2)*(subl[2]+2);
//    buffer = new int[bufsize];

    if (rank == 0) {
        // TODO: Serialize and send flagfields to respective domains
        for (auto k = 0; k < dimensions[2]; ++k) {
            for (auto j = 0; j < dimensions[1]; ++j) {
                for (auto i = 0; i < dimensions[0]; ++i) {
                    // Weird fancy loop
                    auto send_subl = lbm::mpi::compute_sublengths(length, { i, j, k}, cfg);
                    int index = 0;
                    size_t bufsize = (send_subl[0]+2)*(send_subl[1]+2)*(send_subl[2]+2);
                    int* buffer = new int[bufsize];
                    int xrmd = 0, yrmd = 0, zrmd = 0;
                    if (i >= length[0] % dimensions[0]) {
                    	xrmd = length[0] % dimensions[0];
                    }
                    if (j >= length[1] % dimensions[1]) {
                    	yrmd = length[1] % dimensions[1];
                    }
                    if (k >= length[2] % dimensions[2]) {
                    	zrmd = length[2] % dimensions[2];
                    }
                    for (auto z = 0u; z < send_subl[2] + 2; ++z) {
                        for (auto y = 0u; y < send_subl[1] + 2; ++y) {
                            for (auto x = 0u; x < send_subl[0] + 2; ++x) {
//                                buffer.push_back(domain->cell(x, y, z).get_type());
                                buffer[index++] = static_cast<int>(
                                        *domain->cell(
                                                x+i*(send_subl[0])+xrmd,
                                                y+j*(send_subl[1])+yrmd,
                                                z+k*(send_subl[2])+zrmd)
                                        .get_collision_handler());
                            }
                        }
                    }
                    int dest;
                    std::array<int, 3> ijk {i, j, k};
                    MPI_Cart_rank(comm, ijk.data(), &dest);
//                    bufsize = (subl_tmp[0]+2)*(subl_tmp[1]+2)*(subl_tmp[2]+2);
                    if (dest != 0) {
                        // TODO: Also compute displacment of origin and send it!
//                        std::cout << rank << " sending " << bufsize
//                                << " bytes of data to " << dest << std::endl;
                        MPI_Send(buffer, bufsize, MPI_INT, dest, 0, comm);
                    }
                    else {
                        subdomain = lbm::mpi::create_subdomain_from_buffer(send_subl, origin, spacing, buffer, bufsize, *collision);
//                        std::cout << "0 created subdomain with lengths: "
//                                << subdomain->xlength() << ' ' << subdomain->ylength()
//                                << ' ' << subdomain->zlength()<< std::endl;
                    }
                    delete [] buffer;
                }
            }
        }
    }
    else {
        // Compute sublengths to create domains
//        auto subl = lbm::mpi::compute_sublengths(length, cart_rank, cfg);
        size_t bufsize = (subl[0]+2)*(subl[1]+2)*(subl[2]+2);
        int* buffer = new int[bufsize];
        // Receive flag fields from rank 0 and create domain
        MPI_Recv(buffer, bufsize, MPI_INT, 0, 0, comm, MPI_STATUS_IGNORE);
//        std::cout << rank << "Received " << bufsize <<
//                ", creating with (" << subl[0] << ' ' << subl[1] << ' ' << subl[2] <<  ")" << std::endl;
        for (auto i=0u; i<origin.size(); ++i) {
        	int rmd = 0;
            if (cart_rank[i] >= length[i] % dimensions[i]) {
            	rmd = length[i] % dimensions[i];
            }
            origin[i] += (cart_rank[i]*subl[i]+rmd)*spacing[i];
        }
        subdomain = lbm::mpi::create_subdomain_from_buffer(subl, origin, spacing, buffer, bufsize, *collision);
//        std::cout << rank << " created subdomain with lengths: "
//            << subdomain->xlength() << ' ' << subdomain->ylength()
//            << ' ' << subdomain->zlength()<< std::endl;
        delete [] buffer;
    }
    subdomain->replace_fluid_by_parallel_boundary();
    domain = std::move(subdomain);
//    if (cart_rank[0] == 0 && cart_rank[1] == 0 && cart_rank[2] == 1) {
//        for (auto z = 0u; z < domain->zlength()+2; ++z) {
//            for (auto y = 0u; y < domain->ylength()+2; ++y) {
////                cout << '-' << static_cast<size_t>(*domain->cell(domain->xlength()+1, y, z).get_collision_handler());
//            }
//        }
//    }

//    std::cout << domain->xlength() << " " << domain->ylength() << " " << domain->zlength() << std::endl;
    domain->set_nonfluid_cells_nullcollide();
#if 1

    std::ostringstream os;
    os << cfg.output_filename() << '_' << rank ;
    // Run simulation
//    double mlups_duration = 0.0;
//    double walltime = omp_get_wtime();

    size_t d = *std::max_element(subl.begin(), subl.end()) + 2;
    d *= d*9;
//    size_t d = std::max(subl.begin(), subl.end()) + 2;
    double* sendBuffer = new double[d];
    double* recvBuffer = new double[d];

    for (auto t = 1u; t <= cfg.timesteps(); ++t) {
//        double start = omp_get_wtime();
        lbm::mpi::exchangePdfs<model>(sendBuffer, recvBuffer, *domain, comm);
        domain->stream();
        domain->swap();
        domain->collide();
//        mlups_duration += omp_get_wtime() - start;

        // Output file only when timesteps-per-plot was set to != 0
        if (cfg.timesteps_per_plot() && t % cfg.timesteps_per_plot() == 0)
            lbm::io::write_vtk_file(*domain, cfg.output_dir(), os.str(), t);
        // Print percents completed
        if (rank == 0) {
            std::cout << "\r" << (int) ((double) t / cfg.timesteps() * 100) << " %";
            std::cout.flush();
        }
    }
    if (rank == 0) {
        std::cout << "\nFinished!" << std::endl;
//        std::cout << "Total runtime: " << omp_get_wtime()-walltime << " seconds." << std::endl;
//        double mlups = (2+domain->xlength())*(2+domain->ylength())*(2+domain->zlength())
//        *cfg.timesteps()/(mlups_duration*1e6);
//        std::cout << "MLUPS: " << mlups << std::endl;
    }
    delete [] sendBuffer;
    delete [] recvBuffer;
#endif
//    delete [] buffer;
    MPI_Finalize();
    return EXIT_SUCCESS;
}
