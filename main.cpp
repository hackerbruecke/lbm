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

    lbm::io::Config cfg;
	lbm::int3D cart_rank;
    lbm::int3D length  {0,0,0};
    lbm::int3D periods {0,0,0};
    lbm::double_array<3> origin  {0,0,0};
    lbm::double_array<3> spacing {0,0,0};
    lbm::BGKCollision<model>* collision = nullptr;
    lbm::Domain_ptr<model> domain;
    lbm::Domain_ptr<model> subdomain;

    try {
        if (rank == 0) {
            // Read configuration file
            cfg = lbm::io::Config(argc, argv);

            if (cfg.iproc()*cfg.jproc()*cfg.kproc() != number_of_ranks)
                throw std::logic_error("iproc*jproc*kproc != number of ranks!");
            // Allocate collision operator
            collision = new lbm::BGKCollision<model>(cfg.tau());
            domain = lbm::io::parse_scenario_file<model>(cfg.scenario_xml(), cfg, *collision);
            length[0] = domain->xlength();
            length[1] = domain->ylength();
            length[2] = domain->zlength();
            // Print configuration file values and domain lengths
            std::cout << cfg << std::endl;
            std::cout << "> Domain lengths (x,y,z): " << domain->xlength()
                      << ", " << domain->ylength()
                      << ", " << domain->zlength() << std::endl;
            // Copy origin and spacing for broadcasting
            origin = { domain->xorigin(), domain->yorigin(), domain->zorigin() };
            spacing = { domain->xspacing(), domain->yspacing(), domain->zspacing() };
        }
        // Broadcast config values values and create cartesian grid
        lbm::mpi::broadcast(cfg, length, origin, spacing, MPI_COMM_WORLD);
        if (rank != 0)
            collision = new lbm::BGKCollision<model>(cfg.tau());

        // Create cartesian grid
        lbm::int3D dimensions { (int)cfg.iproc(), (int)cfg.jproc(), (int)cfg.kproc() };
        MPI_Cart_create(MPI_COMM_WORLD, 3, dimensions.data(), periods.data(), 0, &comm);
        MPI_Cart_coords(comm, rank, 3, cart_rank.data());

        // Create sudomain for each rank
        domain = lbm::mpi::create_subdomain(rank, dimensions, cfg,
                cart_rank, length, domain, comm,
                origin, spacing, collision);
        // Fluid cells on the boundary of a subdomain are transformed into parallel cells
        domain->replace_fluid_by_parallel_boundary();
        domain->set_nonfluid_cells_nullcollide();

        // Filename is contains rank number for each process that is writing vtk files.
        std::ostringstream os;
        os << cfg.output_filename() << '_' << rank ;

        // Create buffer for sending and receiving messages
        size_t d = std::max(domain->xlength(), std::max(domain->ylength(), domain->zlength()));
        d *= d*9;
        double* sendBuffer = new double[d];
        double* recvBuffer = new double[d];

        // Run simulation
    //    double mlups_duration = 0.0;
    //    double walltime = omp_get_wtime();
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
    }
    catch (const std::exception& ex) {
        std::cerr << "An error occured: " << ex.what() << ". Will exit now." << std::endl;
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}
