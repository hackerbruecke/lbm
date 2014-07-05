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
#include <iostream>
#include <memory>

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

int main(int argc, char** argv)
{
    try {
        // Read configuration file
        lbm::io::Config cfg(argc, argv);
        // Allocate collision operator
        auto collision = lbm::BGKCollision<model>(cfg.tau());
        auto domain = lbm::io::parse_scenario_file<model>(cfg.scenario_xml(), cfg, collision);
        // Print configuration file values and domain lengths
        std::cout << cfg << std::endl;
        std::cout << "> Domain lengths (x,y,z): " << domain->xlength() << ", " << domain->ylength()
                  << ", " << domain->zlength() << std::endl;
        std::cout << "Starting simulation..." << std::endl;

        domain->set_nonfluid_cells_nullcollide();
        // Run simulation
        double mlups_duration = 0.0;
        double walltime = omp_get_wtime();
        for (auto t = 1u; t <= cfg.timesteps(); ++t) {
            double start = omp_get_wtime();
            domain->stream();
            domain->swap();
            domain->collide();
            mlups_duration += omp_get_wtime() - start;

            // Output file only when timesteps-per-plot was set to != 0
            if (cfg.timesteps_per_plot() && t % cfg.timesteps_per_plot() == 0)
                lbm::io::write_vtk_file(*domain, cfg.output_dir(), cfg.output_filename(), t);
            // Print percents completed
            std::cout << "\r" << (int) ((double) t / cfg.timesteps() * 100) << " %";
            std::cout.flush();
        }
        std::cout << "\nFinished!" << std::endl;
        std::cout << "Total runtime: " << omp_get_wtime()-walltime << " seconds." << std::endl;
        double mlups = (2+domain->xlength())*(2+domain->ylength())*(2+domain->zlength())
                *cfg.timesteps()/(mlups_duration*1e6);
        std::cout << "MLUPS: " << mlups << std::endl;
    }
    catch (const std::exception& ex) {
        std::cerr << "An error occured: " << ex.what() << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
