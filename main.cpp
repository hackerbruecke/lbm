#include "model.h"
#define D3Q 19

#if D3Q == 19
    using model = lbm::model::d3q19;
#elif D3Q == 15
    using model = lbm::model::d3q15;
#elif D3Q == 27
    using model = lbm::model::d3q27;
#else
    #error "Set lattice model via -DD3Q15, -DD3Q19 or -DD3Q27"
#endif

#include <mpi.h>
#include <iostream>
#include <memory>

#include "parallel.h"
#include "configuration.h"
#include "lbmdefinitions.h"
#include "helper.h"
#include "collision.h"
#include "boundary.h"
#include "cell.h"
#include "domain.h"
#include "visual.h"

void setBoundaryConditions(lbm::Domain<model>& domain, lbm::NoSlipBoundary<model>& no_slip,
        lbm::InflowBoundary<model>& inflow, lbm::OutflowBoundary<model>& outflow)
{
    auto xl = domain.xlength();
    auto yl = domain.ylength();
    auto zl = domain.zlength();
    // Boundary conditions to be applied on lattice cells.
    // Bottom XY plane
    domain.setBoundaryCondition(no_slip, 0, xl + 1, 0, yl + 1, 0, 0);
    // Top XY plane
    domain.setBoundaryCondition(no_slip, 0, xl + 1, 0, yl + 1, zl + 1, zl + 1);

    // Front XZ plane
    domain.setBoundaryCondition(no_slip, 0, xl + 1, 0, 0, 0, zl);
    // Rear XZ plane
    domain.setBoundaryCondition(no_slip, 0, xl + 1, yl + 1, yl + 1, 0, zl);

    // Left YZ plane
    domain.setBoundaryCondition(inflow, 0, 0, 0, yl + 1, 0, zl);
    // Right YZ plane
    domain.setBoundaryCondition(outflow, xl + 1, xl + 1, 0, yl + 1, 0, zl);
}


int main(int argc, char** argv)
{
    try {
        // Read configuration file
        lbm::io::Config cfg(argc, argv);
        // Allocate collision operator
        auto collision = lbm::BGKCollision<model>(cfg.tau());
        // Create domain from VTK point file
        auto domain = lbm::io::read_vtk_point_file<model,
                lbm::NoSlipBoundary<model>>(cfg.input_vtk(), collision);
        // Create boundary conditions and apply them
        lbm::double_array<model::D> inflowVelocity = { 0.003496, 0.0, 0.0 };
        auto no_slip = lbm::NoSlipBoundary<model>(*domain);
        auto outflow = lbm::OutflowBoundary<model>(*domain);
        auto inflow =  lbm::InflowBoundary<model>(*domain, inflowVelocity);
        setBoundaryConditions(*domain, no_slip, inflow, outflow);

        // Print configuration file values and domain lengths
        std::cout << cfg << std::endl;
        std::cout << "Domain lengths (x,y,z): " << domain->xlength() << ", " << domain->ylength()
                  << ", " << domain->zlength() << std::endl;
        std::cout << "Starting simulation..." << std::endl;

        double duration = 0.0;
        // Run simulation
        for (auto t = 1u; t <= cfg.timesteps(); ++t) {
            double start = omp_get_wtime();
            domain->stream();
            domain->swap();
            domain->collide();
            duration += omp_get_wtime() - start;

            if (t % cfg.timesteps_per_plot() == 0)
                lbm::io::write_vtk_file(*domain, cfg.output_dir(), cfg.output_filename(), t);
            // Print percents completed
            std::cout << "\r" << (int) ((double) t / cfg.timesteps() * 100) << " %";
            std::cout.flush();
        }
        std::cout << "\nFinished!" << std::endl;
        double mlups = domain->xlength()*domain->ylength()*domain->zlength()
                *cfg.timesteps()/(duration*1e6);
        std::cout << "MLUPS: " << mlups << std::endl;
    }
    catch (const std::exception& ex) {
        std::cerr << "An error occured: " << ex.what() << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
