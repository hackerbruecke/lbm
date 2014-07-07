#pragma once

#include <mpi.h>
#include <boost/mpi.hpp>
#include "io/configuration.h"
#include "lbmdefinitions.h"

namespace lbm
{
namespace mpi
{

void broadcast(lbm::io::Config& cfg, std::array<size_t, 3>& domain_length,
        lbm::double_array<3>& origin, lbm::double_array<3>& spacing,
        MPI_Comm comm)
{
    boost::mpi::communicator c(comm, boost::mpi::comm_attach);
	// Broadcast strings
    boost::mpi::broadcast(c, cfg._collision_model, 0);
    boost::mpi::broadcast(c, cfg._input_file, 0);
    boost::mpi::broadcast(c, cfg._output_dir, 0);
    boost::mpi::broadcast(c, cfg._output_filename, 0);
    boost::mpi::broadcast(c, cfg._scenario_xml, 0);

    // Broadcast numerical values
    boost::mpi::broadcast(c, cfg._timesteps, 0);
    boost::mpi::broadcast(c, cfg._timesteps_per_plot, 0);
    boost::mpi::broadcast(c, cfg._tau, 0);
    boost::mpi::broadcast(c, cfg._omp_threads, 0);
    // TODO: Invoke this automatically
    omp_set_num_threads(cfg._omp_threads);

    boost::mpi::broadcast(c, cfg._iproc, 0);
    boost::mpi::broadcast(c, cfg._jproc, 0);
    boost::mpi::broadcast(c, cfg._kproc, 0);

    // broadcast domain lengths
//    boost::mpi::broadcast(c, domain_length, 0);
    boost::mpi::broadcast(c, domain_length[0], 0);
    boost::mpi::broadcast(c, domain_length[1], 0);
    boost::mpi::broadcast(c, domain_length[2], 0);

    boost::mpi::broadcast(c, origin[0], 0);
    boost::mpi::broadcast(c, origin[1], 0);
    boost::mpi::broadcast(c, origin[2], 0);

    boost::mpi::broadcast(c, spacing[0], 0);
    boost::mpi::broadcast(c, spacing[1], 0);
    boost::mpi::broadcast(c, spacing[2], 0);
}

inline auto compute_sublengths(
		const std::array<size_t, 3>& length,
		const std::array<int, 3>& cart_rank,
		const lbm::io::Config& cfg) -> std::array<size_t, 3>
{
	std::array<size_t, 3> subl = {
			length[0] / cfg.iproc(),
			length[1] / cfg.jproc(),
			length[2] / cfg.kproc()
	};
	if (cart_rank[0] < (int)length[0] % (int)cfg.iproc()) ++subl[0];
	if (cart_rank[1] < (int)length[1] % (int)cfg.jproc()) ++subl[1];
	if (cart_rank[2] < (int)length[2] % (int)cfg.kproc()) ++subl[2];
	return subl;
}

template <typename lattice_model>
auto get_collision_from_type(Domain<lattice_model>& domain, size_t type) -> Collision<lattice_model>&
{
    static BGKCollision<lattice_model> bgk(0.6);//tau
    static NoSlipBoundary<lattice_model> noslip(domain);
    static MovingWallBoundary<lattice_model> movingwall(domain, {0.05, 0, 0});//Wall velocity
    static FreeSlipBoundary<lattice_model> freeslip(domain);
    static OutflowBoundary<lattice_model> outflow(domain);
    static InflowBoundary<lattice_model> inflow(domain, {0.0, 0.0, 0.05});//inflow velocity
    static PressureBoundary<lattice_model> pressure(domain, 1.005);//input density
    static parallel::ParallelBoundary<lattice_model> parallel(domain);

    if (type == bgk){
        return bgk;
    }
    else if (type == noslip) {
        return noslip;
    }
    else if (type == movingwall) {
        return movingwall;
    }
    else if (type == freeslip) {
        return freeslip;
    }
    else if (type == outflow) {
        return outflow;
    }
    else if (type == inflow) {
        return inflow;
    }
    else if (type == pressure) {
        return pressure;
    }
    else if (type == parallel) {
        return parallel;
    }
    else {
        throw std::exception();
    }
}

template <typename lattice_model>
auto create_subdomain_from_buffer(
		const uint_array<3>& subl,
		const double_array<3>& origin,
		const double_array<3>& spacing,
		const int* const buffer,
		size_t bufsize,
		FluidCollision<lattice_model>& collision) -> Domain_ptr<lattice_model>
{
	auto domain = make_unique<Domain<lattice_model>>(
	        subl[0], subl[1], subl[2], collision,
	        origin[0], origin[1], origin[2],
	        spacing[0], spacing[1], spacing[2]);
//	std::cout << "Collide field size: " << domain->collide_field.size()
//	        << " - buffer size " << bufsize<< std::endl;
	assert(domain->collide_field.size() == bufsize);

	for (auto i = 0u; i < bufsize; ++i) {
	    auto& collision = get_collision_from_type<lattice_model>(*domain,
	            static_cast<size_t>(buffer[i]));
		domain->collide_field[i].set_collision_handler(&collision);
	}
	return std::move(domain);
}

template <typename lattice_model>
void exchange(double *sendBuffer, double *readBuffer, Domain<lattice_model>& flagField,
        MPI_Comm comm, int dir, int displ)
{
    int source, dest;
    int x0[3], x1[3], nSend=0, nRecv=0;
    int Qi[9], Qn = 0;
    const std::array<int, 3> sublength {
        (int)flagField.xlength(),
        (int)flagField.ylength(),
        (int)flagField.zlength()
    };

    for (int i = 0; i < (int)lattice_model::Q; ++i) {
        if (lattice_model::velocities[i][dir] == displ) {
            Qi[Qn++] = i;
        }
    }
    for (int i = 0; i < 3; ++i) {
        x0[i] = 0;
        x1[i] = sublength[i] + 2;
    }
    x0[dir] = (displ + 1) / 2 * (x1[dir] - 3) + 1;
    x1[dir] = x0[dir] + 1;
    MPI_Cart_shift(comm, dir, displ, &source, &dest);

    if (dest != MPI_PROC_NULL) {
        for (int z = x0[2]; z < x1[2]; ++z) {
            for (int y = x0[1]; y < x1[1]; ++y) {
                for (int x = x0[0]; x < x1[0]; ++x) {
                    auto flag = flagField.cell(x, y, z).get_collision_handler()->get_type();
                    if (flag == CollisionType::bgk ||
                        flag == CollisionType::parallel ||
                        flag == CollisionType::freeslip) {
                        for (int i = 0; i < Qn; ++i) {
                            sendBuffer[nSend++] = flagField.cell(x, y, z)[Qi[i]];
                        }
                    }
                }
            }
        }
    }
    int* recvInd = new int[9 * (sublength[(dir+1)%3] + 2) * (sublength[(dir+2)%3] + 2)];
    if (source != MPI_PROC_NULL) {
        x0[dir] = x0[dir] - sublength[dir] * displ;
        x1[dir] = x0[dir] + 1;
        for (int z = x0[2]; z < x1[2]; ++z) {
            for (int y = x0[1]; y < x1[1]; ++y) {
                for (int x = x0[0]; x < x1[0]; ++x) {
                    auto flag = flagField.cell(x, y, z).get_collision_handler()->get_type();
                    if (flag == CollisionType::parallel || flag == CollisionType::freeslip) {
                        for (int i = 0; i < Qn; ++i) {
                            recvInd[nRecv++] = flagField.idx(x, y, z)*lattice_model::Q + Qi[i];
                        }
                    }
                }
            }
        }
    }
    MPI_Sendrecv(sendBuffer, nSend, MPI_DOUBLE, dest, 0,
                 readBuffer, nRecv, MPI_DOUBLE, source, 0, comm, MPI_STATUS_IGNORE);

    const auto xl = flagField.xlength()+2;
    const auto yl = flagField.ylength()+2;

    for (int k = 0; k < nRecv; ++k) {
        const auto idx = recvInd[k];
        const auto x = idx/lattice_model::Q % xl;
        const auto y = idx/(lattice_model::Q * xl) % yl;
        const auto z = idx/(lattice_model::Q * xl * yl);
        flagField.cell(x, y, z)[idx % lattice_model::Q] = readBuffer[k];
    }
    delete [] recvInd;
}

template <typename lattice_model>
void exchangePdfs(double* sendBuffer, double* readBuffer, Domain<lattice_model>& flagField, MPI_Comm& comm)
{
    for (int dir = 0; dir < 3; ++dir) {
        exchange(sendBuffer, readBuffer, flagField, comm, dir, -1);
        exchange(sendBuffer, readBuffer, flagField, comm, dir, 1);
    }
}


}//namespace mpi
}//namespace lbm

