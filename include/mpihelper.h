#pragma once

#include <mpi.h>
#include <boost/mpi.hpp>
#include "io/configuration.h"
#include "lbmdefinitions.h"

namespace lbm
{
namespace mpi
{

template<typename lattice_model>
class ParallelBoundary : public NonFluidCollision<lattice_model>
{
public:
    ParallelBoundary(Domain<lattice_model>& domain)
        : NonFluidCollision<lattice_model>(CollisionType::parallel, domain)
    {}

    void collide(Cell<lattice_model>& cell,
            const uint_array<lattice_model::D>& position) const override
    {}
};

void broadcast(lbm::io::Config& cfg, int3D& domain_length,
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
		const int3D& length,
		const int3D& cart_rank,
		const lbm::io::Config& cfg) -> int3D
{
	int3D subl = {
			length[0] / (int)cfg.iproc(),
			length[1] / (int)cfg.jproc(),
			length[2] / (int)cfg.kproc()
	};
	// If the dimprocs are not a divisor of the rescpective length,
	// distribute the remainders evenly to the sublengths
	if (cart_rank[0] < (int)length[0] % (int)cfg.iproc()) ++subl[0];
	if (cart_rank[1] < (int)length[1] % (int)cfg.jproc()) ++subl[1];
	if (cart_rank[2] < (int)length[2] % (int)cfg.kproc()) ++subl[2];
	return subl;
}

template <typename lattice_model>
auto get_collision_from_type(Domain<lattice_model>& domain, const io::Config& cfg,
        size_t type) -> Collision<lattice_model>&
{
    // TODO: Fix those constant values!!!
    static BGKCollision<lattice_model> bgk(cfg.tau());//tau
    static NoSlipBoundary<lattice_model> noslip(domain);
    static MovingWallBoundary<lattice_model> movingwall(domain, {0.05, 0, 0});//Wall velocity
    static FreeSlipBoundary<lattice_model> freeslip(domain);
    static OutflowBoundary<lattice_model> outflow(domain);
    static InflowBoundary<lattice_model> inflow(domain, {0.003496177, 0.0, 0.0});//inflow velocity
    static PressureBoundary<lattice_model> pressure(domain, 1.005);//input density
    static ParallelBoundary<lattice_model> parallel(domain);

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
		const int3D& subl,
		const double_array<3>& origin,
		const double_array<3>& spacing,
		const int* const buffer,
		size_t bufsize,
		FluidCollision<lattice_model>& collision,
		const io::Config& cfg) -> Domain_ptr<lattice_model>
{
	auto domain = make_unique<Domain<lattice_model>>(
	        subl[0], subl[1], subl[2], collision,
	        origin[0], origin[1], origin[2],
	        spacing[0], spacing[1], spacing[2]);

	assert(domain->collide_field.size() == bufsize);
	// Copy whole buffer into the domain's collide field
	for (auto i = 0u; i < bufsize; ++i) {
	    auto& collision = get_collision_from_type<lattice_model>(*domain, cfg,
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
    const int3D sublength {
        (int)flagField.xlength(),
        (int)flagField.ylength(),
        (int)flagField.zlength()
    };
    // Copy relevant pdfs into array
    for (int i = 0; i < (int)lattice_model::Q; ++i) {
        if (lattice_model::velocities[i][dir] == displ) {
            Qi[Qn++] = i;
        }
    }
    // Compute first and last plane of current direction
    for (int i = 0; i < 3; ++i) {
        x0[i] = 0;
        x1[i] = sublength[i] + 2;
    }
    x0[dir] = (displ + 1) / 2 * (x1[dir] - 3) + 1;
    x1[dir] = x0[dir] + 1;
    MPI_Cart_shift(comm, dir, displ, &source, &dest);

    // Put relevant pdfs into the send buffer (no-slip etc. are ignored)
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
    // Put indices of receiving pdfs into an index buffer (we skip cells that won't receive anything)
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
    // We must put the received cells into their proper position, according to the index buffer
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
    // Exchange data in all three directions
    for (int dir = 0; dir < 3; ++dir) {
        exchange(sendBuffer, readBuffer, flagField, comm, dir, -1);
        exchange(sendBuffer, readBuffer, flagField, comm, dir, 1);
    }
}

template <typename lattice_model>
Domain_ptr<lattice_model> create_subdomain(int rank, const int3D& dimensions,
        const io::Config& cfg, const int3D& cart_rank,
        const int3D& length, Domain_ptr<lattice_model>& domain, MPI_Comm comm,
        lbm::double_array<3>& origin, const lbm::double_array<3>& spacing,
        FluidCollision<lattice_model>* collision)
{

    Domain_ptr<lattice_model> subdomain;
    if (rank == 0) {
        // Serialize and send flagfields to respective domains
        for (auto k = 0; k < dimensions[2]; ++k) {
            for (auto j = 0; j < dimensions[1]; ++j) {
                for (auto i = 0; i < dimensions[0]; ++i) {
                    // Compute dimensions of destination subdomain
                    auto send_subl = lbm::mpi::compute_sublengths(length, { i, j, k}, cfg);
                    int index = 0;
                    // Compute buffer size for dimensions of the subdomain that has top be send.
                    size_t bufsize = (send_subl[0]+2)*(send_subl[1]+2)*(send_subl[2]+2);
                    int* buffer = new int[bufsize];
                    // If dimprocs % length != 0 we need to treat the remainders
                    int xrmd = 0, yrmd = 0, zrmd = 0;
                    if (i >= (int)length[0] % dimensions[0]) {
                        xrmd = length[0] % dimensions[0];
                    }
                    if (j >= (int)length[1] % dimensions[1]) {
                        yrmd = length[1] % dimensions[1];
                    }
                    if (k >= (int)length[2] % dimensions[2]) {
                        zrmd = length[2] % dimensions[2];
                    }
                    // For the whole subdomain, copy respective parts of the original domain
                    // into the buffer.
                    for (auto z = 0; z < send_subl[2] + 2; ++z) {
                        for (auto y = 0; y < send_subl[1] + 2; ++y) {
                            for (auto x = 0; x < send_subl[0] + 2; ++x) {
                                buffer[index++] = static_cast<int>(
                                        *domain->cell(
                                                x+i*(send_subl[0])+xrmd,
                                                y+j*(send_subl[1])+yrmd,
                                                z+k*(send_subl[2])+zrmd)
                                        .get_collision_handler());
                            }
                        }
                    }
                    // Compute destination rank which shall receive buffer data.
                    int dest;
                    int3D ijk {i, j, k};
                    MPI_Cart_rank(comm, ijk.data(), &dest);
                    // Either send data for the respective process to create the domain
                    if (dest != 0)
                        MPI_Send(buffer, bufsize, MPI_INT, dest, 0, comm);
                    else    // or directly create the domain for the rank 0 process
                        subdomain = lbm::mpi::create_subdomain_from_buffer(
                                send_subl, origin, spacing, buffer, bufsize, *collision, cfg);
                    delete [] buffer;
                }
            }
        }
    }
    else {
        // Compute sublengths to create domains
        auto subl = lbm::mpi::compute_sublengths(length, cart_rank, cfg);
        // Compute sublengths to create domains
        size_t bufsize = (subl[0]+2)*(subl[1]+2)*(subl[2]+2);
        int* buffer = new int[bufsize];
        // Receive flag fields from rank 0 and create domain
        MPI_Recv(buffer, bufsize, MPI_INT, 0, 0, comm, MPI_STATUS_IGNORE);
        // If there are any remainders of dimproc % length, they must be treated here.
        for (auto i=0u; i<origin.size(); ++i) {
            int rmd = 0;
            if (cart_rank[i] >= (int)length[i] % dimensions[i]) {
                rmd = length[i] % dimensions[i];
            }
            origin[i] += (cart_rank[i]*subl[i]+rmd)*spacing[i];
        }
        // Create subdomain from the data received.
        subdomain = lbm::mpi::create_subdomain_from_buffer(subl, origin, spacing,
                buffer, bufsize,*collision, cfg);
        delete [] buffer;
    }
    return std::move(subdomain);
}


}//namespace mpi
}//namespace lbm

