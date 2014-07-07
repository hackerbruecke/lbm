#pragma once

#include "collision.h"
#include "io/configuration.h"
#include "io/scenario.h"
#include "lbmdefinitions.h"

#include <mpi.h>

namespace lbm
{
namespace mpi
{

void init(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
}

template <typename lattice_model>
class sim
{
    int rank;
    int number_of_ranks;
    MPI_Comm comm;
    std::array<int, 3> cart_rank;
    std::array<int, 3> periods { 0, 0, 0 };
    std::array<size_t, 3> length { 0, 0, 0 };
    lbm::io::Config cfg;
    lbm::Domain_ptr<lattice_model> domain;
    lbm::Domain_ptr<lattice_model> subdomain;
    lbm::BGKCollision<lattice_model>* collision;
    lbm::double_array<3> origin{0,0,0};
    lbm::double_array<3> spacing{0,0,0};
    std::array<int, 3> dimensions;
    double* sendBuffer;
    double* recvBuffer;


private:

    sim(int argc, char** argv)
    : periods{0,0,0},
      length{0,0,0},
      origin{0,0,0},
      spacing{0,0,0}
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &number_of_ranks);
        configure(argc, argv);
        create_cartesian_grid();
    }
public:
    ~sim()
    {
        delete [] sendBuffer;
        delete [] recvBuffer;
    }
    static sim* inst;
public:
    io::Config& config() { return cfg; }
    int get_rank() const { return rank; }

    static sim& instance(int argc, char** argv) {
        if (inst == nullptr) {
            inst = new sim(argc, argv);
        }
        return *inst;
    }
private:
    void configure(int argc, char** argv)
    {
        if (rank == 0) {
              // Read configuration file
              cfg = lbm::io::Config(argc, argv);

              // Allocate collision operator
              collision = new lbm::BGKCollision<lattice_model>(cfg.tau());
              domain = lbm::io::parse_scenario_file<lattice_model>(cfg.scenario_xml(), cfg, *collision);
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
              collision = new lbm::BGKCollision<lattice_model>(cfg.tau());
    }

    void create_cartesian_grid()
    {
        dimensions = { (int)cfg.iproc(), (int)cfg.jproc(), (int)cfg.kproc() };
        MPI_Cart_create(MPI_COMM_WORLD, 3, dimensions.data(), periods.data(), 0, &comm);
        MPI_Cart_coords(comm, rank, 3, cart_rank.data());
    }
public:
    Domain_ptr<lattice_model> create_subdomain()
    {
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

                        for (auto z = 0u; z < send_subl[2] + 2; ++z) {
                            for (auto y = 0u; y < send_subl[1] + 2; ++y) {
                                for (auto x = 0u; x < send_subl[0] + 2; ++x) {
    //                                buffer.push_back(domain->cell(x, y, z).get_type());
                                    buffer[index++] = static_cast<int>(
                                            *domain->cell(
                                                    x+i*(send_subl[0]),
                                                    y+j*(send_subl[1]),
                                                    z+k*(send_subl[2]))
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
                origin[i] += cart_rank[i]*subl[i]*spacing[i];
            }
            subdomain = lbm::mpi::create_subdomain_from_buffer(subl, origin, spacing, buffer, bufsize, *collision);
    //        std::cout << rank << " created subdomain with lengths: "
    //            << subdomain->xlength() << ' ' << subdomain->ylength()
    //            << ' ' << subdomain->zlength()<< std::endl;
            delete [] buffer;
        }
        return std::move(subdomain);
    }

    void prepare_buffers()
    {
        size_t d = std::max(domain->xlength(), std::max(domain->ylength(), domain->zlength()));
        d *= d*9;
        sendBuffer = new double[d];
        recvBuffer = new double[d];
    }
private:
    void exchange(int dir, int displ)
    {
        auto& flagField = *this->domain;
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
                     recvBuffer, nRecv, MPI_DOUBLE, source, 0, comm, MPI_STATUS_IGNORE);

        const auto xl = flagField.xlength()+2;
        const auto yl = flagField.ylength()+2;

        for (int k = 0; k < nRecv; ++k) {
            const auto idx = recvInd[k];
            const auto x = idx/lattice_model::Q % xl;
            const auto y = idx/(lattice_model::Q * xl) % yl;
            const auto z = idx/(lattice_model::Q * xl * yl);
            flagField.cell(x, y, z)[idx % lattice_model::Q] = recvBuffer[k];
        }
        delete [] recvInd;
    }
public:
    void exchangePdfs()
    {
        for (int dir = 0; dir < 3; ++dir) {
            exchange(dir, -1);
            exchange(dir, 1);
        }
    }
};

template <typename lattice_model>
sim<lattice_model>* sim<lattice_model>::inst = nullptr;

}//namespace mpi
}//namespace lbm
