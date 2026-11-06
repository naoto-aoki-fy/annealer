#include <annealer/replica_exchange.hpp>
#include "demo_models.hpp"
#include <iomanip>
#include <iostream>

// No input files: every rank constructs the same existing nonlinear demo model.
// This example doubles as a small benchmark; all normal output belongs to rank 0.
int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank = 0, size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    try {
        // Coordinate failures outside the solver too, before entering its collectives.
        std::optional<annealer::PreparedModel> prepared;
        int failed = 0, any_failed = 0;
        try { prepared.emplace(annealer::prepare(demos::nonlinear())); }
        catch (...) { failed = 1; }
        MPI_Allreduce(&failed, &any_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if (any_failed) throw std::runtime_error("model preparation failed");

        annealer::SolverConfig config;
        config.seed = 42;
        config.temperature_levels = 16;
        config.sweeps_per_level = 2000;
        constexpr std::uint64_t starts = 7;
        MPI_Barrier(MPI_COMM_WORLD);
        double begin = MPI_Wtime();
        const auto multi = annealer::solve_multistart(MPI_COMM_WORLD, *prepared, starts, config);
        double elapsed = MPI_Wtime() - begin, maximum = 0;
        MPI_Reduce(&elapsed, &maximum, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (!rank) {
            std::cout << std::setprecision(17)
                << "mode,processes,base_seed,exchange_seed,starts,replicas,levels,sweeps_per_level,total_sweeps_per_replica,exchange_interval,cold,hot,swap_proposals,attempted_local_moves,seconds_max,best_energy,exchange_attempts,exchange_accepts,preparation_included\n"
                << "multistart," << size << ",42,0,7,0,16,2000,0,0,0.1,10,0,"
                << multi.attempted_moves << ',' << maximum << ',' << multi.best_energy << ",0,0,false\n";
        }

        annealer::ReplicaExchangeConfig exchange;
        exchange.solver = config;
        exchange.temperatures = annealer::geometric_ladder(0.1, 10.0, size);
        exchange.total_sweeps = 32000;
        exchange.exchange_interval = 100;
        exchange.exchange_seed = 900;
        MPI_Barrier(MPI_COMM_WORLD);
        begin = MPI_Wtime();
        const auto replicas = annealer::solve_replica_exchange(MPI_COMM_WORLD, *prepared, exchange);
        elapsed = MPI_Wtime() - begin;
        MPI_Reduce(&elapsed, &maximum, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        const std::uint64_t local_moves = replicas.local.attempted_moves;
        std::uint64_t total_moves = 0;
        MPI_Reduce(&local_moves, &total_moves, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
        if (!rank) {
            std::cout << "replica_exchange," << size << ",42,900,0," << size
                << ",0,0,32000,100,0.1,10,0," << total_moves << ',' << maximum << ','
                << replicas.best_energy << ',' << replicas.exchange_attempts << ','
                << replicas.exchange_accepts << ",false\n";
        }
    } catch (const std::exception& error) {
        std::cerr << "rank " << rank << ": " << error.what() << '\n';
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
    MPI_Finalize();
}
