#include <annealer/replica_exchange.hpp>
#include <iostream>
using namespace annealer;
void require(bool ok) { if (!ok) throw std::runtime_error("replica exchange assertion failed"); }
bool same(const SolverResult& a, const SolverResult& b) {
    return a.best_assignment == b.best_assignment && a.final_assignment == b.final_assignment &&
        a.best_energy == b.best_energy && a.final_energy == b.final_energy && a.seed == b.seed &&
        a.attempted_moves == b.attempted_moves && a.accepted_moves == b.accepted_moves;
}
void formula() {
    require(replica_exchange_accept(1, 2, 2, 0, .99));
    require(replica_exchange_accept(1, 2, 0, 2, .3));
    require(!replica_exchange_accept(1, 2, 0, 2, .4)); // Detect reversed sign.
    require(replica_exchange_accept(2, 1, 2, 0, .3));
    require(!replica_exchange_accept(2, 1, 2, 0, .4));
    require(replica_exchange_accept(1, 1, -9, 4, .99));
    require(replica_exchange_accept(1, 2, 4, 4, .99));
    require(replica_exchange_accept(1, 2, -4, 4, 0));
    const auto tiny = std::numeric_limits<double>::denorm_min();
    const auto huge = std::numeric_limits<double>::max();
    require(!replica_exchange_accept(tiny, huge, -huge, huge, .1));
    require(replica_exchange_accept(tiny, huge, huge, -huge, .99));
    for (double bad : {0., -1., double(INFINITY), double(NAN)}) {
        bool caught = false;
        try { replica_exchange_accept(bad, 2, 0, 1, .5); } catch (const std::invalid_argument&) { caught = true; }
        require(caught);
    }
    require(std::abs(geometric_ladder(1, 100, 3)[1] - 10) < 1e-12);
    require(geometric_ladder(1, 100, 1)[0] == 1);
    for (int n : {1, 2, 3, 4, 5}) for (int phase : {0, 1}) for (int r = 0; r < n; ++r) {
        auto partner = mpi_detail::exchange_partner(r, n, phase);
        if (partner == MPI_PROC_NULL) continue;
        require(std::abs(partner-r) == 1 && std::min(partner,r) % 2 == phase);
        require(mpi_detail::exchange_partner(partner, n, phase) == r);
    }
}
void compare(MPI_Comm comm, const PreparedModel& prepared, ReplicaExchangeConfig c) {
    int rank, size; MPI_Comm_rank(comm, &rank); MPI_Comm_size(comm, &size);
    auto result = solve_replica_exchange(comm, prepared, c);
    std::vector<std::unique_ptr<detail::ExecutionState<>>> chains;
    std::vector<std::mt19937_64> random;
    for (int r = 0; r < size; ++r) {
        auto local = c.solver; local.seed += r;
        local.initial_temperature = local.final_temperature = c.temperatures[r];
        chains.push_back(std::make_unique<detail::ExecutionState<>>(prepared, local));
        random.emplace_back(c.exchange_seed + r);
    }
    std::uint64_t attempts = 0, accepts = 0, round = 0, local_attempts = 0, local_accepts = 0;
    for (std::size_t done = 0; done < c.total_sweeps; ++round) {
        auto block = std::min(c.exchange_interval, c.total_sweeps - done); done += block;
        for (int r = 0; r < size; ++r) chains[r]->advance_sweeps(block, c.temperatures[r]);
        if (!c.exchanges) continue;
        for (int r = round % 2; r + 1 < size; r += 2) {
            ++attempts; if (rank == r || rank == r+1) ++local_attempts;
            const auto a = chains[r]->snapshot(), b = chains[r+1]->snapshot();
            const auto u = static_cast<double>(random[r]() >> 11) * 0x1.0p-53;
            if (replica_exchange_accept(c.temperatures[r], c.temperatures[r+1], a.final_energy, b.final_energy, u)) {
                ++accepts; if (rank == r || rank == r+1) ++local_accepts;
                chains[r]->install_assignment(b.final_assignment); chains[r+1]->install_assignment(a.final_assignment);
            }
        }
    }
    require(same(result.local, chains[rank]->finish()));
    require(result.exchange_attempts == attempts && result.exchange_accepts == accepts);
    if (size > 1 && c.exchanges && c.total_sweeps == 31) {
        if (c.temperatures.front() == c.temperatures.back()) require(accepts == attempts && accepts > 0);
        else if (!c.solver.swap_proposals) require(accepts < attempts);
    }
    require(result.local_exchange_attempts == local_attempts && result.local_exchange_accepts == local_accepts);
    auto best = chains[0]->finish(); int best_rank = 0;
    for (int r = 1; r < size; ++r) if (chains[r]->snapshot().best_energy < best.best_energy) {
        best = chains[r]->finish(); best_rank = r;
    }
    require(result.best_rank == best_rank && result.best_energy == best.best_energy && result.best_assignment == best.best_assignment);
    require(full_energy(prepared.source(), result.local.final_assignment) == result.local.final_energy);
    if (!c.exchanges) {
        auto serial = c.solver; serial.seed += rank;
        serial.initial_temperature = serial.final_temperature = c.temperatures[rank];
        serial.temperature_levels = 1; serial.sweeps_per_level = c.total_sweeps;
        require(same(result.local, solve(prepared, serial)));
    }
}
int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    try {
        formula();
        int rank, size; MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD, &size);
        ModelBuilder b;
        for (int i = 0; i < 4; ++i) b.add_variable(Domain{-2, 0, 3});
        b.add_term(1, {{0, [](Value x) { return double(x); }}, {1, [](Value x) { return double(x); }}});
        const auto prepared = prepare(b.build());
        ReplicaExchangeConfig c; c.temperatures = geometric_ladder(.01, 100, size);
        c.total_sweeps = 31; c.exchange_interval = 2; c.solver.seed = 42;
        for (bool swaps : {false, true}) {
            c.solver.swap_proposals = swaps;
            c.solver.initial_assignment = Assignment{-2, 0, 3, -2};
            for (bool enabled : {false, true}) { c.exchanges = enabled; compare(MPI_COMM_WORLD, prepared, c); }
            c.solver.initial_assignment.reset(); compare(MPI_COMM_WORLD, prepared, c);
        }
        c.temperatures.assign(size, 2); compare(MPI_COMM_WORLD, prepared, c); // All exchanges accepted.
        c.total_sweeps = 0; compare(MPI_COMM_WORLD, prepared, c);
        c.exchange_interval = rank == size - 1 ? 0 : 2;
        int caught = 0;
        try { solve_replica_exchange(MPI_COMM_WORLD, prepared, c); } catch (const std::runtime_error&) { caught = 1; }
        int all = 0; MPI_Allreduce(&caught, &all, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD); require(all);
        c.exchange_interval = 2;
        if (size > 1) {
            c.temperatures[0] += rank;
            caught = 0;
            try { solve_replica_exchange(MPI_COMM_WORLD, prepared, c); } catch (const std::runtime_error&) { caught = 1; }
            MPI_Allreduce(&caught, &all, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD); require(all);
        }
        int calls = 0;
        ModelBuilder failing;
        failing.add_variable(Domain{0, 1});
        failing.add_term(1, {{0, [&calls, rank](Value x) -> double {
            if (++calls > 20 && rank == 0) throw std::runtime_error("injected callback failure");
            return double(x);
        }}});
        const auto failing_prepared = prepare(failing.build());
        c.temperatures.assign(size, 1); c.total_sweeps = 100; c.exchange_interval = 1;
        caught = 0;
        try { solve_replica_exchange(MPI_COMM_WORLD, failing_prepared, c); } catch (const std::runtime_error&) { caught = 1; }
        MPI_Allreduce(&caught, &all, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD); require(all);
        MPI_Comm sub; MPI_Comm_split(MPI_COMM_WORLD, rank % 2, -rank, &sub);
        int sub_size; MPI_Comm_size(sub, &sub_size);
        c.temperatures = geometric_ladder(1, 10, sub_size); c.total_sweeps = 5;
        compare(sub, prepared, c); MPI_Comm_free(&sub);
        if (!rank) std::cout << "Replica exchange tests passed\n";
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; MPI_Abort(MPI_COMM_WORLD, 1); }
    MPI_Finalize();
}
