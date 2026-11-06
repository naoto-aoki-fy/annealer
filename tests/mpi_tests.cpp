#include <annealer/mpi.hpp>
#include <iostream>

using namespace annealer;
namespace {
void require(bool ok) {
    if (!ok) throw std::runtime_error("MPI test assertion failed");
}
bool equal(const SolverResult& a, const SolverResult& b) {
    return a.best_assignment == b.best_assignment && a.best_energy == b.best_energy &&
        a.final_assignment == b.final_assignment && a.final_energy == b.final_energy &&
        a.attempted_moves == b.attempted_moves && a.accepted_moves == b.accepted_moves && a.seed == b.seed;
}
void compare(MPI_Comm comm, const Model& model, SolverConfig config, std::uint64_t starts) {
    const auto prepared = prepare(model);
    const auto distributed = solve_multistart(comm, prepared, starts, config);
    const auto unprepared = solve_multistart(comm, model, starts, config);
    require(distributed.best_assignment == unprepared.best_assignment);
    std::optional<SolverResult> best;
    std::uint64_t best_id = 0, attempted = 0, accepted = 0;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    std::size_t local_index = 0;
    for (std::uint64_t id = 0; id < starts; ++id) {
        auto trial = config;
        trial.seed += id;
        const auto serial = solve(prepared, trial);
        if (!best || serial.best_energy < best->best_energy) { best = serial; best_id = id; }
        attempted += serial.attempted_moves;
        accepted += serial.accepted_moves;
        if (id % size == static_cast<unsigned>(rank)) {
            require(local_index < distributed.local_runs.size());
            const auto& local = distributed.local_runs[local_index++];
            require(local.start_id == id && equal(local.result, serial));
        }
    }
    require(local_index == distributed.local_runs.size());
    require(distributed.start_id == best_id && distributed.seed == config.seed + best_id);
    require(distributed.best_energy == best->best_energy && distributed.best_assignment == best->best_assignment);
    require(distributed.attempted_moves == attempted && distributed.accepted_moves == accepted);
}
template<class F> void fails(MPI_Comm comm, F action) {
    int caught = 0;
    try { action(); } catch (const std::runtime_error&) { caught = 1; }
    int all = 0;
    MPI_Allreduce(&caught, &all, 1, MPI_INT, MPI_MIN, comm);
    require(all == 1);
}
void tests(MPI_Comm comm) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    ModelBuilder builder;
    for (int i = 0; i < 3; ++i) builder.add_variable(Domain{-3, 0, 7});
    builder.add_term(1, {{0, [](Value x) { return double(x); }},
                         {1, [](Value x) { return double(x); }}});
    const auto model = builder.build();
    SolverConfig config;
    config.seed = 35;
    config.temperature_levels = 4;
    config.sweeps_per_level = 7;
    for (auto starts : {1u, 2u, 5u}) compare(comm, model, config, starts);
    config.initial_assignment = Assignment{-3, 0, 7};
    compare(comm, model, config, 5);
    config.swap_proposals = true;
    compare(comm, model, config, 5);
    config.temperature_levels = 0;
    compare(comm, model, config, 2);
    config.temperature_levels = 4;
    config.sweeps_per_level = 0;
    compare(comm, model, config, 5);
    ModelBuilder flat;
    flat.add_variable(Domain{-2, 4});
    config = {};
    config.temperature_levels = 1;
    compare(comm, flat.build(), config, 5); // Every start ties: ID zero wins.
    compare(comm, ModelBuilder{}.build(), config, 2);
    fails(comm, [&] { solve_multistart(comm, model, 0, config); });
    config.seed = UINT64_MAX;
    fails(comm, [&] { solve_multistart(comm, model, 2, config); });
    compare(comm, model, config, 1); // UINT64_MAX is valid for one start.
    config = {};
    config.temperature_levels = 0;
    if (rank == size - 1) config.initial_assignment = Assignment{};
    fails(comm, [&] { solve_multistart(comm, model, 1, config); }); // Includes idle ranks.
    config = {};
    if (size > 1) {
        fails(comm, [&] { solve_multistart(comm, model, rank ? 2 : 3, config); });
        config.seed = rank;
        fails(comm, [&] { solve_multistart(comm, model, 2, config); });
        config = {};
        config.initial_assignment = Assignment{rank ? 0 : -3, 0, 7};
        fails(comm, [&] { solve_multistart(comm, model, 2, config); });
    }
    ModelBuilder throwing;
    throwing.add_variable(Domain{0, 1});
    throwing.add_term(1, {{0, [rank](Value x) -> double {
        if (rank == 0) throw std::runtime_error("injected callback failure");
        return double(x);
    }}});
    fails(comm, [&] { solve_multistart(comm, throwing.build(), 2); });
    bool fail_callback = false;
    ModelBuilder late_throwing;
    late_throwing.add_variable(Domain{0, 1});
    late_throwing.add_term(1, {{0, [&fail_callback, rank](Value x) -> double {
        if (fail_callback && rank == 0) throw std::runtime_error("solver callback failure");
        return double(x);
    }}});
    const auto prepared_throwing = prepare(late_throwing.build());
    fail_callback = true;
    fails(comm, [&] { solve_multistart(comm, prepared_throwing, 2); });
    // A successful call after failure checks that all participants returned together.
    compare(comm, model, SolverConfig{}, 1);
}
}
int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    try {
        tests(MPI_COMM_WORLD);
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm sub;
        MPI_Comm_split(MPI_COMM_WORLD, rank % 2, -rank, &sub);
        tests(sub);
        MPI_Comm_free(&sub);
        if (!rank) std::cout << "MPI multi-start tests passed\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Finalize();
}
