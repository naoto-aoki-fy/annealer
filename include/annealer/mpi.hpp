#pragma once

#include <annealer/solver.hpp>
#include <mpi.h>

#include <array>
#include <bit>
#include <climits>

namespace annealer {

struct MultiStartRun {
    std::uint64_t start_id;
    SolverResult result;
};

struct MultiStartResult {
    Assignment best_assignment;
    double best_energy = 0;
    std::uint64_t start_id = 0;
    std::uint64_t seed = 0;
    // Global sums over all starts; overflow is a collective error.
    std::uint64_t attempted_moves = 0;
    std::uint64_t accepted_moves = 0;
    // Only starts executed by this rank, in increasing start ID order.
    std::vector<MultiStartRun> local_runs;
};

namespace mpi_detail {
// MPI transport/process failures follow the communicator's MPI error handler.
inline void agree(MPI_Comm comm, bool failed) {
    int local = failed, any = 0;
    MPI_Allreduce(&local, &any, 1, MPI_INT, MPI_MAX, comm);
    if (any) throw std::runtime_error("MPI multi-start failed on one or more ranks");
}

template<class F>
inline void phase(MPI_Comm comm, F&& action) {
    bool failed = false;
    try { action(); } catch (...) { failed = true; }
    agree(comm, failed);
}

inline void same(MPI_Comm comm, std::uint64_t value) {
    auto root = value;
    MPI_Bcast(&root, 1, MPI_UINT64_T, 0, comm);
    agree(comm, root != value);
}

inline void validate(MPI_Comm comm, const Model& model, const SolverConfig& config,
                     std::uint64_t num_starts) {
    phase(comm, [&] {
        if (!num_starts || num_starts - 1 > UINT64_MAX - config.seed)
            throw std::invalid_argument("invalid start count or seed overflow");
        if (model.domains().size() > INT_MAX)
            throw std::invalid_argument("assignment exceeds MPI count range");
        detail::validate_temperatures(config.initial_temperature, config.final_temperature);
        detail::validate_swap_domains(model, config.swap_proposals);
        if (config.initial_assignment) full_energy(model, *config.initial_assignment);
    });
    // Compare fields explicitly, avoiding padding, hashes, and floating conversions.
    const std::array<std::uint64_t, 10> fields{
        num_starts, config.seed, std::bit_cast<std::uint64_t>(config.initial_temperature),
        std::bit_cast<std::uint64_t>(config.final_temperature), config.temperature_levels,
        config.sweeps_per_level, config.swap_proposals, config.initial_assignment.has_value(),
        config.initial_assignment ? config.initial_assignment->size() : 0, model.domains().size()};
    for (const auto field : fields) same(comm, field);
    if (config.initial_assignment)
        for (const auto value : *config.initial_assignment)
            same(comm, std::bit_cast<std::uint64_t>(value));
}
} // namespace mpi_detail

// Collective on a valid intracommunicator; all ranks must call in the same order.
inline MultiStartResult solve_multistart(MPI_Comm comm, const PreparedModel& prepared,
        std::uint64_t num_starts, const SolverConfig& config = {}) {
    int inter = 0;
    MPI_Comm_test_inter(comm, &inter);
    if (inter) throw std::invalid_argument("multi-start requires an intracommunicator");
    int rank = 0, size = 0;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    const auto& model = prepared.source();
    mpi_detail::validate(comm, model, config, num_starts);
    MultiStartResult output;
    std::size_t local_best = 0;
    mpi_detail::phase(comm, [&] {
        for (std::uint64_t id = rank; id < num_starts;) {
            auto trial = config;
            trial.seed = config.seed + id;
            output.local_runs.push_back({id, solve(prepared, trial)});
            const auto index = output.local_runs.size() - 1;
            if (output.local_runs[index].result.best_energy <
                    output.local_runs[local_best].result.best_energy) local_best = index;
            // Subtraction avoids overflow even at the uint64_t boundary.
            if (num_starts - id <= static_cast<std::uint64_t>(size)) break;
            id += size;
        }
        output.best_assignment.resize(model.domains().size());
    });
    int winner = 0;
    bool have_best = false;
    // Fixed-size broadcasts also let us accumulate counters with checked addition.
    for (int owner = 0; owner < size; ++owner) {
        std::uint64_t id = 0, attempted = 0, accepted = 0;
        double energy = 0;
        int present = 0;
        mpi_detail::phase(comm, [&] {
            if (rank != owner) return;
            present = !output.local_runs.empty();
            if (present) {
                id = output.local_runs[local_best].start_id;
                energy = output.local_runs[local_best].result.best_energy;
            }
            for (const auto& run : output.local_runs) {
                if (run.result.attempted_moves > UINT64_MAX - attempted ||
                    run.result.accepted_moves > UINT64_MAX - accepted)
                    throw std::overflow_error("multi-start counters overflow");
                attempted += run.result.attempted_moves;
                accepted += run.result.accepted_moves;
            }
        });
        MPI_Bcast(&present, 1, MPI_INT, owner, comm);
        MPI_Bcast(&id, 1, MPI_UINT64_T, owner, comm);
        MPI_Bcast(&energy, 1, MPI_DOUBLE, owner, comm);
        MPI_Bcast(&attempted, 1, MPI_UINT64_T, owner, comm);
        MPI_Bcast(&accepted, 1, MPI_UINT64_T, owner, comm);
        mpi_detail::phase(comm, [&] {
            if (attempted > UINT64_MAX - output.attempted_moves ||
                accepted > UINT64_MAX - output.accepted_moves)
                throw std::overflow_error("multi-start counters overflow");
            output.attempted_moves += attempted;
            output.accepted_moves += accepted;
        });
        if (present && (!have_best || energy < output.best_energy ||
                (energy == output.best_energy && id < output.start_id))) {
            have_best = true;
            winner = owner;
            output.best_energy = energy;
            output.start_id = id;
        }
    }
    mpi_detail::phase(comm, [&] {
        if (rank == winner)
            std::copy(output.local_runs[local_best].result.best_assignment.begin(),
                      output.local_runs[local_best].result.best_assignment.end(),
                      output.best_assignment.begin());
    });
    MPI_Bcast(output.best_assignment.data(), static_cast<int>(output.best_assignment.size()),
              MPI_INT64_T, winner, comm);
    mpi_detail::phase(comm, [&] {
        if (full_energy(model, output.best_assignment) != output.best_energy)
            throw std::runtime_error("winning energy differs on local model");
    });
    output.seed = config.seed + output.start_id;
    return output;
}

inline MultiStartResult solve_multistart(MPI_Comm comm, const Model& model,
        std::uint64_t num_starts, const SolverConfig& config = {}) {
    std::optional<PreparedModel> prepared;
    mpi_detail::phase(comm, [&] { prepared.emplace(prepare(model)); });
    return solve_multistart(comm, *prepared, num_starts, config);
}

} // namespace annealer
