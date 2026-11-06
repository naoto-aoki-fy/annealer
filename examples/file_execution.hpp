#pragma once
#include "basic_execution.hpp"
#include <sstream>
#include <vector>
#include <chrono>

namespace file_example {
template<class F> auto phase_value(const basic::Session& session, F action) {
    std::optional<decltype(action())> value;
    session.phase([&] { value.emplace(action()); });
    return std::move(*value);
}
// Parse only execution options here; each example retains its own validation.
inline bool option(basic::Options& o, const std::string& key, const std::string& value,
                   bool& exchange) {
    if (key == "--mode") o.mode = value;
    else if (key == "--exchange-interval") { exchange = true; o.exchange_interval = basic::size_value(value); }
    else if (key == "--exchange-seed") { exchange = true; o.exchange_seed = basic::integer(value); }
    else return false;
    return true;
}
inline void validate(const basic::Session& s, basic::Options& o, bool exchange) {
    s.phase([&] {
        if (o.mode != "serial" && o.mode != "mpi-multistart" && o.mode != "mpi-replica-exchange")
            throw std::invalid_argument("invalid --mode");
#ifndef ANNEALER_ENABLE_MPI
        if (o.mode != "serial") throw std::invalid_argument("MPI mode requires ANNEALER_ENABLE_MPI=ON");
#endif
        if (!o.exchange_interval || (exchange && o.mode != "mpi-replica-exchange"))
            throw std::invalid_argument("invalid exchange options");
        const auto stride = o.mode == "mpi-replica-exchange" ? static_cast<std::uint64_t>(s.size) : 1;
        if (!o.trials || o.trials - 1 > (UINT64_MAX - (stride - 1)) / stride) throw std::overflow_error("trial offset overflow");
        const auto last = (o.trials - 1) * stride + (stride - 1);
        if (last > UINT64_MAX - o.solver.seed ||
            (o.mode == "mpi-replica-exchange" && last > UINT64_MAX - o.exchange_seed))
            throw std::overflow_error("trial seed overflow");
    });
#ifdef ANNEALER_ENABLE_MPI
    annealer::mpi_detail::same(MPI_COMM_WORLD, o.mode == "serial" ? 0 : o.mode == "mpi-multistart" ? 1 : 2);
    annealer::mpi_detail::same(MPI_COMM_WORLD, o.trials);
#endif
}
template<class Verify> std::vector<annealer::SolverResult> solve(const basic::Session& s,
        const annealer::PreparedModel& prepared, const basic::Options& o, std::vector<double>& seconds, std::ostream& output, Verify verify) {
    std::vector<annealer::SolverResult> results;
    auto started = std::chrono::steady_clock::now();
    auto elapsed = [&] { seconds.push_back(std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count()); };
    auto save = [&](const auto& assignment, double energy, std::uint64_t seed) {
        annealer::SolverResult r;
        r.best_assignment = r.final_assignment = assignment;
        r.best_energy = r.final_energy = energy; r.seed = seed;
        verify(r); elapsed(); results.push_back(std::move(r));
    };
#ifdef ANNEALER_ENABLE_MPI
    if (o.mode == "mpi-multistart") {
        const auto r = annealer::solve_multistart(MPI_COMM_WORLD, prepared, o.trials, o.solver);
        s.phase([&] {
            for (const auto& run : r.local_runs) verify(run.result);
            save(r.best_assignment, r.best_energy, r.seed);
            output << "Global best energy: " << r.best_energy
                   << "\nGlobal moves: " << r.accepted_moves << '/' << r.attempted_moves
                   << "; start: " << r.start_id << "; seed: " << r.seed << '\n';
        });
        return results;
    }
    if (o.mode == "mpi-replica-exchange") {
        auto config = phase_value(s, [&] {
            annealer::ReplicaExchangeConfig c;
            c.solver = o.solver;
            c.temperatures = annealer::geometric_ladder(o.solver.final_temperature, o.solver.initial_temperature, s.size);
            if (o.solver.sweeps_per_level && o.solver.temperature_levels > SIZE_MAX / o.solver.sweeps_per_level)
                throw std::overflow_error("replica sweep budget overflow");
            c.total_sweeps = o.solver.temperature_levels * o.solver.sweeps_per_level;
            c.exchange_interval = o.exchange_interval;
            return c;
        });
        for (std::uint64_t t = 0; t < o.trials; ++t) {
            started = std::chrono::steady_clock::now();
            config.solver.seed = o.solver.seed + t * s.size;
            config.exchange_seed = o.exchange_seed + t * s.size;
            const auto r = annealer::solve_replica_exchange(MPI_COMM_WORLD, prepared, config);
            s.phase([&] {
                verify(r.local);
                save(r.best_assignment, r.best_energy, config.solver.seed + r.best_rank);
                output << "Run " << t << " global best energy: " << r.best_energy
                       << "\nGlobal exchanges: " << r.exchange_accepts << '/' << r.exchange_attempts << '\n';
            });
        }
        return results;
    }
#endif
    s.phase([&] {
        for (std::uint64_t t = 0; t < o.trials; ++t) {
            started = std::chrono::steady_clock::now();
            auto config = o.solver; config.seed += t;
            auto r = annealer::solve(prepared, config);
            verify(r); elapsed(); results.push_back(std::move(r));
        }
    });
    return results;
}
} // namespace file_example
