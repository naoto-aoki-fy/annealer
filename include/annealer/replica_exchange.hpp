#pragma once

#include <annealer/mpi.hpp>

namespace annealer {

inline std::vector<double> geometric_ladder(double cold, double hot, std::size_t count) {
    detail::validate_temperatures(hot, cold);
    if (!count) throw std::invalid_argument("temperature ladder must not be empty");
    std::vector<double> ladder(count);
    for (std::size_t i = 0; i < count; ++i)
        ladder[i] = geometric_temperature(hot, cold, count, count - 1 - i);
    // A single replica uses the cold endpoint.
    if (count == 1) ladder[0] = cold;
    return ladder;
}

// Log-domain test; scaling differences avoids overflow at extreme finite inputs.
inline bool replica_exchange_accept(double ti, double tj, double ei, double ej, double uniform) {
    if (!std::isfinite(ti) || !std::isfinite(tj) || ti <= 0 || tj <= 0 ||
        !std::isfinite(ei) || !std::isfinite(ej) || !std::isfinite(uniform) || uniform < 0 || uniform >= 1)
        throw std::invalid_argument("invalid replica exchange inputs");
    if (ti == tj || ei == ej || uniform == 0 || ((tj > ti) == (ei > ej))) return true;
    auto log_difference = [](double a, double b) {
        const long double scale = std::max(std::abs(a), std::abs(b));
        return std::log(scale) + std::log(std::abs(static_cast<long double>(a) / scale -
                                                                 static_cast<long double>(b) / scale));
    };
    const long double log_magnitude = log_difference(tj, ti) - std::log(static_cast<long double>(ti)) -
        std::log(static_cast<long double>(tj)) + log_difference(ei, ej);
    return std::log(-std::log(static_cast<long double>(uniform))) > log_magnitude;
}

namespace mpi_detail {
inline int exchange_partner(int rank, int size, std::uint64_t round) {
    const int parity = static_cast<int>(round % 2);
    const int partner = (rank % 2 == parity) ? rank + 1 : rank - 1;
    return partner < parity || partner >= size ? MPI_PROC_NULL : partner;
}
}

struct ReplicaExchangeConfig {
    SolverConfig solver;
    std::vector<double> temperatures;
    std::size_t total_sweeps = 100;
    std::size_t exchange_interval = 10;
    bool exchanges = true;
    std::uint64_t exchange_seed = 0;
};
struct ReplicaExchangeResult {
    SolverResult local;
    Assignment best_assignment;
    double best_energy = 0;
    int best_rank = 0;
    // Global pair counts, each pair counted once by its lower rank.
    std::uint64_t exchange_attempts = 0, exchange_accepts = 0;
    std::uint64_t local_exchange_attempts = 0, local_exchange_accepts = 0;
};

inline ReplicaExchangeResult solve_replica_exchange(MPI_Comm comm, const PreparedModel& prepared,
                                                    const ReplicaExchangeConfig& config) {
    int inter = 0, rank = 0, size = 0;
    MPI_Comm_test_inter(comm, &inter);
    if (inter) throw std::invalid_argument("replica exchange requires an intracommunicator");
    MPI_Comm_rank(comm, &rank); MPI_Comm_size(comm, &size);
    mpi_detail::validate(comm, prepared.source(), config.solver, size);
    mpi_detail::phase(comm, [&] {
        if (config.temperatures.size() != static_cast<std::size_t>(size) || !config.exchange_interval ||
            static_cast<std::uint64_t>(size - 1) > UINT64_MAX - config.exchange_seed)
            throw std::invalid_argument("invalid replica exchange configuration");
        double previous = 0;
        for (auto t : config.temperatures) {
            if (!std::isfinite(t) || t <= 0 || t < previous)
                throw std::invalid_argument("invalid temperature ladder");
            previous = t;
        }
    });
    mpi_detail::same(comm, config.total_sweeps);
    mpi_detail::same(comm, config.exchange_interval);
    mpi_detail::same(comm, config.exchanges);
    mpi_detail::same(comm, config.exchange_seed);
    for (auto t : config.temperatures) mpi_detail::same(comm, std::bit_cast<std::uint64_t>(t));
    std::optional<detail::ExecutionState<>> chain;
    Assignment received;
    ReplicaExchangeResult output;
    std::mt19937_64 exchange_random(config.exchange_seed + rank);
    mpi_detail::phase(comm, [&] {
        auto local = config.solver;
        local.seed += rank;
        local.initial_temperature = local.final_temperature = config.temperatures[rank];
        chain.emplace(prepared, local);
        detail::checked_count(chain->attempts_per_sweep(), config.total_sweeps);
        received.resize(prepared.source().domains().size());
        output.best_assignment.resize(received.size());
    });
    std::uint64_t owned_attempts = 0, owned_accepts = 0, round = 0;
    for (std::size_t done = 0; done < config.total_sweeps; ++round) {
        const auto block = std::min(config.exchange_interval, config.total_sweeps - done);
        mpi_detail::phase(comm, [&] { chain->advance_sweeps(block, config.temperatures[rank]); });
        done += block;
        if (!config.exchanges) continue;
        const int partner = mpi_detail::exchange_partner(rank, size, round);
        double energy = 0, other = 0;
        mpi_detail::phase(comm, [&] { energy = full_energy(prepared.source(), chain->assignment()); });
        MPI_Sendrecv(&energy, 1, MPI_DOUBLE, partner, 0, &other, 1, MPI_DOUBLE, partner, 0, comm, MPI_STATUS_IGNORE);
        int accept = 0;
        mpi_detail::phase(comm, [&] {
            if (partner != MPI_PROC_NULL && rank < partner) {
                const double uniform = static_cast<double>(exchange_random() >> 11) * 0x1.0p-53;
                accept = replica_exchange_accept(config.temperatures[rank], config.temperatures[partner], energy, other, uniform);
                ++owned_attempts; owned_accepts += accept;
            }
        });
        int decision = 0;
        MPI_Sendrecv(&accept, 1, MPI_INT, partner, 1, &decision, 1, MPI_INT, partner, 1, comm, MPI_STATUS_IGNORE);
        if (partner != MPI_PROC_NULL) {
            if (rank > partner) accept = decision;
            ++output.local_exchange_attempts; output.local_exchange_accepts += accept;
        }
        if (accept) MPI_Sendrecv(chain->assignment().data(), static_cast<int>(received.size()), MPI_INT64_T, partner, 2,
            received.data(), static_cast<int>(received.size()), MPI_INT64_T, partner, 2, comm, MPI_STATUS_IGNORE);
        mpi_detail::phase(comm, [&] { if (accept) chain->install_assignment(received); });
    }
    mpi_detail::phase(comm, [&] { output.local = chain->finish(); });
    for (int owner = 0; owner < size; ++owner) {
        double energy = output.local.best_energy;
        std::uint64_t attempts = owned_attempts, accepts = owned_accepts;
        MPI_Bcast(&energy, 1, MPI_DOUBLE, owner, comm);
        MPI_Bcast(&attempts, 1, MPI_UINT64_T, owner, comm);
        MPI_Bcast(&accepts, 1, MPI_UINT64_T, owner, comm);
        mpi_detail::phase(comm, [&] {
            if (attempts > UINT64_MAX - output.exchange_attempts || accepts > UINT64_MAX - output.exchange_accepts)
                throw std::overflow_error("exchange counters overflow");
            output.exchange_attempts += attempts; output.exchange_accepts += accepts;
        });
        if (owner == 0 || energy < output.best_energy) { output.best_energy = energy; output.best_rank = owner; }
    }
    if (rank == output.best_rank) std::copy(output.local.best_assignment.begin(), output.local.best_assignment.end(), output.best_assignment.begin());
    MPI_Bcast(output.best_assignment.data(), static_cast<int>(received.size()), MPI_INT64_T, output.best_rank, comm);
    mpi_detail::phase(comm, [&] {
        if (full_energy(prepared.source(), output.best_assignment) != output.best_energy)
            throw std::runtime_error("global replica energy mismatch");
    });
    return output;
}
inline ReplicaExchangeResult solve_replica_exchange(MPI_Comm comm, const Model& model,
                                                    const ReplicaExchangeConfig& config) {
    std::optional<PreparedModel> prepared;
    mpi_detail::phase(comm, [&] { prepared.emplace(prepare(model)); });
    return solve_replica_exchange(comm, *prepared, config);
}
} // namespace annealer
