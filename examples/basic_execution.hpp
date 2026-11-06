#pragma once

#include <annealer/solver.hpp>
#ifdef ANNEALER_ENABLE_MPI
#include <annealer/replica_exchange.hpp>
#endif
#include <charconv>
#include <iostream>
#include <cmath>
#include <limits>
#include <optional>
#include <string>
#include <string_view>

namespace basic {
struct Options {
    annealer::SolverConfig solver;
    std::string mode = "serial";
    std::uint64_t trials = 1;
    std::size_t exchange_interval = 10;
    std::uint64_t exchange_seed = 0;
    bool help = false;
};

inline std::uint64_t integer(std::string_view value) {
    std::uint64_t result = 0;
    const auto parsed = std::from_chars(value.data(), value.data() + value.size(), result);
    if (parsed.ec != std::errc{} || parsed.ptr != value.data() + value.size())
        throw std::invalid_argument("invalid unsigned integer: " + std::string(value));
    return result;
}
inline std::size_t size_value(std::string_view value) {
    const auto result = integer(value);
    if (result > std::numeric_limits<std::size_t>::max())
        throw std::overflow_error("option exceeds size_t range");
    return static_cast<std::size_t>(result);
}
inline double temperature(const std::string& value) {
    std::size_t used = 0;
    const auto result = std::stod(value, &used);
    if (used != value.size() || !std::isfinite(result) || result <= 0)
        throw std::invalid_argument("temperature must be finite and positive");
    return result;
}
inline void parse(Options& options, int argc, char** argv) {
    bool exchange_option = false;
    for (int i = 1; i < argc; ++i) {
        const std::string key = argv[i];
        if (key == "--help") { options.help = true; continue; }
        if (++i == argc) throw std::invalid_argument("missing value for " + key);
        const std::string value = argv[i];
        if (key == "--mode") options.mode = value;
        else if (key == "--seed") options.solver.seed = integer(value);
        else if (key == "--trials") options.trials = integer(value);
        else if (key == "--initial-temperature") options.solver.initial_temperature = temperature(value);
        else if (key == "--final-temperature") options.solver.final_temperature = temperature(value);
        else if (key == "--levels") options.solver.temperature_levels = size_value(value);
        else if (key == "--sweeps") options.solver.sweeps_per_level = size_value(value);
        else if (key == "--exchange-interval") {
            exchange_option = true;
            options.exchange_interval = size_value(value);
        } else if (key == "--exchange-seed") {
            exchange_option = true;
            options.exchange_seed = integer(value);
        } else throw std::invalid_argument("unknown option: " + key);
    }
    if (options.mode != "serial" && options.mode != "mpi-multistart" && options.mode != "mpi-replica-exchange")
        throw std::invalid_argument("invalid --mode");
#ifndef ANNEALER_ENABLE_MPI
    if (options.mode != "serial")
        throw std::invalid_argument("MPI mode requires a build with ANNEALER_ENABLE_MPI=ON");
#endif
    if (exchange_option && options.mode != "mpi-replica-exchange")
        throw std::invalid_argument("exchange options require mpi-replica-exchange mode");
    if (!options.trials || !options.exchange_interval) throw std::invalid_argument("trials and exchange interval must be positive");
    annealer::detail::validate_temperatures(options.solver.initial_temperature, options.solver.final_temperature);
}

// All ranks enter each phase, including phases that only print on rank zero.
class Session {
    bool owns_mpi_ = false;
public:
    int rank = 0, size = 1;
    Session(int& argc, char**& argv) {
#ifdef ANNEALER_ENABLE_MPI
        int initialized = 0;
        MPI_Initialized(&initialized);
        if (!initialized) {
            if (MPI_Init(&argc, &argv) != MPI_SUCCESS) throw std::runtime_error("MPI_Init failed");
            owns_mpi_ = true;
        }
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
    }
    ~Session() {
#ifdef ANNEALER_ENABLE_MPI
        if (owns_mpi_) MPI_Finalize();
#endif
    }
    template<class F> void phase(F action) const {
        std::string error;
        bool failed = false;
        try { action(); }
        catch (const std::exception& e) { failed = true; error = e.what(); }
        catch (...) { failed = true; }
#ifdef ANNEALER_ENABLE_MPI
        int local = failed, any = 0;
        MPI_Allreduce(&local, &any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        failed = any;
#endif
        if (failed) throw std::runtime_error(error.empty() ? "example phase failed on another rank" : error);
    }
};

template<class MakeModel, class Before, class Report>
int run(int argc, char** argv, Options options, MakeModel make_model, Before before, Report report,
        annealer::PreparationConfig preparation = {}) {
    Session session(argc, argv);
    try {
        session.phase([&] { parse(options, argc, argv); });
#ifdef ANNEALER_ENABLE_MPI
        // MPMD invocations must agree before any mode-dependent collectives.
        annealer::mpi_detail::same(MPI_COMM_WORLD, options.help);
        annealer::mpi_detail::same(MPI_COMM_WORLD, options.mode == "serial" ? 0 : options.mode == "mpi-multistart" ? 1 : 2);
        annealer::mpi_detail::same(MPI_COMM_WORLD, options.trials);
#endif
        if (options.help) {
            session.phase([&] { if (!session.rank) std::cout
                << "--mode serial|mpi-multistart|mpi-replica-exchange\n"
                << "--seed N --trials N --initial-temperature T --final-temperature T --levels N --sweeps N\n"
                << "--exchange-interval N (default 10) --exchange-seed N (default 0), replica exchange only\n"
                << "Multi-start trials count global starts. Replica-exchange trials count collective runs;\n"
                << "levels * sweeps is the per-replica total budget, with a cold-to-hot geometric ladder.\n"; });
            return 0;
        }
        std::optional<annealer::Model> model;
        std::optional<annealer::PreparedModel> prepared;
        session.phase([&] { model.emplace(make_model()); });
        session.phase([&] { prepared.emplace(annealer::prepare(*model, preparation)); });
        session.phase([&] { if (!session.rank) before(*model); });
        session.phase([&] {
            const auto stride = options.mode == "mpi-replica-exchange" ? static_cast<std::uint64_t>(session.size) : 1;
            if (options.trials - 1 > (UINT64_MAX - (stride - 1)) / stride)
                throw std::overflow_error("trial seed offset overflow");
            const auto last = (options.trials - 1) * stride + (stride - 1);
            if (last > UINT64_MAX - options.solver.seed ||
                (options.mode == "mpi-replica-exchange" && last > UINT64_MAX - options.exchange_seed))
                throw std::overflow_error("trial seed overflow");
        });
#ifdef ANNEALER_ENABLE_MPI
        if (options.mode == "mpi-multistart") {
            const auto result = annealer::solve_multistart(MPI_COMM_WORLD, *prepared, options.trials, options.solver);
            session.phase([&] { if (!session.rank) {
                std::cout << "Global best energy: " << result.best_energy << " at";
                for (auto value : result.best_assignment) std::cout << ' ' << value;
                std::cout << "\nGlobal moves: " << result.accepted_moves << '/' << result.attempted_moves
                          << "; start: " << result.start_id << "; seed: " << result.seed << '\n';
            } });
            return 0;
        }
        if (options.mode == "mpi-replica-exchange") {
            annealer::ReplicaExchangeConfig config;
            session.phase([&] {
                config.solver = options.solver;
                config.temperatures = annealer::geometric_ladder(options.solver.final_temperature,
                    options.solver.initial_temperature, session.size);
                if (options.solver.sweeps_per_level && options.solver.temperature_levels >
                    std::numeric_limits<std::size_t>::max() / options.solver.sweeps_per_level)
                    throw std::overflow_error("replica sweep budget overflow");
                config.total_sweeps = options.solver.temperature_levels * options.solver.sweeps_per_level;
                config.exchange_interval = options.exchange_interval;
            });
            for (std::uint64_t t = 0; t < options.trials; ++t) {
                config.solver.seed = options.solver.seed + t * session.size;
                config.exchange_seed = options.exchange_seed + t * session.size;
                const auto result = annealer::solve_replica_exchange(MPI_COMM_WORLD, *prepared, config);
                session.phase([&] { if (!session.rank) {
                    std::cout << "Run " << t << " global best energy: " << result.best_energy << " at";
                    for (auto value : result.best_assignment) std::cout << ' ' << value;
                    std::cout << "\nGlobal exchanges: " << result.exchange_accepts << '/' << result.exchange_attempts
                              << "\nRank 0 local moves: " << result.local.accepted_moves << '/' << result.local.attempted_moves << '\n';
                } });
            }
            return 0;
        }
#endif
        for (std::uint64_t t = 0; t < options.trials; ++t) {
            std::optional<annealer::SolverResult> result;
            session.phase([&] {
                auto config = options.solver;
                config.seed += t;
                result.emplace(annealer::solve(*prepared, config));
            });
            session.phase([&] { if (!session.rank) report(*result); });
        }
        return 0;
    } catch (const std::exception& error) {
        if (!session.rank) std::cerr << error.what() << '\n';
        return 1;
    }
}
} // namespace basic
