#include "tsplib.hpp"
#include "file_execution.hpp"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>

int main(int argc, char** argv) {
    basic::Session session(argc, argv);
    std::ostringstream output;
    try {
#ifdef ANNEALER_ENABLE_MPI
        annealer::mpi_detail::same(MPI_COMM_WORLD, argc < 2);
        annealer::mpi_detail::same(MPI_COMM_WORLD, argc == 2 && std::string(argv[1]) == "--help");
#endif
        if (argc < 2 || (argc == 2 && std::string(argv[1]) == "--help")) {
            output << "Usage: tsplib_example FILE [--seed N] [--trials N]\n"
                         "  [--initial-temperature T] [--final-temperature T] [--levels N] [--sweeps N] [--reference-cost N]\n"
                         "Defaults: seed=0 trials=10 temperatures=30,0.05 levels=200 sweeps=100\n";
            if (!session.rank) std::cout << output.str()
                << "--mode serial|mpi-multistart|mpi-replica-exchange --exchange-interval N --exchange-seed N\n";
            return argc < 2 ? 1 : 0;
        }
        annealer::SolverConfig config;
        config.initial_temperature = 30;
        config.final_temperature = 0.05;
        config.temperature_levels = 200;
        config.sweeps_per_level = 100;
        std::size_t trials = 10;
        std::optional<tsp::Cost> reference;
        config.swap_proposals = true;
        basic::Options execution;
        bool exchange = false;
        session.phase([&] {
            for (int i = 2; i < argc; ++i) {
                const std::string option = argv[i];
                if (++i == argc) throw std::invalid_argument("Missing value for " + option);
                const std::string value = argv[i];
                if (file_example::option(execution, option, value, exchange)) continue;
                if (option == "--reference-cost") reference = gset::unsigned_number<tsp::Cost>(value);
                else if (option == "--seed") config.seed = gset::unsigned_number<std::uint64_t>(value);
                else if (option == "--trials") trials = gset::unsigned_number<std::size_t>(value);
                else if (option == "--levels") config.temperature_levels = gset::unsigned_number<std::size_t>(value);
                else if (option == "--sweeps") config.sweeps_per_level = gset::unsigned_number<std::size_t>(value);
                else if (option == "--initial-temperature" || option == "--final-temperature") {
                    std::size_t end = 0;
                    const double number = std::stod(value, &end);
                    if (end != value.size() || !std::isfinite(number) || number <= 0)
                        throw std::invalid_argument("Temperature must be positive and finite");
                    if (option == "--initial-temperature") config.initial_temperature = number;
                    else config.final_temperature = number;
                } else throw std::invalid_argument("Unknown option: " + option);
            }
            if (!trials || trials - 1 > std::numeric_limits<std::uint64_t>::max() - config.seed)
                throw std::invalid_argument("Trials must be positive without seed overflow");
            if (!config.temperature_levels || !config.sweeps_per_level)
                throw std::invalid_argument("Levels and sweeps must be positive");
            if (config.final_temperature > config.initial_temperature)
                throw std::invalid_argument("Final temperature exceeds initial temperature");
            if (reference && *reference > tsp::limit) throw std::invalid_argument("Reference exceeds exact range");
        });
        execution.solver = config; execution.trials = trials;
        file_example::validate(session, execution, exchange);
        const auto instance = file_example::phase_value(session, [&] {
            std::ifstream file(argv[1]);
            if (!file) throw std::runtime_error("Cannot open TSPLIB file");
            return tsp::read(file);
        });
        const auto start = std::chrono::steady_clock::now();
        const auto model = file_example::phase_value(session, [&] { return tsp::make_model(instance); });
        const auto prepared = file_example::phase_value(session, [&] { return annealer::prepare(model, {.lookup_bytes = 1024 * 1024}); });
        double best = std::numeric_limits<double>::infinity();
        auto best_seed = config.seed;
        output << std::setprecision(17) << "Variables: " << instance.size << " Terms: " << model.terms().size()
                  << "\nInitial temperature: " << config.initial_temperature
                  << " Final temperature: " << config.final_temperature << " Levels: " << config.temperature_levels
                  << " Sweeps: " << config.sweeps_per_level << " Trials: " << trials << '\n';
        output << "Initialization base seed: " << config.seed << "; mode: " << execution.mode << '\n';
        std::vector<double> solve_seconds;
        const auto results = file_example::solve(session, prepared, execution, solve_seconds, output, [&](const auto& result) {
            const auto objective = tsp::evaluate(instance, result.best_assignment);
            if (tsp::energy(instance, objective) != result.best_energy ||
                tsp::energy(instance, tsp::evaluate(instance, result.final_assignment)) != result.final_energy)
                throw std::runtime_error("Independent objective/energy verification failed");
        }, [&](auto& run) {
            run.initial_assignment = tsp::permutation(instance.size, run.seed);
        });
        session.phase([&] {
            std::size_t record = 0;
            for (const auto& result : results) {
                const double seconds = solve_seconds[record++];
                if (execution.mode != "serial") output << "Global record: " << record - 1 << " ";
                const auto objective = tsp::evaluate(instance, result.best_assignment);
                if (result.best_energy < best) { best = result.best_energy; best_seed = result.seed; }
                output << "Seed: " << result.seed << " Tour length: " << objective.cost
                          << " Collisions: " << objective.collisions
                          << " Total energy: " << result.best_energy << " Solve seconds: " << seconds
                          << " Valid tour: " << (objective.collisions == 0 ? "yes" : "no")
                          << " Status: " << (objective.collisions == 0 ? "feasible candidate (optimality unproven)" : "infeasible assignment") << '\n';
                if (reference && objective.feasible)
                    output << "Reference cost: " << *reference << " Objective gap: "
                              << static_cast<double>(objective.cost) - static_cast<double>(*reference) << '\n';

            }
        });
        output << "Lowest energy seed: " << best_seed << " Lowest energy: " << best
                  << "\nTotal seconds (model, preparation, trials, verification, output): "
                  << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() << '\n';
        if (!session.rank) std::cout << output.str();
    } catch (const std::exception& e) { if (!session.rank) std::cerr << "Error: " << e.what() << '\n'; return 1; }
}
