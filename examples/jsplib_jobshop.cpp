#include "jobshop_model.hpp"
#include "gset_maxcut.hpp"
#include <annealer/solver.hpp>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>

int main(int argc, char** argv) {
    try {
        if (argc < 2 || (argc == 2 && std::string(argv[1]) == "--help")) {
            std::cout << "Usage: jsplib_jobshop_example FILE [--seed N] [--trials N]\n"
                         "  [--initial-temperature T] [--final-temperature T] [--levels N] [--sweeps N] [--memory-mib N]\n"
                         "Defaults: seed=0 trials=10 temperatures=30,0.05 levels=200 sweeps=100\n";
            return argc < 2 ? 1 : 0;
        }
        annealer::SolverConfig config;
        config.initial_temperature = 30;
        config.final_temperature = 0.05;
        config.temperature_levels = 200;
        config.sweeps_per_level = 100;
        std::size_t trials = 10;
        jobshop::ModelConfig model_config;
        for (int i = 2; i < argc; ++i) {
            const std::string option = argv[i];
            if (++i == argc) throw std::invalid_argument("Missing value for " + option);
            const std::string value = argv[i];
            if (option == "--memory-mib") model_config.memory_bytes = jobshop::checked_multiply(gset::unsigned_number<std::size_t>(value), 1024*1024);
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
        std::ifstream file(argv[1]);
        if (!file) throw std::runtime_error("Cannot open JSPLib file");
        const auto instance = jobshop::read(file);
        const auto start = std::chrono::steady_clock::now();
        const auto plan = jobshop::estimate(instance);
        std::cout << "Horizon: " << plan.horizon << " Penalty: " << plan.penalty << " Estimated bytes: " << plan.memory_bytes
                  << " Variables: " << plan.variables << " Terms: " << plan.terms << std::endl;
        const auto model = jobshop::make_model(instance, model_config);
        config.initial_assignment = jobshop::initial_assignment(instance);
        const auto prepared = annealer::prepare(model, {.lookup_bytes = 1024 * 1024});
        double best = std::numeric_limits<double>::infinity();
        auto best_seed = config.seed;
        std::cout << std::setprecision(17) << "Variables: " << plan.variables << " Terms: " << plan.terms
                  << "\nInitial temperature: " << config.initial_temperature
                  << " Final temperature: " << config.final_temperature << " Levels: " << config.temperature_levels
                  << " Sweeps: " << config.sweeps_per_level << " Trials: " << trials << '\n';
        for (std::size_t trial = 0; trial < trials; ++trial) {
            const auto solve_start = std::chrono::steady_clock::now();
            const auto result = annealer::solve(prepared, config);
            const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - solve_start).count();
            const auto objective = jobshop::verify_model(instance, result.best_assignment);
            if (objective.energy != result.best_energy ||
                jobshop::verify_model(instance, result.final_assignment).energy != result.final_energy)
                throw std::runtime_error("Independent objective/energy verification failed");
            if (result.best_energy < best) { best = result.best_energy; best_seed = result.seed; }
            std::cout << "Seed: " << result.seed << " Actual makespan: " << objective.schedule.makespan
                      << " Constraint violations: " << (objective.schedule.precedence_violations + objective.schedule.machine_conflicts + objective.makespan_violations)
                      << " Precedence: " << objective.schedule.precedence_violations
                      << " Machine conflicts: " << objective.schedule.machine_conflicts
                      << " Makespan bounds: " << objective.makespan_violations
                      << " Declared makespan: " << objective.declared_makespan
                      << " Total energy: " << result.best_energy << " Solve seconds: " << seconds
                      << " Feasible: " << ((objective.schedule.precedence_violations + objective.schedule.machine_conflicts + objective.makespan_violations) == 0 ? "yes" : "no")
                      << " Status: " << ((objective.schedule.precedence_violations + objective.schedule.machine_conflicts + objective.makespan_violations) == 0 ? "feasible candidate (optimality unproven)" : "infeasible assignment") << '\n';
            if (trial + 1 < trials) ++config.seed;
        }
        std::cout << "Lowest energy seed: " << best_seed << " Lowest energy: " << best
                  << "\nTotal seconds (model, preparation, trials, verification, output): "
                  << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() << '\n';
    } catch (const std::exception& e) { std::cerr << "Error: " << e.what() << '\n'; return 1; }
}
