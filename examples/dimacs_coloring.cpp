#include "dimacs_coloring.hpp"
#include <annealer/solver.hpp>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>

int main(int argc, char** argv) {
    try {
        if (argc < 2 || (argc == 2 && std::string(argv[1]) == "--help")) {
            std::cout << "Usage: dimacs_coloring_example FILE [--colors K] [--seed N] [--trials N]\n"
                         "  [--initial-temperature T] [--final-temperature T] [--levels N] [--sweeps N]\n"
                         "Defaults: colors=5 seed=0 trials=10 temperatures=30,0.05 levels=200 sweeps=100\n";
            return argc < 2 ? 1 : 0;
        }
        annealer::SolverConfig config;
        config.initial_temperature = 30;
        config.final_temperature = 0.05;
        config.temperature_levels = 200;
        config.sweeps_per_level = 100;
        std::size_t colors = 5, trials = 10;
        for (int i = 2; i < argc; ++i) {
            const std::string option = argv[i];
            if (++i == argc) throw std::invalid_argument("Missing value for " + option);
            const std::string value = argv[i];
            if (option == "--colors") colors = gset::unsigned_number<std::size_t>(value);
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
        coloring::validate_colors(colors);
        if (!trials || trials - 1 > std::numeric_limits<std::uint64_t>::max() - config.seed)
            throw std::invalid_argument("Trials must be positive without seed overflow");
        if (config.final_temperature > config.initial_temperature)
            throw std::invalid_argument("Final temperature exceeds initial temperature");
        std::ifstream file(argv[1]);
        if (!file) throw std::runtime_error("Cannot open graph file");
        const auto graph = coloring::read(file);
        const auto start = std::chrono::steady_clock::now();
        const auto model = coloring::make_model(graph, colors);
        const auto prepared = annealer::prepare(model, {.lookup_bytes = 1024 * 1024});
        double best = std::numeric_limits<double>::infinity();
        auto best_seed = config.seed;
        std::cout << std::setprecision(17) << "Vertices: " << graph.vertices << " Edges: " << graph.edges.size()
                  << " Colors: " << colors << "\nInitial temperature: " << config.initial_temperature
                  << " Final temperature: " << config.final_temperature << " Levels: " << config.temperature_levels
                  << " Sweeps: " << config.sweeps_per_level << " Trials: " << trials << '\n';
        for (std::size_t trial = 0; trial < trials; ++trial) {
            const auto solve_start = std::chrono::steady_clock::now();
            const auto result = annealer::solve(prepared, config);
            const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - solve_start).count();
            const auto count = coloring::conflicts(graph, colors, result.best_assignment);
            if (static_cast<double>(count) != result.best_energy ||
                static_cast<double>(coloring::conflicts(graph, colors, result.final_assignment)) != result.final_energy)
                throw std::runtime_error("Independent conflict/energy verification failed");
            if (result.best_energy < best) { best = result.best_energy; best_seed = result.seed; }
            std::cout << "Seed: " << result.seed << " Conflicting edges: " << count
                      << " Best energy: " << result.best_energy << " Solve seconds: " << seconds
                      << " Status: " << (count == 0 ? "valid K-coloring" : "assignment contains conflicts") << '\n';
            if (trial + 1 < trials) ++config.seed;
        }
        std::cout << "Best seed: " << best_seed << " Best energy: " << best << " Conflicting edges: " << best
                  << " Status: " << (best == 0 ? "valid K-coloring" : "assignment contains conflicts")
                  << "\nTotal seconds (model, preparation, trials, verification, output): "
                  << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() << '\n';
    } catch (const std::exception& e) { std::cerr << "Error: " << e.what() << '\n'; return 1; }
}
