#include "gset_maxcut.hpp"
#include <annealer/solver.hpp>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>

namespace {
double real_number(const std::string& text) {
    std::size_t end = 0;
    const double value = std::stod(text, &end);
    if (end != text.size() || !std::isfinite(value)) throw std::invalid_argument("Invalid number: " + text);
    return value;
}
void usage() {
    std::cout << "Usage: gset_g1_example FILE [--seed N] [--trials N]\n"
                 "  [--initial-temperature T] [--final-temperature T]\n"
                 "  [--levels N] [--sweeps N] [--reference-cut W]\n"
                 "Defaults: seed=0 trials=10 initial-temperature=30 final-temperature=0.05\n"
                 "          levels=200 sweeps=100; no reference cut by default.\n";
}
}

int main(int argc, char** argv) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--help") { usage(); return 0; }
        if (argc < 2) { usage(); return 1; }
        annealer::SolverConfig config;
        config.initial_temperature = 30;
        config.final_temperature = 0.05;
        config.temperature_levels = 200;
        config.sweeps_per_level = 100;
        std::size_t trials = 10;
        std::optional<double> reference;
        for (int i = 2; i < argc; i += 2) {
            const std::string option = argv[i];
            if (i + 1 == argc) throw std::invalid_argument("Missing value for " + option);
            const std::string value = argv[i + 1];
            if (option == "--seed") config.seed = gset::unsigned_number<std::uint64_t>(value);
            else if (option == "--trials") trials = gset::unsigned_number<std::size_t>(value);
            else if (option == "--levels") config.temperature_levels = gset::unsigned_number<std::size_t>(value);
            else if (option == "--sweeps") config.sweeps_per_level = gset::unsigned_number<std::size_t>(value);
            else if (option == "--initial-temperature") config.initial_temperature = real_number(value);
            else if (option == "--final-temperature") config.final_temperature = real_number(value);
            else if (option == "--reference-cut") reference = real_number(value);
            else throw std::invalid_argument("Unknown option: " + option);
        }
        if (!trials || trials - 1 > std::numeric_limits<std::uint64_t>::max() - config.seed) {
            throw std::invalid_argument("Trials must be positive and seed range must not overflow");
        }
        std::ifstream file(argv[1]);
        if (!file) throw std::runtime_error(std::string("Cannot open: ") + argv[1]);
        const auto graph = gset::read(file);
        const auto start = std::chrono::steady_clock::now();
        const auto model = gset::make_model(graph);
        const auto prepared = annealer::prepare(model, {.lookup_bytes = 1024 * 1024});
        std::cout << std::setprecision(17) << "Vertices: " << graph.vertices << " Edges: " << graph.edges.size()
                  << "\nInitial temperature: " << config.initial_temperature << " Final temperature: "
                  << config.final_temperature << " Levels: " << config.temperature_levels
                  << " Sweeps: " << config.sweeps_per_level << '\n';
        double best = std::numeric_limits<double>::infinity();
        std::uint64_t best_seed = config.seed;
        for (std::size_t trial = 0; trial < trials; ++trial, ++config.seed) {
            const auto solve_start = std::chrono::steady_clock::now();
            const auto result = annealer::solve(prepared, config);
            const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - solve_start).count();
            const double cut = gset::cut_weight(graph, result.best_assignment);
            if (std::abs(cut + result.best_energy) > 1e-10 * (1 + std::max(std::abs(cut), std::abs(result.best_energy)))) {
                throw std::runtime_error("Independent cut/energy verification failed");
            }
            if (result.best_energy < best) { best = result.best_energy; best_seed = result.seed; }
            std::cout << "Seed: " << result.seed << " Cut: " << cut << " Best energy: " << result.best_energy
                      << " Attempted: " << result.attempted_moves << " Accepted: " << result.accepted_moves
                      << " Solve seconds: " << seconds << " Verified: yes\n";
        }
        std::cout << "Best seed: " << best_seed << " Best cut: " << -best << " Best energy: " << best;
        if (reference) std::cout << " Reference gap: " << *reference + best;
        std::cout << "\nTotal seconds (model, preparation, trials, verification, output): "
                  << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() << '\n';
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
}
