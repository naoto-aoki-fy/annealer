#include "gset_measurement.hpp"
#include "file_execution.hpp"
#include <annealer/solver.hpp>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <numeric>
#include <algorithm>

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
                 "  [--target-cut W] [--trace-csv PATH] [--profile]\n"
                 "  [--mode serial|mpi-multistart|mpi-replica-exchange]\n"
                 "  [--exchange-interval N] [--exchange-seed N] (replica exchange only)\n"
                 "Trace CSV, profile, and target-cut currently require serial mode.\n"
                 "Multi-start trials are global starts; replica trials are collective runs,\n"
                 "with levels*sweeps per replica and a cold-to-hot geometric ladder.\n"
                 "Defaults: seed=0 trials=10 initial-temperature=30 final-temperature=0.05\n"
                 "          levels=200 sweeps=100; no reference cut by default.\n";
}
}

int main(int argc, char** argv) {
    basic::Session session(argc, argv);
    try {
        if (argc == 2 && std::string(argv[1]) == "--help") { if (!session.rank) usage(); return 0; }
        if (argc < 2) { if (!session.rank) usage(); return 1; }
        basic::Options execution;
        auto& config = execution.solver;
        bool exchange = false;
        config.initial_temperature = 30;
        config.final_temperature = 0.05;
        config.temperature_levels = 200;
        config.sweeps_per_level = 100;
        execution.trials = 10;
        auto& trials = execution.trials;
        std::optional<double> reference, target;
        std::string trace;
        bool profile = false, trace_requested = false;
        session.phase([&] {
            for (int i = 2; i < argc; ++i) {
                const std::string option = argv[i];
                if (option == "--profile") { profile = true; continue; }
                if (i + 1 == argc) throw std::invalid_argument("Missing value for " + option);
                const std::string value = argv[++i];
                if (file_example::option(execution, option, value, exchange)) continue;
                if (option == "--seed") config.seed = gset::unsigned_number<std::uint64_t>(value);
                else if (option == "--trials") trials = gset::unsigned_number<std::size_t>(value);
                else if (option == "--levels") config.temperature_levels = gset::unsigned_number<std::size_t>(value);
                else if (option == "--sweeps") config.sweeps_per_level = gset::unsigned_number<std::size_t>(value);
                else if (option == "--initial-temperature") config.initial_temperature = real_number(value);
                else if (option == "--final-temperature") config.final_temperature = real_number(value);
                else if (option == "--reference-cut") reference = real_number(value);
                else if (option == "--target-cut") target = real_number(value);
                else if (option == "--trace-csv") { trace = value; trace_requested = true; }
                else throw std::invalid_argument("Unknown option: " + option);
            }
            if (execution.mode != "serial") {
                if (trace_requested) throw std::invalid_argument("--trace-csv currently requires serial mode");
                if (profile) throw std::invalid_argument("--profile currently requires serial mode");
                if (target) throw std::invalid_argument("--target-cut currently requires serial mode");
            }
        });
        file_example::validate(session, execution, exchange);
        const auto graph = file_example::phase_value(session, [&] {
            std::ifstream file(argv[1]);
            if (!file) throw std::runtime_error(std::string("Cannot open: ") + argv[1]);
            return gset::read(file);
        });
        const auto start = std::chrono::steady_clock::now();
        const auto model = file_example::phase_value(session, [&] { return gset::make_model(graph); });
        const auto prepared = file_example::phase_value(session, [&] {
            return annealer::prepare(model, {.lookup_bytes = 1024 * 1024});
        });
        if (execution.mode != "serial") {
            std::vector<double> seconds;
            std::ostringstream output;
            const auto results = file_example::solve(session, prepared, execution, seconds, output,
                [&](const auto& result) { gset::verify(graph, result); });
            session.phase([&] { if (!session.rank) {
                const auto& best = *std::min_element(results.begin(), results.end(),
                    [](const auto& a, const auto& b) { return a.best_energy < b.best_energy; });
                const double cut = gset::cut_weight(graph, best.best_assignment);
                std::cout << std::setprecision(17)
                          << "Best seed: " << best.seed << " Best cut: " << cut
                          << " Best energy: " << best.best_energy << " Verified: yes";
                if (reference) std::cout << " Reference gap: " << *reference - cut;
                std::cout << "\nMPI wall seconds (rank 0; solving and verification): "
                          << std::accumulate(seconds.begin(), seconds.end(), 0.0) << '\n';
            } });
            return 0;
        }
        session.phase([&] { if (!session.rank) {
            std::cout << std::setprecision(17) << "Vertices: " << graph.vertices << " Edges: " << graph.edges.size()
                      << "\nInitial temperature: " << config.initial_temperature << " Final temperature: "
                      << config.final_temperature << " Levels: " << config.temperature_levels
                      << " Sweeps: " << config.sweeps_per_level << '\n';
            std::ofstream csv;
            if (!trace.empty()) {
                csv.open(trace);
                if (!csv) throw std::runtime_error("Cannot open trace CSV");
                csv.exceptions(std::ios::badbit | std::ios::failbit);
                gset::csv_header(csv);
            }
            double best = std::numeric_limits<double>::infinity();
            std::uint64_t best_seed = config.seed;
            for (std::size_t trial = 0; trial < trials; ++trial, ++config.seed) {
                auto measured = annealer::solve_instrumented(prepared, config, {.profile = profile});
                const auto& result = measured.result;
                if (profile) {
                    const auto normal_start = std::chrono::steady_clock::now();
                    const auto normal = annealer::solve(prepared, config);
                    const double normal_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-normal_start).count();
                    const auto history = annealer::solve_instrumented(prepared, config);
                    if (!gset::same_result(result, normal) || !gset::same_result(result, history.result))
                        throw std::runtime_error("Instrumentation changed solver result");
                    gset::verify(graph, normal);
                    gset::verify(graph, history.result);
                    std::cout << "Normal seconds: " << normal_seconds << " History seconds: " << history.seconds
                              << " Profile seconds: " << measured.seconds
                              << " History overhead seconds: " << history.seconds-normal_seconds
                              << " Profile overhead seconds: " << measured.seconds-normal_seconds << '\n';
                    double accounted = 0;
                    for (std::size_t j = 0; j < measured.profile.size(); ++j) {
                        const auto& item = measured.profile[j];
                        accounted += item.seconds;
                        std::cout << annealer::profile_names[j] << " seconds: " << item.seconds << " Calls: " << item.calls << '\n';
                    }
                    std::cout << "Initialization/finalization/logging/clock residual seconds: " << measured.seconds-accounted << '\n';
                }
                gset::verify(graph, result);
                if (target) {
                    const auto hit = gset::target_hit(measured.history, *target);
                    if (hit) std::cout << "TTS seconds: " << hit->seconds << " Updates: " << hit->attempted_moves << '\n';
                    else std::cout << "Target not reached\n";
                }
                if (csv.is_open()) gset::history_csv(csv, config, measured.history);
                const double seconds = measured.seconds;
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
        } });
    } catch (const std::exception& e) {
        if (!session.rank) std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
}
