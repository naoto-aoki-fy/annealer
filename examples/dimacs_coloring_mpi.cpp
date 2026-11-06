#include "dimacs_coloring.hpp"
#include <annealer/mpi.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>

int main(int argc, char** argv) {
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 1;
    int rank = 0, size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int status = 0;
    try {
        if (argc < 2 || (argc == 2 && std::string(argv[1]) == "--help")) {
            if (!rank)
                std::cout << "Usage: dimacs_coloring_mpi_example FILE [--colors K] [--seed N] [--trials N]\n"
                             "  [--initial-temperature T] [--final-temperature T] [--levels N] [--sweeps N]\n"
                             "Defaults: colors=5 seed=0 trials=10 temperatures=30,0.05 levels=200 sweeps=100\n"
                             "Trials is the global number of independent starts.\n";
            MPI_Finalize();
            return argc < 2 ? 1 : 0;
        }
        annealer::SolverConfig config;
        std::size_t colors = 5, trials = 10;
        coloring::Graph graph;
        std::optional<annealer::PreparedModel> prepared;
        int failed = 0, any_failed = 0;
        // Agree on local input/preparation failures before entering solver collectives.
        try {
            config.initial_temperature = 30;
            config.final_temperature = 0.05;
            config.temperature_levels = 200;
            config.sweeps_per_level = 100;
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
            graph = coloring::read(file);
            prepared.emplace(annealer::prepare(coloring::make_model(graph, colors),
                                               {.lookup_bytes = 1024 * 1024}));
        } catch (const std::exception& error) {
            std::cerr << "rank " << rank << ": " << error.what() << '\n';
            failed = 1;
        } catch (...) { failed = 1; }
        MPI_Allreduce(&failed, &any_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if (any_failed) throw std::runtime_error("Input or model preparation failed on one or more ranks");

        const auto result = annealer::solve_multistart(MPI_COMM_WORLD, *prepared, trials, config);
        std::size_t count = 0;
        failed = 0;
        try {
            count = coloring::conflicts(graph, colors, result.best_assignment);
            if (static_cast<double>(count) != result.best_energy)
                throw std::runtime_error("Independent conflict/energy verification failed");
            for (const auto& run : result.local_runs) {
                if (static_cast<double>(coloring::conflicts(graph, colors, run.result.best_assignment)) != run.result.best_energy ||
                    static_cast<double>(coloring::conflicts(graph, colors, run.result.final_assignment)) != run.result.final_energy)
                    throw std::runtime_error("Independent local conflict/energy verification failed");
            }
        } catch (...) { failed = 1; }
        MPI_Allreduce(&failed, &any_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if (any_failed) throw std::runtime_error("Independent conflict/energy verification failed on one or more ranks");
        if (!rank) {
            std::cout << std::setprecision(17) << "Vertices: " << graph.vertices
                      << " Edges: " << graph.edges.size() << " Colors: " << colors
                      << " Processes: " << size << "\nInitial temperature: " << config.initial_temperature
                      << " Final temperature: " << config.final_temperature << " Levels: " << config.temperature_levels
                      << " Sweeps: " << config.sweeps_per_level << " Trials: " << trials
                      << "\nBest seed: " << result.seed << " Best energy: " << result.best_energy
                      << " Conflicting edges: " << count
                      << " Status: " << (count == 0 ? "valid K-coloring" : "assignment contains conflicts") << '\n';
        }
    } catch (const std::exception& error) {
        if (!rank) std::cerr << "Error: " << error.what() << '\n';
        status = 1;
    } catch (...) { status = 1; }
    int global_status = 0;
    MPI_Allreduce(&status, &global_status, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Finalize();
    return global_status;
}
