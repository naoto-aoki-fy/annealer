#include "../examples/dimacs_coloring.hpp"
#include <annealer/solver.hpp>
#include <iostream>
#include <sstream>

void require(bool ok) { if (!ok) throw std::runtime_error("Coloring check failed"); }
template<class F> void rejects(F action) {
    bool rejected = false;
    try { action(); } catch (const std::invalid_argument&) { rejected = true; }
    require(rejected);
}
void enumerate(const std::string& text, std::size_t colors, std::size_t minimum) {
    std::istringstream input(text);
    const auto graph = coloring::read(input);
    const auto model = coloring::make_model(graph, colors);
    for (auto budget : {std::size_t{0}, std::size_t{1048576}}) {
        const auto prepared = annealer::prepare(model, {.lookup_bytes = budget});
        std::size_t states = 1, best = graph.edges.size();
        for (std::size_t i = 0; i < graph.vertices; ++i) states *= colors;
        for (std::size_t state = 0; state < states; ++state) {
            annealer::Assignment x(graph.vertices);
            auto digits = state;
            for (auto& value : x) { value = static_cast<annealer::Value>(digits % colors); digits /= colors; }
            const auto count = coloring::conflicts(graph, colors, x);
            require(annealer::full_energy(model, x) == static_cast<double>(count));
            require(annealer::prepared_energy(prepared, x) == static_cast<double>(count));
            best = std::min(best, count);
        }
        require(best == minimum);
        annealer::SolverConfig config;
        config.seed = 123;
        config.temperature_levels = 10;
        config.sweeps_per_level = 5;
        const auto a = annealer::solve(prepared, config);
        const auto b = annealer::solve(prepared, config);
        require(a.best_assignment == b.best_assignment && a.final_assignment == b.final_assignment);
        require(a.best_energy == b.best_energy && a.final_energy == b.final_energy && a.seed == b.seed);
        require(a.attempted_moves == b.attempted_moves && a.accepted_moves == b.accepted_moves);
        require(static_cast<double>(coloring::conflicts(graph, colors, a.best_assignment)) == a.best_energy);
        require(static_cast<double>(coloring::conflicts(graph, colors, a.final_assignment)) == a.final_energy);
    }
}
int main() {
    try {
        const std::string triangle = "c triangle\np edge 3 3\ne 1 2\ne 2 3\ne 3 1\n";
        enumerate(triangle, 1, 3);
        enumerate(triangle, 2, 1);
        enumerate(triangle, 3, 0);
        enumerate("p edge 4 3\ne 1 2\ne 2 3\ne 3 4\n", 2, 0);
        enumerate("p edge 3 0\n", 1, 0);
        enumerate("\n c comment\r\np edge 3 1\r\ne 3 1\r\nc end\n", 2, 0);
        for (const auto* text : {"", "p edge 0 0", "p edge -1 0", "p edge 2 -1", "p edge 2x 0",
             "p edges 2 0", "p edge 2 0 extra", "p edge 2 0\np edge 2 0", "e 1 2\np edge 2 1",
             "p edge 2 1", "p edge 2 1\ne 0 2", "p edge 2 1\ne 1 3", "p edge 2 1\ne 1 1",
             "p edge 2 2\ne 1 2\ne 2 1", "p edge 2 0\ne 1 2", "p edge 2 1\ne 1",
             "p edge 2 1\ne 1 2 extra", "p edge 2 0\nx", "123\np edge 2 0"})
            rejects([&] { std::istringstream input(text); (void)coloring::read(input); });
        const coloring::Graph graph{2, {{0, 1}}};
        rejects([&] { (void)coloring::make_model(graph, 0); });
        rejects([&] { coloring::validate_colors(std::numeric_limits<std::size_t>::max()); });
        rejects([&] { (void)coloring::conflicts(graph, 2, {0}); });
        rejects([&] { (void)coloring::conflicts(graph, 2, {-1, 0}); });
        rejects([&] { (void)coloring::conflicts(graph, 2, {0, 2}); });
        std::cout << "Coloring exhaustive energy, parser, verifier and replay tests passed\n";
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
