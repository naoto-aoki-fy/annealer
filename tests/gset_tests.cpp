#include "../examples/gset_maxcut.hpp"
#include <annealer/prepared.hpp>
#include <iostream>
#include <sstream>

void require(bool condition) {
    if (!condition) throw std::runtime_error("G-set correctness check failed");
}
void enumerate(const std::string& text, double maximum) {
    std::istringstream input(text);
    const auto graph = gset::read(input);
    const auto model = gset::make_model(graph);
    const auto prepared = annealer::prepare(model);
    double best = -std::numeric_limits<double>::infinity();
    for (std::size_t mask = 0; mask < (std::size_t{1} << graph.vertices); ++mask) {
        annealer::Assignment x(graph.vertices);
        for (std::size_t i = 0; i < x.size(); ++i) x[i] = (mask >> i) & 1;
        const double cut = gset::cut_weight(graph, x);
        require(std::abs(annealer::full_energy(model, x) + cut) < 1e-12);
        require(std::abs(annealer::prepared_energy(prepared, x) + cut) < 1e-12);
        best = std::max(best, cut);
    }
    require(best == maximum);
}
int main() {
    try {
        enumerate("2 1\n1 2 3\n", 3);
        enumerate("3 3\n1 2 1\n2 3 1\n1 3 1\n", 2);
        enumerate("3 3\n1 2 2\n2 3 -3\n1 3 1\n", 3);
        enumerate("3 1\n1 3 0.5\n", 0.5);
        enumerate("2 0\n", 0);
        std::istringstream indices("3 1\n1 3 7\n");
        const auto graph = gset::read(indices);
        require(graph.edges[0].u == 0 && graph.edges[0].v == 2);
        require(gset::cut_weight(graph, {0, 0, 1}) == 7);
        for (const auto* invalid : {"", "0 0", "-2 0", "2 -1", "2 1", "2 1\n0 2 1",
             "2 1\n1 3 1", "2 1\n1 1 1", "2 1\n1 2 nan", "2 1\n1 2 inf",
             "2 1\n1 2 1e999", "2 0\nextra", "2x 0", "2 1\n1 2 1\n1 2 2"}) {
            bool rejected = false;
            try { std::istringstream input(invalid); (void)gset::read(input); }
            catch (const std::invalid_argument&) { rejected = true; }
            require(rejected);
        }
        std::cout << "G-set exhaustive formulation and parser tests passed\n";
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
