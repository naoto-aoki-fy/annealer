#include "../examples/tsplib.hpp"
#include <annealer/solver.hpp>
#include <sstream>
#include <iostream>
void require(bool ok) { if (!ok) throw std::runtime_error("TSP check failed"); }
tsp::Instance parse(const std::string& text) { std::istringstream in(text); return tsp::read(in); }
const std::string header = "NAME: tiny\nTYPE : TSP\nDIMENSION: 4\nEDGE_WEIGHT_TYPE: EUC_2D\nNODE_COORD_SECTION\n";
template<class F> void rejects(F f) { bool bad = false; try { f(); } catch (const std::exception&) { bad = true; } require(bad); }
int main() {
    try {
        require(tsp::distance({0, 0}, {0, 0.5}) == 1);
        require(tsp::distance({0, 0}, {0, 0.49}) == 0);
        require(tsp::distance({0, 0}, {1, 1}) == 1);
        require(tsp::distance({-1, -2}, {2, 2}) == 5);
        const auto instance = parse(header + "3 1 1\n1 0 0\n4 0 1\n2 1 0\nEOF\n");
        require(instance.points[0].x == 0 && instance.points[1].x == 1);
        const auto model = tsp::make_model(instance);
        require(model.domains().size() == 4 && model.terms().size() == 40);
        for (auto budget : {std::size_t{0}, std::size_t{32}, std::size_t{1048576}}) {
            const auto prepared = annealer::prepare(model, {.lookup_bytes = budget});
            tsp::Cost minimum = tsp::limit;
            for (unsigned state = 0; state < 256; ++state) {
                annealer::Assignment x(4); auto digits = state;
                for (auto& v : x) { v = digits % 4; digits /= 4; }
                const auto objective = tsp::evaluate(instance, x);
                require(tsp::energy(instance, objective) == annealer::full_energy(model, x));
                require(tsp::energy(instance, objective) == annealer::prepared_energy(prepared, x));
                if (objective.feasible) { require(objective.cost == 4); minimum = std::min(minimum, objective.cost); }
                else require(tsp::energy(instance, objective) >= instance.penalty);
            }
            require(minimum == 4);
            annealer::SolverConfig config;
            config.swap_proposals = true; config.seed = 42; config.initial_assignment = tsp::permutation(4, 42);
            config.temperature_levels = 4; config.sweeps_per_level = 3;
            const auto a = annealer::solve(prepared, config), b = annealer::solve(prepared, config);
            const auto reference = annealer::solve_reference(model, config);
            require(a.best_assignment == b.best_assignment && a.final_assignment == b.final_assignment);
            require(a.final_assignment == reference.final_assignment && a.best_energy == reference.best_energy);
            require(tsp::evaluate(instance, a.best_assignment).feasible && tsp::evaluate(instance, a.final_assignment).feasible);
        }
        const auto invalid = tsp::evaluate(instance, {0, 0, 0, 0});
        require(!invalid.feasible && invalid.collisions == 6 && invalid.cost == 0 && invalid.penalty_cost == 30);
        rejects([&] { (void)tsp::evaluate(instance, {0, 1}); });
        rejects([&] { (void)tsp::evaluate(instance, {-1, 1, 2, 3}); });
        rejects([&] { (void)tsp::evaluate(instance, {0, 1, 2, 4}); });
        require(parse(header + "1 0 0\n2 1 0\n3 1 1\n4 0 1\n").size == 4);
        for (const auto* records : {"", "1 0 0", "0 0 0", "5 0 0", "1 nan 0", "1 0 inf", "1 0x 0",
             "1 0 0 extra", "1 0 0\n1 1 1", "EOF\n1 0 0", "1 1e308 0\n2 -1e308 0\n3 0 0\n4 0 1"})
            rejects([&] { (void)parse(header + records); });
        for (const auto* text : {"", "TYPE: ATSP\nDIMENSION: 2\nEDGE_WEIGHT_TYPE: EUC_2D\nNODE_COORD_SECTION",
             "TYPE: TSP\nDIMENSION: 2\nEDGE_WEIGHT_TYPE: GEO\nNODE_COORD_SECTION",
             "TYPE: TSP\nDIMENSION: -1\nEDGE_WEIGHT_TYPE: EUC_2D\nNODE_COORD_SECTION",
             "TYPE: TSP\nDIMENSION: 1\nEDGE_WEIGHT_TYPE: EUC_2D\nNODE_COORD_SECTION",
             "TYPE: TSP\nTYPE: TSP", "EDGE_WEIGHT_SECTION", "DIMENSION: 18446744073709551616"})
            rejects([&] { (void)parse(text); });
        const auto two = parse("TYPE TSP\nDIMENSION 2\nEDGE_WEIGHT_TYPE EUC_2D\nNODE_COORD_SECTION\n1 0 0\n2 3 4");
        require(tsp::evaluate(two, {0, 1}).cost == 10);
        require(annealer::full_energy(tsp::make_model(two), {0, 1}) == 10);
        std::cout << "TSPLIB parser, rounding, exhaustive energy and replay tests passed\n";
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
