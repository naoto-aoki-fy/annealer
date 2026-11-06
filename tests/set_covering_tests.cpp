#include "../examples/set_covering.hpp"
#include <annealer/solver.hpp>
#include <iostream>
#include <sstream>

void require(bool ok) { if (!ok) throw std::runtime_error("Set covering check failed"); }
covering::Instance parse(const std::string& text) { std::istringstream in(text); return covering::read(in); }
template<class F> void rejects(F f) {
    bool rejected = false;
    try { f(); } catch (const std::invalid_argument&) { rejected = true; }
    require(rejected);
}
void enumerate(const std::string& text, covering::Cost optimum) {
    const auto instance = parse(text);
    const auto model = covering::make_model(instance);
    require(model.domains().size() == instance.costs.size());
    require(model.terms().size() == instance.costs.size() + instance.rows.size());
    for (std::size_t i = 0; i < instance.rows.size(); ++i)
        require(model.terms()[instance.costs.size() + i].factors().size() == instance.rows[i].size());
    for (auto budget : {std::size_t{0}, std::size_t{16}, std::size_t{1048576}}) {
        const auto prepared = annealer::prepare(model, {.lookup_bytes = budget});
        covering::Cost best = instance.cost_sum;
        double best_energy = std::numeric_limits<double>::infinity();
        for (std::size_t state = 0; state < (std::size_t{1} << instance.costs.size()); ++state) {
            annealer::Assignment x(instance.costs.size());
            for (std::size_t j = 0; j < x.size(); ++j) x[j] = (state >> j) & 1;
            const auto objective = covering::evaluate(instance, x);
            const auto energy = covering::energy(instance, objective);
            require(annealer::full_energy(model, x) == energy);
            require(annealer::prepared_energy(prepared, x) == energy);
            if (!objective.uncovered_rows) best = std::min(best, objective.selected_cost);
            else require(energy > static_cast<double>(instance.cost_sum));
            best_energy = std::min(best_energy, energy);
        }
        require(best == optimum && best_energy == static_cast<double>(optimum));
        annealer::SolverConfig config;
        config.seed = 42; config.temperature_levels = 5; config.sweeps_per_level = 3;
        const auto result = annealer::solve(prepared, config);
        require(result.best_energy == covering::energy(instance, covering::evaluate(instance, result.best_assignment)));
        require(result.final_energy == covering::energy(instance, covering::evaluate(instance, result.final_assignment)));
    }
}
int main() {
    try {
        enumerate("3 5\n2 7 0 4 3\n3 1 2 3\n5 1 2 3 4 5\n2 4 5\n", 3);
        enumerate("2 3 5 2 8 2 1 2 2 2 3", 2);
        enumerate("1 1 0 1 1", 0);
        enumerate("2 2 3 4 1 1 1 2", 7);
        const auto instance = parse("2 3 5 2 8 2 1 2 2 2 3");
        auto value = covering::evaluate(instance, {1, 0, 0});
        require(value.selected_cost == 5 && value.uncovered_rows == 1 && instance.penalty == 16);
        for (const auto* text : {"", "0 2", "2 0", "-1 2", "1x 2", "1 2 3", "1 1 -1 1 1",
            "1 1 1.5 1 1", "1 1 nan 1 1", "1 1 1 0", "1 1 1 2 1 1",
            "1 2 1 2 2 1 1", "1 2 1 2 1 0", "1 2 1 2 1 3", "1 2 1 2 1 -1",
            "1 2 1 2 1", "1 1 1 1 1 extra", "2 1 1 1 1", "18446744073709551616 1",
            "1 1 18446744073709551616", "1 1 9007199254740992 1 1",
            "1 1 4503599627370496 1 1", "1 2 9007199254740991 1 1 1"})
            rejects([&] { (void)parse(text); });
        rejects([&] { (void)covering::evaluate(instance, {}); });
        rejects([&] { (void)covering::evaluate(instance, {0, 2, 0}); });
        std::cout << "Set covering exhaustive energy, optimum, parser and verifier tests passed\n";
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
