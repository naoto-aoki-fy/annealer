#include "../examples/maxsat.hpp"
#include <annealer/solver.hpp>
#include <sstream>
#include <iostream>

void require(bool ok) { if (!ok) throw std::runtime_error("MaxSAT check failed"); }
maxsat::Instance parse(const std::string& text) { std::istringstream in(text); return maxsat::read(in); }
template<class F> void rejects(F f) {
    bool rejected = false;
    try { f(); } catch (const std::invalid_argument&) { rejected = true; }
    require(rejected);
}
void enumerate(const std::string& text) {
    const auto instance = parse(text);
    const auto model = maxsat::make_model(instance);
    require(model.terms().size() == instance.clauses.size());
    for (std::size_t i = 0; i < instance.clauses.size(); ++i)
        require(model.terms()[i].factors().size() == instance.clauses[i].literals.size());
    for (auto budget : {std::size_t{0}, std::size_t{1048576}}) {
        const auto prepared = annealer::prepare(model, {.lookup_bytes = budget});
        for (std::size_t state = 0; state < (std::size_t{1} << instance.variables); ++state) {
            annealer::Assignment x(instance.variables);
            for (std::size_t i = 0; i < x.size(); ++i) x[i] = (state >> i) & 1;
            const auto objective = maxsat::evaluate(instance, x);
            const auto expected = maxsat::energy(instance, objective);
            require(annealer::full_energy(model, x) == expected);
            require(annealer::prepared_energy(prepared, x) == expected);
        }
        annealer::SolverConfig config;
        config.seed = 42; config.temperature_levels = 5; config.sweeps_per_level = 3;
        const auto a = annealer::solve(prepared, config);
        const auto b = annealer::solve(prepared, config);
        require(a.best_assignment == b.best_assignment && a.best_energy == b.best_energy);
        require(a.best_energy == maxsat::energy(instance, maxsat::evaluate(instance, a.best_assignment)));
        require(a.final_energy == maxsat::energy(instance, maxsat::evaluate(instance, a.final_assignment)));
    }
}
int main() {
    try {
        enumerate("c mixed\np wcnf 5 4 100\n100 1 -2 3 0\n3 -1 2 -3 4 -5 0\n7 1 0\n2 -1 0\n");
        enumerate("h 1 -2 3 0\n3 -1 2 -3 4 -5 0\n7 1 0\n2 -1 0\n");
        enumerate("p wcnf 5 2\nh 1 -2 3 0\n9 1 2 3 -4 -5 0\n");
        enumerate("p wcnf 2 3 9\n10 1 1 0\n2 1 -1 0\n3 0\n");
        enumerate("p wcnf 1 2\nh 0\nh 1 0\n");
        enumerate("p wcnf 1 0\n");
        enumerate("p wcnf 2 2\n2 1\nc continuation\n-2 0 3 -1 0\n");
        const auto known = parse("p wcnf 1 3 10\n10 1 0\n3 -1 0\n7 0\n");
        auto objective = maxsat::evaluate(known, {1});
        require(objective.soft_cost == 10 && objective.hard_violations == 0 && known.penalty == 11);
        objective = maxsat::evaluate(known, {0});
        require(objective.soft_cost == 7 && objective.hard_violations == 1);
        for (const auto* text : {"", "p cnf 1 0", "p wcnf -1 0", "p wcnf 1 -1", "p wcnf 1 0 0",
            "p wcnf 1 0 extra", "p wcnf 1 0\np wcnf 1 0", "p wcnf 1 1\n1 2 0",
            "p wcnf 1 1\n1 -2 0", "p wcnf 1 1\n1 -0 0", "p wcnf 1 1\n1 1", "p wcnf 1 2\n1 1 0",
            "p wcnf 1 0\n1 1 0", "0 1 0", "-2 1 0", "1.5 1 0", "h +1 0", "h 1x 0",
            "p wcnf 1 1 10\nh 1 0", "18446744073709551616 1 0", "9007199254740993 1 0",
            "9007199254740992 1 0\n1 -1 0", "9007199254740992 1 0\nh 1 0",
            "4503599627370496 1 0\nh 1 0", "h 18446744073709551616 0", "h 1 0\np wcnf 1 1"})
            rejects([&] { (void)parse(text); });
        rejects([&] { (void)maxsat::evaluate(known, {}); });
        rejects([&] { (void)maxsat::evaluate(known, {2}); });
        std::cout << "MaxSAT exhaustive energy, parser, verifier and replay tests passed\n";
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
