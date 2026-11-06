#include "../examples/qaplib.hpp"
#include <annealer/solver.hpp>
#include <iostream>
#include <sstream>
void require(bool ok) { if (!ok) throw std::runtime_error("QAP check failed"); }
qap::Instance parse(const std::string& text) { std::istringstream in(text); return qap::read(in); }
template<class F> void rejects(F f) { bool rejected = false; try { f(); } catch (const std::invalid_argument&) { rejected = true; } require(rejected); }
void enumerate(const std::string& text, qap::Cost optimum) {
    const auto instance = parse(text);
    const auto model = qap::make_model(instance);
    const auto n = instance.size;
    for (auto budget : {std::size_t{0}, std::size_t{32}, std::size_t{1048576}}) {
        const auto prepared = annealer::prepare(model, {.lookup_bytes = budget});
        std::size_t states = 1; for (std::size_t i = 0; i < n; ++i) states *= n;
        qap::Cost minimum = qap::limit;
        for (std::size_t state = 0; state < states; ++state) {
            annealer::Assignment x(n); auto digits = state;
            for (auto& value : x) { value = digits % n; digits /= n; }
            const auto objective = qap::evaluate(instance, x);
            const auto energy = qap::energy(instance, objective);
            require(energy == annealer::full_energy(model, x));
            require(energy == annealer::prepared_energy(prepared, x));
            if (objective.feasible) minimum = std::min(minimum, objective.cost);
            else require(energy >= instance.penalty);
        }
        require(minimum == optimum);
        auto x = qap::permutation(n, 42); std::sort(x.begin(), x.end());
        do { require(qap::evaluate(instance, x).feasible); } while (std::next_permutation(x.begin(), x.end()));
        annealer::SolverConfig config;
        config.seed = 42; config.swap_proposals = true; config.initial_assignment = qap::permutation(n, config.seed);
        config.temperature_levels = 5; config.sweeps_per_level = 4;
        const auto a = annealer::solve(prepared, config), b = annealer::solve(prepared, config);
        const auto reference = annealer::solve_reference(model, config);
        require(a.best_assignment == b.best_assignment && a.final_assignment == b.final_assignment);
        require(a.best_assignment == reference.best_assignment && a.final_assignment == reference.final_assignment);
        require(qap::evaluate(instance, a.best_assignment).feasible && qap::evaluate(instance, a.final_assignment).feasible);
        require(a.best_energy == qap::energy(instance, qap::evaluate(instance, a.best_assignment)));
        require(a.final_energy == qap::energy(instance, qap::evaluate(instance, a.final_assignment)));
    }
}
int main() {
    try {
        enumerate("2 0 1 2 0 0 3 4 0", 10); // asymmetric B, directed objective
        enumerate("2 1 2 3 4 5 6 6 7", 57); // diagonal and symmetric B
        enumerate("3 0 1 0 1 0 2 0 2 0 0 3 4 3 0 5 4 5 0", 20);
        enumerate("1 3 4", 12);
        enumerate("2 0 0 0 0 0 0 0 0", 0);
        const auto known = parse("2 0 1 2 0 0 3 4 0");
        const auto collision = qap::evaluate(known, {0, 0});
        require(!collision.feasible && collision.collisions == 1 && collision.penalty_cost == 13 && collision.cost == 0);
        for (const auto* text : {"", "0", "-1", "2x", "2 0 0", "1 1", "1 1 2 extra", "1 -1 2",
            "1 1.5 2", "1 nan 2", "18446744073709551616", "9007199254740992",
            "1 9007199254740992 1", "1 9007199254740991 2", "2 1 1 1 1 9007199254740991 0 0 0"})
            rejects([&] { (void)parse(text); });
        rejects([&] { (void)qap::evaluate(known, {}); });
        rejects([&] { (void)qap::evaluate(known, {-1, 1}); });
        rejects([&] { (void)qap::evaluate(known, {0, 2}); });
        require(qap::permutation(20, 19) == qap::permutation(20, 19));
        std::cout << "QAP exhaustive energy, permutation, parser and replay checks passed\n";
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
