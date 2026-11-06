#include "../examples/jobshop_model.hpp"
#include <annealer/solver.hpp>
#include <iostream>
#include <sstream>
void require(bool ok) { if (!ok) throw std::runtime_error("Job-shop model check failed"); }
jobshop::Instance parse(const std::string& text) { std::istringstream in(text); return jobshop::read(in); }
void exhaustive(const std::string& text, double optimum) {
    const auto instance = parse(text); const auto plan = jobshop::estimate(instance);
    const auto model = jobshop::make_model(instance);
    require(model.terms().size() == plan.terms && model.domains().size() == plan.variables);
    for (auto budget : {std::size_t{0}, std::size_t{32}, std::size_t{1048576}}) {
        const auto prepared = annealer::prepare(model, {.lookup_bytes = budget});
        std::size_t states = 1; for (const auto& d : model.domains()) states *= d.values().size();
        double best = 1e100;
        for (std::size_t state = 0; state < states; ++state) {
            auto digits = state; annealer::Assignment x;
            for (const auto& d : model.domains()) { x.push_back(d.values()[digits % d.values().size()]); digits /= d.values().size(); }
            const auto direct = jobshop::verify_model(instance, x);
            require(direct.energy == annealer::full_energy(model, x));
            require(direct.energy == annealer::prepared_energy(prepared, x));
            if (!direct.feasible()) require(direct.energy > plan.horizon);
            best = std::min(best, direct.energy);
        }
        require(best == optimum);
        annealer::SolverConfig config;
        config.initial_assignment = jobshop::initial_assignment(instance); config.seed = 42;
        config.temperature_levels = 3; config.sweeps_per_level = 2;
        const auto a = annealer::solve(prepared, config), b = annealer::solve(prepared, config);
        require(a.best_assignment == b.best_assignment && a.final_assignment == b.final_assignment && a.best_energy == b.best_energy);
        require(jobshop::verify_model(instance, a.best_assignment).feasible());
        require(a.best_energy == jobshop::verify_model(instance, a.best_assignment).energy);
        require(a.final_energy == jobshop::verify_model(instance, a.final_assignment).energy);
    }
}
int main() {
    try {
        exhaustive("1 1 0 2", 2);
        exhaustive("2 1 0 1 0 2", 3);
        exhaustive("1 2 0 1 1 1", 2);
        exhaustive("2 2 0 1 1 1 1 1 0 1", 2);
        for (const auto [jobs, machines] : {std::pair{6,6}, std::pair{10,5}}) {
            std::ostringstream text; text << jobs << ' ' << machines << '\n';
            for (int j = 0; j < jobs; ++j) for (int s = 0; s < machines; ++s) text << (s+j)%machines << " 1 ";
            const auto instance = parse(text.str());
            const auto model = jobshop::make_model(instance);
            annealer::SolverConfig config; config.initial_assignment = jobshop::initial_assignment(instance);
            config.temperature_levels = 1; config.sweeps_per_level = 1;
            const auto result = annealer::solve(annealer::prepare(model, {.lookup_bytes = 0}), config);
            require(jobshop::verify_model(instance, result.best_assignment).feasible());
        }
        bool rejected = false;
        try { (void)jobshop::make_model(parse("1 1 0 1"), {.memory_bytes = 0}); } catch (const std::invalid_argument&) { rejected = true; }
        require(rejected);
        rejected = false;
        try { (void)jobshop::estimate(parse("1 1 0 9007199254740992")); } catch (const std::invalid_argument&) { rejected = true; }
        require(rejected);
        std::cout << "Job-shop exhaustive energy, penalty, budget and replay tests passed\n";
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
