#include <annealer/solver.hpp>
#include <iostream>

void require(bool ok) { if (!ok) throw std::runtime_error("Swap check failed"); }
template<class F> void rejects(F f) {
    bool rejected = false;
    try { f(); } catch (const std::invalid_argument&) { rejected = true; }
    require(rejected);
}
bool same(const annealer::SolverResult& a, const annealer::SolverResult& b) {
    return a.best_assignment == b.best_assignment && a.final_assignment == b.final_assignment &&
        a.best_energy == b.best_energy && a.final_energy == b.final_energy &&
        a.attempted_moves == b.attempted_moves && a.accepted_moves == b.accepted_moves && a.seed == b.seed;
}
int main() {
    try {
        std::mt19937_64 random(17);
        for (int trial = 0; trial < 12; ++trial) {
            annealer::ModelBuilder builder;
            for (int i = 0; i < 4; ++i) builder.add_variable(annealer::Domain{-3, 0, 2, 5});
            auto identity = [](annealer::Value x) { return static_cast<double>(x); };
            builder.add_term(3, {{0, identity}, {1, identity}, {0, identity}, {2, identity}});
            for (int t = 0; t < 15; ++t) {
                std::vector<annealer::Factor> factors;
                for (unsigned k = 0, n = random() % 6; k < n; ++k)
                    factors.emplace_back(random() % 4, identity);
                builder.add_term(static_cast<double>(static_cast<int>(random() % 11) - 5), factors);
            }
            const auto model = builder.build();
            for (auto budget : {std::size_t{0}, std::size_t{32}, std::size_t{1048576}}) {
                const auto prepared = annealer::prepare(model, {.lookup_bytes = budget});
                annealer::PreparedState state(prepared, {-3, 0, 2, 5});
                for (int move = 0; move < 300; ++move) {
                    const auto pair = annealer::detail::draw_swap(4, random);
                    require(pair.variable != pair.index && pair.index < 4);
                    const auto before = state.assignment();
                    const auto before_energy = state.energy();
                    auto after = before;
                    std::swap(after[pair.variable], after[pair.index]);
                    const auto delta = state.propose_swap(pair.variable, pair.index);
                    require(delta == annealer::full_energy(model, after) - annealer::full_energy(model, before));
                    require(state.assignment() == before && state.energy() == before_energy);
                    // Alternate acceptance with deliberate rejection (no accept call).
                    if (move % 2) {
                        state.accept();
                        require(state.assignment() == after && state.energy() == annealer::full_energy(model, after));
                    }
                    if (move % 31 == 0) state.rebase();
                }
                annealer::SolverConfig config;
                config.swap_proposals = true; config.seed = 42;
                config.initial_assignment = annealer::Assignment{-3, 0, 2, 5};
                const auto result = annealer::solve(prepared, config);
                require(same(result, annealer::solve(prepared, config)));
                require(same(result, annealer::solve_reference(model, config)));
                require(same(result, annealer::solve_instrumented(prepared, config, {.profile = true}).result));
                auto sorted = result.final_assignment; std::sort(sorted.begin(), sorted.end());
                require(sorted == *config.initial_assignment);
                require(result.attempted_moves == 4000);
                config.temperature_levels = 0;
                require(annealer::solve(prepared, config).attempted_moves == 0);
                require(annealer::solve_reference(model, config).final_assignment == *config.initial_assignment);
                config.temperature_levels = 1; config.sweeps_per_level = 0;
                require(annealer::solve(prepared, config).attempted_moves == 0);
            }
        }
        // Each ordered pair has one equally likely rank, including its reverse.
        std::size_t counts[4][4]{};
        for (std::size_t a = 0; a < 4; ++a)
            for (std::size_t r = 0; r < 3; ++r) ++counts[a][annealer::detail::replacement_index(a, r)];
        for (std::size_t a = 0; a < 4; ++a)
            for (std::size_t b = 0; b < 4; ++b) require(counts[a][b] == (a == b ? 0u : 1u));
        annealer::ModelBuilder invalid;
        invalid.add_variable(annealer::Domain{0, 1}); invalid.add_variable(annealer::Domain{0, 2});
        const auto model = invalid.build();
        annealer::SolverConfig config; config.swap_proposals = true; config.temperature_levels = 0;
        rejects([&] { (void)annealer::solve(model, config); });
        rejects([&] { (void)annealer::solve_reference(model, config); });
        const auto prepared = annealer::prepare(model);
        annealer::PreparedState valid_lifetime(prepared, {0, 0});
        rejects([&] { valid_lifetime.propose_swap(0, 1); });
        rejects([&] { valid_lifetime.propose_swap(0, 0); });
        annealer::ModelBuilder large;
        large.add_variable(annealer::Domain{-1, 1}); large.add_variable(annealer::Domain{-1, 1});
        large.add_term(1e308, {{0, [](annealer::Value x) { return static_cast<double>(x); }}});
        const auto large_prepared = annealer::prepare(large.build());
        annealer::PreparedState large_state(large_prepared, {-1, 1});
        bool overflow = false;
        try { large_state.propose_swap(0, 1); } catch (const std::overflow_error&) { overflow = true; }
        require(overflow && large_state.assignment() == annealer::Assignment({-1, 1}) && large_state.energy() == -1e308);
        config.temperature_levels = 2; config.sweeps_per_level = 3;
        annealer::ModelBuilder identical;
        identical.add_variable(annealer::Domain{7}); identical.add_variable(annealer::Domain{7});
        const auto unchanged = annealer::solve(identical.build(), config);
        require(unchanged.attempted_moves == 12 && unchanged.accepted_moves == 12);
        for (int count : {0, 1}) {
            annealer::ModelBuilder empty;
            if (count) empty.add_variable(annealer::Domain{7});
            require(annealer::solve(empty.build(), config).attempted_moves == 0);
            require(annealer::solve_reference(empty.build(), config).attempted_moves == 0);
        }
        std::cout << "Swap tests passed\n";
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
