#include <annealer/solver.hpp>
#include <iostream>
#include <sstream>

using namespace annealer;
void require(bool ok) { if (!ok) throw std::runtime_error("execution state regression"); }
bool same(const SolverResult& a, const SolverResult& b) {
    return a.best_assignment == b.best_assignment && a.final_assignment == b.final_assignment &&
        a.best_energy == b.best_energy && a.final_energy == b.final_energy &&
        a.attempted_moves == b.attempted_moves && a.accepted_moves == b.accepted_moves && a.seed == b.seed;
}
void history_same(const InstrumentedResult& a, const InstrumentedResult& b) {
    require(a.history.size() == b.history.size());
    for (std::size_t i = 0; i < a.history.size(); ++i) {
        require(a.history[i].best_energy == b.history[i].best_energy);
        require(a.history[i].attempted_moves == b.history[i].attempted_moves);
        require(a.history[i].accepted_moves == b.history[i].accepted_moves);
    }
    for (std::size_t i = 0; i < 6; ++i) require(a.profile[i].calls == b.profile[i].calls);
}
int main() {
 try {
    ModelBuilder b;
    for (int i = 0; i < 4; ++i) b.add_variable(Domain{-2, 0, 3});
    b.add_term(1, {{0, [](Value x) { return double(x); }}, {1, [](Value x) { return double(x); }}});
    b.add_term(-2, {{2, [](Value x) { return double(x); }}});
    const auto p = prepare(b.build());
    // Captured from the serial solver before extraction, including history and profiling calls.
    const char* golden[] = {
        "-12 -12 1800 783 -2 3 3 -2 -2 3 3 -2 h 0 0 -6 h 19 18 -12 p 1800 p 1800 p 1800 p 783 p 3 p 1",
        "-12 -12 1800 757 3 -2 3 -2 -2 3 3 -2 h 0 0 -6 h 13 11 -12 p 1800 p 1800 p 1800 p 757 p 2 p 1",
        "-6 -6 1800 1059 -2 0 3 -2 0 -2 3 -2 h 0 0 -6 p 1800 p 1800 p 1800 p 1059 p 4 p 0",
        "-6 -6 1800 1071 -2 0 3 -2 -2 3 0 -2 h 0 0 -6 p 1800 p 1800 p 1800 p 1071 p 4 p 0"};
    std::size_t fixture = 0;
    for (bool swaps : {false, true}) for (unsigned seed : {0, 42}) {
        SolverConfig c;
        c.seed = seed; c.swap_proposals = swaps; c.initial_assignment = Assignment{-2, 0, 3, -2};
        c.temperature_levels = 3; c.sweeps_per_level = 150;
        auto serial = solve_instrumented(p, c, {.profile = true});
        const auto& r = serial.result;
        std::ostringstream text;
        text << r.best_energy << ' ' << r.final_energy << ' ' << r.attempted_moves << ' ' << r.accepted_moves;
        for (auto x : r.best_assignment) text << ' ' << x;
        for (auto x : r.final_assignment) text << ' ' << x;
        for (auto h : serial.history) text << " h " << h.attempted_moves << ' ' << h.accepted_moves << ' ' << h.best_energy;
        for (auto v : serial.profile) text << " p " << v.calls;
        require(text.str() == golden[fixture++]);
        require(same(r, solve_reference(p.source(), c)));
        InstrumentedResult split;
        detail::ExecutionState<true, true> chain(p, c, &split);
        for (std::size_t level = 0; level < 3; ++level) {
            const auto t = geometric_temperature(c.initial_temperature, c.final_temperature, 3, level);
            chain.advance_attempts(0, t);
            // One-attempt blocks cross every accepted-move rebase boundary.
            for (int i = 0; i < 600; ++i) chain.advance_attempts(1, t);
        }
        split.result = chain.finish();
        require(same(r, split.result)); history_same(serial, split);
        detail::ExecutionState<> sweeps(p, c), imported(p, c);
        for (std::size_t level = 0; level < 3; ++level) {
            const auto t = geometric_temperature(c.initial_temperature, c.final_temperature, 3, level);
            sweeps.advance_sweeps(150, t);
            imported.advance_attempts(257, t);
            imported.install_assignment(imported.assignment());
            imported.advance_attempts(343, t);
        }
        require(same(r, sweeps.finish()) && same(r, imported.finish()));
        auto before = imported.snapshot();
        for (const Assignment& bad : {Assignment{}, Assignment{999, 0, 3, -2}}) {
            bool caught = false;
            try { imported.install_assignment(bad); } catch (const std::invalid_argument&) { caught = true; }
            require(caught && same(before, imported.snapshot()));
        }
    }
    SolverConfig c;
    c.initial_assignment = Assignment{-2, 3, 3, -2}; // energy -12
    detail::ExecutionState<> chain(p, c);
    chain.install_assignment(Assignment{3, 3, -2, 0}); // worse, changes all relevant caches
    auto installed = chain.snapshot();
    require(installed.best_energy == -12 && installed.final_energy == 13);
    require(installed.attempted_moves == 0 && installed.accepted_moves == 0);
    c.initial_assignment = installed.final_assignment;
    c.temperature_levels = 1; c.sweeps_per_level = 300;
    chain.advance_sweeps(300, c.initial_temperature);
    const auto oracle = solve_reference(p.source(), c);
    const auto continued = chain.finish();
    require(continued.final_assignment == oracle.final_assignment && continued.final_energy == oracle.final_energy);
    require(continued.accepted_moves == oracle.accepted_moves && continued.best_energy <= -12);
    c.initial_assignment = Assignment{3, 3, -2, 0};
    detail::ExecutionState<> improving(p, c);
    improving.install_assignment(Assignment{-2, 3, 3, -2});
    require(improving.snapshot().best_energy == -12);
    for (auto budget : {0u, 1u}) {
        c.temperature_levels = budget; c.sweeps_per_level = 0;
        detail::ExecutionState<> zero(p, c);
        zero.advance_attempts(0, 1); zero.advance_sweeps(0, 1);
        require(same(zero.finish(), solve(p, c)));
    }
    const auto empty = prepare(ModelBuilder{}.build());
    detail::ExecutionState<> zero(empty, SolverConfig{});
    zero.advance_sweeps(100, 1);
    require(zero.finish().attempted_moves == 0);
    std::cout << "Execution state regression tests passed\n";
 } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
