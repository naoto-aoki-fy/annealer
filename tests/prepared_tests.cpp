#include <annealer/solver.hpp>

#include <iostream>

namespace {
int failures = 0;
void check(bool ok, const char* message) {
    if (!ok) {
        ++failures;
        std::cerr << "FAIL: " << message << '\n';
    }
}
bool close(double a, double b) {
    // Mixed absolute/relative tolerance for reordered difference arithmetic.
    return std::abs(a - b) <= 1e-10 * (1 + std::max(std::abs(a), std::abs(b)));
}
template <class Exception, class F>
void throws(F action, const char* message) {
    try {
        action();
        check(false, message);
    } catch (const Exception&) {
    } catch (...) {
        check(false, message);
    }
}
double identity(annealer::Value value) { return static_cast<double>(value); }
double nonlinear(annealer::Value value) { return std::sin(static_cast<double>(value)); }

annealer::Model model() {
    using namespace annealer;
    ModelBuilder builder;
    builder.set_constant(0.25);
    builder.add_term(2.0, {});
    builder.add_term(-0.5, {});
    builder.add_variable(Domain{-3, 0, 7});
    builder.add_variable(Domain::cardinality(2));
    builder.add_variable(Domain{-1, 2, 8, 11});
    builder.add_variable(Domain{4});
    builder.add_term(-0.75, {native_factor<identity>(0), native_factor<nonlinear>(0),
                            native_factor<identity>(1), native_factor<nonlinear>(2)});
    builder.add_term(1.2, {{2, identity}, {3, identity}});
    builder.add_term(0, {{0, identity}});
    builder.add_term(0.5, {});
    return builder.build();
}

void differential() {
    using namespace annealer;
    const auto source = model();
    for (const std::size_t budget : {std::size_t{0}, std::size_t{40}, std::size_t{4096}}) {
        const auto prepared = prepare(source, {budget});
        check(prepared.lookup_bytes() <= budget, "table budget respected");
        check(prepared.term_count() == 4, "only leading constants folded");
        check(prepared.incident(0).size() == 2, "repeated variable has one incidence per term");
        check(prepared.incident(1).size() == 1, "incident list correct");
        for (const auto a : source.domains()[0].values()) {
            for (const auto b : source.domains()[1].values()) {
                for (const auto c : source.domains()[2].values()) {
                    Assignment assignment{a, b, c, 4};
                    check(close(prepared_energy(prepared, assignment), full_energy(source, assignment)),
                          "exhaustive prepared energies match reference");
                    for (Variable variable = 0; variable < 4; ++variable) {
                        const auto& values = source.domains()[variable].values();
                        for (std::size_t i = 0; i < values.size(); ++i) {
                            PreparedState state(prepared, assignment);
                            auto proposed = assignment;
                            proposed[variable] = values[i];
                            const double delta = state.propose(variable, i);
                            check(close(delta, full_energy(source, proposed) - full_energy(source, assignment)),
                                  "exhaustive incremental delta matches reference");
                            check(state.assignment() == assignment &&
                                  close(state.energy(), full_energy(source, assignment)),
                                  "unaccepted proposal leaves state unchanged");
                            state.accept();
                            check(state.assignment() == proposed && close(state.energy(), full_energy(source, proposed)),
                                  "accepted proposal commits assignment and caches");
                        }
                    }
                }
            }
        }
        PreparedState state(prepared, {-3, 0, -1, 4});
        std::mt19937_64 random(19);
        for (std::size_t move = 0; move < 10000; ++move) {
            const auto variable = static_cast<Variable>(random() % 3);
            const auto replacement = static_cast<std::size_t>(random() % source.domains()[variable].values().size());
            auto proposed = state.assignment();
            proposed[variable] = source.domains()[variable].values()[replacement];
            check(close(state.propose(variable, replacement),
                        full_energy(source, proposed) - full_energy(source, state.assignment())),
                  "long mixed accept/reject sequence delta");
            if (random() % 2) {
                state.accept();
            }
            if (move % 256 == 0) {
                state.rebase();
            }
            check(close(state.energy(), full_energy(source, state.assignment())), "bounded energy drift");
        }
        throws<std::invalid_argument>([&] { prepared_energy(prepared, {}); }, "assignment length checked");
        throws<std::invalid_argument>([&] { prepared_energy(prepared, {1, 0, -1, 4}); }, "actual values validated");
        throws<std::out_of_range>([&] { state.propose(0, 3); }, "bad position checked");
        throws<std::logic_error>([&] { state.accept(); }, "failed proposal cannot commit stale cache");
    }
}

void solver_comparison() {
    using namespace annealer;
    ModelBuilder builder;
    builder.add_variable(Domain{-2, 0, 3});
    builder.add_variable(Domain::cardinality(2));
    builder.add_term(-2, {{0, identity}, {0, identity}, {1, identity}});
    const auto source = builder.build();
    const auto prepared = prepare(source);
    SolverConfig config;
    config.seed = 812;
    config.temperature_levels = 20;
    config.sweeps_per_level = 80;
    const auto reference = solve_reference(source, config);
    const auto optimized = solve(prepared, config);
    const auto direct = solve(prepare(source, {0}), config);
    check(reference.best_assignment == optimized.best_assignment &&
          reference.final_assignment == optimized.final_assignment &&
          reference.accepted_moves == optimized.accepted_moves &&
          reference.attempted_moves == optimized.attempted_moves, "exact integer trajectory matches reference");
    check(direct.final_assignment == optimized.final_assignment && direct.accepted_moves == optimized.accepted_moves,
          "fallback and tables yield same trajectory");
    check(optimized.best_energy == full_energy(source, optimized.best_assignment) &&
          optimized.final_energy == full_energy(source, optimized.final_assignment), "returned energies recomputed");
    const auto fractional = model();
    config.sweeps_per_level = 100;
    const auto fractional_result = solve(prepare(fractional), config);
    check(fractional_result.best_energy == full_energy(fractional, fractional_result.best_assignment) &&
          fractional_result.final_energy == full_energy(fractional, fractional_result.final_assignment),
          "fractional solver returns reference energies");
}

void numerical_semantics() {
    using namespace annealer;
    ModelBuilder bad;
    bad.add_variable(Domain{0, 1});
    bad.add_term(0, {{0, [](Value x) { return x ? INFINITY : 0.0; }}});
    throws<std::domain_error>([&] { prepare(bad.build()); }, "zero term still validates nonfinite callback");
    const auto fallback = prepare(bad.build(), {0});
    check(prepared_energy(fallback, {0}) == 0, "fallback does not eagerly evaluate unused domain values");
    throws<std::domain_error>([&] { prepared_energy(fallback, {1}); }, "fallback checks callback result");
    ModelBuilder overflow;
    overflow.add_variable(Domain{1});
    const auto huge = [](Value) { return std::numeric_limits<double>::max(); };
    overflow.add_term(0, {{0, huge}, {0, huge}});
    const auto prepared = prepare(overflow.build());
    throws<std::overflow_error>([&] { prepared_energy(prepared, {1}); }, "zero term product overflow retained");
    ModelBuilder constants;
    constants.set_constant(std::numeric_limits<double>::max());
    constants.add_term(std::numeric_limits<double>::max(), {});
    throws<std::overflow_error>([&] { prepare(constants.build()); }, "constant folding checks arithmetic");
}
} // namespace

int main() {
    try {
        differential();
        solver_comparison();
        numerical_semantics();
    } catch (const std::exception& error) {
        std::cerr << "Unexpected exception: " << error.what() << '\n';
        return 1;
    }
    if (failures) {
        return 1;
    }
    std::cout << "All prepared model tests passed\n";
}
