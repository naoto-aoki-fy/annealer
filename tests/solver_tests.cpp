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

void proposals() {
    // Exhaustively enumerate each equally likely rank for heterogeneous domains.
    // Every directed edge must occur exactly once, including its reverse.
    for (const std::vector<annealer::Value>& values : {
            std::vector<annealer::Value>{-7, 9}, {-11, 0, 4, 90}, {0, 1, 2}}) {
        const auto count = values.size();
        std::vector<std::vector<int>> edges(count, std::vector<int>(count));
        for (std::size_t current = 0; current < count; ++current) {
            for (std::size_t rank = 0; rank < count - 1; ++rank) {
                const auto proposed = annealer::detail::replacement_index(current, rank);
                check(proposed < count && proposed != current, "valid replacement without self moves");
                if (proposed < count) {
                    ++edges[current][proposed];
                }
            }
        }
        for (std::size_t a = 0; a < count; ++a) {
            for (std::size_t b = 0; b < count; ++b) {
                check(edges[a][b] == (a == b ? 0 : 1), "uniform alternatives");
                check(edges[a][b] == edges[b][a], "forward and reverse proposal probabilities match");
            }
        }
    }
    // Exercise the actual random scan with unequal domain sizes and a singleton.
    annealer::ModelBuilder builder;
    builder.add_variable(annealer::Domain{-90});
    builder.add_variable(annealer::Domain{-7, 9});
    builder.add_variable(annealer::Domain{-11, 0, 4, 90});
    const auto model = builder.build();
    std::mt19937_64 random(4321);
    std::vector<std::size_t> indices{0, 1, 2};
    std::vector<std::vector<int>> counts{{0}, {0, 0}, {0, 0, 0, 0}};
    for (int attempt = 0; attempt < 24000; ++attempt) {
        const auto proposal = annealer::detail::draw_proposal(model, {1, 2}, indices, random);
        check(proposal.variable == 1 || proposal.variable == 2, "draw selects a mutable variable");
        const auto count = model.domains()[proposal.variable].values().size();
        check(proposal.index < count && proposal.index != indices[proposal.variable],
              "draw gives a valid different value");
        ++counts[proposal.variable][proposal.index];
    }
    check(counts[1][0] > 11000 && counts[1][0] < 13000, "uniform variable selection sanity check");
    for (const auto index : {0, 1, 3}) {
        check(counts[2][index] > 3500 && counts[2][index] < 4500,
              "uniform replacement selection sanity check");
    }
}

void acceptance() {
    using annealer::metropolis_accept;
    check(metropolis_accept(-2, 1, 0.999), "downhill always accepted");
    check(metropolis_accept(0, 1, 0.999), "equal energy always accepted");
    check(metropolis_accept(1, 1, 0.3), "uphill below threshold accepted");
    check(!metropolis_accept(1, 1, 0.4), "uphill above threshold rejected");
    check(!metropolis_accept(1, 1, std::exp(-1.0)), "threshold equality rejects");
    check(metropolis_accept(1, 1, 0), "zero uniform accepted");
    check(!metropolis_accept(1000, 1, 0.1), "large exponent rejects without exp underflow");
    check(!metropolis_accept(std::numeric_limits<double>::max(),
                            std::numeric_limits<double>::min(), 0.5), "ratio overflow rejects");
    for (double u : {-0.1, 1.0, std::numeric_limits<double>::quiet_NaN()}) {
        throws<std::invalid_argument>([&] { metropolis_accept(1, 1, u); }, "bad uniform rejected");
    }
    throws<std::invalid_argument>([] { metropolis_accept(INFINITY, 1, 0.5); }, "nonfinite delta rejected");
}

void schedules() {
    using annealer::geometric_temperature;
    check(geometric_temperature(100, 1, 3, 0) == 100, "initial endpoint exact");
    check(std::abs(geometric_temperature(100, 1, 3, 1) - 10) < 1e-12, "geometric midpoint");
    check(geometric_temperature(100, 1, 3, 2) == 1, "final endpoint exact");
    check(geometric_temperature(100, 1, 1, 0) == 100, "one level uses initial temperature");
    check(geometric_temperature(4, 4, 3, 1) == 4, "constant schedule");
    const double tiny = std::numeric_limits<double>::denorm_min();
    const double huge = std::numeric_limits<double>::max();
    double previous = huge;
    for (std::size_t i = 0; i < 31; ++i) {
        const double t = geometric_temperature(huge, tiny, 31, i);
        check(std::isfinite(t) && t > 0 && t <= previous, "extreme schedule remains finite and monotonic");
        previous = t;
    }
    throws<std::out_of_range>([] { geometric_temperature(2, 1, 0, 0); }, "zero levels have no endpoint");
    throws<std::out_of_range>([] { geometric_temperature(2, 1, 2, 2); }, "invalid schedule index");
    throws<std::invalid_argument>([] { geometric_temperature(1, 2, 2, 0); }, "heating rejected");
}

void runs() {
    using namespace annealer;
    ModelBuilder builder;
    builder.add_variable(Domain{-8});
    builder.add_variable(Domain{-3, 2, 9});
    builder.add_variable(Domain::cardinality(2));
    builder.add_term(1, {{1, identity}, {1, identity}, {2, identity}});
    const auto model = builder.build();
    SolverConfig config;
    config.seed = 734;
    config.temperature_levels = 7;
    config.sweeps_per_level = 11;
    const auto a = solve(model, config);
    const auto b = solve(model, config);
    check(a.best_assignment == b.best_assignment && a.final_assignment == b.final_assignment &&
          a.best_energy == b.best_energy && a.final_energy == b.final_energy &&
          a.attempted_moves == b.attempted_moves && a.accepted_moves == b.accepted_moves &&
          a.seed == config.seed && b.seed == config.seed, "seed reproducibility with random initialization");
    check(a.attempted_moves == 2 * 7 * 11 && a.accepted_moves <= a.attempted_moves,
          "sweep budget counts mutable variables only");
    check(a.best_energy == full_energy(model, a.best_assignment) &&
          a.final_energy == full_energy(model, a.final_assignment), "reported assignments match energies");
    check(a.best_energy <= a.final_energy, "best energy no worse than final");
    check(a.best_assignment[0] == -8 && a.final_assignment[0] == -8, "singleton never changed");

    config.initial_assignment = Assignment{-8, 9, 1};
    const auto supplied = solve(model, config);
    check(supplied.best_energy <= 81, "best includes supplied initial state");
    const auto replay = solve(model, config);
    check(supplied.final_assignment == replay.final_assignment &&
          supplied.best_assignment == replay.best_assignment &&
          supplied.accepted_moves == replay.accepted_moves, "seed reproducibility with supplied assignment");
    config.temperature_levels = 0;
    auto zero = solve(model, config);
    check(zero.attempted_moves == 0 && zero.accepted_moves == 0 &&
          zero.best_assignment == *config.initial_assignment &&
          zero.final_assignment == zero.best_assignment && zero.best_energy == 81 && zero.final_energy == 81,
          "zero levels return initial state");
    config.temperature_levels = 5;
    config.sweeps_per_level = 0;
    check(solve(model, config).attempted_moves == 0, "zero sweeps make no moves");
    config.initial_assignment.reset();
    zero = solve(model, config);
    check(zero.best_energy == full_energy(model, zero.best_assignment), "zero budget random initialization valid");
    config.initial_assignment = Assignment{0, 0, 0};
    throws<std::invalid_argument>([&] { solve(model, config); }, "zero budget still validates assignment");
    config.initial_assignment = Assignment{};
    throws<std::invalid_argument>([&] { solve(model, config); }, "bad assignment size");
    config.initial_assignment.reset();
    for (const double bad : {0.0, -1.0, std::numeric_limits<double>::infinity(),
                            -std::numeric_limits<double>::infinity(),
                            std::numeric_limits<double>::quiet_NaN()}) {
        SolverConfig invalid;
        invalid.initial_temperature = bad;
        throws<std::invalid_argument>([&] { solve(model, invalid); }, "bad initial temperature");
        invalid.initial_temperature = 10;
        invalid.final_temperature = bad;
        throws<std::invalid_argument>([&] { solve(model, invalid); }, "bad final temperature");
    }
    config.temperature_levels = std::numeric_limits<std::size_t>::max();
    config.sweeps_per_level = 2;
    throws<std::overflow_error>([&] { solve(model, config); }, "budget overflow rejected before moves");

    ModelBuilder fixed;
    fixed.add_variable(Domain{-5});
    fixed.set_constant(7);
    const auto fixed_result = solve(fixed.build());
    check(fixed_result.final_assignment == Assignment{-5} && fixed_result.best_energy == 7 &&
          fixed_result.final_energy == 7 && fixed_result.attempted_moves == 0, "all singleton model");
    const auto empty = solve(ModelBuilder{}.build());
    check(empty.final_assignment.empty() && empty.best_assignment.empty() &&
          empty.final_energy == 0 && empty.accepted_moves == 0 && empty.attempted_moves == 0, "zero-variable model");
}

void state_tracking() {
    using namespace annealer;
    ModelBuilder builder;
    builder.add_variable(Domain{0, 1});
    builder.add_term(1, {{0, identity}});
    const auto model = builder.build();
    SolverConfig config;
    config.initial_assignment = Assignment{1};
    config.initial_temperature = config.final_temperature = std::numeric_limits<double>::max();
    config.temperature_levels = 1;
    config.sweeps_per_level = 2;
    auto result = solve(model, config);
    check(result.best_assignment == Assignment{0} && result.best_energy == 0 &&
          result.final_assignment == Assignment{1} && result.final_energy == 1 &&
          result.accepted_moves == 2, "best survives later accepted uphill move");
    config.initial_assignment = Assignment{0};
    config.sweeps_per_level = 1;
    result = solve(model, config);
    check(result.best_energy == 0 && result.final_energy == 1, "initial best retained");
    config.initial_temperature = config.final_temperature = std::numeric_limits<double>::min();
    result = solve(model, config);
    check(result.accepted_moves == 0 && result.final_energy == 0, "uphill rejection preserves state");
    ModelBuilder flat;
    flat.add_variable(Domain{-2, 10});
    config.initial_assignment = Assignment{-2};
    result = solve(flat.build(), config);
    check(result.accepted_moves == 1 && result.final_assignment == Assignment{10} &&
          result.best_assignment == Assignment{-2}, "equal move accepted, first best wins tie");

    ModelBuilder extremes;
    extremes.add_variable(Domain{-1, 1});
    extremes.add_term(std::numeric_limits<double>::max(), {{0, identity}});
    config.initial_assignment = Assignment{-1};
    throws<std::overflow_error>([&] { solve(extremes.build(), config); }, "delta overflow is explicit");
}
} // namespace

int main() {
    try {
        proposals();
        acceptance();
        schedules();
        runs();
        state_tracking();
    } catch (const std::exception& error) {
        std::cerr << "Unexpected exception: " << error.what() << '\n';
        return 1;
    }
    if (failures) {
        return 1;
    }
    std::cout << "All solver tests passed\n";
}
