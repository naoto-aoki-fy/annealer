#include <annealer/solver.hpp>
#include <iostream>

double square(std::int64_t value) {
    const double x = static_cast<double>(value);
    return x * x;
}

int main() {
    annealer::ModelBuilder builder;
    const auto x = builder.add_variable(annealer::Domain{-2, 1, 4});
    const auto y = builder.add_variable(annealer::Domain::cardinality(3));
    builder.set_constant(1.0);
    builder.add_term(-0.5, {{x, square}, {y, [](std::int64_t v) { return double(v) + 1; }}});
    const auto model = builder.build();
    const auto prepared = annealer::prepare(model, {.lookup_bytes = 1024});
    annealer::SolverConfig config;
    config.seed = 42;
    config.initial_assignment = annealer::Assignment{-2, 0};
    config.initial_temperature = 10;
    config.final_temperature = 0.1;
    config.temperature_levels = 100;
    config.sweeps_per_level = 10;
    const auto result = annealer::solve(prepared, config);
    std::cout << "Best energy: " << result.best_energy << '\n';
}
