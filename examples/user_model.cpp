#include "user_model.hpp"
#include "c_factors.h"

annealer::Model make_user_model() {
    annealer::ModelBuilder builder;
    const auto x = builder.add_variable(annealer::Domain{-2, 1, 4});
    const auto y = builder.add_variable(annealer::Domain::cardinality(3));
    const auto z = builder.add_variable(annealer::Domain{3, 7});
    builder.set_constant(1.5);
    const auto square = [](std::int64_t value) {
        const double real = static_cast<double>(value);
        return real * real;
    };
    const auto identity = [](std::int64_t value) { return static_cast<double>(value); };
    builder.add_term(-0.5, {{x, square}, {y, example_shift}, {z, identity}});
    builder.add_term(2.0, {{x, identity}});
    builder.add_term(2.0, {});
    return builder.build();
}
