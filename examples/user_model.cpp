#include "user_model.hpp"
#include "c_factors.h"
#include <annealer/prepared.hpp>

annealer::Model make_user_model() {
    annealer::ModelBuilder builder;
    const auto x = builder.add_variable(annealer::Domain(std::vector<annealer::Value>(
        user_model::x_domain.begin(), user_model::x_domain.end())));
    const auto y = builder.add_variable(annealer::Domain::cardinality(3));
    const auto z = builder.add_variable(annealer::Domain{3, 7});
    builder.set_constant(1.5);
    builder.add_term(user_model::interaction_coefficient,
        {annealer::native_factor<user_model::square>(x), {y, example_shift},
         annealer::native_factor<user_model::identity>(z)});
    builder.add_term(2.0, {annealer::native_factor<user_model::identity>(x)});
    builder.add_term(2.0, {});
    return builder.build();
}
