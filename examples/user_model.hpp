#pragma once

#include <annealer/model.hpp>
#include <array>

namespace user_model {
inline constexpr std::array<annealer::Value, 3> x_domain{-2, 1, 4};
inline constexpr double interaction_coefficient = -0.5;
inline double square(annealer::Value value) {
    const double real = static_cast<double>(value);
    return real * real;
}
inline double identity(annealer::Value value) { return static_cast<double>(value); }
} // namespace user_model

annealer::Model make_user_model();
