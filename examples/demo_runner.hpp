#pragma once
#include "demo_models.hpp"
#include <iostream>

namespace demos {
// For tiny examples only; enumerate actual domain values in registration order.
template <class Visit>
void enumerate(const annealer::Model& model, Visit visit) {
    std::size_t count = 1;
    for (const auto& domain : model.domains()) {
        if (domain.values().size() > 4096 / count) {
            throw std::invalid_argument("demo enumeration is limited to 4096 states");
        }
        count *= domain.values().size();
    }
    for (std::size_t code = 0; code < count; ++code) {
        auto remaining = code;
        annealer::Assignment state;
        for (const auto& domain : model.domains()) {
            state.push_back(domain.values()[remaining % domain.values().size()]);
            remaining /= domain.values().size();
        }
        visit(state);
    }
}

inline int run(const annealer::Model& model) {
    double minimum = std::numeric_limits<double>::infinity();
    enumerate(model, [&](const auto& state) { minimum = std::min(minimum, annealer::full_energy(model, state)); });
    std::cout << "Exact reference minimum: " << minimum << '\n';
    const auto prepared = annealer::prepare(model);
    // Optimization-quality demonstration, not a correctness assertion.
    for (const std::uint64_t seed : {0, 1, 2, 3, 4}) {
        annealer::SolverConfig config;
        config.seed = seed;
        const auto result = annealer::solve(prepared, config);
        std::cout << "Seed " << seed << ": best " << result.best_energy
                  << ", gap " << result.best_energy - minimum << ", final " << result.final_energy << '\n';
    }
    return 0;
}
} // namespace demos
