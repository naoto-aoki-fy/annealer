#pragma once
#include "demo_models.hpp"
#include "basic_execution.hpp"
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

template<class MakeModel>
inline int run(int argc, char** argv, MakeModel make_model) {
    double minimum = std::numeric_limits<double>::infinity();
    basic::Options options;
    options.solver.seed = 0;
    options.trials = 5;
    return basic::run(argc, argv, options, make_model, [&](const auto& model) {
        enumerate(model, [&](const auto& state) { minimum = std::min(minimum, annealer::full_energy(model, state)); });
        std::cout << "Exact reference minimum: " << minimum << '\n';
    }, [&](const auto& result) {
        std::cout << "Seed " << result.seed << ": best " << result.best_energy
                  << ", gap " << result.best_energy - minimum << ", final " << result.final_energy << '\n';
    });
}
} // namespace demos
