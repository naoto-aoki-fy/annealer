#include "../examples/demo_runner.hpp"
#include "../examples/user_model.hpp"

namespace {
void require(bool ok) { if (!ok) { throw std::runtime_error("integration check failed"); } }
bool close(double a, double b) { return std::abs(a-b) <= 1e-10 * (1 + std::max(std::abs(a), std::abs(b))); }

template <class Formula>
void verify(const annealer::Model& model, double expected_minimum, Formula formula) {
    double minimum = std::numeric_limits<double>::infinity();
    for (const std::size_t budget : {std::size_t{0}, std::size_t{24}, std::size_t{65536}}) {
        const auto prepared = annealer::prepare(model, {budget});
        demos::enumerate(model, [&](const auto& assignment) {
            const double energy = annealer::full_energy(model, assignment);
            require(close(energy, formula(assignment)));
            minimum = std::min(minimum, energy);
            require(close(annealer::prepared_energy(prepared, assignment), energy));
            for (annealer::Variable v = 0; v < assignment.size(); ++v) {
                const auto& values = model.domains()[v].values();
                for (std::size_t index = 0; index < values.size(); ++index) {
                    if (values[index] == assignment[v]) { continue; }
                    annealer::PreparedState state(prepared, assignment);
                    auto proposed = assignment;
                    proposed[v] = values[index];
                    const double expected = annealer::full_energy(model, proposed);
                    const double delta = state.propose(v, index);
                    require(close(delta, expected - energy));
                    require(state.assignment() == assignment && close(state.energy(), energy));
                    // Reject, repeat from the unchanged cache, then accept.
                    require(close(state.propose(v, index), delta));
                    state.accept();
                    require(state.assignment() == proposed && close(state.energy(), expected));
                    // A reverse proposal verifies that every cached contribution was committed.
                    const auto old_index = static_cast<std::size_t>(std::find(values.begin(), values.end(), assignment[v]) - values.begin());
                    require(close(state.propose(v, old_index), -delta));
                }
            }
        });
        annealer::SolverConfig config;
        config.seed = 77;
        const auto result = annealer::solve(prepared, config);
        require(result.best_energy == annealer::full_energy(model, result.best_assignment));
        require(result.final_energy == annealer::full_energy(model, result.final_assignment));
        require(result.best_energy >= expected_minimum); // No stochastic success requirement.
    }
    require(minimum == expected_minimum);
}
} // namespace

int main() {
    try {
        verify(demos::binary_hubo(), -3, [](const auto& x) {
            return 1.0 - 2*x[0] - x[1] - 3*x[0]*x[1]*x[2] + 2*x[0]*x[2];
        });
        verify(make_user_model(), -156.5, [](const auto& x) {
            return 3.5 - 0.5*x[0]*x[0]*(x[1]+1)*x[2] + 2*x[0];
        });
        verify(demos::nonlinear(), -1710, [](const auto& x) {
            return 2.0 - double(x[0])*x[0]*x[0]*(x[1]+1)*x[2]*x[2] + double(x[0])*x[0];
        });
        annealer::ModelBuilder constants;
        constants.set_constant(-2);
        constants.add_term(0, {});
        constants.add_term(3, {});
        verify(constants.build(), 1, [](const auto&) { return 1; });
        std::cout << "All integration checks passed\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
