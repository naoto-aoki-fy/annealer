#include "user_model.hpp"
#include <annealer/solver.hpp>

#include <exception>
#include <iostream>

int main() {
    try {
        const auto model = make_user_model();
        // 1.5 - 0.5 * 4 * 3 * 3 + 2 * (-2) + 2 = -18.5
        std::cout << "Energy: " << annealer::full_energy(model, {-2, 2, 3}) << '\n';
        annealer::SolverConfig config;
        config.seed = 42;
        config.initial_assignment = annealer::Assignment{-2, 2, 3};
        const auto result = annealer::solve(model, config);
        std::cout << "Best energy: " << result.best_energy << " at";
        for (const auto value : result.best_assignment) {
            std::cout << ' ' << value;
        }
        std::cout << "\nFinal energy: " << result.final_energy
                  << "\nMoves: " << result.accepted_moves << '/' << result.attempted_moves
                  << " accepted; seed: " << result.seed << '\n';
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
