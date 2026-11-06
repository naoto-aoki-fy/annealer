#include "user_model.hpp"
#include "basic_execution.hpp"

int main(int argc, char** argv) {
    basic::Options options;
    options.solver.seed = 42;
    options.solver.initial_assignment = annealer::Assignment{-2, 2, 3};
    return basic::run(argc, argv, options, make_user_model, [](const auto& model) {
        std::cout << "Energy: " << annealer::full_energy(model, {-2, 2, 3}) << '\n';
    }, [](const auto& result) {
        std::cout << "Best energy: " << result.best_energy << " at";
        for (const auto value : result.best_assignment) {
            std::cout << ' ' << value;
        }
        std::cout << "\nFinal energy: " << result.final_energy
                  << "\nMoves: " << result.accepted_moves << '/' << result.attempted_moves
                  << " accepted; seed: " << result.seed << '\n';
    });
}
