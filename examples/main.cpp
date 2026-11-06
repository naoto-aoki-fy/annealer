#include "user_model.hpp"

#include <exception>
#include <iostream>

int main() {
    try {
        const auto model = make_user_model();
        // 1.5 - 0.5 * 4 * 3 * 3 + 2 * (-2) + 2 = -18.5
        std::cout << "Energy: " << annealer::full_energy(model, {-2, 2, 3}) << '\n';
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
