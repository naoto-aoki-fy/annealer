#pragma once
#include "gset_maxcut.hpp"
#include <cstdint>
#include <set>

namespace covering {
using Cost = std::uint64_t;
inline constexpr Cost exact_limit = Cost{1} << 53;
struct Instance {
    std::vector<Cost> costs;
    std::vector<std::vector<std::size_t>> rows;
    Cost cost_sum = 0, penalty = 0;
};
inline Instance read(std::istream& input) {
    auto number = [&]() {
        std::string token;
        if (!(input >> token)) throw std::invalid_argument("Missing set covering field");
        return gset::unsigned_number<Cost>(token);
    };
    const auto m = number(), n = number();
    if (!m || !n || m > std::vector<std::vector<std::size_t>>{}.max_size() ||
        n > std::vector<Cost>{}.max_size())
        throw std::invalid_argument("Invalid row or column count");
    Instance instance;
    // Grow only as records arrive, rather than allocating from untrusted counts.
    for (Cost j = 0; j < n; ++j) {
        const auto cost = number();
        if (cost >= exact_limit - instance.cost_sum)
            throw std::invalid_argument("Cost sum leaves no exact integer penalty");
        instance.cost_sum += cost;
        instance.costs.push_back(cost);
    }
    instance.penalty = instance.cost_sum + 1;
    if (m > (exact_limit - instance.cost_sum) / instance.penalty)
        throw std::invalid_argument("Worst-case energy exceeds exact integer range");
    for (Cost i = 0; i < m; ++i) {
        const auto count = number();
        if (!count || count > n) throw std::invalid_argument("Row must have 1..N covering columns");
        std::vector<std::size_t> row;
        std::set<Cost> seen;
        for (Cost k = 0; k < count; ++k) {
            const auto column = number();
            if (!column || column > n || !seen.insert(column).second)
                throw std::invalid_argument("Invalid or duplicate covering column");
            row.push_back(static_cast<std::size_t>(column - 1));
        }
        instance.rows.push_back(std::move(row));
    }
    input >> std::ws;
    if (input.bad()) throw std::runtime_error("Set covering input read failed");
    if (!input.eof()) throw std::invalid_argument("Unexpected trailing set covering data");
    return instance;
}
inline annealer::Model make_model(const Instance& instance) {
    annealer::ModelBuilder builder;
    for (auto cost : instance.costs) {
        const auto variable = builder.add_variable(annealer::Domain::cardinality(2));
        builder.add_term(static_cast<double>(cost), {{variable, [](annealer::Value x) { return static_cast<double>(x); }}});
    }
    for (const auto& row : instance.rows) {
        std::vector<annealer::Factor> factors;
        for (auto column : row)
            factors.emplace_back(column, [](annealer::Value x) { return 1.0 - static_cast<double>(x); });
        builder.add_term(static_cast<double>(instance.penalty), std::move(factors));
    }
    return builder.build();
}
struct Objective { Cost selected_cost = 0; std::size_t uncovered_rows = 0; };
// Independent of the factor model: inspect selected columns and row coverage.
inline Objective evaluate(const Instance& instance, const annealer::Assignment& x) {
    if (x.size() != instance.costs.size()) throw std::invalid_argument("Assignment size mismatch");
    Objective result;
    for (std::size_t j = 0; j < x.size(); ++j) {
        if (x[j] != 0 && x[j] != 1) throw std::invalid_argument("Assignment must be Boolean");
        if (x[j]) result.selected_cost += instance.costs[j];
    }
    for (const auto& row : instance.rows) {
        bool covered = false;
        for (auto j : row) covered |= x.at(j) == 1;
        if (!covered) ++result.uncovered_rows;
    }
    return result;
}
inline double energy(const Instance& instance, Objective objective) {
    return static_cast<double>(objective.selected_cost + instance.penalty * objective.uncovered_rows);
}
} // namespace covering
