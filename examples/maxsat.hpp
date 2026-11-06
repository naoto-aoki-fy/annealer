#pragma once
#include "gset_maxcut.hpp"
#include <sstream>
#include <optional>

namespace maxsat {
using Weight = std::uint64_t;
inline constexpr Weight exact_limit = Weight{1} << 53;
struct Literal { std::size_t variable; bool negative; };
struct Clause { Weight weight; bool hard; std::vector<Literal> literals; };
struct Instance {
    std::size_t variables = 0;
    std::vector<Clause> clauses;
    Weight soft_sum = 0, penalty = 0;
};

// Integer arithmetic is bounded before conversion to the model's double coefficients.
inline Instance read(std::istream& input) {
    Instance instance;
    bool header = false, records = false;
    std::optional<std::size_t> expected;
    std::optional<Weight> top;
    std::optional<Clause> current;
    std::string line;
    while (std::getline(input, line)) {
        std::istringstream row(line);
        std::string token;
        if (!(row >> token) || token == "c") continue;
        if (token == "p") {
            std::string format, n, m, t, extra;
            if (header || records || !(row >> format >> n >> m) || format != "wcnf")
                throw std::invalid_argument("Expected one initial p wcnf N M [TOP] header");
            header = true;
            instance.variables = gset::unsigned_number<std::size_t>(n);
            expected = gset::unsigned_number<std::size_t>(m);
            if (row >> t) {
                top = gset::unsigned_number<Weight>(t);
                if (!*top || (row >> extra)) throw std::invalid_argument("Invalid top weight/header");
            }
            continue;
        }
        do {
            records = true;
            if (!current) {
                bool hard = token == "h";
                if (hard && top) throw std::invalid_argument("Explicit h requires modern WCNF");
                Weight weight = hard ? 0 : gset::unsigned_number<Weight>(token);
                if (!hard && !weight) throw std::invalid_argument("Weight must be positive");
                if (top && weight >= *top) hard = true;
                current = Clause{weight, hard, {}};
            } else if (token == "0") {
                if (!current->hard) {
                    if (current->weight > exact_limit - instance.soft_sum)
                        throw std::invalid_argument("Soft weight sum exceeds exact integer range");
                    instance.soft_sum += current->weight;
                }
                instance.clauses.push_back(std::move(*current));
                current.reset();
                if (expected && instance.clauses.size() > *expected)
                    throw std::invalid_argument("Too many clauses");
            } else {
                const bool negative = token.front() == '-';
                const auto index = gset::unsigned_number<std::size_t>(negative ? token.substr(1) : token);
                if (!index || index == std::numeric_limits<std::size_t>::max() || (header && index > instance.variables))
                    throw std::invalid_argument("Literal index out of range");
                if (!header) instance.variables = std::max(instance.variables, index);
                current->literals.push_back({index - 1, negative});
            }
        } while (row >> token);
    }
    if (input.bad()) throw std::runtime_error("WCNF input read failed");
    if (current || (!header && !records) || (expected && instance.clauses.size() != *expected))
        throw std::invalid_argument("Unterminated clause, missing input or clause count mismatch");
    std::size_t hard_count = 0;
    for (const auto& clause : instance.clauses) hard_count += clause.hard;
    if (hard_count) {
        if (instance.soft_sum == exact_limit) throw std::invalid_argument("Penalty exceeds exact integer range");
        instance.penalty = instance.soft_sum + 1;
        if (hard_count > (exact_limit - instance.soft_sum) / instance.penalty)
            throw std::invalid_argument("Worst-case energy exceeds exact integer range");
    }
    return instance;
}

inline annealer::Model make_model(const Instance& instance) {
    annealer::ModelBuilder builder;
    for (std::size_t i = 0; i < instance.variables; ++i)
        builder.add_variable(annealer::Domain::cardinality(2));
    for (const auto& clause : instance.clauses) {
        std::vector<annealer::Factor> factors;
        for (const auto literal : clause.literals)
            factors.emplace_back(literal.variable, [negative = literal.negative](annealer::Value x) {
                return negative ? static_cast<double>(x) : 1.0 - static_cast<double>(x);
            });
        builder.add_term(static_cast<double>(clause.hard ? instance.penalty : clause.weight), std::move(factors));
    }
    return builder.build();
}
struct Objective { Weight soft_cost = 0; std::size_t hard_violations = 0; };
// Direct Boolean clause evaluation, independent of model factors and preparation.
inline Objective evaluate(const Instance& instance, const annealer::Assignment& x) {
    if (x.size() != instance.variables) throw std::invalid_argument("Assignment size mismatch");
    for (auto value : x) if (value != 0 && value != 1) throw std::invalid_argument("Assignment must be Boolean");
    Objective result;
    for (const auto& clause : instance.clauses) {
        bool satisfied = false;
        for (const auto literal : clause.literals)
            satisfied |= literal.negative ? x.at(literal.variable) == 0 : x.at(literal.variable) == 1;
        if (!satisfied) {
            if (clause.hard) ++result.hard_violations;
            else result.soft_cost += clause.weight;
        }
    }
    return result;
}
inline double energy(const Instance& instance, Objective objective) {
    return static_cast<double>(objective.soft_cost + instance.penalty * objective.hard_violations);
}
} // namespace maxsat
