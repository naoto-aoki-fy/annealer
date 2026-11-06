#pragma once
#include "gset_maxcut.hpp"
#include <memory>
#include <numeric>
#include <random>

namespace qap {
using Cost = std::uint64_t;
inline constexpr Cost limit = Cost{1} << 53;
inline Cost multiply(Cost a, Cost b) {
    if (b && a > limit / b) throw std::invalid_argument("QAP arithmetic exceeds exact integer range");
    return a * b;
}
inline Cost add(Cost a, Cost b) {
    if (b > limit - a) throw std::invalid_argument("QAP arithmetic exceeds exact integer range");
    return a + b;
}
struct Instance { std::size_t size; std::vector<Cost> a, b; Cost penalty; };
inline Instance read(std::istream& input) {
    auto number = [&]() {
        std::string token;
        if (!(input >> token)) throw std::invalid_argument("Incomplete QAPLIB input");
        const auto value = gset::unsigned_number<Cost>(token);
        if (value > limit) throw std::invalid_argument("QAP integer exceeds exact range");
        return value;
    };
    const auto n = number();
    if (!n || n > static_cast<Cost>(std::numeric_limits<annealer::Value>::max()) ||
        n > std::vector<Cost>{}.max_size() / n)
        throw std::invalid_argument("Invalid QAP dimensions");
    Instance instance{static_cast<std::size_t>(n), {}, {}, 0};
    Cost sum = 0, maximum = 0;
    for (Cost i = 0; i < n * n; ++i) { auto v = number(); sum = add(sum, v); instance.a.push_back(v); }
    for (Cost i = 0; i < n * n; ++i) { auto v = number(); maximum = std::max(maximum, v); instance.b.push_back(v); }
    const auto upper = multiply(sum, maximum);
    instance.penalty = add(upper, 1);
    const auto pairs = multiply(n, n - 1) / 2;
    (void)add(upper, multiply(instance.penalty, pairs));
    input >> std::ws;
    if (input.bad()) throw std::runtime_error("QAPLIB read failed");
    if (!input.eof()) throw std::invalid_argument("Unexpected trailing QAPLIB data");
    return instance;
}
inline annealer::Model make_model(const Instance& instance) {
    annealer::ModelBuilder builder;
    const auto n = instance.size;
    for (std::size_t i = 0; i < n; ++i) builder.add_variable(annealer::Domain::cardinality(n));
    // One immutable matrix snapshot shared by the pure unary callbacks.
    const auto b = std::make_shared<const std::vector<Cost>>(instance.b);
    bool symmetric = true;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) symmetric &= (*b)[i*n+j] == (*b)[j*n+i];
    auto indicator = [](std::size_t a) { return [a](annealer::Value x) { return x == static_cast<annealer::Value>(a) ? 1.0 : 0.0; }; };
    for (std::size_t i = 0; i < n; ++i) {
        if (instance.a[i*n+i])
            builder.add_term(static_cast<double>(instance.a[i*n+i]), {{i, [b, n](annealer::Value x) { return static_cast<double>((*b)[x*n+x]); }}});
        for (std::size_t j = 0; j < n; ++j) {
            if (i == j || (symmetric && j < i)) continue;
            auto coefficient = instance.a[i*n+j];
            if (symmetric) coefficient = add(coefficient, instance.a[j*n+i]);
            if (!coefficient) continue;
            for (std::size_t a = 0; a < n; ++a) {
                bool nonzero = false;
                for (std::size_t k = 0; k < n; ++k) nonzero |= (*b)[a*n+k] != 0;
                if (nonzero) builder.add_term(static_cast<double>(coefficient), {{i, indicator(a)},
                    {j, [b, a, n](annealer::Value x) { return static_cast<double>((*b)[a*n+x]); }}});
            }
        }
        for (std::size_t j = i + 1; j < n; ++j)
            for (std::size_t a = 0; a < n; ++a)
                builder.add_term(static_cast<double>(instance.penalty), {{i, indicator(a)}, {j, indicator(a)}});
    }
    return builder.build();
}
struct Objective { Cost cost = 0, collisions = 0, penalty_cost = 0; bool feasible = false; };
inline Objective evaluate(const Instance& instance, const annealer::Assignment& x) {
    const auto n = instance.size;
    if (x.size() != n) throw std::invalid_argument("QAP assignment size mismatch");
    for (auto value : x) if (value < 0 || static_cast<std::uint64_t>(value) >= n)
        throw std::invalid_argument("QAP assignment outside domain");
    Objective result;
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) result.cost += instance.a[i*n+j] * instance.b[x[i]*n+x[j]];
        for (std::size_t j = i + 1; j < n; ++j) result.collisions += x[i] == x[j];
    }
    result.penalty_cost = instance.penalty * result.collisions;
    result.feasible = result.collisions == 0;
    return result;
}
inline double energy(const Instance&, Objective objective) { return static_cast<double>(objective.cost + objective.penalty_cost); }
inline annealer::Assignment permutation(std::size_t n, std::uint64_t seed) {
    annealer::Assignment x(n); std::iota(x.begin(), x.end(), annealer::Value{0});
    std::mt19937_64 random(seed); std::shuffle(x.begin(), x.end(), random); return x;
}
} // namespace qap
