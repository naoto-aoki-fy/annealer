#pragma once

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace annealer {

using Value = std::int64_t;
using Variable = std::size_t;
using Assignment = std::vector<Value>;

class Domain {
public:
    explicit Domain(std::vector<Value> values) : values_(std::move(values)) {
        if (values_.empty()) {
            throw std::invalid_argument("a domain must not be empty");
        }
        auto sorted = values_;
        std::sort(sorted.begin(), sorted.end());
        if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end()) {
            throw std::invalid_argument("domain values must be distinct");
        }
    }

    Domain(std::initializer_list<Value> values) : Domain(std::vector<Value>(values)) {}

    static Domain cardinality(std::size_t count) {
        if (count == 0) {
            throw std::invalid_argument("domain cardinality must be positive");
        }
        if (static_cast<std::uintmax_t>(count - 1) >
            static_cast<std::uintmax_t>(std::numeric_limits<Value>::max())) {
            throw std::invalid_argument("domain cardinality exceeds the value range");
        }
        std::vector<Value> values;
        values.reserve(count);
        for (std::size_t i = 0; i < count; ++i) {
            values.push_back(static_cast<Value>(i));
        }
        return Domain(std::move(values));
    }

    const std::vector<Value>& values() const noexcept { return values_; }
    bool contains(Value value) const noexcept {
        return std::find(values_.begin(), values_.end(), value) != values_.end();
    }

private:
    std::vector<Value> values_;
};

// A callback must be deterministic and side-effect-free. Its argument is an
// actual domain value, never an index. Fixed captured constants are permitted.
class Factor {
public:
    template <class F>
        requires std::copy_constructible<std::decay_t<F>> &&
                 std::is_invocable_r_v<double, const std::decay_t<F>&, Value>
    Factor(Variable variable, F function)
        : variable_(variable), function_(std::move(function)) {
        if constexpr (std::is_pointer_v<std::decay_t<F>>) {
            if (function == nullptr) {
                throw std::invalid_argument("factor function must not be null");
            }
        }
        if (!function_) {
            throw std::invalid_argument("factor function must not be empty");
        }
    }

    Variable variable() const noexcept { return variable_; }
    double evaluate(Value value) const { return function_(value); }

private:
    friend class ModelBuilder;
    Variable variable_;
    std::function<double(Value)> function_;
};

class Term {
public:
    double coefficient() const noexcept { return coefficient_; }
    const std::vector<Factor>& factors() const noexcept { return factors_; }

private:
    friend class ModelBuilder;
    Term(double coefficient, std::vector<Factor> factors)
        : coefficient_(coefficient), factors_(std::move(factors)) {}

    double coefficient_;
    std::vector<Factor> factors_;
};

class Model {
public:
    double constant() const noexcept { return constant_; }
    const std::vector<Domain>& domains() const noexcept { return domains_; }
    const std::vector<Term>& terms() const noexcept { return terms_; }

private:
    friend class ModelBuilder;
    Model(double constant, std::vector<Domain> domains, std::vector<Term> terms)
        : constant_(constant), domains_(std::move(domains)), terms_(std::move(terms)) {}

    double constant_;
    std::vector<Domain> domains_;
    std::vector<Term> terms_;
};

class ModelBuilder {
public:
    Variable add_variable(Domain domain) {
        if (domain.values().empty()) {
            throw std::invalid_argument("a domain must not be empty (including moved-from domains)");
        }
        domains_.push_back(std::move(domain));
        return domains_.size() - 1;
    }

    void set_constant(double constant) {
        require_finite(constant);
        constant_ = constant;
    }

    void add_term(double coefficient, std::vector<Factor> factors) {
        require_finite(coefficient);
        for (const auto& factor : factors) {
            if (!factor.function_) {
                throw std::invalid_argument("factor function must not be empty");
            }
            if (factor.variable() >= domains_.size()) {
                throw std::out_of_range("factor references an unknown variable");
            }
        }
        terms_.push_back(Term(coefficient, std::move(factors)));
    }

    // The result owns a snapshot; later builder changes cannot alter it.
    Model build() const { return Model(constant_, domains_, terms_); }

private:
    static void require_finite(double value) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("model constants and coefficients must be finite");
        }
    }

    double constant_ = 0.0;
    std::vector<Domain> domains_;
    std::vector<Term> terms_;
};

// Deliberately direct: no tables, cached products, or incremental updates.
inline double full_energy(const Model& model, const Assignment& assignment) {
    if (assignment.size() != model.domains().size()) {
        throw std::invalid_argument("assignment size does not match the model");
    }
    for (std::size_t i = 0; i < assignment.size(); ++i) {
        if (!model.domains()[i].contains(assignment[i])) {
            throw std::invalid_argument("assignment value is outside its domain");
        }
    }

    const auto checked = [](double value) {
        if (!std::isfinite(value)) {
            throw std::overflow_error("nonfinite intermediate energy arithmetic");
        }
        return value;
    };
    double energy = model.constant();
    for (const auto& term : model.terms()) {
        double product = 1.0;
        for (const auto& factor : term.factors()) {
            const double value = factor.evaluate(assignment[factor.variable()]);
            if (!std::isfinite(value)) {
                throw std::domain_error("factor returned a nonfinite value");
            }
            product = checked(product * value);
        }
        const double contribution = checked(term.coefficient() * product);
        energy = checked(energy + contribution);
    }
    return energy;
}

} // namespace annealer
