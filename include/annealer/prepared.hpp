#pragma once

#include <annealer/model.hpp>
#include <span>

namespace annealer {

struct PreparationConfig {
    // Limit on stored double payload; metadata and the owned source are separate.
    std::size_t lookup_bytes = 64 * 1024 * 1024;
};

// The invoker is specialized on Function, exposing its body to normal native
// optimization when its definition is visible. No constexpr function is required.
template <auto Function>
Factor native_factor(Variable variable) {
    static_assert(Function != nullptr);
    return Factor(variable, [](Value value) { return Function(value); });
}

namespace detail {
inline double finite_arithmetic(double value) {
    if (!std::isfinite(value)) {
        throw std::overflow_error("nonfinite prepared energy arithmetic");
    }
    return value;
}
inline double finite_factor(double value) {
    if (!std::isfinite(value)) {
        throw std::domain_error("factor returned a nonfinite value during prepared evaluation");
    }
    return value;
}
} // namespace detail

class PreparedModel {
public:
    const Model& source() const noexcept { return source_; }
    std::size_t lookup_bytes() const noexcept { return table_.size() * sizeof(double); }
    std::size_t term_count() const noexcept { return terms_.size(); }
    std::span<const std::size_t> incident(Variable variable) const {
        if (variable >= source_.domains().size()) {
            throw std::out_of_range("unknown variable");
        }
        return std::span<const std::size_t>(incidence_).subspan(offsets_[variable], offsets_[variable + 1] - offsets_[variable]);
    }
    double constant() const noexcept { return constant_; }

    double contribution(std::size_t term_index, const std::vector<std::size_t>& indices,
                        Variable changed = std::numeric_limits<Variable>::max(),
                        std::size_t replacement = 0) const {
        const auto& term = terms_.at(term_index);
        double product = 1.0;
        for (std::size_t i = term.begin; i < term.end; ++i) {
            const auto& factor = factors_[i];
            const auto index = factor.variable == changed ? replacement : indices.at(factor.variable);
            // at() also validates manually supplied internal positions.
            const auto value = source_.domains()[factor.variable].values().at(index);
            const double unary = factor.table == no_table
                ? detail::finite_factor(source_.terms()[factor.source_term].factors()[factor.source_factor].evaluate(value))
                : table_[factor.table + index];
            product = detail::finite_arithmetic(product * unary);
        }
        return detail::finite_arithmetic(term.coefficient * product);
    }

private:
    friend PreparedModel prepare(const Model&, const PreparationConfig&);
    struct FlatTerm { double coefficient; std::size_t begin; std::size_t end; };
    struct FlatFactor {
        Variable variable;
        std::size_t table;
        std::size_t source_term;
        std::size_t source_factor;
    };
    static constexpr auto no_table = std::numeric_limits<std::size_t>::max();
    explicit PreparedModel(const Model& model) : source_(model), constant_(model.constant()) {}
    Model source_;
    double constant_;
    std::vector<FlatTerm> terms_;
    std::vector<FlatFactor> factors_;
    std::vector<double> table_;
    std::vector<std::size_t> offsets_;
    std::vector<std::size_t> incidence_;
};

inline PreparedModel prepare(const Model& model, const PreparationConfig& config = {}) {
    PreparedModel result(model);
    std::vector<std::vector<std::size_t>> incidence(model.domains().size());
    std::size_t entries = 0;
    bool leading_constants = true;
    for (std::size_t t = 0; t < model.terms().size(); ++t) {
        const auto& term = model.terms()[t];
        if (leading_constants && term.factors().empty()) {
            result.constant_ = detail::finite_arithmetic(result.constant_ + term.coefficient());
            continue;
        }
        leading_constants = false;
        const auto prepared_term = result.terms_.size();
        const auto begin = result.factors_.size();
        std::vector<Variable> variables;
        for (std::size_t f = 0; f < term.factors().size(); ++f) {
            const auto variable = term.factors()[f].variable();
            const auto size = model.domains()[variable].values().size();
            auto offset = PreparedModel::no_table;
            if (size <= config.lookup_bytes / sizeof(double) - entries) {
                offset = entries;
                entries += size;
            }
            result.factors_.push_back({variable, offset, t, f});
            variables.push_back(variable);
        }
        result.terms_.push_back({term.coefficient(), begin, result.factors_.size()});
        std::sort(variables.begin(), variables.end());
        variables.erase(std::unique(variables.begin(), variables.end()), variables.end());
        for (const auto variable : variables) {
            incidence[variable].push_back(prepared_term);
        }
    }
    result.table_.reserve(entries);
    for (const auto& factor : result.factors_) {
        if (factor.table != PreparedModel::no_table) {
            const auto& callback = model.terms()[factor.source_term].factors()[factor.source_factor];
            for (const auto value : model.domains()[factor.variable].values()) {
                result.table_.push_back(detail::finite_factor(callback.evaluate(value)));
            }
        }
    }
    result.offsets_.push_back(0);
    for (const auto& incident : incidence) {
        result.incidence_.insert(result.incidence_.end(), incident.begin(), incident.end());
        result.offsets_.push_back(result.incidence_.size());
    }
    return result;
}

// Mutable state is separate from the immutable prepared model, which must outlive it.
class PreparedState {
public:
    PreparedState(const PreparedModel& model, Assignment assignment)
        : model_(model), assignment_(std::move(assignment)) {
        if (assignment_.size() != model_.source().domains().size()) {
            throw std::invalid_argument("assignment size does not match model");
        }
        for (Variable i = 0; i < assignment_.size(); ++i) {
            const auto& values = model_.source().domains()[i].values();
            if (!model_.source().domains()[i].contains(assignment_[i])) {
                throw std::invalid_argument("assignment outside domain");
            }
            indices_.push_back(static_cast<std::size_t>(std::find(values.begin(), values.end(), assignment_[i]) - values.begin()));
        }
        contributions_.resize(model_.term_count());
        std::size_t largest = 0;
        for (Variable i = 0; i < assignment_.size(); ++i) {
            largest = std::max(largest, model_.incident(i).size());
        }
        pending_.resize(largest);
        energy_ = model_.constant();
        for (std::size_t t = 0; t < contributions_.size(); ++t) {
            contributions_[t] = model_.contribution(t, indices_);
            energy_ = detail::finite_arithmetic(energy_ + contributions_[t]);
        }
    }
    const Assignment& assignment() const noexcept { return assignment_; }
    const std::vector<std::size_t>& indices() const noexcept { return indices_; }
    double energy() const noexcept { return energy_; }

    double propose(Variable variable, std::size_t replacement) {
        pending_valid_ = false;
        model_.source().domains().at(variable).values().at(replacement);
        double delta = 0;
        std::size_t position = 0;
        for (const auto term : model_.incident(variable)) {
            const double value = model_.contribution(term, indices_, variable, replacement);
            pending_[position++] = value;
            delta = detail::finite_arithmetic(delta + detail::finite_arithmetic(value - contributions_[term]));
        }
        pending_energy_ = detail::finite_arithmetic(energy_ + delta);
        variable_ = variable;
        replacement_ = replacement;
        pending_valid_ = true;
        return delta;
    }

    void accept() {
        if (!pending_valid_) {
            throw std::logic_error("no pending proposal");
        }
        std::size_t position = 0;
        for (const auto term : model_.incident(variable_)) {
            contributions_[term] = pending_[position++];
        }
        assignment_[variable_] = model_.source().domains()[variable_].values()[replacement_];
        indices_[variable_] = replacement_;
        energy_ = pending_energy_;
        pending_valid_ = false;
    }

    // Refresh all caches and anchor the accumulated energy to the reference.
    void rebase() {
        pending_valid_ = false;
        for (std::size_t term = 0; term < contributions_.size(); ++term) {
            contributions_[term] = model_.contribution(term, indices_);
        }
        energy_ = full_energy(model_.source(), assignment_);
    }

private:
    const PreparedModel& model_;
    Assignment assignment_;
    std::vector<std::size_t> indices_;
    std::vector<double> contributions_;
    std::vector<double> pending_;
    double energy_ = 0;
    double pending_energy_ = 0;
    Variable variable_ = 0;
    std::size_t replacement_ = 0;
    bool pending_valid_ = false;
};

inline double prepared_energy(const PreparedModel& model, const Assignment& assignment) {
    PreparedState state(model, assignment);
    return state.energy();
}

} // namespace annealer
