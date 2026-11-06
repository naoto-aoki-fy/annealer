#pragma once

#include <annealer/model.hpp>
#include <annealer/prepared.hpp>

#include <optional>
#include <memory>
#include <array>
#include <chrono>
#include <random>

namespace annealer {

struct SolverConfig {
    std::uint64_t seed = 0;
    std::optional<Assignment> initial_assignment;
    double initial_temperature = 10.0;
    double final_temperature = 0.1;
    std::size_t temperature_levels = 100;
    std::size_t sweeps_per_level = 10;
    bool swap_proposals = false;
};

struct SolverResult {
    Assignment best_assignment;
    double best_energy;
    Assignment final_assignment;
    double final_energy;
    std::size_t attempted_moves;
    std::size_t accepted_moves;
    std::uint64_t seed;
};

struct BestUpdate {
    double seconds;
    std::size_t attempted_moves, accepted_moves;
    double best_energy;
};
struct ProfileInterval {
    double seconds = 0;
    std::size_t calls = 0;
};
inline constexpr std::array<const char*, 6> profile_names{
    "selection", "propose", "metropolis", "accept", "rebase", "best_candidate"};
struct InstrumentationConfig { bool profile = false; };
struct InstrumentedResult {
    SolverResult result;
    std::vector<BestUpdate> history;
    double seconds = 0;
    std::array<ProfileInterval, 6> profile{};
};

namespace detail {
inline void validate_temperatures(double initial, double final) {
    if (!std::isfinite(initial) || !std::isfinite(final) || initial <= 0.0 || final <= 0.0) {
        throw std::invalid_argument("temperatures must be positive and finite");
    }
    if (final > initial) {
        throw std::invalid_argument("final temperature must not exceed initial temperature");
    }
}

// Map a uniform rank in [0, count-2] bijectively to all indices except current.
inline std::size_t replacement_index(std::size_t current, std::size_t rank) {
    return rank >= current ? rank + 1 : rank;
}

struct Proposal {
    Variable variable;
    std::size_t index;
};

inline Proposal draw_proposal(const Model& model, const std::vector<Variable>& mutable_variables,
                              const std::vector<std::size_t>& indices, std::mt19937_64& random) {
    std::uniform_int_distribution<std::size_t> choose_variable(0, mutable_variables.size() - 1);
    const auto variable = mutable_variables[choose_variable(random)];
    const auto count = model.domains()[variable].values().size();
    std::uniform_int_distribution<std::size_t> choose_replacement(0, count - 2);
    return {variable, replacement_index(indices[variable], choose_replacement(random))};
}

inline void validate_swap_domains(const Model& model, bool swaps) {
    if (!swaps || model.domains().empty()) return;
    for (const auto& domain : model.domains())
        if (domain.values() != model.domains().front().values())
            throw std::invalid_argument("swap proposals require identical variable domains");
}

inline Proposal draw_swap(std::size_t count, std::mt19937_64& random) {
    std::uniform_int_distribution<std::size_t> first(0, count - 1), rank(0, count - 2);
    const auto variable = first(random);
    return {variable, replacement_index(variable, rank(random))};
}

inline double temperature_unchecked(double initial, double final,
                                    std::size_t levels, std::size_t level) {
    if (level == 0) {
        return initial;
    }
    if (level == levels - 1) {
        return final;
    }
    const double fraction = static_cast<double>(level) / static_cast<double>(levels - 1);
    // Log interpolation avoids overflow/underflow in final / initial.
    const double temperature = std::exp(std::lerp(std::log(initial), std::log(final), fraction));
    return std::clamp(temperature, final, initial);
}

inline std::size_t checked_count(std::size_t a, std::size_t b) {
    if (b != 0 && a > std::numeric_limits<std::size_t>::max() / b) {
        throw std::overflow_error("proposal budget exceeds the move counter range");
    }
    return a * b;
}
} // namespace detail

inline double geometric_temperature(double initial, double final,
                                    std::size_t levels, std::size_t level) {
    detail::validate_temperatures(initial, final);
    if (level >= levels) {
        throw std::out_of_range("temperature level is outside the schedule");
    }
    return detail::temperature_unchecked(initial, final, levels, level);
}

// uniform must be in [0, 1). Equality with the uphill threshold rejects.
inline bool metropolis_accept(double delta_energy, double temperature, double uniform) {
    if (!std::isfinite(delta_energy) || !std::isfinite(temperature) || temperature <= 0.0 ||
        !std::isfinite(uniform) || uniform < 0.0 || uniform >= 1.0) {
        throw std::invalid_argument("Metropolis inputs require finite delta, positive temperature, and uniform in [0,1)");
    }
    if (delta_energy <= 0.0 || uniform == 0.0) {
        return true;
    }
    // Equivalent to uniform < exp(-delta/T), without exponential underflow.
    return std::log(uniform) < -(delta_energy / temperature);
}

inline SolverResult solve_reference(const Model& model, const SolverConfig& config = {}) {
    detail::validate_temperatures(config.initial_temperature, config.final_temperature);
    detail::validate_swap_domains(model, config.swap_proposals);
    std::mt19937_64 random(config.seed);
    Assignment current;
    if (config.initial_assignment) {
        current = *config.initial_assignment;
    } else {
        current.reserve(model.domains().size());
        for (const auto& domain : model.domains()) {
            const auto& values = domain.values();
            std::uniform_int_distribution<std::size_t> choose(0, values.size() - 1);
            current.push_back(values[choose(random)]);
        }
    }
    double current_energy = full_energy(model, current);
    SolverResult result{current, current_energy, current, current_energy, 0, 0, config.seed};

    std::vector<Variable> mutable_variables;
    std::vector<std::size_t> indices;
    indices.reserve(current.size());
    for (Variable variable = 0; variable < current.size(); ++variable) {
        const auto& values = model.domains()[variable].values();
        indices.push_back(static_cast<std::size_t>(std::find(values.begin(), values.end(), current[variable]) -
                                                   values.begin()));
        if (config.swap_proposals ? current.size() >= 2 : values.size() > 1) {
            mutable_variables.push_back(variable);
        }
    }
    if (mutable_variables.empty() || config.temperature_levels == 0 || config.sweeps_per_level == 0) {
        return result;
    }
    const auto moves_per_level = detail::checked_count(mutable_variables.size(), config.sweeps_per_level);
    detail::checked_count(moves_per_level, config.temperature_levels);
    for (std::size_t level = 0; level < config.temperature_levels; ++level) {
        const double temperature = detail::temperature_unchecked(config.initial_temperature,
            config.final_temperature, config.temperature_levels, level);
        for (std::size_t move = 0; move < moves_per_level; ++move) {
            const auto [variable, proposed_index] = config.swap_proposals
                ? detail::draw_swap(current.size(), random)
                : detail::draw_proposal(model, mutable_variables, indices, random);
            const auto& values = model.domains()[variable].values();
            Assignment proposed = current;
            if (config.swap_proposals) std::swap(proposed[variable], proposed[proposed_index]);
            else proposed[variable] = values[proposed_index];
            const double proposed_energy = full_energy(model, proposed);
            const double delta = proposed_energy - current_energy;
            if (!std::isfinite(delta)) {
                throw std::overflow_error("energy difference is nonfinite");
            }
            ++result.attempted_moves;
            // The top 53 bits give an exactly representable uniform in [0,1).
            const double uniform = delta > 0.0 ? static_cast<double>(random() >> 11) * 0x1.0p-53 : 0.0;
            if (metropolis_accept(delta, temperature, uniform)) {
                current = std::move(proposed);
                current_energy = proposed_energy;
                if (config.swap_proposals) std::swap(indices[variable], indices[proposed_index]);
                else indices[variable] = proposed_index;
                ++result.accepted_moves;
                if (current_energy < result.best_energy) {
                    result.best_assignment = current;
                    result.best_energy = current_energy;
                }
            }
        }
    }
    result.final_assignment = std::move(current);
    result.final_energy = current_energy;
    return result;
}

namespace detail {
template<bool History = false, bool Profile = false>
class ExecutionState {
    using Clock = std::chrono::steady_clock;
    const PreparedModel& prepared_;
    bool swaps_;
    InstrumentedResult* instrumentation;
    Clock::time_point start{};
    std::mt19937_64 random;
    std::unique_ptr<PreparedState> state_;
    SolverResult result;
    std::vector<Variable> mutable_variables;

    template<class F>
    decltype(auto) timed(std::size_t interval, F&& action) {
        if constexpr (Profile) {
            struct Guard {
                ProfileInterval& item;
                Clock::time_point begin;
                ~Guard() { item.seconds += std::chrono::duration<double>(Clock::now()-begin).count(); ++item.calls; }
            } guard{instrumentation->profile[interval], Clock::now()};
            return action();
        } else {
            (void)interval;
            return action();
        }
    }
    void record(const SolverResult& value) {
        if constexpr (History) instrumentation->history.push_back({
            std::chrono::duration<double>(Clock::now()-start).count(),
            value.attempted_moves, value.accepted_moves, value.best_energy});
    }
public:
    // The prepared model and optional instrumentation must outlive this chain.
    ExecutionState(const PreparedModel& prepared, const SolverConfig& config,
                   InstrumentedResult* output = nullptr)
        : prepared_(prepared), swaps_(config.swap_proposals), instrumentation(output), random(config.seed) {
        if constexpr (History || Profile) {
            if (!instrumentation) throw std::invalid_argument("instrumentation output is required");
        }
        if constexpr (History) start = Clock::now();
        detail::validate_temperatures(config.initial_temperature, config.final_temperature);
        const auto& model = prepared.source();
        detail::validate_swap_domains(model, config.swap_proposals);

        Assignment initial;
        if (config.initial_assignment) {
            initial = *config.initial_assignment;
        } else {
            initial.reserve(model.domains().size());
            for (const auto& domain : model.domains()) {
                std::uniform_int_distribution<std::size_t> choose(0, domain.values().size() - 1);
                initial.push_back(domain.values()[choose(random)]);
            }
        }
        state_ = std::make_unique<PreparedState>(prepared, initial);
        const double initial_energy = full_energy(model, initial);
        result = SolverResult{initial, initial_energy, initial, initial_energy, 0, 0, config.seed};
        record(result);

        for (Variable i = 0; i < model.domains().size(); ++i) {
            if (config.swap_proposals ? model.domains().size() >= 2 : model.domains()[i].values().size() > 1) {
                mutable_variables.push_back(i);
            }
        }
    }

    std::size_t attempts_per_sweep() const noexcept { return mutable_variables.size(); }
    const Assignment& assignment() const noexcept { return state_->assignment(); }
    double energy() const noexcept { return state_->energy(); }

    void advance_attempts(std::size_t attempts, double temperature) {
        detail::validate_temperatures(temperature, temperature);
        if (mutable_variables.empty()) return;
        if (attempts > std::numeric_limits<std::size_t>::max() - result.attempted_moves)
            throw std::overflow_error("proposal budget exceeds the move counter range");
        auto& state = *state_;
        const auto& model = prepared_.source();
        for (std::size_t move = 0; move < attempts; ++move) {
            const auto proposal = timed(0, [&] { return swaps_ ? detail::draw_swap(model.domains().size(), random)
                    : detail::draw_proposal(model, mutable_variables, state.indices(), random); });
            const double delta = timed(1, [&] { return swaps_ ? state.propose_swap(proposal.variable, proposal.index)
                    : state.propose(proposal.variable, proposal.index); });
            ++result.attempted_moves;
            const bool accepted = timed(2, [&] {
                const double uniform = delta > 0 ? static_cast<double>(random() >> 11) * 0x1.0p-53 : 0.0;
                return metropolis_accept(delta, temperature, uniform);
            });
            if (accepted) {
                timed(3, [&] { state.accept(); });
                ++result.accepted_moves;
                // Correct drift on a fixed accepted-move cadence.
                if (result.accepted_moves % 256 == 0) {
                    timed(4, [&] { state.rebase(); });
                }
                if (state.energy() < result.best_energy) {
                    bool improved = false;
                    timed(5, [&] {
                        const double exact = full_energy(model, state.assignment());
                        if (exact < result.best_energy) {
                            std::copy(state.assignment().begin(), state.assignment().end(), result.best_assignment.begin());
                            result.best_energy = exact;
                            improved = true;
                        }
                    });
                    if (improved) record(result);
                }
            }
        }
    }
    void advance_sweeps(std::size_t sweeps, double temperature) {
        advance_attempts(detail::checked_count(sweeps, attempts_per_sweep()), temperature);
    }

    // Construct and verify before committing, so failed imports preserve the chain.
    void install_assignment(const Assignment& assignment) {
        auto replacement = std::make_unique<PreparedState>(prepared_, assignment);
        const double exact = full_energy(prepared_.source(), assignment);
        if (exact < result.best_energy) {
            // Allocate a history entry before changing the chain, if requested.
            if constexpr (History) instrumentation->history.push_back({
                std::chrono::duration<double>(Clock::now()-start).count(),
                result.attempted_moves, result.accepted_moves, exact});
            std::copy(assignment.begin(), assignment.end(), result.best_assignment.begin());
            result.best_energy = exact;
        }
        state_.swap(replacement); // Also discards any rejected pending proposal.
    }

    SolverResult snapshot() const {
        auto output = result;
        output.final_assignment = state_->assignment();
        output.final_energy = full_energy(prepared_.source(), output.final_assignment);
        output.best_energy = full_energy(prepared_.source(), output.best_assignment);
        return output;
    }
    SolverResult finish() const {
        auto output = snapshot();
        if constexpr (History) instrumentation->seconds = std::chrono::duration<double>(Clock::now()-start).count();
        return output;
    }
};

template<bool History, bool Profile>
inline SolverResult solve_loop(const PreparedModel& prepared, const SolverConfig& config,
                               InstrumentedResult* instrumentation = nullptr) {
    ExecutionState<History, Profile> chain(prepared, config, instrumentation);
    if (chain.attempts_per_sweep() && config.temperature_levels && config.sweeps_per_level) {
        const auto moves = detail::checked_count(chain.attempts_per_sweep(), config.sweeps_per_level);
        detail::checked_count(moves, config.temperature_levels);
        for (std::size_t level = 0; level < config.temperature_levels; ++level)
            chain.advance_attempts(moves, detail::temperature_unchecked(config.initial_temperature,
                config.final_temperature, config.temperature_levels, level));
    }
    return chain.finish();
}

} // namespace detail

inline SolverResult solve(const PreparedModel& prepared, const SolverConfig& config = {}) {
    return detail::solve_loop<false, false>(prepared, config);
}
inline InstrumentedResult solve_instrumented(const PreparedModel& prepared,
        const SolverConfig& config = {}, const InstrumentationConfig& instrumentation_config = {}) {
    InstrumentedResult output;
    output.result = instrumentation_config.profile
        ? detail::solve_loop<true, true>(prepared, config, &output)
        : detail::solve_loop<true, false>(prepared, config, &output);
    return output;
}

inline SolverResult solve(const Model& model, const SolverConfig& config = {}) {
    return solve(prepare(model), config);
}

} // namespace annealer
