#pragma once
#include "gset_maxcut.hpp"
#include <annealer/solver.hpp>
#include <iomanip>
#include <ostream>

namespace gset {
inline double real_number(const std::string& text) {
    std::size_t end = 0;
    const double value = std::stod(text, &end);
    if (end != text.size() || !std::isfinite(value)) throw std::invalid_argument("Invalid number: " + text);
    return value;
}
inline void verify(const Graph& graph, const annealer::SolverResult& result) {
    for (const auto& pair : {std::pair{&result.best_assignment, result.best_energy},
                            std::pair{&result.final_assignment, result.final_energy}}) {
        const double cut = cut_weight(graph, *pair.first);
        if (std::abs(cut + pair.second) > 1e-10 * (1 + std::max(std::abs(cut), std::abs(pair.second))))
            throw std::runtime_error("Independent cut/energy verification failed");
    }
}
inline std::optional<annealer::BestUpdate> target_hit(
        const std::vector<annealer::BestUpdate>& history, double target,
        double limit = std::numeric_limits<double>::infinity()) {
    for (const auto& entry : history)
        if (entry.seconds <= limit && -entry.best_energy >= target) return entry;
    return std::nullopt;
}
inline double budget_best(const std::vector<annealer::BestUpdate>& history, double limit) {
    double best = -std::numeric_limits<double>::infinity();
    for (const auto& entry : history)
        if (entry.seconds <= limit) best = std::max(best, -entry.best_energy);
    return best;
}
inline void csv_header(std::ostream& out) {
    out << "kind,levels,sweeps,initial_temperature,final_temperature,seed,seconds,attempted,accepted,energy,cut,offset_seconds,within_budget\n";
}
inline void history_csv(std::ostream& out, const annealer::SolverConfig& config,
        const std::vector<annealer::BestUpdate>& history, double offset = 0,
        double budget = std::numeric_limits<double>::infinity()) {
    out << std::setprecision(17);
    for (const auto& h : history)
        out << "history," << config.temperature_levels << ',' << config.sweeps_per_level << ','
            << config.initial_temperature << ',' << config.final_temperature << ',' << config.seed << ','
            << h.seconds << ',' << h.attempted_moves << ',' << h.accepted_moves << ','
            << h.best_energy << ',' << -h.best_energy << ',' << offset << ','
            << (offset + h.seconds <= budget) << '\n';
}
inline bool same_result(const annealer::SolverResult& a, const annealer::SolverResult& b) {
    return a.best_assignment == b.best_assignment && a.final_assignment == b.final_assignment &&
        a.best_energy == b.best_energy && a.final_energy == b.final_energy &&
        a.attempted_moves == b.attempted_moves && a.accepted_moves == b.accepted_moves && a.seed == b.seed;
}
}
