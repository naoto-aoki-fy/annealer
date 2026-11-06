#pragma once
#include "jsplib_jobshop.hpp"
#include <annealer/model.hpp>

namespace jobshop {
struct ModelConfig { std::size_t memory_bytes = std::size_t{1024} * 1024 * 1024; };
struct Estimate { Time horizon, penalty; std::size_t variables, terms, memory_bytes; };
inline std::size_t checked_add(std::size_t a, std::size_t b) {
    if (b > std::numeric_limits<std::size_t>::max() - a) throw std::overflow_error("Model size overflow");
    return a+b;
}
inline std::size_t checked_multiply(std::size_t a, std::size_t b) {
    if (b && a > std::numeric_limits<std::size_t>::max()/b) throw std::overflow_error("Model size overflow");
    return a*b;
}
inline Estimate estimate(const Instance& instance) {
    validate(instance);
    Time h = 0; for (const auto& op : instance.operations) h = completion(h, op.duration);
    constexpr Time exact = Time{1} << 53;
    if (h >= exact || static_cast<std::uintmax_t>(h) >= std::numeric_limits<std::size_t>::max())
        throw std::invalid_argument("Job-shop horizon exceeds supported range");
    std::size_t terms = 1, constraints = 0, domain_entries = static_cast<std::size_t>(h)+1;
    for (std::size_t i = 0; i < instance.operations.size(); ++i) {
        const auto count = static_cast<std::size_t>(h-instance.operations[i].duration)+1;
        domain_entries = checked_add(domain_entries, count);
        terms = checked_add(terms, count); ++constraints; // Makespan bound.
        if (instance.operations[i].step) {
            terms = checked_add(terms, static_cast<std::size_t>(h-instance.operations[i-1].duration)+1); ++constraints;
        }
        for (std::size_t j = i+1; j < instance.operations.size(); ++j)
            if (instance.operations[i].machine == instance.operations[j].machine) {
                terms = checked_add(terms, count); ++constraints;
            }
    }
    if (constraints > static_cast<std::uint64_t>(exact-h)/static_cast<std::uint64_t>(h+1))
        throw std::invalid_argument("Penalized energy exceeds exact integer range");
    // Conservative planning allowance for source + prepared ownership, factors,
    // callback captures, incidence, term caches, domains and proposal scratch.
    const auto bytes = checked_add(checked_multiply(terms, 2048), checked_multiply(domain_entries, 32));
    return {h, h+1, checked_add(instance.operations.size(), 1), terms, bytes};
}
inline annealer::Model make_model(const Instance& instance, const ModelConfig& config = {}) {
    const auto plan = estimate(instance);
    if (plan.memory_bytes > config.memory_bytes) throw std::invalid_argument("Estimated model memory exceeds limit; increase --memory-mib explicitly");
    annealer::ModelBuilder builder;
    for (const auto& op : instance.operations)
        builder.add_variable(annealer::Domain::cardinality(static_cast<std::size_t>(plan.horizon-op.duration)+1));
    const auto makespan = builder.add_variable(annealer::Domain::cardinality(static_cast<std::size_t>(plan.horizon)+1));
    builder.add_term(1, {{makespan, [](annealer::Value x) { return static_cast<double>(x); }}});
    auto add_relation = [&](std::size_t a, std::size_t b, auto violation) {
        const auto last = plan.horizon-instance.operations[a].duration;
        for (Time start = 0; start <= last; ++start)
            builder.add_term(static_cast<double>(plan.penalty), {
                {a, [start](annealer::Value x) { return x == start ? 1.0 : 0.0; }},
                {b, [start, violation](annealer::Value x) { return violation(start, x) ? 1.0 : 0.0; }}});
    };
    for (std::size_t i = 0; i < instance.operations.size(); ++i) {
        const auto duration = instance.operations[i].duration;
        add_relation(i, makespan, [duration](Time start, Time m) { return start+duration > m; });
        if (instance.operations[i].step) {
            const auto previous_duration = instance.operations[i-1].duration;
            add_relation(i-1, i, [previous_duration](Time start, Time next) { return start+previous_duration > next; });
        }
        for (std::size_t j = i+1; j < instance.operations.size(); ++j)
            if (instance.operations[i].machine == instance.operations[j].machine) {
                const auto other_duration = instance.operations[j].duration;
                add_relation(i, j, [duration, other_duration](Time start, Time other) {
                    return start < other+other_duration && other < start+duration;
                });
            }
    }
    return builder.build();
}
struct ModelVerification { Verification schedule; std::size_t makespan_violations; Time declared_makespan; double energy; bool feasible() const { return schedule.feasible() && !makespan_violations; } };
inline ModelVerification verify_model(const Instance& instance, const annealer::Assignment& x) {
    const auto plan = estimate(instance);
    if (x.size() != plan.variables) throw std::invalid_argument("Model assignment size mismatch");
    Schedule starts(x.begin(), x.end()-1);
    for (std::size_t i = 0; i < starts.size(); ++i)
        if (starts[i] < 0 || starts[i] > plan.horizon-instance.operations[i].duration) throw std::invalid_argument("Start outside model horizon");
    if (x.back() < 0 || x.back() > plan.horizon) throw std::invalid_argument("Makespan outside horizon");
    const auto result = evaluate(instance, starts);
    std::size_t bounds = 0;
    for (std::size_t i = 0; i < starts.size(); ++i) bounds += completion(starts[i], instance.operations[i].duration) > x.back();
    return {result, bounds, x.back(), static_cast<double>(x.back()) + static_cast<double>(plan.penalty)*
        static_cast<double>(result.precedence_violations+result.machine_conflicts+bounds)};
}
inline annealer::Assignment initial_assignment(const Instance& instance) {
    auto x = serial_schedule(instance);
    x.push_back(evaluate(instance, x).makespan);
    return x;
}
} // namespace jobshop
