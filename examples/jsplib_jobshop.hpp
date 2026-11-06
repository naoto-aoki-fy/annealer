#pragma once
#include <algorithm>
#include <charconv>
#include <cstdint>
#include <istream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace jobshop {
using Time = std::int64_t;
struct Operation { std::size_t job, step, machine; Time duration; };
struct Instance { std::size_t jobs, machines; std::vector<Operation> operations; };
using Schedule = std::vector<Time>; // Job-major operation order.
inline Time completion(Time start, Time duration) {
    if (start < 0 || duration <= 0) throw std::invalid_argument("Start must be nonnegative and duration positive");
    if (start > std::numeric_limits<Time>::max() - duration)
        throw std::overflow_error("Operation completion exceeds time range");
    return start + duration;
}
inline void validate(const Instance& instance) {
    if (!instance.jobs || !instance.machines || instance.jobs > instance.operations.max_size() / instance.machines ||
        instance.operations.size() != instance.jobs * instance.machines)
        throw std::invalid_argument("Invalid job-shop dimensions or operation count");
    Time total = 0;
    for (std::size_t i = 0; i < instance.operations.size(); ++i) {
        const auto& op = instance.operations[i];
        if (op.job != i / instance.machines || op.step != i % instance.machines || op.machine >= instance.machines)
            throw std::invalid_argument("Invalid job-shop operation metadata");
        total = completion(total, op.duration);
    }
}
template<class T> inline T unsigned_number(const std::string& text) {
    T value{};
    const auto parsed = std::from_chars(text.data(), text.data()+text.size(), value);
    if (text.empty() || text.front() == '-' || text.front() == '+' || parsed.ec != std::errc{} || parsed.ptr != text.data()+text.size())
        throw std::invalid_argument("Invalid nonnegative integer: " + text);
    return value;
}
inline Instance read(std::istream& input) {
    auto token = [&]() {
        std::string value;
        if (!(input >> value)) throw std::invalid_argument("Incomplete JSPLib input");
        return value;
    };
    const auto jobs = unsigned_number<std::size_t>(token());
    const auto machines = unsigned_number<std::size_t>(token());
    Instance instance{jobs, machines, {}};
    if (!jobs || !machines || jobs > instance.operations.max_size() / machines)
        throw std::invalid_argument("Invalid job-shop dimensions");
    Time total = 0;
    for (std::size_t j = 0; j < jobs; ++j) {
        for (std::size_t step = 0; step < machines; ++step) {
            const auto machine = unsigned_number<std::size_t>(token());
            const auto duration = unsigned_number<Time>(token());
            if (machine >= machines || duration <= 0) throw std::invalid_argument("Invalid machine or duration");
            total = completion(total, duration);
            instance.operations.push_back({j, step, machine, duration});
        }
    }
    input >> std::ws;
    if (input.bad()) throw std::runtime_error("JSPLib input read failed");
    if (!input.eof()) throw std::invalid_argument("Unexpected trailing JSPLib data");
    return instance;
}
struct Verification {
    std::size_t precedence_violations = 0, machine_conflicts = 0;
    Time makespan = 0;
    bool feasible() const noexcept { return precedence_violations == 0 && machine_conflicts == 0; }
};
// Independent schedule validation: no annealer or energy model is involved.
inline Verification evaluate(const Instance& instance, const Schedule& starts) {
    validate(instance);
    if (starts.size() != instance.operations.size()) throw std::invalid_argument("Schedule size mismatch");
    std::vector<Time> ends;
    Verification result;
    for (std::size_t i = 0; i < starts.size(); ++i) {
        ends.push_back(completion(starts[i], instance.operations[i].duration));
        result.makespan = std::max(result.makespan, ends.back());
        if (instance.operations[i].step && ends[i-1] > starts[i]) ++result.precedence_violations;
    }
    for (std::size_t i = 0; i < starts.size(); ++i)
        for (std::size_t j = i + 1; j < starts.size(); ++j)
            if (instance.operations[i].machine == instance.operations[j].machine &&
                starts[i] < ends[j] && starts[j] < ends[i]) {
                if (result.machine_conflicts == std::numeric_limits<std::size_t>::max())
                    throw std::overflow_error("Machine conflict counter overflow");
                ++result.machine_conflicts;
            }
    return result;
}
// Run all operations consecutively in job-major order: feasible, not optimized.
inline Schedule serial_schedule(const Instance& instance) {
    validate(instance);
    Schedule starts;
    Time clock = 0;
    for (const auto& op : instance.operations) {
        starts.push_back(clock);
        clock = completion(clock, op.duration);
    }
    return starts;
}
} // namespace jobshop
