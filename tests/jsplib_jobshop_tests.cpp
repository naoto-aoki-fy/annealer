#include "../examples/jsplib_jobshop.hpp"
#include <iostream>
#include <sstream>
void require(bool ok) { if (!ok) throw std::runtime_error("Job-shop check failed"); }
jobshop::Instance parse(const std::string& text) { std::istringstream in(text); return jobshop::read(in); }
template<class F> void rejects(F f) { bool bad = false; try { f(); } catch (const std::exception&) { bad = true; } require(bad); }
void dimensions(std::size_t jobs, std::size_t machines) {
    std::ostringstream text; text << jobs << ' ' << machines << '\n';
    for (std::size_t j = 0; j < jobs; ++j) {
        for (std::size_t s = 0; s < machines; ++s) text << (s+j)%machines << " 1 ";
        text << '\n';
    }
    const auto instance = parse(text.str());
    require(instance.jobs == jobs && instance.machines == machines && instance.operations.size() == jobs*machines);
    for (std::size_t i = 0; i < instance.operations.size(); ++i) {
        const auto& op = instance.operations[i];
        require(op.job == i/machines && op.step == i%machines && op.machine == (op.job+op.step)%machines && op.duration == 1);
    }
    const auto starts = jobshop::serial_schedule(instance);
    require(starts == jobshop::serial_schedule(instance));
    const auto result = jobshop::evaluate(instance, starts);
    require(result.feasible() && result.makespan == static_cast<jobshop::Time>(jobs*machines));
}
int main() {
    try {
        const auto instance = parse("2 2\n0 3 1 2\n1 4 0 1\n");
        auto result = jobshop::evaluate(instance, {0, 4, 0, 4});
        require(result.feasible() && result.makespan == 6);
        result = jobshop::evaluate(instance, {0, 3, 0, 4});
        require(result.precedence_violations == 0 && result.machine_conflicts == 1 && result.makespan == 5);
        result = jobshop::evaluate(instance, {0, 1, 3, 7});
        require(result.precedence_violations == 1 && result.machine_conflicts == 0 && result.makespan == 8);
        result = jobshop::evaluate(instance, {0, 1, 0, 0});
        require(result.precedence_violations == 2 && result.machine_conflicts == 2 && result.makespan == 4);
        const auto serial = jobshop::serial_schedule(instance);
        require(serial == jobshop::Schedule({0, 3, 5, 9}));
        require(jobshop::evaluate(instance, serial).feasible() && jobshop::evaluate(instance, serial).makespan == 10);
        dimensions(6, 6); dimensions(10, 5); // ft06/la01 shapes; synthetic data only.
        for (std::size_t j = 1; j <= 8; ++j) for (std::size_t m = 1; m <= 7; ++m) dimensions(j, m);
        const auto repeat = parse("1 2 0 2 0 3");
        require(jobshop::evaluate(repeat, jobshop::serial_schedule(repeat)).feasible());
        const auto boundary = parse("1 1 0 9223372036854775807");
        require(jobshop::evaluate(boundary, {0}).makespan == std::numeric_limits<jobshop::Time>::max());
        rejects([&] { (void)jobshop::evaluate(boundary, {1}); });
        rejects([&] { (void)jobshop::evaluate(instance, {-1, 3, 5, 9}); });
        rejects([&] { (void)jobshop::evaluate(instance, {}); });
        rejects([&] { (void)jobshop::serial_schedule({0, 0, {}}); });
        auto broken = instance; broken.operations[0].step = 1;
        rejects([&] { (void)jobshop::evaluate(broken, serial); });
        for (const auto* text : {"", "0 1", "1 0", "-1 2", "1x 2", "1 1", "1 1 0", "1 1 1 2",
            "1 1 -1 2", "1 1 0 0", "1 1 0 -1", "1 1 0 1.5", "1 1 0 nan", "1 1 0 2 extra",
            "2 1 0 2", "18446744073709551616 1", "18446744073709551615 2",
            "1 1 0 9223372036854775808", "1 2 0 9223372036854775807 1 1"})
            rejects([&] { (void)parse(text); });
        std::cout << "JSPLib parser, schedule validation and serial constructor tests passed\n";
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
