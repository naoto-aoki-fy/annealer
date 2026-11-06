#include <annealer/solver.hpp>
#include "build_info.hpp"

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

namespace {
using Clock = std::chrono::steady_clock;
template <class F> double seconds(F action) {
    const auto start = Clock::now();
    action();
    return std::chrono::duration<double>(Clock::now() - start).count();
}
double unary(annealer::Value x) { return std::sin(static_cast<double>(x) * 0.2) + 0.25; }
bool close(double a, double b) { return std::abs(a-b) <= 1e-8 * (1 + std::max(std::abs(a), std::abs(b))); }

annealer::Model make_model(std::size_t variables, std::size_t terms, std::size_t order, bool binary) {
    annealer::ModelBuilder builder;
    std::mt19937_64 random(314159);
    for (std::size_t i = 0; i < variables; ++i) {
        if (binary) { builder.add_variable(annealer::Domain::cardinality(2)); }
        else { builder.add_variable(i % 2 ? annealer::Domain{-3, -1, 1, 3, 5} : annealer::Domain{-3, 0, 4}); }
    }
    for (std::size_t term = 0; term < terms; ++term) {
        std::vector<annealer::Factor> factors;
        for (std::size_t i = 0; i < order; ++i) {
            factors.push_back(annealer::native_factor<unary>(random() % variables));
        }
        builder.add_term(term % 2 ? -0.75 : 1.25, std::move(factors));
    }
    return builder.build();
}

void run(const char* name, std::size_t n, std::size_t terms, std::size_t order, bool binary, std::size_t budget) {
    const auto model = make_model(n, terms, order, binary);
    // All inputs are generated before timing and reused for both implementations.
    std::mt19937_64 random(271828);
    std::vector<annealer::Assignment> assignments;
    for (std::size_t sample = 0; sample < 1000; ++sample) {
        annealer::Assignment assignment;
        for (const auto& domain : model.domains()) {
            assignment.push_back(domain.values()[random() % domain.values().size()]);
        }
        assignments.push_back(std::move(assignment));
    }
    std::vector<annealer::detail::Proposal> moves;
    std::vector<std::size_t> indices(n, 0);
    for (std::size_t attempt = 0; attempt < 10000; ++attempt) {
        const auto variable = random() % n;
        const auto count = model.domains()[variable].values().size();
        const auto replacement = annealer::detail::replacement_index(indices[variable], random() % (count-1));
        moves.push_back({variable, replacement});
        indices[variable] = replacement;
    }
    const auto start = Clock::now();
    const auto prepared = annealer::prepare(model, {budget});
    const double preparation = std::chrono::duration<double>(Clock::now() - start).count();
    double ref_sum = 0, opt_sum = 0;
    const double ref_eval = seconds([&] { for (const auto& x : assignments) { ref_sum += annealer::full_energy(model, x); } });
    const double opt_eval = seconds([&] { for (const auto& x : assignments) { opt_sum += annealer::prepared_energy(prepared, x); } });
    if (!close(ref_sum, opt_sum)) { throw std::runtime_error("evaluation checksum mismatch"); }
    annealer::Assignment initial;
    for (const auto& domain : model.domains()) { initial.push_back(domain.values().front()); }
    auto current = initial;
    double energy = annealer::full_energy(model, current);
    double ref_delta_sum = 0, opt_delta_sum = 0;
    const double ref_moves = seconds([&] {
        for (const auto& move : moves) {
            current[move.variable] = model.domains()[move.variable].values()[move.index];
            const double next = annealer::full_energy(model, current);
            ref_delta_sum += next - energy;
            energy = next;
        }
    });
    annealer::PreparedState state(prepared, initial);
    const double opt_moves = seconds([&] {
        for (const auto& move : moves) {
            opt_delta_sum += state.propose(move.variable, move.index);
            state.accept(); // Predetermined all-accept stream, not an annealing run.
        }
    });
    if (!close(ref_delta_sum, opt_delta_sum) || state.assignment() != current || !close(state.energy(), energy)) {
        throw std::runtime_error("proposal checksum mismatch");
    }
    annealer::SolverConfig config;
    config.seed = 42;
    config.initial_assignment = initial;
    config.temperature_levels = 20;
    config.sweeps_per_level = 2;
    double result_sum = 0;
    const double ref_solve = seconds([&] {
        const auto result = annealer::solve_reference(model, config);
        result_sum += result.best_energy + result.final_energy + result.accepted_moves;
    });
    const double opt_solve = seconds([&] {
        const auto result = annealer::solve(prepared, config);
        result_sum += result.best_energy + result.final_energy + result.accepted_moves;
    });
    std::size_t incidences = 0;
    for (std::size_t i = 0; i < n; ++i) { incidences += prepared.incident(i).size(); }
    std::cout << name << ',' << (budget ? "tables" : "direct") << ',' << n << ',' << terms << ','
              << (binary ? "2" : "3/5") << ',' << order << ',' << double(incidences)/n << ','
              << preparation*1000 << ',' << 1000/ref_eval << ',' << 1000/opt_eval << ','
              << 10000/ref_moves << ',' << 10000/opt_moves << ',' << ref_solve*1000 << ','
              << opt_solve*1000 << ',' << (preparation+opt_solve)*1000 << ','
              << prepared.lookup_bytes() << ',' << state.cache_bytes() << ','
              << ref_sum+opt_sum+ref_delta_sum+opt_delta_sum+result_sum << '\n';
}
} // namespace

int main() {
    try {
        std::cout << "# Compiler: " << BENCH_COMPILER << "\n# Flags: " << BENCH_FLAGS
                  << "\n# Build: " << BENCH_BUILD << "\n# CPU: ";
        std::ifstream cpu("/proc/cpuinfo");
        std::string line;
        bool found = false;
        while (std::getline(cpu, line)) {
            if (line.starts_with("model name")) { std::cout << line.substr(line.find(':')+1); found = true; break; }
        }
        if (!found) { std::cout << "unavailable (record CPU externally)"; }
        std::cout << "\n# Seeds: model=314159, inputs=271828, solver=42; serial steady_clock timings\n"
                  << "workload,path,variables,terms,domain_sizes,order,mean_incident_terms,prep_ms,ref_eval_s,opt_eval_s,ref_proposal_s,opt_proposal_s,ref_solve_ms,opt_solve_ms,prep_plus_solve_ms,table_bytes,cache_bytes,checksum\n"
                  << std::setprecision(9);
        for (const auto budget : {std::size_t{64*1024*1024}, std::size_t{0}}) {
            run("tiny", 8, 12, 3, true, budget);
            run("sparse", 128, 256, 5, false, budget);
            run("dense", 32, 1024, 8, false, budget);
        }
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
