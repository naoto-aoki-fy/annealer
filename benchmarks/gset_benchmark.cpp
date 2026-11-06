#include "../examples/gset_measurement.hpp"
#include "build_info.hpp"
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
using Clock = std::chrono::steady_clock;
double elapsed(Clock::time_point start) { return std::chrono::duration<double>(Clock::now()-start).count(); }
std::string quote(const std::string& text) {
    std::string out = "\"";
    for (char c : text) { out += c; if (c == '"') out += c; }
    return out + '"';
}
int main(int argc, char** argv) {
 try {
    if (argc == 2 && std::string(argv[1]) == "--help") {
        std::cout << "Usage: gset_benchmark FILE [--levels 50,100,200] [--budget 240] [--seed 0]\n"
            " [--sweeps 100] [--initial-temperature 30] [--final-temperature 0.05]\n"
            " [--target-cut W] [--csv PATH]\nOne comparison per setting; no detailed profiling.\n";
        return 0;
    }
    if (argc < 2) throw std::invalid_argument("Missing input; use --help");
    annealer::SolverConfig config;
    config.initial_temperature=30; config.final_temperature=0.05; config.sweeps_per_level=100;
    std::vector<std::size_t> levels{50,100,200};
    double budget=240;
    std::optional<double> target;
    std::string path="gset-results.csv";
    for (int i=2; i<argc; ++i) {
        const std::string option=argv[i];
        if (++i==argc) throw std::invalid_argument("Missing value for "+option);
        const std::string value=argv[i];
        if (option=="--budget") budget=gset::real_number(value);
        else if (option=="--target-cut") target=gset::real_number(value);
        else if (option=="--seed") config.seed=gset::unsigned_number<std::uint64_t>(value);
        else if (option=="--sweeps") config.sweeps_per_level=gset::unsigned_number<std::size_t>(value);
        else if (option=="--initial-temperature") config.initial_temperature=gset::real_number(value);
        else if (option=="--final-temperature") config.final_temperature=gset::real_number(value);
        else if (option=="--csv") path=value;
        else if (option=="--levels") {
            levels.clear(); std::istringstream list(value); std::string token;
            while (std::getline(list,token,',')) levels.push_back(gset::unsigned_number<std::size_t>(token));
            if (levels.empty() || value.back()==',') throw std::invalid_argument("Empty level list item");
        } else throw std::invalid_argument("Unknown option: "+option);
    }
    if (budget<=0) throw std::invalid_argument("Budget must be positive");
    annealer::detail::validate_temperatures(config.initial_temperature,config.final_temperature);
    const auto load_start=Clock::now();
    std::ifstream input(argv[1]);
    if (!input) throw std::runtime_error("Cannot open input");
    const auto graph=gset::read(input);
    const auto model=gset::make_model(graph);
    const auto prepared=annealer::prepare(model,{.lookup_bytes=1024*1024});
    const double preparation=elapsed(load_start);
    std::ofstream raw(path), iterations(path+".iterations.csv"), summary(path+".summary.csv"), metadata(path+".metadata.csv");
    for (auto* stream : {&raw,&iterations,&summary,&metadata}) {
        if (!*stream) throw std::runtime_error("Cannot open result file");
        stream->exceptions(std::ios::badbit|std::ios::failbit);
        *stream << std::setprecision(17);
    }
    gset::csv_header(raw);
    iterations << "levels,seed,start_seconds,end_seconds,solve_seconds,completed,tts_seconds,tts_moves,best_cut,attempted,accepted,moves_per_second\n";
    summary << "levels,budget_seconds,elapsed_seconds,overrun_seconds,best_cut_within_budget,first_target_seconds,completed,successful,success_percent\n";
    metadata << "key,value\ncompiler," << quote(BENCH_COMPILER) << "\nbuild," << quote(BENCH_BUILD)
        << "\nflags," << quote(BENCH_FLAGS) << "\ninput," << quote(argv[1]) << "\nvertices," << graph.vertices
        << "\nedges," << graph.edges.size() << "\nloading_preparation_seconds," << preparation
        << "\nbudget," << budget << "\nsweeps," << config.sweeps_per_level << "\ninitial_temperature," << config.initial_temperature
        << "\nfinal_temperature," << config.final_temperature << "\nfirst_seed," << config.seed << "\ncomparisons,1\ntarget_cut,";
    if (target) metadata << *target;
    metadata << '\n';
    std::ifstream cpu("/proc/cpuinfo"); std::string line;
    while (std::getline(cpu,line)) if (line.starts_with("model name")) { metadata << "cpu," << quote(line) << '\n'; break; }
    struct Iteration { double start,end; annealer::SolverConfig config; annealer::InstrumentedResult measured; };
    double output_seconds=0;
    const auto first_seed=config.seed;
    for (auto level : levels) {
        config.temperature_levels=level; config.seed=first_seed;
        std::vector<Iteration> runs;
        double best=-std::numeric_limits<double>::infinity();
        std::optional<double> first;
        std::size_t completed=0,successes=0;
        const auto begin=Clock::now();
        while (elapsed(begin)<budget) {
            const double offset=elapsed(begin);
            auto measured=annealer::solve_instrumented(prepared,config);
            gset::verify(graph,measured.result);
            const double end=elapsed(begin);
            best=std::max(best,gset::budget_best(measured.history,budget-offset));
            const auto hit=target ? gset::target_hit(measured.history,*target,budget-offset) : std::nullopt;
            if (hit && !first) first=offset+hit->seconds;
            if (end<=budget) { ++completed; if (hit) ++successes; }
            runs.push_back({offset,end,config,std::move(measured)});
            if (config.seed==std::numeric_limits<std::uint64_t>::max()) {
                if (elapsed(begin)<budget) throw std::overflow_error("Seed range exhausted");
                break;
            }
            ++config.seed;
        }
        const double total=elapsed(begin);
        const auto output_start=Clock::now();
        for (const auto& run : runs) {
            const auto& m=run.measured; const auto& r=m.result;
            gset::history_csv(raw,run.config,m.history,run.start,budget);
            const auto hit=target ? gset::target_hit(m.history,*target) : std::nullopt;
            iterations << level << ',' << r.seed << ',' << run.start << ',' << run.end << ',' << m.seconds << ',' << (run.end<=budget) << ',';
            if (hit) iterations << hit->seconds;
            iterations << ','; if (hit) iterations << hit->attempted_moves;
            iterations << ',' << -r.best_energy << ',' << r.attempted_moves << ',' << r.accepted_moves << ',' << static_cast<double>(r.attempted_moves)/m.seconds << '\n';
            std::cout << "Levels: " << level << " Seed: " << r.seed << " Completed: " << (run.end<=budget)
                << " Best cut: " << -r.best_energy << " Moves/s: " << static_cast<double>(r.attempted_moves)/m.seconds << " TTS: ";
            if (hit) std::cout << hit->seconds; else std::cout << "not reached";
            std::cout << " Verified: yes\n";
        }
        summary << level << ',' << budget << ',' << total << ',' << std::max(0.0,total-budget) << ',';
        if (std::isfinite(best)) summary << best;
        summary << ','; if (first) summary << *first;
        summary << ',' << completed << ',' << successes << ',';
        if (completed) summary << 100.0*static_cast<double>(successes)/static_cast<double>(completed);
        summary << '\n';
        std::cout << "Summary levels: " << level << " Best within budget: " << best << " First target: ";
        if (first) std::cout << *first; else std::cout << "not reached";
        std::cout << " Completed: " << completed << " Successful: " << successes << " Success percent: ";
        if (completed) std::cout << 100.0*static_cast<double>(successes)/static_cast<double>(completed); else std::cout << "n/a";
        std::cout << " Overrun seconds: " << std::max(0.0,total-budget) << '\n';
        raw.flush(); iterations.flush(); summary.flush();
        output_seconds+=elapsed(output_start);
    }
    metadata << "result_output_seconds," << output_seconds << '\n';
    std::cout << "Loading/preparation seconds: " << preparation << " Result output seconds: " << output_seconds << '\n';
 } catch (const std::exception& e) { std::cerr << "Error: " << e.what() << '\n'; return 1; }
}
