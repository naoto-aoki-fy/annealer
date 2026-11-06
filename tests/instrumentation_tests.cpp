#include "../examples/gset_measurement.hpp"
#include <sstream>
#include <iostream>
void require(bool ok) { if (!ok) throw std::runtime_error("Instrumentation check failed"); }
int main() {
 try {
    std::istringstream input("2 1\n1 2 3\n");
    const auto graph=gset::read(input);
    const auto model=gset::make_model(graph);
    const auto prepared=annealer::prepare(model);
    annealer::SolverConfig c;
    c.initial_assignment=annealer::Assignment{0,0};
    c.temperature_levels=1; c.sweeps_per_level=300;
    c.initial_temperature=c.final_temperature=0.001;
    const auto normal=annealer::solve(prepared,c);
    const auto history=annealer::solve_instrumented(prepared,c);
    const auto profile=annealer::solve_instrumented(prepared,c,{.profile=true});
    require(gset::same_result(normal,history.result) && gset::same_result(normal,profile.result));
    gset::verify(graph,profile.result);
    require(history.history.size()==2);
    require(history.history[0].attempted_moves==0 && history.history[0].accepted_moves==0);
    require(history.history[0].seconds>0);
    const auto hit=gset::target_hit(history.history,3);
    require(hit && hit->attempted_moves==1 && hit->accepted_moves==1);
    require(!gset::target_hit(history.history,4));
    for (std::size_t i=1;i<history.history.size();++i) {
        require(history.history[i].seconds>=history.history[i-1].seconds);
        require(history.history[i].attempted_moves>history.history[i-1].attempted_moves);
        require(history.history[i].accepted_moves>history.history[i-1].accepted_moves);
        require(history.history[i].best_energy<history.history[i-1].best_energy);
    }
    for (int i=0;i<3;++i) require(profile.profile[i].calls==normal.attempted_moves);
    require(profile.profile[3].calls==normal.accepted_moves);
    require(profile.profile[4].calls==normal.accepted_moves/256);
    require(profile.profile[5].calls==1);
    double sum=0; for (const auto& p : profile.profile) sum+=p.seconds;
    require(sum<=profile.seconds);
    for (const auto& p : history.profile) require(p.calls==0 && p.seconds==0);
    c.sweeps_per_level=0; c.initial_assignment=annealer::Assignment{0,1};
    const auto zero=annealer::solve_instrumented(prepared,c);
    require(zero.history.size()==1 && zero.result.attempted_moves==0);
    require(gset::target_hit(zero.history,3)->attempted_moves==0);
    std::istringstream flat_input("2 0\n");
    const auto flat_graph=gset::read(flat_input);
    const auto flat_model=gset::make_model(flat_graph);
    const auto flat_prepared=annealer::prepare(flat_model);
    c.sweeps_per_level=300;
    const auto flat=annealer::solve_instrumented(flat_prepared,c,{.profile=true});
    require(flat.history.size()==1 && flat.result.accepted_moves==600);
    require(flat.profile[4].calls==2 && flat.profile[5].calls==0);
    std::vector<annealer::BestUpdate> artificial{{0.5,0,0,0},{1,1,1,-2},{1.5,2,2,-3}};
    require(gset::budget_best(artificial,1)==2);
    require(gset::target_hit(artificial,2,1)->attempted_moves==1);
    require(!gset::target_hit(artificial,3,1));
    require(!gset::target_hit(artificial,0,0.49));
    std::ostringstream csv;
    gset::csv_header(csv); gset::history_csv(csv,c,artificial,239,240);
    require(csv.str().find("239,1\n")!=std::string::npos);
    require(csv.str().find("239,0\n")!=std::string::npos);
    // Exercise several seeds and a nontrivial schedule for trajectory parity.
    c.initial_assignment.reset(); c.temperature_levels=5; c.sweeps_per_level=100;
    for (c.seed=0;c.seed<10;++c.seed) {
        const auto a=annealer::solve(prepared,c);
        require(gset::same_result(a,annealer::solve_instrumented(prepared,c).result));
        require(gset::same_result(a,annealer::solve_instrumented(prepared,c,{.profile=true}).result));
    }
    std::cout << "Instrumentation, TTS, profile, CSV and budget boundary tests passed\n";
 } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
