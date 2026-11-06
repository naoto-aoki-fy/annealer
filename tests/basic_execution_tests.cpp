#include "../examples/basic_execution.hpp"

int main(int argc, char** argv) {
    basic::Session session(argc, argv);
    for (int stage = 0; stage < 5; ++stage) {
        const bool fail_here = session.rank == session.size - 1;
        basic::Options options;
        options.solver.temperature_levels = 1;
        options.solver.sweeps_per_level = 1;
        char invalid[] = "--unknown";
        char value[] = "value";
        char* arguments[] = {argv[0], invalid, value};
        const int count = stage == 0 && fail_here ? 3 : 1;
        const auto make_model = [&] {
            if (stage == 1 && fail_here) throw std::runtime_error("construction failure");
            annealer::ModelBuilder builder;
            const auto x = builder.add_variable(annealer::Domain::cardinality(2));
            builder.add_term(1, {{x, [=](annealer::Value v) {
                if ((stage == 2 || stage == 3) && fail_here)
                    throw std::runtime_error("factor failure");
                return static_cast<double>(v);
            }}});
            return builder.build();
        };
        const int status = basic::run(count, arguments, options, make_model,
            [](const auto&) {}, [&](const auto&) {
                if (stage == 4) throw std::runtime_error("result processing failure");
            }, {.lookup_bytes = stage == 3 ? 0u : 1024u});
        int local = status == 1, total = 0;
        MPI_Allreduce(&local, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if (total != session.size) return 1;
    }
}
