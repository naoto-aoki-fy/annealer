# Getting started

The library requires C++20, CMake 3.20+, and a C compiler for the C interoperability
example. The default serial build has no third-party runtime dependencies.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DANNEALER_ENABLE_MPI=OFF
cmake --build build --parallel 1
ctest --test-dir build --output-on-failure --parallel 1
./build/minimal_example
./build/binary_hubo_example
./build/multivalued_example
./build/nonlinear_example
```

The minimal example defines a model in one file. The binary example has a cubic
HUBO interaction. The multivalued example compiles a C factor together with a C++
model, using three heterogeneous integer domains. The nonlinear example has a
fifth-order factorized interaction, square functions, a repeated variable, and
a singleton domain. The binary and nonlinear programs enumerate their tiny state
spaces first, then print the best-energy gaps for five sequential seeded runs.
These are optimization-quality demonstrations, not guarantees of optimality.

## Complete minimal example

Save this standalone snippet as `minimal.cpp`. The repository
[CLI example](../examples/minimal.cpp) additionally uses shared execution helpers:

```cpp
#include <annealer/solver.hpp>
#include <iostream>

double square(std::int64_t value) {
    const double x = static_cast<double>(value);
    return x * x;
}

int main() {
    annealer::ModelBuilder builder;
    const auto x = builder.add_variable(annealer::Domain{-2, 1, 4});
    const auto y = builder.add_variable(annealer::Domain::cardinality(3));
    builder.set_constant(1.0);
    builder.add_term(-0.5, {{x, square}, {y, [](std::int64_t v) { return double(v) + 1; }}});
    const auto model = builder.build();
    const auto prepared = annealer::prepare(model, {.lookup_bytes = 1024});
    annealer::SolverConfig config;
    config.seed = 42;
    config.initial_assignment = annealer::Assignment{-2, 0};
    config.initial_temperature = 10;
    config.final_temperature = 0.1;
    config.temperature_levels = 100;
    config.sweeps_per_level = 10;
    const auto result = annealer::solve(prepared, config);
    std::cout << "Best energy: " << result.best_energy << '\n';
}
```

From the repository root, compile and run the saved file:

```sh
c++ -std=c++20 -O2 -Iinclude minimal.cpp -o minimal
./minimal
```

The checked local run prints `Best energy: -23`. This observed result is not a
general solver success guarantee. Callbacks receive actual `std::int64_t` domain
values, and coefficients and results use `double`. Write callbacks as deterministic,
side-effect-free functions. No variable binary encoding is required.

For another CMake project, use `add_subdirectory(path/to/annealer)` and link
`annealer::annealer`. See [user-model compilation](user-models.md) for separate
C/C++ source files and normal `extern "C"` interoperability.

## Configuration and limitations

- [Model API](model.md): domains, coefficients, arbitrary interaction orders,
  validation, and numerical errors.
- [Solver](solver.md): temperatures, sweeps, optional initial assignment, seeds,
  zero budgets, and reported best/final states.
- [Compilation and preparation](preparation.md): native function specialization,
  runtime tabulation, memory limits, direct fallback, and drift correction.
- [Benchmarks](benchmarks.md): reproducible serial workloads and measured results.

Set `PreparationConfig::lookup_bytes` to zero to disable unary tables, or provide
a payload budget. Native compilation and optional IPO optimize user functions;
preparation constructs tables and indices at runtime. Arbitrary callbacks need
not be constexpr. Local execution is CPU-only. Optional MPI distributes independent starts or
exchanges replicas; see [MPI execution](mpi.md).

Use normal floating-point optimization, without fast-math. Finite inputs can still
overflow during products, sums, or differences; the library throws as documented.
Underflow and rounding follow double arithmetic. Incremental rounding can affect
near-tie decisions; reference checkpoints correct drift and returned energies are
recomputed. Runs reproduce within the same model, configuration, seed, and
build/runtime environment; cross-toolchain bitwise reproducibility is not promised.

## Selectable execution

For an MPI build (requires an MPI installation):

```sh
cmake -S . -B build-mpi -DCMAKE_BUILD_TYPE=Release -DANNEALER_ENABLE_MPI=ON
cmake --build build-mpi --parallel 2
ctest --test-dir build-mpi --output-on-failure
```

An MPI-enabled executable still defaults to serial without a launcher. The
[complete command list and mode semantics](mpi.md#ordinary-example-command-list)
cover all eleven examples. MPI-disabled binaries reject requested MPI modes.
The saved snippet above needs only `-Iinclude`; compiling repository CLI sources
also requires their sibling helper headers. Keep `examples/` intact; use CMake
for MPI linkage and the `ANNEALER_ENABLE_MPI` definition.
