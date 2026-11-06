# annealer

A C++20 library for finite-domain optimization with factorized higher-order
interactions, with a serial CPU Metropolis simulated annealing solver and a
reference full-energy evaluator.

[Getting started and examples](docs/getting-started.md) · [Serial benchmarks](docs/benchmarks.md)

Requires CMake 3.20+, a C++20 compiler, and a C compiler for the interoperability
example. No third-party libraries are required.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 1
ctest --test-dir build --output-on-failure --parallel 1
./build/multivalued_example
```

The example evaluates an assignment (`Energy: -18.5`) and runs the solver. Link the CMake target
`annealer::annealer` and include `<annealer/model.hpp>` in your own C++ code.

- [Model API, numerical behavior, and assumptions](docs/model.md)
- [Compiling a user model and C factor functions](docs/user-models.md)
- [Solver configuration, schedule, and results](docs/solver.md)
- [Native compilation, runtime preparation, and incremental updates](docs/preparation.md)
- [G-set G1 Max-Cut tutorial and benchmark](docs/examples/gset-g1.md)
- [DIMACS graph coloring tutorial and DSJC125.1 benchmark](docs/examples/dimacs-coloring.md)
- [Weighted partial MaxSAT WCNF example and benchmarks](docs/examples/maxsat.md)
- [OR-Library set covering scp41 benchmark](docs/examples/set-covering.md)
- [QAPLIB nug12 and nug20 benchmark](docs/examples/qaplib.md)
- [TSPLIB95 eil51 and berlin52 benchmark](docs/examples/tsplib.md)
- [JSPLib job-shop model and schedule validation](docs/examples/jsplib-jobshop.md)
