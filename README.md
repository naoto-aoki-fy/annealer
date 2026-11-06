# annealer

A C++20 library for finite-domain optimization with factorized higher-order
interactions. This first stage provides model construction and a serial reference
energy evaluator. The simulated annealing solver is not implemented yet.

Requires CMake 3.20+, a C++20 compiler, and a C compiler for the interoperability
example. No third-party libraries are required.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 1
ctest --test-dir build --output-on-failure --parallel 1
./build/multivalued_example
```

The example prints `Energy: -18.5`. Link the CMake target
`annealer::annealer` and include `<annealer/model.hpp>` in your own C++ code.

- [Model API, numerical behavior, and assumptions](docs/model.md)
- [Compiling a user model and C factor functions](docs/user-models.md)
