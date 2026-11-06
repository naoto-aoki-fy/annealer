# Compiling a user model

The [example model header](../examples/user_model.hpp) declares a factory returning
an owning `annealer::Model`. Its [C++ implementation](../examples/user_model.cpp)
registers three differently sized domains, a nonlinear square function, a C factor,
a cubic interaction, a unary term, and a constant term. The
[driver](../examples/main.cpp) evaluates `{-2, 2, 3}` and prints `Energy: -18.5`.

The [C header](../examples/c_factors.h) wraps its declaration in `extern "C"` for
C++ callers. Its [implementation](../examples/c_factors.c) compiles as C. The
library's model-building API is C++; ordinary C functions with compatible
`double(int64_t)` signatures can supply factors through normal interoperability.
There is no C-only model-building API.

## In another CMake project

```cmake
cmake_minimum_required(VERSION 3.20)
project(my_model LANGUAGES C CXX)
add_subdirectory(path/to/annealer)
add_executable(my_model main.cpp user_model.cpp c_factors.c)
target_link_libraries(my_model PRIVATE annealer::annealer)
```

The example and tests in the annealer subdirectory are also available. Set
`BUILD_TESTING=OFF` if you do not need the tests. There is no installed CMake
package at this stage.

## Direct compilation

From the annealer repository root, a separate C compilation preserves the
language boundary:

```sh
mkdir -p build-direct
cc -std=c99 -O2 -c examples/c_factors.c -o build-direct/c_factors.o
c++ -std=c++20 -O2 -Iinclude examples/main.cpp examples/user_model.cpp \
    build-direct/c_factors.o -o build-direct/multivalued_example
./build-direct/multivalued_example
```

Replace the example sources with your own factory and driver. Native compilation
lets users write ordinary functions rather than source strings or expressions.
See the [model contract](model.md) for purity, lifetime, and numerical rules.
