# DIMACS graph coloring benchmark

[Executable](../../examples/dimacs_coloring.cpp),
[parser, formulation and verifier](../../examples/dimacs_coloring.hpp), and
[tests](../../tests/coloring_tests.cpp) use only public annealer APIs.

## Formulation

For a requested positive color count K, vertex v has one integer variable
`x_v` with `Domain::cardinality(K)`, containing `{0, ..., K-1}`. Define the pure
single-variable indicator `f_c(x) = (x == c ? 1 : 0)`.

```text
E(x) = sum_(u,v in edges) sum_(c=0..K-1) f_c(x_u) * f_c(x_v)
```

Each edge adds K coefficient-one products of two single-variable factors,
with the color captured by value. There are no auxiliary variables or opaque
multivariable callbacks. `ModelBuilder`, `prepare`, and `solve` are unchanged.
The independent verifier compares endpoint colors directly and checks assignment
length and color range. Energy equals the number of monochromatic edges.
Zero conflicts certifies a valid K-coloring. Positive conflicts do not prove
that the graph is not K-colorable.

## Obtain DSJC125.1 and convert the binary archive

The [official DIMACS archive](https://archive.dimacs.rutgers.edu/pub/challenge/graph/benchmarks/color/)
provides `DSJC125.1.col.b`, a **binary** adjacency bitmap with a text preamble.
Renaming it to `.col` does not convert it. This example accepts text only and
rejects the binary length record. Use the original `bin2asc` translator from
[this source mirror](https://github.com/BartMassey/instances/tree/83b4e11b5d83bcf44feb628d1e4f7a994de0a0cf/binformat).
The pinned source includes the binary format specification and translator.
The following procedure was executed successfully on Linux with GCC 11.5:

```sh
git clone https://github.com/BartMassey/instances /tmp/annealer-dimacs-instances
git -C /tmp/annealer-dimacs-instances checkout --detach 83b4e11b5d83bcf44feb628d1e4f7a994de0a0cf
cc -std=gnu89 /tmp/annealer-dimacs-instances/binformat/bin2asc.c \
  -o /tmp/annealer-dimacs-instances/binformat/bin2asc
curl --fail --location \
  https://archive.dimacs.rutgers.edu/pub/challenge/graph/benchmarks/color/DSJC125.1.col.b \
  -o /tmp/DSJC125.1.col.b
/tmp/annealer-dimacs-instances/binformat/bin2asc \
  /tmp/DSJC125.1.col.b /tmp/DSJC125.1.col
sha256sum /tmp/DSJC125.1.col.b /tmp/DSJC125.1.col
# Correct the archive preamble's doubled edge count for a simple graph.
sed 's/^p edge 125 1472$/p edge 125 736/' /tmp/DSJC125.1.col \
  > /tmp/DSJC125.1-simple.col
sha256sum /tmp/DSJC125.1-simple.col
```

Expected SHA-256 values, binary then converted text:

```text
d1fa485ae23e795500adde49e5904e3b9641c4d80351be538e873001a28b9346
5f796d4e7c395d5817e144e9f8a638c132688f0acea7418875d9070e0c0a7d7e
```

The translator preserves the archive header `p edge 125 1472`, but emits
736 distinct undirected edges. An independent bitmap decode confirms 736:
the preamble declares twice the actual undirected edge count. The explicit `sed` command above corrects
only this known header; the strict parser does not silently accept mismatches.
The normalized file has SHA-256:

```text
391a413b3eda8dc7efeedbcca75d4670c6d68e9846737941012abcb5f0782ba7
```

The normalized instance has 125 vertices and 736 edges. Downloads stay outside
the repository and are not needed by CTest. The legacy translator is an external
utility, not a project dependency; use it with the specified archive and verify
the hashes. Benchmark data is not redistributed here.

## Build and run

From the repository root, with CMake 3.20+, C++20 and a C compiler:

```sh
cmake -S . -B /tmp/annealer-coloring-build -DCMAKE_BUILD_TYPE=Release \
  -DANNEALER_BUILD_BENCHMARKS=ON \
  -DCMAKE_CXX_FLAGS='-Wall -Wextra -Wpedantic -Wconversion'
cmake --build /tmp/annealer-coloring-build --parallel 2
ctest --test-dir /tmp/annealer-coloring-build --output-on-failure
/tmp/annealer-coloring-build/dimacs_coloring_example --help
/tmp/annealer-coloring-build/dimacs_coloring_example /tmp/DSJC125.1-simple.col \
  --colors 5 --seed 0 --trials 10 --initial-temperature 30 \
  --final-temperature 0.05 --levels 200 --sweeps 100
```

The example is built by default, independently of `ANNEALER_BUILD_BENCHMARKS`.
For a quick smoke run, use `--trials 1 --levels 10 --sweeps 10`.

| Option | Default | Meaning |
| --- | --- | --- |
| `--colors` | 5 | Positive K fitting the integer domain; no optimality claim |
| `--seed` | 0 | First unsigned 64-bit seed |
| `--trials` | 10 | Positive number of sequential trials, using seed + trial |
| `--initial-temperature` | 30 | Positive finite starting temperature |
| `--final-temperature` | 0.05 | Positive finite ending temperature, no greater than initial |
| `--levels` | 200 | Geometric temperature levels; zero evaluates initialization only |
| `--sweeps` | 100 | Sweeps per level; zero evaluates initialization only |

One level uses the initial temperature. One sweep makes `vertices` random
single-variable attempts with replacement. Seed range overflow is rejected.
Preparation occurs once with a 1 MiB unary-table payload budget. Each trial
starts with a fresh random assignment. See [solver details](../solver.md).

## Input and output

Text input allows blank lines and `c` comment lines, exactly one `p edge N M`
header with positive N, and exactly M `e U V` records using 1-based endpoints.
Zero edges and isolated vertices are supported. Malformed/unknown records,
extra fields, missing or excess edges, negative or out-of-range indices,
self-loops, and duplicate undirected edges (including reversed duplicates) are
rejected. Other DIMACS variants, weights and binary input are unsupported.

Output records vertices, edges, colors, schedule and trial count. Each trial
reports seed, conflicting edges, best energy, solve runtime and either
`valid K-coloring` or `assignment contains conflicts`. Both best and final
assignments are independently verified against their returned energies; a
mismatch fails the run. The summary retains the first best seed on ties.
A successful process exit means the computation and verification succeeded;
inspect the coloring status to determine feasibility.

Solve seconds include initialization and solving, excluding independent
verification. Total seconds include model construction, preparation, trials,
verification and output, excluding parsing. Timings use `steady_clock`.
Record input hashes, compiler, flags, machine, parameters and seeds when
comparing runs. Fixed seeds replay assignments and move counts in the same
build/runtime environment; timings and cross-platform bitwise replay are not
promised.

## Tests and limitations

CTest requires no downloads. Exhaustive small-graph tests include triangles
with K=1,2,3, a bipartite path, empty edges, isolated vertices, and both cached
and uncached preparation. Every assignment compares independent conflict
counts with `full_energy` and `prepared_energy`. Tests also cover malformed
input, invalid colors/assignments, fixed-seed replay and CLI error/status output.

This is serial CPU Metropolis annealing, without a feasibility or optimality
guarantee or a specialized graph-coloring search. Increasing K, trials or the
schedule can change quality and cost. The model has N variables and K*M terms;
large K and dense graphs can consume substantial memory and time. Domains are
materialized, so very large representable K can still exceed available memory.
Use normal floating-point optimization without fast-math. Unit conflict
energies are exactly representable for this benchmark.

A validation run with GCC 11.5, Release, seed 0, K=5, one trial, 200 levels
and 100 sweeps returned best energy 1 and one independently verified
conflicting edge (2.583 solve seconds on the validation host). This assignment
contains conflicts; it is not a valid 5-coloring. This is a reproducibility
example, not a reference optimum or a runtime threshold.
