# G-set G1 Maximum Cut example

## Overview

Maximum Cut partitions an undirected graph into two sets to maximize the sum
of weights on edges crossing the partition. The
[G-set collection hosted by Stanford](https://web.stanford.edu/~yyye/yyye/Gset/)
provides widely used graph optimization benchmarks. Its G1 instance has 800
vertices and 19,176 edges, all with weight one. This is a useful sparse,
nontrivial workload for comparing solution quality and incident-term evaluation
cost across fixed solver configurations.

The [example source](../../examples/gset_g1.cpp) uses only public library APIs;
its [shared helper](../../examples/gset_maxcut.hpp) implements parsing,
formulation, and independent cut evaluation. No core solver changes are needed.

## Mathematical Formulation

Let `x_i` be a binary variable indicating a vertex's partition. Define the pure
single-variable function `s_i = spin(x_i) = 1 - 2*x_i`. Equal partitions give
` s_i*s_j = 1`; different partitions give `s_i*s_j = -1`. Thus

```text
C(x) = sum_(i,j) w_ij * (1 - s_i*s_j)/2
E(x) = -C(x) = -0.5 * sum_(i,j) w_ij
                   + sum_(i,j) (0.5*w_ij) * spin(x_i)*spin(x_j)
```

The builder creates one `Domain::cardinality(2)` variable per vertex, sets the
constant offset, and adds one factorized pairwise term per edge:

```cpp
builder.add_term(0.5 * e.weight, {{e.u, spin}, {e.v, spin}});
```

There are no auxiliary variables. Each term is a coefficient times a product of
two single-variable functions, as required by the [model API](../model.md).
Minimizing energy maximizes cut weight, and minimum energy equals the negative
maximum cut weight up to floating-point rounding.

## Downloading the Benchmark

Run all commands from the repository root. Download outside the repository;
redistribution terms have not been verified and the data is not included.

```sh
curl --fail --location https://web.stanford.edu/~yyye/yyye/Gset/G1 -o /tmp/annealer-G1
```

The parser expects whitespace-separated `vertices edges`, followed by exactly
that many `u v weight` triples with **1-based** endpoints. It rejects zero vertex
counts, invalid numbers, missing edges, out-of-range endpoints, self-loops,
nonfinite weights, and trailing data. Zero edges are allowed. Repeated edges
are retained as separate weighted contributions; negative weights are supported.
Comments are not part of this input format.

## Building the Example

Requires CMake 3.20+, C++20, and a C compiler for the existing examples.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 1
ctest --test-dir build --output-on-failure --parallel 1
ctest --test-dir build -R '^gset_tests$' --output-on-failure
```

CMake links `gset_g1_example` to `annealer::annealer`. Use normal floating-point
optimization without fast-math.

## Running the Example

```sh
./build/gset_g1_example --help
./build/gset_g1_example /tmp/annealer-G1 --seed 0 --trials 10 --initial-temperature 30 --final-temperature 0.05 --levels 200 --sweeps 100
```

Syntax: `gset_g1_example FILE [OPTION VALUE ...]`.

| Option | Default | Meaning |
| --- | --- | --- |
| `--seed` | 0 | First unsigned 64-bit seed; trial k uses seed+k |
| `--trials` | 10 | Positive number of independent sequential trials |
| `--initial-temperature` | 30 | Positive finite initial temperature |
| `--final-temperature` | 0.05 | Positive finite final temperature, no greater than initial |
| `--levels` | 200 | Number of geometric temperature levels |
| `--sweeps` | 100 | Sweeps per level |
| `--reference-cut` | absent | Optional finite comparison value for this input |
| `--help` | — | Print syntax and defaults, used alone |

Zero levels or sweeps evaluate random initial states without proposing moves.
One level uses the initial temperature. Seed range overflow is rejected.
Preparation runs once with a 1 MiB unary-table payload budget, reused by all
trials. Each trial initializes a fresh random assignment.

## Understanding the Results

Each trial prints its seed, independently computed cut weight, best (minimum
observed) energy, attempted moves, accepted moves, solve seconds, and successful
energy verification. The check uses
`abs(cut + energy) <= 1e-10 * (1 + max(abs(cut), abs(energy)))`.
A mismatch terminates with an error. Small synthetic tests exhaustively verify
both reference and prepared energies without requiring solver optimality.

The final summary reports the best cut and energy across trials and the seed
that produced them; ties retain the first trial. A best-known solution is a
published feasible cut, which need not be a mathematically proven optimum.
The program embeds no best-known value. If you provide a separately verified
value using `--reference-cut`, the printed reference gap is that value minus the
best observed cut. It is a best-known gap only if the supplied reference actually
is a best-known solution for the same graph. Negative gaps indicate an improved
cut or an inappropriate reference. No global-optimality claim follows.

One sweep means `vertices` random single-variable attempts with replacement,
not one guaranteed visit per vertex. Attempt counts are
`levels * sweeps * vertices` per trial; accepted counts include equal-energy
moves. Counts describe the entire trial, not just the path to its saved best.

Solve seconds include initialization, Metropolis moves, reference checkpoints,
and result construction, excluding preparation and independent cut verification.
Total seconds include model construction, preparation, all trials, verification,
and reporting, excluding file parsing. Timings use `steady_clock`, vary with
host load and build settings, and have no correctness threshold.

Identical input, configuration, seed, and build/runtime environment reproduce
a trial's solution and counts, not its elapsed time. Cross-compiler or
cross-standard-library bitwise reproducibility is not promised. Record the input,
compiler, build flags, machine, and all printed settings when comparing runs.

## Performance and Limitations

The solver uses serial CPU-based Metropolis simulated annealing with geometric
cooling. There is no CUDA or GPU support. Higher initial temperatures permit more
uphill moves; lower final temperatures reduce them. More levels and sweeps cost
more time and may improve solution quality. Multiple seeds expose stochastic
variability but provide no guarantee of global optimality.

[Incremental evaluation](../preparation.md) updates only incident terms and uses
unary lookup tables within the preparation budget. It periodically rebases and
checks best states with the reference evaluator; returned energies are fully
recomputed. Dense graphs increase incident-term costs. Floating-point rounding
and poorly scaled signed weights can affect trajectories and verification;
finite input arithmetic can still overflow. G1's unit weights are exactly
representable at this scale. See the [solver documentation](../solver.md) for
schedule, budget, and reproducibility details.
