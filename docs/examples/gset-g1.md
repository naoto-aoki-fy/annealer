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
formulation, and independent cut evaluation.

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

Set `ANNEALER_WORK_DIR` to an absolute directory of your choice before running
these commands. It defaults to `TMPDIR` when set, otherwise `/tmp`. Run the
following blocks in the same shell.

```sh
export ANNEALER_WORK_DIR="${ANNEALER_WORK_DIR:-${TMPDIR:-/tmp}}"
mkdir -p "${ANNEALER_WORK_DIR}"
curl --fail --location https://web.stanford.edu/~yyye/yyye/Gset/G1 -o "${ANNEALER_WORK_DIR}/annealer-G1"
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
cmake --build build --parallel "$(nproc)"
ctest --test-dir build --output-on-failure --parallel 1
ctest --test-dir build -R '^gset_tests$' --output-on-failure
```

CMake links `gset_g1_example` to `annealer::annealer` (and `annealer::mpi` in MPI builds). Use normal floating-point
optimization without fast-math.

## Running the Example

```sh
./build/gset_g1_example --help
./build/gset_g1_example "${ANNEALER_WORK_DIR}/annealer-G1" --seed 0 --trials 10 --initial-temperature 30 --final-temperature 0.05 --levels 200 --sweeps 100
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

## MPI execution

Build with `-DANNEALER_ENABLE_MPI=ON`, then run, for example:

```sh
mpiexec -n 2 ./build/gset_g1_example graph.gset --mode mpi-multistart --trials 10
mpiexec -n 2 ./build/gset_g1_example graph.gset --mode mpi-replica-exchange --trials 2 --exchange-interval 10 --exchange-seed 0
```

The default mode is `serial`. Multi-start distributes global trials with seeds
`seed + trial`. Replica exchange runs collective trials with solver base seed
`seed + trial * ranks` and exchange seed `exchange-seed + trial * ranks`;
replicas add their rank to the solver base seed. It uses a cold-to-hot geometric
temperature ladder and `levels * sweeps` total sweeps per replica per trial.
Exchange options apply only to replica exchange.

MPI modes independently verify local results and global best assignments.
Rank 0 prints one best-cut summary across all runs and MPI wall timing including
solving and verification. `--reference-cut` uses this globally selected cut.
`--trace-csv`, `--profile`, and `--target-cut` currently require serial mode:
the MPI APIs do not return the histories or profiling data these options need.
Requests combining any of these options with MPI fail before opening output
files or entering solver collectives. Distributed instrumentation is not provided.

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

## TTS and detailed measurement

`--target-cut W` reports the first history record with `-best_energy >= W`;
it prints elapsed solve seconds and attempted moves, or explicitly reports
non-achievement. A target present initially has zero attempted moves and the
actual initialization time. It does not shorten a trial. `--reference-cut`
continues to control only the displayed gap.

`--trace-csv PATH` writes initial and strict best updates after each trial,
including levels, sweeps, temperatures, seed, seconds, attempted and accepted
moves, energy, and cut. History is buffered during exploration. Existing files
at output paths are overwritten.

`--profile` compares detailed profiling, normal execution, and history-only
execution for each trial with identical configuration and seed. It checks
assignments, energies, and move counts for equality and independently validates
returned best and final cuts. It reports times, absolute overhead seconds, six
disjoint interval totals/call counts, and the residual described in the
[solver API](../solver.md#instrumented-prepared-solves).

```sh
cmake -S . -B "${ANNEALER_WORK_DIR}/annealer-build" -DCMAKE_BUILD_TYPE=Release \
  -DANNEALER_BUILD_BENCHMARKS=ON \
  -DCMAKE_CXX_FLAGS='-Wall -Wextra -Wpedantic -Wconversion'
cmake --build "${ANNEALER_WORK_DIR}/annealer-build" --parallel "$(nproc)"
ctest --test-dir "${ANNEALER_WORK_DIR}/annealer-build" --output-on-failure
"${ANNEALER_WORK_DIR}/annealer-build/gset_g1_example" "${ANNEALER_WORK_DIR}/annealer-G1" --trials 1 --levels 200 --sweeps 100 \
  --target-cut 11624 --profile --trace-csv docs/measurements/g1-profile.csv
"${ANNEALER_WORK_DIR}/annealer-build/gset_benchmark" "${ANNEALER_WORK_DIR}/annealer-G1" --levels 50,100,200 --budget 240 \
  --target-cut 11624 --csv docs/measurements/g1-budget.csv
```

## Equal-budget multi-start comparison

The optional `gset_benchmark` target shares parsing, model construction, cut
validation, and history helpers with the example. Defaults are levels 50,100,200,
100 sweeps, temperatures 30 to 0.05, 240 seconds per setting, and seeds increasing
from zero independently for each setting. One comparison is performed per
setting. There is no default target; explicitly pass 11624 for this experiment.
CLI overrides include `--levels` (comma-separated), `--budget` (positive finite
seconds), `--seed`, `--sweeps`, both temperatures, `--target-cut`, and `--csv`.

The common input loading/model preparation time is recorded separately.
Per-setting wall time includes restart initialization, history management,
exploration, and independent validation of best and final assignments. No result
file or console output occurs inside the restart loop. Each started iteration
runs its full update budget; an iteration that finishes after the wall deadline
is marked incomplete. Its history through the deadline still contributes to
the best cut and first target time, but it contributes neither to completed
iteration counts nor to their success rate. Updates exactly at the deadline
are included. No target hit means an empty CSV time cell, never zero.

The main CSV contains raw histories, including post-budget updates flagged
`within_budget=0`. Add `offset_seconds` to a history's `seconds` to locate it
within a setting. Sibling files append these suffixes:

- `.iterations.csv`: per-seed start/end, completion flag, solve time, TTS,
  attempted moves at target, final best cut, move counts, and attempted moves/s.
  Per-iteration TTS and best cut describe the entire iteration, including an
  overrun; use the summary for budget-limited results.
- `.summary.csv`: budget-limited best cut, first overall target time, completed
  count, successes, descriptive success percentage, elapsed time, and overrun.
  Success percentage is empty when no iteration completes.
- `.metadata.csv`: compiler, flags, build, CPU (on Linux), input path/dimensions,
  seed, target, schedule, loading/preparation time, and result output time.

The success fraction from this single comparison is descriptive. First overall
target time comes directly from timestamped histories; it is not estimated
from the average TTS among successful iterations. Detailed per-move profiling
is deliberately disabled for these comparisons. See the
[recorded measurement report](../measurements/g1-results.md) and its CSV links.

## Common execution commands

Build with MPI disabled or enabled as shown in [getting started](../getting-started.md).
Replace `graph.gset` with a local input in the format described above; MPI inputs
must be identical and accessible on every rank.

```sh
./build-mpi/gset_g1_example graph.gset --mode serial --seed 42 --trials 2 --levels 2 --sweeps 1
mpiexec -n 3 ./build-mpi/gset_g1_example graph.gset --mode mpi-multistart --seed 42 --trials 2 --levels 2 --sweeps 1
mpiexec -n 3 ./build-mpi/gset_g1_example graph.gset --mode mpi-replica-exchange --seed 42 --trials 2 --levels 2 --sweeps 1 --exchange-interval 1 --exchange-seed 17
```

These short budgets demonstrate the interface, without a solution-quality
promise. See [mode semantics and reporting](../mpi.md#ordinary-example-command-list)
for global trials, seed offsets, temperatures, sweep budgets and timing scope.

`--trace-csv`, `--profile`, and `--target-cut` require serial mode; MPI requests
are rejected before trace files are created or truncated.
