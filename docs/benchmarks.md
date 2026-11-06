# Serial benchmarks and verification

## Reproduce the measurements

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DANNEALER_BUILD_BENCHMARKS=ON \
    -DCMAKE_CXX_FLAGS='-Wall -Wextra -Wpedantic -Werror' \
    -DCMAKE_C_FLAGS='-Wall -Wextra -Wpedantic -Werror'
cmake --build build --parallel 1
ctest --test-dir build --output-on-failure --parallel 1
./build/serial_benchmark
```

Run the executable three times sequentially to reproduce the sampling procedure.
It prints compiler/build/CPU metadata and CSV measurements. No benchmark is a
CTest test, and timings have no pass/fail performance threshold. Optional IPO is
described in [preparation](preparation.md); it was disabled for the results below.

The workload generator uses model seed 314159, input seed 271828, and solver seed
42. Terms select variables with replacement, allowing repeated occurrences.
Coefficients alternate 1.25 and -0.75; factors use `sin(0.2*x) + 0.25`.

| Workload | Variables | Terms | Domain sizes | Term order | Mean incident terms per variable |
| --- | ---: | ---: | --- | ---: | ---: |
| tiny | 8 | 12 | all 2 (`{0,1}`) | 3 | 3.875 |
| sparse | 128 | 256 | alternating 3 and 5 | 5 | 9.859375 |
| dense | 32 | 1024 | alternating 3 and 5 | 8 | 229.46875 |

The multivalued domains are `{-3,0,4}` and `{-3,-1,1,3,5}`. Sparsity is reported
as actual deduplicated incidence counts, rather than a fraction of all possible
interaction tuples. No quadratic encoding or Cartesian tabulation is used.

## Measurement boundaries

- Model construction and deterministic input generation occur outside timings.
- Preparation is one `prepare` call, including source copying, indexing, allocation,
  and tabulation. Both the default table budget and zero-budget fallback are run.
- Full evaluation uses the same 1,000 predetermined assignments for both paths.
  The optimized measurement calls the public `prepared_energy` convenience API,
  which constructs a temporary state and allocates its caches each time.
- Proposal throughput uses the same 10,000 predetermined legal moves from the
  same initial assignment. Both paths accept every move for comparability. The
  reference path performs a full recomputation; the optimized path uses a reused
  `PreparedState`. Input generation and state initialization are excluded. This
  stream is an evaluator comparison, not an annealing experiment; it omits periodic
  rebasing and best tracking. Checksums and final states are checked for agreement.
- Annealing timings use actual Metropolis acceptance, 20 levels, two sweeps per
  level, default temperatures 10 to 0.1, and a supplied initial assignment. They
  include solver initialization, drift correction, best tracking, and returned
  energy recomputation. Optimized annealing reuses the prepared model. The
  preparation-plus-solve column adds preparation to this time for one solve.
  Reference and optimized floating-point trajectories can differ.
- `table_bytes` is logical unary double payload. `cache_bytes` is current term
  contributions plus largest incident-term proposal scratch, in bytes. It excludes
  assignment/index arrays, model copies, vector capacity, and allocator overhead.
  Neither is a total resident-memory measurement.

Every timing uses `std::chrono::steady_clock`. Results are consumed in a printed
checksum so the measured work cannot be discarded as unused. Executions and
workloads run sequentially on one application thread.

## Measured results, 2026-11-06

Environment: Intel Core i7-6700K @ 4.00GHz; GNU C++ 11.5.0; C++20;
`-Wall -Wextra -Wpedantic -Werror -O3 -DNDEBUG`; Release, IPO off, no fast-math.
The following are medians from three sequential local runs. Raw output with all
metadata and checksums: [run 1](benchmarks/run1.csv), [run 2](benchmarks/run2.csv),
[run 3](benchmarks/run3.csv). Measurements are illustrative and depend on host
load, clocks, compiler, and callback cost; no CPU affinity was imposed.

| Workload/path | Prep ms | Reference eval/s | Prepared eval/s | Reference proposals/s | Prepared proposals/s | Proposal speedup |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| tiny/tables | 0.021 | 1,707,866 | 2,563,544 | 1,840,542 | 13,482,886 | 7.33x |
| sparse/tables | 0.174 | 46,848 | 131,801 | 53,274 | 3,333,633 | 62.58x |
| dense/tables | 1.148 | 8,051 | 32,361 | 8,121 | 93,520 | 11.52x |
| tiny/direct | 0.006 | 1,682,723 | 1,149,922 | 1,855,077 | 4,170,450 | 2.25x |
| sparse/direct | 0.110 | 46,637 | 36,527 | 53,104 | 879,242 | 16.56x |
| dense/direct | 0.502 | 8,058 | 6,592 | 8,104 | 26,650 | 3.29x |

| Workload/path | Reference solve ms | Prepared solve ms | Prep + solve ms | Table bytes | Cache bytes |
| --- | ---: | ---: | ---: | ---: | ---: |
| tiny/tables | 0.196 | 0.042 | 0.063 | 576 | 144 |
| sparse/tables | 97.589 | 5.957 | 6.130 | 40,608 | 2,192 |
| dense/tables | 157.353 | 20.028 | 21.161 | 261,968 | 10,200 |
| tiny/direct | 0.195 | 0.097 | 0.103 | 0 | 144 |
| sparse/direct | 97.511 | 10.561 | 10.677 | 0 | 2,192 |
| dense/direct | 157.561 | 54.112 | 54.609 | 0 | 10,200 |

Tables improved full evaluation throughput by 1.50–4.02x in these workloads.
Direct fallback full evaluations regressed to 0.68–0.82x reference throughput:
the convenience API adds indexing/cache allocation and prepared metadata access
without removing callback cost. Incident updates still gained 2.25–16.56x without
tables. Dense interactions touch many terms, reducing the gain from sparsity.
Table preparation costs more time and memory; these runs do not establish that
tabulation pays off for every model, short run, or very large domain.

## Integration and execution review

The tests enumerate the tiny binary, heterogeneous multivalued, nonlinear, and
constant-only examples and their legal single-variable moves, across full, partial,
and absent tables. They check independent formulas, reference energies, deltas,
accepted/rejected caches, and returned assignment energies. Exact reference minima
are -3, -156.5, -1710, and 1 respectively. The existing tests cover invalid inputs,
nonfinite callbacks/arithmetic, schedule boundaries, reproducibility, singleton
domains, zero factors/coefficients, repeated variables, and long update sequences.
Differential tolerance is documented in [preparation](preparation.md).

Correctness assertions do not require the annealer to find an optimum. Separate
binary/nonlinear example experiments used seeds 0–4 and observed zero best-energy
gap on all five runs for each model. These tiny results are not an optimization
quality guarantee for larger models.

Source/CMake dependency review found ordinary serial loops and standard C/C++
library dependencies, with no thread creation, parallel execution policies, GPU,
OpenMP, MPI, or worker-process execution paths in preparation or solving. The
local benchmark links only the normal C/C++ runtime and math libraries. CPU-only,
single-thread execution remains a library requirement.
