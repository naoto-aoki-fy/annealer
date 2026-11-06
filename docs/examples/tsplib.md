# TSPLIB95 eil51 and berlin52

`tsplib_example` represents a tour with one integer variable per tour position,
`x_t` in `{0,...,N-1}`. Its original cyclic length includes the closing edge:

`L(x) = sum_t d(x_t, x_((t+1) mod N))`.

Each distance uses unary-factor products:
`d(x_t,x_next) = sum_a [x_t=a] * d(a,x_next)`. The distance-row callbacks
share an owned immutable matrix; no callback receives two variable values.
There are `N^2` distance terms and `N*N*(N-1)/2` all-different terms:

`E(x) = L(x) + P * sum_(t<u) sum_a [x_t=a]*[x_u=a]`.

Let `D` be the maximum rounded distance and `U=N*D`. The penalty is `P=U+1`.
All lengths lie in `[0,U]`, so every infeasible assignment has energy greater
than every valid tour. Global energy minima solve the TSP, but simulated
annealing does not guarantee finding them. An independent evaluator recomputes
cyclic length from coordinates and checks all city IDs and duplicate pairs.

## Supported input and numerical limits

The [TSPLIB95 specification](https://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/tsp/tsp95.pdf)
defines EUC_2D distances as `floor(sqrt(dx*dx+dy*dy)+0.5)`, computed separately
for each edge in double precision. Rounding the total Euclidean length instead
would give a different objective. Near rounding boundaries, floating point
coordinate parsing and arithmetic can affect the result.

The parser requires `TYPE: TSP`, `DIMENSION: N` with `N>=2`,
`EDGE_WEIGHT_TYPE: EUC_2D`, then `NODE_COORD_SECTION` and one `id x y` record
per city. IDs are 1-based, unique and within the dimension; record order may
vary. Coordinates must be finite. `NAME` and repeatable `COMMENT` headers are
accepted. Header colons and the final `EOF` marker are optional. Blank lines
are ignored. Unknown headers/sections, unsupported types (including ATT, GEO,
CEIL_2D and explicit distance matrices), duplicate headers, extra fields,
missing cities and data after EOF are rejected.

Dimension arithmetic is checked. Nonfinite distance computations and conservative
bounds exceeding `U + P*N*(N-1)/2 <= 2^53` are rejected so integer model energies
are exact in `double` and cost arithmetic cannot overflow. Huge finite coordinates
may overflow intermediate squared distances and are rejected too. The supported
subset covers eil51 and berlin52; it is not a general TSPLIB parser.

## Download, build and run

The [official TSPLIB95 site](https://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/)
links compressed instance files. Download and decompress the targets locally:

Set `ANNEALER_WORK_DIR` to an absolute directory of your choice before running
these commands. It defaults to `TMPDIR` when set, otherwise `/tmp`. Run the
following blocks in the same shell.

```sh
export ANNEALER_WORK_DIR="${ANNEALER_WORK_DIR:-${TMPDIR:-/tmp}}"
mkdir -p "${ANNEALER_WORK_DIR}"
curl -fL https://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/tsp/eil51.tsp.gz -o "${ANNEALER_WORK_DIR}/eil51.tsp.gz"
curl -fL https://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/tsp/berlin52.tsp.gz -o "${ANNEALER_WORK_DIR}/berlin52.tsp.gz"
gzip -dk "${ANNEALER_WORK_DIR}/eil51.tsp.gz" "${ANNEALER_WORK_DIR}/berlin52.tsp.gz"
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel "$(nproc)"
./build/tsplib_example "${ANNEALER_WORK_DIR}/eil51.tsp" --reference-cost 426 --seed 42 --trials 10 \
  --initial-temperature 30 --final-temperature 0.05 --levels 200 --sweeps 100
./build/tsplib_example "${ANNEALER_WORK_DIR}/berlin52.tsp" --reference-cost 7542 --seed 42 --trials 10 \
  --initial-temperature 30 --final-temperature 0.05 --levels 200 --sweeps 100
ctest --test-dir build --output-on-failure -R tsplib
```

The [official optimum listing](https://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/tsp/TSP-BEST.html)
reports **426** for eil51 and **7542** for berlin52. Supply `--reference-cost`
only after verifying that the input is the corresponding original instance;
no reference is inferred from its filename. The absolute gap is
`tour_length-reference_cost`, printed only for a valid tour. Tests contain only
synthetic data and do not download real instances or require stochastic success.

## Execution and scaling

Each trial initializes an identity permutation and shuffles it with a separate
`mt19937_64` generator using that trial's seed. The solver uses
`swap_proposals=true` and the shuffled initial assignment. Seeds increment
consecutively, with overflow rejected. Replay is deterministic within the same
standard-library implementation. Swap proposals preserve the initial multiset;
there are `N` attempts per sweep. See [solver documentation](../solver.md).

Output reports variables, terms, schedule, seed, tour length, collisions, energy,
validity, optional reference gap and solve seconds. Best and final assignments
are checked independently. Feasible results are candidates with optimality
unproven. The summary identifies the lowest-energy trial. Exit status zero means
execution and verification succeeded. Total time includes model construction,
preparation, initialization, solves, verification and output, excluding parsing.
Record input checksums, build settings, schedules and seeds with measurements.

The penalty causes O(N^3) model terms even though the tour objective has O(N^2)
terms. Each swap reevaluates the union of incident terms, O(N^2) for the penalty.
The CLI caps factor lookup payload at 1 MiB; uncached factors use callbacks.
Model metadata, incidence lists, distance storage and mutable caches require
additional memory. This is a public factor-API example, not a specialized TSP
implementation with 2-opt or constant-time swap deltas. Start with short schedules
and tune temperatures for each instance's distance scale. Rotations and reversals
are not fixed, and all positions remain swappable.

## MPI execution

With `ANNEALER_ENABLE_MPI=ON`, select `--mode mpi-multistart` or
`--mode mpi-replica-exchange` (default: `serial`). For example:

```sh
mpiexec -n 2 ./build-mpi/tsplib_example instance.dat --mode mpi-multistart --seed 42 --trials 3
mpiexec -n 2 ./build-mpi/tsplib_example instance.dat --mode mpi-replica-exchange --seed 42 --trials 3 --exchange-interval 10 --exchange-seed 17
```

Multi-start counts global logical starts and uses one common permutation generated
from the requested base seed on every rank and for every start. Solver seeds are
base seed plus logical start ID. Serial trials retain their own permutation from
each trial seed; serial and MPI trajectories therefore need not match.

Replica exchange counts collective runs. Run `t` uses solver base seed
`base + t * ranks`, also its common permutation initialization seed, and exchange
seed `exchange_base + t * ranks`. The MPI API offsets solver seeds by rank;
input configurations and initial permutations are identical across ranks.
The cold-to-hot geometric ladder spans final to initial temperature, with
`levels * sweeps` sweeps per replica. All modes use swap proposals.

Only rank zero reports. The initialization base seed, winning solver seed and
multi-start ID or replica run and winning rank identify results. Both locally
returned best/final assignments and the global best are independently verified.
Reference gaps are printed only for feasible global candidates. Reproducibility
requires the same rank count, configuration, and standard-library implementation.
See [MPI execution](../mpi.md) for build and launcher details.

## Common execution commands

Build with MPI disabled or enabled as shown in [getting started](../getting-started.md).
Replace `instance.tsp` with a local input in the format described above; MPI inputs
must be identical and accessible on every rank.

```sh
./build-mpi/tsplib_example instance.tsp --mode serial --seed 42 --trials 2 --levels 2 --sweeps 1
mpiexec -n 3 ./build-mpi/tsplib_example instance.tsp --mode mpi-multistart --seed 42 --trials 2 --levels 2 --sweeps 1
mpiexec -n 3 ./build-mpi/tsplib_example instance.tsp --mode mpi-replica-exchange --seed 42 --trials 2 --levels 2 --sweeps 1 --exchange-interval 1 --exchange-seed 17
```

These short budgets demonstrate the interface, without a solution-quality
promise. See [mode semantics and reporting](../mpi.md#ordinary-example-command-list)
for global trials, seed offsets, temperatures, sweep budgets and timing scope.
