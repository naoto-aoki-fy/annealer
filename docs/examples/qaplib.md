# QAPLIB nug12 and nug20

`qaplib_example` uses integer facility variables `x_i` with domain `{0,...,N-1}`.
A feasible assignment is a permutation of locations. The original objective is
`Q(x) = sum_i sum_j A[i][j] * B[x_i][x_j]`, including diagonal entries and both
ordered directions. The independent evaluator computes this double sum directly,
counts equal-location pairs, and reports permutation validity and penalty cost.

## Factor decomposition and penalty

The model uses only pure unary callbacks:

`B[x_i][x_j] = sum_a [x_i=a] * B[a][x_j]`.

Each nonzero off-diagonal flow entry produces up to `N` two-factor terms.
Zero flow entries and zero rows of B are skipped. If B is symmetric, the two
ordered flow entries combine as `A[i][j]+A[j][i]` for `i<j`. Diagonal flow entries
use a single unary function `B[x_i][x_i]`. Callbacks share an owned immutable B
matrix snapshot and otherwise capture only fixed indices.

The all-different penalty is
`P * sum_{i<j} sum_a [x_i=a]*[x_j=a]`.
Let `U = sum(A)*max(B)` and `P = U+1`. For nonnegative matrices every QAP cost
is in `[0,U]`, so any infeasible assignment has energy at least `P`, greater
than every feasible energy. Global energy minima therefore solve the original
QAP. Swap optimization starts feasible and never activates the penalty, but the
model and independent evaluator also define energies for infeasible assignments.

The supported format is a positive size followed by exactly two row-major square
matrices of nonnegative decimal integers, separated by arbitrary whitespace.
Negative, fractional or nonfinite entries, incomplete matrices and trailing data
are rejected. Dimensions and integer operations are checked. The conservative
bound `U + P*N*(N-1)/2 <= 2^53`, together with checked input and sum bounds,
keeps model energies exact in `double` and avoids integer overflow. Large weights
or sizes can be rejected even when an individual entry fits in 64 bits. Signed
and floating point QAP variants are outside this example's supported subset.

## Acquire and run

The [official QAPLIB listing](https://coral.ise.lehigh.edu/data-sets/qaplib/qaplib-problem-instances-and-solutions/)
links the original Nug12 and Nug20 `.dat` files and their solution files. Download
those individual files through the listing into a local directory. The listing
reports optimal objectives **578** for nug12 and **2570** for nug20. These apply
only to the exact original data, matrix orientation and objective convention.
No real instance is bundled or required by CTest.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel "$(nproc)"
./build/qaplib_example /path/to/nug12.dat --reference-cost 578 --seed 42 --trials 10 \
  --initial-temperature 30 --final-temperature 0.05 --levels 200 --sweeps 100
./build/qaplib_example /path/to/nug20.dat --reference-cost 2570 --seed 42 --trials 10 \
  --initial-temperature 30 --final-temperature 0.05 --levels 200 --sweeps 100
ctest --test-dir build --output-on-failure -R qaplib
```

Each trial shuffles the identity permutation using `mt19937_64` seeded with the
printed trial seed, then calls the solver with `swap_proposals=true` and that
initial assignment. Seeds increase consecutively, with overflow rejected. Replay
is reproducible within the same standard-library implementation. One sweep
attempts `N` swaps. Schedules may need tuning; there is no guarantee of reaching
an optimum.

Output includes variables, model terms, schedule, seed, original QAP cost,
collisions, total energy, feasibility and solve seconds, followed by the
lowest-energy trial and total time. Best and final assignments are independently
verified. `--reference-cost` is optional and must be supplied only after checking
the reference against the exact input. Its absolute objective gap is
`candidate_cost-reference_cost`, printed only for a feasible candidate; a gap is
not an optimality proof. No reference is inferred from a filename. Exit status
zero means execution and verification succeeded. Total time includes model
construction, preparation, permutation initialization, solves, verification and
output, excluding parsing. Record input checksums, build settings and schedules.

## Computational limits

The decomposition and collision penalty use O(N^3) terms in dense cases, each
with at most two factors. Lookup payload can grow as O(N^4) when all factors
are tabulated; this CLI limits lookup tables to 1 MiB, with the remainder using
callbacks. Model ownership, term metadata and incidence/cache storage are
additional. Each swap recomputes the union of terms incident to either facility,
which is O(N^2) terms for dense instances. This is a demonstration of the public
factor API, not a specialized QAP delta implementation. nug20 costs more than
nug12 in both memory and time; begin with short schedules. See
[swap proposals](../solver.md) and [preparation](../preparation.md).

## MPI execution

With `ANNEALER_ENABLE_MPI=ON`, select `--mode mpi-multistart` or
`--mode mpi-replica-exchange` (default: `serial`). For example:

```sh
mpiexec -n 2 ./build-mpi/qaplib_example instance.dat --mode mpi-multistart --seed 42 --trials 3
mpiexec -n 2 ./build-mpi/qaplib_example instance.dat --mode mpi-replica-exchange --seed 42 --trials 3 --exchange-interval 10 --exchange-seed 17
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
Replace `instance.dat` with a local input in the format described above; MPI inputs
must be identical and accessible on every rank.

```sh
./build-mpi/qaplib_example instance.dat --mode serial --seed 42 --trials 2 --levels 2 --sweeps 1
mpiexec -n 3 ./build-mpi/qaplib_example instance.dat --mode mpi-multistart --seed 42 --trials 2 --levels 2 --sweeps 1
mpiexec -n 3 ./build-mpi/qaplib_example instance.dat --mode mpi-replica-exchange --seed 42 --trials 2 --levels 2 --sweeps 1 --exchange-interval 1 --exchange-seed 17
```

These short budgets demonstrate the interface, without a solution-quality
promise. See [mode semantics and reporting](../mpi.md#ordinary-example-command-list)
for global trials, seed offsets, temperatures, sweep budgets and timing scope.
