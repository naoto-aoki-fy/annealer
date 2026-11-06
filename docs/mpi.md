# MPI multi-start and replica exchange

MPI support is optional and disabled by default. Serial headers and the
`annealer::annealer` target do not depend on MPI.

## Build, run, and test

From the repository root, with CMake 3.20+, C/C++20 compilers, and (for
MPI-enabled builds) an MPI implementation providing C++ discovery and `mpiexec`:

```sh
# Serial build: no MPI installation required.
cmake -S . -B build-serial -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DANNEALER_ENABLE_MPI=OFF
cmake --build build-serial --parallel 2
ctest --test-dir build-serial --output-on-failure
./build-serial/multivalued_example

# MPI build: includes the serial targets and the self-contained MPI example.
cmake -S . -B build-mpi -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DANNEALER_ENABLE_MPI=ON
cmake --build build-mpi --parallel 2
mpiexec -n 1 ./build-mpi/mpi_example
mpiexec -n 3 ./build-mpi/mpi_example
ctest --test-dir build-mpi --output-on-failure
# Run only communication integration tests:
ctest --test-dir build-mpi -R '^(mpi_tests|replica_exchange_tests)_' --output-on-failure
```

CMake launches tests through `MPIEXEC_EXECUTABLE`, `MPIEXEC_NUMPROC_FLAG`,
`MPIEXEC_PREFLAGS`, and `MPIEXEC_POSTFLAGS`; set these at configure time for
your launcher or scheduler. Both integration suites run with 1, 2, 3, and 4
processes, finite 60-second timeouts, and matching CTest `PROCESSORS` counts.
Reserve at least four process slots or configure launcher oversubscription as
appropriate for your MPI installation. The suites check serial equivalence,
odd process counts, uneven start distribution, idle ranks, subcommunicators,
and coordinated rank-local validation/callback failures with recovery.

`examples/mpi.cpp` calls both APIs and reuses `demos::nonlinear()` from
`examples/demo_models.hpp`. It has no external inputs. Every rank constructs
and prepares the same model; preparation failures are agreed before entering
the solver. Normal CSV output comes only from communicator rank 0. No files
are written by ranks. For applications that write local states, include the
communicator rank (and a communicator/group identifier when needed) in filenames.

Enabling the option requires the MPI C++ component and creates the interface
target `annealer::mpi`, linked to `MPI::MPI_CXX`. Include `<annealer/mpi.hpp>`.
The caller initializes MPI before calling and finalizes it afterward:

```cpp
MPI_Init(&argc, &argv);
// Construct an equivalent Model or PreparedModel locally on every rank.
annealer::SolverConfig config;
config.seed = 42; // base seed
const auto result = annealer::solve_multistart(comm, prepared, 10, config);
// result.best_assignment, best_energy, start_id, seed exist on every rank.
MPI_Finalize();
```

`comm` is a caller-owned valid intracommunicator, including a subcommunicator.
Every member must call the same overload in the same collective order. There
are overloads for `const Model&` and `const PreparedModel&`. Model preparation
is local; prepared inputs are reused across starts. Model objects and callbacks
are never serialized. Model construction outside this API remains the caller's
responsibility, including coordinating construction failures before entry.

`num_starts` is a positive `uint64_t` total across the communicator. Logical
start IDs are `[0, num_starts)`. Rank `r` executes IDs `r, r + size, ...` exactly
once, using serial `solve` with seed `config.seed + start_id`. Seed overflow is
rejected before any start. Idle ranks still participate in all collectives.
Initial assignments, replacement proposals, swap proposals, and zero budgets
have the same semantics as the [serial solver](solver.md).

The winner is selected from every start's best-so-far state, ordered by lowest
energy then lowest start ID. Its assignment is broadcast and its full energy
verified on every rank. Changing process count preserves individual runs and
the winner within the serial solver's reproducibility guarantees.

`MultiStartResult` contains the global `best_assignment`, `best_energy`,
`start_id`, and `seed`. `attempted_moves` and `accepted_moves` are checked
`uint64_t` sums over all starts, including losing starts. Counter overflow is
an error. `local_runs` contains only this rank's starts, ordered by ID, each
with `start_id` and the complete serial `SolverResult` in `result`. Its final
states and counters are local; there is no global final state. Idle ranks have
an empty vector. Retaining these results costs memory proportional to local
start count times assignment size.

## Collective requirements and errors

All ranks must provide identical start count, base seed, schedule, proposal
mode, and optional initial assignment (including presence and value order).
These fields and model variable counts are checked collectively. Local models
must be equivalent, with identical domain order, factors, callback behavior,
and preparation semantics; arbitrary callbacks cannot be compared. Winner
energy verification detects some, but not all, model inconsistencies. Callbacks
must obey the serial purity requirements and must not invoke MPI collectives.
Assignments larger than `INT_MAX` are rejected.

Validation, allocation, preparation, solver, callback, and verification exceptions
inside the API are caught at phase boundaries. If any rank fails, all ranks
throw `std::runtime_error("MPI multi-start failed on one or more ranks")` after
an error agreement collective; no partial result is returned. Original exception
messages/types are intentionally not transported. A failing rank may wait for
other ranks to finish their local runs before the error is reported. A subsequent
collective call is allowed once all callers have handled the exception.

MPI transport errors use the communicator's MPI error handler. Process loss,
invalid communicators, mismatched collective calls, exceptions before entry,
nonterminating callbacks, and MPI failures are outside this exception policy;
there is no fault-tolerant MPI recovery. The caller owns communicator lifetime
and MPI thread-level compliance.

## Temperature replica exchange

Include `<annealer/replica_exchange.hpp>` and link `annealer::mpi`:

```cpp
annealer::ReplicaExchangeConfig config;
config.temperatures = annealer::geometric_ladder(0.1, 30.0, communicator_size);
config.solver.seed = 42;
config.total_sweeps = 1000;
config.exchange_interval = 10;
config.exchange_seed = 900;
const auto result = annealer::solve_replica_exchange(comm, prepared, config);
```

The Model and PreparedModel overloads use the same MPI ownership, equivalent-model
requirements, and collective exception policy as multi-start. Each rank keeps
one fixed temperature; the explicit ladder must have exactly communicator size
entries in nondecreasing order, all positive and finite. Equal temperatures are
allowed. `geometric_ladder(cold, hot, count)` includes both endpoints for two or
more replicas; one replica uses `cold`. A zero count is invalid.

All configuration fields must agree across ranks. `solver.seed + rank` seeds the
local chain, with checked overflow. The optional initial assignment and proposal
mode are inherited from `solver`. The solver schedule fields are checked for
consistency but replaced by the constant rank temperature and `total_sweeps`.
`exchange_interval` must be positive, even when exchanges are disabled. Total
sweeps may be zero: initialization and global best selection still happen, with
no moves or exchanges. A single rank performs constant-temperature local work
and has no exchange partner. `exchanges = false` produces the same local result
as serial solving for one temperature level with `total_sweeps` sweeps, rank seed,
and the same initialization/proposal settings.

After each interval, including a final partial interval, round zero pairs
(0,1), (2,3), ...; round one pairs (1,2), (3,4), ...; these patterns alternate.
Unpaired ranks still join error agreement collectives. Current energies are fully
evaluated before communication. For a pair i,j, acceptance is
`min(1, exp((1/T_i - 1/T_j) * (E_i - E_j)))`. The lower rank makes exactly one
decision and sends it to its partner. It draws a uniform from the top 53 bits of
a separate `mt19937_64(exchange_seed + rank)` stream, once per owned pair attempt.
Local proposal RNGs are never used for exchanges. Seed overflow is rejected.
Equal temperatures or energies accept; zero uniform accepts. Other decisions
compare logarithms, using scaled differences to avoid overflow for extreme finite
temperatures and energies. `replica_exchange_accept` exposes this numerical test;
invalid inputs throw `invalid_argument`.

Communication uses matched `MPI_Sendrecv` calls for energies, decisions, and
accepted `MPI_INT64_T` assignment payloads. No C++ model or cache objects are
transmitted. Received assignments rebuild the local caches through the internal
execution interface. Temperatures, RNG streams, move counters, rebase cadence,
and historical bests stay on their ranks. This temperature exchange is distinct
from `solver.swap_proposals`, which swaps two values within a local assignment.
The two features can be used together.

`ReplicaExchangeResult.local` is the rank's complete `SolverResult`.
`best_assignment`, `best_energy`, and `best_rank` describe the best historical
state observed anywhere, available and energy-verified everywhere. Equal energies
choose the lowest rank; a rank retains its first strict best. Imported states
participate in this best tracking. `exchange_attempts` and `exchange_accepts` are
global counts of pairs, each counted once. `local_exchange_attempts` and
`local_exchange_accepts` count participation on this rank, so their communicator
sums are twice the global counts. All exchange counts are separate from local
Metropolis move counts; global exchange counter overflow is a collective error.

All members must call the same API in the same order. Use a communicator dedicated
to this operation, without concurrent traffic using its internal tags 0, 1, or 2.
Preparation, initialization, local updates, decision-making, assignment installation,
and final verification have collective error boundaries before peers proceed to
communication. There is no partial result or process-failure recovery. This first
implementation has one replica per rank, a fixed ladder, and blocking exchanges;
it does not adapt temperatures or move replicas between rank temperatures.

## Reproducible small benchmark

The MPI example is also a compact benchmark. Run the same executable for each
process count; preserve separate output files for each run:

```sh
for repeat in 1 2 3; do
    for ranks in 1 2 3 4; do
        mpiexec -n "$ranks" ./build-mpi/mpi_example > "mpi-${ranks}-${repeat}.csv"
    done
done
```

Record your CPU, compiler/version, MPI implementation/version, build type, process
placement, and machine load alongside the CSVs. Use the same build and placement
policy for comparisons. These are measurements, with no runtime or quality
thresholds. This four-variable demonstration (three mutable variables, one
singleton) is deliberately small: communication overhead can dominate, and it
is not evidence of application-scale speedup.

Each CSV reports process count, base and exchange seeds, starts or replicas,
schedule parameters, proposal mode, total attempted local moves, maximum elapsed
solve time across ranks, best energy, and global exchange pair attempts/accepts.
Timing uses `MPI_Wtime` around each prepared-model solver call after a barrier,
then `MPI_MAX` reduction. Model construction and preparation, the initial barrier,
reporting, and reporting reductions are excluded; solver initialization, internal
communication, and final verification are included. No warmup is performed.

Multi-start always executes exactly the same seven start IDs `[0,7)` with seeds
42 through 48, independent of process count. Each start uses 16 geometric
annealing levels from 10 to 0.1 with 2,000 sweeps per level and replacement
proposals: 96,000 attempted moves per start, 672,000 globally. Compare its
maximum solve time to the one-process multi-start run for a fixed-work scaling
measurement. Process counts above seven leave some ranks idle.

Replica exchange uses one replica per process, a fixed geometric ladder from
0.1 to 10 (one process uses only 0.1), 32,000 sweeps per replica, and exchange
interval 100. Chain seeds are `42 + rank`; decision seeds are `900 + rank`.
Total attempted local work is `96,000 * processes`, reported by summing the
returned local counters. There are 320 exchange rounds with alternating adjacent
pairs; exchange counters count pairs once, separately from local work. Increasing
process count changes both the ladder and total work, so this is a quality/work
comparison, not fixed-work strong scaling. The two modes also have different
budgets and temperature policies: their raw runtimes do not establish an
algorithmic speedup. To compare fixed local work for replica exchange, edit
`total_sweeps` per process count and report the resulting actual counters and
ladder. Repeat runs and examine energy as well as runtime; neither better energy
nor speedup is guaranteed.

### Example measurement (2026-11-06)

Local x86-64 Intel Core i7-6700K, GCC 11.5.0, Release (`-O3 -DNDEBUG`),
MPICH 5.0.2, default launcher placement, no explicit affinity. One measured run
per process count, no warmup, with a concurrent build and uncontrolled background
load; these numbers illustrate the output and should not establish scaling claims.
All configurations and seeds are those above. Preparation is excluded.

| Processes | Multi-start moves | Multi-start seconds (max) | Replica moves | Replica seconds (max) | Exchange accepts/attempts |
| --- | --- | --- | --- | --- | --- |
| 1 | 672000 | 0.043580 | 96000 | 0.006308 | 0/0 |
| 2 | 672000 | 0.025133 | 192000 | 0.008050 | 160/160 |
| 3 | 672000 | 0.021336 | 288000 | 0.014067 | 320/320 |
| 4 | 672000 | 0.023691 | 384000 | 0.015851 | 480/480 |

Every row in both modes obtained best energy -1710. All exchanges accepted on
this small model; the integration tests separately exercise rejection behavior.
The MPI-disabled Release build passed all 24 serial tests. Both MPI integration
suites passed at 1, 2, 3, and 4 ranks (8/8 tests), including the existing
subcommunicator, idle-rank, uneven-distribution, and controlled failure checks.

## Basic example execution

`minimal_example`, `binary_hubo_example`, `nonlinear_example`, and
`multivalued_example` accept `--mode serial|mpi-multistart|mpi-replica-exchange`.
The default is serial, including in MPI-enabled builds without `mpiexec`.
MPI modes require `-DANNEALER_ENABLE_MPI=ON`. For example,

```sh
mpiexec -n 4 ./build-mpi/binary_hubo_example --mode mpi-multistart --trials 5 --seed 0
mpiexec -n 4 ./build-mpi/minimal_example --mode mpi-replica-exchange --trials 2 \
  --initial-temperature 10 --final-temperature 0.1 --levels 100 --sweeps 10 \
  --exchange-interval 10 --exchange-seed 0
```

All four accept `--seed`, `--trials`, `--initial-temperature`,
`--final-temperature`, `--levels`, `--sweeps`, and `--help`. No arguments retain
the original model, initial assignment, solver defaults and output. Binary HUBO
and nonlinear default to five trials with seeds 0 through 4. Minimal and
multivalued default to one trial with seed 42. Temperature and sweep defaults
remain those defined by each example and `SolverConfig`. Serial trials use
consecutive seeds from the requested seed.

In multi-start mode, trials are the total logical starts across all ranks,
including when there are fewer starts than ranks. The MPI API receives the base
seed unchanged and adds the global start ID. Output gives the global best
assignment, winning start and seed, and global move counters.

In replica-exchange mode, trials are independent collective runs with one replica
per rank. The geometric ladder runs from the final (cold) temperature to the
initial (hot) temperature. A single rank uses the cold endpoint. Levels times
sweeps determines the total per-replica sweep budget, rather than a cooling
schedule. `--exchange-interval` defaults to 10 and `--exchange-seed` defaults to
0. These two options are accepted only in replica-exchange mode. Collective run
`t` adds `t * communicator_size` to both base seeds, after which the MPI API adds
the rank. Seed offsets and sweep-budget products are checked for overflow.
Output reports the global best assignment and global exchange pair counters for
each run, plus explicitly labeled rank-zero local move counters.

Only rank zero prints normal output. Argument handling, construction,
preparation, serial solving and result processing use collective failure
boundaries in MPI builds. MPI solvers retain their own collective error handling.
MPI initialization and finalization are owned by the example when necessary.

### File-based examples

`gset_g1_example`, `dimacs_coloring_example`, `maxsat_example`,
`set_covering_example`, `qaplib_example`, `tsplib_example`, and `jsplib_jobshop_example` accept the same `--mode serial|mpi-multistart|mpi-replica-exchange`,
`--exchange-interval` (default 10), and `--exchange-seed` (default 0) options as
the basic examples. Their existing file argument, defaults, preparation settings,
and problem-specific validation remain in effect. For example:

```sh
mpiexec -n 2 build-mpi/maxsat_example instance.wcnf --mode mpi-multistart --trials 4
mpiexec -n 2 build-mpi/dimacs_coloring_example graph.col --colors 3 \
  --mode mpi-replica-exchange --trials 2 --exchange-interval 10
```

All ranks must receive identical arguments and construct the same instance and
model. Each rank reads the input: the file must be accessible and identical on
all ranks. Input, parsing, construction, preparation, and independent result
verification errors are coordinated before further collectives.

Multi-start trials count global starts, with seed `base + start_id`.
Replica-exchange trials count collective runs, with rank seed
`base + trial * ranks + rank`; the exchange seed advances by `trial * ranks`.
The cold-to-hot geometric ladder spans the final and initial temperatures;
`levels * sweeps` is each replica's total sweep budget.

Rank 0 emits one report. MPI records are numbered `Global record` (one winner
for multi-start, one winner per replica-exchange run). Every local solver result
is independently checked, as is each global winner. Set covering's separate
best-feasible summary includes local candidates from every rank. An infeasible
candidate is a successful search execution, with infeasibility in its status.
Solve timings include result verification and, in MPI modes, collective work.

Job shop retains its initial schedule and model-memory check. `--memory-mib`
applies separately to each process's model construction in MPI mode; it is not
a total memory budget for the MPI job.

## Ordinary example command list

Use the build commands above. These commands use MPI-enabled binaries, including
serial execution without a launcher. For an MPI-disabled build use the same serial
commands with `build-serial` instead of `build-mpi`.

```sh
./build-mpi/minimal_example --mode serial --seed 42 --trials 2 --levels 2 --sweeps 1
mpiexec -n 3 ./build-mpi/minimal_example --mode mpi-multistart --seed 42 --trials 2 --levels 2 --sweeps 1
mpiexec -n 3 ./build-mpi/minimal_example --mode mpi-replica-exchange --seed 42 --trials 2 --levels 2 --sweeps 1 --exchange-interval 1 --exchange-seed 17
./build-mpi/binary_hubo_example --mode serial --seed 42 --trials 2 --levels 2 --sweeps 1
mpiexec -n 3 ./build-mpi/binary_hubo_example --mode mpi-multistart --seed 42 --trials 2 --levels 2 --sweeps 1
mpiexec -n 3 ./build-mpi/binary_hubo_example --mode mpi-replica-exchange --seed 42 --trials 2 --levels 2 --sweeps 1 --exchange-interval 1 --exchange-seed 17
./build-mpi/nonlinear_example --mode serial --seed 42 --trials 2 --levels 2 --sweeps 1
mpiexec -n 3 ./build-mpi/nonlinear_example --mode mpi-multistart --seed 42 --trials 2 --levels 2 --sweeps 1
mpiexec -n 3 ./build-mpi/nonlinear_example --mode mpi-replica-exchange --seed 42 --trials 2 --levels 2 --sweeps 1 --exchange-interval 1 --exchange-seed 17
./build-mpi/multivalued_example --mode serial --seed 42 --trials 2 --levels 2 --sweeps 1
mpiexec -n 3 ./build-mpi/multivalued_example --mode mpi-multistart --seed 42 --trials 2 --levels 2 --sweeps 1
mpiexec -n 3 ./build-mpi/multivalued_example --mode mpi-replica-exchange --seed 42 --trials 2 --levels 2 --sweeps 1 --exchange-interval 1 --exchange-seed 17
```

File-example commands: [G-set](examples/gset-g1.md#common-execution-commands),
[coloring](examples/dimacs-coloring.md#common-execution-commands),
[MaxSAT](examples/maxsat.md#common-execution-commands),
[set covering](examples/set-covering.md#common-execution-commands),
[QAP](examples/qaplib.md#common-execution-commands),
[TSP](examples/tsplib.md#common-execution-commands), and
[job shop](examples/jsplib-jobshop.md#common-execution-commands).

Serial trials are sequential starts with seeds `base + trial`. Multi-start trials
are global starts, distributed across ranks, with seeds `base + start_id`; ranks
may have no starts. Both use a geometric cooling schedule with `levels`
temperatures and `sweeps` sweeps at each level. Replica-exchange trials are
collective runs, each with one replica per rank, `levels * sweeps` sweeps per
replica, and a fixed cold-to-hot ladder from final to initial temperature. One
rank uses only the cold endpoint. Changing rank count changes this ladder and
the search. Modes do not imply equal computational budgets or identical results.
Replica run `t` uses chain seed `base + t * ranks + rank` and exchange decision
seed `exchange_base + t * ranks + rank`. All seed arithmetic and the sweep
product are checked. Exchange interval must be positive and exchange options
are valid only for replica exchange.

QAP/TSP use swap proposals: MPI multi-start uses one common permutation generated
from the base seed for all starts; replica run `t` uses a common permutation from
`base + t * ranks`. Serial starts shuffle independently from their trial seed.
Thus even multi-start need not reproduce serial QAP/TSP trajectories.

Normal reports come only from rank zero. A report may contain multiple trial or
collective-run records; these are not duplicate rank reports. Returned candidates
are independently evaluated, including permutation/constraint checks. Search
infeasibility is reported and is not necessarily an execution error. Basic
examples report global winners and global multi-start move counts; replica move
counts are explicitly local to rank zero. File-example solve timings include
verification and MPI collective work, measured on rank zero, not a maximum across
ranks. Total timings include model preparation and reporting as labeled, excluding
input parsing. G-set labels its MPI solving/verification wall time separately.
The dedicated `mpi_example` instead reports the maximum solver time across ranks,
as detailed above. Timing fields are excluded from deterministic replay checks.

The CLI audit uses tiny generated local fixtures, one and three ranks in addition
to the existing two-rank checks, idle multi-start ranks, malformed/missing inputs,
checked arithmetic, deterministic replay and root-only summaries. Existing MPI
algorithm suites retain responsibility for communication and algorithm invariants;
no stochastic optimum, speedup or improvement is required.
