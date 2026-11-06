# MPI multi-start solver

MPI support is optional and disabled by default. Serial headers and the
`annealer::annealer` target do not depend on MPI.

```sh
cmake -S . -B build-mpi -DANNEALER_ENABLE_MPI=ON
cmake --build build-mpi --parallel 2
ctest --test-dir build-mpi -R mpi_tests --output-on-failure
```

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
