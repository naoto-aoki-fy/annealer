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
