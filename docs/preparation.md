# Native compilation and runtime preparation

## Native compilation

User functions compile ahead of time as ordinary C/C++ code. The example exposes
function bodies, a `constexpr` domain array, and a fixed interaction coefficient in
[user_model.hpp](../examples/user_model.hpp). `native_factor<function>(variable)`
specializes a small native invoker on that function. Its body and any fixed
template parameters are available for compiler inlining and simplification. The
function need not be constexpr. There is no interaction unrolling, code generator,
parser, JIT, symbolic analysis, or Cartesian interaction table.

The general `Factor(variable, callable)` path supports runtime model information
and fixed captured parameters. Domain vectors, term storage, coefficients copied
into the builder, and incidence data remain runtime objects. Direct fallback uses
the owning type-erased callback; a specialized invoker can inline the known body
but may still have an outer indirect call. This does not promise elimination of
all dispatch from the general path.

Use a Release build. Optional supported link-time optimization:

```sh
cmake -S . -B build-ipo -DCMAKE_BUILD_TYPE=Release -DANNEALER_ENABLE_IPO=ON
cmake --build build-ipo --parallel 1
ctest --test-dir build-ipo --output-on-failure --parallel 1
```

CMake rejects IPO when unsupported. The project never enables fast-math, which
would invalidate the documented finite checks.

## Runtime preparation and memory

`prepare(model, PreparationConfig{})` returns an owning immutable `PreparedModel`.
This is runtime preprocessing, including when the original constants are visible
to the compiler. The source model is retained for fallback and reference checks.
`solve(prepared, config)` reuses preparation; `solve(model, config)` prepares
automatically. `solve_reference(model, config)` retains the original solver.

Preparation flattens terms and factors into contiguous vectors and builds a
contiguous variable-to-term index with offsets. Sorting each term's variable list
deduplicates incidences: repeated variables still have only one incidence per term.
Storage is proportional to variables, terms, factors, unique incidences, and unary
table entries, plus the owned source model. Work is linear in these sizes, plus
per-term variable sorting, callback copying, and tabulation callback costs.

The default `PreparationConfig::lookup_bytes` is 64 MiB of stored double payload.
Zero disables tables. Each factor table is complete or absent. If it does not fit
the remaining budget, direct compiled evaluation is used; later smaller tables
may still fit. `lookup_bytes()` reports payload size. Vector capacity, allocator
overhead, metadata, temporary preparation storage, and source copies are outside
this limit, so it is not a total process memory cap.

Tables use internal domain positions; callbacks receive actual integer values.
Tables are not shared, because arbitrary function identity and captured parameter
equality are not assumed. This also prevents sharing across different domains.
Fully tabulated contribution evaluation reads doubles without callback dispatch.
Partial tables invoke callbacks only for missing factors. No interaction Cartesian
tables are built.

Every tabulated result is checked, including values never visited by a run. A
nonfinite value throws `std::domain_error` during preparation; fallback checks it
when evaluated instead. This earlier error detection is intentional. Leading
empty-factor terms are folded into the offset in original addition order with
checks at each step. Other constants retain their positions. Zero-coefficient
terms with factors are retained: deleting them would hide invalid callbacks or
overflowing products rejected by the reference contract. Leading zero constant
terms are folded using the same ordered rule.

## Incremental state and numerical policy

`PreparedState(prepared, assignment)` owns assignment values, domain positions,
current term contributions, and reusable proposal storage. The prepared model
must outlive the state. `propose(variable, replacement_index)` recomputes products
for incident terms, updating every occurrence of the changed variable, and sums
contribution differences. It uses no division, so zero factors are safe.
`accept()` commits the pending assignment and cache values. Ignoring a proposal
rejects it; the next proposal overwrites scratch without changing current caches.
Invalid or failed proposals invalidate the pending commit.

Ordinary proposal evaluation and acceptance allocate no library storage. Their
cost is proportional to factors in affected terms plus the incident-term count;
direct fallback also pays callback costs. Updating the best assignment copies into
preallocated storage. Preparation, initialization, and result construction may
allocate. Every operation remains serial on the CPU.

`prepared_energy` sums prepared contributions in reference term order. Incremental
deltas sum local differences instead of subtracting two full sums, so rounding
and trajectories can differ near thresholds or nearly tied energies. Factor
values, products, contributions, deltas, and accumulated energies must remain
finite. Arithmetic failures throw `std::overflow_error`. Global full-sum
intermediate overflow is checked at reference checkpoints rather than on every
proposed state. Use `solve_reference` when exact reference summation/error timing
is required.

Every 256 accepted moves, `rebase()` refreshes term caches and anchors accumulated
energy to `full_energy`. Potential strict best improvements are also checked with
the reference evaluator. These validation calls may dispatch callbacks even when
all factors are tabulated, outside ordinary term updates. A near tie within
accumulated rounding error can affect which assignment is saved as best. Both
returned energies are always recomputed from their corresponding assignments by
the reference evaluator. Cross-platform bitwise trajectories are not promised.

Differential tests exhaustively check a small multivalued model and 10,000-step
mixed accept/reject sequences for full, partial, and disabled tables. The tolerance
is `abs(a-b) <= 1e-10 * (1 + max(abs(a), abs(b)))`. Integer-valued solver trajectories
are compared exactly. These tests cover their exercised numerical scales and do
not establish a universal drift bound for ill-conditioned models. The reference
path remains the correctness oracle.
