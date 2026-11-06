# Serial Metropolis solver

Include `<annealer/solver.hpp>` and call `annealer::solve(model, config)`:

```cpp
annealer::SolverConfig config;
config.seed = 42;
config.initial_assignment = annealer::Assignment{-2, 2, 3}; // optional
config.initial_temperature = 10.0;
config.final_temperature = 0.1;
config.temperature_levels = 100;
config.sweeps_per_level = 10;
const auto result = annealer::solve(model, config);
```

These are also the defaults, except the default seed is zero and the initial
assignment is absent. Without an assignment, each variable is independently
initialized uniformly over its domain. Supplied assignments are validated by
`full_energy`, including when the move budget is zero. Model construction,
initialization, and execution run serially on the CPU. The solver takes a const
model and keeps all mutable state locally.

## Proposals and acceptance

Each attempt selects uniformly among variables with at least two values, with
replacement. It then selects uniformly among that variable's other values.
Singletons are skipped. With `m` mutable variables and selected domain size `q`,
the probability of either direction of a particular move is `1 / (m * (q - 1))`.
The set of mutable variables and their cardinalities do not depend on the state,
so heterogeneous and nonconsecutive domains preserve symmetry. The implementation
samples an alternative index and passes the corresponding actual value to model
callbacks. One sweep is exactly `m` attempts, rather than a visit to each variable.

The default solver uses [prepared incident-term updates](preparation.md).
`solve_reference(model, config)` evaluates every proposed assignment with
the reference `full_energy` function.
For `delta = proposed_energy - current_energy`, downhill and equal moves always
accept. Uphill moves accept with probability `exp(-delta / T)`. The testable
`metropolis_accept(delta, T, u)` helper requires a finite delta, positive finite T,
and `u` in `[0,1)`. It uses `log(u) < -delta/T`; zero accepts and threshold equality
rejects. This avoids forming an exponential that can underflow. The solver uses
the top 53 random bits of `std::mt19937_64` to obtain an exactly representable
uniform double in `[0,1)`, and draws an acceptance variate only for uphill moves.
Discrete indices use `std::uniform_int_distribution<std::size_t>`.

The evaluator's numerical errors and callback exceptions propagate. A nonfinite
energy difference raises `std::overflow_error`, even if both energies are finite;
no partial result is returned. Overflow in the positive `delta/T` ratio means
rejection for positive `u`. The callback purity and floating-point assumptions in
the [model documentation](model.md) still apply.

## Cooling and move budgets

Both endpoint temperatures must be positive and finite, even for zero moves.
The final temperature must not exceed the initial temperature. Violations throw
`std::invalid_argument`. Equal endpoints give a constant schedule.

For `L >= 2`, level `j` uses geometric cooling:

`T(j) = exp((1 - j/(L-1)) * log(T_initial) + j/(L-1) * log(T_final))`.

The endpoints are returned exactly; intermediate results are clamped to the
endpoint range to contain rounding error. Log interpolation avoids forming an
extreme endpoint ratio. `geometric_temperature(initial, final, L, j)` exposes
this schedule and throws `std::out_of_range` if `j >= L`.

- One level uses the initial temperature only. The final setting is still validated.
- Zero levels or zero sweeps return the valid initialized assignment without moves.
- No mutable variables also returns immediately without moves.
- Otherwise attempts equal `L * sweeps_per_level * m`. A count exceeding
  `std::size_t` raises `std::overflow_error` before attempting moves.

There is no schedule allocation or iteration when no moves are requested.

## Results and reproducibility

`SolverResult` contains `best_assignment`, `best_energy`, `final_assignment`,
`final_energy`, `attempted_moves`, `accepted_moves`, and `seed`. The initial state
participates in best tracking. Strict improvement replaces the best; ties retain
the first observed best. The current state can accept an uphill move without
changing the saved best. Reported energies come from full evaluations of their
corresponding assignments. See the [drift policy](preparation.md) for floating-point
best-state comparisons in the optimized path.

The same model, configuration, seed, and build/runtime environment reproduce a
run. The library does not promise bitwise reproducibility across compilers,
standard libraries, or floating-point environments. Optimization is stochastic;
finding a global optimum is not a correctness guarantee.

The [runnable example](../examples/main.cpp) prints its best assignment, best and
final energy, move counts, and seed after evaluating the original assignment.

## Instrumented prepared solves

`solve_instrumented(prepared, config, InstrumentationConfig{.profile = false})`
returns an `InstrumentedResult`: the ordinary `result`, `history`, total
`seconds`, and six `profile` intervals. The prepared model owns its source model and can be reused across solves.

History begins with the validated initial energy at attempted/accepted count
zero, then records only strict improvements confirmed by full energy evaluation.
Each `BestUpdate` contains `seconds`, `attempted_moves`, `accepted_moves`, and
`best_energy`. Time starts before initialization, excluding model construction
and preparation. Equal-energy accepted moves do not create records. Zero levels,
zero sweeps, and immutable models still produce an initial record. Instrumentation
never stops a solve on target achievement.

The normal, history, and profile paths instantiate the same internal loop.
History and timing operations are removed with `if constexpr` in normal solves.
Random draws, acceptance, best-candidate validation, and the rebase cadence of
256 accepted moves are shared.

With `.profile = true`, `profile_names` identifies six disjoint intervals:
proposal selection, incremental proposal evaluation, acceptance RNG and criterion,
state acceptance, rebase, and full evaluation/storage of best candidates.
Each interval exposes seconds and calls. Best-candidate calls include candidates
that fail exact validation. Initial and final full evaluations are outside that
interval. Total time minus the interval sum contains initialization, finalization,
history logging, scheduling, counters, and clock/instrumentation overhead; it
cannot isolate each of these costs. Per-move clock calls perturb performance.
Use history-only solves for time-budget comparisons and compare all three paths
with the same seed to assess measurement cost. No timing threshold defines test
success. See the [G1 guide](examples/gset-g1.md) for executable commands and CSVs.

## Optional swap proposals

Set `config.swap_proposals = true` to exchange the values at two distinct
variable positions. The default is false and retains the existing replacement
proposal and random draws. Swap mode is supported by `solve_reference`, prepared
`solve`, and `solve_instrumented`.

All variables must have identical domains, including the same stored value
order. This is checked even when the move budget is zero. Models with fewer
than two variables perform no swaps. Otherwise one sweep is exactly `n`
attempts for `n` variables, including singleton domains and equal-valued swaps.
Each attempt draws a uniform first position and a uniform distinct second
position. Each ordered pair has probability `1/(n*(n-1))`; a particular
unordered pair has probability `2/(n*(n-1))`, giving a symmetric proposal.
Equal-valued swaps are allowed and count as accepted zero-energy moves.

Supply a permutation through `initial_assignment` for permutation problems.
Every accepted swap preserves its multiset, including repeated values if present.
Default independent initialization does not generate a permutation. Swaps cannot
change the initial multiset or enforce other constraints; they only explore
reorderings. This mode provides a useful neighborhood for permutation models
such as QAP and TSP without implementing those models.

Acceptance uses the complete simultaneous swap energy difference. Prepared
updates evaluate the union of incident terms, including terms involving both
positions only once. Finite arithmetic checks, the accepted-move rebase cadence,
and reference checks for best candidates apply as in replacement mode. Swap
replay is deterministic for a fixed seed and configuration under the same
standard-library implementation. It uses its own sequence of proposals.
