# G1 instrumentation and multi-start measurements

Measured on 2026-11-06 using a Release build with GNU C++ 11.5.0,
`-Wall -Wextra -Wpedantic -Wconversion -O3 -DNDEBUG`, IPO disabled.
The [metadata CSV](g1-budget.csv.metadata.csv) records the CPU, input dimensions,
configuration, and separately measured preparation/output times. The local G1
input has SHA-256
`73bf704d8ffc55ba42260ab4cb659e3dcb6e729be70404d2cf476ba4e46d1665`.
The dataset and pre-existing executable are not included.

Commands and timing scopes are in the [G1 guide](../examples/gset-g1.md).
All returned best and final assignments were checked using independent edge-cut
calculations. Release validation passed all ten CTests. The final build used the
warning flags above. No warnings came from the changed code; existing
`-Wconversion` warnings remain in unchanged `benchmarks/serial_benchmark.cpp`
and `tests/integration_tests.cpp`.

## Measurement cost

One seed (0), 200 levels, 100 sweeps per level, temperatures 30 to 0.05.
Execution order was detailed profiling, normal, history-only. All three paths
produced identical best/final assignments, energies, and move counts.
They attempted 16,000,000 moves, accepted 4,297,572, and found cut 11624.

| Path | Seconds | Increase over normal |
| --- | ---: | ---: |
| Normal | 21.7324 | — |
| History only | 21.8100 | 0.0776 s (0.36%) |
| Detailed profiling | 24.5741 | 2.8417 s (13.08%) |

These are single measurements, not confidence intervals. Host load and clock
instrumentation can affect the differences. The detailed-path target TTS was
17.1982 s at attempted move 9,594,886. Use history-only results below for the
budget comparison.

| Detailed interval | Seconds | Calls | Fraction of total |
| --- | ---: | ---: | ---: |
| Random variable/value selection | 0.7101 | 16,000,000 | 2.89% |
| Incremental proposal | 15.5991 | 16,000,000 | 63.48% |
| Acceptance RNG and Metropolis | 0.7846 | 16,000,000 | 3.19% |
| State/cache acceptance | 0.2431 | 4,297,572 | 0.99% |
| Rebase | 5.8110 | 16,787 | 23.65% |
| Full best-candidate evaluation/storage | 0.0895 | 558 | 0.36% |
| Initialization/finalization/logging and other residual | 1.3367 | — | 5.44% |

Proposal evaluation dominates this profile, followed by rebasing. These
measurements identify where to investigate; this change introduces no
Max-Cut-specific optimization. The residual includes clock overhead and loop
bookkeeping and is not a pure initialization measurement.

Raw evidence: [profile stdout](g1-profile.txt),
[profile best-update history](g1-profile.csv).

## Equal-budget comparison

Each setting receives 240 seconds, starting at seed zero, using history-only
instrumentation. A started iteration always finishes all its updates. Only
records at or before the setting deadline contribute to budget-limited results.
Completed-iteration success rates exclude the final overrun iteration, even
when its early history contributes a target hit.

Raw evidence: [best-update histories](g1-budget.csv),
[per-iteration results](g1-budget.csv.iterations.csv),
[comparison summary](g1-budget.csv.summary.csv),
[metadata](g1-budget.csv.metadata.csv), and [stdout](g1-budget.txt).

| Levels | Best cut within 240 s | First target (s) | Completed | Successful | Success rate | Overrun (s) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 50 | 11624 | 32.9006 | 41 | 14 | 34.15% | 1.5191 |
| 100 | 11624 | 7.8628 | 21 | 10 | 47.62% | 9.6878 |
| 200 | 11624 | 15.3102 | 11 | 10 | 90.91% | 19.7818 |

The input was loaded and prepared once in 0.0137 s. Buffered result reporting
took 0.1369 s outside the setting budgets. Each setting has one additional
overrun iteration, excluded from completed counts. The CSV aggregate was
independently recomputed from all raw histories and iteration rows, including
deadline filtering and success/completion counts.

For this run, 100 levels reached the target first, while 200 levels had the
largest fraction of successful completed iterations. All settings had the same
best cut within budget. These results alone do not establish a statistically
reliable ranking or a general expected time to solution.
