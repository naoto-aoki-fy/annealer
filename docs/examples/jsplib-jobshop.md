# JSPLib job-shop parsing and schedule validation

This first job-shop component provides an instance parser, an independent
schedule evaluator and a deterministic feasible initialization. It does not
construct an energy model or run an optimizer.

A job is an ordered chain of nonpreemptive operations. Each operation requires
one machine exclusively for its positive integer processing duration. Operations
from different jobs may run concurrently on different machines. The objective
is to minimize makespan, the latest operation completion time.

## Supported input and representation

[JSPLib's standard text format](https://scheduleopt.github.io/benchmarks/jsplib#standard-format)
contains positive job and machine counts `J M`, followed by `M` ordered
`machine duration` pairs for each of the `J` jobs. Machine IDs are **zero-based**,
in `[0,M)`. All fields are decimal integers, separated by arbitrary whitespace.
For example:

```text
2 2
0 3 1 2
1 4 0 1
```

The first job runs on machine 0 for three time units, then on machine 1 for two.
The second job runs on machine 1 for four, then machine 0 for one.

`examples/jsplib_jobshop.hpp` defines `jobshop::Instance` and `Operation`.
Every operation explicitly stores zero-based `job`, `step`, `machine`, and
`duration`. The flat vector is in job-major order; operation index is
`job*M+step`. `read(std::istream&)` validates dimensions and allocation bounds,
positive durations, machine IDs, completeness, trailing data, and integer ranges.
It grows the operation vector as records arrive. Repeated machine visits within
a job are accepted; the format still requires exactly `M` operations per job.

`Time` is signed 64-bit. Total duration must fit in that type, providing a safe
serial horizon. This conservative requirement may reject an instance whose
parallel makespan alone would fit. Overflowing sums throw `std::overflow_error`;
invalid fields or structure throw `std::invalid_argument`. The parser accepts
neither comments, JSON, nor DaColTeppan variable-length `-1 -1` job records.
Use the source repository's standard `.txt` files for ft06 and la01.

## Schedule verification

A `jobshop::Schedule` contains one nonnegative integer start time per operation,
in the same job-major order. `evaluate(instance, starts)` returns a structured
`Verification` with:

- `precedence_violations`: adjacent job edges where the predecessor completes
  after the successor starts.
- `machine_conflicts`: unordered operation pairs on the same machine whose
  half-open intervals `[start,start+duration)` overlap.
- `makespan`: maximum completion time, computed even for infeasible schedules.
- `feasible()`: true exactly when both violation counts are zero.

Touching endpoints are allowed. Every conflicting pair is counted, including
same-job pairs and nested intervals; three mutually overlapping operations give
three conflicts. A pair may contribute to both precedence and machine counts.
Wrong schedule sizes, negative starts, malformed instance metadata and completion
overflow throw exceptions rather than producing a misleading verification result.
An infeasible schedule's makespan is not a valid scheduling upper bound.

For the example above, starts `{0,4,0,4}` are feasible with makespan 6.
Starts `{0,3,0,4}` have one machine conflict and makespan 5.
`serial_schedule(instance)` runs all operations consecutively in job-major order,
returning `{0,3,5,9}` with makespan 10. It always respects job order and avoids
all machine overlap on supported instances. This initialization is reproducible
and feasible, but intentionally makes no effort to optimize parallelism.

Validation takes O(K^2) time for K operations to count every machine-conflict
pair, and O(K) temporary storage for completion times. The serial constructor is
O(K). Both validate the instance structure and duration horizon first.

## Benchmarks and acquisition

[JSPLib](https://scheduleopt.github.io/benchmarks/jsplib) includes the classic
Fisher–Thompson `ft06` (6 jobs, 6 machines, 36 operations) and Lawrence `la01`
(10 jobs, 5 machines, 50 operations). Download standard instances locally:

```sh
curl -fL https://raw.githubusercontent.com/ScheduleOpt/benchmarks/main/jobshop/instances/text/FisherThompson1963/ft06.txt -o /tmp/ft06.txt
curl -fL https://raw.githubusercontent.com/ScheduleOpt/benchmarks/main/jobshop/instances/text/Lawrence1984/la01.txt -o /tmp/la01.txt
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 1
ctest --test-dir build --output-on-failure -R jsplib_jobshop
```

Record the source revision and file checksums when reproducing results from the
[ScheduleOpt repository](https://github.com/ScheduleOpt/benchmarks). Public data
is not bundled and is not required or downloaded by tests. Regression tests use
synthetic ft06/la01 dimensions, hand-calculated schedules and malformed input.

The helper is header-only C++20 and depends only on the standard library:

```cpp
#include "examples/jsplib_jobshop.hpp"
#include <fstream>

std::ifstream input("/tmp/ft06.txt");
const auto instance = jobshop::read(input);
const auto starts = jobshop::serial_schedule(instance);
const auto result = jobshop::evaluate(instance, starts);
// result.feasible() is true; result.makespan is a valid serial upper bound.
```

A subsequent change can build a factorized energy model and reuse this independent
validator; this component has no annealer or solver dependency.
