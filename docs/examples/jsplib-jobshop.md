# JSPLib job-shop parsing and schedule validation

The job-shop example provides an instance parser, an independent schedule
evaluator, deterministic feasible initialization and a factorized annealing model.

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

Set `ANNEALER_WORK_DIR` to an absolute directory of your choice before running
these commands. It defaults to `TMPDIR` when set, otherwise `/tmp`. Run the
following blocks in the same shell.

```sh
export ANNEALER_WORK_DIR="${ANNEALER_WORK_DIR:-${TMPDIR:-/tmp}}"
mkdir -p "${ANNEALER_WORK_DIR}"
curl -fL https://raw.githubusercontent.com/ScheduleOpt/benchmarks/main/jobshop/instances/text/FisherThompson1963/ft06.txt -o "${ANNEALER_WORK_DIR}/ft06.txt"
curl -fL https://raw.githubusercontent.com/ScheduleOpt/benchmarks/main/jobshop/instances/text/Lawrence1984/la01.txt -o "${ANNEALER_WORK_DIR}/la01.txt"
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
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <string>

const char* work_dir = std::getenv("ANNEALER_WORK_DIR");
if (work_dir == nullptr || *work_dir == '\0') {
    throw std::runtime_error("Set ANNEALER_WORK_DIR before running this example");
}
std::ifstream input(std::string(work_dir) + "/ft06.txt");
const auto instance = jobshop::read(input);
const auto starts = jobshop::serial_schedule(instance);
const auto result = jobshop::evaluate(instance, starts);
// result.feasible() is true; result.makespan is a valid serial upper bound.
```

The parser/validator header remains independent of the annealer. The model and
runner below use it to verify solver results.


## Factorized finite-horizon model

`examples/jobshop_model.hpp` uses operation starts `s_i` in `[0,H-duration_i]`
and a declared makespan `M` in `[0,H]`, with `H=sum(duration)`. The serial
initial schedule fits this horizon, and some optimum fits it as well. The energy is

`E = M + P*(precedence violations + machine conflicts + makespan-bound violations)`.

Each violated adjacent precedence, overlapping same-machine pair, or completion
beyond M counts once (not once per overlapping time unit). With `P=H+1`, every
infeasible energy exceeds H, while the serial feasible initialization has energy
H. Thus global minima are feasible and minimize makespan. All arithmetic and the
bound `H+P*constraint_count <= 2^53` are checked for exact integer energies in
double; this is more restrictive than the standalone parser's time range.

For any binary constraint predicate V, the exact factorization is
`V(s_i,y) = sum_a [s_i=a]*V(a,y)`. In each term `a` is an immutable captured
constant, so `V(a,y)` is a pure unary function. Precedence uses
`[a+duration_i > s_next]`; makespan uses `[a+duration_i > M]`; overlap uses
`[a < s_j+duration_j && s_j < a+duration_i]`. No model factor takes two variable
values, and no occupancy array or auxiliary Boolean variables are introduced.
Although often called a HUBO formulation, these finite-domain terms have at most
two factors, using the same public factorized higher-order API.

For K operations, R adjacent precedence edges and C same-machine pairs, term
count is bounded by `1+(K+R+C)*(H+1)`. Start domains are shortened by duration.
Domain storage is O(KH), model/cache storage O((K+R+C)H), and full evaluation
has the same term complexity. An incremental proposal touches all incident
terms; the makespan variable touches all makespan-bound expansions.

`estimate(instance)` computes exact term counts before model allocation. The
planning allowance is 2048 bytes per term plus 32 bytes per domain entry, covering
source and prepared copies, callback captures, metadata, caches and scratch.
It is a conservative estimate, not a portable allocator guarantee or a process
RSS limit. Preparation uses an additional lookup budget of 1 MiB. The CLI prints
the estimate before construction and enforces `--memory-mib` (default 1024).
Allocation can still fail on a constrained host.

For the original files checked during development:

| Instance | H | Variables | Terms | Planning bytes |
| --- | ---: | ---: | ---: | ---: |
| ft06 | 197 | 37 | 30,008 | 61,684,512 |
| la01 | 2849 | 51 | 880,759 | 1,808,354,464 |

la01 exceeds the default limit: constructing it requires an explicit larger
budget and enough host memory. These sizes make long runs expensive. Start with
short schedules rather than assuming this generic formulation is competitive
with a specialized scheduling solver.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 1
./build/jsplib_jobshop_example "${ANNEALER_WORK_DIR}/ft06.txt" --seed 42 --trials 1 \
  --initial-temperature 30 --final-temperature 0.05 --levels 2 --sweeps 1
./build/jsplib_jobshop_example "${ANNEALER_WORK_DIR}/la01.txt" --memory-mib 2048 --seed 42 --trials 1 \
  --initial-temperature 30 --final-temperature 0.05 --levels 2 --sweeps 1
ctest --test-dir build --output-on-failure -R 'jobshop'
```

The runner uses prepared evaluation and the default single-variable replacement
neighborhood. Every trial starts from the same feasible serial schedule; seeds
increment consecutively. Moving one start independently often breaks tightly
coupled constraints, so progress may be poor. Swap mode is unsuitable for these
heterogeneous time domains. No optimality claim is made. The
[JSPLib reference listing](https://scheduleopt.github.io/benchmarks/jsplib)
reports optimal makespans 55 for ft06 and 666 for la01; tests do not require
annealing to reach them.

Each trial prints seed, independently computed actual makespan, individual
constraint counts, declared makespan, energy, feasibility and solve time. The
model estimate includes variables and terms. Best and final assignments are
checked independently; a feasible declared M may exceed the actual makespan,
so both are reported. The total time includes estimation, construction,
preparation, trials, verification and output, excluding parsing. Exit status zero
means the run and verification completed. `--help` lists schedule/trial options.
