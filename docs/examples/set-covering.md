# OR-Library set covering: scp41

This serial CPU example creates one binary integer variable `x_j` per column.
Selecting a column incurs its nonnegative integer cost. Each row contributes
one higher-order interaction with a `(1-x_j)` factor for every column covering
that row:

`E(x) = sum_j cost_j*x_j + P*sum_i product_{j in C_i}(1-x_j)`.

A row with five covering columns produces one five-factor term. There are no
auxiliary variables or pairwise QUBO conversions. All factors are deterministic,
side-effect-free functions of one variable. A separate Boolean evaluator sums
selected-column costs and counts uncovered rows, independently of model factors.

## Penalty and supported input

Let `S` be the sum of all column costs. The example uses `P = S+1`.
Every supported row has at least one covering column, so selecting all columns
is feasible and costs `S`. Nonnegative costs imply every feasible assignment has
energy at most `S`, while any infeasible assignment has energy at least `P`.
Thus every global energy minimum is a minimum-cost feasible cover. Simulated
annealing does not guarantee that it will reach that minimum or find a cover.

Costs are decimal nonnegative integers, including zero. Floating point and
negative costs are rejected. Checked unsigned parsing and the conservative bound
`S + P*row_count <= 2^53` ensure exact integer energies in `double` and avoid
integer overflow; some large instances or weights are therefore rejected.

The [official OR-Library format](https://people.brunel.ac.uk/~mastjjb/jeb/orlib/scpinfo.html)
is whitespace-separated and row-oriented:

1. Row count `m` and column count `n`, both positive.
2. The `n` column costs.
3. For each of the `m` rows, its covering-column count followed by that many
   1-based column indices.

Lists may span lines. The parser rejects missing or extra fields, invalid sizes
or indices, duplicate indices within a row, and empty rows (which cannot be
covered). This is the set covering format; do not supply column-oriented set
partitioning data. Repeated coverage lists for different rows are allowed.

## Download and reproduce

Download the [official scp41 instance](https://people.brunel.ac.uk/~mastjjb/jeb/orlib/files/scp41.txt)
locally. It has 200 rows and 1,000 columns. Benchmark data is not bundled and is
not downloaded by the build or tests; consult the provider about redistribution
before sharing data copies.

Set `ANNEALER_WORK_DIR` to an absolute directory of your choice before running
these commands. It defaults to `TMPDIR` when set, otherwise `/tmp`. Run the
following blocks in the same shell.

```sh
export ANNEALER_WORK_DIR="${ANNEALER_WORK_DIR:-${TMPDIR:-/tmp}}"
mkdir -p "${ANNEALER_WORK_DIR}"
curl -fL https://people.brunel.ac.uk/~mastjjb/jeb/orlib/files/scp41.txt -o "${ANNEALER_WORK_DIR}/scp41.txt"
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 1
./build/set_covering_example "${ANNEALER_WORK_DIR}/scp41.txt" --seed 42 --trials 10 \
  --initial-temperature 30 --final-temperature 0.05 --levels 200 --sweeps 100
ctest --test-dir build --output-on-failure -R set_covering
```

`--help` lists defaults. Seeds increase consecutively across trials; overflowing
seed ranges are rejected. Schedule parameters are configurable and should be
tuned for the cost and penalty scale. No reference optimum is assumed.

Each trial reports this structure (values depend on the run):

```text
Seed: ... Selected-column cost: ... Uncovered rows: ... Total energy: ... Solve seconds: ... Feasible: yes/no Status: ...
Best feasible cost: ... Seed: ...
Lowest energy seed: ... Lowest energy: ...
Total seconds (model, preparation, trials, verification, output): ...
```

If none of the trial-best assignments is feasible, the best feasible cost is
reported as unavailable. An infeasible assignment's selected-column cost is not
a valid covering cost. Feasible candidates are labeled with optimality unproven.
Both best and final assignments are checked against the independent objective.
Exit status zero means execution and verification succeeded, not optimality.
Solve time covers each solver call; total time includes model construction,
preparation, trials, verification and output, excluding parsing. Record the input
checksum, compiler/build configuration, schedules and seeds with measurements.

## Memory and runtime

Model size grows with column count plus the total number of row incidences.
Long rows remain long product terms; evaluation and updates touch their factors
and incident terms, so dense rows can increase work substantially. Preparation
uses a 1 MiB factor lookup budget in this executable, not an exponential table
over all subsets of a row. Runtime also grows with columns, temperature levels,
sweeps and trials. Large rail instances may require substantial memory and time;
start with scp41 and a short schedule. See [preparation](../preparation.md) for
lookup behavior. Tests use exhaustive small synthetic instances with several
lookup budgets and do not require annealing to find the scp41 optimum.
