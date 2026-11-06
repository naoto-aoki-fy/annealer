# Weighted partial MaxSAT with WCNF

Build `maxsat_example` using the README CMake commands. This serial CPU example
uses the existing factorized higher-order model without auxiliary variables.

## Formulation and exactness

Each Boolean variable is an integer variable with domain `{0,1}`. A positive
literal's falsity factor is `1-x`; a negative literal's is `x`. Each clause is
one term whose factors are multiplied, with coefficient equal to its soft
weight or the hard penalty `P`. Thus a five-literal clause has five factors.
Callbacks capture only an immutable sign and are deterministic and side-effect-free.
Repeated literals remain repeated factors; opposite literals make a tautology's
violation zero. Empty clauses use the empty product (one), so are always violated.

The energy is `soft_cost + P * hard_violations`, where `P = S + 1` and `S`
is the sum of all positive integer soft weights. If the hard clauses are
jointly satisfiable, every global energy minimum is feasible and minimizes soft
cost among feasible assignments: all feasible energies are at most `S`, while
any infeasible energy is at least `P`. If hard clauses are inconsistent there
is no valid MaxSAT solution. Annealing does not establish feasibility or
optimality in general; a feasible candidate is not an optimality certificate.

The parser uses checked unsigned integer conversions and rejects instances whose
soft sum, penalty, or conservative bound `S + P * number_of_hard_clauses`
exceeds `2^53`. This ensures all possible model energies are exact integers in
`double` (and prevents integer overflow), including incremental energy updates.
Large weights can therefore be rejected even if they fit in 64 bits. Hard
weights in traditional input are used only to classify clauses, not as penalties.
The independent evaluator checks Boolean assignments and computes soft cost and
hard violations directly by testing clause satisfaction.

## Supported input

Plain text WCNF supports:

- Traditional `p wcnf N M TOP`, followed by `weight literals... 0` clauses.
  Weights at least `TOP` mark hard clauses; `TOP` must be positive.
- Modern headerless clauses: `h literals... 0` for hard clauses and
  `weight literals... 0` for soft clauses. Variable count is the largest index.
- An optional `p wcnf N M` header for explicit-hard or all-soft input.

Blank lines and lines starting with `c` after whitespace are ignored. Clauses
may span lines, and multiple terminated clauses may share a line. Every clause
must terminate with `0`. Weights must be positive decimal integers, indices must
be nonzero and within any declared variable count, and declared clause counts
must match. Mixed `h` and traditional TOP notation, duplicate headers, extra
records, signed weights, floating weights, and invalid tokens are rejected.
Compressed input must be decompressed first. No third-party parser is used.

## Reproduction

Create a small instance locally:

Set `ANNEALER_WORK_DIR` to an absolute directory of your choice before running
these commands. It defaults to `TMPDIR` when set, otherwise `/tmp`. Run the
following blocks in the same shell.

```sh
export ANNEALER_WORK_DIR="${ANNEALER_WORK_DIR:-${TMPDIR:-/tmp}}"
mkdir -p "${ANNEALER_WORK_DIR}"
cat > "${ANNEALER_WORK_DIR}/small.wcnf" <<'WCNF'
p wcnf 5 4 100
100 1 -2 3 0
3 -1 2 -3 4 -5 0
7 1 0
2 -1 0
WCNF
./build/maxsat_example "${ANNEALER_WORK_DIR}/small.wcnf" --seed 42 --trials 10 \
  --initial-temperature 30 --final-temperature 0.05 --levels 200 --sweeps 100
ctest --test-dir build --output-on-failure -R maxsat
```

`--help` lists defaults and options. Trials use consecutive seeds, with seed
range overflow rejected. Each trial prints total energy, independent soft cost,
hard violations, feasibility, solve seconds, and seed. The summary identifies the
lowest-energy trial; infeasible assignments are explicitly labeled. Exit status
zero means the run completed and verification succeeded, not that a feasible or
optimal solution was found. Total time includes model construction, preparation,
trials, verification and output, excluding input parsing. Temperature and sweep
settings need tuning for each instance and weight scale.

## Benchmark acquisition

The [MaxSAT Evaluation complete benchmark collection](https://zenodo.org/records/21283148)
provides year-grouped ZIP archives containing individually xz-compressed instances.
Choose a suitable small weighted partial instance from an archive, check the
published archive checksum, extract it locally, then decompress the chosen file:

```sh
xz -dk chosen.wcnf.xz
./build/maxsat_example chosen.wcnf --seed 42 --trials 10 \
  --initial-temperature 30 --final-temperature 0.05 --levels 200 --sweeps 100
```

Record the archive year, instance filename, checksum, build configuration, and
printed schedule/seeds when reporting results. The collection is large; no
archive download is part of the build or tests. Tests use small self-contained
instances and exhaustive enumeration rather than requiring stochastic optimality.
