# Model API and reference evaluation

The library minimizes objectives of the form

`E(x) = c0 + sum_t (c_t * product_r f_t,r(x_i(t,r)))`.

`annealer::Value` is `std::int64_t`; coefficients, function results, and energy
are `double`. `Variable` is a zero-based `std::size_t` identifier returned by
`ModelBuilder::add_variable`. It identifies a variable, not a value in a domain.

## Construction

- `Domain{-2, 1, 4}` or `Domain(std::vector<Value>{...})` supplies distinct
  integer values. Input order is retained. A domain must be nonempty.
- `Domain::cardinality(q)` constructs `{0, ..., q - 1}`. Zero or a cardinality
  whose largest value cannot fit in `Value` is rejected. Domains are stored
  explicitly, so available memory also limits cardinality.
- `ModelBuilder::add_variable(domain)` registers a variable. Domains can differ.
- `set_constant(c0)` sets the offset, initially zero.
- `add_term(coefficient, factors)` appends a term. Each `Factor(variable, function)`
  references an already registered variable. Functions may differ per occurrence.
- `build()` returns an owning model snapshot with read-only accessors. Subsequent
  builder edits do not change that snapshot.
- `full_energy(model, assignment)` evaluates the objective. `Assignment` is a
  `std::vector<Value>` in variable registration order, containing **actual domain
  values**, not domain indices.

An empty factor list has product one, so `add_term(c, {})` adds a constant.
One factor gives a unary term; arbitrary factor counts are supported subject to
memory limits. Repeated occurrences of a variable are evaluated individually,
even if they use different functions. Binary HUBO uses cardinality-two domains
and identity functions without a separate encoding. Zero-variable models and
empty assignments are valid.

Models own their domains, coefficients, terms, and callback wrappers. No mutable
assignment or future solver state is stored in a model. Keep the model unchanged
for the duration of any evaluation or future solve. Do not replace or move it
while an operation uses it. Builder preparation and evaluation are serial CPU
operations; there are no worker threads or parallel execution facilities.

## Function contract

A factor is conceptually `double f(std::int64_t value)`. Function pointers,
stateless lambdas, and copyable function objects callable through a const reference
are accepted. Immutable captured constants are also allowed. Return types must
be convertible to `double`. Null pointers and empty `std::function` wrappers are
rejected.

Functions must be deterministic and side-effect-free. They must not depend on
call history, randomness, clocks, I/O, mutable global data, or other assignment
variables. Purity and captured-resource lifetimes are the user's responsibility;
the library cannot prove them. Avoid captures by reference unless the referenced
immutable data outlives every model copy. User callback exceptions propagate.

Callbacks compile as native C/C++ code. The header-only API enables ordinary
compiler optimization; `std::function` supplies ownership and type erasure and may
retain an indirect call. There is no interpreter, JIT, symbolic restriction,
precomputation, or requirement for `constexpr` functions. Use a Release build for
optimization. Do not use fast-math or other options that assume all results are
finite: they can defeat the validation below. Integer overflow inside a callback
is the user's responsibility; convert to double before potentially overflowing
integer arithmetic, as in the square function in the example.

## Numerical and validation behavior

The reference evaluator validates the whole assignment before invoking callbacks.
It starts with `c0`, visits terms and factors in insertion order, multiplies all
factor values starting at one, multiplies the product by the coefficient, then
adds the contribution to the energy. It evaluates zero-coefficient terms too.
There are no cached products, tables, or incremental updates; this implementation
is intended as an oracle for future evaluators.

| Condition | Behavior |
| --- | --- |
| Empty/duplicate domain, invalid cardinality | `std::invalid_argument` |
| Nonfinite offset/coefficient or empty callback | `std::invalid_argument` |
| Unknown variable at term insertion | `std::out_of_range` |
| Wrong assignment length or value outside its domain | `std::invalid_argument` |
| Callback returns NaN or infinity | `std::domain_error` |
| Nonfinite product, scaled contribution, or sum | `std::overflow_error` |

Every intermediate must be finite, even if a later multiplication by zero or
cancellation would produce a finite mathematical answer. No reordering or
compensated summation is performed. Floating-point rounding, loss of integer
precision when converting large values to double, and finite underflow (including
rounding to zero) follow normal double arithmetic. Allocation failures propagate.
The library does not change the floating-point environment; use the compiler's
normal nontrapping environment for exception-based overflow reporting.

The [serial Metropolis solver](solver.md) uses this evaluator for every proposal,
keeping its mutable state separate from the model.
