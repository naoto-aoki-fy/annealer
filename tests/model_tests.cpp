#include <annealer/model.hpp>

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace {
int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}

template <class Exception, class F>
void throws(F action, const std::string& message) {
    try {
        action();
        check(false, message + " (no exception)");
    } catch (const Exception&) {
    } catch (...) {
        check(false, message + " (wrong exception)");
    }
}

double identity(std::int64_t value) { return static_cast<double>(value); }
struct Square {
    double operator()(std::int64_t value) const {
        const double real = static_cast<double>(value);
        return real * real;
    }
};

void model_and_evaluation() {
    using namespace annealer;
    ModelBuilder builder;
    check(full_energy(builder.build(), {}) == 0.0, "empty model");
    builder.set_constant(1.5);
    builder.add_term(2.0, {});
    check(full_energy(builder.build(), {}) == 3.5, "constant-only model");

    const auto x = builder.add_variable(Domain{-2, 1, 4});
    const auto y = builder.add_variable(Domain::cardinality(3));
    const auto z = builder.add_variable(Domain{3, 7});
    builder.add_term(-0.5, {{x, Square{}}, {y, [](Value v) { return double(v) + 1.0; }},
                            {z, identity}});
    builder.add_term(2.0, {{x, identity}});
    const auto snapshot = builder.build();
    check(full_energy(snapshot, {-2, 2, 3}) == -18.5,
          "mixed domains, nonlinear function and cubic term");
    for (const Value a : {-2, 1, 4}) {
        for (const Value b : {0, 1, 2}) {
            for (const Value c : {3, 7}) {
                const double expected = 3.5 - 0.5 * double(a) * double(a) *
                                                  (double(b) + 1.0) * double(c) + 2.0 * double(a);
                check(full_energy(snapshot, {a, b, c}) == expected,
                      "exhaustive multivalued objective matches independent formula");
            }
        }
    }
    builder.set_constant(100.0);
    builder.add_term(50.0, {});
    builder.add_variable(Domain{9});
    check(full_energy(snapshot, {-2, 2, 3}) == -18.5, "builder snapshot ownership");
    check(snapshot.domains()[0].values() == std::vector<Value>({-2, 1, 4}),
          "domain order is preserved");
    check(snapshot.domains()[1].values() == std::vector<Value>({0, 1, 2}),
          "cardinality represents zero through q minus one");

    ModelBuilder repeated;
    const auto r = repeated.add_variable(Domain{-2, 3});
    repeated.add_term(0.5, {{r, identity}, {r, Square{}}, {r, identity}});
    check(full_energy(repeated.build(), {-2}) == 8.0, "all repeated occurrences evaluated");

    ModelBuilder high_order;
    high_order.add_variable(Domain{2});
    std::vector<Factor> factors;
    for (int i = 0; i < 64; ++i) {
        factors.emplace_back(0, identity);
    }
    high_order.add_term(1.0, std::move(factors));
    check(full_energy(high_order.build(), {2}) == std::ldexp(1.0, 64),
          "interaction order is not restricted to a small bound");

    ModelBuilder hubo;
    for (int i = 0; i < 3; ++i) {
        hubo.add_variable(Domain::cardinality(2));
    }
    hubo.add_term(-3, {{0, identity}, {1, identity}, {2, identity}});
    for (Value a = 0; a < 2; ++a) {
        for (Value b = 0; b < 2; ++b) {
            for (Value c = 0; c < 2; ++c) {
                check(full_energy(hubo.build(), {a, b, c}) == -3.0 * a * b * c,
                      "binary HUBO truth table");
            }
        }
    }

    // Integer extremes reach callbacks intact, without index substitution.
    const auto low = std::numeric_limits<Value>::min();
    const auto high = std::numeric_limits<Value>::max();
    ModelBuilder extremes;
    extremes.add_variable(Domain{high, low});
    extremes.add_term(1, {{0, [low](Value v) { return v == low ? 7.0 : 9.0; }}});
    check(full_energy(extremes.build(), {low}) == 7.0, "minimum int64 actual value");
    check(full_energy(extremes.build(), {high}) == 9.0, "maximum int64 actual value");
}

void validation() {
    using namespace annealer;
    throws<std::invalid_argument>([] { Domain empty{}; }, "empty explicit domain");
    throws<std::invalid_argument>([] { Domain::cardinality(0); }, "zero cardinality");
    throws<std::invalid_argument>([] { Domain duplicate{-1, 2, -1}; }, "duplicate domain");
    if (std::numeric_limits<std::size_t>::max() >
        static_cast<std::uintmax_t>(std::numeric_limits<Value>::max())) {
        throws<std::invalid_argument>([] { Domain::cardinality(std::numeric_limits<std::size_t>::max()); },
                                      "unrepresentable cardinality");
    }
    ModelBuilder builder;
    builder.add_variable(Domain{-5, 8});
    throws<std::out_of_range>([&] { builder.add_term(1, {{1, identity}}); }, "unknown variable");
    throws<std::out_of_range>([&] {
        builder.add_term(1, {{std::numeric_limits<Variable>::max(), identity}});
    }, "maximum variable index");
    throws<std::invalid_argument>([] {
        Factor factor(0, static_cast<double (*)(Value)>(nullptr));
    }, "null function pointer");
    throws<std::invalid_argument>([] {
        Factor factor(0, std::function<double(Value)>{});
    }, "empty function wrapper");
    Domain source{1, 2};
    Domain destination(std::move(source));
    // A moved-from domain is valid but can be empty. Never accept an empty one.
    if (source.values().empty()) {
        throws<std::invalid_argument>([&] { builder.add_variable(source); }, "moved-from empty domain");
    }
    Factor factor_source(0, identity);
    Factor factor_destination(std::move(factor_source));
    // std::function's moved-from state is unspecified; rejection is required if empty.
    try {
        factor_source.evaluate(1);
    } catch (const std::bad_function_call&) {
        throws<std::invalid_argument>([&] { builder.add_term(1.0, {factor_source}); },
                                      "moved-from empty factor");
    }
    const auto model = builder.build();
    throws<std::invalid_argument>([&] { full_energy(model, {}); }, "short assignment");
    throws<std::invalid_argument>([&] { full_energy(model, {-5, 8}); }, "long assignment");
    throws<std::invalid_argument>([&] { full_energy(model, {0}); }, "index is not domain value");
    throws<std::invalid_argument>([&] { full_energy(model, {9}); }, "out-of-domain value");
    for (const double value : {std::numeric_limits<double>::infinity(),
                               -std::numeric_limits<double>::infinity(),
                               std::numeric_limits<double>::quiet_NaN()}) {
        throws<std::invalid_argument>([&] { builder.set_constant(value); }, "nonfinite constant");
        throws<std::invalid_argument>([&] { builder.add_term(value, {}); }, "nonfinite coefficient");
        ModelBuilder bad_factor;
        bad_factor.add_variable(Domain{0});
        bad_factor.add_term(0.0, {{0, [value](Value) { return value; }}});
        throws<std::domain_error>([&] { full_energy(bad_factor.build(), {0}); },
                                  "nonfinite factor even with zero coefficient");
    }
    check(full_energy(builder.build(), {-5}) == 0.0, "failed builder operations preserve model");

    ModelBuilder throwing;
    throwing.add_variable(Domain{0});
    throwing.add_term(1, {{0, [](Value) -> double { throw std::logic_error("user failure"); }}});
    throws<std::logic_error>([&] { full_energy(throwing.build(), {0}); }, "callback exceptions propagate");
    throws<std::invalid_argument>([&] { full_energy(throwing.build(), {1}); },
                                  "assignment validated before callbacks");
}

void arithmetic() {
    using namespace annealer;
    const double max = std::numeric_limits<double>::max();
    ModelBuilder product;
    product.add_variable(Domain{0});
    const auto huge = [max](Value) { return max; };
    product.add_term(0.0, {{0, huge}, {0, huge}});
    throws<std::overflow_error>([&] { full_energy(product.build(), {0}); }, "product overflow");
    ModelBuilder scaled;
    scaled.add_variable(Domain{0});
    scaled.add_term(2.0, {{0, huge}});
    throws<std::overflow_error>([&] { full_energy(scaled.build(), {0}); }, "coefficient overflow");
    ModelBuilder sum;
    sum.set_constant(-max);
    sum.add_term(-max, {});
    throws<std::overflow_error>([&] { full_energy(sum.build(), {}); }, "negative sum overflow");
    ModelBuilder positive_sum;
    positive_sum.set_constant(max);
    positive_sum.add_term(max, {});
    throws<std::overflow_error>([&] { full_energy(positive_sum.build(), {}); }, "positive sum overflow");
    ModelBuilder underflow;
    underflow.add_variable(Domain{0});
    underflow.add_term(1e-300, {{0, [](Value) { return 1e-300; }}});
    check(full_energy(underflow.build(), {0}) == 0.0, "finite underflow is permitted");
}
} // namespace

int main() {
    try {
        model_and_evaluation();
        validation();
        arithmetic();
    } catch (const std::exception& error) {
        std::cerr << "Unexpected exception: " << error.what() << '\n';
        return 1;
    }
    if (failures != 0) {
        return 1;
    }
    std::cout << "All model tests passed\n";
}
