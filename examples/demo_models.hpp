#pragma once
#include <annealer/solver.hpp>

namespace demos {
inline double identity(annealer::Value x) { return static_cast<double>(x); }
inline double square(annealer::Value x) { const double v = static_cast<double>(x); return v * v; }
inline double shift(annealer::Value x) { return static_cast<double>(x) + 1; }

inline annealer::Model binary_hubo() {
    annealer::ModelBuilder b;
    for (int i = 0; i < 3; ++i) { b.add_variable(annealer::Domain::cardinality(2)); }
    b.set_constant(1);
    b.add_term(-2, {{0, identity}});
    b.add_term(-1, {{1, identity}});
    b.add_term(-3, {{0, identity}, {1, identity}, {2, identity}});
    b.add_term(2, {{0, identity}, {2, identity}});
    return b.build();
}

inline annealer::Model nonlinear() {
    annealer::ModelBuilder b;
    b.add_variable(annealer::Domain{-2, 1, 4});
    b.add_variable(annealer::Domain::cardinality(3));
    b.add_variable(annealer::Domain{-1, 3});
    b.add_variable(annealer::Domain{2});
    b.set_constant(2);
    b.add_term(-0.5, {annealer::native_factor<square>(0), {1, shift},
                       {2, square}, {0, identity}, {3, identity}});
    b.add_term(1, {{0, square}});
    b.add_term(0, {{2, identity}});
    return b.build();
}
} // namespace demos
