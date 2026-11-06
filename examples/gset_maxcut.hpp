#pragma once

#include <annealer/model.hpp>
#include <charconv>
#include <istream>
#include <string>

namespace gset {
struct Edge { std::size_t u, v; double weight; };
struct Graph { std::size_t vertices; std::vector<Edge> edges; };

// Reject signs, partial tokens, and values outside the destination range.
template<class T> T unsigned_number(const std::string& text) {
    T value{};
    const auto result = std::from_chars(text.data(), text.data() + text.size(), value);
    if (result.ec != std::errc{} || result.ptr != text.data() + text.size()) {
        throw std::invalid_argument("Invalid unsigned integer: " + text);
    }
    return value;
}

inline Graph read(std::istream& input) {
    std::string n, m;
    if (!(input >> n >> m)) throw std::invalid_argument("Missing graph header (vertices edges)");
    Graph graph{unsigned_number<std::size_t>(n), {}};
    const auto count = unsigned_number<std::size_t>(m);
    if (!graph.vertices) throw std::invalid_argument("Vertex count must be positive");
    for (std::size_t k = 0; k < count; ++k) {
        std::string a, b;
        double weight;
        if (!(input >> a >> b >> weight)) {
            throw std::invalid_argument("Missing or invalid edge " + std::to_string(k + 1));
        }
        const auto u = unsigned_number<std::size_t>(a);
        const auto v = unsigned_number<std::size_t>(b);
        if (!u || !v || u > graph.vertices || v > graph.vertices || u == v || !std::isfinite(weight)) {
            throw std::invalid_argument("Invalid endpoints or weight at edge " + std::to_string(k + 1));
        }
        graph.edges.push_back({u - 1, v - 1, weight});
    }
    input >> std::ws;
    if (!input.eof()) throw std::invalid_argument("Unexpected data after declared edges");
    return graph;
}

inline double spin(std::int64_t x) {
    return 1.0 - 2.0 * static_cast<double>(x);
}

inline annealer::Model make_model(const Graph& graph) {
    annealer::ModelBuilder builder;
    for (std::size_t i = 0; i < graph.vertices; ++i) {
        builder.add_variable(annealer::Domain::cardinality(2));
    }
    double total = 0;
    for (const auto& e : graph.edges) {
        total += e.weight;
        if (!std::isfinite(total)) throw std::overflow_error("Total edge weight overflow");
        builder.add_term(0.5 * e.weight, {{e.u, spin}, {e.v, spin}});
    }
    builder.set_constant(-0.5 * total);
    return builder.build();
}

inline double cut_weight(const Graph& graph, const annealer::Assignment& x) {
    if (x.size() != graph.vertices) throw std::invalid_argument("Assignment size mismatch");
    for (auto value : x) {
        if (value != 0 && value != 1) throw std::invalid_argument("Assignment must be binary");
    }
    double cut = 0;
    for (const auto& e : graph.edges) {
        if (x[e.u] != x[e.v]) cut += e.weight;
        if (!std::isfinite(cut)) throw std::overflow_error("Cut weight overflow");
    }
    return cut;
}
} // namespace gset
