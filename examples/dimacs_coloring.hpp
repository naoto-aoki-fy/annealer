#pragma once
#include "gset_maxcut.hpp"
#include <set>
#include <sstream>

namespace coloring {
struct Edge { std::size_t u, v; };
struct Graph { std::size_t vertices = 0; std::vector<Edge> edges; };

// Strict simple, undirected DIMACS text: comments, one header, then edges.
inline Graph read(std::istream& input) {
    Graph graph;
    std::size_t expected = 0;
    std::set<std::pair<std::size_t, std::size_t>> seen;
    std::string line;
    while (std::getline(input, line)) {
        std::istringstream row(line);
        std::string kind, a, b, extra;
        if (!(row >> kind) || kind == "c") continue;
        if (kind == "p") {
            std::string format;
            if (graph.vertices || !(row >> format >> a >> b) || format != "edge" || (row >> extra))
                throw std::invalid_argument("Expected one p edge N M header");
            graph.vertices = gset::unsigned_number<std::size_t>(a);
            expected = gset::unsigned_number<std::size_t>(b);
            if (!graph.vertices) throw std::invalid_argument("Vertex count must be positive");
        } else if (kind == "e") {
            if (!graph.vertices || !(row >> a >> b) || (row >> extra))
                throw std::invalid_argument("Invalid edge record or missing header");
            auto u = gset::unsigned_number<std::size_t>(a);
            auto v = gset::unsigned_number<std::size_t>(b);
            if (!u || !v || u > graph.vertices || v > graph.vertices || u == v)
                throw std::invalid_argument("Invalid vertex index or self-loop");
            if (u > v) std::swap(u, v);
            if (!seen.emplace(u, v).second) throw std::invalid_argument("Duplicate undirected edge");
            graph.edges.push_back({u - 1, v - 1});
            if (graph.edges.size() > expected) throw std::invalid_argument("Too many edges");
        } else throw std::invalid_argument("Unknown DIMACS text record (binary .col.b requires conversion)");
    }
    if (input.bad()) throw std::runtime_error("Graph input read failed");
    if (!graph.vertices || graph.edges.size() != expected) throw std::invalid_argument("Missing header or edge count mismatch");
    return graph;
}
inline void validate_colors(std::size_t colors) {
    if (!colors || static_cast<std::uintmax_t>(colors - 1) >
        static_cast<std::uintmax_t>(std::numeric_limits<annealer::Value>::max()))
        throw std::invalid_argument("Colors must be positive and fit the integer domain");
}
inline annealer::Model make_model(const Graph& graph, std::size_t colors) {
    validate_colors(colors);
    annealer::ModelBuilder builder;
    for (std::size_t i = 0; i < graph.vertices; ++i)
        builder.add_variable(annealer::Domain::cardinality(colors));
    for (const auto& edge : graph.edges) {
        for (std::size_t c = 0; c < colors; ++c) {
            const auto indicator = [color = static_cast<annealer::Value>(c)](annealer::Value x) {
                return x == color ? 1.0 : 0.0;
            };
            builder.add_term(1.0, {{edge.u, indicator}, {edge.v, indicator}});
        }
    }
    return builder.build();
}
// Independent of factor evaluation and preparation.
inline std::size_t conflicts(const Graph& graph, std::size_t colors, const annealer::Assignment& x) {
    validate_colors(colors);
    if (x.size() != graph.vertices) throw std::invalid_argument("Assignment size mismatch");
    for (auto value : x)
        if (value < 0 || static_cast<std::uintmax_t>(value) >= colors)
            throw std::invalid_argument("Assignment color out of range");
    std::size_t count = 0;
    for (const auto& edge : graph.edges) count += x.at(edge.u) == x.at(edge.v);
    return count;
}
} // namespace coloring
