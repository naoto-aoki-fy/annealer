#pragma once
#include "gset_maxcut.hpp"
#include <memory>
#include <numeric>
#include <random>
#include <sstream>
#include <map>

namespace tsp {
using Cost = std::uint64_t;
inline constexpr Cost limit = Cost{1} << 53;
struct Point { double x, y; };
inline Cost distance(Point a, Point b) {
    const double dx = a.x - b.x, dy = a.y - b.y;
    const double rounded = std::floor(std::sqrt(dx*dx + dy*dy) + 0.5);
    if (!std::isfinite(rounded) || rounded < 0 || rounded > static_cast<double>(limit))
        throw std::invalid_argument("TSPLIB distance exceeds exact integer range");
    return static_cast<Cost>(rounded);
}
struct Instance { std::size_t size; std::vector<Point> points; Cost penalty; };
inline std::string trim(const std::string& s) {
    const auto begin = s.find_first_not_of(" \t\r\n");
    return begin == std::string::npos ? "" : s.substr(begin, s.find_last_not_of(" \t\r\n") - begin + 1);
}
inline Instance read(std::istream& input) {
    std::map<std::string, std::string> headers;
    Instance instance{0, {}, 0};
    bool section = false, ended = false;
    std::vector<bool> seen;
    std::string line;
    while (std::getline(input, line)) {
        line = trim(line);
        if (line.empty()) continue;
        if (ended) throw std::invalid_argument("Trailing TSPLIB data");
        if (!section) {
            if (line == "NODE_COORD_SECTION") {
                if (headers["TYPE"] != "TSP" || headers["EDGE_WEIGHT_TYPE"] != "EUC_2D")
                    throw std::invalid_argument("Only TYPE TSP with EUC_2D is supported");
                instance.size = gset::unsigned_number<std::size_t>(headers["DIMENSION"]);
                if (instance.size < 2 || instance.size > std::vector<Point>{}.max_size() ||
                    instance.size > std::vector<Cost>{}.max_size() / instance.size ||
                    instance.size > static_cast<std::uintmax_t>(std::numeric_limits<annealer::Value>::max()))
                    throw std::invalid_argument("Invalid TSPLIB dimension");
                section = true;
            } else {
                const auto colon = line.find(':');
                std::string key, value;
                if (colon != std::string::npos) { key = trim(line.substr(0, colon)); value = trim(line.substr(colon + 1)); }
                else { std::istringstream row(line); row >> key; std::getline(row, value); value = trim(value); }
                if (value.empty() || (key != "NAME" && key != "COMMENT" && key != "TYPE" && key != "DIMENSION" && key != "EDGE_WEIGHT_TYPE"))
                    throw std::invalid_argument("Unsupported TSPLIB header or section");
                if (key != "COMMENT" && headers.contains(key)) throw std::invalid_argument("Duplicate TSPLIB header");
                headers[key] = value;
            }
        } else if (line == "EOF") ended = true;
        else {
            std::istringstream row(line);
            std::string id, sx, sy, extra;
            if (!(row >> id >> sx >> sy) || (row >> extra)) throw std::invalid_argument("Invalid coordinate record");
            const auto index = gset::unsigned_number<std::size_t>(id);
            if (!index || index > instance.size) throw std::invalid_argument("City ID outside dimension");
            auto coordinate = [](const std::string& token) {
                std::size_t end = 0; const double value = std::stod(token, &end);
                if (end != token.size() || !std::isfinite(value)) throw std::invalid_argument("Invalid coordinate");
                return value;
            };
            const Point point{coordinate(sx), coordinate(sy)};
            // Allocate from the dimension only after a valid record is present.
            if (seen.empty()) { instance.points.resize(instance.size); seen.resize(instance.size); }
            if (seen[index - 1]) throw std::invalid_argument("Duplicate city ID");
            seen[index - 1] = true; instance.points[index - 1] = point;
        }
    }
    if (input.bad()) throw std::runtime_error("TSPLIB input read failed");
    if (!section || seen.empty() || std::find(seen.begin(), seen.end(), false) != seen.end())
        throw std::invalid_argument("Missing coordinate section or cities");
    Cost maximum = 0;
    for (auto a : instance.points) for (auto b : instance.points) maximum = std::max(maximum, distance(a, b));
    const auto n = instance.size;
    if (maximum > (limit - 1) / n) throw std::invalid_argument("Tour bound exceeds exact range");
    const auto upper = maximum * n;
    instance.penalty = upper + 1;
    if (n > limit / (n - 1)) throw std::invalid_argument("Pair count exceeds exact range");
    const auto pairs = n * (n - 1) / 2;
    if (pairs > (limit - upper) / instance.penalty) throw std::invalid_argument("Penalized energy exceeds exact range");
    return instance;
}
inline annealer::Model make_model(const Instance& instance) {
    const auto n = instance.size;
    annealer::ModelBuilder builder;
    for (std::size_t t = 0; t < n; ++t) builder.add_variable(annealer::Domain::cardinality(n));
    auto values = std::make_shared<std::vector<Cost>>();
    for (auto a : instance.points) for (auto b : instance.points) values->push_back(distance(a, b));
    const std::shared_ptr<const std::vector<Cost>> distances = values;
    auto indicator = [](std::size_t a) { return [a](annealer::Value x) { return x == static_cast<annealer::Value>(a) ? 1.0 : 0.0; }; };
    for (std::size_t t = 0; t < n; ++t) {
        for (std::size_t a = 0; a < n; ++a)
            builder.add_term(1, {{t, indicator(a)}, {(t+1)%n, [distances, a, n](annealer::Value x) { return static_cast<double>((*distances)[a*n+x]); }}});
        for (std::size_t u = t + 1; u < n; ++u)
            for (std::size_t a = 0; a < n; ++a)
                builder.add_term(static_cast<double>(instance.penalty), {{t, indicator(a)}, {u, indicator(a)}});
    }
    return builder.build();
}
struct Objective { Cost cost = 0, collisions = 0, penalty_cost = 0; bool feasible = false; };
inline Objective evaluate(const Instance& instance, const annealer::Assignment& x) {
    if (x.size() != instance.size) throw std::invalid_argument("Tour size mismatch");
    for (auto v : x) if (v < 0 || static_cast<std::uint64_t>(v) >= instance.size) throw std::invalid_argument("City outside domain");
    Objective result;
    for (std::size_t t = 0; t < x.size(); ++t) {
        result.cost += distance(instance.points[x[t]], instance.points[x[(t+1)%x.size()]]);
        for (std::size_t u = t + 1; u < x.size(); ++u) result.collisions += x[t] == x[u];
    }
    result.penalty_cost = instance.penalty * result.collisions;
    result.feasible = result.collisions == 0;
    return result;
}
inline double energy(const Instance&, Objective objective) { return static_cast<double>(objective.cost + objective.penalty_cost); }
inline annealer::Assignment permutation(std::size_t n, std::uint64_t seed) {
    annealer::Assignment x(n); std::iota(x.begin(), x.end(), annealer::Value{0});
    std::mt19937_64 random(seed); std::shuffle(x.begin(), x.end(), random); return x;
}
} // namespace tsp
