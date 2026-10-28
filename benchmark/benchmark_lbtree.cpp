/**
 * @file benchmark_lbtree.cpp
 * @brief Benchmarks LBTree insertion order and containment classification.
 */

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "benchmark_common.hpp"
#include "local_bounds/bound_set.hpp"
#include "local_bounds/neighborhood_bound_set.hpp"
#include "structures/lb_tree.hpp"

using namespace local_bounds;

template <Objective Sense>
void benchmark_insertion_orders(const std::vector<LocalBound<int64_t>>& snapshot,
                                InstanceType type, unsigned shuffles) {
  std::vector<std::vector<int64_t>> native;
  for (const auto& bound : snapshot) native.push_back(bound.coordinates);
  auto expected_bounds = native;
  std::sort(expected_bounds.begin(), expected_bounds.end());
  const auto dimensions = native.front().size();
  constexpr bool minimize = Sense == Objective::MINIMIZE;

  std::vector<std::vector<int64_t>> queries;
  const auto samples = std::min(native.size(), std::size_t{128});
  for (std::size_t i = 0; i < samples; ++i) {
    auto query = native[i * native.size() / samples];
    queries.push_back(query);
    for (auto& value : query) value += minimize ? -1 : 1;
    queries.push_back(std::move(query));
  }
  std::vector<bool> expected_hits;
  for (const auto& query : queries)
    expected_hits.push_back(std::any_of(native.begin(), native.end(),
        [&](const auto& bound) { return strictly_dominates<int64_t, Sense>(query, bound); }));

  for (unsigned order = 0; order < shuffles + 2; ++order) {
    const unsigned seed = order < 2 ? 0 : 42 + order - 2;
    const std::string name = order == 0 ? "native" : order == 1 ? "sorted"
                                                               : "shuffled_" + std::to_string(seed);
    auto coordinates = native;
    if (order == 1) std::sort(coordinates.begin(), coordinates.end());
    else if (order >= 2) {
      std::mt19937_64 rng(seed);
      std::shuffle(coordinates.begin(), coordinates.end(), rng);
    }
    const double build_ms = median_time_ms([&] {
      LBTree<int64_t, Sense> index(32, 8, dimensions);
      for (const auto& point : coordinates) index.Insert(point);
      return index.Size();
    });
    LBTree<int64_t, Sense> index(32, 8, dimensions);
    for (const auto& point : coordinates) index.Insert(point);
    auto exported = index.GetAllBounds();
    std::sort(exported.begin(), exported.end());
    const auto context = std::to_string(dimensions) + "D " + (minimize ? "MIN " : "MAX ")
                       + instance_type_str(type) + " " + name;
    if (exported != expected_bounds || index.Size() != snapshot.size())
      throw std::runtime_error("LBTree bound mismatch: " + context);
    for (std::size_t i = 0; i < queries.size(); ++i)
      if ((index.FindStrictlyDominated(queries[i]) != nullptr) != expected_hits[i])
        throw std::runtime_error("LBTree query mismatch: " + context);
    const double query_ms = median_time_ms([&] {
      std::size_t hits = 0;
      for (const auto& query : queries) hits += index.FindStrictlyDominated(query) != nullptr;
      return hits;
    });
    std::cout << dimensions << ',' << (minimize ? "MIN" : "MAX") << ',' << instance_type_str(type) << ','
              << snapshot.size() << ',' << name << ',' << index.Depth() << ','
              << build_ms << ',' << query_ms << '\n' << std::flush;
  }
}

template <Objective Sense>
void run_benchmark(std::size_t dimensions, std::size_t target_points,
                   InstanceType type, unsigned shuffles) {
  const auto points = generate_stable_set<Sense>(target_points, dimensions, 10000, type, 42);
  constexpr bool minimize = Sense == Objective::MINIMIZE;
  BoundSet<int64_t, Sense> bounds(std::vector<int64_t>(dimensions, minimize ? 10001 : 0),
                                 std::vector<int64_t>(dimensions, minimize ? 0 : 10001));
  const std::vector<std::size_t> targets = {500, 1000, 3000, 5000};
  std::size_t next = 0;
  for (std::size_t i = 0; i < points.size(); ++i) {
    bounds.update_re_enhanced(points[i]);
    if (bounds.size() < targets[next] && i + 1 != points.size()) continue;
    benchmark_insertion_orders<Sense>(bounds.bounds(), type, shuffles);
    if (++next == targets.size()) break;
  }
}

template <Objective Sense>
std::vector<bool> classify_snapshot(const std::vector<LocalBound<int64_t>>& raw, bool indexed) {
  std::vector<bool> quasi(raw.size(), false);
  if (!indexed) {
    for (std::size_t i = 0; i < raw.size(); ++i)
      for (std::size_t j = 0; j < raw.size(); ++j)
        if (i != j && weakly_dominates<int64_t, Sense>(raw[i].coordinates, raw[j].coordinates) &&
            (raw[i].coordinates != raw[j].coordinates || j < i)) {
          quasi[i] = true;
          break;
        }
    return quasi;
  }
  LBTree<int64_t, Sense> index(32, 8, raw.front().coordinates.size());
  for (const auto& bound : raw) index.Insert(bound.coordinates);
  std::set<std::vector<int64_t>> seen;
  for (std::size_t i = 0; i < raw.size(); ++i)
    quasi[i] = !seen.insert(raw[i].coordinates).second || index.HasProperContainingBound(raw[i].coordinates);
  return quasi;
}

template <Objective Sense>
void run_classification_benchmark(std::size_t dimensions, std::size_t target_points,
                                  InstanceType type, unsigned seed) {
  const auto points = generate_stable_set<Sense>(target_points, dimensions, 10000, type, seed);
  constexpr bool minimize = Sense == Objective::MINIMIZE;
  NeighborhoodBoundSet<int64_t, Sense> bounds(std::vector<int64_t>(dimensions, minimize ? 10001 : 0),
                                             std::vector<int64_t>(dimensions, minimize ? 0 : 10001));
  const std::vector<std::size_t> targets = {500, 1000, 1500, 2000, 2500, 3000, 4000, 5000};
  std::size_t next = 0;
  for (const auto& point : points) {
    bounds.update(point);
    if (bounds.size() < targets[next]) continue;
    const auto raw = bounds.bounds();
    const auto expected = classify_snapshot<Sense>(raw, false);
    if (expected != classify_snapshot<Sense>(raw, true) || bounds.get_adjacency_graph(true).quasi != expected)
      throw std::runtime_error("Classification benchmark correctness failure");
    const auto expected_count = std::count(expected.begin(), expected.end(), false);
    if (bounds.nonredundant_size() != static_cast<std::size_t>(expected_count))
      throw std::runtime_error("Classification benchmark count mismatch");
    const double scan_ms = median_time_ms([&] {
      const auto flags = classify_snapshot<Sense>(raw, false);
      return std::count(flags.begin(), flags.end(), false);
    });
    const double index_ms = median_time_ms([&] {
      const auto flags = classify_snapshot<Sense>(raw, true);
      return std::count(flags.begin(), flags.end(), false);
    });
    const double auto_ms = median_time_ms([&] { return bounds.nonredundant_size(); });
    const double graph_ms = median_time_ms([&] { return bounds.get_adjacency_graph().nodes.size(); });
    std::cout << dimensions << ',' << (minimize ? "MIN" : "MAX") << ',' << instance_type_str(type) << ','
              << seed << ',' << targets[next] << ',' << bounds.size() << ',' << bounds.max_node_index() << ','
              << scan_ms << ',' << index_ms << ',' << auto_ms << ',' << graph_ms << '\n' << std::flush;
    if (++next == targets.size()) break;
  }
}

int main(int argc, char* argv[]) try {
  const bool classification = argc > 1 && std::string(argv[1]) == "--classification";
  const int count_arg = classification ? 2 : 1;
  if (argc > count_arg + 1)
    throw std::invalid_argument("Usage: benchmark_lbtree [shuffles | --classification [runs]]");
  std::size_t consumed = 0;
  const int count = argc > count_arg ? std::stoi(argv[count_arg], &consumed) : 3;
  if (count < 1 || (argc > count_arg && consumed != std::string(argv[count_arg]).size()))
    throw std::invalid_argument(classification ? "Classification runs must be a positive integer"
                                               : "Shuffle count must be a positive integer");
  if (classification)
    std::cout << "dimensions,sense,instance,seed,target_active,active,slots,scan_ms,index_ms,auto_count_ms,graph_ms\n";
  else
    std::cout << "dimensions,sense,instance,bounds,order,depth,build_ms,query_ms\n";
  for (const auto& [dimensions, points] : std::vector<std::pair<std::size_t, std::size_t>>{
           {3, 5000}, {4, 1000}, {5, 400}, {6, 200}})
    for (auto type : {InstanceType::GENERAL_POSITION, InstanceType::GENERAL_CASE}) {
      if (classification) {
        for (int run = 0; run < count; ++run) {
          run_classification_benchmark<Objective::MINIMIZE>(dimensions, points, type, 42 + run);
          run_classification_benchmark<Objective::MAXIMIZE>(dimensions, points, type, 42 + run);
        }
      } else {
        run_benchmark<Objective::MINIMIZE>(dimensions, points, type, count);
        run_benchmark<Objective::MAXIMIZE>(dimensions, points, type, count);
      }
    }
  return 0;
} catch (const std::exception& error) {
  std::cerr << error.what() << '\n';
  return 1;
}
