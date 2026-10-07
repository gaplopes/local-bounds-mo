#include <algorithm>
#include <stdexcept>
#define CHECK(condition) do { if (!(condition)) throw std::runtime_error(#condition); } while (false)
#include <cstdint>
#include <iostream>
#include <set>
#include <vector>

#include "local_bounds.hpp"

using namespace local_bounds;

using IntBound = LocalBound<int64_t>;
using IntPoint = Point<int64_t>;

std::vector<std::vector<int64_t>> canonicalize_bounds(
    const std::vector<IntBound>& bounds) {
  std::vector<std::vector<int64_t>> out;
  out.reserve(bounds.size());
  for (const auto& b : bounds) {
    out.push_back(b.coordinates);
  }
  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

bool compare_bounds(
    const std::vector<IntBound>& a,
    const std::vector<IntBound>& b) {
  return canonicalize_bounds(a) == canonicalize_bounds(b);
}

void assert_bounds_equal(
    const std::vector<IntBound>& got,
    const std::vector<std::vector<int64_t>>& expected_coords) {
  auto got_c = canonicalize_bounds(got);
  auto exp_c = expected_coords;
  std::sort(exp_c.begin(), exp_c.end());
  exp_c.erase(std::unique(exp_c.begin(), exp_c.end()), exp_c.end());
  CHECK(got_c == exp_c);
}

void test_example_2_8_from_paper() {
  std::cout << "--- Example 2.8 (paper, GP) ---\n";

  const std::vector<int64_t> M = {10, 10, 10};
  const std::vector<int64_t> m = {0, 0, 0};

  NeighborhoodBoundSet<int64_t> nbs(M, m);

  // (1) N = {z1}, z1 = (4,0,4)
  nbs.update(IntPoint("z1", {4, 0, 4}));
  assert_bounds_equal(nbs.bounds(), {{4, 10, 10}, {10, 0, 10}, {10, 10, 4}});

  // (2) N = {z1,z2}, z2 = (3,3,1)
  nbs.update(IntPoint("z2", {3, 3, 1}));
  assert_bounds_equal(
      nbs.bounds(),
      {
          {3, 10, 10},
          {4, 3, 10},
          {10, 0, 10},
          {10, 3, 4},
          {10, 10, 1},
      });

  // (3) N = {z1,z2,z3}, z3 = (2,2,2)
  nbs.update(IntPoint("z3", {2, 2, 2}));
  assert_bounds_equal(
      nbs.bounds(),
      {
          {2, 10, 10},
          {3, 10, 2},
          {4, 2, 10},
          {10, 0, 10},
          {10, 2, 4},
          {10, 3, 2},
          {10, 10, 1},
      });

  std::cout << "Example 2.8 passed.\n\n";
}

void test_example_4_2_from_paper() {
  std::cout << "--- Example 4.2 (paper) ---\n";

  const std::vector<int64_t> M = {10, 10, 10};
  const std::vector<int64_t> m = {0, 0, 0};

  NeighborhoodBoundSet<int64_t> nbs(M, m);

  // Start from N={z1} then insert z2 exactly as Example 4.2.
  const IntPoint z1("z1", {4, 0, 4});
  const IntPoint z2("z2", {3, 3, 1});

  nbs.update(z1);
  nbs.update(z2);

  assert_bounds_equal(
      nbs.bounds(),
      {
          {3, 10, 10},
          {4, 3, 10},
          {10, 0, 10},
          {10, 3, 4},
          {10, 10, 1},
      });

  std::cout << "Example 4.2 passed.\n\n";
}

void test_example_4_13_ngp_from_paper() {
  std::cout << "--- Example 4.13 (paper, NGP) ---\n";

  const std::vector<int64_t> M = {10, 10, 10};
  const std::vector<int64_t> m = {0, 0, 0};

  NeighborhoodBoundSet<int64_t> nbs(M, m);

  // z1=(4,0,4), z2=(4,3,1), z3=(2,3,2)
  nbs.update(IntPoint("z1", {4, 0, 4}));
  nbs.update(IntPoint("z2", {4, 3, 1}));
  nbs.update(IntPoint("z3", {2, 3, 2}));

  // Quasi-upper bound set from paper includes 7 vectors.
  assert_bounds_equal(
      nbs.bounds(),
      {
          {2, 10, 10},
          {4, 3, 10},
          {10, 0, 10},
          {4, 3, 4},
          {4, 10, 2},
          {10, 3, 4},
          {10, 10, 1},
      });

  std::cout << "Example 4.13 passed.\n\n";
}

void test_zbar_dominates_several_bounds() {
  std::cout << "--- Multi-dominated update (z_bar dominates several bounds) ---\n";

  const std::vector<int64_t> M = {10, 10, 10};
  const std::vector<int64_t> m = {0, 0, 0};

  NeighborhoodBoundSet<int64_t> nbs(M, m);
  BoundSet<int64_t> reference(M, m);

  const std::vector<IntPoint> seed_points = {
      IntPoint("a", {8, 8, 8}),
      IntPoint("b", {7, 9, 9}),
      IntPoint("c", {9, 7, 9}),
      IntPoint("d", {9, 9, 7}),
  };

  for (const auto& p : seed_points) {
    nbs.update(p);
    reference.update_ra(p);
  }

  // z_bar that dominates multiple current bounds.
  const IntPoint z_bar("z_bar", {2, 2, 2});
  nbs.update(z_bar);
  reference.update_ra(z_bar);

  CHECK(compare_bounds(nbs.bounds(), reference.bounds()));
  std::cout << "Multi-dominated update passed.\n\n";
}

void test_fixed_cross_validation() {
  std::cout << "--- Fixed cross-validation (Algorithm 1 vs Algorithm 5) ---\n";

  const std::vector<int64_t> M = {100, 100, 100, 100};
  const std::vector<int64_t> m = {0, 0, 0, 0};

  NeighborhoodBoundSet<int64_t> nbs(M, m);
  BoundSet<int64_t> ra(M, m);

  const std::vector<IntPoint> points = {
      IntPoint("p1", {50, 50, 50, 50}),
      IntPoint("p2", {40, 60, 60, 60}),
      IntPoint("p3", {60, 40, 60, 60}),
      IntPoint("p4", {60, 60, 40, 60}),
      IntPoint("p5", {60, 60, 60, 40}),
      IntPoint("p6", {30, 70, 70, 70}),
      IntPoint("p7", {70, 30, 70, 70}),
      IntPoint("p8", {70, 70, 30, 70}),
      IntPoint("p9", {70, 70, 70, 30}),
      IntPoint("p10", {20, 20, 80, 80}),
      IntPoint("p11", {80, 80, 20, 20}),
      IntPoint("p12", {10, 90, 90, 90}),
      IntPoint("p13", {90, 10, 90, 90}),
      IntPoint("p14", {90, 90, 10, 90}),
      IntPoint("p15", {90, 90, 90, 10}),
  };

  for (const auto& p : points) {
    nbs.update(p);
    ra.update_ra(p);
    CHECK(compare_bounds(nbs.nonredundant_bounds(), ra.bounds()));
  }

  std::cout << "Fixed cross-validation passed.\n\n";
}

void test_update_return_value() {
  std::cout << "--- Test update() return value ---\n";

  const std::vector<int64_t> M = {10, 10, 10};
  const std::vector<int64_t> m = {0, 0, 0};

  NeighborhoodBoundSet<int64_t> nbs(M, m);

  // Valid point in search region -> returns true
  bool r1 = nbs.update(IntPoint("z1", {4, 0, 4}));
  CHECK(r1);
  (void)r1;

  // Dominated point (5, 5, 5) is outside search region -> returns false, bounds unchanged
  size_t size_before = nbs.size();
  bool r_dom = nbs.update(IntPoint("z_dom", {5, 5, 5}));
  CHECK(!r_dom);
  CHECK(nbs.size() == size_before);
  (void)size_before;
  (void)r_dom;

  std::cout << "update() return value tests passed.\n\n";
}

template <Objective Sense>
void test_tied_graph_export() {
  using NBS = NeighborhoodBoundSet<int64_t, Sense>;
  constexpr bool minimize = Sense == Objective::MINIMIZE;
  NBS nbs(std::vector<int64_t>(3, minimize ? 6 : 0),
          std::vector<int64_t>(3, minimize ? 0 : 6));
  const std::vector<std::vector<int64_t>> points = {
      {5, 2, 5}, {2, 4, 3}, {4, 4, 1}, {4, 3, 2}, {3, 3, 3}};
  for (std::size_t i = 0; i < points.size(); ++i) {
    auto coordinates = points[i];
    if (!minimize) for (auto& coordinate : coordinates) coordinate = 6 - coordinate;
    CHECK(nbs.update(IntPoint("z" + std::to_string(i), coordinates)));
    // Export between updates to catch classifications retained across mutations.
    CHECK(nbs.get_adjacency_graph().nodes.size() == nbs.nonredundant_size());
  }

  const auto raw = nbs.get_adjacency_graph(true);
  const std::vector<std::string> ids = {
      "u011", "u01213", "u02", "u0331", "u01211", "u01212",
      "u03321", "u0333", "u03323", "u0322", "u0122"};
  const std::vector<std::vector<int64_t>> coordinates = {
      {2, 6, 6}, {4, 4, 3}, {6, 2, 6}, {4, 6, 3}, {3, 4, 6}, {4, 3, 6},
      {4, 4, 3}, {6, 6, 1}, {6, 4, 2}, {6, 3, 5}, {5, 3, 6}};
  const std::vector<std::vector<std::string>> defining_ids = {
      {"z1", "z_hat2", "z_hat3"}, {"z3", "z1", "z4"},
      {"z_hat1", "z0", "z_hat3"}, {"z2", "z_hat2", "z1"},
      {"z4", "z1", "z_hat3"}, {"z3", "z4", "z_hat3"},
      {"z3", "z2", "z1"}, {"z_hat1", "z_hat2", "z2"},
      {"z_hat1", "z2", "z3"}, {"z_hat1", "z3", "z0"},
      {"z0", "z3", "z_hat3"}};
  CHECK(raw.nodes.size() == ids.size());
  CHECK(raw.quasi == std::vector<bool>({false, true, false, false, false, true,
                                       true, false, false, false, false}));
  CHECK(raw.adjacency_list == std::vector<std::vector<std::size_t>>({
      {4, 3}, {4, 5, 6}, {10, 9}, {0, 6, 7}, {0, 5, 1}, {4, 10, 1},
      {3, 1, 8}, {3, 8}, {6, 9, 7}, {10, 2, 8}, {5, 2, 9}}));
  const auto missing = NBS::npos;
  CHECK(raw.k_neighbors == std::vector<std::vector<std::size_t>>({
      {missing, 4, 3}, {4, 5, 6}, {10, missing, 9}, {0, 6, 7}, {0, 5, 1},
      {4, 10, 1}, {3, 1, 8}, {3, 8, missing}, {6, 9, 7}, {10, 2, 8}, {5, 2, 9}}));
  for (std::size_t i = 0; i < ids.size(); ++i) {
    CHECK(raw.nodes[i].id == ids[i]);
    auto expected = coordinates[i];
    if (!minimize) for (auto& coordinate : expected) coordinate = 6 - coordinate;
    CHECK(raw.nodes[i].coordinates == expected);
    CHECK(raw.nodes[i].defining_points.size() == 3);
    CHECK(raw.nodes[i].defining_point_sets.size() == 3);
    for (std::size_t k = 0; k < 3; ++k) {
      CHECK(raw.nodes[i].defining_points[k].id == defining_ids[i][k]);
      CHECK(raw.nodes[i].defining_point_sets[k].size() == 1);
      CHECK(raw.nodes[i].defining_point_sets[k][0].id == defining_ids[i][k]);
    }
  }

  const auto clean = nbs.get_adjacency_graph();
  const std::vector<std::size_t> retained = {0, 2, 3, 4, 7, 8, 9, 10};
  CHECK(clean.nodes.size() == retained.size());
  CHECK(clean.quasi == std::vector<bool>(retained.size(), false));
  for (std::size_t i = 0; i < retained.size(); ++i) CHECK(clean.nodes[i].id == ids[retained[i]]);
  // Complete contraction includes edges absent from the component-neighbor array.
  CHECK(clean.adjacency_list == std::vector<std::vector<std::size_t>>({
      {3, 2}, {7, 6}, {0, 4, 5, 3, 7}, {0, 7, 2, 5},
      {2, 5}, {6, 4, 2, 3, 7}, {7, 1, 5}, {1, 6, 3, 2, 5}}));
  CHECK(clean.k_neighbors == std::vector<std::vector<std::size_t>>({
      {missing, 3, 2}, {7, missing, 6}, {0, 7, 4}, {0, 7, 5},
      {2, 5, missing}, {2, 6, 4}, {7, 1, 5}, {3, 1, 6}}));
  CHECK(nbs.get_clean_neighbors(3) == std::vector<std::size_t>({0, 7, 8, 4, 11}));
  for (auto inactive : {std::size_t{10}, nbs.max_node_index()}) {
    bool rejected = false;
    try { (void)nbs.get_clean_neighbors(inactive); }
    catch (const std::out_of_range& error) { rejected = std::string(error.what()) == "Node must be active"; }
    CHECK(rejected);
  }
}

int main() {
  test_example_2_8_from_paper();
  test_example_4_2_from_paper();
  test_example_4_13_ngp_from_paper();
  test_zbar_dominates_several_bounds();
  test_fixed_cross_validation();
  test_update_return_value();
  test_tied_graph_export<Objective::MINIMIZE>();
  test_tied_graph_export<Objective::MAXIMIZE>();

  std::cout << "All NeighborhoodBoundSet tests passed successfully!\n";
  return 0;
}
