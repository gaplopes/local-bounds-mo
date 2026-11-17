#include <algorithm>
#include <stdexcept>
#define CHECK(condition) do { if (!(condition)) throw std::runtime_error(#condition); } while (false)
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <set>
#include <string>
#include <type_traits>
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

  // Valid point in search region returns true and the removed initial bound.
  auto r1 = nbs.update(IntPoint("z1", {4, 0, 4}));
  CHECK(r1.first);
  CHECK(r1.second.size() == 1 && r1.second[0].id == "u0");
  (void)r1;

  // Dominated point (5, 5, 5) returns {false, {}}, leaving bounds unchanged.
  size_t size_before = nbs.size();
  auto r_dom = nbs.update(IntPoint("z_dom", {5, 5, 5}));
  CHECK(!r_dom.first && r_dom.second.empty());
  CHECK(nbs.size() == size_before);
  (void)size_before;
  (void)r_dom;

  std::cout << "update() return value tests passed.\n\n";
}

// Paper 2, Example 4.13: with ties, the defining points follow the insertion
// order of eq. (6). z2 defines the first component of the quasi-nonredundant
// bound (4,3,4) although z2 ties with it in the second component.
void test_defining_points_match_paper_example_4_13() {
  std::cout << "--- Defining points of Example 4.13 (paper, NGP) ---\n";

  const int64_t M = 10;
  NeighborhoodBoundSet<int64_t> nbs({M, M, M}, {0, 0, 0});
  nbs.update(IntPoint("z1", {4, 0, 4}));
  nbs.update(IntPoint("z2", {4, 3, 1}));
  nbs.update(IntPoint("z3", {2, 3, 2}));

  const std::vector<std::pair<std::vector<int64_t>, std::vector<std::string>>> table = {
      {{2, M, M}, {"z3", "z_hat2", "z_hat3"}},
      {{4, 3, M}, {"z1", "z3", "z_hat3"}},
      {{M, 0, M}, {"z_hat1", "z1", "z_hat3"}},
      {{4, 3, 4}, {"z2", "z3", "z1"}},
      {{4, M, 2}, {"z2", "z_hat2", "z3"}},
      {{M, 3, 4}, {"z_hat1", "z2", "z1"}},
      {{M, M, 1}, {"z_hat1", "z_hat2", "z2"}},
  };
  const auto bounds = nbs.bounds();
  CHECK(bounds.size() == table.size());
  for (const auto& row : table) {
    const auto u = std::find_if(bounds.begin(), bounds.end(),
                                [&](const IntBound& b) { return b.coordinates == row.first; });
    CHECK(u != bounds.end());
    std::vector<std::string> got;
    for (const auto& z : u->defining_points) got.push_back(z.id);
    CHECK(got == row.second);
    CHECK(u->defining_point_sets.empty());  // Algorithm 1 does not maintain the sets Z^j(u)
  }

  std::cout << "Example 4.13 defining points passed.\n\n";
}

// Returns how many quasi-nonredundant bounds the run produced.
template <Objective Sense>
std::size_t check_exports_match_pairwise_scan(std::size_t p, int64_t K, unsigned seed) {
  constexpr bool minimize = Sense == Objective::MINIMIZE;
  const std::vector<int64_t> ref(p, minimize ? K + 1 : 0), anti(p, minimize ? 0 : K + 1);
  NeighborhoodBoundSet<int64_t, Sense> nbs(ref, anti);
  BoundSet<int64_t, Sense> re(ref);
  std::mt19937 rng(seed);
  // Includes the anti-reference boundary, which updates accept.
  std::uniform_int_distribution<int64_t> value(minimize ? 0 : 1, minimize ? K : K + 1);
  std::vector<std::vector<int64_t>> points;
  std::size_t quasi_seen = 0;

  for (int attempt = 0; attempt < 60; ++attempt) {
    std::vector<int64_t> z(p);
    for (auto& c : z) c = value(rng);
    // Keep the input stable: skip points comparable with an earlier one.
    if (std::any_of(points.begin(), points.end(), [&](const auto& q) {
          return weakly_dominates<int64_t, Sense>(q, z) || weakly_dominates<int64_t, Sense>(z, q);
        }))
      continue;
    points.push_back(z);
    nbs.update(IntPoint(z));
    re.update_re(IntPoint(z));

    // Independent classification: compare every active bound with every other,
    // keeping the first slot among equal coordinates. bounds() is in slot order.
    const auto raw = nbs.bounds();
    std::vector<bool> quasi(raw.size(), false);
    std::vector<std::vector<int64_t>> expected;
    for (std::size_t i = 0; i < raw.size(); ++i) {
      for (std::size_t j = 0; j < raw.size(); ++j)
        if (i != j && weakly_dominates<int64_t, Sense>(raw[i].coordinates, raw[j].coordinates) &&
            (raw[i].coordinates != raw[j].coordinates || j < i))
          quasi[i] = true;
      if (!quasi[i]) expected.push_back(raw[i].coordinates);
    }
    quasi_seen += raw.size() - expected.size();

    std::vector<std::vector<int64_t>> got;
    for (const auto& b : nbs.nonredundant_bounds()) got.push_back(b.coordinates);
    CHECK(got == expected);
    CHECK(nbs.nonredundant_size() == expected.size());
    CHECK(nbs.get_adjacency_graph(true).quasi == quasi);

    // The contracted graph must agree with the standalone per-node query.
    const auto graph = nbs.get_adjacency_graph();
    CHECK(graph.nodes.size() == expected.size());
    std::vector<std::size_t> internal, dense(nbs.max_node_index(), nbs.npos);
    for (std::size_t slot = 0, active = 0; slot < nbs.max_node_index(); ++slot) {
      if (!nbs.is_active(slot)) continue;
      if (!quasi[active++]) { dense[slot] = internal.size(); internal.push_back(slot); }
    }
    for (std::size_t u = 0; u < internal.size(); ++u) {
      std::vector<std::size_t> want;
      for (std::size_t slot : nbs.get_clean_neighbors(internal[u])) want.push_back(dense[slot]);
      CHECK(graph.adjacency_list[u] == want);
    }

    // And the filtered set is exactly U(N), as computed by Algorithm 2.
    std::sort(expected.begin(), expected.end());
    CHECK(expected == canonicalize_bounds(re.bounds()));
    CHECK(expected.size() == re.size());
  }
  return quasi_seen;
}

void test_nonredundant_exports_match_pairwise_scan() {
  std::cout << "--- Nonredundant exports vs pairwise scan (random ties) ---\n";

  std::size_t quasi_seen = 0;
  for (std::size_t p : {2, 3, 4})
    for (int64_t K : {4, 8, 30})
      for (unsigned seed = 0; seed < 15; ++seed) {
        quasi_seen += check_exports_match_pairwise_scan<Objective::MINIMIZE>(p, K, seed);
        quasi_seen += check_exports_match_pairwise_scan<Objective::MAXIMIZE>(p, K, seed);
      }
  CHECK(quasi_seen > 0);  // otherwise the inputs never exercised the filter

  std::cout << "Export classification passed (" << quasi_seen << " quasi bounds seen).\n\n";
}

int main() {
  test_defining_points_match_paper_example_4_13();
  test_nonredundant_exports_match_pairwise_scan();
  test_example_2_8_from_paper();
  test_example_4_2_from_paper();
  test_example_4_13_ngp_from_paper();
  test_zbar_dominates_several_bounds();
  test_fixed_cross_validation();
  test_update_return_value();

  std::cout << "All NeighborhoodBoundSet tests passed successfully!\n";
  return 0;
}
