#include <algorithm>
#include <cassert>
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
  assert(got_c == exp_c);
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

  assert(compare_bounds(nbs.bounds(), reference.bounds()));
  std::cout << "Multi-dominated update passed.\n\n";
}

void test_random_cross_validation() {
  std::cout << "--- Random cross-validation (Algorithm 1 vs Algorithm 5) ---\n";

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
    assert(compare_bounds(nbs.bounds(), ra.bounds()));
  }

  std::cout << "Random cross-validation passed.\n\n";
}

int main() {
  test_example_2_8_from_paper();
  test_example_4_2_from_paper();
  test_example_4_13_ngp_from_paper();
  test_zbar_dominates_several_bounds();
  test_random_cross_validation();

  std::cout << "All NeighborhoodBoundSet tests passed successfully!\n";
  return 0;
}
