/**
 * @file test_bound_set.cpp
 * @brief Tests for BoundSet class.
 *
 * Uses examples from the paper "On the representation of the search region
 * in multiobjective optimization" to validate the implementation.
 */

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "local_bounds.hpp"

using namespace local_bounds;

int tests_passed = 0;
int tests_failed = 0;

#define TEST(name, condition)                                        \
    do {                                                             \
        if (condition) {                                             \
            std::cout << "[PASS] " << name << std::endl;             \
            ++tests_passed;                                          \
        } else {                                                     \
            std::cout << "[FAIL] " << name << std::endl;             \
            ++tests_failed;                                          \
        }                                                            \
    } while (0)

void print_bounds(const std::vector<LocalBound<int64_t>>& bounds) {
    for (const auto& b : bounds) {
        std::cout << "  " << b.to_string() << std::endl;
    }
}

/**
 * Example 1 from the paper (p = 2, bi-objective).
 * 
 * For a stable set N = {z^1, z^2, z^3} where:
 * z^1 = (3, 7), z^2 = (5, 5), z^3 = (7, 3)
 * 
 * The upper bound set should be:
 * U(N) = {(3, M), (5, 7), (7, 5), (M, 3)}
 * 
 * With reference point M = (10, 10).
 */
void test_2d_example_minimization() {
    std::cout << "\n=== Example 1: 2D Minimization ===" << std::endl;

    const int64_t M = 10;
    BoundSet<int64_t, Objective::MINIMIZE> bound_set({M, M});

    // Initial state: single bound at reference point
    TEST("Initial bound set has 1 bound", bound_set.size() == 1);

    // Add z^1 = (3, 7)
    Point<int64_t> z1("z1", {3, 7});
    bound_set.update_re(z1);
    std::cout << "After adding z1 = (3, 7): " << bound_set.size() << " bounds" << std::endl;
    print_bounds(bound_set.bounds());
    TEST("After z1: 2 bounds", bound_set.size() == 2);

    // Add z^2 = (5, 5)
    Point<int64_t> z2("z2", {5, 5});
    bound_set.update_re(z2);
    std::cout << "After adding z2 = (5, 5): " << bound_set.size() << " bounds" << std::endl;
    print_bounds(bound_set.bounds());
    TEST("After z2: 3 bounds", bound_set.size() == 3);

    // Add z^3 = (7, 3)
    Point<int64_t> z3("z3", {7, 3});
    bound_set.update_re(z3);
    std::cout << "After adding z3 = (7, 3): " << bound_set.size() << " bounds" << std::endl;
    print_bounds(bound_set.bounds());
    
    // For 2D: |U(N)| = |N| + 1 = 4
    TEST("After z3: 4 bounds (|N| + 1 for 2D)", bound_set.size() == 4);
}

/**
 * Example 2 from the paper (p = 3, tri-objective).
 * 
 * For N = {z^1} where z^1 = (3, 5, 7):
 * U(N) = {(3, M, M), (M, 5, M), (M, M, 7)}
 * 
 * Adding z^2 = (6, 2, 4):
 * z^2 strictly dominates (M, 5, M) and (M, M, 7)
 */
void test_3d_example_minimization() {
    std::cout << "\n=== Example 2: 3D Minimization ===" << std::endl;

    const int64_t M = 100;
    BoundSet<int64_t, Objective::MINIMIZE> bound_set({M, M, M});

    TEST("Initial bound set has 1 bound", bound_set.size() == 1);

    // Add z^1 = (3, 5, 7)
    Point<int64_t> z1("z1", {3, 5, 7});
    bound_set.update_re(z1);
    std::cout << "After adding z1 = (3, 5, 7): " << bound_set.size() << " bounds" << std::endl;
    print_bounds(bound_set.bounds());
    TEST("After z1: 3 bounds", bound_set.size() == 3);

    // Add z^2 = (6, 2, 4)
    Point<int64_t> z2("z2", {6, 2, 4});
    bound_set.update_re(z2);
    std::cout << "After adding z2 = (6, 2, 4): " << bound_set.size() << " bounds" << std::endl;
    print_bounds(bound_set.bounds());
    
    // For 3D: |U(N)| <= 2|N| + 1 = 5
    TEST("After z2: 5 bounds", bound_set.size() == 5);
}

/**
 * Test maximization - same examples but inverted.
 */
void test_2d_example_maximization() {
    std::cout << "\n=== Example: 2D Maximization ===" << std::endl;

    const int64_t m = 0;  // Lower bound for maximization
    BoundSet<int64_t, Objective::MAXIMIZE> bound_set({m, m});

    TEST("Initial bound set has 1 bound", bound_set.size() == 1);

    // Add points (higher is better now)
    Point<int64_t> z1("z1", {7, 3});
    bound_set.update_re(z1);
    std::cout << "After adding z1 = (7, 3): " << bound_set.size() << " bounds" << std::endl;
    print_bounds(bound_set.bounds());
    TEST("After z1: 2 bounds", bound_set.size() == 2);

    Point<int64_t> z2("z2", {5, 5});
    bound_set.update_re(z2);
    std::cout << "After adding z2 = (5, 5): " << bound_set.size() << " bounds" << std::endl;
    print_bounds(bound_set.bounds());
    TEST("After z2: 3 bounds", bound_set.size() == 3);

    Point<int64_t> z3("z3", {3, 7});
    bound_set.update_re(z3);
    std::cout << "After adding z3 = (3, 7): " << bound_set.size() << " bounds" << std::endl;
    print_bounds(bound_set.bounds());
    TEST("After z3: 4 bounds", bound_set.size() == 4);
}

/**
 * Test search region membership.
 */
void test_search_region() {
    std::cout << "\n=== Testing Search Region Membership ===" << std::endl;

    const int64_t M = 10;
    BoundSet<int64_t, Objective::MINIMIZE> bound_set({M, M});

    // Add z^1 = (3, 7)
    Point<int64_t> z1("z1", {3, 7});
    bound_set.update_re(z1);

    // Point (2, 6) should be in search region (dominates some bound)
    TEST("(2, 6) is in search region", bound_set.is_in_search_region({2, 6}));

    // Point (4, 8) should NOT be in search region (dominated by z1)
    TEST("(4, 8) is NOT in search region", !bound_set.is_in_search_region({4, 8}));

    // Point (5, 5) should be in search region
    TEST("(5, 5) is in search region", bound_set.is_in_search_region({5, 5}));
}

/**
 * Test with double coordinates.
 */
void test_double_coordinates() {
    std::cout << "\n=== Testing Double Coordinates ===" << std::endl;

    BoundSet<double, Objective::MINIMIZE> bound_set({100.0, 100.0});

    Point<double> z1("z1", {3.5, 7.5});
    bound_set.update_re(z1);
    TEST("After z1: 2 bounds (double)", bound_set.size() == 2);

    Point<double> z2("z2", {5.5, 5.5});
    bound_set.update_re(z2);
    TEST("After z2: 3 bounds (double)", bound_set.size() == 3);
}

/**
 * Test the enhanced update algorithm.
 */
void test_enhanced_algorithm() {
    std::cout << "\n=== Testing Enhanced Algorithm ===" << std::endl;

    const int64_t M = 10;
    BoundSet<int64_t, Objective::MINIMIZE> bound_set({M, M});

    Point<int64_t> z1("z1", {3, 7});
    bound_set.update_re_enhanced(z1);
    TEST("Enhanced: After z1, 2 bounds", bound_set.size() == 2);

    Point<int64_t> z2("z2", {5, 5});
    bound_set.update_re_enhanced(z2);
    TEST("Enhanced: After z2, 3 bounds", bound_set.size() == 3);

    Point<int64_t> z3("z3", {7, 3});
    bound_set.update_re_enhanced(z3);
    TEST("Enhanced: After z3, 4 bounds", bound_set.size() == 4);
}

/**
 * Test Algorithm 5 (redundancy avoidance).
 */
void test_avoidance_algorithm() {
    std::cout << "\n=== Testing Algorithm 5 (Avoidance) ===" << std::endl;

    const int64_t M = 10;
    const int64_t m = 0;
    BoundSet<int64_t, Objective::MINIMIZE> bound_set({M, M}, {m, m});

    Point<int64_t> z1("z1", {3, 7});
    bound_set.update_ra(z1);
    TEST("Avoidance: After z1, 2 bounds", bound_set.size() == 2);

    Point<int64_t> z2("z2", {5, 5});
    bound_set.update_ra(z2);
    TEST("Avoidance: After z2, 3 bounds", bound_set.size() == 3);

    Point<int64_t> z3("z3", {7, 3});
    bound_set.update_ra(z3);
    TEST("Avoidance: After z3, 4 bounds", bound_set.size() == 4);

    std::cout << "Final bounds (Algorithm 5):" << std::endl;
    print_bounds(bound_set.bounds());
}

/**
 * Compare all algorithms produce the same number of bounds.
 */
void test_algorithm_consistency() {
    std::cout << "\n=== Testing Algorithm Consistency ===" << std::endl;

    const int64_t M = 100;
    const int64_t m = 0;
    std::vector<Point<int64_t>> points = {
        Point<int64_t>("z1", {3, 5, 7}),
        Point<int64_t>("z2", {6, 2, 4}),
        Point<int64_t>("z3", {4, 4, 5}),
        Point<int64_t>("z4", {8, 3, 2})
    };

    // Naive
    BoundSet<int64_t, Objective::MINIMIZE> bs_naive({M, M, M});
    for (const auto& p : points) bs_naive.update_naive(p);

    // Algorithm 2
    BoundSet<int64_t, Objective::MINIMIZE> bs_alg2({M, M, M});
    for (const auto& p : points) bs_alg2.update_re(p);

    // Algorithm 3
    BoundSet<int64_t, Objective::MINIMIZE> bs_alg3({M, M, M});
    for (const auto& p : points) bs_alg3.update_re_enhanced(p);

    // Algorithm 4 (SA)
    BoundSet<int64_t, Objective::MINIMIZE> bs_alg4({M, M, M}, {m, m, m});
    for (const auto& p : points) bs_alg4.update_ra_sa(p);

    // Algorithm 5
    BoundSet<int64_t, Objective::MINIMIZE> bs_alg5({M, M, M}, {m, m, m});
    for (const auto& p : points) bs_alg5.update_ra(p);

    std::cout << "Naive:       " << bs_naive.size() << " bounds" << std::endl;
    std::cout << "Algorithm 2: " << bs_alg2.size() << " bounds" << std::endl;
    std::cout << "Algorithm 3: " << bs_alg3.size() << " bounds" << std::endl;
    std::cout << "Algorithm 4: " << bs_alg4.size() << " bounds" << std::endl;
    std::cout << "Algorithm 5: " << bs_alg5.size() << " bounds" << std::endl;

    TEST("Naive and Algorithm 2 produce same number of bounds", 
         bs_naive.size() == bs_alg2.size());
    TEST("Algorithms 2 and 3 produce same number of bounds", 
         bs_alg2.size() == bs_alg3.size());
    TEST("Algorithms 2 and 4 produce same number of bounds", 
         bs_alg2.size() == bs_alg4.size());
    TEST("Algorithms 2 and 5 produce same number of bounds", 
         bs_alg2.size() == bs_alg5.size());
}

// --------------------------------------------------------------------------
// Helper: collect sorted bound coordinates for value-level comparison.
// --------------------------------------------------------------------------
std::vector<std::vector<int64_t>> sorted_coords(
    const std::vector<LocalBound<int64_t>>& bounds) {
    std::vector<std::vector<int64_t>> result;
    result.reserve(bounds.size());
    for (const auto& b : bounds) result.push_back(b.coordinates);
    std::sort(result.begin(), result.end());
    return result;
}

/**
 * Test the naive algorithm individually (was never tested).
 */
void test_naive_algorithm() {
    std::cout << "\n=== Testing Naive Algorithm ===" << std::endl;

    const int64_t M = 10;
    BoundSet<int64_t, Objective::MINIMIZE> bound_set({M, M});

    Point<int64_t> z1("z1", {3, 7});
    bound_set.update_naive(z1);
    TEST("Naive: After z1, 2 bounds", bound_set.size() == 2);

    Point<int64_t> z2("z2", {5, 5});
    bound_set.update_naive(z2);
    TEST("Naive: After z2, 3 bounds", bound_set.size() == 3);

    Point<int64_t> z3("z3", {7, 3});
    bound_set.update_naive(z3);
    TEST("Naive: After z3, 4 bounds", bound_set.size() == 4);
}

/**
 * Test Algorithm 4 (RA, General Position) individually.
 */
void test_avoidance_sa_algorithm() {
    std::cout << "\n=== Testing Algorithm 4 (RA-SA) ===" << std::endl;

    const int64_t M = 10;
    const int64_t m = 0;
    BoundSet<int64_t, Objective::MINIMIZE> bound_set({M, M}, {m, m});

    Point<int64_t> z1("z1", {3, 7});
    bound_set.update_ra_sa(z1);
    TEST("RA-SA: After z1, 2 bounds", bound_set.size() == 2);

    Point<int64_t> z2("z2", {5, 5});
    bound_set.update_ra_sa(z2);
    TEST("RA-SA: After z2, 3 bounds", bound_set.size() == 3);

    Point<int64_t> z3("z3", {7, 3});
    bound_set.update_ra_sa(z3);
    TEST("RA-SA: After z3, 4 bounds", bound_set.size() == 4);
}

/**
 * Test update_auto dispatch logic.
 */
void test_auto_algorithm() {
    std::cout << "\n=== Testing update_auto ===" << std::endl;

    // 2D (p < 6): should dispatch to RE Enhanced
    const int64_t M = 10;
    const int64_t m = 0;
    BoundSet<int64_t, Objective::MINIMIZE> bs_2d({M, M}, {m, m});
    bs_2d.update_auto(Point<int64_t>("z1", {3, 7}));
    bs_2d.update_auto(Point<int64_t>("z2", {7, 3}));
    TEST("Auto 2D: correct bound count", bs_2d.size() == 3);

    // 3D (p < 6): should dispatch to RE Enhanced
    BoundSet<int64_t, Objective::MINIMIZE> bs_3d({M, M, M}, {m, m, m});
    bs_3d.update_auto(Point<int64_t>("z1", {3, 5, 7}));
    bs_3d.update_auto(Point<int64_t>("z2", {6, 2, 4}));
    TEST("Auto 3D: correct bound count", bs_3d.size() == 5);

    // 6D (p >= 6): should dispatch to RA (Alg 5) when anti-reference provided
    std::vector<int64_t> ref6(6, 100);
    std::vector<int64_t> anti6(6, 0);
    BoundSet<int64_t, Objective::MINIMIZE> bs_6d(ref6, anti6);
    bs_6d.update_auto(Point<int64_t>("z1", {10, 20, 30, 40, 50, 60}));
    TEST("Auto 6D: bound set updated", bs_6d.size() == 6);

    // 6D without anti-reference: should dispatch to RE Enhanced
    BoundSet<int64_t, Objective::MINIMIZE> bs_6d_no_anti(ref6);
    bs_6d_no_anti.update_auto(Point<int64_t>("z1", {10, 20, 30, 40, 50, 60}));
    TEST("Auto 6D no anti-ref: bound set updated", bs_6d_no_anti.size() == 6);
}

/**
 * Test adding a dominated point (outside all search zones).
 * The bound set should remain unchanged.
 */
void test_dominated_point() {
    std::cout << "\n=== Testing Dominated Point (No-Op) ===" << std::endl;

    const int64_t M = 10;
    BoundSet<int64_t, Objective::MINIMIZE> bound_set({M, M});

    Point<int64_t> z1("z1", {3, 7});
    bound_set.update_re(z1);
    TEST("Before: 2 bounds", bound_set.size() == 2);

    // (4, 8) is dominated by z1 = (3,7) — outside all search zones
    Point<int64_t> z_dom("z_dom", {4, 8});
    bound_set.update_re(z_dom);
    TEST("After dominated point: still 2 bounds", bound_set.size() == 2);

    // Also test with naive
    BoundSet<int64_t, Objective::MINIMIZE> bs_naive({M, M});
    bs_naive.update_naive(Point<int64_t>("z1", {3, 7}));
    bs_naive.update_naive(Point<int64_t>("z_dom", {4, 8}));
    TEST("Naive: After dominated point, still 2 bounds", bs_naive.size() == 2);
}

/**
 * Verify actual bound coordinates, not just sizes.
 * Uses the 2D example where the exact bound set is known.
 */
void test_bound_coordinates() {
    std::cout << "\n=== Testing Bound Coordinates ===" << std::endl;

    const int64_t M = 10;
    BoundSet<int64_t, Objective::MINIMIZE> bound_set({M, M});

    bound_set.update_re(Point<int64_t>("z1", {3, 7}));
    bound_set.update_re(Point<int64_t>("z2", {5, 5}));
    bound_set.update_re(Point<int64_t>("z3", {7, 3}));

    // Expected: {(3,10), (5,7), (7,5), (10,3)}
    auto coords = sorted_coords(bound_set.bounds());
    std::vector<std::vector<int64_t>> expected = {
        {3, 10}, {5, 7}, {7, 5}, {10, 3}
    };
    std::sort(expected.begin(), expected.end());

    TEST("2D bound coordinates match expected", coords == expected);

    // Same with Naive
    BoundSet<int64_t, Objective::MINIMIZE> bs_naive({M, M});
    bs_naive.update_naive(Point<int64_t>("z1", {3, 7}));
    bs_naive.update_naive(Point<int64_t>("z2", {5, 5}));
    bs_naive.update_naive(Point<int64_t>("z3", {7, 3}));
    TEST("Naive: same coordinates", sorted_coords(bs_naive.bounds()) == expected);

    // Same with Algorithm 5
    BoundSet<int64_t, Objective::MINIMIZE> bs_alg5({M, M}, {0, 0});
    bs_alg5.update_ra(Point<int64_t>("z1", {3, 7}));
    bs_alg5.update_ra(Point<int64_t>("z2", {5, 5}));
    bs_alg5.update_ra(Point<int64_t>("z3", {7, 3}));
    TEST("Alg 5: same coordinates", sorted_coords(bs_alg5.bounds()) == expected);
}

/**
 * Test with non-general-position points (duplicate component values).
 * Algorithm 5 should handle these correctly via Z^k(u) sets.
 */
void test_non_general_position() {
    std::cout << "\n=== Testing Non-General-Position Points ===" << std::endl;

    const int64_t M = 100;
    const int64_t m = 0;

    // Points with shared component values (violates general position)
    std::vector<Point<int64_t>> points = {
        Point<int64_t>("z1", {3, 5, 7}),
        Point<int64_t>("z2", {3, 7, 4}),  // shares z1[0] = 3
        Point<int64_t>("z3", {6, 5, 2}),  // shares z1[1] = 5
    };

    // Algorithm 5 (handles general case)
    BoundSet<int64_t, Objective::MINIMIZE> bs_alg5({M, M, M}, {m, m, m});
    for (const auto& p : points) bs_alg5.update_ra(p);

    // Naive as ground truth
    BoundSet<int64_t, Objective::MINIMIZE> bs_naive({M, M, M});
    for (const auto& p : points) bs_naive.update_naive(p);

    // Algorithm 2 (RE)
    BoundSet<int64_t, Objective::MINIMIZE> bs_alg2({M, M, M});
    for (const auto& p : points) bs_alg2.update_re(p);

    std::cout << "Non-GP — Naive: " << bs_naive.size()
              << ", Alg 2: " << bs_alg2.size()
              << ", Alg 5: " << bs_alg5.size() << std::endl;

    TEST("Non-GP: Naive and Alg 2 agree", bs_naive.size() == bs_alg2.size());
    TEST("Non-GP: Naive and Alg 5 agree", bs_naive.size() == bs_alg5.size());
    TEST("Non-GP: Coordinate-level Naive vs Alg 5",
         sorted_coords(bs_naive.bounds()) == sorted_coords(bs_alg5.bounds()));
}

/**
 * Test maximization with all algorithms.
 */
void test_maximization_all_algorithms() {
    std::cout << "\n=== Testing Maximization (All Algorithms) ===" << std::endl;

    const int64_t m = 0;    // lower bound (reference for MAX)
    const int64_t M = 100;  // upper bound (anti-reference for MAX)

    std::vector<Point<int64_t>> points = {
        Point<int64_t>("z1", {7, 3}),
        Point<int64_t>("z2", {5, 5}),
        Point<int64_t>("z3", {3, 7}),
    };

    // Expected: {(7,0), (5,3), (3,5), (0,7)}
    std::vector<std::vector<int64_t>> expected = {
        {0, 7}, {3, 5}, {5, 3}, {7, 0}
    };
    std::sort(expected.begin(), expected.end());

    // Naive
    BoundSet<int64_t, Objective::MAXIMIZE> bs_naive({m, m});
    for (const auto& p : points) bs_naive.update_naive(p);
    TEST("MAX Naive: 4 bounds", bs_naive.size() == 4);
    TEST("MAX Naive: correct coordinates",
         sorted_coords(bs_naive.bounds()) == expected);

    // Algorithm 3 (Enhanced RE)
    BoundSet<int64_t, Objective::MAXIMIZE> bs_alg3({m, m});
    for (const auto& p : points) bs_alg3.update_re_enhanced(p);
    TEST("MAX Alg 3: 4 bounds", bs_alg3.size() == 4);
    TEST("MAX Alg 3: correct coordinates",
         sorted_coords(bs_alg3.bounds()) == expected);

    // Algorithm 4 (RA-SA)
    BoundSet<int64_t, Objective::MAXIMIZE> bs_alg4({m, m}, {M, M});
    for (const auto& p : points) bs_alg4.update_ra_sa(p);
    TEST("MAX Alg 4: 4 bounds", bs_alg4.size() == 4);
    TEST("MAX Alg 4: correct coordinates",
         sorted_coords(bs_alg4.bounds()) == expected);

    // Algorithm 5 (RA)
    BoundSet<int64_t, Objective::MAXIMIZE> bs_alg5({m, m}, {M, M});
    for (const auto& p : points) bs_alg5.update_ra(p);
    TEST("MAX Alg 5: 4 bounds", bs_alg5.size() == 4);
    TEST("MAX Alg 5: correct coordinates",
         sorted_coords(bs_alg5.bounds()) == expected);
}

/**
 * Test with higher-dimensional points (4D).
 */
void test_higher_dimensions() {
    std::cout << "\n=== Testing Higher Dimensions (4D) ===" << std::endl;

    const int64_t M = 100;
    const int64_t m = 0;

    std::vector<Point<int64_t>> points = {
        Point<int64_t>("z1", {10, 20, 30, 40}),
        Point<int64_t>("z2", {30, 10, 40, 20}),
        Point<int64_t>("z3", {20, 40, 10, 30}),
    };

    BoundSet<int64_t, Objective::MINIMIZE> bs_naive({M, M, M, M});
    for (const auto& p : points) bs_naive.update_naive(p);

    BoundSet<int64_t, Objective::MINIMIZE> bs_alg3({M, M, M, M});
    for (const auto& p : points) bs_alg3.update_re_enhanced(p);

    BoundSet<int64_t, Objective::MINIMIZE> bs_alg5({M, M, M, M}, {m, m, m, m});
    for (const auto& p : points) bs_alg5.update_ra(p);

    std::cout << "4D — Naive: " << bs_naive.size()
              << ", Alg 3: " << bs_alg3.size()
              << ", Alg 5: " << bs_alg5.size() << std::endl;

    TEST("4D: Naive and Alg 3 agree", bs_naive.size() == bs_alg3.size());
    TEST("4D: Naive and Alg 5 agree", bs_naive.size() == bs_alg5.size());
    TEST("4D: Coordinates Naive vs Alg 3",
         sorted_coords(bs_naive.bounds()) == sorted_coords(bs_alg3.bounds()));
    TEST("4D: Coordinates Naive vs Alg 5",
         sorted_coords(bs_naive.bounds()) == sorted_coords(bs_alg5.bounds()));
}

int main() {
    std::cout << "===========================================\n";
    std::cout << "  BoundSet Test Suite\n";
    std::cout << "===========================================\n";

    test_2d_example_minimization();
    test_3d_example_minimization();
    test_2d_example_maximization();
    test_search_region();
    test_double_coordinates();
    test_enhanced_algorithm();
    test_avoidance_algorithm();
    test_naive_algorithm();
    test_avoidance_sa_algorithm();
    test_auto_algorithm();
    test_dominated_point();
    test_bound_coordinates();
    test_non_general_position();
    test_maximization_all_algorithms();
    test_higher_dimensions();
    test_algorithm_consistency();

    std::cout << "\n===========================================\n";
    std::cout << "  RESULTS: " << tests_passed << " passed, " 
              << tests_failed << " failed\n";
    std::cout << "===========================================\n";

    return tests_failed > 0 ? 1 : 0;
}

