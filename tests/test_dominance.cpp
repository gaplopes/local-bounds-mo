/**
 * @file test_dominance.cpp
 * @brief Tests for dominance relation functions.
 */

#include <iostream>
#include <string>
#include <vector>

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

// Helper functions to avoid macro parsing issues with template commas
template <typename T>
bool weakly_dominates_min(const std::vector<T>& v1, const std::vector<T>& v2) {
    return weakly_dominates<T, Objective::MINIMIZE>(v1, v2);
}

template <typename T>
bool weakly_dominates_max(const std::vector<T>& v1, const std::vector<T>& v2) {
    return weakly_dominates<T, Objective::MAXIMIZE>(v1, v2);
}

template <typename T>
bool strictly_dominates_min(const std::vector<T>& v1, const std::vector<T>& v2) {
    return strictly_dominates<T, Objective::MINIMIZE>(v1, v2);
}

template <typename T>
bool strictly_dominates_max(const std::vector<T>& v1, const std::vector<T>& v2) {
    return strictly_dominates<T, Objective::MAXIMIZE>(v1, v2);
}

template <typename T>
bool dominates_min(const std::vector<T>& v1, const std::vector<T>& v2) {
    return dominates<T, Objective::MINIMIZE>(v1, v2);
}

template <typename T>
bool dominates_max(const std::vector<T>& v1, const std::vector<T>& v2) {
    return dominates<T, Objective::MAXIMIZE>(v1, v2);
}

template <typename T>
bool incomparable_min(const std::vector<T>& v1, const std::vector<T>& v2) {
    return incomparable<T, Objective::MINIMIZE>(v1, v2);
}

template <typename T>
bool incomparable_max(const std::vector<T>& v1, const std::vector<T>& v2) {
    return incomparable<T, Objective::MAXIMIZE>(v1, v2);
}

void test_weakly_dominates_minimize() {
    std::cout << "\n=== Testing weakly_dominates (MINIMIZE) ===" << std::endl;

    std::vector<int> v1 = {1, 2, 3};
    std::vector<int> v2 = {2, 3, 4};
    std::vector<int> v3 = {1, 2, 3};
    std::vector<int> v4 = {1, 3, 2};

    TEST("v1 weakly dominates v2 (all smaller)", weakly_dominates_min(v1, v2));
    TEST("v2 does NOT weakly dominate v1", !weakly_dominates_min(v2, v1));
    TEST("v1 weakly dominates v3 (equal)", weakly_dominates_min(v1, v3));
    TEST("v1 does NOT weakly dominate v4 (v1[2]=3 > v4[2]=2)", !weakly_dominates_min(v1, v4));
    TEST("v4 does NOT weakly dominate v1 (v4[1]=3 > v1[1]=2)", !weakly_dominates_min(v4, v1));
}

void test_weakly_dominates_maximize() {
    std::cout << "\n=== Testing weakly_dominates (MAXIMIZE) ===" << std::endl;

    std::vector<int> v1 = {4, 5, 6};
    std::vector<int> v2 = {3, 4, 5};
    std::vector<int> v3 = {4, 5, 6};

    TEST("v1 weakly dominates v2 (all larger)", weakly_dominates_max(v1, v2));
    TEST("v2 does NOT weakly dominate v1", !weakly_dominates_max(v2, v1));
    TEST("v1 weakly dominates v3 (equal)", weakly_dominates_max(v1, v3));
}

void test_strictly_dominates_minimize() {
    std::cout << "\n=== Testing strictly_dominates (MINIMIZE) ===" << std::endl;

    std::vector<int> v1 = {1, 2, 3};
    std::vector<int> v2 = {2, 3, 4};
    std::vector<int> v3 = {1, 2, 3};
    std::vector<int> v4 = {2, 2, 4};

    TEST("v1 strictly dominates v2", strictly_dominates_min(v1, v2));
    TEST("v1 does NOT strictly dominate v3 (equal)", !strictly_dominates_min(v1, v3));
    TEST("v1 does NOT strictly dominate v4 (2 == 2)", !strictly_dominates_min(v1, v4));
}

void test_strictly_dominates_maximize() {
    std::cout << "\n=== Testing strictly_dominates (MAXIMIZE) ===" << std::endl;

    std::vector<int> v1 = {5, 6, 7};
    std::vector<int> v2 = {4, 5, 6};
    std::vector<int> v3 = {5, 6, 7};

    TEST("v1 strictly dominates v2", strictly_dominates_max(v1, v2));
    TEST("v1 does NOT strictly dominate v3 (equal)", !strictly_dominates_max(v1, v3));
}

void test_dominates_minimize() {
    std::cout << "\n=== Testing dominates (MINIMIZE) ===" << std::endl;

    std::vector<int> v1 = {1, 2, 3};
    std::vector<int> v2 = {2, 3, 4};
    std::vector<int> v3 = {1, 2, 3};
    std::vector<int> v4 = {1, 2, 4};

    TEST("v1 dominates v2", dominates_min(v1, v2));
    TEST("v1 does NOT dominate v3 (equal vectors)", !dominates_min(v1, v3));
    TEST("v1 dominates v4 (1<=1, 2<=2, 3<4)", dominates_min(v1, v4));
}

void test_dominates_maximize() {
    std::cout << "\n=== Testing dominates (MAXIMIZE) ===" << std::endl;

    std::vector<int> v1 = {5, 6, 7};
    std::vector<int> v2 = {4, 5, 6};
    std::vector<int> v3 = {5, 6, 7};
    std::vector<int> v4 = {5, 6, 6};

    TEST("v1 dominates v2", dominates_max(v1, v2));
    TEST("v1 does NOT dominate v3 (equal)", !dominates_max(v1, v3));
    TEST("v1 dominates v4", dominates_max(v1, v4));
}

void test_incomparable() {
    std::cout << "\n=== Testing incomparable ===" << std::endl;

    std::vector<int> v1 = {1, 3};
    std::vector<int> v2 = {2, 2};

    TEST("v1 and v2 are incomparable (MINIMIZE)", incomparable_min(v1, v2));
    TEST("v1 and v2 are incomparable (MAXIMIZE)", incomparable_max(v1, v2));

    std::vector<int> v3 = {1, 1};
    std::vector<int> v4 = {2, 2};

    TEST("v3 and v4 are NOT incomparable (v3 dominates v4 in MIN)", !incomparable_min(v3, v4));
}

void test_self_dominance() {
    std::cout << "\n=== Testing Self-Dominance (Reflexivity) ===" << std::endl;

    std::vector<int> v = {3, 5, 7};

    TEST("v weakly dominates itself (MIN)", weakly_dominates_min(v, v));
    TEST("v weakly dominates itself (MAX)", weakly_dominates_max(v, v));
    TEST("v does NOT strictly dominate itself (MIN)", !strictly_dominates_min(v, v));
    TEST("v does NOT strictly dominate itself (MAX)", !strictly_dominates_max(v, v));
    TEST("v does NOT dominate itself (MIN)", !dominates_min(v, v));
    TEST("v does NOT dominate itself (MAX)", !dominates_max(v, v));
    TEST("v is incomparable with itself (MIN)", incomparable_min(v, v));
    TEST("v is incomparable with itself (MAX)", incomparable_max(v, v));
}

void test_double_coordinates() {
    std::cout << "\n=== Testing with double coordinates ===" << std::endl;

    std::vector<double> v1 = {1.5, 2.5, 3.5};
    std::vector<double> v2 = {1.6, 2.6, 3.6};

    TEST("v1 strictly dominates v2 (double, MINIMIZE)", strictly_dominates_min(v1, v2));
    TEST("v2 strictly dominates v1 (double, MAXIMIZE)", strictly_dominates_max(v2, v1));
}

int main() {
    std::cout << "===========================================\n";
    std::cout << "  Dominance Functions Test Suite\n";
    std::cout << "===========================================\n";

    test_weakly_dominates_minimize();
    test_weakly_dominates_maximize();
    test_strictly_dominates_minimize();
    test_strictly_dominates_maximize();
    test_dominates_minimize();
    test_dominates_maximize();
    test_incomparable();
    test_self_dominance();
    test_double_coordinates();

    std::cout << "\n===========================================\n";
    std::cout << "  RESULTS: " << tests_passed << " passed, " 
              << tests_failed << " failed\n";
    std::cout << "===========================================\n";

    return tests_failed > 0 ? 1 : 0;
}
