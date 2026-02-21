/**
 * @file basic_usage.cpp
 * @brief Example demonstrating the Local Bounds library.
 *
 * Shows how to use the library for both minimization and maximization problems.
 */

#include <iostream>
#include <vector>

#include "local_bounds.hpp"

using namespace local_bounds;

void example_minimization() {
    std::cout << "=== Minimization Example (Bi-objective) ===\n\n";

    // Reference point (nadir) - upper bound of the search space
    std::vector<int64_t> nadir = {100, 100};
    
    // Create a bound set for minimization
    BoundSet<int64_t, Objective::MINIMIZE> bounds(nadir);
    std::cout << "Initial bounds (reference point): " << bounds.size() << "\n";

    // Simulate finding nondominated points during optimization
    std::vector<Point<int64_t>> solutions = {
        Point<int64_t>("z1", {30, 70}),
        Point<int64_t>("z2", {50, 50}),
        Point<int64_t>("z3", {70, 30})
    };

    for (const auto& sol : solutions) {
        // Check if point is in search region (could find new nondominated points)
        bool in_region = bounds.is_in_search_region(sol.coordinates);
        std::cout << "\nPoint " << sol.to_string() 
                  << " is in search region: " << (in_region ? "yes" : "no") << "\n";

        // Update bounds with new nondominated point
        bounds.update_re(sol);
        std::cout << "After update: " << bounds.size() << " bounds\n";
    }

    // Display final bounds
    std::cout << "\nFinal upper bound set:\n";
    for (const auto& b : bounds.bounds()) {
        std::cout << "  " << b.to_string() << "\n";
    }
    std::cout << "\n";
}

void example_maximization() {
    std::cout << "=== Maximization Example (Tri-objective) ===\n\n";

    // Reference point (ideal) - lower bound of the search space
    std::vector<double> ideal = {0.0, 0.0, 0.0};
    
    // Create a bound set for maximization
    BoundSet<double, Objective::MAXIMIZE> bounds(ideal);
    std::cout << "Initial bounds: " << bounds.size() << "\n";

    // Add some solutions
    std::vector<Point<double>> solutions = {
        Point<double>("z1", {0.8, 0.3, 0.5}),
        Point<double>("z2", {0.4, 0.7, 0.6}),
        Point<double>("z3", {0.5, 0.5, 0.9})
    };

    for (const auto& sol : solutions) {
        bounds.update_re(sol);
        std::cout << "After " << sol.id << ": " << bounds.size() << " bounds\n";
    }

    // Display final bounds
    std::cout << "\nFinal lower bound set:\n";
    for (const auto& b : bounds.bounds()) {
        std::cout << "  " << b.to_string() << "\n";
    }
    std::cout << "\n";
}

void example_dominance_checks() {
    std::cout << "=== Dominance Check Examples ===\n\n";

    std::vector<int> a = {1, 2, 3};
    std::vector<int> b = {2, 3, 4};
    std::vector<int> c = {1, 3, 2};

    std::cout << "a = (1, 2, 3), b = (2, 3, 4), c = (1, 3, 2)\n\n";

    // Minimization
    std::cout << "For MINIMIZE:\n";
    std::cout << "  a dominates b: " 
              << (dominates<int, Objective::MINIMIZE>(a, b) ? "yes" : "no") << "\n";
    std::cout << "  a strictly dominates b: " 
              << (strictly_dominates<int, Objective::MINIMIZE>(a, b) ? "yes" : "no") << "\n";
    std::cout << "  a and c incomparable: " 
              << (incomparable<int, Objective::MINIMIZE>(a, c) ? "yes" : "no") << "\n";

    // Maximization
    std::cout << "\nFor MAXIMIZE:\n";
    std::cout << "  b dominates a: " 
              << (dominates<int, Objective::MAXIMIZE>(b, a) ? "yes" : "no") << "\n";
    std::cout << "  b strictly dominates a: " 
              << (strictly_dominates<int, Objective::MAXIMIZE>(b, a) ? "yes" : "no") << "\n";
}

int main() {
    std::cout << "=============================================\n";
    std::cout << "  Local Bounds Library - Usage Examples\n";
    std::cout << "=============================================\n\n";

    example_minimization();
    example_maximization();
    example_dominance_checks();

    std::cout << "\n=============================================\n";
    return 0;
}
