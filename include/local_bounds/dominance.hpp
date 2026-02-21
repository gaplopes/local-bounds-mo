/**
 * @file dominance.hpp
 * @brief Dominance relation functions for multiobjective optimization.
 *
 * Provides templated dominance check functions that work for both
 * minimization and maximization problems.
 *
 * Reference: "On the representation of the search region in multiobjective optimization"
 * http://dx.doi.org/10.1016/j.ejor.2015.03.031
 */

#ifndef LOCAL_BOUNDS_DOMINANCE_HPP
#define LOCAL_BOUNDS_DOMINANCE_HPP

#include <cassert>
#include <vector>

#include "types.hpp"

namespace local_bounds {

namespace detail {

/**
 * @brief Compares two values based on optimization sense.
 *
 * For MINIMIZE: returns true if a < b
 * For MAXIMIZE: returns true if a > b
 */
template <typename T, Objective Sense>
constexpr bool is_better(const T& a, const T& b) {
    if constexpr (Sense == Objective::MINIMIZE) {
        return a < b;
    } else {
        return a > b;
    }
}

/**
 * @brief Compares two values for weak improvement based on optimization sense.
 *
 * For MINIMIZE: returns true if a <= b
 * For MAXIMIZE: returns true if a >= b
 */
template <typename T, Objective Sense>
constexpr bool is_at_least_as_good(const T& a, const T& b) {
    if constexpr (Sense == Objective::MINIMIZE) {
        return a <= b;
    } else {
        return a >= b;
    }
}

}  // namespace detail

/**
 * @brief Checks if v1 weakly dominates v2.
 *
 * For MINIMIZE: v1 weakly dominates v2 if v1[i] <= v2[i] for all i.
 * For MAXIMIZE: v1 weakly dominates v2 if v1[i] >= v2[i] for all i.
 *
 * @tparam T Numeric type for coordinates.
 * @tparam Sense Optimization sense (MINIMIZE or MAXIMIZE).
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return true if v1 weakly dominates v2.
 */
template <typename T, Objective Sense = Objective::MINIMIZE>
bool weakly_dominates(const std::vector<T>& v1, const std::vector<T>& v2) {
    assert(v1.size() == v2.size() && "Vectors must have same dimensions");
    for (std::size_t i = 0; i < v1.size(); ++i) {
        if (!detail::is_at_least_as_good<T, Sense>(v1[i], v2[i])) {
            return false;
        }
    }
    return true;
}

/**
 * @brief Checks if v1 strictly dominates v2.
 *
 * For MINIMIZE: v1 strictly dominates v2 if v1[i] < v2[i] for all i.
 * For MAXIMIZE: v1 strictly dominates v2 if v1[i] > v2[i] for all i.
 *
 * @tparam T Numeric type for coordinates.
 * @tparam Sense Optimization sense (MINIMIZE or MAXIMIZE).
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return true if v1 strictly dominates v2.
 */
template <typename T, Objective Sense = Objective::MINIMIZE>
bool strictly_dominates(const std::vector<T>& v1, const std::vector<T>& v2) {
    assert(v1.size() == v2.size() && "Vectors must have same dimensions");
    for (std::size_t i = 0; i < v1.size(); ++i) {
        if (!detail::is_better<T, Sense>(v1[i], v2[i])) {
            return false;
        }
    }
    return true;
}

/**
 * @brief Checks if v1 dominates v2 (Pareto dominance).
 *
 * v1 dominates v2 if v1 weakly dominates v2 and v1 != v2.
 *
 * @tparam T Numeric type for coordinates.
 * @tparam Sense Optimization sense (MINIMIZE or MAXIMIZE).
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return true if v1 dominates v2.
 */
template <typename T, Objective Sense = Objective::MINIMIZE>
bool dominates(const std::vector<T>& v1, const std::vector<T>& v2) {
    if (v1 == v2) return false;
    return weakly_dominates<T, Sense>(v1, v2);
}

/**
 * @brief Checks if two points are incomparable (mutually nondominating).
 *
 * @tparam T Numeric type for coordinates.
 * @tparam Sense Optimization sense (MINIMIZE or MAXIMIZE).
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return true if neither v1 dominates v2 nor v2 dominates v1.
 */
template <typename T, Objective Sense = Objective::MINIMIZE>
bool incomparable(const std::vector<T>& v1, const std::vector<T>& v2) {
    return !dominates<T, Sense>(v1, v2) && !dominates<T, Sense>(v2, v1);
}

}  // namespace local_bounds

#endif  // LOCAL_BOUNDS_DOMINANCE_HPP
