/**
 * @file types.hpp
 * @brief Core data structures for the Local Bounds library.
 *
 * This header defines the fundamental types used for representing points
 * and local bounds in multiobjective optimization problems.
 *
 * Reference: "On the representation of the search region in multiobjective optimization"
 * http://dx.doi.org/10.1016/j.ejor.2015.03.031
 */

#ifndef LOCAL_BOUNDS_TYPES_HPP
#define LOCAL_BOUNDS_TYPES_HPP

#include <algorithm>
#include <cassert>
#include <sstream>
#include <string>
#include <vector>

namespace local_bounds {

/**
 * @brief Optimization sense (minimization or maximization).
 */
enum class Objective { MINIMIZE,
                       MAXIMIZE };

/**
 * @brief Represents a point in the objective space.
 *
 * @tparam T The numeric type for coordinates (e.g., int64_t, double).
 */
template <typename T = double>
struct Point {
  std::string id;              ///< Unique identifier for the point
  std::vector<T> coordinates;  ///< Objective values

  Point() = default;

  Point(std::string id, std::vector<T> coordinates)
      : id(std::move(id)), coordinates(std::move(coordinates)) {}

  explicit Point(std::vector<T> coordinates)
      : id(""), coordinates(std::move(coordinates)) {}

  bool operator==(const Point& other) const {
    return coordinates == other.coordinates;
  }

  bool operator!=(const Point& other) const {
    return !(*this == other);
  }

  /**
   * @brief Returns the dimensionality of the point.
   */
  [[nodiscard]] std::size_t dimensions() const {
    return coordinates.size();
  }

  /**
   * @brief Returns a string representation of the point.
   */
  [[nodiscard]] std::string to_string() const {
    std::ostringstream oss;
    oss << id << " (";
    for (std::size_t i = 0; i < coordinates.size(); ++i) {
      oss << coordinates[i];
      if (i != coordinates.size() - 1) oss << ", ";
    }
    oss << ")";
    return oss.str();
  }
};

/**
 * @brief Represents a local bound in the objective space.
 *
 * A local bound defines a search zone in multiobjective optimization.
 * For minimization, these are "upper bounds" (points to be dominated).
 * For maximization, these are "lower bounds" (points that dominate).
 *
 * @tparam T The numeric type for coordinates (e.g., int64_t, double).
 */
template <typename T = double>
struct LocalBound {
  std::string id;              ///< Unique identifier
  std::vector<T> coordinates;  ///< Bound coordinates

  /// Defining points z^j(u) for each component j (Algorithm 4, SA case).
  /// Size = p. Entry j is the single point z in N_hat such that z_j = u_j and z_{-j} < u_{-j}.
  std::vector<Point<T>> defining_points;

  /// Defining point sets Z^j(u) for each component j (Algorithm 5, general case).
  /// Size = p. Entry j is the set {z in N_hat : z_j = u_j and z_{-j} < u_{-j}}.
  std::vector<std::vector<Point<T>>> defining_point_sets;

  LocalBound() = default;

  LocalBound(std::string id, std::vector<T> coordinates,
             std::vector<Point<T>> defining_points = {},
             std::vector<std::vector<Point<T>>> defining_point_sets = {})
      : id(std::move(id)),
        coordinates(std::move(coordinates)),
        defining_points(std::move(defining_points)),
        defining_point_sets(std::move(defining_point_sets)) {}

  bool operator==(const LocalBound& other) const {
    return coordinates == other.coordinates;
  }

  bool operator!=(const LocalBound& other) const {
    return !(*this == other);
  }

  /**
   * @brief Returns the dimensionality of the bound.
   */
  [[nodiscard]] std::size_t dimensions() const {
    return coordinates.size();
  }

  /**
   * @brief Creates dummy points ẑ^j = (M_j, m_{-j}) for j = 1..p.
   *
   * These are used to ensure that every component of every local upper bound
   * is defined by some point in N_hat = N ∪ {ẑ^1, ..., ẑ^p}.
   *
   * @param reference_point The reference point M (nadir for min, ideal for max).
   * @param anti_reference  The anti-reference point m (ideal for min, nadir for max).
   * @return Vector of p dummy points.
   */
  static std::vector<Point<T>> make_dummy_points(
      const std::vector<T>& reference_point,
      const std::vector<T>& anti_reference) {
    const std::size_t p = reference_point.size();
    std::vector<Point<T>> dummies;
    dummies.reserve(p);
    for (std::size_t j = 0; j < p; ++j) {
      std::vector<T> coords = anti_reference;  // m_{-j}
      coords[j] = reference_point[j];          // M_j
      dummies.emplace_back("z_hat" + std::to_string(j + 1), std::move(coords));
    }
    return dummies;
  }

  /**
   * @brief Creates the initial bound U(∅) = {M}.
   *
   * For Algorithm 4 (SA): defining_points[j] = ẑ^j (the dummy point for component j).
   * For Algorithm 5 (general): defining_point_sets[j] = {ẑ^j}.
   *
   * @param reference_point The reference point M (nadir for min, ideal for max).
   * @param anti_reference  The anti-reference point m (ideal for min, nadir for max).
   * @param id Optional identifier for the bound.
   * @return The initial local bound with properly initialized defining sets.
   */
  static LocalBound initial(const std::vector<T>& reference_point,
                            const std::vector<T>& anti_reference,
                            const std::string& id = "u0") {
    const std::size_t p = reference_point.size();
    auto dummies = make_dummy_points(reference_point, anti_reference);

    // SA: z^j(M) = ẑ^j
    std::vector<Point<T>> def_pts(dummies.begin(), dummies.end());

    // General: Z^j(M) = {ẑ^j}
    std::vector<std::vector<Point<T>>> def_sets(p);
    for (std::size_t j = 0; j < p; ++j) {
      def_sets[j] = {dummies[j]};
    }

    return LocalBound(id, reference_point, std::move(def_pts), std::move(def_sets));
  }

  /**
   * @brief Convenience overload: creates the initial bound without dummy-point tracking.
   *
   * Used by Algorithms 2/3 (Naive, RE, Enhanced RE) which don't need defining sets.
   */
  static LocalBound initial(const std::vector<T>& reference_point,
                            const std::string& id = "u0") {
    return LocalBound(id, reference_point);
  }

  /**
   * @brief Returns a string representation of the local bound.
   */
  [[nodiscard]] std::string to_string() const {
    std::ostringstream oss;
    oss << id << " (";
    for (std::size_t i = 0; i < coordinates.size(); ++i) {
      oss << coordinates[i];
      if (i != coordinates.size() - 1) oss << ", ";
    }
    oss << ")";
    if (!defining_points.empty()) {
      oss << " def_pts={";
      for (std::size_t i = 0; i < defining_points.size(); ++i) {
        oss << defining_points[i].id;
        if (i != defining_points.size() - 1) oss << ", ";
      }
      oss << "}";
    }
    if (!defining_point_sets.empty()) {
      oss << " def_sets={";
      for (std::size_t i = 0; i < defining_point_sets.size(); ++i) {
        oss << "[";
        for (std::size_t j = 0; j < defining_point_sets[i].size(); ++j) {
          oss << defining_point_sets[i][j].id;
          if (j != defining_point_sets[i].size() - 1) oss << ", ";
        }
        oss << "]";
        if (i != defining_point_sets.size() - 1) oss << ", ";
      }
      oss << "}";
    }
    return oss.str();
  }
};

}  // namespace local_bounds

#endif  // LOCAL_BOUNDS_TYPES_HPP
