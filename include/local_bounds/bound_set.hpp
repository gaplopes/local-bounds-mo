/**
 * @file bound_set.hpp
 * @brief Core class for maintaining local bounds in multiobjective optimization.
 *
 * Implements the algorithms from the paper for updating the bound set
 * when new nondominated points are discovered.
 *
 * Reference: "On the representation of the search region in multiobjective optimization"
 * http://dx.doi.org/10.1016/j.ejor.2015.03.031
 */

#ifndef LOCAL_BOUNDS_BOUND_SET_HPP
#define LOCAL_BOUNDS_BOUND_SET_HPP

#include <algorithm>
#include <cassert>
#include <iterator>
#include <limits>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "dominance.hpp"
#include "types.hpp"

namespace local_bounds {

/**
 * @brief Maintains a set of local bounds for multiobjective optimization.
 *
 * This class implements algorithms for updating local bounds when new
 * nondominated points are found during optimization. It supports both
 * minimization and maximization problems.
 *
 * For MINIMIZE: bounds represent "upper bounds" - the search region consists
 *               of points that are strictly dominated by at least one bound.
 * For MAXIMIZE: bounds represent "lower bounds" - the search region consists
 *               of points that strictly dominate at least one bound.
 *
 * @tparam T The numeric type for coordinates (e.g., int64_t, double).
 * @tparam Sense Objective sense (MINIMIZE or MAXIMIZE).
 */
template <typename T = double, Objective Sense = Objective::MINIMIZE>
class BoundSet {
 public:
  /**
   * @brief Constructs a BoundSet with the given reference point.
   *
   * @param reference_point For MINIMIZE: the nadir point (upper bound of search space).
   *                        For MAXIMIZE: the ideal point (lower bound of search space).
   */
  /**
   * @brief Constructs a BoundSet with only the reference point.
   *
   * Used by Algorithms 2/3/Naive which don't need defining-point tracking.
   * The anti-reference defaults to 0 for MINIMIZE, max for MAXIMIZE.
   *
   * @param reference_point For MINIMIZE: the nadir point M.
   *                        For MAXIMIZE: the ideal point m.
   */
  explicit BoundSet(const std::vector<T>& reference_point)
      : dimensions_(reference_point.size()),
        reference_point_(reference_point),
        next_bound_id_(1) {
    if (reference_point.empty()) {
      throw std::invalid_argument("Reference point must have at least one dimension");
    }
    bounds_.push_back(LocalBound<T>::initial(reference_point, "u0"));
  }

  /**
   * @brief Constructs a BoundSet with both reference and anti-reference points.
   *
   * Required by Algorithms 4/5 (Redundancy Avoidance) which need dummy points
   * ẑ^j = (M_j, m_{-j}) to track defining sets.
   *
   * @param reference_point For MINIMIZE: the nadir point M.
   *                        For MAXIMIZE: the ideal point m.
   * @param anti_reference  For MINIMIZE: the ideal point m (lower bound of search space).
   *                        For MAXIMIZE: the nadir point M (upper bound of search space).
   */
  BoundSet(const std::vector<T>& reference_point,
           const std::vector<T>& anti_reference)
      : dimensions_(reference_point.size()),
        reference_point_(reference_point),
        anti_reference_(anti_reference),
        next_bound_id_(1) {
    if (reference_point.empty()) {
      throw std::invalid_argument("Reference point must have at least one dimension");
    }
    if (anti_reference.size() != reference_point.size()) {
      throw std::invalid_argument("Anti-reference must have same dimensions as reference");
    }
    bounds_.push_back(LocalBound<T>::initial(reference_point, anti_reference, "u0"));
  }

  /**
   * @brief Updates the bound set with a new nondominated point.
   *
   * Implements Algorithm 2 (Redundancy Elimination) from the paper.
   *
   * @param point The new nondominated point.
   */
  void update_re(const Point<T>& point) {
    if (point.dimensions() != dimensions_) {
      throw std::invalid_argument("Point dimensions must match bound set dimensions");
    }
    update_redundancy_elimination(bounds_, point);
  }

  /**
   * @brief Updates the bound set using the enhanced filtering algorithm.
   *
   * Implements Algorithm 3 (Enhanced Redundancy Elimination) from the paper.
   *
   * @param point The new nondominated point.
   */
  void update_re_enhanced(const Point<T>& point) {
    if (point.dimensions() != dimensions_) {
      throw std::invalid_argument("Point dimensions must match bound set dimensions");
    }
    update_redundancy_elimination_enhanced(bounds_, point);
  }

  /**
   * @brief Updates the bound set using the redundancy avoidance algorithm.
   *
   * Implements Algorithm 4 (Redundancy Avoidance - General Position) from the paper.
   * This algorithm avoids the filtering step by leveraging the structure of bounds
   * under the general position assumption (SA). It maintains additional information to
   * directly generate only non-redundant bounds.
   *
   * @param point The new nondominated point.
   */
  void update_ra_sa(const Point<T>& point) {
    if (point.dimensions() != dimensions_) {
      throw std::invalid_argument("Point dimensions must match bound set dimensions");
    }
    update_redundancy_avoidance_sa(bounds_, point);
  }

  /**
   * @brief Updates the bound set using redundancy avoidance.
   *
   * Implements Algorithm 5 (Redundancy Avoidance - General Case) from the paper.
   * This algorithm theoretically avoids the filtering step entirely by tracking
   * the sets of points Z^j(u) that define each component of each local bound.
   *
   * @param point The new nondominated point.
   */
  void update_ra(const Point<T>& point) {
    if (point.dimensions() != dimensions_) {
      throw std::invalid_argument("Point dimensions must match bound set dimensions");
    }
    update_redundancy_avoidance(bounds_, point);
  }

  /**
   * @brief Automatically selects the best algorithm based on the number of
   *        objectives p, following the empirical analysis from the paper.
   *
   * Selection rule:
   *   - p <= 5: Algorithm 3 (Enhanced RE) — better cache locality, small |A|
   *   - p >= 6: Algorithm 5 (RA, general case) — avoids expensive filtering
   *
   * Algorithm 5 is preferred over Algorithm 4 because it handles both the
   * general position assumption (SA) and general case instances. If the anti-
   * reference point was not provided (single-argument constructor), Algorithm 3
   * is used as fallback regardless of p.
   *
   * @param point The new nondominated point.
   */
  void update_auto(const Point<T>& point) {
    if (point.dimensions() != dimensions_) {
      throw std::invalid_argument("Point dimensions must match bound set dimensions");
    }
    if (dimensions_ >= 6 && has_anti_reference()) {
      update_redundancy_avoidance(bounds_, point);
    } else {
      update_redundancy_elimination_enhanced(bounds_, point);
    }
  }

  /**
   * @brief Updates the bound set using the naive approach.
   *
   * This is a brute-force algorithm that generates all candidate bounds
   * and filters by checking dominance against the entire set. Useful for
   * correctness verification of the more efficient algorithms.
   *
   * @param point The new nondominated point.
   */
  void update_naive(const Point<T>& point) {
    if (point.dimensions() != dimensions_) {
      throw std::invalid_argument("Point dimensions must match bound set dimensions");
    }
    _update_naive(bounds_, point);
  }

  /**
   * @brief Returns the current set of local bounds.
   */
  [[nodiscard]] const std::vector<LocalBound<T>>& bounds() const {
    return bounds_;
  }

  /**
   * @brief Returns the number of local bounds.
   */
  [[nodiscard]] std::size_t size() const {
    return bounds_.size();
  }

  /**
   * @brief Returns the dimensionality of the objective space.
   */
  [[nodiscard]] std::size_t dimensions() const {
    return dimensions_;
  }

  /**
   * @brief Checks if a point is in the search region.
   *
   * For MINIMIZE: returns true if the point is strictly dominated by some bound.
   * For MAXIMIZE: returns true if the point strictly dominates some bound.
   *
   * @param point The point to check.
   * @return true if the point is in the search region.
   */
  [[nodiscard]] bool is_in_search_region(const std::vector<T>& point) const {
    if (point.size() != dimensions_) {
      throw std::invalid_argument("Point dimensions must match bound set dimensions");
    }
    for (const auto& bound : bounds_) {
      if (strictly_dominates<T, Sense>(point, bound.coordinates)) {
        return true;
      }
    }
    return false;
  }

  /**
   * @brief Finds the bound whose search zone contains the given point.
   *
   * @param point The point to locate.
   * @return Optional containing the bound if found, empty otherwise.
   */
  [[nodiscard]] std::optional<LocalBound<T>> find_containing_bound(
      const std::vector<T>& point) const {
    for (const auto& bound : bounds_) {
      if (strictly_dominates<T, Sense>(point, bound.coordinates)) {
        return bound;
      }
    }
    return std::nullopt;
  }

 private:
  /**
   * @brief Returns true if the anti-reference point was provided at construction.
   */
  [[nodiscard]] bool has_anti_reference() const {
    return !anti_reference_.empty();
  }

  std::size_t dimensions_;
  std::vector<T> reference_point_;
  std::vector<T> anti_reference_;
  std::vector<LocalBound<T>> bounds_;
  std::size_t next_bound_id_;

  /**
   * @brief Generates a unique ID for a new bound.
   */
  std::string generate_bound_id() {
    return "u" + std::to_string(next_bound_id_++);
  }

  /**
   * @brief Algorithm 2: Update procedure based on redundancy elimination.
   *
   * This is an in-place update that modifies current_bounds.
   * Optimized version: uses indices instead of copies, pre-allocates vectors.
   *
   * @param current_bounds Current set of local bounds (modified in place).
   * @param point New point z̄.
   */
  void update_redundancy_elimination(
      std::vector<LocalBound<T>>& current_bounds, const Point<T>& point) {
    // Step 1: Find strongly dominated bounds (Set A) and weakly dominated bounds (Set B)
    // A contains bounds strictly dominated by the point (z < u)
    // B contains bounds weakly dominated by the point (z <= u, but not z < u)
    std::vector<std::size_t> A_idx;
    std::vector<std::size_t> B_idx;

    // Reserving small capacity to prevent minor reallocations
    A_idx.reserve(current_bounds.size());
    B_idx.reserve(current_bounds.size());

    for (std::size_t i = 0; i < current_bounds.size(); ++i) {
      if (strictly_dominates<T, Sense>(point.coordinates, current_bounds[i].coordinates)) {
        A_idx.push_back(i);  // u in A
      } else if (weakly_dominates<T, Sense>(point.coordinates, current_bounds[i].coordinates)) {
        B_idx.push_back(i);  // u in B
      }
    }

    // If no search zones contain the new point, the bound set is unaffected
    if (A_idx.empty()) return;

    // Step 2: Generate all candidate bounds (projections) into Set P
    std::vector<LocalBound<T>> P;
    P.reserve(A_idx.size() * dimensions_);

    for (std::size_t idx : A_idx) {
      const auto& u = current_bounds[idx];
      for (std::size_t j = 0; j < dimensions_; ++j) {
        std::string new_id = u.id + std::to_string(j + 1);
        std::vector<T> new_coords = u.coordinates;
        new_coords[j] = point.coordinates[j];
        P.emplace_back(new_id, std::move(new_coords));
      }
    }

    std::vector<std::size_t> redundant_idx;
    redundant_idx.reserve(P.size());

    // Step 3: Filter out candidates in P that are weakly dominated by other candidates in P
    for (std::size_t i = 0; i < P.size(); ++i) {
      for (std::size_t j = i + 1; j < P.size(); ++j) {
        if (weakly_dominates<T, Sense>(P[i].coordinates, P[j].coordinates)) {
          redundant_idx.push_back(i);  // P[i] is redundant
        } else if (weakly_dominates<T, Sense>(P[j].coordinates, P[i].coordinates)) {
          redundant_idx.push_back(j);  // P[j] is redundant
        }
      }
    }

    // Step 4: Filter out candidates in P that are weakly dominated by existing bounds in B
    // By only checking against B, we avoid useless dominance checks against all of U(N)\A
    for (std::size_t i = 0; i < P.size(); ++i) {
      for (std::size_t b_idx : B_idx) {
        if (weakly_dominates<T, Sense>(P[i].coordinates, current_bounds[b_idx].coordinates)) {
          redundant_idx.push_back(i);
          break;  // P[i] is already marked for deletion, move to the next candidate
        }
      }
    }

    // Step 5: Execute deletions on P using the reverse sorted swap-and-pop technique
    std::sort(redundant_idx.begin(), redundant_idx.end());
    redundant_idx.erase(std::unique(redundant_idx.begin(), redundant_idx.end()), redundant_idx.end());

    for (auto it = redundant_idx.rbegin(); it != redundant_idx.rend(); ++it) {
      std::swap(P[*it], P.back());
      P.pop_back();
    }

    // Step 6: Remove the original dominated bounds (Set A) from current_bounds
    // Because A_idx was populated in ascending order, reverse iteration is perfectly safe here.
    for (auto it = A_idx.rbegin(); it != A_idx.rend(); ++it) {
      std::swap(current_bounds[*it], current_bounds.back());
      current_bounds.pop_back();
    }

    // Step 7: Append the filtered new bounds (Set P) into current_bounds
    current_bounds.insert(current_bounds.end(),
                          std::make_move_iterator(P.begin()),
                          std::make_move_iterator(P.end()));
  }

  /**
   * @brief Algorithm 3: update with enhanced filtering.
   *
   * This implementation follows Algorithm 3 in Section 3 of the paper:
   * 1) A = {u in U(N): z̄ < u}
   * 2) B_j = {u in U(N): z̄_j = u_j and z̄_{-j} < u_{-j}}
   * 3) P_j = {(z̄_j, u_{-j}) : u in A}
   * 4) P_j <- {p in P_j : p ≰ u0, for all u0 in P_j ∪ B_j}
   * 5) U(N ∪ {z̄}) = (U(N) \ A) ∪ (⋃_j P_j)
   *
   * @param current_bounds Current set of local bounds (modified in place).
   * @param point New point z̄.
   */
  void update_redundancy_elimination_enhanced(
      std::vector<LocalBound<T>>& current_bounds, const Point<T>& point) {
    // Step 1: Find A (strictly dominated bounds)
    std::vector<std::size_t> A_idx;
    A_idx.reserve(current_bounds.size());
    for (std::size_t i = 0; i < current_bounds.size(); ++i) {
      if (strictly_dominates<T, Sense>(point.coordinates, current_bounds[i].coordinates)) {
        A_idx.push_back(i);
      }
    }

    if (A_idx.empty()) return;

    std::sort(A_idx.begin(), A_idx.end());

    // Step 2: Build B_j with strict inequality on all non-j dimensions
    std::vector<std::vector<LocalBound<T>>> B(dimensions_);
    for (std::size_t j = 0; j < dimensions_; ++j) {
      for (const auto& u : current_bounds) {
        if (point.coordinates[j] != u.coordinates[j]) {
          continue;
        }

        bool strict_all_other = true;
        for (std::size_t k = 0; k < dimensions_; ++k) {
          if (k == j) continue;
          if (!detail::is_better<T, Sense>(point.coordinates[k], u.coordinates[k])) {
            strict_all_other = false;
            break;
          }
        }

        if (strict_all_other) {
          B[j].push_back(u);
        }
      }
    }

    // Step 3: Generate P_j from A
    std::vector<std::vector<LocalBound<T>>> P(dimensions_);
    for (std::size_t j = 0; j < dimensions_; ++j) {
      P[j].reserve(A_idx.size());
    }

    for (std::size_t idx : A_idx) {
      const auto& u = current_bounds[idx];
      for (std::size_t j = 0; j < dimensions_; ++j) {
        std::string new_id = u.id + std::to_string(j + 1);
        std::vector<T> new_coords = u.coordinates;
        new_coords[j] = point.coordinates[j];
        P[j].emplace_back(new_id, std::move(new_coords));
      }
    }

    // Step 4: Per-dimension filtering against P_j and B_j
    for (std::size_t j = 0; j < dimensions_; ++j) {
      std::vector<LocalBound<T>> filtered;
      filtered.reserve(P[j].size());

      for (const auto& p : P[j]) {
        bool redundant = false;

        // Check against other candidates in same P_j
        for (const auto& u0 : P[j]) {
          if (p != u0 && weakly_dominates<T, Sense>(p.coordinates, u0.coordinates)) {
            redundant = true;
            break;
          }
        }

        // Check against B_j
        if (!redundant) {
          for (const auto& u0 : B[j]) {
            if (weakly_dominates<T, Sense>(p.coordinates, u0.coordinates)) {
              redundant = true;
              break;
            }
          }
        }

        if (!redundant) {
          filtered.push_back(p);
        }
      }

      P[j] = std::move(filtered);
    }

    // Step 5: U' = (U \ A) ∪ (⋃_j P_j)
    for (auto it = A_idx.rbegin(); it != A_idx.rend(); ++it) {
      std::swap(current_bounds[*it], current_bounds.back());
      current_bounds.pop_back();
    }

    std::size_t total_candidates = 0;
    for (const auto& pj : P) {
      total_candidates += pj.size();
    }
    current_bounds.reserve(current_bounds.size() + total_candidates);

    for (auto& pj : P) {
      current_bounds.insert(current_bounds.end(),
                            std::make_move_iterator(pj.begin()),
                            std::make_move_iterator(pj.end()));
    }
  }

  /**
   * @brief Algorithm 4: Update procedure based on redundancy avoidance (SA).
   *
   * This algorithm avoids redundancy checks by keeping track of the sets Z^j(u)
   * of points that define each component of each local bound.
   *
   * @param current_bounds Current set of local bounds (modified in place).
   * @param def_sets Current defining sets (modified in place).
   * @param point New point z̄.
   */
  void update_redundancy_avoidance_sa(
      std::vector<LocalBound<T>>& current_bounds, const Point<T>& z_bar) {
    // Step 1: A ← {u ∈ U(N) : z̄ < u}
    std::vector<std::size_t> A_idx;
    for (std::size_t i = 0; i < current_bounds.size(); ++i) {
      if (strictly_dominates<T, Sense>(z_bar.coordinates, current_bounds[i].coordinates)) {
        A_idx.push_back(i);
      }
    }

    if (A_idx.empty()) return;

    // Step 2: P ← ∅
    std::vector<LocalBound<T>> P;

    // Steps 3–10: For each u ∈ A, for each j ∈ {1,...,p}
    for (std::size_t idx : A_idx) {
      const auto& u = current_bounds[idx];

      for (std::size_t j = 0; j < dimensions_; ++j) {
        // Step 5: z_j^max(u) = max_{k≠j} { z^k(u)_j }
        // z^k(u) is the defining point for component k of bound u.
        // We look at its j-th coordinate.
        T z_max_j;
        if constexpr (Sense == Objective::MINIMIZE) {
          z_max_j = std::numeric_limits<T>::lowest();
          for (std::size_t k = 0; k < dimensions_; ++k) {
            if (k == j) continue;
            T val = u.defining_points[k].coordinates[j];
            if (val > z_max_j) z_max_j = val;
          }
        } else {
          // For MAXIMIZE: z_j^max becomes z_j^min (we flip the comparison)
          z_max_j = std::numeric_limits<T>::max();
          for (std::size_t k = 0; k < dimensions_; ++k) {
            if (k == j) continue;
            T val = u.defining_points[k].coordinates[j];
            if (val < z_max_j) z_max_j = val;
          }
        }

        // Step 6: Check condition of Theorem 4.3
        // MIN: z̄_j > z_j^max(u)   MAX: z̄_j < z_j^min(u)
        bool keep;
        if constexpr (Sense == Objective::MINIMIZE) {
          keep = z_bar.coordinates[j] > z_max_j;
        } else {
          keep = z_bar.coordinates[j] < z_max_j;
        }

        if (keep) {
          // Let u_j = (z̄_j, u_{-j})
          LocalBound<T> u_j;
          u_j.id = u.id + std::to_string(j + 1);
          u_j.coordinates = u.coordinates;
          u_j.coordinates[j] = z_bar.coordinates[j];

          // Step 8: z^j(u_j) ← z̄
          u_j.defining_points = u.defining_points;
          u_j.defining_points[j] = z_bar;
          // Steps 9-10: For k ≠ j, z^k(u_j) ← z^k(u) (inherited)

          P.push_back(std::move(u_j));
        }
      }
    }

    // Step 11: U(N ∪ {z̄}) ← (U(N) \ A) ∪ P
    for (auto it = A_idx.rbegin(); it != A_idx.rend(); ++it) {
      std::swap(current_bounds[*it], current_bounds.back());
      current_bounds.pop_back();
    }
    current_bounds.insert(current_bounds.end(),
                          std::make_move_iterator(P.begin()),
                          std::make_move_iterator(P.end()));
  }

  /**
   * @brief Algorithm 5: Update procedure based on redundancy avoidance (general case).
   *
   * This algorithm avoids redundancy checks by keeping track of the sets Z^j(u)
   * of points that define each component of each local bound.
   *
   * @param current_bounds Current set of local bounds (modified in place).
   * @param def_sets Current defining sets (modified in place).
   * @param point New point z̄.
   */
  void update_redundancy_avoidance(
      std::vector<LocalBound<T>>& current_bounds,
      const Point<T>& z_bar) {
    // Step 1: A ← {u ∈ U(N) : z̄ < u}
    std::vector<std::size_t> A_idx;
    for (std::size_t i = 0; i < current_bounds.size(); ++i) {
      if (strictly_dominates<T, Sense>(z_bar.coordinates, current_bounds[i].coordinates)) {
        A_idx.push_back(i);
      }
    }

    // Steps 3–4: Update Z^j(u) for ALL bounds u ∈ U(N)
    // Per Proposition 4.1: if z̄_j = u_j and z̄_{-j} < u_{-j}, add z̄ to Z^j(u)
    for (auto& u : current_bounds) {
      for (std::size_t j = 0; j < dimensions_; ++j) {
        if (z_bar.coordinates[j] == u.coordinates[j]) {
          bool all_strict = true;
          for (std::size_t k = 0; k < dimensions_; ++k) {
            if (k == j) continue;
            if (!detail::is_better<T, Sense>(z_bar.coordinates[k], u.coordinates[k])) {
              all_strict = false;
              break;
            }
          }
          if (all_strict) {
            u.defining_point_sets[j].push_back(z_bar);
          }
        }
      }
    }

    if (A_idx.empty()) return;

    // Step 2: P ← ∅
    std::vector<LocalBound<T>> P;

    // Steps 5–12: For each u ∈ A, for each j ∈ {1,...,p}
    for (std::size_t idx : A_idx) {
      const auto& u = current_bounds[idx];

      for (std::size_t j = 0; j < dimensions_; ++j) {
        // Step 7: z_j^max(u) = max_{k≠j} min{z_j : z ∈ Z^k(u)}
        T max_of_mins;
        if constexpr (Sense == Objective::MINIMIZE) {
          max_of_mins = std::numeric_limits<T>::lowest();
          for (std::size_t k = 0; k < dimensions_; ++k) {
            if (k == j || u.defining_point_sets[k].empty()) continue;
            T min_val = std::numeric_limits<T>::max();
            for (const auto& z : u.defining_point_sets[k]) {
              if (z.coordinates[j] < min_val) min_val = z.coordinates[j];
            }
            if (min_val > max_of_mins) max_of_mins = min_val;
          }
        } else {
          // MAXIMIZE: z_j^max(u) = min_{k≠j} max{z_j : z ∈ Z^k(u)}
          max_of_mins = std::numeric_limits<T>::max();
          for (std::size_t k = 0; k < dimensions_; ++k) {
            if (k == j || u.defining_point_sets[k].empty()) continue;
            T max_val = std::numeric_limits<T>::lowest();
            for (const auto& z : u.defining_point_sets[k]) {
              if (z.coordinates[j] > max_val) max_val = z.coordinates[j];
            }
            if (max_val < max_of_mins) max_of_mins = max_val;
          }
        }

        // Step 8: Check condition of Theorem 4.4
        // MIN: z̄_j > z_j^max(u)   MAX: z̄_j < z_j^max(u)
        bool keep;
        if constexpr (Sense == Objective::MINIMIZE) {
          keep = z_bar.coordinates[j] > max_of_mins;
        } else {
          keep = z_bar.coordinates[j] < max_of_mins;
        }

        if (keep) {
          // Let u_j = (z̄_j, u_{-j})
          LocalBound<T> u_j;
          u_j.id = u.id + std::to_string(j + 1);
          u_j.coordinates = u.coordinates;
          u_j.coordinates[j] = z_bar.coordinates[j];

          // Step 10: Z^j(u_j) ← {z̄}
          u_j.defining_point_sets.resize(dimensions_);
          u_j.defining_point_sets[j] = {z_bar};

          // Steps 11–12: Z^k(u_j) ← {z ∈ Z^k(u) : z_j < z̄_j}  (MIN)
          //               Z^k(u_j) ← {z ∈ Z^k(u) : z_j > z̄_j}  (MAX)
          for (std::size_t k = 0; k < dimensions_; ++k) {
            if (k == j) continue;
            for (const auto& z : u.defining_point_sets[k]) {
              if constexpr (Sense == Objective::MINIMIZE) {
                if (z.coordinates[j] < z_bar.coordinates[j]) {
                  u_j.defining_point_sets[k].push_back(z);
                }
              } else {
                if (z.coordinates[j] > z_bar.coordinates[j]) {
                  u_j.defining_point_sets[k].push_back(z);
                }
              }
            }
          }

          P.push_back(std::move(u_j));
        }
      }
    }

    // Step 13: U(N ∪ {z̄}) ← (U(N) \ A) ∪ P
    for (auto it = A_idx.rbegin(); it != A_idx.rend(); ++it) {
      std::swap(current_bounds[*it], current_bounds.back());
      current_bounds.pop_back();
    }
    current_bounds.insert(current_bounds.end(),
                          std::make_move_iterator(P.begin()),
                          std::make_move_iterator(P.end()));
  }

  /**
   * @brief Naive O(n²) update algorithm for correctness verification.
   *
   * Generates all candidate bounds and explicitly filters redundant ones by
   * checking weak dominance against all other bounds.
   *
   * @param current_bounds Current set of local bounds (modified in place).
   * @param point New point z̄.
   */
  void _update_naive(std::vector<LocalBound<T>>& current_bounds, const Point<T>& point) {
    // Step 1: Find strongly dominated bounds and their indices
    std::vector<LocalBound<T>> A;
    std::vector<std::size_t> A_idx;
    for (std::size_t i = 0; i < current_bounds.size(); ++i) {
      if (strictly_dominates<T, Sense>(point.coordinates, current_bounds[i].coordinates)) {
        A.push_back(current_bounds[i]);
        A_idx.push_back(i);
      }
    }

    if (A.empty()) return;

    // Step 2: Generate all candidate bounds (projections)
    std::vector<LocalBound<T>> new_bounds;
    for (const auto& u : A) {
      for (std::size_t j = 0; j < dimensions_; ++j) {
        std::string new_id = u.id + std::to_string(j + 1);
        std::vector<T> new_coords = u.coordinates;
        new_coords[j] = point.coordinates[j];
        new_bounds.emplace_back(new_id, std::move(new_coords));
      }
    }

    // Step 3: Remove dominated bounds from current_bounds
    for (auto it = A_idx.rbegin(); it != A_idx.rend(); ++it) {
      std::swap(current_bounds[*it], current_bounds.back());
      current_bounds.pop_back();
    }

    // Step 4.1: Find new bounds weakly dominated by other new bounds
    // A bound is redundant when its search zone is a subset of another bound's zone.
    //
    // - Minimization (upper bounds): the search zone of u is {z : z[i] < u[i] ∀i}.
    //   If u_a ≤ u_b component-wise, then zone(u_a) ⊆ zone(u_b), so u_a is redundant.
    //   We keep the larger bound (u_b) and remove the smaller one (u_a).
    //
    // - Maximization (lower bounds): the search zone of u is {z : z[i] > u[i] ∀i}.
    //   If u_a ≥ u_b component-wise, then zone(u_a) ⊆ zone(u_b), so u_a is redundant.
    //   We keep the smaller bound (u_b) and remove the larger one (u_a).
    //
    // In both cases, weak dominance of u_a over u_b implies u_a is redundant.
    std::vector<std::size_t> redundant_idx;
    for (std::size_t i = 0; i < new_bounds.size(); ++i) {
      for (std::size_t j = i + 1; j < new_bounds.size(); ++j) {
        if (weakly_dominates<T, Sense>(new_bounds[i].coordinates, new_bounds[j].coordinates)) {
          redundant_idx.push_back(i);
        } else if (weakly_dominates<T, Sense>(new_bounds[j].coordinates, new_bounds[i].coordinates)) {
          redundant_idx.push_back(j);
        }
      }
    }

    // Step 4.2: Find new bounds weakly dominated by remaining bounds
    for (std::size_t i = 0; i < new_bounds.size(); ++i) {
      for (const auto& b : current_bounds) {
        if (weakly_dominates<T, Sense>(new_bounds[i].coordinates, b.coordinates)) {
          redundant_idx.push_back(i);
          break;
        }
      }
    }

    // Step 4.3: Remove redundant new bounds (in reverse order for safe removal)
    std::sort(redundant_idx.begin(), redundant_idx.end());
    redundant_idx.erase(std::unique(redundant_idx.begin(), redundant_idx.end()), redundant_idx.end());
    for (auto it = redundant_idx.rbegin(); it != redundant_idx.rend(); ++it) {
      std::swap(new_bounds[*it], new_bounds.back());
      new_bounds.pop_back();
    }

    // Step 5: Add new bounds
    current_bounds.insert(current_bounds.end(), new_bounds.begin(), new_bounds.end());
  }
};

}  // namespace local_bounds

#endif  // LOCAL_BOUNDS_BOUND_SET_HPP
