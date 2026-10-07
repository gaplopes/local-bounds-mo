#ifndef LOCAL_BOUNDS_BOUND_SET_TREE_HPP
#define LOCAL_BOUNDS_BOUND_SET_TREE_HPP

#include <optional>
#include <string>
#include <vector>

#include "../structures/lb_tree.hpp"
#include "dominance.hpp"
#include "types.hpp"

namespace local_bounds {

/**
 * @brief Tree-accelerated bound set using a LUBTree spatial index.
 *
 * Provides the same API as BoundSet but uses a LBTree to accelerate the
 * filtering steps of RE. Enhanced and naive methods are API aliases.
 *
 * @note Algorithms 4/5 (Redundancy Avoidance) require tracking defining point
 *       sets (Z^j(u)) which a spatial index does not support. Use BoundSet
 *       directly for those algorithms.
 *
 * @tparam T     Coordinate type (e.g., double, int64_t).
 * @tparam Sense Optimization objective (MINIMIZE or MAXIMIZE).
 */
template <typename T = double, Objective Sense = Objective::MINIMIZE>
class BoundSetTree {
public:
  /**
   * @brief Constructs a BoundSetTree with the given reference point.
   *
   * @param reference_point For MINIMIZE: the nadir point (upper bound of search
   *                        space). For MAXIMIZE: the ideal point (lower bound
   *                        of search space).
   * @param max_leaf_size   Maximum points per LBTree leaf before splitting
   *                        (default: 32).
   * @param num_children    Number of children created on split (default: 8).
   */
  explicit BoundSetTree(const std::vector<T> &reference_point,
                        size_t max_leaf_size = 32, size_t num_children = 8)
      : dimensions_(reference_point.size()),
        tree_(max_leaf_size, num_children, reference_point.size()),
        reference_point_(reference_point) {
    detail::validate_coordinates(reference_point, dimensions_);
    tree_.Insert(reference_point);
  }

  /**
   * @brief Constructs a BoundSetTree with reference and anti-reference points.
   *
   * @param reference_point For MINIMIZE: the nadir point M.
   *                        For MAXIMIZE: the ideal point m.
   * @param anti_reference  For MINIMIZE: the ideal point m.
   *                        For MAXIMIZE: the nadir point M.
   * @param max_leaf_size   Maximum points per LBTree leaf before splitting.
   * @param num_children    Number of children created on split.
   */
  BoundSetTree(const std::vector<T> &reference_point,
               const std::vector<T> &anti_reference,
               size_t max_leaf_size = 32, size_t num_children = 8)
      : dimensions_(reference_point.size()),
        tree_(max_leaf_size, num_children, reference_point.size()),
        reference_point_(reference_point), anti_reference_(anti_reference) {
    detail::validate_interval<T, Sense>(reference_point, anti_reference);
    tree_.Insert(reference_point);
  }

  /**
   * @brief Updates using Algorithm 2 (Redundancy Elimination).
   *
   * @param point The new nondominated point.
   * @return true if the bound set was updated, false if the point did not
   * dominate any local bound.
   */
  bool update_re(const Point<T> &point) { return update_re_impl(point); }

  /**
   * @brief API-compatible alias of this class's indexed RE implementation.
   *
   * @param point The new nondominated point.
   * @return true if the bound set was updated, false if the point did not
   * dominate any local bound.
   */
  bool update_re_enhanced(const Point<T> &point) {
    return update_re_impl(point);
  }

  /**
   * @brief Updates using the naive algorithm.
   *
   * Delegates to update_re since the tree-based implementation already
   * performs filtering during generation.
   *
   * @param point The new nondominated point.
   * @return true if the bound set was updated, false if the point did not
   * dominate any local bound.
   */
  bool update_naive(const Point<T> &point) { return update_re(point); }

  /**
   * @brief Automatically selects the best update algorithm.
   *
   * Always dispatches to update_re_enhanced for the tree-based implementation.
   *
   * @param point The new nondominated point.
   * @return true if the bound set was updated, false if the point did not
   * dominate any local bound.
   */
  bool update_auto(const Point<T> &point) { return update_re_enhanced(point); }

  /**
   * @brief Returns the current set of local bounds.
   */
  [[nodiscard]] std::vector<LocalBound<T>> bounds() const {
    auto all_lubs = tree_.GetAllBounds();
    std::vector<LocalBound<T>> result;
    result.reserve(all_lubs.size());
    for (size_t i = 0; i < all_lubs.size(); ++i) {
      result.emplace_back("u" + std::to_string(i), std::move(all_lubs[i]));
    }
    return result;
  }

  /**
   * @brief Returns the number of local bounds.
   */
  [[nodiscard]] std::size_t size() const { return tree_.Size(); }

  /**
   * @brief Returns the dimensionality of the objective space.
   */
  [[nodiscard]] std::size_t dimensions() const { return dimensions_; }

  /**
   * @brief Checks if a point is in the search region.
   *
   * For MINIMIZE: returns true if the point strictly dominates some bound.
   * For MAXIMIZE: returns true if the point strictly dominates some bound.
   *
   * @param point The point to check.
   * @return true if the point is in the search region.
   */
  [[nodiscard]] bool is_in_search_region(const std::vector<T> &point) const {
    detail::validate_coordinates(point, dimensions_);
    if (!detail::in_interval<T, Sense>(point, reference_point_, anti_reference_)) return false;
    return tree_.FindStrictlyDominated(point) != nullptr;
  }

  /**
   * @brief Finds a bound whose search zone contains the given point.
   *
   * @param point The point to locate.
   * @return Optional containing the bound if found, empty otherwise.
   */
  [[nodiscard]] std::optional<LocalBound<T>>
  find_containing_bound(const std::vector<T> &point) const {
    detail::validate_coordinates(point, dimensions_);
    if (!detail::in_interval<T, Sense>(point, reference_point_, anti_reference_)) return std::nullopt;
    size_t ordinal = 0;
    if (const auto *bound = tree_.FindStrictlyDominated(point, &ordinal))
      return LocalBound<T>("u" + std::to_string(ordinal), *bound);
    return std::nullopt;
  }

private:
  std::size_t dimensions_;
  LBTree<T, Sense> tree_;
  std::vector<T> reference_point_;
  std::vector<T> anti_reference_;

  /**
   * @brief Shared implementation for update_re and update_re_enhanced.
   *
   * Both algorithms share the same tree-based implementation since the
   * LBTree handles the spatial indexing uniformly.
   */
  bool update_re_impl(const Point<T> &point) {
    detail::validate_update<T, Sense>(point.coordinates, reference_point_, anti_reference_);
    const auto &z = point.coordinates;
    std::vector<std::vector<T>> A = tree_.ExtractStrictlyDominated(z);

    if (A.empty())
      return false;

    // Extract B_j
    std::vector<std::vector<std::vector<T>>> B(dimensions_);
    for (std::size_t j = 0; j < dimensions_; ++j) {
      tree_.FindBoundsWithEqualComponent(z, j, B[j]);
    }

    // Generate candidate bounds and filter out redundant ones
    std::vector<std::vector<T>> P;

    for (std::size_t i = 0; i < A.size(); ++i) {
      const auto &u = A[i];
      for (std::size_t j = 0; j < dimensions_; ++j) {
        if (!detail::is_better<T, Sense>(z[j], u[j]))
          continue;

        bool dominated = false;

        // Filter against A
        for (std::size_t k = 0; k < A.size(); ++k) {
          if (i != k) {
            bool cur_dom = true;
            for (std::size_t dim = 0; dim < dimensions_; ++dim) {
              if (!detail::is_at_least_as_good<T, Sense>(dim == j ? z[j] : u[dim],
                                                         A[k][dim])) {
                cur_dom = false;
                break;
              }
            }
            if (cur_dom) {
              dominated = true;
              break;
            }
          }
        }

        // Filter against B_j
        if (!dominated) {
          for (const auto &w : B[j]) {
            bool p_le_w = true;
            for (std::size_t k = 0; k < dimensions_; ++k) {
              if (!detail::is_at_least_as_good<T, Sense>(k == j ? z[j] : u[k], w[k])) {
                p_le_w = false;
                break;
              }
            }
            if (p_le_w) {
              dominated = true;
              break;
            }
          }
        }

        // A already rejects candidates covered by another source's projection.
        if (!dominated) {
          std::vector<T> p_cand = u;
          p_cand[j] = z[j];
          P.push_back(std::move(p_cand));
        }
      }
    }

    // Insert filtered new bounds into tree
    for (const auto &p_cand : P) {
      tree_.Insert(p_cand);
    }
    return true;
  }
};

} // namespace local_bounds

#endif // LOCAL_BOUNDS_BOUND_SET_TREE_HPP
