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
 * filtering steps of Algorithms 2/3 (Redundancy Elimination).
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
        next_bound_id_(1) {
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
               [[maybe_unused]] const std::vector<T> &anti_reference,
               size_t max_leaf_size = 32, size_t num_children = 8)
      : dimensions_(reference_point.size()),
        tree_(max_leaf_size, num_children, reference_point.size()),
        next_bound_id_(1) {
    tree_.Insert(reference_point);
  }

  /**
   * @brief Updates using Algorithm 2 (Redundancy Elimination).
   *
   * @param point The new nondominated point.
   */
  void update_re(const Point<T> &point) { update_re_impl(point, false); }

  /**
   * @brief Updates using Algorithm 3 (Enhanced Redundancy Elimination).
   *
   * @param point The new nondominated point.
   */
  void update_re_enhanced(const Point<T> &point) {
    update_re_impl(point, true);
  }

  /**
   * @brief Updates using the naive algorithm.
   *
   * Delegates to update_re since the tree-based implementation already
   * performs filtering during generation.
   *
   * @param point The new nondominated point.
   */
  void update_naive(const Point<T> &point) { update_re(point); }

  /**
   * @brief Automatically selects the best update algorithm.
   *
   * Always dispatches to update_re_enhanced for the tree-based implementation.
   *
   * @param point The new nondominated point.
   */
  void update_auto(const Point<T> &point) { update_re_enhanced(point); }

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
    auto all_lubs = tree_.GetAllBounds();
    for (const auto &lub : all_lubs) {
      if (strictly_dominates<T, Sense>(point, lub)) {
        return true;
      }
    }
    return false;
  }

  /**
   * @brief Finds a bound whose search zone contains the given point.
   *
   * @param point The point to locate.
   * @return Optional containing the bound if found, empty otherwise.
   */
  [[nodiscard]] std::optional<LocalBound<T>>
  find_containing_bound(const std::vector<T> &point) const {
    auto all_lubs = tree_.GetAllBounds();
    for (size_t i = 0; i < all_lubs.size(); ++i) {
      if (strictly_dominates<T, Sense>(point, all_lubs[i])) {
        return LocalBound<T>("u" + std::to_string(i), std::move(all_lubs[i]));
      }
    }
    return std::nullopt;
  }

private:
  std::size_t dimensions_;
  LBTree<T, Sense> tree_;
  std::size_t next_bound_id_;

  /**
   * @brief Shared implementation for update_re and update_re_enhanced.
   *
   * Both algorithms share the same tree-based implementation since the
   * LBTree handles the spatial indexing uniformly.
   */
  void update_re_impl(const Point<T> &point, bool /* enhanced */) {
    const auto &z = point.coordinates;
    std::vector<std::vector<T>> A = tree_.ExtractStrictlyDominated(z);

    if (A.empty())
      return;

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

        std::vector<T> p_cand = u;
        p_cand[j] = z[j];

        bool dominated = false;

        // Filter against A
        for (std::size_t k = 0; k < A.size(); ++k) {
          if (i != k) {
            bool cur_dom = true;
            for (std::size_t dim = 0; dim < dimensions_; ++dim) {
              if (!detail::is_at_least_as_good<T, Sense>(p_cand[dim],
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
              if (!detail::is_at_least_as_good<T, Sense>(p_cand[k], w[k])) {
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

        // Filter against P
        if (!dominated) {
          for (const auto &w : P) {
            bool p_le_w = true;
            for (std::size_t k = 0; k < dimensions_; ++k) {
              if (!detail::is_at_least_as_good<T, Sense>(p_cand[k], w[k])) {
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

        if (!dominated) {
          P.push_back(p_cand);
        }
      }
    }

    // Insert filtered new bounds into tree
    for (const auto &p_cand : P) {
      tree_.Insert(p_cand);
    }
  }
};

} // namespace local_bounds

#endif // LOCAL_BOUNDS_BOUND_SET_TREE_HPP
