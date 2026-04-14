#ifndef LINEAR_LIST_HPP
#define LINEAR_LIST_HPP

#include "../local_bounds/dominance.hpp"
#include "../local_bounds/types.hpp"
#include <cstddef> // for size_t
#include <vector>

namespace local_bounds {

/**
 * @brief O(N) baseline for non-dominance checking (Pareto archive).
 *
 * Maintains a set of mutually non-dominated points using a simple linear list.
 * Useful as a correctness reference for tree-based structures.
 *
 * @tparam T    Coordinate type (e.g., double, int64_t).
 * @tparam Sense Optimization objective (MINIMIZE or MAXIMIZE).
 */
template <typename T = double, Objective Sense = Objective::MINIMIZE>
class LinearList {
public:
  using Point = std::vector<T>;

  /**
   * @brief Construct a new Linear List object.
   *
   * @param p The number of objectives (dimensions) for each point.
   */
  explicit LinearList(size_t p) : p_(p) {}

  /**
   * @brief Updates the linear list with a new candidate point `y`.
   *
   * This is an O(N) iterative process that compares `y` against every existing
   * point:
   * 1. If `y` is covered by any existing point, it is instantly rejected.
   * 2. If `y` enters the archive, any existing point dominated by `y` is
   * removed.
   *
   * @param y A candidate point to insert.
   * @param pruned Optional pointer to a vector that will be populated with all
   *               points that were pruned/dominated by `y`. When nullptr
   *               (default), no collection overhead is introduced.
   * @return true if the point was inserted successfully, false if it was
   * rejected.
   */
  bool Update(const Point &y, std::vector<Point> *pruned = nullptr) {
    // check if y is covered (weakly dominated by an existing point)
    for (const auto &p : archive_) {
      if (weakly_dominates<T, Sense>(p, y))
        return false;
    }

    // y is not covered, add it and remove dominated
    auto it = archive_.begin();
    while (it != archive_.end()) {
      if (dominates<T, Sense>(y, *it)) {
        if (pruned)
          pruned->push_back(*it);
        it = archive_.erase(it);
      } else {
        ++it;
      }
    }
    archive_.push_back(y);
    return true;
  }

  /**
   * @brief Returns a copy of the non-dominated Pareto archive.
   *
   * @return std::vector<Point> The points residing in the list.
   */
  std::vector<Point> GetPoints() const { return archive_; }

private:
  std::vector<Point> archive_;
  size_t p_;
};

} // namespace local_bounds

#endif // LINEAR_LIST_HPP
