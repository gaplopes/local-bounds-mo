#ifndef LB_TREE_HPP
#define LB_TREE_HPP

#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <vector>

#include "../local_bounds/dominance.hpp"
#include "../local_bounds/types.hpp"

namespace local_bounds {

template <typename T = double, Objective Sense = Objective::MINIMIZE>
class LBTreeNode {
public:
  using Point = std::vector<T>;

  LBTreeNode() = default;

  bool is_leaf() const { return children.empty(); }

  std::vector<Point> L;
  std::vector<LBTreeNode *> children;
  Point best_point;  // Represents the ideal bound condition (min for MINIMIZE,
                     // max for MAXIMIZE)
  Point worst_point; // Represents the nadir bound condition (max for MINIMIZE,
                     // min for MAXIMIZE)
  LBTreeNode *parent = nullptr;
};

/**
 * @brief A generic spatial tree structure for indexing local bounds
 *
 * This tree indexes a set of bounds based on their coordinates to allow
 * efficient spatial queries (such as extracting strictly dominated bounds).
 * It can be used for both lower bounds and upper bounds.
 */
template <typename T = double, Objective Sense = Objective::MINIMIZE>
class LBTree {
public:
  using Point = std::vector<T>;
  using Node = LBTreeNode<T, Sense>;

  LBTree(size_t max_leaf_size, size_t num_children, size_t num_objectives)
      : max_leaf_size_(max_leaf_size), num_children_(num_children),
        p_(num_objectives) {
    root_ = allocate_node();
  }

  ~LBTree() {}

  /**
   * @brief Extracts and removes all bounds strictly dominated by the given point.
   *
   * @param z_bar The new point that determines which bounds to invalidate.
   * @return A vector of bounds that were strictly dominated and removed.
   */
  std::vector<Point> ExtractStrictlyDominated(const Point &z_bar) {
    std::vector<Point> strictly_dominated;
    if (!root_)
      return strictly_dominated;

    bool delete_root = ExtractRecursive(root_, z_bar, strictly_dominated);
    if (delete_root) {
      free_node(root_);
      root_ = allocate_node();
    }
    item_count_ -= strictly_dominated.size();
    return strictly_dominated;
  }

  /**
   * @brief Finds all bounds `u` that share the exact value at component `j` with `z`,
   *        and are strictly dominated by `z` in all other components.
   *
   * In mathematical terms, finds u such that u_j == z_j and for all k != j, z_k is strictly better than u_k.
   *
   * @param z The reference point.
   * @param j The component index to test for equality.
   * @param out A vector where the matching bounds will be appended.
   */
  void FindBoundsWithEqualComponent(const Point &z, size_t j,
                                    std::vector<Point> &out) const {
    if (root_)
      FindBoundsWithEqualComponentRecursive(root_, z, j, out);
  }

  /**
   * @brief Inserts a new bound into the tree.
   *
   * @param u_new The coordinates of the new bound.
   */
  void Insert(const Point &u_new) {
    if (!root_)
      root_ = allocate_node();
    InsertRecursive(root_, u_new);
    item_count_++;
  }

  /**
   * @brief Returns the total number of bounds currently in the tree.
   */
  [[nodiscard]] size_t Size() const { return item_count_; }

  /**
   * @brief Returns a copy of all bounds currently stored in the tree.
   */
  std::vector<Point> GetAllBounds() const {
    std::vector<Point> result;
    if (root_)
      GetPointsRecursive(root_, result);
    return result;
  }

private:
  size_t max_leaf_size_;
  size_t num_children_;
  size_t p_;
  size_t item_count_ = 0;
  Node *root_;

  std::deque<Node> arena_;
  std::vector<Node *> free_list_;

  /**
   * @brief Compares two coordinate values depending on the optimization sense.
   * @param a The first coordinate.
   * @param b The second coordinate.
   * @return True if 'a' is strictly better than 'b' (e.g., a < b for MINIMIZE).
   */
  bool is_better(const T &a, const T &b) const {
    return local_bounds::detail::is_better<T, Sense>(a, b);
  }

  /**
   * @brief Compares two coordinate values for weak superiority depending on the optimization sense.
   * @param a The first coordinate.
   * @param b The second coordinate.
   * @return True if 'a' is better than or equal to 'b'.
   */
  bool is_at_least_as_good(const T &a, const T &b) const {
    return local_bounds::detail::is_at_least_as_good<T, Sense>(a, b);
  }

  /**
   * @brief Allocates a new node, using the free list if available, or the arena otherwise.
   * @return A pointer to the newly allocated and initialized LBTreeNode.
   */
  Node *allocate_node() {
    if (!free_list_.empty()) {
      Node *node = free_list_.back();
      free_list_.pop_back();
      node->parent = nullptr;
      node->L.clear();
      node->children.clear();
      node->best_point.clear();
      node->worst_point.clear();
      return node;
    }
    arena_.emplace_back();
    return &arena_.back();
  }

  /**
   * @brief Recursively frees a node and its children, returning them to the free list.
   * @param node The root node of the subtree to free.
   */
  void free_node(Node *node) {
    if (!node)
      return;
    for (auto *child : node->children)
      free_node(child);
    node->children.clear();
    free_list_.push_back(node);
  }

  /**
   * @brief Recursively finds and extracts all bounds strictly dominated by z_bar.
   * @param n The current node being processed.
   * @param z_bar The dominating point.
   * @param extracted Active vector storing the extracted bounds.
   * @return True if the subtree rooted at 'n' is entirely extracted/empty and should be deleted.
   */
  bool ExtractRecursive(Node *n, const Point &z_bar,
                        std::vector<Point> &extracted) {
    if (n->best_point.empty() || n->worst_point.empty())
      return false;

    // Prune if z_bar cannot dominate the worst possible bound in this node.
    for (size_t i = 0; i < p_; ++i) {
      if (!is_better(z_bar[i], n->worst_point[i]))
        return false;
    }

    // Check if the entire node is completely dominated.
    bool mass_hit = true;
    for (size_t i = 0; i < p_; ++i) {
      if (!is_better(z_bar[i], n->best_point[i])) {
        mass_hit = false;
        break;
      }
    }

    if (mass_hit) {
      GetPointsRecursive(n, extracted);
      return true;
    }

    if (n->is_leaf()) {
      size_t i = 0;
      while (i < n->L.size()) {
        if (strictly_inside(z_bar, n->L[i])) {
          extracted.push_back(n->L[i]);
          std::swap(n->L[i], n->L.back());
          n->L.pop_back();
        } else {
          ++i;
        }
      }
      if (n->L.empty())
        return true;
    } else {
      size_t i = 0;
      while (i < n->children.size()) {
        if (ExtractRecursive(n->children[i], z_bar, extracted)) {
          free_node(n->children[i]);
          n->children.erase(n->children.begin() + i);
        } else {
          ++i;
        }
      }
      if (n->children.empty())
        return true;
      
      // Compress the tree path if only one child remains.
      if (n->children.size() == 1) {
        auto *remaining = n->children[0];
        n->L = std::move(remaining->L);
        n->best_point = std::move(remaining->best_point);
        n->worst_point = std::move(remaining->worst_point);
        n->children = std::move(remaining->children);
        for (auto *c : n->children)
          c->parent = n;

        remaining->children.clear();
        free_node(remaining);
        return false;
      }
    }

    RecalculateBounds(n);
    return false;
  }

  /**
   * @brief Recursively finds the appropriate bounding box and inserts a new bound into the tree.
   * @param n The current node being inspected.
   * @param u_new The new bound coordinate to insert.
   */
  void InsertRecursive(Node *n, const Point &u_new) {
    if (n->is_leaf()) {
      // Ignore if it's already entirely covered by an existing bound
      for (const auto &u_existing : n->L) {
        if (covers_bound(u_existing, u_new))
          return;
      }
      n->L.push_back(u_new);
      UpdateBestWorst(n, u_new);
      if (n->L.size() > max_leaf_size_)
        Split(n);
    } else {
      // Find the child with the minimal expected distance metric to the new point
      T min_dist = std::numeric_limits<T>::max();
      Node *closest_child = nullptr;
      for (auto *child : n->children) {
        double d = expected_distance(u_new, child);
        if (d < min_dist) {
          min_dist = static_cast<T>(d);
          closest_child = child;
        }
      }
      UpdateBestWorst(n, u_new);
      if (closest_child) {
        InsertRecursive(closest_child, u_new);
      } else {
        // Fallback if no valid children (unexpected case)
        n->L.push_back(u_new);
        n->children.clear();
      }
    }
  }

  /**
   * @brief Recalculates the exact limits (best/worst bounding boxes) of a node.
   * @param n The node to recalculate.
   */
  void RecalculateBounds(Node *n) {
    if (n->is_leaf()) {
      if (n->L.empty())
        return;
      n->best_point = n->L[0];
      n->worst_point = n->L[0];
      for (size_t k = 1; k < n->L.size(); ++k) {
        for (size_t i = 0; i < p_; ++i) {
          if (is_better(n->L[k][i], n->best_point[i]))
            n->best_point[i] = n->L[k][i];
          if (is_better(n->worst_point[i], n->L[k][i]))
            n->worst_point[i] = n->L[k][i];
        }
      }
    } else {
      if (n->children.empty())
        return;
      n->best_point = n->children[0]->best_point;
      n->worst_point = n->children[0]->worst_point;
      for (size_t k = 1; k < n->children.size(); ++k) {
        for (size_t i = 0; i < p_; ++i) {
          if (is_better(n->children[k]->best_point[i], n->best_point[i]))
            n->best_point[i] = n->children[k]->best_point[i];
          if (is_better(n->worst_point[i], n->children[k]->worst_point[i]))
            n->worst_point[i] = n->children[k]->worst_point[i];
        }
      }
    }
  }

  /**
   * @brief Splits an over-filled leaf node into multiple children using a k-means style clustering distance approach.
   * @param n The leaf node to be split.
   */
  void Split(Node *n) {
    if (n->L.size() <= 1)
      return;

    size_t z_idx = 0;
    double max_avg_dist = -1.0;
    for (size_t i = 0; i < n->L.size(); ++i) {
      double dist_sum = 0.0;
      for (size_t j = 0; j < n->L.size(); ++j) {
        if (i != j)
          dist_sum += distance(n->L[i], n->L[j]);
      }
      double avg_dist = dist_sum / (n->L.size() - 1);
      if (avg_dist > max_avg_dist) {
        max_avg_dist = avg_dist;
        z_idx = i;
      }
    }

    Point z = n->L[z_idx];
    std::swap(n->L[z_idx], n->L.back());
    n->L.pop_back();

    Node *first_child = allocate_node();
    first_child->parent = n;
    first_child->L.push_back(z);
    first_child->best_point = z;
    first_child->worst_point = z;

    std::vector<Node *> new_children;
    new_children.reserve(num_children_);
    new_children.push_back(first_child);

    size_t target_children = std::min(num_children_, n->L.size() + 1);
    while (new_children.size() < target_children) {
      double best_avg_dist = -1.0;
      size_t best_z_idx = 0;

      for (size_t i = 0; i < n->L.size(); ++i) {
        double dist_sum = 0.0;
        size_t child_points = 0;
        for (auto *child : new_children) {
          for (const auto &cp : child->L) {
            dist_sum += distance(n->L[i], cp);
            child_points++;
          }
        }
        double avg_dist = (child_points > 0) ? (dist_sum / child_points) : 0.0;
        if (avg_dist > best_avg_dist) {
          best_avg_dist = avg_dist;
          best_z_idx = i;
        }
      }

      Point next_z = n->L[best_z_idx];
      std::swap(n->L[best_z_idx], n->L.back());
      n->L.pop_back();

      Node *child_prime = allocate_node();
      child_prime->parent = n;
      child_prime->L.push_back(next_z);
      child_prime->best_point = next_z;
      child_prime->worst_point = next_z;
      new_children.push_back(child_prime);
    }

    while (!n->L.empty()) {
      Point p_pt = std::move(n->L.back());
      n->L.pop_back();

      double min_d = std::numeric_limits<double>::max();
      Node *closest = nullptr;
      for (auto *child : new_children) {
        double d = expected_distance(p_pt, child);
        if (d < min_d) {
          min_d = d;
          closest = child;
        }
      }

      if (closest) {
        closest->L.push_back(p_pt);
        for (size_t i = 0; i < p_; ++i) {
          if (is_better(p_pt[i], closest->best_point[i]))
            closest->best_point[i] = p_pt[i];
          if (is_better(closest->worst_point[i], p_pt[i]))
            closest->worst_point[i] = p_pt[i];
        }
      }
    }

    n->children = std::move(new_children);
    RecalculateBounds(n);
  }

  /**
   * @brief Dynamically updates the best and worst aggregate bounding-box constraints with a new point.
   * @param n The node whose boundaries should be updated.
   * @param y The incoming bound coordinate point.
   */
  void UpdateBestWorst(Node *n, const Point &y) {
    if (n->best_point.empty()) {
      n->best_point = y;
      n->worst_point = y;
      return;
    }
    for (size_t i = 0; i < p_; ++i) {
      if (is_better(y[i], n->best_point[i]))
        n->best_point[i] = y[i];
      if (is_better(n->worst_point[i], y[i]))
        n->worst_point[i] = y[i];
    }
  }

  /**
   * @brief Helper function to perform recursive extraction for tied bounds queries.
   * @param n Current node.
   * @param z The query point.
   * @param j The objective component to check for equivalence.
   * @param out The appendable container to return the successfully located bounds.
   */
  void FindBoundsWithEqualComponentRecursive(const Node *n, const Point &z,
                                             size_t j,
                                             std::vector<Point> &out) const {
    if (n->best_point.empty() || n->worst_point.empty())
      return;

    // Reject region if it's strictly impossible to match or be dominated
    if (!is_at_least_as_good(z[j], n->worst_point[j]) ||
        !is_at_least_as_good(n->best_point[j], z[j]))
      return;

    for (size_t k = 0; k < p_; ++k) {
      if (k != j && !is_better(z[k], n->worst_point[k]))
        return;
    }

    if (n->is_leaf()) {
      for (const auto &u : n->L) {
        if (u[j] == z[j]) {
          bool strictly_better_other = true;
          for (size_t k = 0; k < p_; ++k) {
            if (k != j && !is_better(z[k], u[k])) {
              strictly_better_other = false;
              break;
            }
          }
          if (strictly_better_other)
            out.push_back(u);
        }
      }
    } else {
      for (auto *child : n->children)
        FindBoundsWithEqualComponentRecursive(child, z, j, out);
    }
  }

  /**
   * @brief Verifies if u is strictly dominated by bounds z_bar in all directions.
   * @param z_bar The query reference coordinate block.
   * @param u The candidate bound to test.
   * @return True if z_bar dominates u.
   */
  bool strictly_inside(const Point &z_bar, const Point &u) const {
    for (size_t i = 0; i < p_; ++i) {
      if (!is_better(z_bar[i], u[i]))
        return false;
    }
    return true;
  }

  /**
   * @brief Verifies if local bound a weakly-dominates b across all bounds.
   * @param a Bound coordinate 1.
   * @param b Bound coordinate 2.
   * @return True if a dominates or is equal to b.
   */
  bool covers_bound(const Point &a, const Point &b) const {
    for (size_t i = 0; i < p_; ++i) {
      if (!is_at_least_as_good(b[i], a[i]))
        return false;
    }
    return true;
  }

  /**
   * @brief Determines Euclidian distance equivalent between two abstract bounds.
   * @param a First bound.
   * @param b Second bound.
   * @return Euclidian distance.
   */
  double distance(const Point &a, const Point &b) const {
    double dist = 0.0;
    for (size_t i = 0; i < p_; ++i)
      dist += static_cast<double>((a[i] - b[i]) * (a[i] - b[i]));
    return std::sqrt(dist);
  }

  /**
   * @brief Guesses approximate centroid metrics of the sub-tree from bounding box parameters.
   * @param a Input bound being evaluated for nearest neighbor insertion.
   * @param node Evaluator node candidate.
   * @return Numeric estimate representation of bounding box center distance.
   */
  double expected_distance(const Point &a, const Node *node) const {
    double dist = 0.0;
    for (size_t i = 0; i < p_; ++i) {
      double mid_i = static_cast<double>(node->best_point[i] +
                                         node->worst_point[i]) /
                     2.0;
      double diff = static_cast<double>(a[i]) - mid_i;
      dist += diff * diff;
    }
    return std::sqrt(dist);
  }

  /**
   * @brief Helper to dump all stored coordinates recursively without removal.
   * @param n Starting structural node.
   * @param result Appended collection container.
   */
  void GetPointsRecursive(const Node *n, std::vector<Point> &result) const {
    if (n->is_leaf()) {
      result.insert(result.end(), n->L.begin(), n->L.end());
    } else {
      for (auto *child : n->children)
        GetPointsRecursive(child, result);
    }
  }
};

} // namespace local_bounds

#endif // LB_TREE_HPP
