/**
 * @file neighborhood_bound_set.hpp
 * @brief Core class for maintaining local bounds using a neighborhood graph.
 *
 * Implements Algorithm 1 from the paper "Efficient computation of the search region
 * in multi-objective optimization" (Dächert et al., 2017).
 *
 * This algorithm uses a slot map to maintain a neighborhood graph of local bounds,
 * allowing for highly efficient updates without scanning the entire bound set.
 */

#ifndef LOCAL_BOUNDS_NEIGHBORHOOD_BOUND_SET_HPP
#define LOCAL_BOUNDS_NEIGHBORHOOD_BOUND_SET_HPP

#include <algorithm>
#include <cassert>
#include <limits>
#include <optional>
#include <string>
#include <vector>
#include <stdexcept>

#include "dominance.hpp"
#include "types.hpp"

namespace local_bounds {

/**
 * @brief Maintains a set of local bounds using a neighborhood graph.
 *
 * This class implements Algorithm 1 from the paper, which uses a neighborhood
 * relation to efficiently update the search region. It supports both minimization
 * and maximization problems, and includes full handling for General Case position.
 *
 * @tparam T The numeric type for coordinates (e.g., int64_t, double).
 * @tparam Sense Objective sense (MINIMIZE or MAXIMIZE).
 */
template <typename T = double, Objective Sense = Objective::MINIMIZE>
class NeighborhoodBoundSet {
 public:
  static constexpr std::size_t npos = std::numeric_limits<std::size_t>::max();

  /**
   * @brief Constructs a NeighborhoodBoundSet with reference and anti-reference points.
   *
   * @param reference_point For MINIMIZE: the nadir point M.
   *                        For MAXIMIZE: the ideal point m.
   * @param anti_reference  For MINIMIZE: the ideal point m.
   *                        For MAXIMIZE: the nadir point M.
   */
  NeighborhoodBoundSet(const std::vector<T>& reference_point,
                       const std::vector<T>& anti_reference)
      : dimensions_(reference_point.size()),
        reference_point_(reference_point),
        anti_reference_(anti_reference) {
    if (reference_point.empty()) {
      throw std::invalid_argument("Reference point must have at least one dimension");
    }
    if (anti_reference.size() != reference_point.size()) {
      throw std::invalid_argument("Anti-reference must have same dimensions as reference");
    }
    allocate_node(LocalBound<T>::initial(reference_point, "u0"));
  }

  /**
   * @brief Updates the bound set with a new nondominated point.
   *
   * Implements Algorithm 1 from the paper.
   *
   * @param z_bar The new nondominated point.
   */
  void update(const Point<T>& z_bar) {
    if (z_bar.dimensions() != dimensions_) {
      throw std::invalid_argument("Point dimensions must match bound set dimensions");
    }

    // 1. Find an initial bound u_bar such that z_bar < u_bar
    std::size_t u_bar_idx = npos;
    for (std::size_t i = 0; i < nodes_.size(); ++i) {
      if (nodes_[i].is_active && strictly_dominates<T, Sense>(z_bar.coordinates, nodes_[i].bound.coordinates)) {
        u_bar_idx = i;
        break;
      }
    }
    
    // If no search zone contains the new point, the bound set is unaffected
    if (u_bar_idx == npos) return;

    std::size_t original_size = nodes_.size();
    std::vector<std::size_t> O;
    std::vector<bool> visited(original_size, false);
    std::vector<std::size_t> visited_indices;

    O.push_back(u_bar_idx);
    visited[u_bar_idx] = true;
    visited_indices.push_back(u_bar_idx);

    struct NewBound {
      std::size_t id;
      std::size_t parent_id;
      std::size_t i;
      std::vector<std::size_t> I;
    };
    std::vector<NewBound> B;

    // 2. Process O (Find all affected bounds and create children)
    while (!O.empty()) {
      std::size_t u_idx = O.back();
      O.pop_back();

      std::vector<std::size_t> I;
      for (std::size_t i = 0; i < dimensions_; ++i) {
        std::size_t neighbor_idx = nodes_[u_idx].neighbors[i];
        T neighbor_val;
        if (neighbor_idx != npos) {
          neighbor_val = nodes_[neighbor_idx].bound.coordinates[i];
        } else {
          neighbor_val = anti_reference_[i];
        }

        // Condition (C2) with >= for General Case handling
        bool keep;
        if constexpr (Sense == Objective::MINIMIZE) {
          keep = z_bar.coordinates[i] >= neighbor_val;
        } else {
          keep = z_bar.coordinates[i] <= neighbor_val;
        }
        if (keep) I.push_back(i);
      }

      // Create children
      for (std::size_t i : I) {
        LocalBound<T> ui_bound = nodes_[u_idx].bound;
        ui_bound.id = ui_bound.id + std::to_string(i + 1);
        ui_bound.coordinates[i] = z_bar.coordinates[i];
        
        std::size_t ui_idx = allocate_node(std::move(ui_bound));
        nodes_[u_idx].children[i] = ui_idx;
        B.push_back({ui_idx, u_idx, i, I});
      }

      // Link children with each other and with the i-neighbor
      for (std::size_t i : I) {
        std::size_t ui_idx = nodes_[u_idx].children[i];
        
        for (std::size_t j : I) {
          if (j < i) {
            std::size_t uj_idx = nodes_[u_idx].children[j];
            nodes_[ui_idx].neighbors[j] = uj_idx;
            nodes_[uj_idx].neighbors[i] = ui_idx;
          }
        }

        std::size_t i_neighbor = nodes_[u_idx].neighbors[i];
        nodes_[ui_idx].neighbors[i] = i_neighbor;
        if (i_neighbor != npos) {
          std::size_t r = get_reverse_index(i_neighbor, u_idx);
          nodes_[i_neighbor].neighbors[r] = ui_idx;
        }
      }

      // Add unvisited j-neighbors to O
      for (std::size_t j = 0; j < dimensions_; ++j) {
        if (std::find(I.begin(), I.end(), j) == I.end()) {
          std::size_t neighbor_idx = nodes_[u_idx].neighbors[j];
          assert(neighbor_idx != npos); // Guaranteed by Proposition 4.3
          if (!visited[neighbor_idx]) {
            O.push_back(neighbor_idx);
            visited[neighbor_idx] = true;
            visited_indices.push_back(neighbor_idx);
          }
        }
      }
    }

    // 3. Link remaining neighbors
    for (const auto& new_bound : B) {
      std::size_t ui = new_bound.id;
      std::size_t u = new_bound.parent_id;
      std::size_t i = new_bound.i;

      for (std::size_t j = 0; j < dimensions_; ++j) {
        if (std::find(new_bound.I.begin(), new_bound.I.end(), j) == new_bound.I.end()) {
          std::size_t u_hat = nodes_[u].neighbors[j];
          assert(u_hat != npos);
          std::size_t r = get_reverse_index(u_hat, u);
          std::size_t t = j;

          while (nodes_[ui].neighbors[j] == npos) {
            if (r != i) t = i;
            std::size_t u_hat_t = nodes_[u_hat].children[t];
            if (u_hat_t != npos) {
              nodes_[ui].neighbors[j] = u_hat_t;
              nodes_[u_hat_t].neighbors[r] = ui;
            } else {
              std::size_t next_u_hat = nodes_[u_hat].neighbors[t];
              assert(next_u_hat != npos);
              r = get_reverse_index(next_u_hat, u_hat);
              u_hat = next_u_hat;
            }
          }
        }
      }
    }

    // 4. Cleanup affected bounds
    for (std::size_t idx : visited_indices) {
      free_node(idx);
    }
  }

  /**
   * @brief Returns all active local bounds (including quasi-nonredundant ones).
   *
   * Note: In the General Case position, this may include quasi-nonredundant bounds
   * whose search zones are contained within another bound's search zone.
   * These are kept to preserve neighborhood graph connectivity.
   * Use nonredundant_bounds() to get only bounds with distinct search zones.
   */
  [[nodiscard]] std::vector<LocalBound<T>> bounds() const {
    std::vector<LocalBound<T>> active_bounds;
    active_bounds.reserve(nodes_.size() - free_list_.size());
    for (const auto& node : nodes_) {
      if (node.is_active) {
        active_bounds.push_back(node.bound);
      }
    }
    return active_bounds;
  }

  /**
   * @brief Returns only nonredundant local bounds (excludes quasi-nonredundant ones).
   *
   * In the General Case position, Algorithm 1 maintains quasi-nonredundant
   * bounds to preserve neighborhood connectivity (see paper, Section 4.3). A bound u
   * is quasi-nonredundant if one of its i-neighbors weakly dominates it, meaning its
   * search zone is contained within the neighbor's search zone.
   *
   * This method filters those out, returning only bounds with genuinely distinct
   * search zones — matching the output of Algorithms 2/3/5.
   */
  [[nodiscard]] std::vector<LocalBound<T>> nonredundant_bounds() const {
    std::vector<LocalBound<T>> result;
    result.reserve(nodes_.size() - free_list_.size());
    for (std::size_t idx = 0; idx < nodes_.size(); ++idx) {
      if (nodes_[idx].is_active && !is_quasi_nonredundant(idx)) {
        result.push_back(nodes_[idx].bound);
      }
    }
    return result;
  }

  /**
   * @brief Returns the total number of active local bounds (including quasi-nonredundant).
   */
  [[nodiscard]] std::size_t size() const {
    return nodes_.size() - free_list_.size();
  }

  /**
   * @brief Returns the number of nonredundant local bounds (excludes quasi-nonredundant).
   *
   * This count should match the output size of Algorithms 2/3/5 on the same input.
   */
  [[nodiscard]] std::size_t nonredundant_size() const {
    std::size_t count = 0;
    for (std::size_t idx = 0; idx < nodes_.size(); ++idx) {
      if (nodes_[idx].is_active && !is_quasi_nonredundant(idx)) {
        ++count;
      }
    }
    return count;
  }

  /**
   * @brief Returns the dimensionality of the objective space.
   */
  [[nodiscard]] std::size_t dimensions() const {
    return dimensions_;
  }

  /**
   * @brief Checks if a point is in the search region.
   */
  [[nodiscard]] bool is_in_search_region(const std::vector<T>& point) const {
    if (point.size() != dimensions_) {
      throw std::invalid_argument("Point dimensions must match bound set dimensions");
    }
    for (const auto& node : nodes_) {
      if (node.is_active && strictly_dominates<T, Sense>(point, node.bound.coordinates)) {
        return true;
      }
    }
    return false;
  }

  /**
   * @brief Finds the bound whose search zone contains the given point.
   */
  [[nodiscard]] std::optional<LocalBound<T>> find_containing_bound(
      const std::vector<T>& point) const {
    for (const auto& node : nodes_) {
      if (node.is_active && strictly_dominates<T, Sense>(point, node.bound.coordinates)) {
        return node.bound;
      }
    }
    return std::nullopt;
  }

 private:
  struct Node {
    LocalBound<T> bound;
    std::vector<std::size_t> neighbors;
    std::vector<std::size_t> children;
    bool is_active = false;
  };

  std::size_t dimensions_;
  std::vector<T> reference_point_;
  std::vector<T> anti_reference_;

  std::vector<Node> nodes_;
  std::vector<std::size_t> free_list_;

  std::size_t allocate_node(LocalBound<T> bound) {
    if (!free_list_.empty()) {
      std::size_t idx = free_list_.back();
      free_list_.pop_back();
      nodes_[idx].bound = std::move(bound);
      nodes_[idx].is_active = true;
      std::fill(nodes_[idx].children.begin(), nodes_[idx].children.end(), npos);
      std::fill(nodes_[idx].neighbors.begin(), nodes_[idx].neighbors.end(), npos);
      return idx;
    }
    Node node;
    node.bound = std::move(bound);
    node.neighbors.resize(dimensions_, npos);
    node.children.resize(dimensions_, npos);
    node.is_active = true;
    nodes_.push_back(std::move(node));
    return nodes_.size() - 1;
  }

  void free_node(std::size_t idx) {
    nodes_[idx].is_active = false;
    std::fill(nodes_[idx].children.begin(), nodes_[idx].children.end(), npos);
    std::fill(nodes_[idx].neighbors.begin(), nodes_[idx].neighbors.end(), npos);
    free_list_.push_back(idx);
  }

  /**
   * @brief Checks whether a bound is quasi-nonredundant.
   *
   * A bound u is quasi-nonredundant if its search zone C(u) is contained
   * within the search zone of one of its i-neighbors: C(u) ⊆ C(ν_i(u)).
   * 
   * Example in the MINIMIZE case:
   * For MINIMIZE: C(u) ⊆ C(ν_i(u)) iff u_j ≤ ν_i(u)_j for all j,
   * i.e., u weakly dominates ν_i(u) in the optimization sense.
   * This happens in the NGP case when z̄_i = ν_i(u)_i (equality in C2).
   *
   * Note:
   * When two neighbors have identical coordinates (mutual weak domination),
   * their search zones are identical. Only one should be marked redundant;
   * we use the node index as a tiebreaker (higher index is flagged).
   *
   * This is an O(p²) check using only the neighborhood graph pointers.
   */
  bool is_quasi_nonredundant(std::size_t idx) const {
    const auto& u = nodes_[idx];
    for (std::size_t i = 0; i < dimensions_; ++i) {
      std::size_t nb_idx = u.neighbors[i];
      if (nb_idx == npos || !nodes_[nb_idx].is_active) continue;
      if (weakly_dominates<T, Sense>(u.bound.coordinates,
                                     nodes_[nb_idx].bound.coordinates)) {
        // Check if the relation is mutual (identical search zones)
        if (weakly_dominates<T, Sense>(nodes_[nb_idx].bound.coordinates,
                                       u.bound.coordinates)) {
          // Both have identical search zones — only flag the higher-indexed one
          if (idx > nb_idx) return true;
        } else {
          // Strict containment: C(u) ⊊ C(ν_i(u)) — u is genuinely redundant
          return true;
        }
      }
    }
    return false;
  }

  std::size_t get_reverse_index(std::size_t u_hat, std::size_t u) const {
    for (std::size_t r = 0; r < dimensions_; ++r) {
      if (nodes_[u_hat].neighbors[r] == u) return r;
    }
    throw std::runtime_error("Reverse index not found. Graph is corrupted.");
  }
};

}  // namespace local_bounds

#endif  // LOCAL_BOUNDS_NEIGHBORHOOD_BOUND_SET_HPP
