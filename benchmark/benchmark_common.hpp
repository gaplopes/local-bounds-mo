#ifndef LOCAL_BOUNDS_BENCHMARK_COMMON_HPP
#define LOCAL_BOUNDS_BENCHMARK_COMMON_HPP

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <random>
#include <string>
#include <vector>

#include "local_bounds/dominance.hpp"
#include "local_bounds/types.hpp"

/**
 * Instance type controlling whether points satisfy the
 * General Position assumption (SA) from the paper.
 */
enum class InstanceType {
  GENERAL_POSITION,  // SA: no two points share a value in any dimension
  GENERAL_CASE       // General: points may share component values
};

inline std::string instance_type_str(InstanceType type) {
  switch (type) {
    case InstanceType::GENERAL_POSITION:
      return "SA";
    case InstanceType::GENERAL_CASE:
      return "General";
  }
  return "Unknown";
}

/**
 * Generate a random stable set of points (no point dominates another).
 */
template <local_bounds::Objective Sense>
std::vector<local_bounds::Point<int64_t>> generate_stable_set(
    std::size_t num_points,
    std::size_t dimensions,
    int64_t K,
    InstanceType instance_type,
    unsigned seed) {
  std::mt19937_64 rng(seed);
  std::uniform_int_distribution<int64_t> dist(1, K);

  std::vector<local_bounds::Point<int64_t>> points;
  points.reserve(num_points);

  std::size_t attempts = 0;
  const std::size_t max_attempts = num_points * 100;

  while (points.size() < num_points && attempts < max_attempts) {
    ++attempts;

    // Generate random point
    std::vector<int64_t> coords(dimensions);
    for (std::size_t j = 0; j < dimensions; ++j) {
      coords[j] = dist(rng);
    }

    // Check if dominated by or dominates any existing point
    bool is_valid = true;
    for (const auto& existing : points) {
      if (local_bounds::weakly_dominates<int64_t, Sense>(existing.coordinates, coords) ||
          local_bounds::weakly_dominates<int64_t, Sense>(coords, existing.coordinates)) {
        is_valid = false;
        break;
      }
    }

    // Under general position (SA), also reject if any component value
    // is shared with an existing point in any dimension
    if (is_valid && instance_type == InstanceType::GENERAL_POSITION) {
      for (const auto& existing : points) {
        for (std::size_t j = 0; j < dimensions; ++j) {
          if (coords[j] == existing.coordinates[j]) {
            is_valid = false;
            break;
          }
        }
        if (!is_valid) break;
      }
    }

    if (is_valid) {
      points.emplace_back("z" + std::to_string(points.size() + 1), std::move(coords));
    }
  }

  // Randomly reorder the points
  std::shuffle(points.begin(), points.end(), rng);

  return points;
}

template <typename Operation>
double median_time_ms(Operation operation) {
  std::vector<double> samples;
  volatile std::size_t checksum = 0;
  for (int i = 0; i < 3; ++i) {
    const auto start = std::chrono::steady_clock::now();
    checksum += operation();
    samples.push_back(std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count());
  }
  std::sort(samples.begin(), samples.end());
  return samples[1];
}

#endif // LOCAL_BOUNDS_BENCHMARK_COMMON_HPP
