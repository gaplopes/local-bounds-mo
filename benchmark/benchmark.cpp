/**
 * @file benchmark.cpp
 * @brief Benchmark comparing algorithm performance.
 *
 * Replicates the experiments from the paper:
 * "On the representation of the search region in multiobjective optimization"
 *
 * Generates random stable sets of points and measures:
 * - Number of local upper bounds
 * - Computation time for each algorithm
 * - Average |A| (dominated bounds per update)
 */

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <vector>

#include "local_bounds.hpp"

using namespace local_bounds;

/**
 * Instance type controlling whether points satisfy the
 * General Position assumption (SA) from the paper.
 */
enum class InstanceType {
  GENERAL_POSITION,  // SA: no two points share a value in any dimension
  GENERAL_CASE       // General: points may share component values
};

std::string instance_type_str(InstanceType type) {
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
template <Objective Sense>
std::vector<Point<int64_t>> generate_stable_set(
    std::size_t num_points,
    std::size_t dimensions,
    int64_t K,
    InstanceType instance_type,
    unsigned seed) {
  std::mt19937_64 rng(seed);
  std::uniform_int_distribution<int64_t> dist(1, K);

  std::vector<Point<int64_t>> points;
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
      if (dominates<int64_t, Sense>(existing.coordinates, coords) ||
          dominates<int64_t, Sense>(coords, existing.coordinates)) {
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
      points.emplace_back("z" + std::to_string(points.size() + 1), coords);
    }
  }

  // Randomly reorder the points
  std::shuffle(points.begin(), points.end(), rng);

  return points;
}

struct BenchmarkResult {
  std::string algorithm;
  std::string sense;
  std::string instance_type;
  std::size_t dimensions;
  std::size_t num_points;
  std::size_t final_bounds;
  double time_ms;
  double avg_A;
};

std::vector<std::vector<int64_t>> canonicalize_bounds(
    const std::vector<LocalBound<int64_t>>& bounds) {
  std::vector<std::vector<int64_t>> out;
  out.reserve(bounds.size());
  for (const auto& b : bounds) {
    out.push_back(b.coordinates);
  }
  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

std::string format_bound(const std::vector<int64_t>& v) {
  std::ostringstream oss;
  oss << "(";
  for (std::size_t i = 0; i < v.size(); ++i) {
    if (i) oss << ",";
    oss << v[i];
  }
  oss << ")";
  return oss.str();
}

void print_bound_set_delta(
    std::ostream& out,
    const std::vector<std::vector<int64_t>>& a,
    const std::vector<std::vector<int64_t>>& b,
    const std::string& a_name,
    const std::string& b_name,
    std::size_t max_to_print = 15) {
  std::vector<std::vector<int64_t>> only_a;
  std::vector<std::vector<int64_t>> only_b;

  std::set_difference(
      a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(only_a));
  std::set_difference(
      b.begin(), b.end(), a.begin(), a.end(), std::back_inserter(only_b));

  out << "    |" << a_name << "|=" << a.size()
      << " |" << b_name << "|=" << b.size() << "\n";
  out << "    In " << a_name << " and not in " << b_name << ": "
      << only_a.size() << "\n";
  for (std::size_t i = 0; i < std::min(max_to_print, only_a.size()); ++i) {
    out << "      " << format_bound(only_a[i]) << "\n";
  }
  if (only_a.size() > max_to_print) {
    out << "      ... (" << (only_a.size() - max_to_print) << " more)\n";
  }

  out << "    In " << b_name << " and not in " << a_name << ": "
      << only_b.size() << "\n";
  for (std::size_t i = 0; i < std::min(max_to_print, only_b.size()); ++i) {
    out << "      " << format_bound(only_b[i]) << "\n";
  }
  if (only_b.size() > max_to_print) {
    out << "      ... (" << (only_b.size() - max_to_print) << " more)\n";
  }
}

template <
    Objective Sense,
    void (BoundSet<int64_t, Sense>::*UpdateA)(const Point<int64_t>&),
    void (BoundSet<int64_t, Sense>::*UpdateB)(const Point<int64_t>&)>
bool compare_iteration_by_iteration(
    std::ostream& out,
    const std::string& a_name,
    const std::string& b_name,
    const std::vector<Point<int64_t>>& points,
    std::size_t dimensions,
    int64_t ref_val,
    int64_t anti_ref_val,
    std::size_t max_points_to_check = std::numeric_limits<std::size_t>::max()) {
  BoundSet<int64_t, Sense> a(std::vector<int64_t>(dimensions, ref_val),
                             std::vector<int64_t>(dimensions, anti_ref_val));
  BoundSet<int64_t, Sense> b(std::vector<int64_t>(dimensions, ref_val),
                             std::vector<int64_t>(dimensions, anti_ref_val));

  const std::size_t limit = std::min(points.size(), max_points_to_check);
  for (std::size_t i = 0; i < limit; ++i) {
    const auto& point = points[i];
    (a.*UpdateA)(point);
    (b.*UpdateB)(point);

    auto ca = canonicalize_bounds(a.bounds());
    auto cb = canonicalize_bounds(b.bounds());

    if (ca != cb) {
      out << "  [DIFF] " << a_name << " vs " << b_name
          << " at iteration " << i
          << " after inserting " << point.id << "="
          << format_bound(point.coordinates) << "\n";
      print_bound_set_delta(out, ca, cb, a_name, b_name);
      return false;
    }
  }

  out << "  [OK] " << a_name << " == " << b_name
      << " for " << limit << " iterations\n";
  return true;
}

template <Objective Sense,
          void (BoundSet<int64_t, Sense>::*UpdateMethod)(const Point<int64_t>&)>
BenchmarkResult run_benchmark(
    const std::string& name,
    const std::string& sense_str,
    const std::string& inst_type_str,
    const std::vector<Point<int64_t>>& points,
    std::size_t dimensions,
    int64_t ref_val,
    int64_t anti_ref_val) {
  BoundSet<int64_t, Sense> bounds(std::vector<int64_t>(dimensions, ref_val),
                                  std::vector<int64_t>(dimensions, anti_ref_val));

  double total_A = 0;

  auto start = std::chrono::high_resolution_clock::now();

  for (const auto& point : points) {
    // Count |A| before update
    std::size_t A_size = 0;
    for (const auto& b : bounds.bounds()) {
      if (strictly_dominates<int64_t, Sense>(point.coordinates, b.coordinates)) {
        ++A_size;
      }
    }
    total_A += A_size;

    (bounds.*UpdateMethod)(point);
  }

  auto end = std::chrono::high_resolution_clock::now();
  double time_ms = std::chrono::duration<double, std::milli>(end - start).count();

  return {
      name,
      sense_str,
      inst_type_str,
      dimensions,
      points.size(),
      bounds.size(),
      time_ms,
      points.empty() ? 0 : total_A / points.size()};
}

void print_header(std::ostream& out) {
  out << std::left
      << std::setw(24) << "Algorithm"
      << std::setw(6) << "Sense"
      << std::setw(10) << "Instance"
      << std::setw(5) << "p"
      << std::setw(10) << "|N|"
      << std::setw(12) << "|U(N)|"
      << std::setw(14) << "Time (ms)"
      << std::setw(10) << "Avg |A|"
      << std::endl;
  out << std::string(92, '-') << std::endl;
}

void print_result(std::ostream& out, const BenchmarkResult& r) {
  out << std::left
      << std::setw(24) << r.algorithm
      << std::setw(6) << r.sense
      << std::setw(10) << r.instance_type
      << std::setw(5) << r.dimensions
      << std::setw(10) << r.num_points
      << std::setw(12) << r.final_bounds
      << std::setw(14) << std::fixed << std::setprecision(2) << r.time_ms
      << std::setw(10) << std::fixed << std::setprecision(2) << r.avg_A
      << std::endl;
}

template <Objective Sense>
void run_benchmarks_for_sense(
    std::ostream& out,
    const std::string& sense_str,
    const std::vector<std::pair<std::size_t, std::size_t>>& configs,
    InstanceType instance_type,
    int64_t K,
    int64_t ref_val,
    int64_t anti_ref_val,
    unsigned num_runs) {
  const std::string inst_str = instance_type_str(instance_type);

  for (const auto& [dims, target_points] : configs) {
    out << "\n--- " << sense_str << " p = " << dims
        << ", target |N| = " << target_points << " ---\n";

    for (unsigned run = 0; run < num_runs; ++run) {
      unsigned seed = 42 + run;
      auto points = generate_stable_set<Sense>(
          target_points, dims, K, instance_type, seed);

      out << "Generated " << points.size() << " points\n";

      // Algorithm 2: Basic Redundancy Elimination (General Case)
      auto alg2 = run_benchmark<Sense, &BoundSet<int64_t, Sense>::update_re>(
          "Algorithm 2", sense_str, inst_str, points, dims, ref_val, anti_ref_val);
      print_result(out, alg2);

      // Algorithm 3: Enhanced Filtering (General Case)
      auto alg3 = run_benchmark<Sense, &BoundSet<int64_t, Sense>::update_re_enhanced>(
          "Algorithm 3", sense_str, inst_str, points, dims, ref_val, anti_ref_val);
      print_result(out, alg3);
      
      // Algorithm 4: Redundancy Avoidance (General Position only)
      BenchmarkResult alg4{};
      if (instance_type == InstanceType::GENERAL_POSITION) {
        alg4 = run_benchmark<Sense, &BoundSet<int64_t, Sense>::update_ra_sa>(
            "Algorithm 4", sense_str, inst_str, points, dims, ref_val, anti_ref_val);
        print_result(out, alg4);
      }

      // Algorithm 5: Redundancy Avoidance (General Case)
      auto alg5 = run_benchmark<Sense, &BoundSet<int64_t, Sense>::update_ra>(
          "Algorithm 5", sense_str, inst_str, points, dims, ref_val, anti_ref_val);
      print_result(out, alg5);

      // Naive Algorithm: Brute-force generation and filtering of all candidate bounds (General Case)
      auto naive = run_benchmark<Sense, &BoundSet<int64_t, Sense>::update_naive>(
          "Naive", sense_str, inst_str, points, dims, ref_val, anti_ref_val);
      print_result(out, naive);

      // Iteration-by-iteration comparison to identify exact differences.
      // To keep runtime practical, limit detailed checks to a subset.
      if (naive.final_bounds != alg2.final_bounds) {
        out << "  [WARNING] Final bound counts differ between Naive and Algorithm 2\n";
        const std::size_t comparison_limit = std::max(naive.final_bounds, alg2.final_bounds);
        out << "  Iterative set comparison (first " << comparison_limit << " points):\n";
        compare_iteration_by_iteration<
            Sense,
            &BoundSet<int64_t, Sense>::update_naive,
            &BoundSet<int64_t, Sense>::update_re>(
            out, "Naive", "Algorithm 2", points, dims, ref_val, anti_ref_val, comparison_limit);
      }
      if (alg2.final_bounds != alg3.final_bounds) {
        out << "  [WARNING] Final bound counts differ between Algorithm 2 and Algorithm 3\n";
        const std::size_t comparison_limit = std::max(alg2.final_bounds, alg3.final_bounds);
        out << "  Iterative set comparison (first " << comparison_limit << " points):\n";
        compare_iteration_by_iteration<
            Sense,
            &BoundSet<int64_t, Sense>::update_re,
            &BoundSet<int64_t, Sense>::update_re_enhanced>(
            out, "Algorithm 2", "Algorithm 3", points, dims, ref_val, anti_ref_val, comparison_limit);
      }
      if (instance_type == InstanceType::GENERAL_POSITION) {
        if (alg2.final_bounds != alg4.final_bounds) {
          out << "  [WARNING] Final bound counts differ between Algorithm 2 and Algorithm 4\n";
          const std::size_t comparison_limit = std::max(alg2.final_bounds, alg4.final_bounds);
          out << "  Iterative set comparison (first " << comparison_limit << " points):\n";
          compare_iteration_by_iteration<
              Sense,
              &BoundSet<int64_t, Sense>::update_re,
              &BoundSet<int64_t, Sense>::update_ra_sa>(
              out, "Algorithm 2", "Algorithm 4", points, dims, ref_val, anti_ref_val, comparison_limit);
        }
      }
      if (alg2.final_bounds != alg5.final_bounds) {
        out << "  [WARNING] Final bound counts differ between Algorithm 2 and Algorithm 5\n";
        const std::size_t comparison_limit = std::max(alg2.final_bounds, alg5.final_bounds);
        out << "  Iterative set comparison (first " << comparison_limit << " points):\n";
        compare_iteration_by_iteration<
            Sense,
            &BoundSet<int64_t, Sense>::update_re,
            &BoundSet<int64_t, Sense>::update_ra>(
            out, "Algorithm 2", "Algorithm 5", points, dims, ref_val, anti_ref_val, comparison_limit);
      }
    }
  }
}

int main(int /*argc*/, char* /*argv*/[]) {
  std::cout << "=============================================\n";
  std::cout << "  Local Bounds Algorithm Benchmark\n";
  std::cout << "=============================================\n\n";

  // Open results file
  std::ofstream results_file("benchmark_results.txt");
  if (!results_file) {
    std::cerr << "Warning: Could not open results file for writing\n";
  }

  // Write to both console and file
  auto write_both = [&](const std::string& msg) {
    std::cout << msg;
    if (results_file) results_file << msg;
  };

  write_both("=============================================\n");
  write_both("  Local Bounds Algorithm Benchmark\n");
  write_both("  Replicating experiments from the paper:\n");
  write_both("  'On the representation of the search region\n");
  write_both("  in multiobjective optimization'\n");
  write_both("=============================================\n\n");

  // Configuration (reduced sizes for reasonable runtime)
  std::vector<std::pair<std::size_t, std::size_t>> configs = {
      {3, 5000},  // (dimensions, num_points)
      {4, 2500},
      {5, 1000},
      {6, 500}};

  const int64_t K = 10000;
  const unsigned num_runs = 10;
  const std::vector<InstanceType> instance_types = {
      InstanceType::GENERAL_POSITION,
      InstanceType::GENERAL_CASE};

  // Print header
  print_header(std::cout);
  if (results_file) print_header(results_file);

  for (const auto& instance_type : instance_types) {
    // Run benchmarks for both MINIMIZE and MAXIMIZE senses
    // Run MINIMIZATION benchmarks
    // ref_val = M = K+1 (nadir), anti_ref_val = m = 0 (ideal)
    write_both("\n========== MINIMIZATION ==========\n");
    run_benchmarks_for_sense<Objective::MINIMIZE>(
        std::cout, "MIN", configs, instance_type, K, K + 1, 0, num_runs);
    if (results_file) {
      run_benchmarks_for_sense<Objective::MINIMIZE>(
          results_file, "MIN", configs, instance_type, K, K + 1, 0, num_runs);
    }

    // Run MAXIMIZATION benchmarks
    // ref_val = m = 0 (ideal), anti_ref_val = M = K+1 (nadir)
    write_both("\n========== MAXIMIZATION ==========\n");
    run_benchmarks_for_sense<Objective::MAXIMIZE>(
        std::cout, "MAX", configs, instance_type, K, 0, K + 1, num_runs);
    if (results_file) {
      run_benchmarks_for_sense<Objective::MAXIMIZE>(
          results_file, "MAX", configs, instance_type, K, 0, K + 1, num_runs);
    }
  }

  write_both("\n=============================================\n");
  write_both("  Benchmark Complete\n");
  write_both("  Results saved to: benchmark_results.txt\n");
  write_both("=============================================\n");

  return 0;
}
