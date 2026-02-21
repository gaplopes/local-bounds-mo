# Local Bounds Library

A C++ header-only library for maintaining local bounds in multiobjective optimization.

Based on the paper: *"On the representation of the search region in multiobjective optimization"* by Klamroth, Lacour, and Vanderpooten (European Journal of Operational Research, 2015). [DOI: 10.1016/j.ejor.2015.03.031](http://dx.doi.org/10.1016/j.ejor.2015.03.031)

This library is an **independent C++ reimplementation** of the five algorithms described in the paper. The original paper does not provide source code. Test cases are derived from the paper's worked examples (Sections 2–4) and validated against the naive brute-force algorithm for correctness.

## Features

- **Header-only** - Just include and use
- **Min/Max support** - Works for both minimization and maximization problems
- **Template-based** - Supports any numeric type (`int`, `int64_t`, `double`, etc.)
- **Multiple algorithms**:
  - `update_naive()`: Brute-force for correctness verification
  - `update_re()`: Algorithm 2 - Redundancy Elimination (RE)
  - `update_re_enhanced()`: Algorithm 3 - Redundancy Elimination (RE) Enhanced
  - `update_ra_sa()`: Algorithm 4 - Redundancy Avoidance (RA) (General Position)
  - `update_ra()`: Algorithm 5 - Redundancy Avoidance (RA) (General Case)

### Complexity (from the paper)

All algorithms share a common first step: finding the set A of bounds whose search zones contain the new point. This costs O(\|U(N)\|) with a linear scan, or O(log^p \|U(N)\| + \|A\|) with a range tree (see paper Section 5.2). **This library uses the linear scan approach.** 

The complexities below are for the **remaining steps** (candidate generation and filtering/avoidance):

| Algorithm | Remaining-step complexity | Notes |
|-----------|--------------------------|-------|
| Naive | O(p²\|A\|² + p\|A\|·\|U\|) | Brute-force: generates p\|A\| candidates, pairwise + cross filtering |
| Algorithm 2 (RE) | O(p\|A\|·(p\|A\| + \|B\|)) | Filters candidates against P ∪ B |
| Algorithm 3 (RE enhanced) | O(\|A\|²) | Prop. 5.1; reducible to O(\|A\| log \|A\|) for p ∈ {2,3}, O(\|A\| log^(p-3) \|A\| log log \|A\|) for p ≥ 4 |
| Algorithm 4 (RA, general position) | O(\|A\|) | Prop. 5.2; no filtering needed |
| Algorithm 5 (RA, general case) | O(\|N\|·\|A\|) worst case | Due to Z^k(u) sets; in practice much smaller |

Where:
- |U| is the total number of local bounds
- |A| is the number of search zones containing the new point
- |N| is the total number of points
- p is the number of objectives

### Which algorithm should I use?

Based on the empirical analysis from the paper (Section 5.3), the recommended algorithm depends on **p** (the number of objectives):

| Objectives (p) | Recommended algorithm | Rationale |
|---|---|---|
| p ≤ 5 | `update_re_enhanced()` (Algorithm 3, RE) | Better cache locality and smaller \|A\|, making the filtering step cheap |
| p = 5 | Either RE or RA | Performance is close; RE is marginally faster |
| p ≥ 6 | `update_ra_sa()` / `update_ra()` (Algorithms 4/5, RA) | Avoids expensive filtering; advantage grows rapidly with p |

The crossover is driven by **|A|** (the average number of search zones containing the new point), which grows rapidly with p. From the paper's experiments on SA instances (Section 5.3): ≈4 for p=3, ≈22 for p=4, ≈142 for p=5, ≈736 for p=6. The RE approach must filter all |A| dominated bounds each iteration, while the RA approach avoids generating redundant bounds entirely.

> **Note on Algorithm 4 vs 5:** Algorithm 4 (`update_ra_sa`) assumes *general position* (SA) — that no two distinct points share the same value in any objective. If your points may have duplicate component values, use Algorithm 5 (`update_ra`). On pathological instances with many duplicates and small objective ranges, RE may still outperform RA even at p=6.

> **Note:** Algorithms 4 and 5 require the two-argument `BoundSet` constructor (with both reference and anti-reference points). See the Quick Start section below.

## Quick Start

The simplest way is to use `update_auto()`, which selects the best algorithm based on the number of objectives `p`:

```cpp
#include "local_bounds.hpp"
using namespace local_bounds;

// Provide both reference (nadir) and anti-reference (ideal) points
// so that update_auto() can use RA algorithms when p >= 6
BoundSet<double, Objective::MINIMIZE> bounds(
    {100.0, 100.0, 100.0},  // nadir (reference)
    {  0.0,   0.0,   0.0}   // ideal (anti-reference)
);

// update_auto() picks the best algorithm for your number of objectives
bounds.update_auto(Point<double>("z1", {30.0, 70.0, 50.0}));
bounds.update_auto(Point<double>("z2", {50.0, 50.0, 40.0}));

// Check search region
bool in_region = bounds.is_in_search_region({40.0, 60.0, 45.0});

// Get current bounds
for (const auto& b : bounds.bounds()) {
    std::cout << b.to_string() << std::endl;
}
```

You can also choose a specific algorithm explicitly:

```cpp
// For p <= 5: Algorithm 3 (Enhanced RE) — single-argument constructor suffices
BoundSet<double, Objective::MINIMIZE> bounds_re({100.0, 100.0});
bounds_re.update_re_enhanced(Point<double>("z1", {30.0, 70.0}));

// For p >= 6: Algorithm 4 or 5 (RA) — requires both reference and anti-reference
std::vector<double> nadir = {100.0, 100.0, 100.0, 100.0, 100.0, 100.0};
std::vector<double> ideal = {  0.0,   0.0,   0.0,   0.0,   0.0,   0.0};
BoundSet<double, Objective::MINIMIZE> bounds_ra(nadir, ideal);
bounds_ra.update_ra(Point<double>("z1", {10.0, 50.0, 30.0, 70.0, 20.0, 60.0}));
```

## Project Structure

```
local-bounds-mo/
├── include/
│   ├── local_bounds.hpp              # Single-include entry point
│   └── local_bounds/
│       ├── types.hpp                 # Point<T>, LocalBound<T>, Objective enum
│       ├── dominance.hpp             # Dominance relation functions
│       └── bound_set.hpp             # BoundSet<T, Sense> with all 5 algorithms
├── tests/
│   ├── test_dominance.cpp            # Dominance relation tests (32 assertions)
│   └── test_bound_set.cpp            # BoundSet tests (57 assertions)
├── examples/
│   └── basic_usage.cpp               # Min/max/dominance usage examples
├── benchmark/
│   └── benchmark.cpp                 # Performance comparison across algorithms
├── CMakeLists.txt
└── CMakePresets.json                  # debug, release, dev presets
```

## Installation

### Option 1 — CMake FetchContent (recommended)

Add to your `CMakeLists.txt`:

```cmake
include(FetchContent)
FetchContent_Declare(
    local_bounds
    GIT_REPOSITORY https://github.com/gaplopes/local-bounds-mo.git
    GIT_TAG        main   # or a specific tag, e.g. v1.0.0
)
FetchContent_MakeAvailable(local_bounds)

target_link_libraries(your_target PRIVATE local_bounds)
```

### Option 2 — Add as a subdirectory

Clone or add as a git submodule, then:

```cmake
add_subdirectory(external/local-bounds-mo)
target_link_libraries(your_target PRIVATE local_bounds)
```

### Option 3 — System-wide install

```bash
cmake -B build -DBUILD_TESTS=OFF
cmake --install build --prefix /usr/local
```

Then in downstream projects:

```cmake
find_package(local_bounds REQUIRED)
target_link_libraries(your_target PRIVATE local_bounds::local_bounds)
```

### Option 4 — Copy the headers

Since this is a header-only library, you can simply copy the `include/` directory into your project and add it to your include path.

## Building

```bash
# Quick start
cmake -B build -DBUILD_TESTS=ON
cmake --build build
ctest --test-dir build

# Or use CMake presets
cmake --preset dev        # configure (tests + examples + benchmark)
cmake --build --preset dev
ctest --preset dev
```

Available presets: `debug`, `release`, `dev` (see [CMakePresets.json](CMakePresets.json)).

## API

### BoundSet

```cpp
template <typename T = double, Objective Sense = Objective::MINIMIZE>
class BoundSet {
    // For Algorithms 2/3/Naive (no defining-point tracking needed)
    explicit BoundSet(const std::vector<T>& reference_point);
    
    // For Algorithms 4/5 (requires anti-reference for dummy points)
    BoundSet(const std::vector<T>& reference_point,
             const std::vector<T>& anti_reference);
    
    void update_auto(const Point<T>& point);         // Auto-select best algorithm based on p
    void update_naive(const Point<T>& point);        // Naive
    void update_re(const Point<T>& point);              // Algorithm 2 (RE)
    void update_re_enhanced(const Point<T>& point);     // Algorithm 3 (Enhanced RE)
    void update_ra_sa(const Point<T>& point); // Algorithm 4 (RA, General Position)
    void update_ra(const Point<T>& point);    // Algorithm 5 (RA, General Case)
    
    const std::vector<LocalBound<T>>& bounds() const;
    std::size_t size() const;
    std::size_t dimensions() const;
    bool is_in_search_region(const std::vector<T>& point) const;
};
```

### Dominance Functions

```cpp
template <typename T, Objective Sense = Objective::MINIMIZE>
bool weakly_dominates(const std::vector<T>& v1, const std::vector<T>& v2);
bool strictly_dominates(const std::vector<T>& v1, const std::vector<T>& v2);
bool dominates(const std::vector<T>& v1, const std::vector<T>& v2);
bool incomparable(const std::vector<T>& v1, const std::vector<T>& v2);
```

## Contributing

Contributions are welcome!  If you have bug fixes, new problem implementations,
additional strategies, or other improvements, please fork the repository and
open a pull request.  For major changes, consider opening an issue first to
discuss the approach.

## Citation

If you use this library in your research, please cite the original paper:

```bibtex
@software{lopes2026localboundsmo,
  author       = {Lopes, Gon{\c{c}}alo},
  title = {{Local Bounds Library}: C++ header-only implementation of local bounds algorithms for multiobjective optimization},
  year         = {2026},
  url          = {https://github.com/gaplopes/local-bounds-mo}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For any questions or issues, please open an issue on the repository or contact the author at galopes@dei.uc.pt or via GitHub or LinkedIn (see profile for contact information).
