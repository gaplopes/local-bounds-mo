# Local-Bounds-MO Library

[![CI](https://github.com/gaplopes/local-bounds-mo/actions/workflows/ci.yml/badge.svg)](https://github.com/gaplopes/local-bounds-mo/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://isocpp.org/)
[![CMake 3.15+](https://img.shields.io/badge/CMake-3.15%2B-blue.svg)](https://cmake.org/)
[![Header-only](https://img.shields.io/badge/Type-Header--only-success.svg)](#)
[![Python 3.8+](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)

> Header-only C++ library for maintaining local upper/lower bounds in multiobjective optimization.

Based on two papers:
- **Paper 1:** *"On the representation of the search region in multiobjective optimization"* by Klamroth, Lacour, and Vanderpooten (EJOR, 2015). [DOI: 10.1016/j.ejor.2015.03.031](http://dx.doi.org/10.1016/j.ejor.2015.03.031)
- **Paper 2:** *"Efficient computation of the search region in multi-objective optimization"* by Dächert, Klamroth, Lacour, and Vanderpooten (EJOR, 2017). [DOI: 10.1016/j.ejor.2016.05.029](http://dx.doi.org/10.1016/j.ejor.2016.05.029)

This library is an **independent C++ reimplementation** of the algorithms described in both papers. The original papers do not provide source code. Tests cover worked examples, input requirements, tree/list invariants, visualization/API errors, and ordered-output regressions with ties and both optimization senses.

## Features

- **Header-only** - Just include and use
- **Python bindings** - Also available as a `pip`-installable Python package (via [nanobind](https://github.com/wjakob/nanobind))
- **Min/Max support** - Works for both minimization and maximization problems
- **Template-based** - Supports any numeric type (`int`, `int64_t`, `double`, etc.)
- **Redundancy Algorithms from Paper 1** (via `BoundSet`):
  - `update_naive()`: Brute-force for correctness verification
  - `update_re()`: Algorithm 2 - Redundancy Elimination (RE)
  - `update_re_enhanced()`: Algorithm 3 - Redundancy Elimination (RE) Enhanced
  - `update_ra_sa()`: Algorithm 4 - Redundancy Avoidance (RA) (General Position)
  - `update_ra()`: Algorithm 5 - Redundancy Avoidance (RA) (General Case)
- **Neighborhood-based Algorithm from Paper 2** (via `NeighborhoodBoundSet`):
  - `update()`: Algorithm 1 - Neighborhood-Based Update (Paper 2)

### Complexity

**Algorithms 2–5** (Paper 1) share a common first step: finding the set A of bounds whose search zones contain the new point. For fixed dimension p, this costs O(|U(N)|) with a linear scan, or O(log^p |U(N)| + |A|) with a range tree (see Paper 1 Section 5.2). `BoundSet` uses the linear scan approach. `BoundSetTree` uses its `LBTree` spatial index, whose pruning depends on the workload. The complexities below describe **candidate generation and filtering/avoidance for fixed p**:

| Algorithm | Class | Remaining-step complexity | Notes |
|-----------|-------|--------------------------|-------|
| Naive | `BoundSet` | O(\|A\|² + \|A\|·\|U\|) | Brute-force: generates p\|A\| candidates, pairwise + cross filtering |
| Algorithm 2 (RE) | `BoundSet` | O(\|A\|·(\|A\| + \|B\|)) | Filters candidates against P ∪ B |
| Algorithm 3 (RE enhanced) | `BoundSet` | O(\|A\|·(\|A\| + \|B\|)) | GP: O(\|A\|²) (Prop. 5.1), theoretically reducible to O(\|A\| log \|A\|) for p ∈ {2,3}, O(\|A\| log^(p-3) \|A\| log log \|A\|) for p ≥ 4 using specialized filtering not implemented here |
| Algorithm 4 (RA, GP) | `BoundSet` | O(\|A\|) | Prop. 5.2; no filtering needed |
| Algorithm 5 (RA, general) | `BoundSet` | O(\|N\|·\|A\|) worst case | Due to Z^k(u) sets; in practice much smaller |

RE also scans current bounds to build the per-coordinate B_j sets, and general-case RA scans current bounds to update defining sets. These passes are additional to the candidate costs above.

**Algorithm 1** (Paper 2) traverses the affected region through a neighborhood graph. The paper's O(|U_z̄|) bound assumes a containing bound is supplied. This implementation first scans allocated node storage and initializes a visitation array over that storage. For fixed p, an update therefore has O(C + |U_z̄|) overhead, where C is the allocated node capacity (including inactive slots).

`nonredundant_bounds()`, `nonredundant_size()`, and `get_adjacency_graph()` share global pairwise containment classification across allocated storage to handle tied-coordinate aliases correctly, retaining the first active slot for equal-coordinate bounds. The graph and update algorithm retain all active nodes, including quasi-nonredundant ones. Each call owns its classification mask; copies and subsequent updates require no cache invalidation. Standalone `get_clean_neighbors()` performs containment scans for the nodes it visits.

Classification takes O(p·C²) time in the worst case and O(C) temporary storage for the mask, where C includes inactive slots. Graph export also traverses and materializes adjacency, which can be dense. The benchmark reports exports separately from updates.

Where:
- |U(N)| is the total number of local bounds
- |A| is the number of search zones containing the new point (|A| ⊆ |U_z̄|)
- |B| is the total size of the per-coordinate B_j sets
- |U_z̄| is the number of bounds created and destroyed during one update
- |N| is the total number of points
- p is the number of objectives

### Which algorithm should I use?

Choose based on the required metadata, input assumptions, and measurements on your workload. The papers' timings do not establish a universal fastest implementation for this repository.

| Requirement | Candidate |
|---|---|
| Neighborhood relationships | `NeighborhoodBoundSet` |
| Simple bound coordinates, including ties | `BoundSet::update_re()` or `update_re_enhanced()` |
| Defining sets, with ties | `BoundSet::update_ra()` |
| Defining points and general position | `BoundSet::update_ra_sa()` |

The crossover between Algorithms 2–5 is driven by **|A|** (the average number of search zones containing the new point), which grows rapidly with p. From Paper 1's experiments: ≈4 for p=3, ≈22 for p=4, ≈142 for p=5, ≈736 for p=6.

> **Note on Algorithm 1 and General Case:** In the General Case (points sharing component values), Algorithm 1 maintains *quasi-nonredundant* bounds to preserve neighborhood graph connectivity. Use `nonredundant_bounds()` / `nonredundant_size()` to get the filtered set matching Algorithms 2–5. In the General Position case, all bounds are nonredundant, so `bounds()` and `nonredundant_bounds()` return the same set. See the article for more details.

> **Note on Algorithm 4 vs 5:** Algorithm 4 (`update_ra_sa`) assumes *general position* (SA) — that no two distinct points share the same value in any objective. If your points may have duplicate component values, use Algorithm 5 (`update_ra`).

> **Note:** Algorithm 1, 4, and 5 all require both reference and anti-reference points. See the Quick Start section below.

### Input and ownership contracts

All coordinates must be finite and have the configured nonzero dimension. Reference and anti-reference coordinates must form a strictly ordered interval. RE, naive, tree, and neighborhood updates accept the anti-reference boundary and exclude the reference boundary: `anti <= z < reference` for MINIMIZE, `reference < z <= anti` for MAXIMIZE. RA and RA-SA require strict interior points. Queries outside the configured interval return false/no bound; malformed coordinates raise an exception before mutation.

The algorithms assume a stable input set. RA-SA additionally assumes that distinct points have distinct coordinates in every objective. Use a fresh `BoundSet` when switching to RA or RA-SA: RE/naive updates discard defining metadata, and the two RA variants maintain different metadata. Unsupported transitions throw `std::logic_error`. `update_auto()` falls back to RE if RA metadata or interior-point conditions are unavailable.

`LBTree` and `BoundSetTree` are deliberately noncopyable and nonmovable: their arena pointers cannot be safely copied. `BoundSetTree::update_re_enhanced()`, `update_naive()`, and `update_auto()` are aliases of its single RE implementation, not separate implementations of Algorithms 3 or naive filtering.

## Quick Start

### Algorithm 1 — NeighborhoodBoundSet

Uses a neighborhood graph to traverse affected bounds after locating a containing bound:

```cpp
#include <iostream>
#include "local_bounds.hpp"
using namespace local_bounds;

// Requires both reference (nadir) and anti-reference (ideal) points
NeighborhoodBoundSet<double, Objective::MINIMIZE> bounds(
    {100.0, 100.0, 100.0},  // nadir (reference)
    {  0.0,   0.0,   0.0}   // ideal (anti-reference)
);

auto [updated, dominated_bounds] = bounds.update(Point<double>("z1", {30.0, 70.0, 50.0}));
bounds.update(Point<double>("z2", {50.0, 50.0, 40.0}));

// Get bounds (all, including quasi-nonredundant in NGP case)
for (const auto& b : bounds.bounds()) {
    std::cout << b.to_string() << std::endl;
}

// Or get only nonredundant bounds (matches Algorithms 2–5 output)
for (const auto& b : bounds.nonredundant_bounds()) {
    std::cout << b.to_string() << std::endl;
}

std::cout << "Total: " << bounds.size()
          << ", Nonredundant: " << bounds.nonredundant_size() << std::endl;
```

### Algorithms 2–5 — BoundSet

`update_auto()` uses a fixed heuristic: Algorithm 2 (RE) for p ≤ 5, and Algorithm 5 (RA) for p ≥ 6 when defining metadata and strict interior-point requirements are satisfied. Otherwise, it falls back to RE:

```cpp
#include <iostream>
#include "local_bounds.hpp"
using namespace local_bounds;

// Provide both reference (nadir) and anti-reference (ideal) points
// so that update_auto() can use Algorithm 5 (RA) when p >= 6.
// For p <= 5, update_auto() currently dispatches to Algorithm 2 (RE).
BoundSet<double, Objective::MINIMIZE> bounds(
    {100.0, 100.0, 100.0},  // nadir (reference)
    {  0.0,   0.0,   0.0}   // ideal (anti-reference)
);

// update_auto() uses the dimension heuristic and falls back to RE when needed
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
// Algorithm 2 (RE) — single-argument constructor suffices
BoundSet<double, Objective::MINIMIZE> bounds_re({100.0, 100.0});
bounds_re.update_re(Point<double>("z1", {30.0, 70.0}));

// Algorithm 5 (RA, general case) — requires both reference and anti-reference
std::vector<double> nadir = {100.0, 100.0, 100.0, 100.0, 100.0, 100.0};
std::vector<double> ideal = {  0.0,   0.0,   0.0,   0.0,   0.0,   0.0};
BoundSet<double, Objective::MINIMIZE> bounds_ra(nadir, ideal);
bounds_ra.update_ra(Point<double>("z1", {10.0, 50.0, 30.0, 70.0, 20.0, 60.0}));
```

### Python Quick Start

The library is also available as a Python package. Create a virtual environment and install it:

```bash
# 1. Create and activate a virtual environment (recommended)
python3 -m venv venv
source venv/bin/activate

# 2. Install the local_bounds package
pip install -e .     # Editable install (recommended for development)
# OR: pip install .  # Standard install
```

Then use it in Python:

```python
from local_bounds import BoundSet, NeighborhoodBoundSet, Point, Objective

# BoundSet with sense parameter (defaults to MINIMIZE)
bs = BoundSet([100.0, 100.0, 100.0], sense=Objective.MINIMIZE)

updated, dominated_bounds = bs.update_auto(Point("z1", [30.0, 70.0, 50.0]))
bs.update_auto(Point("z2", [50.0, 50.0, 40.0]))

print(f"Number of bounds: {bs.size()}")
for b in bs.bounds:
    print(f"  {b}")

# Check if a point is in the search region
print(bs.is_in_search_region([40.0, 60.0, 45.0]))
```

NeighborhoodBoundSet works similarly:

```python
nbs = NeighborhoodBoundSet(
    [100.0, 100.0, 100.0],  # reference (nadir)
    [0.0, 0.0, 0.0],        # anti-reference (ideal)
    sense=Objective.MINIMIZE
)

nbs.update(Point("z1", [30.0, 70.0, 50.0]))
nbs.update(Point("z2", [50.0, 50.0, 40.0]))

print(f"Total: {nbs.size()}, Nonredundant: {nbs.nonredundant_size()}")
```

> **Note:** The Python bindings use `double` coordinates. The `sense` parameter (`Objective.MINIMIZE` or `Objective.MAXIMIZE`) replaces the C++ template parameter.

## Project Structure

```
local-bounds-mo/
├── include/
│   ├── local_bounds.hpp              # Single-include entry point
│   ├── local_bounds/
│   │   ├── types.hpp                 # Point<T>, LocalBound<T>, Objective enum
│   │   ├── dominance.hpp             # Dominance relation functions
│   │   ├── bound_set.hpp             # BoundSet<T, Sense> with Algorithms 2–5
│   │   ├── bound_set_tree.hpp        # BoundSetTree<T, Sense>: indexed RE and API aliases
│   │   └── neighborhood_bound_set.hpp # NeighborhoodBoundSet<T, Sense> (Algorithm 1)
│   └── structures/
│       ├── lb_tree.hpp               # LBTree (Local Bounds Tree)
│       └── linear_list.hpp           # O(N) internal baseline list for testing
├── python/
│   ├── bindings.cpp                  # nanobind C++ bindings
│   └── local_bounds/
│       └── __init__.py               # Python package with unified API
├── visualization/
│   ├── README.md                     # Dashboard, CLI, and visualization API guide
│   ├── app.py                        # Interactive Flask dashboard
│   ├── cli.py                        # Standalone HTML reports
│   └── engine.py                     # Shared plotting and bounds extraction
├── tests/
│   ├── test_dominance.cpp            # Dominance relation tests
│   ├── test_bound_set.cpp            # BoundSet tests
│   └── test_neighborhood_bound_set.cpp # NeighborhoodBoundSet tests
├── examples/
│   └── basic_usage.cpp               # Min/max/dominance usage examples
├── benchmark/
│   ├── benchmark.cpp                 # Performance comparison across algorithms
│   ├── benchmark_lbtree.cpp          # LBTree insertion-order and containment benchmarks
│   └── benchmark_common.hpp          # Shared instance generation and timing
├── pyproject.toml                    # Python package configuration (scikit-build-core + nanobind)
├── CMakeLists.txt
└── CMakePresets.json                  # debug, release, dev presets
```

## Requirements

- **C++17** or higher
- **CMake 3.15** or higher (CMake 3.21+ required if using the provided `CMakePresets.json`)

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

### Option 5 — Python package (pip)

If you want to use the library from Python, install it directly with `pip`:

```bash
# Basic package (C++ extension & core algorithms)
pip install .

# Or install with visualization dependencies (Flask, Plotly, NetworkX)
pip install ".[vis]"
```

This compiles the C++ code with [nanobind](https://github.com/wjakob/nanobind) and installs `local_bounds` as a regular Python package. Requires a C++17 compiler and CMake 3.15+.

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

Every update method returns `UpdateResult<T>`, an alias for
`std::pair<bool, std::vector<LocalBound<T>>>`. Python returns the equivalent
`(updated, dominated_bounds)` tuple. `updated` is true exactly when the list is
nonempty; a no-op returns `{false, {}}` in C++ or `(False, [])` in Python.

The list contains owned snapshots of all existing local bounds strictly dominated
by the incoming point, with their IDs, coordinates, and available defining-point
metadata before removal. Later updates do not change these snapshots. The list
order is unspecified. Neighborhood results include affected quasi-nonredundant
bounds. Tree results use the IDs from the pre-update `bounds()` export; tree IDs
are assigned on export and do not remain stable across updates.

This changes the previous boolean return type. Unpack the result before checking
success: Python treats even `(False, [])` as truthy.

```cpp
auto [updated, dominated_bounds] = bounds.update_auto(point);
if (updated) {
    std::cout << dominated_bounds.size() << " old bounds removed\n";
}
```

```python
updated, dominated_bounds = bounds.update_auto(point)
if updated:
    print(f"{len(dominated_bounds)} old bounds removed")
```

### NeighborhoodBoundSet (Algorithm 1 — Paper 2)

```cpp
template <typename T = double, Objective Sense = Objective::MINIMIZE>
class NeighborhoodBoundSet {
    // Requires both reference and anti-reference points
    NeighborhoodBoundSet(const std::vector<T>& reference_point,
                         const std::vector<T>& anti_reference);
    
    UpdateResult<T> update(const Point<T>& point);                       // Full-storage lookup + neighborhood traversal
    
    std::vector<LocalBound<T>> bounds() const;                // All bounds (incl. quasi-nonredundant)
    std::vector<LocalBound<T>> nonredundant_bounds() const;   // Excludes quasi-nonredundant
    std::size_t size() const;                                 // Total count (incl. quasi-nonredundant)
    std::size_t nonredundant_size() const;                    // Nonredundant count only
    std::size_t dimensions() const;
    bool is_in_search_region(const std::vector<T>& point) const;
    std::optional<LocalBound<T>> find_containing_bound(const std::vector<T>& point) const;
};
```

### BoundSet (Algorithms 2–5 — Paper 1)

```cpp
template <typename T = double, Objective Sense = Objective::MINIMIZE>
class BoundSet {
    // For Algorithms 2/3/Naive (no defining-point tracking needed)
    explicit BoundSet(const std::vector<T>& reference_point);
    
    // For Algorithms 4/5 (requires anti-reference for dummy points)
    BoundSet(const std::vector<T>& reference_point,
             const std::vector<T>& anti_reference);
    
    UpdateResult<T> update_auto(const Point<T>& point);            // Dimension heuristic with RE fallback
    UpdateResult<T> update_naive(const Point<T>& point);           // Naive
    UpdateResult<T> update_re(const Point<T>& point);              // Algorithm 2 (RE)
    UpdateResult<T> update_re_enhanced(const Point<T>& point);     // Algorithm 3 (Enhanced RE)
    UpdateResult<T> update_ra_sa(const Point<T>& point);           // Algorithm 4 (RA, General Position)
    UpdateResult<T> update_ra(const Point<T>& point);              // Algorithm 5 (RA, General Case)
    
    const std::vector<LocalBound<T>>& bounds() const;
    std::size_t size() const;
    std::size_t dimensions() const;
    bool is_in_search_region(const std::vector<T>& point) const;
    std::optional<LocalBound<T>> find_containing_bound(const std::vector<T>& point) const;
};
```

### BoundSetTree (indexed RE)

> **What is LBTree?** The `LBTree` (Local Bounds Tree) is a spatial index that recursively clusters bounds and prunes multidimensional box queries. Its benefit over a linear scan depends on the workload; benchmark the intended dimensions and input distributions.

```cpp
template <typename T = double, Objective Sense = Objective::MINIMIZE>
class BoundSetTree {
    // Requires reference point (and optionally anti-reference, max_leaf_size, num_children)
    BoundSetTree(const std::vector<T>& reference_point);
    
    UpdateResult<T> update_auto(const Point<T>& point);            // Delegates to update_re_enhanced
    UpdateResult<T> update_naive(const Point<T>& point);           // Delegates to update_re
    UpdateResult<T> update_re(const Point<T>& point);              // Algorithm 2 (RE) accelerated by LBTree
    UpdateResult<T> update_re_enhanced(const Point<T>& point);     // Alias of tree RE
    
    std::vector<LocalBound<T>> bounds() const;
    std::size_t size() const;
    std::size_t dimensions() const;
    bool is_in_search_region(const std::vector<T>& point) const;
    std::optional<LocalBound<T>> find_containing_bound(const std::vector<T>& point) const;
};
```

#### Query order, ordinals, and pruning

An **ordinal** is a bound's zero-based position in the tree's current depth-first export order. `LBTree::GetAllBounds()` uses this order, and `BoundSetTree::bounds()` assigns the corresponding IDs `u0`, `u1`, `u2`, and so on. `find_containing_bound(point)` returns the first matching bound in that same order, with the same ID. These positions and IDs can change after updates; they are not permanent identifiers.

**Pruning** means skipping a subtree whose bounding box proves that none of its bounds can contain the query point in their search zones. Each node records how many bounds its subtree contains (`subtree_size`). If a skipped subtree contains 100 bounds, an ID-preserving query advances the ordinal counter by 100 without inspecting those bounds individually. This preserves the first matching bound and its ID while avoiding unnecessary dominance comparisons.

Boolean queries through `is_in_search_region()` already use bounding-box pruning. Ordinal-aware pruning allows `find_containing_bound()` to skip subtrees while still returning the correct ID. This optimization does not change the local-bound update algorithms or replace their storage.

### Dominance Functions

```cpp
template <typename T, Objective Sense = Objective::MINIMIZE>
bool weakly_dominates(const std::vector<T>& v1, const std::vector<T>& v2);
template <typename T, Objective Sense = Objective::MINIMIZE>
bool strictly_dominates(const std::vector<T>& v1, const std::vector<T>& v2);
template <typename T, Objective Sense = Objective::MINIMIZE>
bool dominates(const std::vector<T>& v1, const std::vector<T>& v2);
template <typename T, Objective Sense = Objective::MINIMIZE>
bool incomparable(const std::vector<T>& v1, const std::vector<T>& v2);
```

### Python API

The Python package exposes the same classes through a unified API where the `sense` parameter replaces the C++ template parameter:

```python
class BoundSet:
    def __init__(self, reference_point, anti_reference=None, sense=Objective.MINIMIZE): ...
    def update_auto(self, point) -> tuple[bool, list[LocalBound]]: ...
    def update_re(self, point) -> tuple[bool, list[LocalBound]]: ...
    def update_re_enhanced(self, point) -> tuple[bool, list[LocalBound]]: ...
    def update_ra_sa(self, point) -> tuple[bool, list[LocalBound]]: ...
    def update_ra(self, point) -> tuple[bool, list[LocalBound]]: ...
    def update_naive(self, point) -> tuple[bool, list[LocalBound]]: ...
    def size(self) -> int: ...
    def dimensions(self) -> int: ...
    def is_in_search_region(self, point: list[float]) -> bool: ...
    def find_containing_bound(self, point: list[float]) -> LocalBound | None: ...
    bounds: list[LocalBound]  # read-only property

class NeighborhoodBoundSet:
    def __init__(self, reference_point, anti_reference, sense=Objective.MINIMIZE): ...
    def update(self, point) -> tuple[bool, list[LocalBound]]: ...
    def size(self) -> int: ...
    def nonredundant_size(self) -> int: ...
    def dimensions(self) -> int: ...
    def is_in_search_region(self, point: list[float]) -> bool: ...
    def find_containing_bound(self, point: list[float]) -> LocalBound | None: ...
    def bounds(self) -> list[LocalBound]: ...
    def nonredundant_bounds(self) -> list[LocalBound]: ...
    def get_adjacency_graph(self, include_quasi: bool = False) -> AdjacencyGraph: ...

class BoundSetTree:
    def __init__(self, reference_point, anti_reference=None, max_leaf_size=32,
                 num_children=8, sense=Objective.MINIMIZE): ...
    def update_auto(self, point) -> tuple[bool, list[LocalBound]]: ...
    def update_re(self, point) -> tuple[bool, list[LocalBound]]: ...
    def update_re_enhanced(self, point) -> tuple[bool, list[LocalBound]]: ...
    def update_naive(self, point) -> tuple[bool, list[LocalBound]]: ...
    def size(self) -> int: ...
    def dimensions(self) -> int: ...
    def is_in_search_region(self, point: list[float]) -> bool: ...
    def find_containing_bound(self, point: list[float]) -> LocalBound | None: ...
    bounds: list[LocalBound]  # read-only property
```

## Visualization & Interactive Dashboard

The [`visualization/`](visualization/) tool provides an interactive dashboard, standalone HTML reports, and a Python plotting API for exploring local bounds in 2D and 3D. It includes spatial views, neighborhood graphs, defining-point tables, and step-by-step point insertions, with presets from both papers.

See the [visualization guide](visualization/README.md) for installation, dashboard usage, CLI options, Python examples, and troubleshooting.

## Benchmarks

The library includes a benchmarking tool (`benchmark/benchmark.cpp`) that evaluates the performance of the various algorithms by replicating the experiments from the original papers. It generates random stable sets of nondominated points and measures:

- The final number of local bounds (`|U(N)|`)
- The computation time (in milliseconds) for each algorithm
- Export time separately from update time; export timing excludes verification sorting
- The average number of bounds updated per point insertion (`|A|`)

To build and run the benchmarks, make sure to compile in Release mode:

```bash
cmake -B build -DBUILD_BENCHMARK=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build
./build/benchmark
```

The tool will run iterations for both minimization and maximization, across both *General Position* and *General Case* instance types. Default comparison results are saved to `benchmark_results.txt`, falling back to standard output if the file cannot be opened.

The `benchmark_lbtree --classification` mode compares pairwise containment scans with an LBTree index:

```bash
./build/benchmark_lbtree --classification > build/classification_benchmark.csv
./build/benchmark_lbtree --classification 1 > build/classification_smoke.csv
```

This mode uses three seeds by default (the optional argument sets the run count), both objectives, both instance types, and dimensions 3–6. It samples snapshots around 500–5,000 active bounds and reports actual active counts and allocated slots. Each timing is the median of three measurements. `scan_ms` and `index_ms` classify the same exported snapshot; `index_ms` includes building and destroying the tree and duplicate-coordinate set. `auto_count_ms` measures the library's actual `nonredundant_size()` call over allocated storage, and `graph_ms` measures complete clean graph export. This mode checks indexed and library classification masks against an independent pairwise scan. The generator can exhaust its attempt budget before reaching every target.

Its default mode measures how insertion order affects tree depth, construction time, and lookup time:

```bash
./build/benchmark_lbtree > build/lbtree.csv
./build/benchmark_lbtree 1 > build/lbtree_smoke.csv
```

It builds the same nonredundant bound snapshots in native, lexicographically sorted, and seeded shuffled orders, using the default leaf size 32 and eight children. Fixtures use seed 42, both objectives, both instance types, and dimensions 3–6, with snapshots near 500, 1,000, 3,000, and 5,000 bounds or the final available snapshot. Three shuffles (seeds 42–44) run by default; the optional argument changes this count. The `order` column includes the shuffle seed.

CSV has eight columns: `dimensions,sense,instance,bounds,order,depth,build_ms,query_ms`. Depth counts nodes from root to deepest leaf. Times are medians of three measurements; `build_ms` includes construction, insertion, and destruction. Sorting, shuffling, depth measurement, and correctness checks are untimed, so this is an order-sensitivity diagnostic rather than a comparison of total preprocessing costs. Every order must preserve the exact stored bounds and count, and match linear-scan answers for the same fixed query batch. Queries sample bound coordinates and points one integer unit inside them in every objective. Detailed query-ID checks remain in the unit tests.

## Contributing

Contributions are welcome!  If you have bug fixes, new problem implementations,
additional strategies, or other improvements, please fork the repository and
open a pull request.  For major changes, consider opening an issue first to
discuss the approach.

## Citation

If you use this library in your research, please cite the original papers and this implementation as follows:

```bibtex
@software{localboundsmo,
  author       = {Lopes, Gon{\c{c}}alo},
  title        = {{Local Bounds Library}: C++ header-only implementation of
                  local bounds algorithms for multiobjective optimization},
  year         = {2026},
  url          = {https://github.com/gaplopes/local-bounds-mo}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For any questions or issues, please open an issue on the repository or contact the author at galopes@dei.uc.pt or via GitHub or LinkedIn (see profile for contact information).
