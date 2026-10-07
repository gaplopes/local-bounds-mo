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

**Algorithms 2–5** (Paper 1) share a common first step: finding the set A of bounds whose search zones contain the new point. This costs O(|U(N)|) with a linear scan, or O(log^p |U(N)| + |A|) with a range tree (see Paper 1 Section 5.2). **This library uses the linear scan approach.** The complexities below are for the **remaining steps** (candidate generation and filtering/avoidance):

| Algorithm | Class | Remaining-step complexity | Notes |
|-----------|-------|--------------------------|-------|
| Naive | `BoundSet` | O(p²\|A\|² + p\|A\|·\|U\|) | Brute-force: generates p\|A\| candidates, pairwise + cross filtering |
| Algorithm 2 (RE) | `BoundSet` | O(p\|A\|·(p\|A\| + \|B\|)) | Filters candidates against P ∪ B |
| Algorithm 3 (RE enhanced) | `BoundSet` | O(\|A\|²) | Prop. 5.1; reducible to O(\|A\| log \|A\|) for p ∈ {2,3}, O(\|A\| log^(p-3) \|A\| log log \|A\|) for p ≥ 4 |
| Algorithm 4 (RA, GP) | `BoundSet` | O(\|A\|) | Prop. 5.2; no filtering needed |
| Algorithm 5 (RA, general) | `BoundSet` | O(\|N\|·\|A\|) worst case | Due to Z^k(u) sets; in practice much smaller |

**Algorithm 1** (Paper 2) traverses the affected region through a neighborhood graph. The paper's O(|U_z̄|) bound assumes a containing bound is supplied. This implementation first scans allocated node storage and initializes a visitation array over that storage. For fixed p, an update therefore has O(C + |U_z̄|) overhead, where C is the allocated node capacity (including inactive slots).

`nonredundant_bounds()` and `nonredundant_size()` check global containment to handle tied-coordinate aliases correctly. These exports take O(p·C²) in the worst case. `get_adjacency_graph()` classifies each active node once per export and reuses the results while contracting redundant nodes. Its cost also includes graph traversal and materializing the complete adjacency, which can be dense. Exports occur outside the benchmark's timed update section.

Where:
- |U(N)| is the total number of local bounds
- |A| is the number of search zones containing the new point (|A| ⊆ |U_z̄|)
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
#include "local_bounds.hpp"
using namespace local_bounds;

// Requires both reference (nadir) and anti-reference (ideal) points
NeighborhoodBoundSet<double, Objective::MINIMIZE> bounds(
    {100.0, 100.0, 100.0},  // nadir (reference)
    {  0.0,   0.0,   0.0}   // ideal (anti-reference)
);

bounds.update(Point<double>("z1", {30.0, 70.0, 50.0}));
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

The simplest way is to use `update_auto()`, which selects the best algorithm based on the number of objectives `p`:

```cpp
#include "local_bounds.hpp"
using namespace local_bounds;

// Provide both reference (nadir) and anti-reference (ideal) points
// so that update_auto() can use RA algorithms when p >= 6.
// For p <= 5, update_auto() currently dispatches to Algorithm 2 (RE).
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
// For p <= 5: Algorithm 2 (RE) — single-argument constructor suffices
BoundSet<double, Objective::MINIMIZE> bounds_re({100.0, 100.0});
bounds_re.update_re(Point<double>("z1", {30.0, 70.0}));

// For p >= 6: Algorithm 4 or 5 (RA) — requires both reference and anti-reference
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

bs.update_auto(Point("z1", [30.0, 70.0, 50.0]))
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
├── tests/
│   ├── test_dominance.cpp            # Dominance relation tests
│   ├── test_bound_set.cpp            # BoundSet tests
│   └── test_neighborhood_bound_set.cpp # NeighborhoodBoundSet tests
├── examples/
│   └── basic_usage.cpp               # Min/max/dominance usage examples
├── benchmark/
│   └── benchmark.cpp                 # Performance comparison across algorithms
├── pyproject.toml                    # Python package configuration (scikit-build-core + nanobind)
├── CMakeLists.txt
└── CMakePresets.json                  # debug, release, dev presets
```

## Requirements

- **C++17** or higher
- **CMake 3.14** or higher (CMake 3.21+ required if using the provided `CMakePresets.json`)

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

### NeighborhoodBoundSet (Algorithm 1 — Paper 2)

```cpp
template <typename T = double, Objective Sense = Objective::MINIMIZE>
class NeighborhoodBoundSet {
    // Requires both reference and anti-reference points
    NeighborhoodBoundSet(const std::vector<T>& reference_point,
                         const std::vector<T>& anti_reference);
    
    bool update(const Point<T>& point);                       // Full-storage lookup + neighborhood traversal
    
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
    
    bool update_auto(const Point<T>& point);            // Auto-select best algorithm based on p
    bool update_naive(const Point<T>& point);           // Naive
    bool update_re(const Point<T>& point);              // Algorithm 2 (RE)
    bool update_re_enhanced(const Point<T>& point);     // Algorithm 3 (Enhanced RE)
    bool update_ra_sa(const Point<T>& point);           // Algorithm 4 (RA, General Position)
    bool update_ra(const Point<T>& point);              // Algorithm 5 (RA, General Case)
    
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
    
    bool update_auto(const Point<T>& point);            // Delegates to update_re_enhanced
    bool update_naive(const Point<T>& point);           // Delegates to update_re
    bool update_re(const Point<T>& point);              // Algorithm 2 (RE) accelerated by LBTree
    bool update_re_enhanced(const Point<T>& point);     // Alias of tree RE
    
    std::vector<LocalBound<T>> bounds() const;
    std::size_t size() const;
    std::size_t dimensions() const;
    bool is_in_search_region(const std::vector<T>& point) const;
    std::optional<LocalBound<T>> find_containing_bound(const std::vector<T>& point) const;
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

### Python API

The Python package exposes the same classes through a unified API where the `sense` parameter replaces the C++ template parameter:

```python
class BoundSet:
    def __init__(self, reference_point, anti_reference=None, sense=Objective.MINIMIZE): ...
    def update_auto(self, point): ...
    def update_re(self, point): ...
    def update_re_enhanced(self, point): ...
    def update_ra_sa(self, point): ...
    def update_ra(self, point): ...
    def update_naive(self, point): ...
    def size(self) -> int: ...
    def dimensions(self) -> int: ...
    def is_in_search_region(self, point: list[float]) -> bool: ...
    def find_containing_bound(self, point: list[float]) -> LocalBound | None: ...
    bounds: list[LocalBound]  # read-only property

class NeighborhoodBoundSet:
    def __init__(self, reference_point, anti_reference, sense=Objective.MINIMIZE): ...
    def update(self, point): ...
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
    def update_auto(self, point): ...
    def update_re(self, point): ...
    def update_re_enhanced(self, point): ...
    def update_naive(self, point): ...
    def size(self) -> int: ...
    def dimensions(self) -> int: ...
    def is_in_search_region(self, point: list[float]) -> bool: ...
    def find_containing_bound(self, point: list[float]) -> LocalBound | None: ...
    bounds: list[LocalBound]  # read-only property
```

## Visualization & Interactive Dashboard

The library provides comprehensive 2D and 3D visualization tools, neighbor graph rendering, local bounds tables, and an interactive dashboard mirroring the figures from Klamroth et al. (2015) and Dächert et al. (2017).

All visualization code is organized in the [`visualization/`](visualization/) package with HTML templates decoupled in `visualization/templates/index.html`.

### Prerequisites & Installation

To run the visualization dashboard or generate reports, make sure your virtual environment has the required dependencies:

```bash
# 1. Activate your virtual environment
source venv/bin/activate

# 2. Install the current library and visualization dependencies
python -m pip install -e ".[vis]"
```

For development, the editable install loads the Python wrapper from this checkout. Rerun the install command after changing C++ headers or bindings to rebuild the native extension. Running the dashboard from the checkout with an older installed library can otherwise produce missing attributes or incompatible graph-export arguments.

`visualization/requirements.txt` installs only visualization dependencies (Flask, Plotly, NetworkX, and NumPy); it does not install or rebuild the library.

### Interactive Web Dashboard

Launch the interactive dashboard:

```bash
python -m visualization.app --port 8050
```

Open `http://127.0.0.1:8050` in your browser. Features include:
- **In-Place Problem Configuration**: Directly switch between 2D and 3D, toggle Minimization ($U(N)$) vs Maximization ($L(N)$), and customize Lower Bound ($LB$) and Upper Bound ($UB$) search spaces without intrusive popups.
- **Custom Instances**: Click **New custom instance** beside the dimension selector to start empty with the displayed dimension, sense, and interval bounds, then enter points. Applying settings or loading a preset starts a new problem and resets the views.
- **Strict Point Validation**: Automatically checks for dimension matching, finite coordinates, search space containment $[LB, UB]$, duplicate points, and Pareto dominance violations (points dominated by $N$ or dominating existing points in $N$) with clean, non-intrusive inline error banners.
- **Generation Timeline Slider**: Step forward and backward through point insertions, tracking destroyed, created, and persistent bounds.
- **3D & 2D Spatial Views**:
  - 3D view: Translucent (30% opacity by default), uniquely color-coded Pareto dominance cones $D(z) = [z, M]$ (matching Fig. 2 in Klamroth et al. 2015), opposite-axis camera perspective (looking from Ideal $m$ towards Nadir $M$), real-time interactive opacity slider, and single unified legend entry for one-click toggling.
  - Centered Pairwise 2D Projections: Centered matrix of projections $(f_1, f_2)$, $(f_1, f_3)$, $(f_2, f_3)$ showing bounding boxes.
  - 2D view: 2D search rectangles and staircase Pareto front.
- **View Controls**: Camera, zoom, and surviving graph node positions persist while stepping or highlighting. Spatial and pairwise views retain their own zoom. Use **Wireframe** to remove filled 3D faces and **Labels** to show all labels by default or hide them. Labels use white boxes with dark 12 px text; the reference labels are **m** and **M**. Coincident pairwise markers are grouped; clicking a bound group offers each underlying bound for selection.
- **Neighbor Graph Visualization**: Choose contracted adjacency or the raw Algorithm 1 graph, retaining quasi-bounds as labeled ellipses. Both views show total, quasi, and nonredundant counts and support component arrows or undirected edges. `get_adjacency_graph(include_quasi=True)` exports the raw topology and per-node `quasi` flags.
- **Membership Probe**: Check a point at the selected step without inserting it. The dashboard identifies and highlights a nonredundant containing bound and displays exact search-zone inequalities, including the strict reference-side boundary.

The dashboard computes history from bound IDs without building unused graph snapshots. It retains up to eight step snapshots and eight serialized figures, keyed by problem settings, ordered point IDs/coordinates, and view options; every API mutation clears these caches. Public step data is copied to keep snapshots independent. The CLI constructs only the final state unless `--step-by-step` is requested.

CLI HTML reports embed Plotly and work offline. The interactive dashboard serves Plotly locally but currently requires internet access for Tailwind, Lucide, and Vis-Network CDN assets; bundling those dashboard dependencies remains a deployment enhancement.
- **Local Bounds Table**: Displays coordinates $u$, defining points $z^j(u)$ (e.g. $z^1, z^2, z^3$), and neighbor pointers $\nu_k(u)$, with cross-highlighting across scenes and graphs.
- **Preset Library**: One-click presets for Paper 2 Example 2.8, Paper 1 Example 2 (SA), Paper 1 Example 3 (NGP Ties), 2D, and 3D Maximization.

### Command-Line Visualization Script

Generate standalone HTML reports from terminal:

```bash
# 3D with Paper 2 Example 2.8 preset and step-by-step terminal output
python -m visualization.cli --dim 3 --preset paper2 --step-by-step --output report_3d.html

# 2D example
python -m visualization.cli --dim 2 --preset 2d --step-by-step --output report_2d.html

# Custom points
python -m visualization.cli --dim 3 --points "[[4,0,4],[3,3,1],[2,2,2]]" --output custom.html

# Or launch the dashboard directly via CLI
python -m visualization.cli --interactive --port 8050
```

### Python API Visualization

```python
import local_bounds as lb
import visualization.engine as vis

# Initialize bound set
nbs = lb.NeighborhoodBoundSet([10.0, 10.0, 10.0], [0.0, 0.0, 0.0])
nbs.update(lb.Point("z1", [4.0, 0.0, 4.0]))
nbs.update(lb.Point("z2", [3.0, 3.0, 1.0]))
nbs.update(lb.Point("z3", [2.0, 2.0, 2.0]))

# Extract bounds, defining points, and neighbor relationships
data = vis.extract_bounds_data(nbs, [10.0, 10.0, 10.0], [0.0, 0.0, 0.0])

# Generate bounds table
table = vis.create_bounds_table(data)

# Generate Plotly figures
fig_3d = vis.plot_3d_bounds(data, show_occupied_boxes=True)
fig_graph = vis.plot_neighbor_graph_plotly(data, mode="combined")

# Show or export
fig_3d.write_html("bounds_3d.html")
fig_graph.write_html("graph.html")
```


## Benchmarks

The library includes a benchmarking tool (`benchmark/benchmark.cpp`) that evaluates the performance of the various algorithms by replicating the experiments from the original papers. It generates random stable sets of nondominated points and measures:

- The final number of local bounds (`|U(N)|`)
- The computation time (in milliseconds) for each algorithm
- The average number of bounds updated per point insertion (`|A|`)

To build and run the benchmarks, make sure to compile in Release mode:

```bash
cmake -B build -DBUILD_BENCHMARK=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build
./build/benchmark
```

The tool will run iterations for both minimization and maximization, across both *General Position* and *General Case* instance types. Results are printed to standard output and automatically saved to `benchmark_results.txt`.

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
