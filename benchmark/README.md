# Benchmarks

Build in Release mode:

```bash
cmake -B build -DBUILD_BENCHMARK=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

## Algorithm comparison

```bash
./build/benchmark
```

`benchmark.cpp` follows the experiments of both papers. It generates random
stable sets of points and runs every update method on them, for minimization
and maximization, in general position and with ties. For each run it reports:

- the final number of local bounds, |U(N)|
- the update time in milliseconds
- the export time, measured separately and without the sorting used for checks
- the average number of bounds replaced per insertion, |A|

Every method must produce the same bounds as `update_naive()`, or the benchmark
stops with an error. Results go to `benchmark_results.txt`, or to standard
output if that file cannot be opened.

## Tree and classification diagnostics

`benchmark_lbtree` has two modes. Both cover both senses, both instance types
and 3 to 6 objectives, and report the median of three measurements.

**Insertion order** (default). Measures how the order of insertion changes the
tree's depth, build time and lookup time.

```bash
./build/benchmark_lbtree > build/lbtree.csv
./build/benchmark_lbtree 1 > build/lbtree_smoke.csv    # one shuffle instead of three
```

The same bound snapshots (instance seed 42; near 500, 1,000, 3,000 and 5,000
bounds) are inserted in native, sorted and shuffled order with the default leaf
size 32 and 8 children. The shuffles use seeds 42 to 44, shown in the `order`
column. Columns: `dimensions,sense,instance,bounds,order,depth,build_ms,query_ms`.
`depth` counts nodes from the root to the deepest leaf. `build_ms` covers
construction, insertion and destruction. Sorting, shuffling and the correctness
checks are not timed, so this shows sensitivity to order, not total cost. Each
order must keep exactly the same bounds and answer a fixed batch of queries
like a linear scan does.

**Classification** (`--classification`). Compares ways of finding the
quasi-nonredundant bounds of a `NeighborhoodBoundSet`.

```bash
./build/benchmark_lbtree --classification > build/classification.csv
./build/benchmark_lbtree --classification 1 > build/classification_smoke.csv    # one seed instead of three
```

It samples snapshots between about 500 and 5,000 active bounds. `scan_ms` is a
pairwise scan of the exported snapshot and `index_ms` the same job with an
`LBTree`, including building it. `auto_count_ms` times the library's
`nonredundant_size()` and `graph_ms` a full `get_adjacency_graph()`. All three
results are checked against the pairwise scan. The generator may run out of
attempts before reaching the largest sizes.
