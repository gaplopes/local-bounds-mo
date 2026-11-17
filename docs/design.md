# Design notes

Details for people who want to know what each update costs or who work on the
library itself. The [README](../README.md) has everything needed to use it.

Notation, as in the papers:

- N is the set of points inserted so far and p the number of objectives.
- U(N) is the set of local bounds.
- A is the set of bounds whose search zones contain the new point z̄. Paper 2
  calls the same set U_z̄.
- B is the set of bounds that have z̄ on the boundary of their search zone.
- *General position* means no two points share a value in any objective
  (Paper 1 calls it SA, Paper 2 calls it GP). Otherwise there are *ties*.

## Cost of one update

**Algorithms 2–5** (Paper 1) start by finding A. For fixed p this costs
O(|U(N)|) with a linear scan, or O(log^p |U(N)| + |A|) with a range tree
(Paper 1, Section 5.2). `BoundSet` scans linearly. `BoundSetTree` uses its
`LBTree` index, whose pruning depends on the workload.

The table gives the cost of what follows, candidate generation and
filtering or avoidance, for fixed p:

| Method | Cost | Notes |
|---|---|---|
| `update_naive()` | O(\|A\|² + \|A\|·\|U\|) | Generates p\|A\| candidates and filters them pairwise and against all bounds |
| `update_re()` | O(\|A\|·(\|A\| + \|B\|)) | See below |
| `update_re_enhanced()` | O(\|A\|·(\|A\| + \|B\|)) | Algorithm 3. General position: O(\|A\|²) (Prop. 5.1) |
| `update_ra_sa()` | O(\|A\|) | Algorithm 4 (Prop. 5.2); no filtering |
| `update_ra()` | O(\|N\|·\|A\|) worst case | Algorithm 5; the sets Z^k(u) are usually small |

`update_re()` and `update_re_enhanced()` return the same bounds. Both filter a
candidate only against candidates of the same component and against the bounds
tied in that component, which is the refinement of Paper 1's Proposition 3.2
and Algorithm 3. They differ in how B is collected: `update_re()` takes every
weakly dominated bound, `update_re_enhanced()` only those with strict
inequality in the other components. The paper's own Algorithm 2 filters each
candidate against all of P ∪ B; this library does not implement that variant.

The paper notes that the filtering of Algorithm 3 can be reduced to
O(|A| log |A|) for p ∈ {2, 3} and O(|A| log^(p-3) |A| log log |A|) for p ≥ 4
with specialized algorithms. Those are not implemented here.

RE also scans the current bounds to build the per-component sets B_j, and
`update_ra()` scans them to update the defining sets. These passes come on top
of the costs above.

|A| grows quickly with p. Paper 1 reports averages of about 4 for p = 3, 22 for
p = 4, 142 for p = 5 and 736 for p = 6 on its random instances.

**Algorithm 1** (Paper 2) walks the affected region through the neighborhood
graph. The paper's O(|A|) bound assumes the caller supplies one bound that
contains the new point (Section 5.1). This implementation finds that bound by
scanning the node storage, and sets up a visited array over the same storage.
For fixed p an update therefore costs O(C + |A|), where C is the number of
allocated node slots, including inactive ones.

## Quasi-nonredundant bounds

With ties, Algorithm 1 keeps *quasi-nonredundant* bounds: bounds whose search
zone is contained in another one. They stay in the graph because removing them
could disconnect it (Paper 2, Section 4.3). `bounds()` returns them;
`nonredundant_bounds()`, `nonredundant_size()` and `get_adjacency_graph()`
filter them out.

A bound is quasi-nonredundant when another active bound contains its search
zone. Among bounds with equal coordinates, the one in the lowest slot is kept.

Only *tie-born* bounds are tested: bounds created where condition (C2) held
with equality, and their descendants. A bound with no tie in its ancestry keeps
strict defining points, so it is a local upper bound by Corollary 2.6 of
Paper 2 and is accepted without a test.
`tests/test_neighborhood_bound_set.cpp` checks this against a full pairwise
scan on random tied inputs.

The test costs O(p·T·C) time and O(C) temporary memory, where T is the number
of tie-born bounds. T is zero in general position. When every bound is
tie-born the cost is O(p·C²). Graph export also walks and copies the adjacency
lists. Each call computes its own mask, so nothing has to be invalidated after
an update. `get_clean_neighbors()` runs the same test for the nodes it visits.

## IDs and order in `BoundSetTree`

`BoundSetTree` does not store IDs. `bounds()` numbers the bounds `u0`, `u1`,
... in the tree's depth-first order, and `find_containing_bound()` returns the
first match in that order with the same ID. The bounds returned by an update
carry the IDs they had in the export before that update. These IDs change
whenever the tree changes; they are positions, not identifiers.

To return an ID without visiting every bound, each tree node stores how many
bounds lie below it (`subtree_size`). When a query skips a subtree because its
bounding box cannot match, it adds that count to the running position.

## Ownership

`LBTree` and `BoundSetTree` cannot be copied or moved: their nodes point into
an arena owned by the tree.
