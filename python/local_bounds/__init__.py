from ._core import (
    Objective,
    Point,
    LocalBound,
    BoundSetMinimize,
    BoundSetMaximize,
    NeighborhoodBoundSetMinimize,
    NeighborhoodBoundSetMaximize,
    BoundSetTreeMinimize,
    BoundSetTreeMaximize
)

class BoundSet:
    def __init__(self, reference_point, anti_reference=None, sense=Objective.MINIMIZE):
        self.sense = sense
        if sense == Objective.MINIMIZE:
            if anti_reference is not None:
                self._impl = BoundSetMinimize(reference_point, anti_reference)
            else:
                self._impl = BoundSetMinimize(reference_point)
        else:
            if anti_reference is not None:
                self._impl = BoundSetMaximize(reference_point, anti_reference)
            else:
                self._impl = BoundSetMaximize(reference_point)

    def update_auto(self, point): return self._impl.update_auto(point)
    def update_re(self, point): return self._impl.update_re(point)
    def update_re_enhanced(self, point): return self._impl.update_re_enhanced(point)
    def update_ra_sa(self, point): return self._impl.update_ra_sa(point)
    def update_ra(self, point): return self._impl.update_ra(point)
    def update_naive(self, point): return self._impl.update_naive(point)
    def size(self): return self._impl.size()
    def dimensions(self): return self._impl.dimensions()
    def is_in_search_region(self, point): return self._impl.is_in_search_region(point)
    def find_containing_bound(self, point): return self._impl.find_containing_bound(point)
    @property
    def bounds(self): return self._impl.bounds

class NeighborhoodBoundSet:
    def __init__(self, reference_point, anti_reference, sense=Objective.MINIMIZE):
        self.sense = sense
        if sense == Objective.MINIMIZE:
            self._impl = NeighborhoodBoundSetMinimize(reference_point, anti_reference)
        else:
            self._impl = NeighborhoodBoundSetMaximize(reference_point, anti_reference)

    def update(self, point): return self._impl.update(point)
    def size(self): return self._impl.size()
    def nonredundant_size(self): return self._impl.nonredundant_size()
    def dimensions(self): return self._impl.dimensions()
    def is_in_search_region(self, point): return self._impl.is_in_search_region(point)
    def find_containing_bound(self, point): return self._impl.find_containing_bound(point)
    def bounds(self): return self._impl.bounds()
    def nonredundant_bounds(self): return self._impl.nonredundant_bounds()
    def get_adjacency_graph(self): return self._impl.get_adjacency_graph()

class BoundSetTree:
    def __init__(self, reference_point, anti_reference=None, max_leaf_size=32, num_children=8, sense=Objective.MINIMIZE):
        self.sense = sense
        if sense == Objective.MINIMIZE:
            if anti_reference is not None:
                self._impl = BoundSetTreeMinimize(reference_point, anti_reference, max_leaf_size, num_children)
            else:
                self._impl = BoundSetTreeMinimize(reference_point, max_leaf_size, num_children)
        else:
            if anti_reference is not None:
                self._impl = BoundSetTreeMaximize(reference_point, anti_reference, max_leaf_size, num_children)
            else:
                self._impl = BoundSetTreeMaximize(reference_point, max_leaf_size, num_children)

    def update_auto(self, point): return self._impl.update_auto(point)
    def update_re(self, point): return self._impl.update_re(point)
    def update_re_enhanced(self, point): return self._impl.update_re_enhanced(point)
    def update_naive(self, point): return self._impl.update_naive(point)
    def size(self): return self._impl.size()
    def dimensions(self): return self._impl.dimensions()
    def is_in_search_region(self, point): return self._impl.is_in_search_region(point)
    def find_containing_bound(self, point): return self._impl.find_containing_bound(point)
    @property
    def bounds(self): return self._impl.bounds

__all__ = [
    "Objective",
    "Point",
    "LocalBound",
    "BoundSet",
    "NeighborhoodBoundSet",
    "BoundSetTree"
]
