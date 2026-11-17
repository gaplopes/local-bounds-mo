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

def _normalize_sense(sense):
    if isinstance(sense, str) and sense in ("MINIMIZE", "MAXIMIZE"):
        sense = getattr(Objective, sense)
    if sense not in (Objective.MINIMIZE, Objective.MAXIMIZE):
        raise ValueError("Sense must be MINIMIZE or MAXIMIZE")
    return sense


def _pick(sense, minimize, maximize):
    return minimize if _normalize_sense(sense) == Objective.MINIMIZE else maximize


# The C++ classes take the sense as a template parameter, so _core has one
# class per sense. These factories return an instance of the matching one.

def BoundSet(reference_point, anti_reference=None, sense=Objective.MINIMIZE):
    cls = _pick(sense, BoundSetMinimize, BoundSetMaximize)
    if anti_reference is None:
        return cls(reference_point)
    return cls(reference_point, anti_reference)


def NeighborhoodBoundSet(reference_point, anti_reference, sense=Objective.MINIMIZE):
    cls = _pick(sense, NeighborhoodBoundSetMinimize, NeighborhoodBoundSetMaximize)
    return cls(reference_point, anti_reference)


def BoundSetTree(reference_point, anti_reference=None, max_leaf_size=32, num_children=8,
                 sense=Objective.MINIMIZE):
    cls = _pick(sense, BoundSetTreeMinimize, BoundSetTreeMaximize)
    if anti_reference is None:
        return cls(reference_point, max_leaf_size, num_children)
    return cls(reference_point, anti_reference, max_leaf_size, num_children)


__all__ = [
    "Objective",
    "Point",
    "LocalBound",
    "BoundSet",
    "NeighborhoodBoundSet",
    "BoundSetTree"
]
