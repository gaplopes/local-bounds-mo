from local_bounds import Point, BoundSet, BoundSetTree, NeighborhoodBoundSet, Objective


def _raises(error, call, *args, **kwargs):
    try:
        call(*args, **kwargs)
    except error:
        return
    raise AssertionError(f"expected {error.__name__}")


def test_install_adds_nothing_outside_the_package():
    # The C++ headers and CMake files belong to the CMake install, not to site-packages.
    from importlib.metadata import files
    top_level = {path.parts[0] for path in files("local_bounds")}
    assert not top_level & {"include", "lib", "lib64", "share"}, sorted(top_level)


def test_basic():
    print("Testing local_bounds python bindings...")
    ref = [10.0, 10.0, 10.0]
    bs = BoundSet(ref, sense=Objective.MINIMIZE)
    
    print(f"Initial size: {bs.size()}")
    assert bs.size() == 1
    
    p = Point("p1", [1.0, 2.0, 3.0])
    updated, dominated = bs.update_auto(p)
    assert updated is True
    assert [b.coordinates for b in dominated] == [ref]
    
    print(f"Size after 1 point: {bs.size()}")
    assert bs.size() == 3

    # Dominated point should return (False, [])
    p_dom = Point("p_dom", [4.0, 5.0, 6.0])
    assert bs.update_auto(p_dom) == (False, [])
    assert bs.size() == 3

def test_bounds_is_a_method_on_every_class():
    for sense, reference, anti, coordinates in (
        (Objective.MINIMIZE, [10.0]*3, [0.0]*3, [1.0, 2.0, 3.0]),
        (Objective.MAXIMIZE, [0.0]*3, [10.0]*3, [9.0, 8.0, 7.0]),
    ):
        for cls, update in ((BoundSet, "update_re"), (BoundSetTree, "update_re"),
                            (NeighborhoodBoundSet, "update")):
            bounds = cls(reference, anti, sense=sense)
            initial = bounds.bounds()
            assert isinstance(initial, list)
            assert [b.coordinates for b in initial] == [reference]
            updated, _ = getattr(bounds, update)(Point("p1", coordinates))
            assert updated
            assert len(bounds.bounds()) == bounds.size() == 3
            assert [b.coordinates for b in initial] == [reference]  # earlier export is a snapshot


def test_sense_is_reported_and_validated():
    assert BoundSet([1.0], sense="MAXIMIZE").sense == Objective.MAXIMIZE
    assert BoundSetTree([1.0]).sense == Objective.MINIMIZE
    assert NeighborhoodBoundSet([1.0], [0.0]).sense == Objective.MINIMIZE
    for cls in (BoundSet, BoundSetTree):
        _raises(ValueError, cls, [1.0], sense="sideways")


def test_interval_boundaries_and_argument_errors():
    # The anti-reference boundary is inside the interval; the reference boundary is not.
    for sense, valid, invalid in ((Objective.MINIMIZE, [0, 5], ([-5e-8, 5], [10, 5])),
                                  (Objective.MAXIMIZE, [10, 5], ([10 + 5e-8, 5], [0, 5]))):
        ref, anti = ([10, 10], [0, 0]) if sense == Objective.MINIMIZE else ([0, 0], [10, 10])
        for cls, name in ((BoundSet, "update_re"), (BoundSetTree, "update_re"),
                          (NeighborhoodBoundSet, "update")):
            bounds = cls(ref, anti, sense=sense)
            for coords in invalid:
                _raises(ValueError, getattr(bounds, name), Point("z", coords))
            assert getattr(bounds, name)(Point("z", valid))[0]
            _raises(ValueError, bounds.find_containing_bound, [1])
    # RA needs the anti-reference point for its dummy defining points.
    _raises(RuntimeError, BoundSet([10, 10]).update_ra, Point("z", [3, 7]))


def test_raw_graph_keeps_quasi_bounds_and_is_symmetric():
    tied = [[5, 2, 5], [2, 4, 3], [4, 4, 1], [4, 3, 2], [3, 3, 3]]
    for sense in (Objective.MINIMIZE, Objective.MAXIMIZE):
        is_min = sense == Objective.MINIMIZE
        bounds = NeighborhoodBoundSet([6]*3 if is_min else [0]*3, [0]*3 if is_min else [6]*3, sense=sense)
        for i, z in enumerate(tied):
            bounds.update(Point(f"z{i}", z if is_min else [6 - c for c in z]))
        raw = bounds.get_adjacency_graph(include_quasi=True)
        assert len(raw.nodes) == bounds.size() > 8
        assert sum(not q for q in raw.quasi) == bounds.nonredundant_size() == 8
        for u, neighbors in enumerate(raw.adjacency_list):
            for v in neighbors:
                assert u in raw.adjacency_list[v]


def test_neighborhood():
    print("Testing NeighborhoodBoundSet...")
    ref = [10.0, 10.0, 10.0]
    anti = [0.0, 0.0, 0.0]
    nbs = NeighborhoodBoundSet(ref, anti, sense=Objective.MINIMIZE)
    
    print(f"Initial size: {nbs.size()}")
    assert nbs.size() == 1
    
    p = Point("p1", [1.0, 2.0, 3.0])
    updated, dominated = nbs.update(p)
    assert updated is True
    assert [b.coordinates for b in dominated] == [ref]
    
    print(f"Size after 1 point: {nbs.nonredundant_size()}")
    assert nbs.nonredundant_size() == 3

    # Dominated point should return (False, [])
    p_dom = Point("p_dom", [4.0, 5.0, 6.0])
    assert nbs.update(p_dom) == (False, [])
    assert nbs.nonredundant_size() == 3
    print("Success!")


def test_update_results():
    def snapshot(bound):
        def points(values):
            return tuple((p.id, tuple(p.coordinates)) for p in values)
        return (bound.id, tuple(bound.coordinates), points(bound.defining_points),
                tuple(points(values) for values in bound.defining_point_sets))

    methods = ["update_re", "update_re_enhanced", "update_ra_sa", "update_ra",
               "update_naive", "update_auto"]
    cases = [(BoundSet, methods), (BoundSetTree, [m for m in methods if "ra" not in m]),
             (NeighborhoodBoundSet, ["update"])]
    for sense in (Objective.MINIMIZE, Objective.MAXIMIZE):
        def orient(coords):
            return coords if sense == Objective.MINIMIZE else [10 - x for x in coords]

        for cls, names in cases:
            for name in names:
                for tied in (False, True):
                    if tied and name == "update_ra_sa":
                        continue
                    options = {"max_leaf_size": 2, "num_children": 2} if cls is BoundSetTree else {}
                    bounds = cls(orient([10]*3), orient([0]*3), sense=sense, **options)
                    update = getattr(bounds, name)
                    export = bounds.bounds
                    data = ([[5,2,5], [2,4,3], [4,4,1], [4,3,2], [3,3,3]] if tied else
                            [[3,5,7], [6,2,4], [2,6,5], [5,4,3]])
                    saved = []
                    saw_multiple = False
                    for i, coords in enumerate(data):
                        before = export()
                        z = orient(coords)
                        expected = sorted(snapshot(b) for b in before if all(
                            x < y if sense == Objective.MINIMIZE else x > y
                            for x, y in zip(z, b.coordinates)))
                        result = update(Point(f"z{i}", z))
                        assert isinstance(result, tuple) and len(result) == 2
                        updated, dominated = result
                        assert updated is bool(expected)
                        assert isinstance(dominated, list)
                        assert sorted(snapshot(b) for b in dominated) == expected
                        saw_multiple |= len(dominated) > 1
                        saved.append((dominated, sorted(snapshot(b) for b in dominated)))
                        assert update(Point(f"duplicate{i}", z)) == (False, [])
                    assert saw_multiple
                    for dominated, expected in saved:
                        assert sorted(snapshot(b) for b in dominated) == expected
                    before = export()
                    expected = sorted(snapshot(b) for b in before)
                    assert update(Point("dominated", orient([9]*3))) == (False, [])
                    assert sorted(snapshot(b) for b in export()) == expected
                    if cls is BoundSetTree:
                        for coords in ([1,1,1], [0,9,9]):
                            z = orient(coords)
                            expected = sorted(snapshot(b) for b in export() if all(
                                x < y if sense == Objective.MINIMIZE else x > y
                                for x, y in zip(z, b.coordinates)))
                            updated, dominated = update(Point("replacement", z))
                            assert updated is bool(expected)
                            assert sorted(snapshot(b) for b in dominated) == expected


def test_high_dimensional_bound_ids():
    # Objective paths "11" and "1,1" collide without separators.
    coordinates = [[10]*11, [5] + [20]*10]
    cases = [(BoundSet, ["update_re", "update_re_enhanced", "update_ra_sa",
                         "update_ra", "update_naive", "update_auto"]),
             (NeighborhoodBoundSet, ["update"])]
    for sense in (Objective.MINIMIZE, Objective.MAXIMIZE):
        def orient(values):
            return values if sense == Objective.MINIMIZE else [100 - x for x in values]

        for cls, methods in cases:
            for method in methods:
                bounds = cls(orient([100]*11), orient([0]*11), sense=sense)
                for i, values in enumerate(coordinates):
                    getattr(bounds, method)(Point(f"p{i}", orient(values)))
                exported = bounds.bounds()
                assert len(exported) == 21
                assert len({b.id for b in exported}) == len(exported), (cls.__name__, method, sense)


def test_indexed_neighborhood_exports():
    for sense in (Objective.MINIMIZE, Objective.MAXIMIZE):
        def orient(coordinates):
            return [(x if sense == Objective.MINIMIZE else 100000 - x) / 4
                    for x in coordinates]

        bounds = NeighborhoodBoundSet(orient([100000]*3), orient([0]*3), sense=sense)
        for i in range(1, 1601):
            bounds.update(Point(f"p{i}", orient([10*i, 100000 - 10*i,
                                                 10*(1 + (37*i) % 9999)])))
        exported = bounds.nonredundant_bounds()
        raw = bounds.get_adjacency_graph(True)
        clean = bounds.get_adjacency_graph()
        assert bounds.size() == bounds.nonredundant_size() == len(exported) == 3201
        assert not any(raw.quasi)
        expected = [(b.id, b.coordinates) for b in exported]
        assert [(b.id, b.coordinates) for b in raw.nodes] == expected
        assert [(b.id, b.coordinates) for b in clean.nodes] == expected
        assert raw.adjacency_list == clean.adjacency_list
        assert raw.k_neighbors == clean.k_neighbors

        bounds.update(Point("collapse", orient([1, 1, 1])))
        assert bounds.size() == bounds.nonredundant_size() == 3
        bounds.update(Point("reuse", orient([0, 2, 2])))
        assert len(bounds.nonredundant_bounds()) == bounds.nonredundant_size()
        assert len(bounds.get_adjacency_graph().nodes) == bounds.nonredundant_size()


if __name__ == "__main__":
    test_install_adds_nothing_outside_the_package()
    test_basic()
    test_bounds_is_a_method_on_every_class()
    test_sense_is_reported_and_validated()
    test_interval_boundaries_and_argument_errors()
    test_raw_graph_keeps_quasi_bounds_and_is_symmetric()
    test_neighborhood()
    test_update_results()
    test_high_dimensional_bound_ids()
    test_indexed_neighborhood_exports()
