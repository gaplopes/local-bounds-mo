from local_bounds import Point, BoundSet, BoundSetTree, NeighborhoodBoundSet, Objective

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
                    def export():
                        values = bounds.bounds
                        return values() if callable(values) else values
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


if __name__ == "__main__":
    test_basic()
    test_neighborhood()
    test_update_results()
