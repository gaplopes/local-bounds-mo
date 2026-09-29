from local_bounds import Point, BoundSet, NeighborhoodBoundSet, Objective

def test_basic():
    print("Testing local_bounds python bindings...")
    ref = [10.0, 10.0, 10.0]
    bs = BoundSet(ref, sense=Objective.MINIMIZE)
    
    print(f"Initial size: {bs.size()}")
    assert bs.size() == 1
    
    p = Point("p1", [1.0, 2.0, 3.0])
    updated = bs.update_auto(p)
    assert updated is True
    
    print(f"Size after 1 point: {bs.size()}")
    assert bs.size() == 3

    # Dominated point should return False
    p_dom = Point("p_dom", [4.0, 5.0, 6.0])
    assert bs.update_auto(p_dom) is False
    assert bs.size() == 3

def test_neighborhood():
    print("Testing NeighborhoodBoundSet...")
    ref = [10.0, 10.0, 10.0]
    anti = [0.0, 0.0, 0.0]
    nbs = NeighborhoodBoundSet(ref, anti, sense=Objective.MINIMIZE)
    
    print(f"Initial size: {nbs.size()}")
    assert nbs.size() == 1
    
    p = Point("p1", [1.0, 2.0, 3.0])
    updated = nbs.update(p)
    assert updated is True
    
    print(f"Size after 1 point: {nbs.nonredundant_size()}")
    assert nbs.nonredundant_size() == 3

    # Dominated point should return False
    p_dom = Point("p_dom", [4.0, 5.0, 6.0])
    assert nbs.update(p_dom) is False
    assert nbs.nonredundant_size() == 3
    print("Success!")

if __name__ == "__main__":
    test_basic()
    test_neighborhood()
