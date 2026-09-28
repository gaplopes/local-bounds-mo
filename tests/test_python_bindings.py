from local_bounds import Point, BoundSet, NeighborhoodBoundSet, Objective

def test_basic():
    print("Testing local_bounds python bindings...")
    ref = [10.0, 10.0, 10.0]
    bs = BoundSet(ref, sense=Objective.MINIMIZE)
    
    print(f"Initial size: {bs.size()}")
    assert bs.size() == 1
    
    p = Point("p1", [1.0, 2.0, 3.0])
    bs.update_auto(p)
    
    print(f"Size after 1 point: {bs.size()}")
    assert bs.size() == 3

def test_neighborhood():
    print("Testing NeighborhoodBoundSet...")
    ref = [10.0, 10.0, 10.0]
    anti = [0.0, 0.0, 0.0]
    nbs = NeighborhoodBoundSet(ref, anti, sense=Objective.MINIMIZE)
    
    print(f"Initial size: {nbs.size()}")
    assert nbs.size() == 1
    
    p = Point("p1", [1.0, 2.0, 3.0])
    nbs.update(p)
    
    print(f"Size after 1 point: {nbs.nonredundant_size()}")
    assert nbs.nonredundant_size() == 3
    print("Success!")

if __name__ == "__main__":
    test_basic()
    test_neighborhood()
