"""
Unit and Integration Tests for Local Bounds Visualization, Table, Neighbor Graphs, and Interactive App.
"""

import unittest
import json
import contextlib
import io
import tempfile
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from unittest.mock import patch
import local_bounds as lb
import visualization.engine as vis
import visualization.app as web_app


class TestVisualization(unittest.TestCase):

    def test_example_2_8_paper2(self):
        """Validates bounds, coordinates, defining points (z^j), and neighbor graphs against Paper 2 Example 2.8."""
        ref = [10.0, 10.0, 10.0]
        anti = [0.0, 0.0, 0.0]
        nbs = lb.NeighborhoodBoundSet(ref, anti)

        # Step 1: z1 = (4, 0, 4)
        z1 = lb.Point("z1", [4.0, 0.0, 4.0])
        nbs.update(z1)
        self.assertEqual(nbs.size(), 3)

        # Step 2: z2 = (3, 3, 1)
        z2 = lb.Point("z2", [3.0, 3.0, 1.0])
        nbs.update(z2)
        self.assertEqual(nbs.size(), 5)

        # Step 3: z3 = (2, 2, 2)
        z3 = lb.Point("z3", [2.0, 2.0, 2.0])
        nbs.update(z3)
        self.assertEqual(nbs.size(), 7)

        # Extract data
        data = vis.extract_bounds_data(nbs, ref, anti, [z1, z2, z3])
        self.assertEqual(len(data["bounds"]), 7)

        # Check extreme bounds: exactly 3 extreme bounds in 3D
        extreme_bounds = [b for b in data["bounds"] if b.is_extreme]
        self.assertEqual(len(extreme_bounds), 3)

        extreme_dims = set(b.extreme_dim for b in extreme_bounds)
        self.assertEqual(extreme_dims, {0, 1, 2})

        # Expected coordinates from Example 2.8 (3)
        expected_coords = {
            (2.0, 10.0, 10.0),
            (3.0, 10.0, 2.0),
            (4.0, 2.0, 10.0),
            (10.0, 0.0, 10.0),
            (10.0, 2.0, 4.0),
            (10.0, 3.0, 2.0),
            (10.0, 10.0, 1.0)
        }
        got_coords = set(tuple(b.coordinates) for b in data["bounds"])
        self.assertEqual(got_coords, expected_coords)

        # Check table formatting and z^3 notation
        table = vis.create_bounds_table(data)
        self.assertEqual(len(table), 7)
        found_z3_dummy = False
        for row in table:
            self.assertIn("ID", row)
            self.assertIn("Coordinates", row)
            self.assertIn("z^1(u)", row)
            self.assertIn("z^2(u)", row)
            self.assertIn("z^3(u)", row)
            self.assertIn("ν_1", row)
            self.assertIn("ν_2", row)
            self.assertIn("ν_3", row)
            # Ensure z_hat3 is NEVER present, replaced with z^3
            self.assertNotIn("z_hat3", row["z^3(u)"])
            if row["z^3(u)"] == "z^3":
                found_z3_dummy = True

        self.assertTrue(found_z3_dummy, "Dummy point for component 3 should be displayed as z^3")

    def test_2d_bounds_and_graph(self):
        """Validates 2D local bounds generation and neighbor relationships."""
        ref = [10.0, 10.0]
        anti = [0.0, 0.0]
        nbs = lb.NeighborhoodBoundSet(ref, anti)

        pts = [
            lb.Point("z1", [3.0, 7.0]),
            lb.Point("z2", [5.0, 4.0]),
            lb.Point("z3", [7.0, 2.0])
        ]
        for p in pts:
            nbs.update(p)

        self.assertEqual(nbs.dimensions(), 2)
        self.assertEqual(nbs.size(), 4)

        data = vis.extract_bounds_data(nbs, ref, anti, pts)
        self.assertEqual(len(data["bounds"]), 4)

        # Exactly 2 extreme bounds in 2D
        extreme_bounds = [b for b in data["bounds"] if b.is_extreme]
        self.assertEqual(len(extreme_bounds), 2)

        # Check graph extraction
        graph_elems = vis.extract_network_graph_elements(data, mode="combined")
        self.assertEqual(len(graph_elems["nodes"]), 4)
        self.assertGreater(len(graph_elems["edges"]), 0)

    def test_figures_generation(self):
        """Tests that 2D, 3D, and pairwise projection figure generators produce valid Plotly objects."""
        ref3 = [10.0, 10.0, 10.0]
        anti3 = [0.0, 0.0, 0.0]
        nbs3 = lb.NeighborhoodBoundSet(ref3, anti3)
        pt3 = lb.Point("z1", [4.0, 2.0, 5.0])
        nbs3.update(pt3)
        data3 = vis.extract_bounds_data(nbs3, ref3, anti3, [pt3])

        # 3D Plot (Figure 2 alignment)
        fig3d = vis.plot_3d_bounds(data3)
        self.assertIsNotNone(fig3d)
        self.assertGreater(len(fig3d.data), 0)
        self.assertEqual(fig3d.layout.scene.camera.eye.x, -1.65)
        self.assertEqual(fig3d.layout.scene.camera.eye.y, -1.65)
        self.assertTrue(any(getattr(t, "legendgroup", None) == "dominated_zones" for t in fig3d.data))
        self.assertFalse(any(getattr(t, "legendgroup", None) == "defining_rays" for t in fig3d.data))

        # Highlighted bound S(u)
        fig3d_hl = vis.plot_3d_bounds(data3, highlight_bound_id=data3["bounds"][0].id)
        self.assertTrue(any(getattr(t, "legendgroup", None) == "highlight_zone" for t in fig3d_hl.data))

        # Pairwise 2D Projections (centered)
        fig_pair = vis.plot_pairwise_projections_2d(data3)
        self.assertIsNotNone(fig_pair)
        self.assertEqual(fig_pair.layout.legend.itemclick, "toggle")
        self.assertTrue(any(t.legendgroup == "local_bounds" and t.showlegend for t in fig_pair.data))
        self.assertTrue(any(t.legendgroup == "points" and t.showlegend for t in fig_pair.data))

        # Graph Plotly
        fig_graph = vis.plot_neighbor_graph_plotly(data3, mode="combined")
        self.assertIsNotNone(fig_graph)

        # 2D Plot
        ref2 = [10.0, 10.0]
        anti2 = [0.0, 0.0]
        nbs2 = lb.NeighborhoodBoundSet(ref2, anti2)
        nbs2.update(lb.Point("z1", [3.0, 5.0]))
        data2 = vis.extract_bounds_data(nbs2, ref2, anti2)
        fig2d = vis.plot_2d_bounds(data2)
        self.assertIsNotNone(fig2d)

    def test_tracker_history(self):
        """Tests the step-by-step history tracking."""
        tracker = vis.LocalBoundsTracker([10.0, 10.0, 10.0], [0.0, 0.0, 0.0])
        self.assertEqual(len(tracker.steps), 1)
        self.assertEqual(tracker.steps[0].step, 0)
        self.assertEqual(len(tracker.steps[0].bounds_data["bounds"]), 1)

        step1 = tracker.add_point(lb.Point("z1", [4.0, 0.0, 4.0]))
        self.assertEqual(step1.step, 1)
        self.assertEqual(step1.destroyed_bounds, ["u0"])
        self.assertEqual(len(step1.bounds_data["bounds"]), 3)

        step2 = tracker.add_point(lb.Point("z2", [3.0, 3.0, 1.0]))
        self.assertEqual(step2.step, 2)
        self.assertEqual(len(step2.bounds_data["bounds"]), 5)

        step3 = tracker.add_point(lb.Point("z3", [2.0, 2.0, 2.0]))
        self.assertEqual(step3.step, 3)
        self.assertEqual(len(step3.bounds_data["bounds"]), 7)

    def test_flask_interactive_app(self):
        """Tests interactive application API endpoints and external template."""
        client = web_app.app.test_client()

        # 1. Main index (rendered via templates/index.html)
        res = client.get("/")
        self.assertEqual(res.status_code, 200)
        self.assertIn(b"Local Bounds Visualization", res.data)

        # 2. Plotly JS asset
        res_js = client.get("/static/plotly.js")
        self.assertEqual(res_js.status_code, 200)

        # 3. Get state
        res_state = client.get("/api/state")
        self.assertEqual(res_state.status_code, 200)
        state_data = json.loads(res_state.data)
        self.assertIn("dimensions", state_data)
        self.assertIn("table", state_data)
        self.assertIn("history", state_data)

        # 4. Load preset
        res_preset = client.post("/api/load_preset", json={"preset": "paper2"})
        self.assertEqual(res_preset.status_code, 200)
        p_data = json.loads(res_preset.data)
        self.assertEqual(p_data["bounds_count"], 7)

        # 5. Add point (mutually nondominated)
        res_add = client.post("/api/add_point", json={"id": "z4", "coordinates": [1.5, 5.0, 1.5]})
        self.assertEqual(res_add.status_code, 200)
        add_data = json.loads(res_add.data)
        self.assertEqual(len(add_data["points"]), 4)

        # 6. Delete point
        res_del = client.post("/api/delete_point", json={"index": 3})
        self.assertEqual(res_del.status_code, 200)
        del_data = json.loads(res_del.data)
        self.assertEqual(len(del_data["points"]), 3)

        # 7. Graph data
        res_graph = client.get("/api/graph?mode=combined")
        self.assertEqual(res_graph.status_code, 200)
        g_data = json.loads(res_graph.data)
        self.assertIn("nodes", g_data)
        self.assertIn("edges", g_data)

        # 8. Figure 3D
        res_fig3d = client.get("/api/figure_3d")
        self.assertEqual(res_fig3d.status_code, 200)

        # 9. Figure Pairwise Projections
        res_pair = client.get("/api/figure_pairwise")
        self.assertEqual(res_pair.status_code, 200)

        # 10. Configure 2D
        res_cfg = client.post("/api/configure", json={
            "dimensions": 2,
            "sense": "MINIMIZE",
            "lower_bound": [0.0, 0.0],
            "upper_bound": [10.0, 10.0]
        })
        self.assertEqual(res_cfg.status_code, 200)
        cfg_data = json.loads(res_cfg.data)
        self.assertEqual(cfg_data["dimensions"], 2)

    def test_point_validation_rules(self):
        """Verifies thorough validation checks for invalid, dominated, and duplicate points."""
        existing_pts = [
            lb.Point("z1", [4.0, 0.0, 4.0]),
            lb.Point("z2", [3.0, 3.0, 1.0]),
            lb.Point("z3", [2.0, 2.0, 2.0])
        ]
        lb_box = [0.0, 0.0, 0.0]
        ub_box = [10.0, 10.0, 10.0]

        # 1. Dimension Mismatch
        with self.assertRaises(vis.PointValidationError) as cm:
            vis.validate_point([1.0, 2.0], "p", existing_pts, lb_box, ub_box)
        self.assertEqual(cm.exception.code, "DIMENSION_MISMATCH")

        with self.assertRaises(vis.PointValidationError) as cm:
            vis.validate_point([1.0, 2.0, 3.0, 4.0], "p", existing_pts, lb_box, ub_box)
        self.assertEqual(cm.exception.code, "DIMENSION_MISMATCH")

        # 2. Non-numeric / NaN / Inf
        with self.assertRaises(vis.PointValidationError) as cm:
            vis.validate_point([float("nan"), 2.0, 3.0], "p", existing_pts, lb_box, ub_box)
        self.assertEqual(cm.exception.code, "NON_FINITE")

        with self.assertRaises(vis.PointValidationError) as cm:
            vis.validate_point([float("inf"), 2.0, 3.0], "p", existing_pts, lb_box, ub_box)
        self.assertEqual(cm.exception.code, "NON_FINITE")

        # 3. Out of bounds (below LB or above UB)
        with self.assertRaises(vis.PointValidationError) as cm:
            vis.validate_point([-0.5, 2.0, 3.0], "p", existing_pts, lb_box, ub_box)
        self.assertEqual(cm.exception.code, "OUT_OF_BOUNDS")

        with self.assertRaises(vis.PointValidationError) as cm:
            vis.validate_point([10.5, 2.0, 3.0], "p", existing_pts, lb_box, ub_box)
        self.assertEqual(cm.exception.code, "OUT_OF_BOUNDS")

        # 4. Duplicate ID
        with self.assertRaises(vis.PointValidationError) as cm:
            vis.validate_point([1.5, 5.0, 1.5], "z1", existing_pts, lb_box, ub_box)
        self.assertEqual(cm.exception.code, "DUPLICATE_ID")

        # 5. Duplicate Coordinates
        with self.assertRaises(vis.PointValidationError) as cm:
            vis.validate_point([4.0, 0.0, 4.0], "z_new", existing_pts, lb_box, ub_box)
        self.assertEqual(cm.exception.code, "DUPLICATE_COORDINATES")

        # 6. Dominated by existing point in N (e.g. z3 = (2,2,2) dominates (5,5,5))
        with self.assertRaises(vis.PointValidationError) as cm:
            vis.validate_point([5.0, 5.0, 5.0], "p_dom", existing_pts, lb_box, ub_box)
        self.assertEqual(cm.exception.code, "DOMINATED_POINT")

        # 7. Dominates existing point in N (e.g. (1,1,1) dominates z3 = (2,2,2))
        with self.assertRaises(vis.PointValidationError) as cm:
            vis.validate_point([1.0, 1.0, 1.0], "p_super", existing_pts, lb_box, ub_box)
        self.assertEqual(cm.exception.code, "DOMINATES_EXISTING")

        # 8. Valid nondominated point succeeds
        vis.validate_point([1.5, 5.0, 1.5], "z4", existing_pts, lb_box, ub_box)

    def test_api_rejection_of_invalid_points(self):
        """Verifies that Flask API returns HTTP 400 with descriptive JSON when adding invalid points."""
        client = web_app.app.test_client()
        client.post("/api/load_preset", json={"preset": "paper2"})

        # Dominated point rejection
        res_dom = client.post("/api/add_point", json={"id": "z_bad", "coordinates": [5.0, 5.0, 5.0]})
        self.assertEqual(res_dom.status_code, 400)
        err = json.loads(res_dom.data)
        self.assertEqual(err["code"], "DOMINATED_POINT")
        self.assertIn("dominated", err["error"].lower())

        # Out of bounds rejection
        res_oob = client.post("/api/add_point", json={"id": "z_oob", "coordinates": [-1.0, 2.0, 2.0]})
        self.assertEqual(res_oob.status_code, 400)
        err = json.loads(res_oob.data)
        self.assertEqual(err["code"], "OUT_OF_BOUNDS")

        # Wrong dimension rejection
        res_dim = client.post("/api/add_point", json={"id": "z_dim", "coordinates": [2.0, 2.0]})
        self.assertEqual(res_dim.status_code, 400)
        err = json.loads(res_dim.data)
        self.assertEqual(err["code"], "DIMENSION_MISMATCH")

    def test_maximization_problem_configuration(self):
        """Verifies that MAXIMIZE problem setting works properly with lower bounds L(N)."""
        client = web_app.app.test_client()

        # Configure MAXIMIZE problem
        res_cfg = client.post("/api/configure", json={
            "dimensions": 3,
            "sense": "MAXIMIZE",
            "lower_bound": [0.0, 0.0, 0.0],
            "upper_bound": [10.0, 10.0, 10.0]
        })
        self.assertEqual(res_cfg.status_code, 200)
        cfg_data = json.loads(res_cfg.data)
        self.assertEqual(cfg_data["sense"], "MAXIMIZE")
        self.assertEqual(cfg_data["bounds_count"], 1)

        # Initial bound for MAXIMIZE is u0 = (0, 0, 0)
        initial_row = cfg_data["table"][0]
        self.assertEqual(initial_row["u_1"], 0.0)
        self.assertEqual(initial_row["u_2"], 0.0)
        self.assertEqual(initial_row["u_3"], 0.0)

        # Insert nondominated point z1 = (6, 10, 6)
        res_add = client.post("/api/add_point", json={"id": "z1", "coordinates": [6.0, 10.0, 6.0]})
        self.assertEqual(res_add.status_code, 200)
        state_after_z1 = json.loads(res_add.data)
        self.assertEqual(state_after_z1["bounds_count"], 3)

        # Check 3D figure generation for MAXIMIZE
        res_fig3d = client.get("/api/figure_3d")
        self.assertEqual(res_fig3d.status_code, 200)

        # Check 2D figure generation for MAXIMIZE
        res_cfg2d = client.post("/api/configure", json={
            "dimensions": 2,
            "sense": "MAXIMIZE",
            "lower_bound": [0.0, 0.0],
            "upper_bound": [10.0, 10.0]
        })
        self.assertEqual(res_cfg2d.status_code, 200)
        res_fig2d = client.get("/api/figure_3d")
        self.assertEqual(res_fig2d.status_code, 200)

    def test_preset_loads_all_settings(self):
        """Verifies that loading presets completely updates sense, dimensions, LB, UB, and points."""
        client = web_app.app.test_client()

        # 1. 2D Preset
        res_2d = client.post("/api/load_preset", json={"preset": "2d"})
        self.assertEqual(res_2d.status_code, 200)
        d_2d = json.loads(res_2d.data)
        self.assertEqual(d_2d["dimensions"], 2)
        self.assertEqual(d_2d["sense"], "MINIMIZE")
        self.assertEqual(d_2d["lower_bound"], [0.0, 0.0])
        self.assertEqual(d_2d["upper_bound"], [10.0, 10.0])
        self.assertEqual(len(d_2d["points"]), 3)

        # 2. 3D Maximization Preset
        res_max = client.post("/api/load_preset", json={"preset": "max_3d"})
        self.assertEqual(res_max.status_code, 200)
        d_max = json.loads(res_max.data)
        self.assertEqual(d_max["dimensions"], 3)
        self.assertEqual(d_max["sense"], "MAXIMIZE")
        self.assertEqual(d_max["lower_bound"], [0.0, 0.0, 0.0])
        self.assertEqual(d_max["upper_bound"], [10.0, 10.0, 10.0])
        self.assertEqual(len(d_max["points"]), 3)

        # 3. Paper 2 Preset
        res_p2 = client.post("/api/load_preset", json={"preset": "paper2"})
        self.assertEqual(res_p2.status_code, 200)
        d_p2 = json.loads(res_p2.data)
        self.assertEqual(d_p2["dimensions"], 3)
        self.assertEqual(d_p2["sense"], "MINIMIZE")
        self.assertEqual(d_p2["lower_bound"], [0.0, 0.0, 0.0])
        self.assertEqual(d_p2["upper_bound"], [10.0, 10.0, 10.0])
        self.assertEqual(len(d_p2["points"]), 3)

        # 4. Paper 1 SA Preset (Example 2 / Figure 2)
        res_p1 = client.post("/api/load_preset", json={"preset": "paper1_sa"})
        self.assertEqual(res_p1.status_code, 200)
        d_p1 = json.loads(res_p1.data)
        self.assertEqual(len(d_p1["points"]), 2)
        self.assertEqual(d_p1["points"][0]["coordinates"], [3.0, 5.0, 7.0])
        self.assertEqual(d_p1["points"][1]["coordinates"], [6.0, 2.0, 4.0])
        self.assertEqual(d_p1["bounds_count"], 5)

    def test_overlapping_configuration_requests(self):
        def configure_and_read(index):
            dims = 2 + index % 2
            with web_app.app.test_client() as client:
                configured = client.post("/api/configure", json={"dimensions": dims})
                self.assertEqual(configured.status_code, 200)
                self.assertEqual(configured.get_json()["dimensions"], dims)
                response = client.get("/api/state")
                self.assertEqual(response.status_code, 200)
                data = response.get_json()
                self.assertEqual(len(data["lower_bound"]), data["dimensions"])
                self.assertEqual(len(data["upper_bound"]), data["dimensions"])
                self.assertEqual(len(data["table"][0]["Coordinates"]), data["dimensions"])
        with ThreadPoolExecutor(max_workers=4) as executor:
            list(executor.map(configure_and_read, range(24)))

    def test_malformed_requests_are_transactional(self):
        client = web_app.app.test_client()
        client.post("/api/load_preset", json={"preset": "paper2"})
        before = client.get("/api/state").get_json()
        invalid = [
            ("configure", {"dimensions": "oops"}),
            ("configure", {"dimensions": 2.5}),
            ("configure", {"dimensions": True}),
            ("configure", {"sense": "TYPO"}),
            ("configure", {"sense": []}),
            ("configure", {"lower_bound": ["nan", 0, 0], "upper_bound": [10, 10, 10]}),
            ("configure", {"lower_bound": [0, 0, 0]}),
            ("configure", {"lower_bound": "000", "upper_bound": [10, 10, 10]}),
            ("configure", {"lower_bound": [True, 0, 0], "upper_bound": [10, 10, 10]}),
            ("add_point", {"id": 42, "coordinates": [1, 5, 1]}),
            ("add_point", {"coordinates": "123"}),
            ("add_point", {"coordinates": [float("inf"), 5, 1]}),
            ("add_point", {"coordinates": [-5e-8, 5, 1]}),
            ("add_point", {"coordinates": [10, 5, 1]}),
            ("delete_point", {"index": "oops"}),
            ("delete_point", {"index": 1.5}),
            ("delete_point", {"index": -1}),
            ("load_preset", {"preset": "unknown"}),
            ("load_preset", {"preset": []}),
        ]
        for route, payload in invalid:
            with self.subTest(route=route, payload=payload):
                response = client.post("/api/" + route, json=payload)
                self.assertEqual(response.status_code, 400)
                self.assertIn("code", response.get_json())
                self.assertEqual(client.get("/api/state").get_json(), before)
        for body in ('[]', 'null', '"text"', '{broken'):
            response = client.post("/api/configure", data=body, content_type="application/json")
            self.assertEqual(response.status_code, 400)
            self.assertIn("code", response.get_json())
        for query in ("/api/state?step=oops", "/api/graph?mode=unknown", "/api/figure_3d?dom_opacity=nan"):
            self.assertEqual(client.get(query).status_code, 400)

    def test_exact_boundaries_and_native_errors(self):
        for sense in (lb.Objective.MINIMIZE, lb.Objective.MAXIMIZE):
            is_min = sense == lb.Objective.MINIMIZE
            valid = [0, 5] if is_min else [10, 5]
            invalid = ([-5e-8, 5], [10, 5]) if is_min else ([10+5e-8, 5], [0, 5])
            ref, anti = ([10, 10], [0, 0]) if is_min else ([0, 0], [10, 10])
            nbs = lb.NeighborhoodBoundSet(ref, anti, sense=sense)
            for coords in invalid:
                with self.assertRaises(vis.PointValidationError):
                    vis.validate_point(coords, "z", [], [0, 0], [10, 10], sense)
                with self.assertRaises(ValueError):
                    nbs.update(lb.Point("z", coords))
            vis.validate_point(valid, "z", [], [0, 0], [10, 10], sense)
            self.assertTrue(nbs.update(lb.Point("z", valid)))
            with self.assertRaises(ValueError):
                nbs.find_containing_bound([1])
        missing = lb.BoundSet([10, 10])
        with self.assertRaises(RuntimeError):
            missing.update_ra(lb.Point("z", [3, 7]))
        with self.assertRaises(ValueError):
            lb.BoundSet([10, 10], sense="TYPO")
        tracker = vis.LocalBoundsTracker([10, 10], [0, 0], sense="MINIMIZE")
        with self.assertRaises(vis.PointValidationError):
            tracker.add_point(lb.Point("z", [-5e-8, 5]))
        self.assertEqual(len(tracker.steps), 1)
        self.assertEqual(tracker.points_history, [])

    def test_ngp_graph_uses_full_contracted_adjacency(self):
        ref, anti = [10]*3, [0]*3
        nbs = lb.NeighborhoodBoundSet(ref, anti)
        points = [lb.Point(f"z{i}", coords) for i, coords in enumerate([[4,0,4], [4,3,1], [2,3,2]])]
        for point in points:
            nbs.update(point)
        graph = nbs.get_adjacency_graph()
        expected = {(min(u,v), max(u,v)) for u, neighbors in enumerate(graph.adjacency_list)
                    for v in neighbors if u != v}
        data = vis.extract_bounds_data(nbs, ref, anti, points)
        self.assertEqual(len(expected), 9)
        self.assertEqual(set(data["edges_undirected"]), expected)
        elements = vis.extract_network_graph_elements(data, mode="undirected")
        self.assertEqual({(e["from"], e["to"]) for e in elements["edges"]}, expected)

    def test_2d_highlight_and_escaped_report(self):
        from visualization.cli import generate_html_report
        payload = '<img src=x onerror="alert(1)"></script>'
        nbs = lb.NeighborhoodBoundSet([10,10], [0,0])
        point = lb.Point(payload, [3,7])
        nbs.update(point)
        data = vis.extract_bounds_data(nbs, [10,10], [0,0], [point])
        highlighted = data["bounds"][0].id
        fig = vis.plot_2d_bounds(data, highlight_bound_id=highlighted)
        bounds = next(trace for trace in fig.data if trace.mode == "markers" and trace.legendgroup == "local_bounds")
        self.assertEqual(bounds.marker.color[0], "#f97316")
        self.assertEqual(bounds.customdata[0], highlighted)
        self.assertNotIn(payload, fig.to_json())
        with tempfile.TemporaryDirectory() as directory, contextlib.redirect_stdout(io.StringIO()):
            output = Path(directory) / "report.html"
            generate_html_report(data, str(output))
            html = output.read_text()
            self.assertNotIn(payload, html)
            self.assertIn("&lt;img", html)
            self.assertNotIn('src="https://cdn.plot.ly', html)
        maxdata = dict(data, sense="MAXIMIZE")
        self.assertTrue(any(trace.name == "Local Bounds L(N)" for trace in vis.plot_2d_bounds(maxdata).data))

    def test_raw_graph_preserves_quasi_nodes(self):
        for sense in (lb.Objective.MINIMIZE, lb.Objective.MAXIMIZE):
            is_min = sense == lb.Objective.MINIMIZE
            ref, anti = ([6]*3, [0]*3) if is_min else ([0]*3, [6]*3)
            nbs = lb.NeighborhoodBoundSet(ref, anti, sense=sense)
            for i, z in enumerate([[5,2,5], [2,4,3], [4,4,1], [4,3,2], [3,3,3]]):
                nbs.update(lb.Point(f"z{i}", z if is_min else [6-c for c in z]))
            graph = nbs.get_adjacency_graph(include_quasi=True)
            self.assertEqual(len(graph.nodes), nbs.size())
            self.assertEqual(sum(not q for q in graph.quasi), 8)
            self.assertGreater(sum(graph.quasi), 0)
            for u, neighbors in enumerate(graph.adjacency_list):
                for v in neighbors:
                    self.assertIn(u, graph.adjacency_list[v])
            data = vis.extract_bounds_data(nbs, ref, anti, sense=sense, include_quasi=True)
            self.assertEqual(data["counts"]["nonredundant"], 8)
            self.assertEqual(data["counts"]["total"], len(data["bounds"]))
            elements = vis.extract_network_graph_elements(data)
            self.assertEqual(sum(n["is_quasi"] for n in elements["nodes"]), sum(graph.quasi))
            self.assertTrue(all(n["shape"] == "ellipse" for n in elements["nodes"] if n["is_quasi"]))
        client = web_app.app.test_client()
        client.post("/api/load_preset", json={"preset": "paper1_ngp"})
        raw = client.get("/api/graph?view=raw").get_json()
        clean = client.get("/api/graph").get_json()
        self.assertEqual(raw["view"], "raw")
        self.assertEqual(raw["counts"], clean["counts"])
        self.assertEqual(len(raw["nodes"]), raw["counts"]["total"])
        self.assertEqual(len(clean["nodes"]), raw["counts"]["nonredundant"])
        self.assertEqual(client.get("/api/graph?view=invalid").status_code, 400)

    def test_membership_probe_exact_boundaries_and_no_mutation(self):
        client = web_app.app.test_client()
        for sense in ("MINIMIZE", "MAXIMIZE"):
            client.post("/api/configure", json={"dimensions": 2, "sense": sense})
            is_min = sense == "MINIMIZE"
            point = [3,7] if is_min else [7,3]
            client.post("/api/add_point", json={"id": "z", "coordinates": point})
            before = client.get("/api/state").get_json()
            outside = [10,5] if is_min else [0,5]
            anti = [0,5] if is_min else [10,5]
            for coords, domain, search in [(point, True, False), (outside, False, False), (anti, True, True)]:
                result = client.post("/api/probe", json={"coordinates": coords})
                self.assertEqual(result.status_code, 200)
                result = result.get_json()
                self.assertEqual(result["in_domain"], domain)
                self.assertEqual(result["in_search_region"], search)
                self.assertEqual(result["containing_bound"] is not None, search)
                self.assertIn("<", "; ".join(result["domain"]))
                if search:
                    self.assertIn(result["containing_bound"]["id"], {r["ID"] for r in before["table"]})
            initial = client.post("/api/probe?step=0", json={"coordinates": point}).get_json()
            self.assertTrue(initial["in_search_region"])
            for coords in ([1], [float("nan"), 1], [True,1], "1,1"):
                self.assertEqual(client.post("/api/probe", json={"coordinates": coords}).status_code, 400)
            self.assertEqual(client.get("/api/state").get_json(), before)

    def test_history_matches_full_tracker_without_graph_exports(self):
        from visualization.presets import PRESETS
        for preset in PRESETS:
            state = web_app.AppState()
            state.load_preset(preset)
            tracker = vis.LocalBoundsTracker(state.reference_point, state.anti_reference,
                                            sense=state.sense)
            for point in state.points:
                tracker.add_point(point)
            expected = [{"step": step.step,
                         "point_added": step.point_added.id if step.point_added else "Initial State",
                         "destroyed_bounds": step.destroyed_bounds,
                         "new_bounds": step.new_bounds,
                         "total_bounds": len(step.bounds_data["bounds"])}
                        for step in tracker.steps]
            with patch.object(vis, "extract_bounds_data", side_effect=AssertionError("unused graph export")):
                self.assertEqual(state.compute_step_history(), expected)
            state.compute_step_history()[0]["new_bounds"].clear()
            self.assertEqual(state.compute_step_history(), expected)

    def test_step_cache_isolation_keys_and_mutation_invalidation(self):
        state = web_app.AppState()
        state.load_preset("paper1_ngp")
        with patch.object(vis, "extract_bounds_data", wraps=vis.extract_bounds_data) as extract:
            first = state.get_bounds_at_step(3)
            expected_table = vis.create_bounds_table(first)
            first["bounds"][0].coordinates[0] = -99
            first["points"][0]["coordinates"][0] = -99
            again = state.get_bounds_at_step(3)
            self.assertEqual(vis.create_bounds_table(again), expected_table)
            self.assertNotEqual(again["points"][0]["coordinates"][0], -99)
            self.assertEqual(extract.call_count, 1)
            state.get_bounds_at_step(3, include_quasi=True)
            state.get_bounds_at_step(0)
            self.assertEqual(extract.call_count, 3)
            # Mutable native point fields must be part of the cache key.
            state.points[0].id = "renamed"
            self.assertEqual(state.get_bounds_at_step(3)["points"][0]["id"], "renamed")
            self.assertEqual(extract.call_count, 4)
            state.load_preset("paper1_ngp")
            state.get_bounds_at_step(3)
            self.assertEqual(extract.call_count, 5)
            state.reset()
            self.assertEqual(state.get_bounds_at_step(0)["counts"]["total"], 1)
            self.assertEqual(extract.call_count, 6)

    def test_step_cache_keeps_native_graph_argument_errors(self):
        state = web_app.AppState()
        state.get_bounds_at_step(0, False)
        state.get_bounds_at_step(0, True)
        for invalid in (0, 1, 0.0, 1.0, [], {}, None, "raw"):
            with self.subTest(include_quasi=invalid):
                with self.assertRaisesRegex(TypeError, "get_adjacency_graph"):
                    state.get_bounds_at_step(0, invalid)

    def test_step_cache_preserves_slice_errors_and_signed_zero(self):
        state = web_app.AppState()
        state.load_preset("2d")
        state.get_bounds_at_step(0)
        with self.assertRaisesRegex(TypeError, "slice indices must be integers"):
            state.get_bounds_at_step(0.0)
        state.lower_bound[0] = -0.0
        self.assertEqual(str(state.get_bounds_at_step(0)["lower_bound"][0]), "-0.0")
        state.lower_bound[0] = 0.0
        self.assertEqual(str(state.get_bounds_at_step(0)["lower_bound"][0]), "0.0")
        state.points[0].coordinates = [3.5, 7]
        self.assertEqual(state.get_bounds_at_step(1)["points"][0]["coordinates"], [3.5, 7])

    def test_figure_cache_preserves_options_and_validation(self):
        client = web_app.app.test_client()
        client.post("/api/load_preset", json={"preset": "paper2"})
        with patch.object(vis, "plot_3d_bounds", wraps=vis.plot_3d_bounds) as spatial:
            original = client.get("/api/figure_3d?step=2").data
            self.assertEqual(client.get("/api/figure_3d?step=2").data, original)
            self.assertEqual(spatial.call_count, 1)
            changed = client.get("/api/figure_3d?step=2&show_boxes=1&dom_opacity=0.5").data
            self.assertNotEqual(changed, original)
            self.assertEqual(spatial.call_count, 2)
            self.assertEqual(client.get("/api/figure_3d?step=2&show_boxes=1&dom_opacity=0.5").data, changed)
            self.assertEqual(spatial.call_count, 2)
            for suffix in ("step=oops", "dom_opacity=nan", "dom_opacity=inf"):
                self.assertEqual(client.get("/api/figure_3d?" + suffix).status_code, 400)
            client.post("/api/delete_point", json={"index": 2})
            self.assertEqual(client.get("/api/figure_3d?step=2").data, original)
            self.assertEqual(spatial.call_count, 3)
        with patch.object(vis, "plot_pairwise_projections_2d", wraps=vis.plot_pairwise_projections_2d) as pairwise:
            original = client.get("/api/figure_pairwise").data
            self.assertEqual(client.get("/api/figure_pairwise").data, original)
            self.assertEqual(pairwise.call_count, 1)
            client.post("/api/configure", json={"dimensions": 2, "sense": "MAXIMIZE"})
            self.assertNotEqual(client.get("/api/figure_pairwise").data, original)
            self.assertEqual(pairwise.call_count, 2)

    def test_shared_presets_and_custom_cli_dimensions(self):
        from visualization import cli
        from visualization.presets import PRESETS
        client = web_app.app.test_client()
        for preset, (sense, coordinates) in PRESETS.items():
            response = client.post("/api/load_preset", json={"preset": preset})
            self.assertEqual(response.status_code, 200)
            self.assertEqual(response.get_json()["preset"], preset)
            self.assertEqual([p["coordinates"] for p in response.get_json()["points"]], coordinates)
            with patch("sys.argv", ["cli", "--preset", preset]), patch.object(cli, "generate_html_report") as report, contextlib.redirect_stdout(io.StringIO()):
                cli.main()
            cli_data = report.call_args.args[0]
            self.assertEqual([p["coordinates"] for p in cli_data["points"]], coordinates)
            self.assertEqual(cli_data["sense"], sense)
        client.post("/api/load_preset", json={"preset": "2d"})
        client.post("/api/delete_point", json={"index": 0})
        custom = client.post("/api/add_point", json={"coordinates": [1,9]}).get_json()
        self.assertIsNone(custom["preset"])
        self.assertEqual(custom["points"][-1]["id"], "z4")
        with patch("sys.argv", ["cli", "--points", "[[3,7]]", "--ref", "20,20", "--anti", "0,0"]), patch.object(cli, "generate_html_report") as report, contextlib.redirect_stdout(io.StringIO()):
            cli.main()
        self.assertEqual(report.call_args.args[0]["dimensions"], 2)
        self.assertEqual(report.call_args.args[0]["reference_point"], [20,20])
        with patch("sys.argv", ["cli", "--points", "[]"]), patch.object(cli, "generate_html_report") as report, contextlib.redirect_stdout(io.StringIO()):
            cli.main()
        self.assertEqual(report.call_args.args[0]["points"], [])
        with patch("sys.argv", ["cli", "--points", "[1]"]), contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as error:
            cli.main()
        self.assertEqual(error.exception.code, 2)

    def test_readability_labels_and_wireframe_preserve_geometry(self):
        state = web_app.AppState()
        data = state.get_bounds_at_step(2)
        selected = data["bounds"][0].id
        solid = vis.plot_3d_bounds(data, show_search_boxes=True, highlight_bound_id=selected)
        meshes = [t for t in solid.data if t.type == "mesh3d" and t.legendgroup == "dominated_zones"]
        self.assertEqual(len(meshes), len(data["points"]))
        self.assertTrue(all(t.opacity == 0.30 for t in meshes))
        reference_labels = [label for label in solid.layout.scene.annotations if label.name == "reference_points"]
        self.assertEqual([label.text for label in reference_labels], ["m", "M"])
        self.assertEqual([label.x for label in reference_labels], [data["lower_bound"][0], data["upper_bound"][0]])
        wire = vis.plot_3d_bounds(data, show_search_boxes=True, highlight_bound_id=selected, wireframe=True)
        self.assertFalse(any(t.type == "mesh3d" for t in wire.data))
        solid_edges = [(t.legendgroup, tuple(t.x), tuple(t.y), tuple(t.z))
                       for t in solid.data if t.type == "scatter3d" and t.mode == "lines"]
        wire_edges = [(t.legendgroup, tuple(t.x), tuple(t.y), tuple(t.z))
                      for t in wire.data if t.type == "scatter3d" and t.mode == "lines"]
        self.assertEqual(wire_edges, solid_edges)
        self.assertEqual({t.legendgroup for t in solid.data if t.showlegend},
                         {t.legendgroup for t in wire.data if t.showlegend})
        for plot, plot_data in ((vis.plot_3d_bounds, data), (vis.plot_pairwise_projections_2d, data),
                                (vis.plot_2d_bounds, dict(data, dimensions=2))):
            for mode in ("all", "none"):
                figure = plot(plot_data, highlight_bound_id=selected, label_mode=mode)
                annotations = figure.layout.scene.annotations if plot == vis.plot_3d_bounds else figure.layout.annotations
                labels = [label for label in annotations if label.name]
                self.assertEqual(bool(labels), mode == "all")
                for label in labels:
                    self.assertEqual(label.bgcolor, "#ffffff")
                    self.assertEqual(label.font.color, "#0f172a")
                    self.assertFalse(label.captureevents)
                    self.assertEqual(label.font.size, 12)
                markers = [trace for trace in figure.data if getattr(trace, "mode", None) == "markers"]
                self.assertTrue(markers)
                for trace in markers:
                    self.assertTrue(all(trace.hovertext))
                    self.assertTrue(all(bool(text) == (mode == "all") for text in trace.text))

    def test_pairwise_coincident_markers_keep_all_members_and_boxes(self):
        records = [vis.BoundRecord(i, name, coords, [], [None] * 3)
                   for i, (name, coords) in enumerate((("u0", [4,4,2]), ("u1", [4,4,7]), ("u2", [8,3,1])))]
        data = dict(dimensions=3, sense="MINIMIZE", lower_bound=[0]*3, upper_bound=[10]*3,
                    reference_point=[10]*3, anti_reference=[0]*3, bounds=records,
                    points=[{"id": "z0", "coordinates": [2,2,5]}, {"id": "z1", "coordinates": [2,2,7]}])
        fig = vis.plot_pairwise_projections_2d(data, highlight_bound_id="u1", label_mode="all")
        bounds = next(t for t in fig.data if t.legendgroup == "local_bounds" and t.mode == "markers")
        self.assertEqual(list(bounds.customdata), [["u0", "u1"], "u2"])
        self.assertEqual(list(bounds.x), [4,8])
        self.assertEqual(bounds.text[0], "u0, u1")
        self.assertEqual(bounds.textfont.color, "#ffffff")
        self.assertTrue(any(a.text == "u0, u1" and a.xref == "x" for a in fig.layout.annotations))
        self.assertEqual(bounds.marker.color[0], "#f97316")
        self.assertIn("u0, u1", bounds.hovertext[0])
        points = next(t for t in fig.data if t.legendgroup == "points")
        self.assertEqual(list(points.customdata), [["z0", "z1"]])
        self.assertEqual(points.text[0], "z0, z1")
        self.assertIn("z0, z1", points.hovertext[0])
        self.assertEqual(len([t for t in fig.data if t.mode == "lines"]), len(records) * 3)
        self.assertEqual([b.coordinates for b in data["bounds"]], [[4,4,2], [4,4,7], [8,3,1]])
        all_labels = vis.plot_pairwise_projections_2d(data, label_mode="all")
        all_bounds = next(t for t in all_labels.data if t.legendgroup == "local_bounds" and t.mode == "markers")
        self.assertEqual(all_bounds.text[0], "u0, u1")

    def test_plot_ranges_and_view_revision_follow_problem_only(self):
        state = web_app.AppState()
        initial, later = state.get_bounds_at_step(0), state.get_bounds_at_step(2)
        for plot in (vis.plot_3d_bounds, vis.plot_2d_bounds, vis.plot_pairwise_projections_2d):
            a = plot(initial)
            b = plot(later, highlight_bound_id=later["bounds"][0].id, label_mode="all")
            self.assertEqual(a.layout.uirevision, b.layout.uirevision)
            self.assertNotEqual(a.layout.uirevision, plot(dict(initial, sense="MAXIMIZE")).layout.uirevision)
            self.assertNotEqual(a.layout.uirevision, plot(dict(initial, lower_bound=[-1]*3)).layout.uirevision)
        data = dict(initial, lower_bound=[-2,-4,1], upper_bound=[8,16,21])
        fig = vis.plot_pairwise_projections_2d(data)
        for suffix, (di, dj) in zip(("", "2", "3"), ((0,1), (0,2), (1,2))):
            self.assertEqual(list(fig.layout["xaxis" + suffix].range), [data["lower_bound"][di]-0.5, data["upper_bound"][di]+0.5])
            self.assertEqual(list(fig.layout["yaxis" + suffix].range), [data["lower_bound"][dj]-0.5, data["upper_bound"][dj]+0.5])
        two = vis.plot_2d_bounds(data)
        self.assertEqual(tuple(two.layout.xaxis.range), (-2.5,8.5))
        self.assertEqual(tuple(two.layout.yaxis.range), (-4.5,16.5))

    def test_readability_api_options_validation_and_cache_keys(self):
        client = web_app.app.test_client()
        client.post("/api/load_preset", json={"preset": "paper2"})
        selected = client.get("/api/state").get_json()["table"][0]["ID"]
        for route, plot_name in (("figure_3d", "plot_3d_bounds"), ("figure_pairwise", "plot_pairwise_projections_2d")):
            with patch.object(vis, plot_name, wraps=getattr(vis, plot_name)) as plot:
                original = client.get(f"/api/{route}?highlight={selected}").get_json()
                self.assertEqual(client.get(f"/api/{route}?highlight={selected}").get_json(), original)
                self.assertEqual(plot.call_count, 1)
                for mode in ("none",):
                    changed = client.get(f"/api/{route}?labels={mode}&highlight={selected}").get_json()
                    self.assertNotEqual(changed, original)
                    self.assertEqual(client.get(f"/api/{route}?labels={mode}&highlight={selected}").get_json(), changed)
                self.assertEqual(plot.call_count, 2)
                for invalid in ("focus", "bad", "", "ALL"):
                    response = client.get(f"/api/{route}?labels={invalid}")
                    self.assertEqual(response.status_code, 400)
                    self.assertEqual(response.get_json()["code"], "INVALID_LABELS")
                self.assertEqual(plot.call_count, 2)
        for route in ("figure_3d", "figure_pairwise"):
            default = client.get(f"/api/{route}").get_json()
            self.assertEqual(default, client.get(f"/api/{route}?labels=all").get_json())
            self.assertTrue(all(all(t["text"]) and t["textfont"]["size"] == 12
                                for t in default["data"] if t.get("mode") == "markers"))
        solid = client.get("/api/figure_3d").get_json()
        self.assertTrue(all(t["opacity"] == 0.30 for t in solid["data"] if t["type"] == "mesh3d" and t["legendgroup"] == "dominated_zones"))
        with patch.object(vis, "plot_3d_bounds", wraps=vis.plot_3d_bounds) as plot:
            wire = client.get("/api/figure_3d?wireframe=1").get_json()
            self.assertFalse(any(t["type"] == "mesh3d" for t in wire["data"]))
            self.assertEqual(client.get("/api/figure_3d?wireframe=1").get_json(), wire)
            self.assertEqual(plot.call_count, 1)
            self.assertEqual(client.get("/api/figure_3d").get_json(), solid)
        client.post("/api/load_preset", json={"preset": "2d"})
        figure = client.get("/api/figure_3d?labels=all").get_json()
        self.assertTrue(all(t["textfont"]["color"] == "#ffffff" and all(t["text"])
                            for t in figure["data"] if t.get("mode") == "markers"))

    def test_final_cli_extracts_only_final_state(self):
        from visualization import cli
        for args, count in (([], 7), (["--preset", "max_3d"], 7), (["--points", "[]"], 1)):
            results = []
            for step_by_step in (False, True):
                argv = ["cli", *args, *(["--step-by-step"] if step_by_step else [])]
                stdout = io.StringIO()
                with patch("sys.argv", argv), patch.object(cli, "generate_html_report") as report, \
                        patch.object(vis, "extract_bounds_data", wraps=vis.extract_bounds_data) as extract, \
                        contextlib.redirect_stdout(stdout):
                    cli.main()
                data = report.call_args.args[0]
                assert len(data["bounds"]) == count
                assert extract.call_count == (len(data["points"]) + 1 if step_by_step else 1)
                assert ("--- Step" in stdout.getvalue()) == (step_by_step and bool(data["points"]))
                results.append(dict(data, bounds=[bound.to_dict() for bound in data["bounds"]]))
            assert results[0] == results[1]

    def test_final_cli_preserves_validation_and_point_construction_order(self):
        from visualization import cli
        for args, message in ((["--points", "[[10,5],[\"bad\",2]]"], "incompatible function arguments"),
                              (["--points", "[[3,7],[3,7]]"], "Duplicate coordinates"),
                              (["--points", "[[3,7]]", "--anti", "10,10"], "Anti-reference must be strictly better")):
            errors = []
            for step_by_step in (False, True):
                stderr = io.StringIO()
                with patch("sys.argv", ["cli", *args, *(["--step-by-step"] if step_by_step else [])]), \
                        contextlib.redirect_stderr(stderr), self.assertRaises(SystemExit) as error:
                    cli.main()
                assert error.exception.code == 2
                errors.append(stderr.getvalue())
            assert errors[0] == errors[1]
            assert message in errors[0]

    def test_tracker_snapshots_survive_new_points_and_replay_mutations(self):
        def snapshot(step):
            return json.dumps(step.bounds_data, default=vis.BoundRecord.to_dict, sort_keys=True)

        tracker = vis.LocalBoundsTracker([10, 10], [0, 0])
        point = lb.Point("z1", [3, 7])
        first = tracker.add_point(point)
        first_snapshot = snapshot(first)
        point.id = "renamed"
        point.coordinates = [2, 8]
        second = tracker.add_point(lb.Point("z2", [7, 2]))
        assert tracker.get_step(1) is first
        assert snapshot(first) == first_snapshot
        assert second.bounds_data["points"][0] == {"id": "renamed", "coordinates": [2, 8]}
        assert [2, 10] in [bound.coordinates for bound in second.bounds_data["bounds"]]
        before = snapshot(second)
        with self.assertRaises(vis.PointValidationError):
            tracker.add_point(lb.Point("bad", [10, 5]))
        assert len(tracker.steps) == 3
        assert snapshot(tracker.get_step(2)) == before


if __name__ == "__main__":
    unittest.main()
