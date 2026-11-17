"""Tests for the visualization engine, the dashboard API and the report CLI."""

import contextlib
import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import local_bounds as lb
import local_bounds.visualization.app as web_app
import local_bounds.visualization.engine as vis
from local_bounds.visualization import cli
from local_bounds.visualization.presets import PRESETS

MIN, MAX = lb.Objective.MINIMIZE, lb.Objective.MAXIMIZE

# Dächert et al. (2017), Example 2.8: bounds destroyed and created by each insertion.
PAPER2_HISTORY = [
    ("Initial State", [], ["u0"]),
    ("z1", ["u0"], ["u01", "u02", "u03"]),
    ("z2", ["u01", "u03"], ["u011", "u012", "u032", "u033"]),
    ("z3", ["u011", "u012", "u032"], ["u0111", "u0113", "u0122", "u0322", "u0323"]),
]
PAPER2_POINTS = [[4, 0, 4], [3, 3, 1], [2, 2, 2]]
TIED_POINTS = [[5, 2, 5], [2, 4, 3], [4, 4, 1], [4, 3, 2], [3, 3, 3]]  # 8 nonredundant bounds


def build(coordinates, sense=MIN, **options):
    """Inserts z1, z2, ... into the box [0, 10]^p; returns the bound set and its extracted data."""
    dims = len(coordinates[0])
    ref, anti = ([10] * dims, [0] * dims) if sense == MIN else ([0] * dims, [10] * dims)
    nbs = lb.NeighborhoodBoundSet(ref, anti, sense=sense)
    points = [lb.Point(f"z{i+1}", c) for i, c in enumerate(coordinates)]
    for point in points:
        nbs.update(point)
    return nbs, vis.extract_bounds_data(nbs, ref, anti, points, sense=sense, **options)


def defining_ids(data):
    return {tuple(b.coordinates): [dp["id"] for dp in b.defining_points] for b in data["bounds"]}


class TestEngine(unittest.TestCase):

    def test_table_matches_paper_example_2_8(self):
        # Dächert et al. (2017), Example 2.8 (3): defining points and extreme bounds.
        _, data = build(PAPER2_POINTS)
        rows = {row["ID"]: [row[f"z^{j}(u)"] for j in (1, 2, 3)] for row in vis.create_bounds_table(data)}
        self.assertEqual(rows, {
            "u0111": ["z3", "z^2", "z^3"],
            "u0113": ["z2", "z^2", "z3"],
            "u0122": ["z1", "z3", "z^3"],
            "u02": ["z^1", "z1", "z^3"],
            "u0322": ["z^1", "z3", "z1"],
            "u0323": ["z^1", "z2", "z3"],
            "u033": ["z^1", "z^2", "z2"],
        })
        self.assertEqual({b.id: b.extreme_dim for b in data["bounds"] if b.is_extreme},
                         {"u0111": 0, "u02": 1, "u033": 2})

    def test_neighbors_match_paper_example_4_2(self):
        # After z1 alone: u1 has no 1-neighbor, its 2-neighbor is u2 and its 3-neighbor is u3.
        _, data = build(PAPER2_POINTS[:1])
        ids = [b.id for b in data["bounds"]]
        u1 = data["bounds"][ids.index("u01")]
        self.assertEqual([n if n is None else ids[n] for n in u1.neighbors], [None, "u02", "u03"])

    def test_defining_points_match_paper_example_4_13(self):
        # Dächert et al. (2017), Example 4.13: the table of defining points under ties.
        _, data = build([[4, 0, 4], [4, 3, 1], [2, 3, 2]], include_quasi=True)
        self.assertEqual(defining_ids(data), {
            (2, 10, 10): ["z3", "z^2", "z^3"],
            (4, 3, 10): ["z1", "z3", "z^3"],
            (10, 0, 10): ["z^1", "z1", "z^3"],
            (4, 3, 4): ["z2", "z3", "z1"],
            (4, 10, 2): ["z2", "z^2", "z3"],
            (10, 3, 4): ["z^1", "z2", "z1"],
            (10, 10, 1): ["z^1", "z^2", "z2"],
        })

    def test_defining_points_are_derived_when_not_stored(self):
        # RE keeps no defining points; the table then shows a point that satisfies
        # z_j = u_j and z_-j < u_-j (Klamroth et al. (2015), Definition 4.2).
        bounds = lb.BoundSet([10] * 3)
        points = [lb.Point("z1", [3, 5, 7]), lb.Point("z2", [6, 2, 4])]
        for point in points:
            bounds.update_re(point)
        data = vis.extract_bounds_data(bounds, [10] * 3, [0] * 3, points)
        self.assertEqual(defining_ids(data)[(6, 5, 10)], ["z2", "z1", "z^3"])

    def test_2d_bounds_form_the_chain_of_paper_1_example_1(self):
        # Klamroth et al. (2015), Example 1: consecutive points define the bounds.
        _, data = build([[3, 7], [5, 4], [7, 2]])
        coords = [tuple(b.coordinates) for b in data["bounds"]]
        self.assertEqual(sorted(coords), [(3, 10), (5, 7), (7, 4), (10, 2)])
        self.assertEqual({c for c, b in zip(coords, data["bounds"]) if b.is_extreme}, {(3, 10), (10, 2)})
        self.assertEqual({frozenset((coords[u], coords[v])) for u, v in data["edges_undirected"]},
                         {frozenset(edge) for edge in (((3, 10), (5, 7)), ((5, 7), (7, 4)), ((7, 4), (10, 2)))})

    def test_graph_edges_under_ties_come_from_the_contracted_adjacency(self):
        nbs, data = build([[4, 0, 4], [4, 3, 1], [2, 3, 2]])
        graph = nbs.get_adjacency_graph()
        expected = {(min(u, v), max(u, v)) for u, neighbors in enumerate(graph.adjacency_list)
                    for v in neighbors if u != v}
        self.assertEqual(len(expected), 9)
        self.assertEqual(set(data["edges_undirected"]), expected)
        elements = vis.extract_network_graph_elements(data, mode="undirected")
        self.assertEqual({(e["from"], e["to"]) for e in elements["edges"]}, expected)

    def test_raw_view_keeps_quasi_bounds(self):
        for sense in (MIN, MAX):
            coordinates = TIED_POINTS if sense == MIN else [[10 - c for c in z] for z in TIED_POINTS]
            nbs, data = build(coordinates, sense, include_quasi=True)
            with self.subTest(sense=sense):
                self.assertGreater(nbs.size(), 8)
                self.assertEqual(data["counts"], {"total": nbs.size(), "nonredundant": 8, "quasi": nbs.size() - 8})
                self.assertEqual(len(data["bounds"]), nbs.size())
                nodes = vis.extract_network_graph_elements(data)["nodes"]
                self.assertEqual(sum(n["is_quasi"] for n in nodes), nbs.size() - 8)

    def test_point_validation_rules(self):
        existing = [lb.Point(f"z{i+1}", c) for i, c in enumerate(PAPER2_POINTS)]
        invalid = [
            ([1, 2], "p", "DIMENSION_MISMATCH"),
            ([1, 2, 3, 4], "p", "DIMENSION_MISMATCH"),
            ([float("nan"), 2, 3], "p", "NON_FINITE"),
            ([float("inf"), 2, 3], "p", "NON_FINITE"),
            ([-0.5, 2, 3], "p", "OUT_OF_BOUNDS"),
            ([10.5, 2, 3], "p", "OUT_OF_BOUNDS"),
            ([1.5, 5, 1.5], "z1", "DUPLICATE_ID"),
            ([4, 0, 4], "z_new", "DUPLICATE_COORDINATES"),
            ([5, 5, 5], "p_dom", "DOMINATED_POINT"),
            ([1, 1, 1], "p_super", "DOMINATES_EXISTING"),
        ]
        for coords, point_id, code in invalid:
            with self.subTest(coords=coords, point_id=point_id):
                with self.assertRaises(vis.PointValidationError) as error:
                    vis.validate_point(coords, point_id, existing, [0] * 3, [10] * 3)
                self.assertEqual(error.exception.code, code)
        vis.validate_point([1.5, 5, 1.5], "z4", existing, [0] * 3, [10] * 3)

    def test_validation_agrees_with_the_library_at_interval_boundaries(self):
        # The anti-reference boundary is inside the interval; the reference boundary is not.
        for sense, valid, invalid in ((MIN, [0, 5], ([-5e-8, 5], [10, 5])),
                                      (MAX, [10, 5], ([10 + 5e-8, 5], [0, 5]))):
            ref, anti = ([10, 10], [0, 0]) if sense == MIN else ([0, 0], [10, 10])
            nbs = lb.NeighborhoodBoundSet(ref, anti, sense=sense)
            for coords in invalid:
                with self.subTest(sense=sense, coords=coords):
                    with self.assertRaises(vis.PointValidationError):
                        vis.validate_point(coords, "z", [], [0, 0], [10, 10], sense)
                    with self.assertRaises(ValueError):
                        nbs.update(lb.Point("z", coords))
            vis.validate_point(valid, "z", [], [0, 0], [10, 10], sense)
            self.assertTrue(nbs.update(lb.Point("z", valid))[0])

    def test_tracker_history_matches_paper_example_2_8(self):
        tracker = vis.LocalBoundsTracker([10] * 3, [0] * 3)
        for i, coords in enumerate(PAPER2_POINTS):
            tracker.add_point(lb.Point(f"z{i+1}", coords))
        history = [(s.point_added.id if s.point_added else "Initial State", s.destroyed_bounds, s.new_bounds)
                   for s in tracker.steps]
        self.assertEqual(history, PAPER2_HISTORY)

    def test_tracker_keeps_its_history_when_a_point_is_rejected(self):
        tracker = vis.LocalBoundsTracker([10, 10], [0, 0])
        tracker.add_point(lb.Point("z1", [3, 7]))
        with self.assertRaises(vis.PointValidationError):
            tracker.add_point(lb.Point("bad", [10, 5]))
        self.assertEqual(len(tracker.steps), 2)
        self.assertEqual([p.id for p in tracker.points_history], ["z1"])


class TestFigures(unittest.TestCase):

    def test_every_plot_builds_for_both_senses(self):
        for sense in (MIN, MAX):
            _, data3 = build([[4, 2, 5]], sense)
            _, data2 = build([[3, 5]], sense)
            for plot, data in ((vis.plot_3d_bounds, data3), (vis.plot_pairwise_projections_2d, data3),
                               (vis.plot_neighbor_graph_plotly, data3), (vis.plot_2d_bounds, data2)):
                with self.subTest(plot=plot.__name__, sense=sense):
                    self.assertGreater(len(plot(data).data), 0)

    def test_highlight_draws_the_selected_search_zone(self):
        _, data = build([[4, 2, 5]])
        plain = vis.plot_3d_bounds(data)
        highlighted = vis.plot_3d_bounds(data, highlight_bound_id=data["bounds"][0].id)
        self.assertGreater(len(highlighted.data), len(plain.data))

    def test_bound_markers_carry_their_ids_for_selection(self):
        _, data = build([[3, 7]])
        markers = next(t for t in vis.plot_2d_bounds(data).data
                       if t.mode == "markers" and t.legendgroup == "local_bounds")
        self.assertEqual(list(markers.customdata), [b.id for b in data["bounds"]])

    def test_wireframe_drops_meshes_and_keeps_edges(self):
        def edges(fig):
            return [(t.legendgroup, tuple(t.x), tuple(t.y), tuple(t.z))
                    for t in fig.data if t.type == "scatter3d" and t.mode == "lines"]

        _, data = build(PAPER2_POINTS[:2])
        solid = vis.plot_3d_bounds(data, show_search_boxes=True)
        wire = vis.plot_3d_bounds(data, show_search_boxes=True, wireframe=True)
        self.assertTrue(any(t.type == "mesh3d" for t in solid.data))
        self.assertFalse(any(t.type == "mesh3d" for t in wire.data))
        self.assertEqual(edges(wire), edges(solid))

    def test_label_mode_none_hides_labels_and_keeps_hover_text(self):
        _, data3 = build(PAPER2_POINTS[:2])
        _, data2 = build([[3, 7], [5, 4]])
        for plot, data in ((vis.plot_3d_bounds, data3), (vis.plot_pairwise_projections_2d, data3),
                           (vis.plot_2d_bounds, data2)):
            for mode in ("all", "none"):
                with self.subTest(plot=plot.__name__, mode=mode):
                    figure = plot(data, label_mode=mode)
                    annotations = (figure.layout.scene.annotations if plot == vis.plot_3d_bounds
                                   else figure.layout.annotations)
                    self.assertEqual(any(label.name for label in annotations), mode == "all")
                    markers = [t for t in figure.data if getattr(t, "mode", None) == "markers"]
                    self.assertTrue(markers)
                    for trace in markers:
                        self.assertTrue(all(trace.hovertext))
                        self.assertTrue(all(bool(text) == (mode == "all") for text in trace.text))

    def test_coincident_projections_share_one_marker(self):
        # u0 and u1 project onto the same spot in the (f1, f2) panel, as do z0 and z1.
        records = [vis.BoundRecord(i, name, coords, [], [None] * 3)
                   for i, (name, coords) in enumerate((("u0", [4, 4, 2]), ("u1", [4, 4, 7]), ("u2", [8, 3, 1])))]
        data = dict(dimensions=3, sense="MINIMIZE", lower_bound=[0] * 3, upper_bound=[10] * 3,
                    reference_point=[10] * 3, anti_reference=[0] * 3, bounds=records,
                    points=[{"id": "z0", "coordinates": [2, 2, 5]}, {"id": "z1", "coordinates": [2, 2, 7]}])
        fig = vis.plot_pairwise_projections_2d(data)
        bounds = next(t for t in fig.data if t.legendgroup == "local_bounds" and t.mode == "markers")
        self.assertEqual(list(bounds.customdata), [["u0", "u1"], "u2"])
        self.assertEqual(list(bounds.x), [4, 8])
        points = next(t for t in fig.data if t.legendgroup == "points")
        self.assertEqual(list(points.customdata), [["z0", "z1"]])
        self.assertEqual([b.coordinates for b in records], [[4, 4, 2], [4, 4, 7], [8, 3, 1]])  # input untouched

    def test_view_revision_changes_only_with_the_problem(self):
        # Plotly keeps the camera while uirevision is unchanged: stepping must not reset the view.
        state = web_app.AppState()
        initial, later = state.get_bounds_at_step(0), state.get_bounds_at_step(2)
        for plot in (vis.plot_3d_bounds, vis.plot_2d_bounds, vis.plot_pairwise_projections_2d):
            with self.subTest(plot=plot.__name__):
                revision = plot(initial).layout.uirevision
                self.assertEqual(plot(later, highlight_bound_id=later["bounds"][0].id).layout.uirevision, revision)
                self.assertNotEqual(plot(dict(initial, sense="MAXIMIZE")).layout.uirevision, revision)
                self.assertNotEqual(plot(dict(initial, lower_bound=[-1] * 3)).layout.uirevision, revision)

    def test_point_ids_are_escaped_in_figures_and_reports(self):
        payload = '<img src=x onerror="alert(1)"></script>'
        nbs = lb.NeighborhoodBoundSet([10, 10], [0, 0])
        point = lb.Point(payload, [3, 7])
        nbs.update(point)
        data = vis.extract_bounds_data(nbs, [10, 10], [0, 0], [point])
        self.assertNotIn(payload, vis.plot_2d_bounds(data).to_json())
        with tempfile.TemporaryDirectory() as directory, contextlib.redirect_stdout(io.StringIO()):
            output = Path(directory) / "report.html"
            cli.generate_html_report(data, str(output))
            html = output.read_text()
        self.assertNotIn(payload, html)
        self.assertIn("&lt;img", html)
        self.assertNotIn('src="https://cdn.plot.ly', html)  # reports work offline


class TestDashboardApi(unittest.TestCase):

    def setUp(self):
        web_app.state.load_preset("paper2")  # the app keeps one shared state; start each test from it
        self.client = web_app.app.test_client()

    def test_walkthrough_of_a_session(self):
        client = self.client
        self.assertEqual(client.get("/").status_code, 200)
        self.assertEqual(client.get("/static/plotly.js").status_code, 200)

        state = client.get("/api/state").get_json()
        self.assertEqual(state["bounds_count"], 7)
        self.assertEqual([(h["point_added"], h["destroyed_bounds"], h["new_bounds"]) for h in state["history"]],
                         PAPER2_HISTORY)
        figure = client.get("/api/figure_3d").get_json()

        rejected = client.post("/api/add_point", json={"id": "bad", "coordinates": [5, 5, 5]})
        self.assertEqual((rejected.status_code, rejected.get_json()["code"]), (400, "DOMINATED_POINT"))

        deleted = client.post("/api/delete_point", json={"index": 2}).get_json()
        self.assertEqual((len(deleted["points"]), deleted["bounds_count"]), (2, 5))  # Example 2.8 (2)

        # Same number of points as at the start, but a different third point:
        # a stale cache would serve the original table and figure again.
        replaced = client.post("/api/add_point", json={"id": "z4", "coordinates": [1.5, 5, 1.5]}).get_json()
        self.assertEqual(len(replaced["points"]), 3)
        self.assertIsNone(replaced["preset"])  # an edited preset is no longer that preset
        self.assertNotEqual(replaced["table"], state["table"])
        self.assertNotEqual(client.get("/api/figure_3d").get_json(), figure)

        self.assertLessEqual({"nodes", "edges"}, set(client.get("/api/graph?mode=combined").get_json()))
        self.assertEqual(client.get("/api/figure_pairwise").status_code, 200)

    def test_maximize_starts_from_the_lower_corner(self):
        configured = self.client.post("/api/configure", json={"dimensions": 3, "sense": "MAXIMIZE"}).get_json()
        self.assertEqual((configured["sense"], configured["bounds_count"]), ("MAXIMIZE", 1))
        self.assertEqual(list(configured["table"][0]["Coordinates"]), [0, 0, 0])
        added = self.client.post("/api/add_point", json={"id": "z1", "coordinates": [6, 10, 6]}).get_json()
        self.assertEqual(added["bounds_count"], 3)
        self.assertEqual(self.client.get("/api/figure_3d").status_code, 200)
        self.assertEqual(self.client.post("/api/configure", json={"dimensions": 2, "sense": "MAXIMIZE"})
                         .get_json()["dimensions"], 2)
        self.assertEqual(self.client.get("/api/figure_3d").status_code, 200)

    def test_every_preset_loads(self):
        for preset, (sense, coordinates) in PRESETS.items():
            with self.subTest(preset=preset):
                data = self.client.post("/api/load_preset", json={"preset": preset}).get_json()
                self.assertEqual((data["preset"], data["sense"], data["dimensions"]),
                                 (preset, sense, len(coordinates[0])))
                self.assertEqual([p["coordinates"] for p in data["points"]], coordinates)
        # Klamroth et al. (2015), Example 2 ends with five local upper bounds.
        loaded = self.client.post("/api/load_preset", json={"preset": "paper1_sa"}).get_json()
        self.assertEqual(loaded["bounds_count"], 5)

    def test_malformed_requests_leave_the_state_unchanged(self):
        client = self.client
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
            ("add_point", {"coordinates": [2.0, 2.0]}),
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
            with self.subTest(body=body):
                response = client.post("/api/configure", data=body, content_type="application/json")
                self.assertEqual(response.status_code, 400)
                self.assertIn("code", response.get_json())
        for query in ("/api/state?step=oops", "/api/graph?mode=unknown", "/api/graph?view=invalid",
                      "/api/figure_3d?dom_opacity=nan", "/api/figure_3d?labels=bad"):
            with self.subTest(query=query):
                self.assertEqual(client.get(query).status_code, 400)

    def test_probe_reports_membership_at_exact_boundaries(self):
        client = self.client
        for sense in ("MINIMIZE", "MAXIMIZE"):
            is_min = sense == "MINIMIZE"
            client.post("/api/configure", json={"dimensions": 2, "sense": sense})
            point = [3, 7] if is_min else [7, 3]
            client.post("/api/add_point", json={"id": "z", "coordinates": point})
            before = client.get("/api/state").get_json()
            outside = [10, 5] if is_min else [0, 5]
            anti = [0, 5] if is_min else [10, 5]
            for coords, in_domain, in_region in ((point, True, False), (outside, False, False), (anti, True, True)):
                with self.subTest(sense=sense, coords=coords):
                    result = client.post("/api/probe", json={"coordinates": coords}).get_json()
                    self.assertEqual((result["in_domain"], result["in_search_region"]), (in_domain, in_region))
                    if in_region:
                        self.assertIn(result["containing_bound"]["id"], {row["ID"] for row in before["table"]})
                    else:
                        self.assertIsNone(result["containing_bound"])
            # Before the point was inserted (step 0) it was still in the search region.
            self.assertTrue(client.post("/api/probe?step=0", json={"coordinates": point})
                            .get_json()["in_search_region"])
            for coords in ([1], [float("nan"), 1], [True, 1], "1,1"):
                self.assertEqual(client.post("/api/probe", json={"coordinates": coords}).status_code, 400)
            self.assertEqual(client.get("/api/state").get_json(), before)  # probing changes nothing

    def test_graph_offers_raw_and_contracted_views(self):
        self.client.post("/api/load_preset", json={"preset": "paper1_ngp"})
        raw = self.client.get("/api/graph?view=raw").get_json()
        clean = self.client.get("/api/graph").get_json()
        self.assertEqual(raw["counts"], clean["counts"])
        self.assertGreater(raw["counts"]["quasi"], 0)
        self.assertEqual(len(raw["nodes"]), raw["counts"]["total"])
        self.assertEqual(len(clean["nodes"]), raw["counts"]["nonredundant"])

    def test_step_data_is_not_shared_between_calls(self):
        state = web_app.AppState()
        first = state.get_bounds_at_step(3)
        table = vis.create_bounds_table(first)
        first["bounds"][0].coordinates[0] = -99
        first["points"][0]["coordinates"][0] = -99
        again = state.get_bounds_at_step(3)
        self.assertEqual(vis.create_bounds_table(again), table)
        self.assertEqual(again["points"][0]["coordinates"], PAPER2_POINTS[0])


def run_cli(*args):
    """Runs the CLI without writing a report; returns the report data and what was printed."""
    stdout = io.StringIO()
    with patch("sys.argv", ["cli", *args]), patch.object(cli, "generate_html_report") as report, \
            contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(io.StringIO()):
        cli.main()
    return report.call_args.args[0], stdout.getvalue()


class TestCli(unittest.TestCase):

    def test_cli_uses_the_shared_presets(self):
        for preset, (sense, coordinates) in PRESETS.items():
            with self.subTest(preset=preset):
                data, _ = run_cli("--preset", preset)
                self.assertEqual(data["sense"], sense)
                self.assertEqual([p["coordinates"] for p in data["points"]], coordinates)

    def test_cli_accepts_custom_points_and_interval(self):
        data, _ = run_cli("--points", "[[3,7]]", "--ref", "20,20", "--anti", "0,0")
        self.assertEqual((data["dimensions"], data["reference_point"]), (2, [20, 20]))
        self.assertEqual(run_cli("--points", "[]")[0]["points"], [])

    def test_step_by_step_reports_the_same_final_state(self):
        for args, count in (([], 7), (["--preset", "max_3d"], 7), (["--points", "[]"], 1)):
            with self.subTest(args=args):
                final, quiet = run_cli(*args)
                stepped, printed = run_cli(*args, "--step-by-step")
                self.assertEqual(len(final["bounds"]), count)
                self.assertEqual([b.to_dict() for b in stepped["bounds"]], [b.to_dict() for b in final["bounds"]])
                self.assertEqual(stepped["points"], final["points"])
                self.assertNotIn("--- Step", quiet)
                self.assertEqual("--- Step" in printed, bool(final["points"]))

    def test_malformed_input_exits_with_code_2(self):
        for args in (["--points", "[1]"],
                     ["--points", '[[10,5],["bad",2]]'],
                     ["--points", "[[3,7],[3,7]]"],
                     ["--points", "[[3,7]]", "--anti", "10,10"]):
            for mode in ([], ["--step-by-step"]):
                with self.subTest(args=args + mode), self.assertRaises(SystemExit) as error:
                    run_cli(*args, *mode)
                self.assertEqual(error.exception.code, 2)


if __name__ == "__main__":
    unittest.main()
