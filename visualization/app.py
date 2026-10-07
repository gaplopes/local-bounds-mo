#!/usr/bin/env python3
"""
Flask Web Application for Local Bounds Interactive Dashboard.
Uses templates/index.html for presentation.
"""

import os
import argparse
import math
from threading import RLock
from typing import List, Dict, Any

from flask import Flask, jsonify, request, Response, render_template, g
import plotly.offline as po

import local_bounds as lb
from . import engine as vis
from .presets import PRESETS

TEMPLATES_DIR = os.path.join(os.path.dirname(__file__), "templates")
STATIC_DIR = os.path.join(os.path.dirname(__file__), "static")

app = Flask(__name__, template_folder=TEMPLATES_DIR, static_folder=STATIC_DIR)

PLOTLY_JS = None

def get_cached_plotly_js():
    global PLOTLY_JS
    if PLOTLY_JS is None:
        PLOTLY_JS = po.get_plotlyjs()
    return PLOTLY_JS


class AppState:
    """Manages problem settings and points for 2D/3D instances."""
    def __init__(self):
        self.dimensions = 3
        self.sense = lb.Objective.MINIMIZE
        self.lower_bound = [0.0, 0.0, 0.0]
        self.upper_bound = [10.0, 10.0, 10.0]
        self.points: List[lb.Point] = []
        self.current_step = 0
        self.load_preset("paper2")

    @property
    def reference_point(self) -> List[float]:
        # For MINIMIZE: reference point is Nadir M (upper bound)
        # For MAXIMIZE: reference point is Nadir m (lower bound)
        return list(self.upper_bound if self.sense == lb.Objective.MINIMIZE else self.lower_bound)

    @property
    def anti_reference(self) -> List[float]:
        # For MINIMIZE: anti-reference is Ideal m (lower bound)
        # For MAXIMIZE: anti-reference is Ideal M (upper bound)
        return list(self.lower_bound if self.sense == lb.Objective.MINIMIZE else self.upper_bound)

    def reset(self):
        self.preset_name = None
        self.points = []
        self.current_step = 0

    def load_preset(self, preset_name: str):
        if not isinstance(preset_name, str) or preset_name not in PRESETS:
            raise vis.PointValidationError("Unknown preset.", "INVALID_PRESET")
        self.preset_name = preset_name
        sense, coordinates = PRESETS[preset_name]
        self.dimensions = len(coordinates[0])
        self.sense = lb._normalize_sense(sense)
        self.lower_bound = [0.0] * self.dimensions
        self.upper_bound = [10.0] * self.dimensions
        self.points = [lb.Point(f"z{i+1}", coords) for i, coords in enumerate(coordinates)]
        self.current_step = len(self.points)

    def compute_step_history(self) -> List[Dict[str, Any]]:
        tracker = vis.LocalBoundsTracker(
            self.reference_point,
            self.anti_reference,
            sense=self.sense,
            lower_bound=self.lower_bound,
            upper_bound=self.upper_bound
        )
        for p in self.points:
            tracker.add_point(p)

        history = []
        for step in tracker.steps:
            history.append({
                "step": step.step,
                "point_added": step.point_added.id if step.point_added else "Initial State",
                "destroyed_bounds": step.destroyed_bounds,
                "new_bounds": step.new_bounds,
                "total_bounds": len(step.bounds_data["bounds"])
            })
        return history

    def get_bound_set_at_step(self, step_idx: int):
        pts_subset = self.points[:step_idx]
        nbs = lb.NeighborhoodBoundSet(
            self.reference_point, self.anti_reference, sense=self.sense
        )
        for p in pts_subset:
            nbs.update(p)
        return nbs

    def get_bounds_at_step(self, step_idx: int, include_quasi=False) -> Dict[str, Any]:
        nbs = self.get_bound_set_at_step(step_idx)
        return vis.extract_bounds_data(
            nbs,
            self.reference_point,
            self.anti_reference,
            self.points[:step_idx],
            lower_bound=self.lower_bound,
            upper_bound=self.upper_bound,
            sense=self.sense,
            include_quasi=include_quasi
        )


state = AppState()
# ponytail: serialize the shared research session; use per-session snapshots if throughput matters.
state_lock = RLock()


@app.before_request
def lock_state():
    if request.path.startswith("/api/"):
        state_lock.acquire()
        g.state_locked = True


@app.teardown_request
def unlock_state(error):
    if getattr(g, "state_locked", False):
        g.state_locked = False
        state_lock.release()


@app.errorhandler(vis.PointValidationError)
def validation_error(error):
    return jsonify({"error": error.message, "code": error.code}), 400


@app.before_request
def validate_json_body():
    if request.method == "POST" and not isinstance(request.get_json(silent=True), dict):
        raise vis.PointValidationError("Request must contain a JSON object.", "INVALID_JSON")


def numeric_list(value, code):
    if not isinstance(value, list) or any(isinstance(x, bool) for x in value):
        raise vis.PointValidationError("Expected a list of finite numbers.", code)
    try:
        values = [float(x) for x in value]
    except (TypeError, ValueError, OverflowError):
        raise vis.PointValidationError("Expected a list of finite numbers.", code)
    if not all(math.isfinite(x) for x in values):
        raise vis.PointValidationError("Numbers must be finite.", code)
    return values


def selected_step():
    raw = request.args.get("step", str(state.current_step))
    try:
        step = int(raw)
    except (TypeError, ValueError):
        raise vis.PointValidationError("Step must be an integer.", "INVALID_STEP")
    return max(0, min(len(state.points), step))


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/static/plotly.js")
def serve_plotly_js():
    return Response(get_cached_plotly_js(), mimetype="application/javascript")


@app.route("/api/state", methods=["GET"])
def get_state():
    step_idx = selected_step()
    
    b_data = state.get_bounds_at_step(step_idx)
    table_rows = vis.create_bounds_table(b_data)
    history = state.compute_step_history()
    
    return jsonify({
        "preset": state.preset_name,
        "dimensions": state.dimensions,
        "sense": "MINIMIZE" if state.sense == lb.Objective.MINIMIZE else "MAXIMIZE",
        "lower_bound": state.lower_bound,
        "upper_bound": state.upper_bound,
        "reference_point": state.reference_point,
        "anti_reference": state.anti_reference,
        "points": [{"id": p.id, "coordinates": list(p.coordinates)} for p in state.points],
        "current_step": step_idx,
        "total_steps": len(state.points),
        "history": history,
        "bounds_count": len(b_data["bounds"]),
        "bounds_counts": b_data["counts"],
        "table": table_rows
    })


@app.route("/api/configure", methods=["POST"])
def configure():
    data = request.json or {}
    dims = data.get("dimensions", state.dimensions)
    if type(dims) is not int or dims not in (2, 3):
        return jsonify({"error": "Dimensions must be 2 or 3", "code": "INVALID_DIMENSIONS"}), 400

    sense_str = data.get("sense", "MINIMIZE")
    if not isinstance(sense_str, str) or sense_str.upper() not in ("MINIMIZE", "MAXIMIZE"):
        raise vis.PointValidationError("Sense must be MINIMIZE or MAXIMIZE.", "INVALID_SENSE")
    sense = lb._normalize_sense(sense_str.upper())

    # Allow setting by lower_bound/upper_bound directly, with fallback to legacy ref/anti
    if ("lower_bound" in data) != ("upper_bound" in data):
        raise vis.PointValidationError("Supply both lower_bound and upper_bound.", "INVALID_BOUNDS")
    if "lower_bound" in data and "upper_bound" in data:
        raw_lb = data.get("lower_bound")
        raw_ub = data.get("upper_bound")
    else:
        # Legacy fallback
        if sense == lb.Objective.MINIMIZE:
            raw_lb = data.get("anti_reference", [0.0] * dims)
            raw_ub = data.get("reference_point", [10.0] * dims)
        else:
            raw_lb = data.get("reference_point", [0.0] * dims)
            raw_ub = data.get("anti_reference", [10.0] * dims)

    lb_vals = numeric_list(raw_lb, "INVALID_BOUNDS")
    ub_vals = numeric_list(raw_ub, "INVALID_BOUNDS")

    if len(lb_vals) != dims or len(ub_vals) != dims:
        return jsonify({"error": f"Bounds must match dimension {dims}.", "code": "DIMENSION_MISMATCH"}), 400

    for i in range(dims):
        if lb_vals[i] >= ub_vals[i]:
            return jsonify({
                "error": f"Lower bound ({lb_vals[i]}) must be strictly less than upper bound ({ub_vals[i]}) at dimension {i+1}.",
                "code": "INVALID_INTERVAL"
            }), 400

    state.dimensions = dims
    state.sense = sense
    state.lower_bound = lb_vals
    state.upper_bound = ub_vals
    state.reset()
    return get_state()


@app.route("/api/load_preset", methods=["POST"])
def load_preset():
    data = request.json or {}
    preset = data.get("preset", "paper2")
    state.load_preset(preset)
    return get_state()


@app.route("/api/add_point", methods=["POST"])
def add_point():
    data = request.json or {}
    p_id = data.get("id", "")
    if not isinstance(p_id, str):
        raise vis.PointValidationError("Point ID must be a string.", "INVALID_ID")
    p_id = p_id.strip()
    if not p_id:
        next_id = len(state.points) + 1
        while any(p.id == f"z{next_id}" for p in state.points):
            next_id += 1
        p_id = f"z{next_id}"
    coords = numeric_list(data.get("coordinates", []), "NON_NUMERIC")

    try:
        vis.validate_point(
            coords=coords,
            point_id=p_id,
            existing_points=state.points,
            lower_bound=state.lower_bound,
            upper_bound=state.upper_bound,
            sense=state.sense
        )
    except vis.PointValidationError as e:
        return jsonify({"error": e.message, "code": e.code}), 400

    state.preset_name = None
    state.points.append(lb.Point(p_id, coords))
    state.current_step = len(state.points)
    return get_state()


@app.route("/api/delete_point", methods=["POST"])
def delete_point():
    data = request.json or {}
    idx = data.get("index")
    if type(idx) is not int or not 0 <= idx < len(state.points):
        raise vis.PointValidationError("Index must identify an existing point.", "INVALID_INDEX")
    state.preset_name = None
    state.points.pop(idx)
    state.current_step = min(state.current_step, len(state.points))
    return get_state()


@app.route("/api/graph", methods=["GET"])
def get_graph():
    step_idx = selected_step()
    mode = request.args.get("mode", default="combined", type=str)
    
    if mode not in ("combined", "undirected"):
        raise vis.PointValidationError("Unknown graph mode.", "INVALID_MODE")
    view = request.args.get("view", "contracted")
    if view not in ("contracted", "raw"):
        raise vis.PointValidationError("Unknown graph view.", "INVALID_VIEW")
    b_data = state.get_bounds_at_step(step_idx, include_quasi=view == "raw")
    elements = vis.extract_network_graph_elements(b_data, mode=mode)
    elements.update(view=view, counts=b_data["counts"])
    return jsonify(elements)


@app.route("/api/probe", methods=["POST"])
def probe():
    coords = numeric_list(request.json.get("coordinates", []), "NON_NUMERIC")
    if len(coords) != state.dimensions:
        raise vis.PointValidationError("Coordinates must match the problem dimension.", "DIMENSION_MISMATCH")
    nbs = state.get_bound_set_at_step(selected_step())
    is_min = state.sense == lb.Objective.MINIMIZE
    in_domain = all(lo <= z < hi if is_min else lo < z <= hi
                    for lo, z, hi in zip(state.lower_bound, coords, state.upper_bound))
    in_search = nbs.is_in_search_region(coords)
    containing = None
    if in_search:
        containing = next(b for b in nbs.nonredundant_bounds()
                          if all(z < u if is_min else z > u for z, u in zip(coords, b.coordinates)))
    inequalities = [f"{lo:g} ≤ f{j+1} < {hi:g}" if is_min else f"{lo:g} < f{j+1} ≤ {hi:g}"
                    for j, (lo, hi) in enumerate(zip(state.lower_bound, state.upper_bound))]
    zone = [] if containing is None else [
        f"{a:g} ≤ f{j+1} < {u:g}" if is_min else f"{u:g} < f{j+1} ≤ {a:g}"
        for j, (a, u) in enumerate(zip(state.anti_reference, containing.coordinates))]
    return jsonify(in_domain=in_domain, in_search_region=in_search,
                   containing_bound=None if containing is None else {"id": containing.id, "coordinates": list(containing.coordinates)},
                   domain=inequalities, search_zone=zone)


@app.route("/api/figure_3d", methods=["GET"])
def get_figure_3d():
    step_idx = selected_step()
    show_occupied = request.args.get("show_occupied", default="1") == "1"
    show_boxes = request.args.get("show_boxes", default="0") == "1"
    hl_bound = request.args.get("highlight", default=None, type=str)
    try:
        dom_opacity = float(request.args.get("dom_opacity", 1.0))
        if not math.isfinite(dom_opacity):
            raise ValueError()
        dom_opacity = max(0.0, min(1.0, dom_opacity))
    except (ValueError, TypeError):
        raise vis.PointValidationError("Opacity must be finite.", "INVALID_OPACITY")

    b_data = state.get_bounds_at_step(step_idx)
    if state.dimensions == 2:
        fig = vis.plot_2d_bounds(b_data, highlight_bound_id=hl_bound)
    else:
        fig = vis.plot_3d_bounds(
            b_data,
            show_occupied_boxes=show_occupied,
            show_search_boxes=show_boxes,
            highlight_bound_id=hl_bound,
            dom_opacity=dom_opacity
        )

    return Response(fig.to_json(), mimetype="application/json")


@app.route("/api/figure_pairwise", methods=["GET"])
def get_figure_pairwise():
    step_idx = selected_step()
    hl_bound = request.args.get("highlight", default=None, type=str)
    b_data = state.get_bounds_at_step(step_idx)
    fig = vis.plot_pairwise_projections_2d(b_data, highlight_bound_id=hl_bound)
    return Response(fig.to_json(), mimetype="application/json")


def main():
    parser = argparse.ArgumentParser(description="Interactive Local Bounds Visualizer")
    parser.add_argument("--port", type=int, default=8050, help="Port to run web dashboard (default: 8050)")
    parser.add_argument("--host", type=str, default="127.0.0.1", help="Host IP address (default: 127.0.0.1)")
    parser.add_argument("--preset", type=str, default="paper2", choices=list(PRESETS), help="Initial preset")
    args = parser.parse_args()

    state.load_preset(args.preset)
    print(f"\n=======================================================")
    print(f"  Local Bounds Multiobjective Visualizer running at:")
    print(f"  --> http://{args.host}:{args.port}")
    print(f"  Preset: {args.preset} ({state.dimensions}D)")
    print(f"=======================================================\n")
    app.run(host=args.host, port=args.port, debug=False)


if __name__ == "__main__":
    main()
