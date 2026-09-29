#!/usr/bin/env python3
"""
Flask Web Application for Local Bounds Interactive Dashboard.
Uses templates/index.html for presentation.
"""

import os
import sys
import argparse
from typing import List, Dict, Any

from flask import Flask, jsonify, request, Response, render_template
import plotly.offline as po

# Ensure library is accessible
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))

import local_bounds as lb
from . import engine as vis

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
        self.points = []
        self.current_step = 0

    def load_preset(self, preset_name: str):
        self.reset()
        if preset_name == "paper2":
            # Example 2.8 from Dächert et al. (2017)
            self.dimensions = 3
            self.sense = lb.Objective.MINIMIZE
            self.lower_bound = [0.0, 0.0, 0.0]
            self.upper_bound = [10.0, 10.0, 10.0]
            self.points = [
                lb.Point("z1", [4.0, 0.0, 4.0]),
                lb.Point("z2", [3.0, 3.0, 1.0]),
                lb.Point("z3", [2.0, 2.0, 2.0])
            ]
        elif preset_name == "paper1_sa":
            # Example 2 from Klamroth et al. (2015) (Figure 2)
            self.dimensions = 3
            self.sense = lb.Objective.MINIMIZE
            self.lower_bound = [0.0, 0.0, 0.0]
            self.upper_bound = [10.0, 10.0, 10.0]
            self.points = [
                lb.Point("z1", [3.0, 5.0, 7.0]),
                lb.Point("z2", [6.0, 2.0, 4.0])
            ]
        elif preset_name == "paper1_ngp":
            # Example 3 from Klamroth et al. (2015) - Non-General Position (ties)
            self.dimensions = 3
            self.sense = lb.Objective.MINIMIZE
            self.lower_bound = [0.0, 0.0, 0.0]
            self.upper_bound = [10.0, 10.0, 10.0]
            self.points = [
                lb.Point("z1", [4.0, 3.0, 7.0]),
                lb.Point("z2", [4.0, 5.0, 4.0]),
                lb.Point("z3", [2.0, 5.0, 7.0])
            ]
        elif preset_name == "2d":
            # 2-Objective Example
            self.dimensions = 2
            self.sense = lb.Objective.MINIMIZE
            self.lower_bound = [0.0, 0.0]
            self.upper_bound = [10.0, 10.0]
            self.points = [
                lb.Point("z1", [3.0, 7.0]),
                lb.Point("z2", [5.0, 4.0]),
                lb.Point("z3", [7.0, 2.0])
            ]
        elif preset_name == "max_3d":
            # 3D Maximization Example
            self.dimensions = 3
            self.sense = lb.Objective.MAXIMIZE
            self.lower_bound = [0.0, 0.0, 0.0]
            self.upper_bound = [10.0, 10.0, 10.0]
            self.points = [
                lb.Point("z1", [6.0, 10.0, 6.0]),
                lb.Point("z2", [7.0, 7.0, 9.0]),
                lb.Point("z3", [8.0, 8.0, 8.0])
            ]
        else:
            self.load_preset("paper2")

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

    def get_bounds_at_step(self, step_idx: int) -> Dict[str, Any]:
        pts_subset = self.points[:step_idx]
        nbs = lb.NeighborhoodBoundSet(
            self.reference_point, self.anti_reference, sense=self.sense
        )
        for p in pts_subset:
            nbs.update(p)
        return vis.extract_bounds_data(
            nbs,
            self.reference_point,
            self.anti_reference,
            pts_subset,
            lower_bound=self.lower_bound,
            upper_bound=self.upper_bound,
            sense=self.sense
        )


state = AppState()


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/static/plotly.js")
def serve_plotly_js():
    return Response(get_cached_plotly_js(), mimetype="application/javascript")


@app.route("/api/state", methods=["GET"])
def get_state():
    step_idx = request.args.get("step", default=state.current_step, type=int)
    step_idx = max(0, min(len(state.points), step_idx))
    
    b_data = state.get_bounds_at_step(step_idx)
    table_rows = vis.create_bounds_table(b_data)
    history = state.compute_step_history()
    
    return jsonify({
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
        "table": table_rows
    })


@app.route("/api/configure", methods=["POST"])
def configure():
    data = request.json or {}
    dims = int(data.get("dimensions", state.dimensions))
    if dims not in (2, 3):
        return jsonify({"error": "Dimensions must be 2 or 3", "code": "INVALID_DIMENSIONS"}), 400

    sense_str = data.get("sense", "MINIMIZE").upper()
    sense = lb.Objective.MINIMIZE if sense_str == "MINIMIZE" else lb.Objective.MAXIMIZE

    # Allow setting by lower_bound/upper_bound directly, with fallback to legacy ref/anti
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

    try:
        lb_vals = [float(x) for x in raw_lb]
        ub_vals = [float(x) for x in raw_ub]
    except (ValueError, TypeError):
        return jsonify({"error": "Bounds must be numeric lists.", "code": "INVALID_BOUNDS"}), 400

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
    state.points = []
    state.current_step = 0
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
    p_id = data.get("id", "").strip() or f"z{len(state.points) + 1}"
    raw_coords = data.get("coordinates", [])

    try:
        coords = [float(x) for x in raw_coords]
    except (ValueError, TypeError):
        return jsonify({"error": "Coordinates must be a list of numbers.", "code": "NON_NUMERIC"}), 400

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

    state.points.append(lb.Point(p_id, coords))
    state.current_step = len(state.points)
    return get_state()


@app.route("/api/delete_point", methods=["POST"])
def delete_point():
    data = request.json or {}
    idx = int(data.get("index", -1))
    if 0 <= idx < len(state.points):
        state.points.pop(idx)
        state.current_step = min(state.current_step, len(state.points))
    return get_state()


@app.route("/api/graph", methods=["GET"])
def get_graph():
    step_idx = request.args.get("step", default=state.current_step, type=int)
    step_idx = max(0, min(len(state.points), step_idx))
    mode = request.args.get("mode", default="combined", type=str)
    
    b_data = state.get_bounds_at_step(step_idx)
    elements = vis.extract_network_graph_elements(b_data, mode=mode)
    return jsonify(elements)


@app.route("/api/figure_3d", methods=["GET"])
def get_figure_3d():
    step_idx = request.args.get("step", default=state.current_step, type=int)
    step_idx = max(0, min(len(state.points), step_idx))
    show_occupied = request.args.get("show_occupied", default="1") == "1"
    show_boxes = request.args.get("show_boxes", default="0") == "1"
    hl_bound = request.args.get("highlight", default=None, type=str)
    try:
        dom_opacity = max(0.0, min(1.0, float(request.args.get("dom_opacity", 1.0))))
    except (ValueError, TypeError):
        dom_opacity = 1.0

    b_data = state.get_bounds_at_step(step_idx)
    if state.dimensions == 2:
        fig = vis.plot_2d_bounds(b_data)
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
    step_idx = request.args.get("step", default=state.current_step, type=int)
    step_idx = max(0, min(len(state.points), step_idx))
    hl_bound = request.args.get("highlight", default=None, type=str)
    b_data = state.get_bounds_at_step(step_idx)
    fig = vis.plot_pairwise_projections_2d(b_data, highlight_bound_id=hl_bound)
    return Response(fig.to_json(), mimetype="application/json")


def main():
    parser = argparse.ArgumentParser(description="Interactive Local Bounds Visualizer")
    parser.add_argument("--port", type=int, default=8050, help="Port to run web dashboard (default: 8050)")
    parser.add_argument("--host", type=str, default="127.0.0.1", help="Host IP address (default: 127.0.0.1)")
    parser.add_argument("--preset", type=str, default="paper2", choices=["paper2", "paper1_sa", "paper1_ngp", "2d"], help="Initial preset")
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
