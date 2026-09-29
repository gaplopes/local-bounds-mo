#!/usr/bin/env python3
"""
Command-Line Script to Visualize Local Bounds Generation in 2D and 3D.
"""

import os
import sys
import argparse
import json
from typing import List

# Ensure library is accessible
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))

import local_bounds as lb
from . import engine as vis


def format_terminal_table(table_data: List[dict], dimensions: int) -> str:
    """Renders a clean ASCII table for the terminal."""
    if not table_data:
        return "No bounds to display."

    headers = ["ID", "Coordinates"]
    for j in range(1, dimensions + 1):
        headers.append(f"z^{j}(u)")
    for k in range(1, dimensions + 1):
        headers.append(f"ν_{k}")
    headers.append("Status")

    rows = []
    for row in table_data:
        r = [row["ID"], str(row["Coordinates"])]
        for j in range(1, dimensions + 1):
            r.append(str(row[f"z^{j}(u)"]))
        for k in range(1, dimensions + 1):
            r.append(str(row[f"ν_{k}"]))
        r.append(str(row["Status"]))
        rows.append(r)

    col_widths = [len(h) for h in headers]
    for r in rows:
        for i, val in enumerate(r):
            col_widths[i] = max(col_widths[i], len(val))

    def make_divider(char="-"):
        return "+" + "+".join(char * (w + 2) for w in col_widths) + "+"

    def format_row(r_vals):
        return "| " + " | ".join(val.ljust(col_widths[i]) for i, val in enumerate(r_vals)) + " |"

    lines = [
        make_divider("="),
        format_row(headers),
        make_divider("="),
    ]
    for r in rows:
        lines.append(format_row(r))
    lines.append(make_divider("-"))

    return "\n".join(lines)


def generate_html_report(bounds_data: dict, output_path: str):
    """Generates an all-in-one standalone HTML report with 3D/2D views, graph, and table."""
    dims = bounds_data["dimensions"]

    if dims == 2:
        fig_spatial = vis.plot_2d_bounds(bounds_data)
        spatial_title = "2D Search Region & Local Bounds"
    else:
        fig_spatial = vis.plot_3d_bounds(bounds_data, show_occupied_boxes=True, show_search_boxes=False)
        spatial_title = "3D Search Region & Local Bounds"

    fig_graph = vis.plot_neighbor_graph_plotly(bounds_data, mode="combined")
    table_rows = vis.create_bounds_table(bounds_data)

    spatial_html = fig_spatial.to_html(full_html=False, include_plotlyjs="cdn")
    graph_html = fig_graph.to_html(full_html=False, include_plotlyjs=False)

    extra_pairwise_html = ""
    if dims == 3:
        fig_pairwise = vis.plot_pairwise_projections_2d(bounds_data)
        extra_pairwise_html = f"""
        <div class="card">
            <h2>Pairwise Projections Matrix</h2>
            {fig_pairwise.to_html(full_html=False, include_plotlyjs=False)}
        </div>
        """

    table_headers = "".join([f"<th>{k}</th>" for k in table_rows[0].keys()])
    table_body = ""
    for r in table_rows:
        tds = "".join([f"<td>{v}</td>" for v in r.values()])
        table_body += f"<tr>{tds}</tr>"

    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>Local Bounds Generation Report ({dims}D)</title>
    <style>
        body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; background: #0f172a; color: #f8fafc; margin: 0; padding: 24px; }}
        h1, h2 {{ color: #e2e8f0; }}
        .header {{ border-bottom: 2px solid #334155; padding-bottom: 16px; margin-bottom: 24px; }}
        .grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 24px; margin-bottom: 24px; }}
        .card {{ background: #1e293b; border-radius: 12px; border: 1px solid #334155; padding: 20px; box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1); margin-bottom: 24px; }}
        table {{ width: 100%; border-collapse: collapse; font-size: 13px; font-family: monospace; }}
        th, td {{ padding: 10px 12px; text-align: left; border-bottom: 1px solid #334155; }}
        th {{ background: #334155; color: #cbd5e1; font-weight: 600; text-transform: uppercase; font-size: 11px; }}
        tr:hover {{ background: #283548; }}
        .badge {{ display: inline-block; padding: 4px 8px; border-radius: 6px; font-size: 12px; background: #4f46e5; color: white; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>Local Bounds Generation Report <span class="badge">{dims}D Space</span></h1>
        <p>Total Points: <b>{len(bounds_data['points'])}</b> | Total Local Bounds: <b>{len(bounds_data['bounds'])}</b></p>
    </div>

    <div class="grid">
        <div class="card">
            <h2>{spatial_title}</h2>
            {spatial_html}
        </div>
        <div class="card">
            <h2>Neighbor Graph G = (U(N), ν)</h2>
            {graph_html}
        </div>
    </div>

    {extra_pairwise_html}

    <div class="card">
        <h2>Local Bounds Table U(N)</h2>
        <table>
            <thead><tr>{table_headers}</tr></thead>
            <tbody>{table_body}</tbody>
        </table>
    </div>
</body>
</html>
"""
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html_content)
    print(f"[OK] Report successfully written to: {os.path.abspath(output_path)}")


def main():
    parser = argparse.ArgumentParser(description="Visualize Local Bounds in 2D and 3D")
    parser.add_argument("--dim", type=int, default=3, choices=[2, 3], help="Dimensionality (2 or 3)")
    parser.add_argument("--preset", type=str, default="paper2", choices=["paper2", "paper1_sa", "paper1_ngp", "2d"], help="Preset dataset")
    parser.add_argument("--points", type=str, default=None, help='Custom JSON points list e.g. "[[4,0,4],[3,3,1],[2,2,2]]"')
    parser.add_argument("--ref", type=str, default=None, help='Reference point M (nadir) e.g. "10,10,10"')
    parser.add_argument("--anti", type=str, default=None, help='Anti-reference point m (ideal) e.g. "0,0,0"')
    parser.add_argument("--output", type=str, default="local_bounds_report.html", help="HTML report output path")
    parser.add_argument("--step-by-step", action="store_true", help="Print step-by-step generation evolution")
    parser.add_argument("--interactive", action="store_true", help="Launch interactive web dashboard")
    parser.add_argument("--port", type=int, default=8050, help="Port for interactive web dashboard")
    args = parser.parse_args()

    if args.interactive:
        from . import app as web_app
        web_app.state.load_preset(args.preset)
        web_app.app.run(port=args.port, host="127.0.0.1")
        return

    dims = args.dim
    ref = [10.0] * dims
    anti = [0.0] * dims

    if args.ref:
        ref = [float(x.strip()) for x in args.ref.split(",")]
    if args.anti:
        anti = [float(x.strip()) for x in args.anti.split(",")]

    points_to_insert = []

    if args.points:
        raw_pts = json.loads(args.points)
        for i, coords in enumerate(raw_pts):
            points_to_insert.append(lb.Point(f"z{i+1}", coords))
        dims = len(raw_pts[0])
    elif args.preset == "paper2":
        dims = 3
        ref = [10.0, 10.0, 10.0]
        anti = [0.0, 0.0, 0.0]
        points_to_insert = [
            lb.Point("z1", [4.0, 0.0, 4.0]),
            lb.Point("z2", [3.0, 3.0, 1.0]),
            lb.Point("z3", [2.0, 2.0, 2.0])
        ]
    elif args.preset == "paper1_sa":
        dims = 3
        ref = [10.0, 10.0, 10.0]
        anti = [0.0, 0.0, 0.0]
        points_to_insert = [
            lb.Point("z1", [3.0, 7.0, 5.0]),
            lb.Point("z2", [5.0, 5.0, 4.0])
        ]
    elif args.preset == "paper1_ngp":
        dims = 3
        ref = [10.0, 10.0, 10.0]
        anti = [0.0, 0.0, 0.0]
        points_to_insert = [
            lb.Point("z1", [4.0, 3.0, 7.0]),
            lb.Point("z2", [4.0, 5.0, 4.0]),
            lb.Point("z3", [2.0, 5.0, 7.0])
        ]
    elif args.preset == "2d":
        dims = 2
        ref = [10.0, 10.0]
        anti = [0.0, 0.0]
        points_to_insert = [
            lb.Point("z1", [3.0, 7.0]),
            lb.Point("z2", [5.0, 4.0]),
            lb.Point("z3", [7.0, 2.0])
        ]

    print("===================================================================")
    print(f"  Local Bounds Generator ({dims}D Space)")
    print(f"  Reference Point M: {ref}")
    print(f"  Anti-Reference m: {anti}")
    print(f"  Points to insert: {len(points_to_insert)}")
    print("===================================================================\n")

    tracker = vis.LocalBoundsTracker(ref, anti)
    for p in points_to_insert:
        step = tracker.add_point(p)
        if args.step_by_step:
            print(f"--- Step {step.step}: Added Point {step.point_added.id} ({step.point_added.coordinates}) ---")
            print(f"  Bounds destroyed: {step.destroyed_bounds if step.destroyed_bounds else 'None'}")
            print(f"  Bounds created:   {step.new_bounds}")
            print(f"  Current |U(N)|:   {len(step.bounds_data['bounds'])}\n")

    final_data = tracker.steps[-1].bounds_data
    table = vis.create_bounds_table(final_data)

    print("Final Local Bounds Table U(N):")
    print(format_terminal_table(table, dims))
    print(f"\nTotal Local Bounds: {len(table)}\n")

    generate_html_report(final_data, args.output)


if __name__ == "__main__":
    main()
