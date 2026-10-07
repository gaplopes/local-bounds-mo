"""
Core visualization engine for Local Bounds in 2D and 3D.
Based on:
- Klamroth, Lacour, Vanderpooten (EJOR 2015)
- Dächert, Klamroth, Lacour, Vanderpooten (EJOR 2017)
"""

from __future__ import annotations
import math
from html import escape
from typing import List, Dict, Any, Optional, Tuple, Union

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False

try:
    import networkx as nx
    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False

import local_bounds as lb


class PointValidationError(ValueError):
    """Exception raised when an invalid point is submitted."""
    def __init__(self, message: str, code: str):
        super().__init__(message)
        self.message = message
        self.code = code


def validate_point(
    coords: List[float],
    point_id: str,
    existing_points: List[lb.Point],
    lower_bound: List[float],
    upper_bound: List[float],
    sense: Union[lb.Objective, str] = lb.Objective.MINIMIZE
) -> None:
    """
    Validates a point prior to insertion into the nondominated set N.

    Checks performed:
    1. Dimension match with expected problem dimension p.
    2. Coordinates are finite numeric values (no NaN, Inf).
    3. Coordinates lie in the anti-inclusive, reference-exclusive interval.
    4. Unique point identifier (no duplicate ID).
    5. Unique coordinates (no duplicate in N).
    6. Pareto dominance checks:
       - Reject if dominated by an existing point in N.
       - Reject if it strictly dominates an existing point in N (violates N stability).
    """
    sense = lb._normalize_sense(sense)
    is_min = sense == lb.Objective.MINIMIZE
    if (len(lower_bound) != len(upper_bound) or not lower_bound or
            any(not math.isfinite(lo) or not math.isfinite(hi) or lo >= hi
                for lo, hi in zip(lower_bound, upper_bound))):
        raise PointValidationError("Bounds must define a finite, nonempty interval.", "INVALID_INTERVAL")
    if not isinstance(point_id, str):
        raise PointValidationError("Point ID must be a string.", "INVALID_ID")
    expected_dim = len(lower_bound)
    if len(coords) != expected_dim:
        raise PointValidationError(
            f"Dimension mismatch: expected {expected_dim} coordinates, received {len(coords)}.",
            "DIMENSION_MISMATCH"
        )

    for i, c in enumerate(coords):
        if isinstance(c, bool) or not isinstance(c, (int, float)) or not math.isfinite(c):
            raise PointValidationError(
                f"Coordinate f_{i+1} = '{c}' is invalid: must be a finite real number.",
                "NON_FINITE"
            )
        if not (lower_bound[i] <= c < upper_bound[i] if is_min else lower_bound[i] < c <= upper_bound[i]):
            raise PointValidationError(
                f"Coordinate f_{i+1} = {c} is out of bounds; must be within the anti-inclusive, reference-exclusive interval ({lower_bound[i]}, {upper_bound[i]}).",
                "OUT_OF_BOUNDS"
            )

    if point_id:
        for p in existing_points:
            if p.id == point_id:
                raise PointValidationError(
                    f"Point ID '{point_id}' is already used by an existing point in N.",
                    "DUPLICATE_ID"
                )

    for p in existing_points:
        p_coords = list(p.coordinates)
        if coords == p_coords:
            raise PointValidationError(
                f"Duplicate coordinates ({', '.join(f'{c:.2f}' for c in coords)}): matches existing point '{p.id}'.",
                "DUPLICATE_COORDINATES"
            )

    is_min = (sense == lb.Objective.MINIMIZE or sense == "MINIMIZE")

    # Pareto dominance checks
    for p in existing_points:
        p_coords = list(p.coordinates)
        if is_min:
            # For MINIMIZE: p dominates coords if p_coords <= coords and p_coords != coords
            p_weakly_dominates_new = all(pc <= c for pc, c in zip(p_coords, coords))
            new_weakly_dominates_p = all(c <= pc for c, pc in zip(coords, p_coords))
        else:
            # For MAXIMIZE: p dominates coords if p_coords >= coords and p_coords != coords
            p_weakly_dominates_new = all(pc >= c for pc, c in zip(p_coords, coords))
            new_weakly_dominates_p = all(c >= pc for c, pc in zip(coords, p_coords))

        if p_weakly_dominates_new:
            raise PointValidationError(
                f"Point is dominated by existing point '{p.id}' ({', '.join(f'{x:.2f}' for x in p_coords)}). "
                f"Only nondominated points can be added to N.",
                "DOMINATED_POINT"
            )
        if new_weakly_dominates_p:
            raise PointValidationError(
                f"Point dominates existing point '{p.id}' ({', '.join(f'{x:.2f}' for x in p_coords)}). "
                f"Nondominated set N must remain mutually nondominating.",
                "DOMINATES_EXISTING"
            )


def format_point_name(name: str) -> str:
    """Formats dummy point names like z_hat1 -> z^1 and z_hat3 -> z^3."""
    if not name:
        return name
    if name.startswith("z_hat"):
        dim_num = name.replace("z_hat", "").strip()
        return f"z^{dim_num}"
    return name


class BoundRecord:
    """Represents a formatted record for a local bound."""
    def __init__(
        self,
        index: int,
        bound_id: str,
        coordinates: List[float],
        defining_points: List[Dict[str, Any]],
        neighbors: List[Optional[int]],
        is_extreme: bool = False,
        extreme_dim: Optional[int] = None,
        is_quasi: bool = False
    ):
        self.index = index
        self.id = bound_id
        self.coordinates = coordinates
        self.defining_points = defining_points
        self.neighbors = neighbors
        self.is_extreme = is_extreme
        self.extreme_dim = extreme_dim
        self.is_quasi = is_quasi

    def to_dict(self) -> Dict[str, Any]:
        return {
            "index": self.index,
            "id": self.id,
            "coordinates": self.coordinates,
            "defining_points": self.defining_points,
            "neighbors": self.neighbors,
            "is_extreme": self.is_extreme,
            "extreme_dim": self.extreme_dim,
            "is_quasi": self.is_quasi
        }


def extract_bounds_data(
    bound_set: Union[lb.NeighborhoodBoundSet, lb.BoundSet],
    reference_point: List[float],
    anti_reference: List[float],
    points: Optional[List[lb.Point]] = None,
    lower_bound: Optional[List[float]] = None,
    upper_bound: Optional[List[float]] = None,
    sense: Union[lb.Objective, str] = lb.Objective.MINIMIZE,
    include_quasi: bool = False
) -> Dict[str, Any]:
    """
    Extracts bounds, defining points, and neighbor relationships for 2D and 3D sets.
    """
    dims = bound_set.dimensions()
    if dims not in (2, 3):
        raise ValueError(f"This visualizer supports 2D and 3D spaces (received {dims}D).")

    # Access C++ adjacency graph if available
    graph = None
    if hasattr(bound_set, "get_adjacency_graph"):
        graph = bound_set.get_adjacency_graph(include_quasi=include_quasi)

    if graph is not None:
        raw_bounds = graph.nodes
        adj_list = graph.adjacency_list
        k_neighbors = graph.k_neighbors
    else:
        if hasattr(bound_set, "nonredundant_bounds"):
            raw_bounds = bound_set.nonredundant_bounds()
        elif hasattr(bound_set, "bounds"):
            raw_bounds = bound_set.bounds
        else:
            raw_bounds = []
        adj_list = [[] for _ in raw_bounds]
        k_neighbors = [[None] * dims for _ in raw_bounds]

    bound_records: List[BoundRecord] = []
    NPOS = 18446744073709551615

    for idx, b in enumerate(raw_bounds):
        # Extract defining points and format dummy names (z_hat3 -> z^3)
        def_pts_info = []
        if hasattr(b, "defining_points") and b.defining_points:
            for p in b.defining_points:
                raw_id = p.id if p.id else "dummy"
                formatted_id = format_point_name(raw_id)
                def_pts_info.append({
                    "id": formatted_id,
                    "raw_id": raw_id,
                    "coordinates": list(p.coordinates)
                })
        else:
            def_pts_info = [{"id": f"z^{j+1}", "raw_id": f"z_{j+1}", "coordinates": []} for j in range(dims)]

        # Extract neighbors
        node_neighbors: List[Optional[int]] = []
        is_extreme = False
        extreme_dim = None

        if idx < len(k_neighbors):
            for k in range(dims):
                nb_val = k_neighbors[idx][k]
                if nb_val == NPOS or nb_val is None or nb_val >= len(raw_bounds):
                    node_neighbors.append(None)
                    is_extreme = True
                    extreme_dim = k
                else:
                    node_neighbors.append(int(nb_val))
        else:
            node_neighbors = [None] * dims

        record = BoundRecord(
            index=idx,
            bound_id=b.id if b.id else f"u_{idx}",
            coordinates=list(b.coordinates),
            defining_points=def_pts_info,
            neighbors=node_neighbors,
            is_extreme=is_extreme,
            extreme_dim=extreme_dim,
            is_quasi=bool(graph.quasi[idx]) if graph is not None else False
        )
        bound_records.append(record)

    edges_directed = [
        (u_idx, v_idx, k)
        for u_idx, record in enumerate(bound_records)
        for k, v_idx in enumerate(record.neighbors)
        if v_idx is not None and v_idx != u_idx and v_idx < len(bound_records)
    ]
    adjacency = ([(u, v) for u, neighbors in enumerate(adj_list) for v in neighbors]
                 if graph is not None else [(u, v) for u, v, _ in edges_directed])
    edges_undirected = sorted({(min(u, v), max(u, v)) for u, v in adjacency
                               if u != v and 0 <= v < len(bound_records)})

    is_min = (sense == lb.Objective.MINIMIZE or sense == "MINIMIZE")
    if lower_bound is None or upper_bound is None:
        if is_min:
            computed_lb = list(anti_reference)
            computed_ub = list(reference_point)
        else:
            computed_lb = list(reference_point)
            computed_ub = list(anti_reference)
    else:
        computed_lb = list(lower_bound)
        computed_ub = list(upper_bound)

    total = bound_set.size()
    nonredundant = sum(not b.is_quasi for b in bound_records)
    return {
        "dimensions": dims,
        "sense": "MINIMIZE" if is_min else "MAXIMIZE",
        "lower_bound": computed_lb,
        "upper_bound": computed_ub,
        "reference_point": reference_point,
        "anti_reference": anti_reference,
        "bounds": bound_records,
        "edges_undirected": edges_undirected,
        "edges_directed": edges_directed,
        "graph_view": "raw" if include_quasi else "contracted",
        "counts": {"total": total, "nonredundant": nonredundant, "quasi": total - nonredundant},
        "points": [{"id": p.id, "coordinates": list(p.coordinates)} for p in (points or [])]
    }


def create_bounds_table(bounds_data: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Creates a clean tabular representation of local bounds.
    Uses z^1, z^2, z^3 notation for defining points.
    """
    rows = []
    dims = bounds_data["dimensions"]
    for b in bounds_data["bounds"]:
        row = {
            "ID": b.id,
            "Coordinates": tuple(round(c, 4) for c in b.coordinates),
        }
        for j in range(dims):
            row[f"u_{j+1}"] = round(b.coordinates[j], 4)
            if j < len(b.defining_points):
                dp = b.defining_points[j]
                row[f"z^{j+1}(u)"] = dp["id"]
            else:
                row[f"z^{j+1}(u)"] = "-"

        for k in range(dims):
            nb_idx = b.neighbors[k]
            if nb_idx is not None and nb_idx < len(bounds_data["bounds"]):
                row[f"ν_{k+1}"] = bounds_data["bounds"][nb_idx].id
            else:
                row[f"ν_{k+1}"] = "∅ (Extreme)"

        row["Status"] = f"Extreme (dim {b.extreme_dim + 1})" if b.is_extreme else "Internal"
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# Plotly 3D Visualization
# ---------------------------------------------------------------------------

# Curated color palette for distinct Pareto dominance cones D(z^i) (100% solid opacity)
DOMINATED_CONE_PALETTE = [
    # (fill_color, wire_color)
    ("rgb(99, 102, 241)", "rgb(224, 231, 255)"),   # 1: Indigo (z1)
    ("rgb(168, 85, 247)", "rgb(243, 232, 255)"),   # 2: Purple (z2)
    ("rgb(16, 185, 129)", "rgb(209, 250, 229)"),   # 3: Emerald (z3)
    ("rgb(245, 158, 11)", "rgb(254, 243, 199)"),   # 4: Amber (z4)
    ("rgb(14, 165, 233)", "rgb(224, 242, 254)"),   # 5: Sky (z5)
    ("rgb(244, 63, 94)",  "rgb(255, 228, 230)"),   # 6: Rose (z6)
    ("rgb(20, 184, 166)", "rgb(204, 251, 241)"),   # 7: Teal (z7)
    ("rgb(217, 70, 239)", "rgb(250, 232, 255)"),   # 8: Fuchsia (z8)
]


def _make_box_mesh(
    p_min: List[float],
    p_max: List[float],
    color: str = "rgb(100, 116, 139)",
    wire_color: str = "rgb(203, 213, 225)",
    wire_width: float = 1.8,
    opacity: float = 1.0,
    name: str = "Zone",
    flatshading: bool = True
) -> Tuple[go.Mesh3d, go.Scatter3d]:
    """Generates 3D mesh faces and wireframe edges for a bounding box."""
    x0, y0, z0 = p_min[:3]
    x1, y1, z1 = p_max[:3]

    vx = [x0, x1, x1, x0, x0, x1, x1, x0]
    vy = [y0, y0, y1, y1, y0, y0, y1, y1]
    vz = [z0, z0, z0, z0, z1, z1, z1, z1]

    # 12 triangles (2 per cube face) with outward-facing normals
    i = [0, 0, 4, 4, 0, 0, 3, 3, 0, 0, 1, 1]
    j = [1, 2, 5, 6, 1, 5, 2, 6, 3, 7, 2, 6]
    k = [2, 3, 6, 7, 5, 4, 6, 7, 7, 4, 6, 5]

    mesh = go.Mesh3d(
        x=vx, y=vy, z=vz,
        i=i, j=j, k=k,
        color=color,
        opacity=opacity,
        flatshading=flatshading,
        lighting=dict(ambient=0.75, diffuse=0.8, specular=0.2, roughness=0.5) if flatshading else None,
        hoverinfo="skip",
        showlegend=False,
        name=name
    )

    edge_x = [
        x0, x1, x1, x0, x0, None,
        x0, x1, x1, x0, x0, None,
        x0, x0, None, x1, x1, None, x1, x1, None, x0, x0, None
    ]
    edge_y = [
        y0, y0, y1, y1, y0, None,
        y0, y0, y1, y1, y0, None,
        y0, y0, None, y0, y0, None, y1, y1, None, y1, y1, None
    ]
    edge_z = [
        z0, z0, z0, z0, z0, None,
        z1, z1, z1, z1, z1, None,
        z0, z1, None, z0, z1, None, z0, z1, None, z0, z1, None
    ]

    wireframe = go.Scatter3d(
        x=edge_x, y=edge_y, z=edge_z,
        mode="lines",
        line=dict(color=wire_color, width=wire_width),
        hoverinfo="skip",
        showlegend=False,
        name=f"{name} Edges"
    )

    return mesh, wireframe


def plot_3d_bounds(
    bounds_data: Dict[str, Any],
    show_occupied_boxes: bool = True,
    show_search_boxes: bool = False,
    show_defining_rays: bool = False,
    highlight_bound_id: Optional[str] = None,
    dom_opacity: float = 1.0
) -> go.Figure:
    """
    Creates an interactive 3D scene in Plotly matching Klamroth et al. (2015) Figure 2.
    - Occupied / Dominated cones D(N) in translucent shaded gray.
    - Local bounds U(N) sit on the reflex corners of the occupied staircase facing Ideal m.
    - Perspective looking from Ideal m towards Nadir M (for minimization) so search zones face the viewer.
    - Clicking a bound in the table highlights its specific search zone S(u).
    """
    if not HAS_PLOTLY:
        raise ImportError("plotly is required for 3D visualization.")

    fig = go.Figure()
    lb_pt = bounds_data.get("lower_bound", bounds_data["anti_reference"])
    ub_pt = bounds_data.get("upper_bound", bounds_data["reference_point"])
    is_max = (bounds_data.get("sense") == "MAXIMIZE")
    pts = bounds_data.get("points", [])
    b_records = bounds_data.get("bounds", [])

    # 1. Bounding overall search space wireframe [LB, UB]
    _, space_wire = _make_box_mesh(
        lb_pt, ub_pt,
        wire_color="rgba(148, 163, 184, 0.45)",
        wire_width=1.5,
        name="Search Space [m, M]"
    )
    space_wire.line.dash = "dot"
    fig.add_trace(space_wire)

    # 2. Occupied / Dominated Zones D(N) (Pareto Dominance Cones from Figure 2)
    # For MINIMIZE: D(z) = [z, M] (dominated region above point z)
    # For MAXIMIZE: D(z) = [m, z] (dominated region below point z)
    if show_occupied_boxes and pts:
        for idx, p in enumerate(pts):
            if not is_max:
                p_min = p["coordinates"][:3]
                p_max = ub_pt[:3]
            else:
                p_min = lb_pt[:3]
                p_max = p["coordinates"][:3]

            c_fill, c_wire = DOMINATED_CONE_PALETTE[idx % len(DOMINATED_CONE_PALETTE)]

            d_mesh, d_wire = _make_box_mesh(
                p_min, p_max,
                color=c_fill,
                wire_color=c_wire,
                wire_width=2.0,
                opacity=dom_opacity,
                name=f"Dominated D({escape(p['id'])})",
                flatshading=True
            )
            d_mesh.legendgroup = "dominated_zones"
            d_wire.legendgroup = "dominated_zones"

            if idx == 0:
                d_mesh.showlegend = True
                d_mesh.name = "Dominated Zones D(N)"
            else:
                d_mesh.showlegend = False

            fig.add_trace(d_mesh)
            fig.add_trace(d_wire)

    # 3. All Search Zones S(u) (Optional)
    # For MINIMIZE: S(u) = [LB, u]
    # For MAXIMIZE: S(u) = [u, UB]
    if show_search_boxes:
        for idx, b in enumerate(b_records):
            if not is_max:
                s_min = lb_pt[:3]
                s_max = b.coordinates[:3]
            else:
                s_min = b.coordinates[:3]
                s_max = ub_pt[:3]

            s_mesh, s_wire = _make_box_mesh(
                s_min, s_max,
                color="rgba(56, 189, 248, 0.10)",
                wire_color="rgba(56, 189, 248, 0.35)",
                wire_width=1.0,
                opacity=0.10,
                name=f"Zone {b.id}",
                flatshading=False
            )
            s_mesh.legendgroup = "search_zones"
            s_wire.legendgroup = "search_zones"
            if idx == 0:
                s_mesh.showlegend = True
                s_mesh.name = "Search Zones S(u)"
            fig.add_trace(s_mesh)
            fig.add_trace(s_wire)

    # 4. Highlighted Specific Search Zone S(u*)
    # When user selects a local bound in the table, highlight its exact search envelope
    if highlight_bound_id:
        target_b = next((b for b in b_records if b.id == highlight_bound_id), None)
        if target_b:
            if not is_max:
                h_min = lb_pt[:3]
                h_max = target_b.coordinates[:3]
            else:
                h_min = target_b.coordinates[:3]
                h_max = ub_pt[:3]

            h_mesh, h_wire = _make_box_mesh(
                h_min, h_max,
                color="rgba(249, 115, 22, 0.30)",
                wire_color="rgba(249, 115, 22, 0.95)",
                wire_width=3.0,
                opacity=0.30,
                name=f"Search Zone S({target_b.id})",
                flatshading=True
            )
            h_mesh.legendgroup = "highlight_zone"
            h_wire.legendgroup = "highlight_zone"
            h_mesh.showlegend = True
            fig.add_trace(h_mesh)
            fig.add_trace(h_wire)

    # 5. Nondominated Points N (Pareto)
    if pts:
        px = [p["coordinates"][0] for p in pts]
        py = [p["coordinates"][1] for p in pts]
        pz = [p["coordinates"][2] for p in pts]
        p_text = [
            f"<b>Point {escape(p['id'])}</b><br>Coords: ({p['coordinates'][0]:.2f}, {p['coordinates'][1]:.2f}, {p['coordinates'][2]:.2f})"
            for p in pts
        ]
        fig.add_trace(go.Scatter3d(
            x=px, y=py, z=pz,
            mode="markers+text",
            marker=dict(size=8, color="#ef4444", symbol="circle", line=dict(color="#ffffff", width=1.5)),
            text=[escape(p["id"]) for p in pts],
            textposition="top right",
            hovertext=p_text,
            hoverinfo="text",
            name="Points N (Pareto)",
            legendgroup="points",
            showlegend=True
        ))

    # 6. Local Bounds U(N) or L(N)
    bx = [b.coordinates[0] for b in b_records]
    by = [b.coordinates[1] for b in b_records]
    bz = [b.coordinates[2] for b in b_records]

    b_colors = []
    b_sizes = []
    for b in b_records:
        if highlight_bound_id == b.id:
            b_colors.append("#f97316")  # bright orange
            b_sizes.append(11)
        elif b.is_extreme:
            b_colors.append("#a855f7")  # purple
            b_sizes.append(8)
        else:
            b_colors.append("#38bdf8")  # cyan
            b_sizes.append(7)

    bounds_trace_name = "Local Bounds L(N)" if is_max else "Local Bounds U(N)"
    b_text = []
    for b in b_records:
        def_pts_str = ", ".join([escape(dp["id"]) for dp in b.defining_points])
        nbs_str = ", ".join([
            f"ν_{k+1}={b_records[nb].id if nb is not None and nb < len(b_records) else '∅'}"
            for k, nb in enumerate(b.neighbors)
        ])
        hover_info = (
            f"<b>Bound {escape(b.id)}</b><br>"
            f"Coords: ({b.coordinates[0]:.2f}, {b.coordinates[1]:.2f}, {b.coordinates[2]:.2f})<br>"
            f"Defining Points: [{def_pts_str}]<br>"
            f"Neighbors: {nbs_str}<br>"
            f"Status: {'Extreme Local Bound' if b.is_extreme else 'Interior Bound'}"
        )
        b_text.append(hover_info)

    fig.add_trace(go.Scatter3d(
        x=bx, y=by, z=bz,
        mode="markers+text",
        marker=dict(size=b_sizes, color=b_colors, symbol="diamond", line=dict(color="#0f172a", width=1.5)),
        text=[escape(b.id) for b in b_records],
        customdata=[b.id for b in b_records],
        textposition="top center",
        textfont=dict(color="#ffffff", size=11),
        hovertext=b_text,
        hoverinfo="text",
        name=bounds_trace_name,
        legendgroup="local_bounds",
        showlegend=True
    ))

    # 7. Reference and Anti-reference points
    ref_labels = ["Reference m", "Anti-reference M"] if is_max else ["Anti-reference m", "Reference M"]
    ref_colors = ["#94a3b8", "#10b981"] if is_max else ["#10b981", "#94a3b8"]
    fig.add_trace(go.Scatter3d(
        x=[lb_pt[0], ub_pt[0]],
        y=[lb_pt[1], ub_pt[1]],
        z=[lb_pt[2], ub_pt[2]],
        mode="markers+text",
        marker=dict(size=7, color=ref_colors, symbol="square"),
        text=ref_labels,
        textposition="top left",
        textfont=dict(color="#e2e8f0", size=11),
        name="Reference Point",
        legendgroup="reference_points",
        showlegend=True
    ))

    # 8. Perspective & Camera (matching Klamroth et al. 2015 Figure 2)
    # In MINIMIZATION, view from Ideal m towards Nadir M:
    # eye at negative coordinates (-1.65, -1.65, 1.35), so m is in foreground,
    # dominated boxes extend into depth, and search zones face towards the viewer.
    if not is_max:
        camera_eye = dict(x=-1.65, y=-1.65, z=1.35)
    else:
        camera_eye = dict(x=1.65, y=1.65, z=1.35)

    camera_settings = dict(
        eye=camera_eye,
        center=dict(x=0, y=0, z=-0.05),
        up=dict(x=0, y=0, z=1)
    )

    # Readable, high-contrast dark-mode legend
    fig.update_layout(
        scene=dict(
            camera=camera_settings,
            xaxis=dict(title="Objective f₁ (z₁)", backgroundcolor="#0f172a", gridcolor="#334155", color="#e2e8f0"),
            yaxis=dict(title="Objective f₂ (z₂)", backgroundcolor="#0f172a", gridcolor="#334155", color="#e2e8f0"),
            zaxis=dict(title="Objective f₃ (z₃)", backgroundcolor="#0f172a", gridcolor="#334155", color="#e2e8f0"),
            aspectmode="data"
        ),
        paper_bgcolor="#090d16",
        plot_bgcolor="#090d16",
        font=dict(color="#f8fafc", size=12),
        margin=dict(l=0, r=0, b=45, t=25),
        legend=dict(
            orientation="h",
            x=0.5,
            xanchor="center",
            y=-0.06,
            yanchor="top",
            bgcolor="rgba(15, 23, 42, 0.85)",
            bordercolor="rgba(51, 65, 85, 0.6)",
            borderwidth=1,
            font=dict(color="#f8fafc", size=11),
            itemclick="toggle",
            itemdoubleclick="toggleothers"
        )
    )

    return fig


# ---------------------------------------------------------------------------
# 2D Visualization
# ---------------------------------------------------------------------------

def plot_2d_bounds(bounds_data: Dict[str, Any], highlight_bound_id: Optional[str] = None) -> go.Figure:
    """
    Creates 2D visualization of points, local bounds, and search zones
    matching Fig. 1 from Paper 1 & Paper 2.
    """
    fig = go.Figure()
    lb_pt = bounds_data.get("lower_bound", bounds_data["anti_reference"])
    ub_pt = bounds_data.get("upper_bound", bounds_data["reference_point"])
    is_max = (bounds_data.get("sense") == "MAXIMIZE")
    b_records = bounds_data["bounds"]
    pts = bounds_data["points"]

    # Search space outline [LB, UB]
    fig.add_shape(
        type="rect",
        x0=lb_pt[0], y0=lb_pt[1], x1=ub_pt[0], y1=ub_pt[1],
        line=dict(color="#64748b", width=2, dash="dash"),
        fillcolor="rgba(30, 41, 59, 0.3)"
    )

    # Shaded search zones S(u): [LB, u] for min, [u, UB] for max
    for b in b_records:
        if not is_max:
            x0, y0 = lb_pt[0], lb_pt[1]
            x1, y1 = b.coordinates[0], b.coordinates[1]
        else:
            x0, y0 = b.coordinates[0], b.coordinates[1]
            x1, y1 = ub_pt[0], ub_pt[1]
        fig.add_trace(go.Scatter(
            x=[x0, x1, x1, x0, x0],
            y=[y0, y0, y1, y1, y0],
            fill="toself",
            fillcolor="rgba(249, 115, 22, 0.25)" if b.id == highlight_bound_id else "rgba(56, 189, 248, 0.15)",
            line=dict(color="#f97316" if b.id == highlight_bound_id else "rgba(56, 189, 248, 0.5)", width=1.5),
            mode="lines",
            hoverinfo="skip",
            showlegend=False,
            legendgroup="local_bounds"
        ))

    # Local bounds
    bx = [b.coordinates[0] for b in b_records]
    by = [b.coordinates[1] for b in b_records]
    b_hover = [
        f"<b>Bound {escape(b.id)}</b><br>({b.coordinates[0]:.2f}, {b.coordinates[1]:.2f})"
        for b in b_records
    ]
    bounds_name = "Local Bounds L(N)" if is_max else "Local Bounds U(N)"
    fig.add_trace(go.Scatter(
        x=bx, y=by,
        mode="markers+text",
        marker=dict(size=[12 if b.id == highlight_bound_id else 10 for b in b_records],
                    color=["#f97316" if b.id == highlight_bound_id else "#38bdf8" for b in b_records], symbol="diamond", line=dict(color="#0f172a", width=1.5)),
        text=[escape(b.id) for b in b_records],
        customdata=[b.id for b in b_records],
        textposition="top center",
        hovertext=b_hover,
        hoverinfo="text",
        name=bounds_name,
        legendgroup="local_bounds"
    ))

    # Points
    if pts:
        px = [p["coordinates"][0] for p in pts]
        py = [p["coordinates"][1] for p in pts]
        p_hover = [
            f"<b>Point {escape(p['id'])}</b><br>({p['coordinates'][0]:.2f}, {p['coordinates'][1]:.2f})"
            for p in pts
        ]
        fig.add_trace(go.Scatter(
            x=px, y=py,
            mode="markers+text",
            marker=dict(size=10, color="#ef4444", symbol="circle", line=dict(color="#ffffff", width=1.5)),
            text=[escape(p["id"]) for p in pts],
            textposition="bottom left",
            hovertext=p_hover,
            hoverinfo="text",
            name="Points N (Pareto)",
            legendgroup="points"
        ))

    fig.update_layout(
        xaxis=dict(title="Objective f₁ (z₁)", range=[lb_pt[0] - 0.5, ub_pt[0] + 0.5], gridcolor="#334155"),
        yaxis=dict(title="Objective f₂ (z₂)", range=[lb_pt[1] - 0.5, ub_pt[1] + 0.5], gridcolor="#334155"),
        paper_bgcolor="#090d16",
        plot_bgcolor="#0c1220",
        font=dict(color="#f8fafc", size=12),
        margin=dict(l=40, r=40, b=40, t=50),
        legend=dict(
            x=0.02, y=0.98,
            bgcolor="rgba(30, 41, 59, 0.92)",
            bordercolor="rgba(100, 116, 139, 0.5)",
            borderwidth=1,
            font=dict(color="#f8fafc")
        )
    )
    return fig


# ---------------------------------------------------------------------------
# Centered Pairwise Projections Matrix (3D)
# ---------------------------------------------------------------------------

def plot_pairwise_projections_2d(
    bounds_data: Dict[str, Any],
    highlight_bound_id: Optional[str] = None
) -> go.Figure:
    """
    Creates a centered matrix of pairwise 2D projections for 3D spaces:
    (f1 vs f2), (f1 vs f3), and (f2 vs f3).
    """
    dims = bounds_data["dimensions"]
    pairs = [(i, j) for i in range(dims) for j in range(i + 1, dims)]
    n_pairs = len(pairs)

    cols = n_pairs
    rows = 1

    subplot_titles = [f"Plane (f_{i+1}, f_{j+1})" for i, j in pairs]
    fig = make_subplots(
        rows=rows, cols=cols,
        subplot_titles=subplot_titles,
        horizontal_spacing=0.08
    )

    lb_pt = bounds_data.get("lower_bound", bounds_data["anti_reference"])
    ub_pt = bounds_data.get("upper_bound", bounds_data["reference_point"])
    is_max = (bounds_data.get("sense") == "MAXIMIZE")
    b_records = bounds_data["bounds"]
    pts = bounds_data["points"]

    for idx, (di, dj) in enumerate(pairs):
        c = idx + 1

        # Draw 2D bounding boxes for each bound as filled traces in legendgroup="local_bounds"
        for b in b_records:
            is_hl = (highlight_bound_id == b.id)
            if not is_max:
                x0, y0 = lb_pt[di], lb_pt[dj]
                x1, y1 = b.coordinates[di], b.coordinates[dj]
            else:
                x0, y0 = b.coordinates[di], b.coordinates[dj]
                x1, y1 = ub_pt[di], ub_pt[dj]

            box_fill = "rgba(249, 115, 22, 0.25)" if is_hl else "rgba(56, 189, 248, 0.12)"
            box_line = "rgba(249, 115, 22, 0.95)" if is_hl else "rgba(56, 189, 248, 0.4)"
            line_w = 2.5 if is_hl else 1

            fig.add_trace(
                go.Scatter(
                    x=[x0, x1, x1, x0, x0],
                    y=[y0, y0, y1, y1, y0],
                    fill="toself",
                    fillcolor=box_fill,
                    line=dict(color=box_line, width=line_w),
                    mode="lines",
                    hoverinfo="skip",
                    showlegend=False,
                    legendgroup="local_bounds"
                ),
                row=1, col=c
            )

        # Plot bounds
        b_colors = ["#f97316" if highlight_bound_id == b.id else "#38bdf8" for b in b_records]
        b_sizes = [11 if highlight_bound_id == b.id else 7 for b in b_records]
        fig.add_trace(
            go.Scatter(
                x=[b.coordinates[di] for b in b_records],
                y=[b.coordinates[dj] for b in b_records],
                mode="markers+text",
                marker=dict(size=b_sizes, color=b_colors, symbol="diamond", line=dict(color="#090d16", width=1.5)),
                text=[escape(b.id) for b in b_records],
                customdata=[b.id for b in b_records],
                textposition="top center",
                showlegend=(idx == 0),
                name="Local Bounds",
                legendgroup="local_bounds"
            ),
            row=1, col=c
        )

        # Plot points
        if pts:
            fig.add_trace(
                go.Scatter(
                    x=[p["coordinates"][di] for p in pts],
                    y=[p["coordinates"][dj] for p in pts],
                    mode="markers+text",
                    marker=dict(size=8, color="#ef4444", symbol="circle", line=dict(color="#ffffff", width=1.2)),
                    text=[escape(p["id"]) for p in pts],
                    textposition="bottom right",
                    showlegend=(idx == 0),
                    name="Points N",
                    legendgroup="points"
                ),
                row=1, col=c
            )

        fig.update_xaxes(
            title_text=f"Objective f_{di+1}",
            gridcolor="#334155",
            zerolinecolor="#475569",
            constrain="domain",
            row=1, col=c
        )
        fig.update_yaxes(
            title_text=f"Objective f_{dj+1}",
            gridcolor="#334155",
            zerolinecolor="#475569",
            scaleanchor=f"x{'' if c == 1 else c}",
            scaleratio=1,
            row=1, col=c
        )

    fig.update_layout(
        paper_bgcolor="#090d16",
        plot_bgcolor="#0c1220",
        font=dict(color="#f8fafc", size=11),
        height=390,
        margin=dict(l=45, r=45, b=65, t=45),
        legend=dict(
            orientation="h",
            xanchor="center",
            x=0.5,
            y=-0.22,
            bgcolor="rgba(15, 23, 42, 0.85)",
            bordercolor="rgba(51, 65, 85, 0.6)",
            borderwidth=1,
            font=dict(color="#f8fafc", size=10),
            itemclick="toggle",
            itemdoubleclick="toggleothers"
        )
    )
    return fig


# ---------------------------------------------------------------------------
# Neighbor Graph Extraction & Plotting (Well-Spaced Layout)
# ---------------------------------------------------------------------------

def extract_network_graph_elements(
    bounds_data: Dict[str, Any],
    mode: str = "combined"  # 'combined' or 'undirected'
) -> Dict[str, Any]:
    """
    Extracts nodes and edges formatted for Vis.js with generous spacing.
    """
    component_colors = [
        "#38bdf8",  # sky blue for k=1
        "#34d399",  # emerald green for k=2
        "#fbbf24"   # amber for k=3
    ]

    nodes = []
    for b in bounds_data["bounds"]:
        coords_str = ", ".join([f"{c:.2f}" for c in b.coordinates])
        label = f"{b.id}\n({coords_str})"
        title = (
            f"<b>Bound {escape(b.id)}</b><br>"
            f"Coordinates: [{coords_str}]<br>"
            f"Status: {'Quasi (redundant search zone)' if b.is_quasi else 'Extreme Local Bound' if b.is_extreme else 'Interior'}"
        )
        nodes.append({
            "id": b.index,
            "label": label,
            "title": title,
            "bound_id": b.id,
            "is_quasi": b.is_quasi,
            "shape": "ellipse" if b.is_quasi else "box",
            "margin": 12,
            "color": {
                "background": "#422006" if b.is_quasi else "#312e81" if b.is_extreme else "#1e293b",
                "border": "#fbbf24" if b.is_quasi else "#818cf8" if b.is_extreme else "#38bdf8",
                "highlight": {"background": "#4338ca", "border": "#f97316"}
            },
            "font": {"color": "#f8fafc", "face": "monospace", "size": 13, "bold": True}
        })

    edges = []
    edge_id = 0

    if mode == "undirected":
        for u, v in bounds_data["edges_undirected"]:
            edges.append({
                "id": f"e_{edge_id}",
                "from": u,
                "to": v,
                "arrows": "",
                "color": {"color": "#94a3b8", "highlight": "#f97316"},
                "width": 2.5
            })
            edge_id += 1
    else:  # combined multi-graph
        for u, v, k in bounds_data["edges_directed"]:
            c = component_colors[k % len(component_colors)]
            edges.append({
                "id": f"e_{edge_id}",
                "from": u,
                "to": v,
                "arrows": "to",
                "label": f"ν_{k+1}",
                "font": {"align": "top", "size": 12, "color": c, "strokeWidth": 0},
                "color": {"color": c, "highlight": "#f97316"},
                "width": 2.2,
                "smooth": {
                    "type": "curvedCW",
                    "roundness": 0.20 + (k * 0.16)
                }
            })
            edge_id += 1

    return {"nodes": nodes, "edges": edges}


def plot_neighbor_graph_plotly(bounds_data: Dict[str, Any], mode: str = "combined") -> go.Figure:
    """
    Renders neighbor graph using Plotly with wide node spacing.
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx is required for graph layout generation.")

    dims = bounds_data["dimensions"]
    component_colors = ["#38bdf8", "#34d399", "#fbbf24"]

    G = nx.DiGraph() if mode != "undirected" else nx.Graph()
    for b in bounds_data["bounds"]:
        G.add_node(b.index, label=b.id, coords=b.coordinates, is_extreme=b.is_extreme)

    if mode == "undirected":
        for u, v in bounds_data["edges_undirected"]:
            G.add_edge(u, v)
    else:
        for u, v, k in bounds_data["edges_directed"]:
            G.add_edge(u, v, component=k)

    # Wide spacing layout: larger k multiplier
    pos = nx.spring_layout(G, seed=42, k=3.2 / math.sqrt(max(1, len(G.nodes))), iterations=70)

    fig = go.Figure()

    if mode == "undirected":
        edge_x, edge_y = [], []
        for u, v in G.edges():
            x0, y0 = pos[u]
            x1, y1 = pos[v]
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])
        fig.add_trace(go.Scatter(
            x=edge_x, y=edge_y,
            mode="lines",
            line=dict(width=2, color="#94a3b8"),
            hoverinfo="none",
            name="Neighbor Relation ν"
        ))
    else:
        for k in range(dims):
            k_edges = [(u, v) for u, v, data in G.edges(data=True) if data.get("component") == k]
            ex, ey = [], []
            for u, v in k_edges:
                x0, y0 = pos[u]
                x1, y1 = pos[v]
                ex.extend([x0, x1, None])
                ey.extend([y0, y1, None])
            if ex:
                c = component_colors[k % len(component_colors)]
                fig.add_trace(go.Scatter(
                    x=ex, y=ey,
                    mode="lines",
                    line=dict(width=2.5, color=c),
                    name=f"k-neighbor ν_{k+1}"
                ))

    node_x = [pos[node][0] for node in G.nodes()]
    node_y = [pos[node][1] for node in G.nodes()]
    node_text = [G.nodes[node]["label"] for node in G.nodes()]
    node_colors = ["#818cf8" if G.nodes[node]["is_extreme"] else "#38bdf8" for node in G.nodes()]

    hover_info = []
    for node in G.nodes():
        b = bounds_data["bounds"][node]
        c_str = ", ".join([f"{c:.2f}" for c in b.coordinates])
        hover_info.append(f"<b>{escape(b.id)}</b><br>Coords: [{c_str}]<br>{'Extreme' if b.is_extreme else 'Internal'}")

    fig.add_trace(go.Scatter(
        x=node_x, y=node_y,
        mode="markers+text",
        marker=dict(size=28, color=node_colors, line=dict(color="#0f172a", width=2)),
        text=node_text,
        textposition="middle center",
        textfont=dict(color="#0f172a", size=12, family="monospace"),
        hovertext=hover_info,
        hoverinfo="text",
        name="Local Bounds"
    ))

    fig.update_layout(
        showlegend=True,
        paper_bgcolor="#0f172a",
        plot_bgcolor="#0f172a",
        font=dict(color="#f8fafc", size=12),
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        margin=dict(l=20, r=20, b=20, t=40),
        legend=dict(
            bgcolor="rgba(30, 41, 59, 0.92)",
            bordercolor="rgba(100, 116, 139, 0.4)",
            borderwidth=1,
            font=dict(color="#f8fafc")
        )
    )
    return fig


# ---------------------------------------------------------------------------
# Step-by-Step History Tracker
# ---------------------------------------------------------------------------

class GenerationStep:
    """Snapshot of bounds state at step t."""
    def __init__(
        self,
        step: int,
        point_added: Optional[lb.Point],
        bounds_data: Dict[str, Any],
        destroyed_bounds: List[str],
        new_bounds: List[str]
    ):
        self.step = step
        self.point_added = point_added
        self.bounds_data = bounds_data
        self.destroyed_bounds = destroyed_bounds
        self.new_bounds = new_bounds


class LocalBoundsTracker:
    """
    Maintains history of local bounds generation as points are added one by one.
    """
    def __init__(
        self,
        reference_point: List[float],
        anti_reference: List[float],
        sense: Union[lb.Objective, str] = lb.Objective.MINIMIZE,
        lower_bound: Optional[List[float]] = None,
        upper_bound: Optional[List[float]] = None
    ):
        sense = lb._normalize_sense(sense)
        self.sense = sense
        self.reference_point = reference_point
        self.anti_reference = anti_reference
        is_min = (sense == lb.Objective.MINIMIZE or sense == "MINIMIZE")
        if lower_bound is None or upper_bound is None:
            if is_min:
                self.lower_bound = list(anti_reference)
                self.upper_bound = list(reference_point)
            else:
                self.lower_bound = list(reference_point)
                self.upper_bound = list(anti_reference)
        else:
            self.lower_bound = list(lower_bound)
            self.upper_bound = list(upper_bound)

        self.points_history: List[lb.Point] = []
        self.steps: List[GenerationStep] = []
        
        # Step 0
        nbs = lb.NeighborhoodBoundSet(reference_point, anti_reference, sense=sense)
        b_data = extract_bounds_data(
            nbs, reference_point, anti_reference, [],
            lower_bound=self.lower_bound, upper_bound=self.upper_bound, sense=self.sense
        )
        self.steps.append(GenerationStep(
            step=0,
            point_added=None,
            bounds_data=b_data,
            destroyed_bounds=[],
            new_bounds=[b.id for b in b_data["bounds"]]
        ))

    def add_point(self, point: lb.Point) -> GenerationStep:
        """Adds a point and computes step changes."""
        validate_point(list(point.coordinates), point.id, self.points_history,
                       self.lower_bound, self.upper_bound, self.sense)
        points = self.points_history + [point]
        nbs = lb.NeighborhoodBoundSet(self.reference_point, self.anti_reference, sense=self.sense)
        
        prev_bound_ids = set(b.id for b in self.steps[-1].bounds_data["bounds"])
        
        for p in points:
            nbs.update(p)
            
        curr_b_data = extract_bounds_data(
            nbs, self.reference_point, self.anti_reference, points,
            lower_bound=self.lower_bound, upper_bound=self.upper_bound, sense=self.sense
        )
        curr_bound_ids = set(b.id for b in curr_b_data["bounds"])
        
        destroyed = sorted(prev_bound_ids - curr_bound_ids)
        created = sorted(curr_bound_ids - prev_bound_ids)
        
        step = GenerationStep(
            step=len(self.steps),
            point_added=point,
            bounds_data=curr_b_data,
            destroyed_bounds=destroyed,
            new_bounds=created
        )
        self.points_history.append(point)
        self.steps.append(step)
        return step

    def get_step(self, step_idx: int) -> GenerationStep:
        if 0 <= step_idx < len(self.steps):
            return self.steps[step_idx]
        raise IndexError(f"Step {step_idx} out of range (0..{len(self.steps)-1})")
