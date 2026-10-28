# Local Bounds Visualization

Explore local bounds in two or three objectives through an interactive web dashboard, standalone HTML reports, or Python figures. The tool uses `NeighborhoodBoundSet` (Algorithm 1 from Dächert et al., 2017) to expose search regions, defining points, and neighborhood relationships. Presets reproduce examples from both papers referenced in the [library README](../README.md).

## Contents

- [Installation](#installation)
- [Interactive dashboard](#interactive-dashboard)
- [Reading the views](#reading-the-views)
- [Input rules](#input-rules)
- [Presets](#presets)
- [Command-line reports](#command-line-reports)
- [Python API](#python-api)
- [Dashboard HTTP API](#dashboard-http-api)
- [Runtime behavior and limitations](#runtime-behavior-and-limitations)
- [Troubleshooting](#troubleshooting)
- [Files](#files)

## Installation

Run the commands in this guide from the **repository root**, which contains `pyproject.toml`, rather than from inside `visualization/`.

The Python package requires Python 3.8+, a C++17 compiler, and CMake 3.15+. Visualization dependencies are Flask, Plotly, NetworkX, and NumPy; their declared versions are in [requirements.txt](requirements.txt) and the `vis` extra in [pyproject.toml](../pyproject.toml).

```bash
python3 -m venv venv
source venv/bin/activate
python -m pip install -e ".[vis]"
```

If you already have a virtual environment, activate it and run the install command. An editable install loads the Python wrapper from this checkout. **Rerun the install command after changing C++ headers or bindings** to rebuild the native extension.

For a regular installation, use `python -m pip install ".[vis]"` instead. Installing `visualization/requirements.txt` alone installs only visualization dependencies; it does not install or rebuild `local_bounds`.

## Interactive dashboard

```bash
python -m visualization.app --port 8050
```

Open [http://127.0.0.1:8050](http://127.0.0.1:8050). The initial problem is the 3D minimization `paper2` preset, with all three points inserted.

The dashboard launcher accepts:

| Option | Default | Purpose |
|---|---|---|
| `--port` | `8050` | Server port |
| `--host` | `127.0.0.1` | Server bind address |
| `--preset` | `paper2` | Initial problem; see [Presets](#presets) |

For example, start with tied coordinates on a different port:

```bash
python -m visualization.app --preset paper1_ngp --port 8051
```

### Explore an example

1. Select a preset. Changing the selection loads it; **Load** reloads the selected preset.
2. Move the timeline to step 0 to see the initial bound before any insertions.
3. Advance one step at a time, or use playback, to inspect created and destroyed bounds.
4. Select a bound in the table, spatial plot, or graph to highlight it across views.
5. Switch graph topology or edge mode to inspect neighborhood relationships.
6. Enter coordinates in **Probe search region** to test membership at the selected step without inserting a point. The result identifies a nonredundant containing bound, when one exists, and shows the domain and search-zone inequalities.

### Create a custom instance

1. Choose **2D** or **3D**, then **Minimization** or **Maximization**.
2. Enter comma-separated lower and upper interval bounds, with one value per objective, and apply them.
3. Click **New custom instance** to start empty with the displayed settings.
4. Enter a point ID and comma-separated coordinates, then add the point. An empty ID receives an unused `z1`, `z2`, … identifier.
5. Add more mutually nondominating points, delete points if needed, and inspect their insertion history.

**Applying settings or loading a preset clears the current points and resets the views.** Changing dimension also resets the interval to 0–10 in each objective. Invalid point submissions appear as inline errors without inserting the point.

## Reading the views

### Spatial and pairwise plots

- **3D view:** local bounds and input points inside the configured box. Dominated zones have distinct colors and 30% opacity by default. For minimization, each point's dominance cone extends toward the upper reference point, as in Figure 2 of Klamroth et al. (2015); for maximization, it extends toward the lower reference point. The default minimization camera looks from the ideal side toward the reference side.
- **Pairwise view:** a centered matrix of the three projections `(f1, f2)`, `(f1, f3)`, and `(f2, f3)` for a 3D problem, with projected search boxes. Coincident markers are grouped; selecting a bound group offers its underlying bounds for selection.
- **2D view:** search rectangles and the staircase Pareto front.

Use **Dominated Zones D(N)** and **Search Zones S(u)** to toggle the 3D regions. The opacity slider adjusts dominated zones; **Wireframe** removes filled 3D faces. Legend entries toggle trace groups, including the dominated zones as a group.

**Labels** shows all point and bound labels by default or hides them. Labels use white boxes with dark 12 px text; `m` and `M` mark the lower and upper reference corners. Camera, zoom, and surviving graph-node positions persist while stepping or highlighting. Spatial and pairwise views keep separate zoom settings.

### Bounds table and neighborhood graph

The table lists each bound's ID, coordinates `u`, defining points `z^j(u)`, component neighbors `ν_k(u)`, and status. Dummy defining points are displayed as `z^1`, `z^2`, and `z^3`. Table selections highlight the corresponding search zone and graph node.

The graph offers two topologies:

- **Contracted:** displays nonredundant bounds, contracting paths through quasi-nonredundant nodes.
- **Raw:** displays Algorithm 1's active topology, including quasi-nonredundant bounds as labeled ellipses.

Both show total, quasi-nonredundant, and nonredundant counts. Use component arrows for the coordinate-specific `ν_k` links, or undirected edges for adjacency. Graph controls also fit the view and toggle layout physics.

Quasi-nonredundant bounds can occur when points share coordinate values. Algorithm 1 retains them to preserve graph connectivity. The default spatial plots and table use the nonredundant set, which matches the bound coordinates computed by Algorithms 2–5. More raw nodes therefore do not imply a larger nonredundant search region.

## Input rules

Each point must have the configured dimension and finite numeric coordinates. Every lower interval bound must be strictly smaller than its upper bound.

| Objective sense | Reference point | Anti-reference point | Allowed point coordinates |
|---|---|---|---|
| `MINIMIZE` | Upper bound | Lower bound | `lower[j] <= z[j] < upper[j]` |
| `MAXIMIZE` | Lower bound | Upper bound | `lower[j] < z[j] <= upper[j]` |

The anti-reference boundary is included; the reference boundary is excluded. In particular, a coordinate equal to 0 can be valid in the default minimization problem, while a coordinate equal to 10 is not.

Point IDs and coordinate vectors must be unique. The input set `N` must remain **stable**, meaning its points are mutually nondominating. The tool rejects both a point dominated by an existing point and a point that dominates an existing point. It does not automatically replace dominated points. Delete conflicting points first if you want to change the set.

Ties in individual coordinates are allowed. They differ from duplicate coordinate vectors, which are rejected. The membership probe checks a query without adding it, so insertion-only restrictions such as unique IDs do not apply to probes.

## Presets

All presets use lower bounds of 0 and upper bounds of 10 in every objective. Loading one replaces the current dimension, objective sense, interval, and points.

| Name | Dimension | Sense | Points, in insertion order |
|---|---|---|---|
| `paper2` | 3 | MINIMIZE | `(4,0,4)`, `(3,3,1)`, `(2,2,2)` — Dächert et al. (2017), Example 2.8 |
| `paper1_sa` | 3 | MINIMIZE | `(3,5,7)`, `(6,2,4)` — Klamroth et al. (2015), Example 2, general position (SA) |
| `paper1_ngp` | 3 | MINIMIZE | `(2,7,7)`, `(5,7,5)`, `(8,7,3)`, `(4,3,7)` — Example 3, including the later insertion, with coordinate ties |
| `2d` | 2 | MINIMIZE | `(3,7)`, `(5,4)`, `(7,2)` |
| `max_3d` | 3 | MAXIMIZE | `(6,10,6)`, `(7,7,9)`, `(8,8,8)` |

The shared definitions live in [presets.py](presets.py), so dashboard and CLI presets use the same data.

## Command-line reports

The CLI prints a final bounds table and writes a standalone HTML file containing the spatial plot, neighborhood graph, and bounds table. A 3D report also includes pairwise projections. Reports embed Plotly and work offline; open the generated file directly in a browser.

```bash
# Paper 2 example, with per-insertion details in the terminal
python -m visualization.cli --preset paper2 --step-by-step --output report_3d.html

# 2D example
python -m visualization.cli --preset 2d --step-by-step --output report_2d.html

# Maximization example
python -m visualization.cli --preset max_3d --output report_max.html

# Custom minimization points and interval
python -m visualization.cli --dim 3 --points '[[4,0,4],[3,3,1],[2,2,2]]' --ref '10,10,10' --anti '0,0,0' --output custom.html

# Launch the dashboard instead of writing a report
python -m visualization.cli --interactive --preset paper1_ngp --port 8050
```

### CLI options

| Option | Default | Meaning |
|---|---|---|
| `--dim {2,3}` | Inferred | Coordinate dimension, inferred from the first point; 3 for an empty list |
| `--preset NAME` | `paper2` | Selects preset points **and objective sense** |
| `--points JSON` | Preset points | Replaces the points with a JSON list of coordinate lists |
| `--ref VALUES` | 10 per objective for MINIMIZE; 0 for MAXIMIZE | Comma-separated reference coordinates |
| `--anti VALUES` | 0 per objective for MINIMIZE; 10 for MAXIMIZE | Comma-separated anti-reference coordinates |
| `--output PATH` | `local_bounds_report.html` | Output HTML path; its parent directory must exist |
| `--step-by-step` | Off | Prints destroyed/created bound IDs and counts for each insertion |
| `--interactive` | Off | Starts the dashboard using the selected preset |
| `--port PORT` | `8050` | Dashboard port when using `--interactive` |

There is no separate `--sense` flag. Use `--preset max_3d` with `--points` for custom maximization data. `--dim` does not convert preset coordinates: use the `2d` preset or supply 2D points for a 2D problem.

`--step-by-step` changes terminal output and retains intermediate snapshots; the HTML report still shows the final state. Without it, the CLI constructs only the final state. `--interactive` uses only the preset and port; configure custom points and interval bounds in the dashboard.

Run `python -m visualization.cli --help` or `python -m visualization.app --help` for launcher help. Existing output files are overwritten.

## Python API

### Extract data and export figures

```python
import local_bounds as lb
import visualization.engine as vis

reference = [10.0, 10.0, 10.0]
anti_reference = [0.0, 0.0, 0.0]
points = [
    lb.Point("z1", [4.0, 0.0, 4.0]),
    lb.Point("z2", [3.0, 3.0, 1.0]),
    lb.Point("z3", [2.0, 2.0, 2.0]),
]
nbs = lb.NeighborhoodBoundSet(reference, anti_reference)
for point in points:
    nbs.update(point)

# Pass the input points so plots can draw them and their dominated zones.
data = vis.extract_bounds_data(nbs, reference, anti_reference, points)
table = vis.create_bounds_table(data)

fig_3d = vis.plot_3d_bounds(data, show_occupied_boxes=True)
fig_pairwise = vis.plot_pairwise_projections_2d(data)
fig_graph = vis.plot_neighbor_graph_plotly(data, mode="combined")

fig_3d.write_html("bounds_3d.html")
fig_pairwise.write_html("pairwise.html")
fig_graph.write_html("graph.html")
```

Each plot function returns a Plotly figure. Use `.show()` to open it interactively or `.write_html(...)` to save it. For 2D data, use `vis.plot_2d_bounds(data)`.

`extract_bounds_data()` includes the dimension, sense, interval, input points, bound records, directed and undirected edges, and bound counts. Bound records expose IDs, coordinates, defining points, neighbor indices, and extreme/quasi status. `create_bounds_table()` turns those records into a list of row dictionaries.

For maximization, construct `NeighborhoodBoundSet` with the lower reference and upper anti-reference points, and pass `sense=lb.Objective.MAXIMIZE` to both the constructor and `extract_bounds_data()`.

To inspect raw topology, reuse the variables above:

```python
raw_data = vis.extract_bounds_data(
    nbs, reference, anti_reference, points, include_quasi=True
)
raw_graph = vis.plot_neighbor_graph_plotly(raw_data, mode="undirected")
raw_graph.write_html("raw_graph.html")

# Native graph export also exposes a per-node quasi flag.
native_graph = nbs.get_adjacency_graph(include_quasi=True)
print(list(native_graph.quasi))
```

Use the `paper1_ngp` points to see quasi-nonredundant nodes. `vis.extract_network_graph_elements(data, mode="combined")` returns node and edge dictionaries for a Vis-Network graph, as used by the dashboard.

### Track insertion history

`LocalBoundsTracker` validates each insertion and stores an independent bounds snapshot for every step, including step 0 before any points are added:

```python
from visualization.engine import LocalBoundsTracker

tracker = LocalBoundsTracker(reference, anti_reference)
for point in points:
    step = tracker.add_point(point)
    print(step.step, step.destroyed_bounds, step.new_bounds)

initial = tracker.get_step(0).bounds_data
final = tracker.get_step(len(points)).bounds_data
```

For direct `NeighborhoodBoundSet` usage, call `vis.validate_point(coords, point_id, existing_points, lower_bound, upper_bound, sense)` before inserting untrusted input if you need the visualizer's stable-set validation. Data extraction itself does not validate the original point sequence.

## Dashboard HTTP API

These routes operate on the same in-memory problem displayed in the browser. POST bodies must be JSON objects. Invalid input returns HTTP 400 with `error` and `code` fields.

| Route | Purpose |
|---|---|
| `GET /api/state` | Settings, points, history, bounds counts, and table for the selected step |
| `POST /api/configure` | Start empty with `dimensions`, `sense`, `lower_bound`, and `upper_bound` |
| `POST /api/load_preset` | Replace the problem using `{"preset": "paper2"}` |
| `POST /api/add_point` | Insert `{"id": "z4", "coordinates": [1.5, 5, 1.5]}`; ID is optional |
| `POST /api/delete_point` | Delete a point using its zero-based `index` |
| `POST /api/probe` | Query membership using `{"coordinates": [1, 1, 1]}` |
| `GET /api/graph` | Graph nodes, edges, and counts; accepts `view=contracted` or `raw`, and `mode=combined` or `undirected` |
| `GET /api/figure_3d` | Plotly JSON for the spatial view; also handles 2D problems |
| `GET /api/figure_pairwise` | Plotly JSON for pairwise projections |

State, graph, figure, and probe requests accept `?step=N`; omitted steps use the current step, and integer values are clamped to the available history. Figure routes also accept `highlight=BOUND_ID` and `labels=all` or `none`. The spatial route accepts `show_occupied=0/1`, `show_boxes=0/1`, `wireframe=0/1`, and a finite `dom_opacity` clamped to `[0, 1]`.

## Runtime behavior and limitations

- **Dimensions:** plotting and data extraction support only 2D and 3D, even though the underlying library supports more objectives.
- **Shared session:** one server process holds one problem in memory. Browser tabs and clients share it; restarting the server reloads the startup preset. There is no session save/load feature.
- **Network access:** reports work offline. The dashboard serves Plotly locally but loads Tailwind, Lucide, and Vis-Network from CDNs, so those dashboard assets require internet access.
- **History and caches:** dashboard history is computed from bound IDs without unused graph snapshots. It caches up to eight step snapshots and eight serialized figures, keyed by problem settings, ordered point IDs/coordinates, and view options. API mutations clear the caches; public step data is copied so callers cannot alter cached snapshots.
- **Local use:** the default launcher binds to `127.0.0.1` and runs Flask with debug disabled. API requests serialize access to the shared state. The dashboard is a research tool with no authentication or per-user isolation.

## Troubleshooting

| Symptom | Check or fix |
|---|---|
| `No module named local_bounds` or a missing visualization dependency | Activate the intended virtual environment and run `python -m pip install -e ".[vis]"` from the repository root. |
| Missing library attributes or incompatible graph-export arguments | Rebuild the native extension with the same install command; an older installed extension can disagree with the checkout. |
| Relative import error when launching a script | Use `python -m visualization.app` or `python -m visualization.cli`, rather than running the `.py` file directly. |
| Port already in use | Choose another port, such as `--port 8051`, and open that address. |
| Dashboard lacks styling, icons, or its graph | Check access to the CDN assets. Standalone CLI reports are the offline option. |
| A point is rejected at an interval boundary | Check the sense-dependent inequalities in [Input rules](#input-rules); the reference boundary is excluded. |
| A new point dominates an existing point | Delete the conflicting point first; insertions must preserve a stable set. |
| `--dim 2` fails with a 3D preset | Select `--preset 2d` or provide 2D `--points`; the flag does not project coordinates. |
| Report output cannot be opened | Ensure its parent directory exists and is writable. |

## Files

| File | Responsibility |
|---|---|
| [app.py](app.py) | Flask server, shared problem state, request validation, and figure caches |
| [cli.py](cli.py) | CLI argument parsing, terminal tables, and standalone HTML reports |
| [engine.py](engine.py) | Point validation, bound extraction, plots, graph conversion, and insertion tracking |
| [presets.py](presets.py) | Shared example coordinates and objective senses |
| [templates/index.html](templates/index.html) | Dashboard layout, controls, Plotly rendering, and Vis-Network interactions |
| [requirements.txt](requirements.txt) | Visualization dependencies, without the core library |
| [../tests/test_visualization.py](../tests/test_visualization.py) | Visualization and dashboard integration tests |

To run the visualization tests in an environment with the library and visualization dependencies installed:

```bash
python -m unittest discover -s tests -p test_visualization.py
```
