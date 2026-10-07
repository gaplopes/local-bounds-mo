"""Shared example coordinates for the web dashboard and CLI."""

# (objective sense, coordinates); intervals are [0, 10] in each objective.
PRESETS = {
    # Dächert et al. (2017), Example 2.8.
    "paper2": ("MINIMIZE", [[4, 0, 4], [3, 3, 1], [2, 2, 2]]),
    # Klamroth et al. (2015), Examples 2 and 3 (including the later insertion).
    "paper1_sa": ("MINIMIZE", [[3, 5, 7], [6, 2, 4]]),
    "paper1_ngp": ("MINIMIZE", [[2, 7, 7], [5, 7, 5], [8, 7, 3], [4, 3, 7]]),
    "2d": ("MINIMIZE", [[3, 7], [5, 4], [7, 2]]),
    "max_3d": ("MAXIMIZE", [[6, 10, 6], [7, 7, 9], [8, 8, 8]]),
}
