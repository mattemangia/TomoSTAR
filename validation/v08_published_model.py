# Copyright 2026 Matteo Mangiagalli
# SPDX-License-Identifier: Apache-2.0
"""
The Vp model of the 2016-2017 central Italy example against the published 3-D Vp model of the
region: Carannante et al. (2013, JGR Solid Earth 118, 5391-5403; SIMULPS14 inversion of 2,099
hand-picked earthquakes, 15 km node spacing), from the authors' data set (Carannante, Cattaneo and
Monachesi 2025, doi:10.5281/zenodo.16535187, files plane.0000.xyz to plane.0012.xyz).

The two models come from different data (other years, other events and stations) and different
parameterisations, so their difference measures the data and the resolution, not the code; it is
reported for information. Compared: the published nodes inside the example's grid with resolution
(diagonal of the resolution matrix) of at least 0.5, where the example's ray coverage (DWS) is at
least 50 km, at 0, 4 and 8 km depth (the depth of 12 km has no node passing both conditions).

Inputs: TOMOSTAR_NORCIA_OUT, and CARANNANTE2025, the folder with the plane files of the data set.
"""
import json
import os
import struct

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from common import record

NORCIA = os.environ["TOMOSTAR_NORCIA_OUT"]
PUBLISHED = os.environ["CARANNANTE2025"]


def volume(path):
    with open(path, "rb") as f:
        b = f.read()
    n = struct.unpack("<i", b[8:12])[0]
    g = json.loads(b[12:12 + n])["Grid"]
    start = ((12 + n + 63) // 64) * 64
    v = np.frombuffer(b[start:start + 4 * g["Nx"] * g["Ny"] * g["Nz"]], "<f4").astype(float).reshape(g["Nz"], g["Ny"], g["Nx"])
    axes = (np.linspace(g["MinDepthKm"], g["MaxDepthKm"], g["Nz"]), np.linspace(g["MinLat"], g["MaxLat"], g["Ny"]),
            np.linspace(g["MinLon"], g["MaxLon"], g["Nx"]))
    return g, RegularGridInterpolator(axes, v, bounds_error=False)


g, vp = volume(os.path.join(NORCIA, "vel", "volumes", "vel_Vp.qvol"))
_, dws = volume(os.path.join(NORCIA, "vel", "volumes", "vel_DWS_P.qvol"))
rel, per_depth = [], {}
for depth in ("0000", "0004", "0008", "0012"):
    a = np.loadtxt(os.path.join(PUBLISHED, f"plane.{depth}.xyz"), comments="#")
    # Columns: x, y, lon, lat, z, vp, ..., presol (index 11). The inversion nodes are every 15 km.
    a = a[(a[:, 0] % 15 == 0) & (a[:, 1] % 15 == 0)]
    a = a[(a[:, 2] >= g["MinLon"]) & (a[:, 2] <= g["MaxLon"]) & (a[:, 3] >= g["MinLat"]) & (a[:, 3] <= g["MaxLat"]) & (a[:, 11] >= 0.5)]
    pts = np.c_[a[:, 4], a[:, 3], a[:, 2]]
    ok = dws(pts) >= 50
    r = (vp(pts[ok]) - a[ok, 5]) / a[ok, 5]
    per_depth[f"{int(depth)} km"] = {"nodes": int(ok.sum()), "median_abs": float(np.median(np.abs(r))) if ok.any() else None}
    rel.extend(r)
rel = np.array(rel)
record("published_model", "Vp of the central Italy 2016-2017 example (TomoSTAR)",
       "published 3-D Vp of Carannante et al. 2013 (SIMULPS14, other data)", "median |Vp - Vp_published| / Vp_published",
       float(np.median(np.abs(rel))), threshold=None,
       details={"nodes": int(rel.size), "rms": float(np.sqrt((rel ** 2).mean())), "p95_abs": float(np.percentile(np.abs(rel), 95)),
                "mean": float(rel.mean()), "by_depth": per_depth})
