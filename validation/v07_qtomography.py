# Copyright 2026 Matteo Mangiagalli
# SPDX-License-Identifier: Apache-2.0
"""
Qp tomography of the 2016-2017 central Italy example (real t* of the M >= 3.5 earthquakes, real
rays in the 3-D velocity model): the forward operator recomputed independently.

TomoSTAR models t* = sum_n L_n s_n q_n + c_station, with L_n the length of the ray on the trilinear
basis function of node n, s_n the slowness and q_n = 1/Q at the node, and c the station term. This
script recomputes the sum from the stored rays (rays.qray, polylines of 48 points), the final Qp
and Vp volumes and the station terms, with its own trilinear weights, and compares it with the t*
TomoSTAR predicts (observed t* minus the final residual). It checks the matrix rows of the Q
inversion and the consistency of the written model, residuals and station terms. The linear solver
of each Gauss-Newton step (unknowns ln(q/q0)) is the LSQR checked in v02, and the t* measurements are checked against AttenTIon in v06.
Metric: |t*_TomoSTAR - t*_recomputed| / t*_TomoSTAR; its 95th percentile must be below 5 %.

Inputs: TOMOSTAR_NORCIA_OUT, the out folder of examples/norcia2016/norcia.tomo.
"""
import csv
import json
import os
import struct

import numpy as np

from common import record

NORCIA = os.environ["TOMOSTAR_NORCIA_OUT"]
R_EARTH = 6371.0


def read_qvol(path):
    """A .qvol volume: (grid, values[depth, lat, lon])."""
    with open(path, "rb") as f:
        b = f.read()
    n = struct.unpack("<i", b[8:12])[0]
    header = json.loads(b[12:12 + n])
    start = ((12 + n + 63) // 64) * 64
    g = header["Grid"]
    v = np.frombuffer(b[start:start + 4 * g["Nx"] * g["Ny"] * g["Nz"]], "<f4").astype(float)
    return g, v.reshape(g["Nz"], g["Ny"], g["Nx"])


def read_rays(path):
    rays = []
    with open(path, "rb") as f:
        assert f.read(4) == b"QRAY"
        (n,) = struct.unpack("<i", f.read(4))
        for _ in range(n):
            arrival, event, station = struct.unpack("<3i", f.read(12))
            f.read(1)
            (m,) = struct.unpack("<i", f.read(4))
            pts = np.frombuffer(f.read(12 * m), "<f4").reshape(3, m)
            rays.append((arrival, pts[0].astype(float), pts[1].astype(float), pts[2].astype(float)))
    return rays


def xyz(lon, lat, dep):
    r = R_EARTH - dep
    return np.stack([r * np.cos(np.radians(lat)) * np.cos(np.radians(lon)), r * np.cos(np.radians(lat)) * np.sin(np.radians(lon)), r * np.sin(np.radians(lat))], axis=-1)


g, qp = read_qvol(os.path.join(NORCIA, "q", "volumes", "q_Qp.qvol"))
gv, vp = read_qvol(os.path.join(NORCIA, "vel", "volumes", "vel_Vp.qvol"))
assert g == gv, "the Q and velocity grids must be the same"
sq = (1.0 / vp) * (1.0 / qp)  # s_n q_n at the nodes
dx = (g["MaxLon"] - g["MinLon"]) / (g["Nx"] - 1)
dy = (g["MaxLat"] - g["MinLat"]) / (g["Ny"] - 1)
dz = (g["MaxDepthKm"] - g["MinDepthKm"]) / (g["Nz"] - 1)


def path_integral(lon, lat, dep, sub=20):
    """sum_n L_n (s q)_n along a polyline: each segment cut in `sub` pieces, trilinear weights at their midpoints."""
    total = 0.0
    p = xyz(lon, lat, dep)
    for k in range(len(lon) - 1):
        seg = np.linalg.norm(p[k + 1] - p[k]) / sub
        f = (np.arange(sub) + 0.5) / sub
        x = (lon[k] + f * (lon[k + 1] - lon[k]) - g["MinLon"]) / dx
        y = (lat[k] + f * (lat[k + 1] - lat[k]) - g["MinLat"]) / dy
        z = (dep[k] + f * (dep[k + 1] - dep[k]) - g["MinDepthKm"]) / dz
        i0 = np.clip(np.floor(x).astype(int), 0, g["Nx"] - 2)
        j0 = np.clip(np.floor(y).astype(int), 0, g["Ny"] - 2)
        k0 = np.clip(np.floor(z).astype(int), 0, g["Nz"] - 2)
        fx, fy, fz = np.clip(x - i0, 0, 1), np.clip(y - j0, 0, 1), np.clip(z - k0, 0, 1)
        val = np.zeros(sub)
        for a in (0, 1):
            for b in (0, 1):
                for c in (0, 1):
                    w = (fx if a else 1 - fx) * (fy if b else 1 - fy) * (fz if c else 1 - fz)
                    val += w * sq[k0 + c, j0 + b, i0 + a]
        total += seg * val.sum()
    return total


with open(os.path.join(NORCIA, "q", "residuals.csv")) as f:
    residuals = list(csv.DictReader(f))
terms = {}
with open(os.path.join(NORCIA, "q", "station_terms.csv")) as f:
    for r in csv.DictReader(f):
        terms[r["station"]] = float(r["tstar_term_s"])
rel, absd, n = [], [], 0
for arrival, lon, lat, dep in read_rays(os.path.join(NORCIA, "q", "rays.qray")):
    r = residuals[arrival]
    if r["rejected"] != "0":
        continue
    predicted = float(r["observed_s"]) - float(r["final_residual_s"])
    recomputed = path_integral(lon, lat, dep) + terms.get(r["station"], 0.0)
    absd.append(abs(predicted - recomputed))
    rel.append(abs(predicted - recomputed) / abs(predicted))
    n += 1
rel, absd = np.array(rel), np.array(absd)
record("qtomography", "Qp tomography forward operator and model (t* = sum L s q + station term), real rays and t*",
       "independent path integration of the written model along the stored rays", "95th percentile of relative t* difference",
       float(np.percentile(rel, 95)),
       details={"t_star": n, "median_relative": float(np.median(rel)), "median_abs_s": float(np.median(absd)),
                "p95_abs_s": float(np.percentile(absd, 95))})
