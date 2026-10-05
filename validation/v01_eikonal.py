# Copyright 2026 Matteo Mangiagalli
# SPDX-License-Identifier: Apache-2.0
"""
Forward problem (fast marching on the sphere, ray tracing, the rows of G) against reference codes
on published models and real data:
  1. TauP (Crotwell et al. 1999) through ObsPy (Krischer et al. 2015): first-arrival P and S times
     in ak135 (Kennett et al. 1995), sources at 5, 15 and 30 km, receivers 10 to 500 km away;
  2. PyKonal (White et al. 2020) in the real 3-D model of the 2016-2017 central Italy example (the
     Vp and Vs volumes of examples/norcia2016, on the forward grid of that inversion), from the eight
     stations with most picks to every relocated earthquake: travel times, ray paths, and the time
     integrated along each ray as a row of the tomography matrix computes it. PyKonal's error is first
     order (several per cent within a few cells of the source on the grid TomoSTAR uses, as its
     authors note), so it is run on the grid refined four times, and pairs closer than 5 km are left
     out.

Inputs: TOMOSTAR_NORCIA_OUT, the out folder of examples/norcia2016/norcia.tomo.
"""
import csv
import json
import os
from collections import Counter

import numpy as np
from obspy.taup import TauPyModel
from scipy.interpolate import RegularGridInterpolator

# PyKonal 0.4.1 still uses np.infty, which NumPy 2 removed: the old name is restored for it.
if not hasattr(np, "infty"):
    np.infty = np.inf
import pykonal  # noqa: E402

from common import record, tomostar_validation, work  # noqa: E402

R_EARTH = 6371.0
NORCIA = os.environ["TOMOSTAR_NORCIA_OUT"]

# 1. ak135, TauP.
out = work("eikonal", "taup.json")
tomostar_validation("eikonal-taup", out)
model = TauPyModel("ak135")
errors = []
for t in json.load(open(out))["Times"]:
    phases = ["p", "P", "Pn", "Pg"] if t["Phase"] == "P" else ["s", "S", "Sn", "Sg"]
    arrivals = model.get_travel_times(source_depth_in_km=t["SourceDepthKm"], distance_in_degree=t["DistanceDeg"], phase_list=phases)
    ref = min(a.time for a in arrivals)
    errors.append((t["Phase"], t["SourceDepthKm"], t["DistanceKm"], t["Time"], ref))
e = np.array([abs(x[3] - x[4]) / x[4] for x in errors])
record("eikonal", "fast marching, ak135 first arrivals (P and S)", "TauP (ObsPy)", "max relative error", e.max(),
       details={"times": len(errors), "rms_relative_error": float(np.sqrt((e ** 2).mean())),
                "max_abs_error_s": float(max(abs(x[3] - x[4]) for x in errors))})


# 2. The real 3-D model of the central Italy example.
def axes3(g, f=1):
    return (np.linspace(g["MinLon"], g["MaxLon"], (g["Nx"] - 1) * f + 1),
            np.linspace(g["MinLat"], g["MaxLat"], (g["Ny"] - 1) * f + 1),
            np.linspace(g["MinDepthKm"], g["MaxDepthKm"], (g["Nz"] - 1) * f + 1))


def xyz(lo, la, de):
    r = R_EARTH - de
    return np.stack([r * np.cos(np.radians(la)) * np.cos(np.radians(lo)), r * np.cos(np.radians(la)) * np.sin(np.radians(lo)), r * np.sin(np.radians(la))])


def pykonal_solver(g, velocity, src, f):
    """PyKonal on the grid refined f times; velocity is an array [depth, lat, lon] on that grid."""
    lon, lat, dep = axes3(g, f)
    solver = pykonal.solver.PointSourceSolver(coord_sys="spherical")
    # PyKonal's spherical axes (rho, colatitude, longitude) all increase; ours are longitude,
    # latitude, depth: radius and colatitude run the other way, so those two axes are reversed.
    solver.velocity.min_coords = R_EARTH - g["MaxDepthKm"], np.radians(90 - g["MaxLat"]), np.radians(g["MinLon"])
    solver.velocity.node_intervals = dep[1] - dep[0], np.radians(lat[1] - lat[0]), np.radians(lon[1] - lon[0])
    solver.velocity.npts = len(dep), len(lat), len(lon)
    solver.velocity.values = np.ascontiguousarray(velocity[::-1, ::-1, :])
    solver.src_loc = np.array([R_EARTH - src["DepthKm"], np.radians(90 - src["Lat"]), np.radians(src["Lon"])])
    solver.solve()
    return solver


# The eight stations with most P picks inside the grid.
counts = Counter()
with open(os.path.join(NORCIA, "vel", "residuals.csv")) as f:
    for r in csv.DictReader(f):
        if r["phase"] == "P" and r["rejected"] == "0":
            counts[r["station"]] += 1
stations = [s for s, _ in counts.most_common(8)]
folder = work("eikonal", "norcia")
if not os.path.exists(os.path.join(folder, "times.json")):
    tomostar_validation("tt-real", os.path.join(NORCIA, "vel", "volumes", "vel_Vp.qvol"), os.path.join(NORCIA, "vel", "volumes", "vel_Vs.qvol"),
                        os.path.join(NORCIA, "vel", "stations.csv"), os.path.join(NORCIA, "vel", "events_relocated.csv"), ",".join(stations), 2, folder)
case = json.load(open(os.path.join(folder, "times.json")))
g = case["Grid"]
nx, ny, nz = g["Nx"], g["Ny"], g["Nz"]
F = 4
lon, lat, dep = axes3(g)
lon_f, lat_f, dep_f = axes3(g, F)
DD, LA, LO = np.meshgrid(dep_f, lat_f, lon_f, indexing="ij")
points = np.stack([DD.ravel(), LA.ravel(), LO.ravel()], axis=1)
del DD, LA, LO
fine = {}
for phase, name in (("P", "vp.f64"), ("S", "vs.f64")):
    v = np.fromfile(os.path.join(folder, name)).reshape(nz, ny, nx)
    # The same model on the refined grid: trilinear interpolation of the node values.
    fine[phase] = RegularGridInterpolator((dep, lat, lon), v)(points).reshape(len(dep_f), len(lat_f), len(lon_f))
del points
rows = {}
for r in case["Rows"]:
    rows.setdefault((r["Station"], r["Phase"]), []).append(r)
rel_t, deviation, kernel = [], [], []
for st in case["Stations"]:
    for phase in ("P", "S"):
        # Each station and phase is saved as it is done, so that an interrupted run resumes.
        part = os.path.join(folder, f"pykonal_{st['Id']}_{phase}.json")
        if not os.path.exists(part):
            solver = pykonal_solver(g, fine[phase], st, F)
            here = [r for r in rows[(st["Id"], phase)]
                    if np.linalg.norm(xyz(r["Lon"], r["Lat"], r["DepthKm"]) - xyz(st["Lon"], st["Lat"], st["DepthKm"])) >= 5]
            pts = np.array([[R_EARTH - r["DepthKm"], np.radians(90 - r["Lat"]), np.radians(r["Lon"])] for r in here])
            ref = solver.tt.resample(pts)
            result = {"relative": [abs(r["Time"] - t) / t for r, t in zip(here, ref)], "deviation": [], "kernel": []}
            for ray in [x for x in case["Rays"] if x["Station"] == st["Id"]] if phase == "P" else []:
                if np.linalg.norm(xyz(ray["Lon"], ray["Lat"], ray["Depth"]) - xyz(st["Lon"], st["Lat"], st["DepthKm"])) < 5:
                    continue
                end = np.array([R_EARTH - ray["Depth"], np.radians(90 - ray["Lat"]), np.radians(ray["Lon"])])
                pk = solver.tt.trace_ray(end)
                pk_xyz = np.stack([pk[:, 0] * np.sin(pk[:, 1]) * np.cos(pk[:, 2]), pk[:, 0] * np.sin(pk[:, 1]) * np.sin(pk[:, 2]), pk[:, 0] * np.cos(pk[:, 1])], axis=1)
                ts_xyz = xyz(np.array(ray["RayLon"]), np.array(ray["RayLat"]), np.array(ray["RayDepth"])).T
                dense = np.concatenate([np.linspace(pk_xyz[i], pk_xyz[i + 1], 10, endpoint=False) for i in range(len(pk_xyz) - 1)] + [pk_xyz[-1:]])
                result["deviation"].append(float(np.array([np.min(np.linalg.norm(dense - p, axis=1)) for p in ts_xyz]).max() / ray["LengthKm"]))
                result["kernel"].append(abs(ray["KernelTime"] - ray["FieldTime"]) / ray["FieldTime"])
            del solver
            json.dump(result, open(part, "w"))
            print(f"PyKonal {st['Id']} {phase}: {len(result['relative'])} times, {len(result['deviation'])} rays", flush=True)
        result = json.load(open(part))
        rel_t += result["relative"]
        deviation += result["deviation"]
        kernel += result["kernel"]
rel_t, deviation, kernel = np.array(rel_t), np.array(deviation), np.array(kernel)
details = {"stations": stations, "pairs": int(rel_t.size), "p99_relative_error": float(np.percentile(rel_t, 99)),
           "rms_relative_error": float(np.sqrt((rel_t ** 2).mean())), "grid": g}
record("eikonal", "fast marching, real 3-D Vp and Vs model (central Italy 2016-2017), station to hypocentre, beyond 5 km",
       "PyKonal (grid refined 4 times)", "max relative error", rel_t.max(), details=details)
record("rays", "ray paths traced back along the travel-time gradient, real 3-D model", "PyKonal", "max deviation / ray length",
       deviation.max(), details={"rays": int(deviation.size), "median": float(np.median(deviation))})
record("rays", "travel time integrated along the ray (row of G), real 3-D model", "eikonal travel time", "max relative error",
       kernel.max(), details={"rays": int(kernel.size), "median": float(np.median(kernel))})
