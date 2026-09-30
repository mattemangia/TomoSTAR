# Copyright 2026 Matteo Mangiagalli
# SPDX-License-Identifier: Apache-2.0
"""
Eikonal solver (fast marching on the sphere) against:
  1. the closed-form travel time in a medium whose velocity increases linearly with depth;
  2. TauP (Crotwell et al. 1999) through ObsPy (Krischer et al. 2015): first-arrival P and S times
     in ak135 (Kennett et al. 1995), sources at 5, 15 and 30 km, receivers 10 to 500 km away;
  3. the exact solution in a uniform medium, for TomoSTAR and PyKonal (White et al. 2020) alike;
  4. PyKonal: the travel-time fields of three sources in a 3-D spherical model (a regional 1-D model
     with an 8 % checkerboard); PyKonal is run on the grid refined four times, since its own error
     near a source is several per cent on the grid TomoSTAR uses (case 3).
"""
import json

import numpy as np
from obspy.taup import TauPyModel

# PyKonal 0.4.1 still uses np.infty, which NumPy 2 removed: the old name is restored for it.
if not hasattr(np, "infty"):
    np.infty = np.inf
import pykonal  # noqa: E402

from common import record, tomostar_validation, work

R_EARTH = 6371.0

# 1. Linear gradient, closed form.
out = work("eikonal", "gradient.json")
tomostar_validation("eikonal-gradient", out)
rows = np.array(json.load(open(out))["Rows"])
rel = np.abs(rows[:, 4] - rows[:, 5]) / rows[:, 5]
record("eikonal", "fast marching, v = 4 + 0.05 z km/s", "closed-form solution", "max relative error", rel.max(),
       details={"nodes": int(len(rel)), "rms_relative_error": float(np.sqrt((rel ** 2).mean()))})

# 2. ak135, TauP.
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

# 3. Near-source accuracy of both solvers: a uniform 6 km/s on the grid of the 3-D case, where the
# exact time is the distance divided by 6 km/s. PyKonal is run with its defaults.
R3 = dict(MinLon=12.8, MaxLon=13.7, MinLat=42.3, MaxLat=43.2, MinDepthKm=-2.0, MaxDepthKm=30.0, Nx=91, Ny=91, Nz=65)


def axes3(g, f=1):
    return (np.linspace(g["MinLon"], g["MaxLon"], (g["Nx"] - 1) * f + 1),
            np.linspace(g["MinLat"], g["MaxLat"], (g["Ny"] - 1) * f + 1),
            np.linspace(g["MinDepthKm"], g["MaxDepthKm"], (g["Nz"] - 1) * f + 1))


def xyz(lo, la, de):
    r = R_EARTH - de
    return np.stack([r * np.cos(np.radians(la)) * np.cos(np.radians(lo)), r * np.cos(np.radians(la)) * np.sin(np.radians(lo)), r * np.sin(np.radians(la))])


def distances(g, src, f=1):
    lon, lat, dep = axes3(g, f)
    LON, LAT, DEP = np.meshgrid(lon, lat, dep, indexing="ij")
    return np.linalg.norm(xyz(LON, LAT, DEP) - xyz(src["Lon"], src["Lat"], src["Depth"])[:, None, None, None], axis=0).transpose(2, 1, 0)


def pykonal_times(g, velocity, src, f=1):
    """PyKonal's travel times on the grid refined f times, as an array [depth, lat, lon]."""
    lon, lat, dep = axes3(g, f)
    solver = pykonal.solver.PointSourceSolver(coord_sys="spherical")
    # PyKonal's spherical axes (rho, colatitude, longitude) all increase; ours are longitude,
    # latitude, depth: radius and colatitude run the other way, so those two axes are reversed.
    solver.velocity.min_coords = R_EARTH - g["MaxDepthKm"], np.radians(90 - g["MaxLat"]), np.radians(g["MinLon"])
    solver.velocity.node_intervals = dep[1] - dep[0], np.radians(lat[1] - lat[0]), np.radians(lon[1] - lon[0])
    solver.velocity.npts = len(dep), len(lat), len(lon)
    solver.velocity.values = np.ascontiguousarray(velocity[::-1, ::-1, :])
    solver.src_loc = np.array([R_EARTH - src["Depth"], np.radians(90 - src["Lat"]), np.radians(src["Lon"])])
    solver.solve()
    return solver, solver.tt.values[::-1, ::-1, :]


src0 = {"Lon": 13.05, "Lat": 42.55, "Depth": 0.0}
out = work("eikonal", "homogeneous.f64")
tomostar_validation("eikonal-homogeneous", out)
d = distances(R3, src0)
exact = d / 6.0
far = d > 3
t_ts = np.fromfile(out).reshape(R3["Nz"], R3["Ny"], R3["Nx"])
_, t_pk = pykonal_times(R3, np.full((R3["Nz"], R3["Ny"], R3["Nx"]), 6.0), src0)
rel_ts = np.abs(t_ts[far] - exact[far]) / exact[far]
rel_pk = np.abs(t_pk[far] - exact[far]) / exact[far]
record("eikonal", "fast marching, uniform medium, nodes beyond 3 km", "exact solution", "max relative error", rel_ts.max(),
       details={"nodes": int(far.sum()), "pykonal_max_relative_error": float(rel_pk.max()),
                "pykonal_median_relative_error": float(np.median(rel_pk)), "tomostar_median_relative_error": float(np.median(rel_ts))})

# 4. 3-D model: PyKonal on the grid refined four times (its error is first order and halves with
# every refinement), sampled at TomoSTAR's nodes, as the reference.
folder = work("eikonal", "3d")
tomostar_validation("eikonal-3d", folder)
case = json.load(open(f"{folder}/case.json"))
g = case["Grid"]
nx, ny, nz = g["Nx"], g["Ny"], g["Nz"]
vp = np.fromfile(f"{folder}/vp.f64").reshape(nz, ny, nx)
F = 4
lon_f, lat_f, dep_f = axes3(g, F)
# The velocity on the refined grid: trilinear interpolation of the same node values.
from scipy.interpolate import RegularGridInterpolator  # noqa: E402
lon, lat, dep = axes3(g)
interp = RegularGridInterpolator((dep, lat, lon), vp)
DD, LA, LO = np.meshgrid(dep_f, lat_f, lon_f, indexing="ij")
vp_f = interp(np.stack([DD.ravel(), LA.ravel(), LO.ravel()], axis=1)).reshape(DD.shape)
del DD, LA, LO
all_rel = []
for q, src in enumerate(case["Sources"]):
    solver, t_ref = pykonal_times(g, vp_f, src, F)
    t_ref = t_ref[::F, ::F, ::F]
    t_ts = np.fromfile(f"{folder}/t{q}.f64").reshape(nz, ny, nx)
    m = distances(g, src) > 5
    all_rel.append(np.abs(t_ts[m] - t_ref[m]) / t_ref[m])
    # Rays: the path of TomoSTAR against PyKonal's (each traced back along the gradient of its own field).
    for ray in [r for r in case["Rays"] if r["Source"] == q]:
        end = np.array([R_EARTH - ray["Depth"], np.radians(90 - ray["Lat"]), np.radians(ray["Lon"])])
        pk = solver.tt.trace_ray(end)
        pk_xyz = np.stack([pk[:, 0] * np.sin(pk[:, 1]) * np.cos(pk[:, 2]), pk[:, 0] * np.sin(pk[:, 1]) * np.sin(pk[:, 2]), pk[:, 0] * np.cos(pk[:, 1])], axis=1)
        ts_xyz = xyz(np.array(ray["RayLon"]), np.array(ray["RayLat"]), np.array(ray["RayDepth"])).T
        dense = np.concatenate([np.linspace(pk_xyz[i], pk_xyz[i + 1], 10, endpoint=False) for i in range(len(pk_xyz) - 1)] + [pk_xyz[-1:]])
        dev = np.array([np.min(np.linalg.norm(dense - p, axis=1)) for p in ts_xyz])
        ray["deviation"] = float(dev.max() / ray["LengthKm"])
        ray["kernel_vs_field"] = abs(ray["KernelTime"] - ray["FieldTime"]) / ray["FieldTime"]
    del solver
rel = np.concatenate(all_rel)
record("eikonal", "fast marching, 3-D spherical model, nodes beyond 5 km", "PyKonal (grid refined 4 times)", "max relative error", rel.max(),
       details={"nodes": int(rel.size), "rms_relative_error": float(np.sqrt((rel ** 2).mean())), "p99_relative_error": float(np.percentile(rel, 99))})
rays = case["Rays"]
record("rays", "ray paths traced back along the travel-time gradient", "PyKonal", "max deviation / ray length",
       max(r["deviation"] for r in rays), details={"rays": len(rays)})
record("rays", "travel time integrated along the ray (row of G)", "eikonal travel time", "max relative error",
       max(r["kernel_vs_field"] for r in rays), details={"rays": len(rays)})
