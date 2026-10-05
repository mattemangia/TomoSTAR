# Copyright 2026 Matteo Mangiagalli
# SPDX-License-Identifier: Apache-2.0
"""
Absolute earthquake location on real, published data: the 2016-2017 central Italy sequence.

Reference case (peer reviewed): Chiaraluce et al. (2022), Scientific Data 9, 710, catalogue CAT1,
located with NonLinLoc (Lomax et al., 2000) from the analyst P and S picks of the INGV bulletin in
the 1-D gradient model of Carannante et al. (2013), with station corrections.

Two comparisons:
  1. code against code on identical input: the picks of the INGV bulletin (ISIDe Working Group,
     2007) for the M >= 2.5 events of the sequence, the Carannante et al. (2013) model (the nodes of
     the authors' data set, Carannante et al., 2025), no station corrections, weighted least squares
     (L2) in both programs, no outlier rejection. TomoSTAR (grid search + Levenberg-Marquardt) is
     compared with NonLinLoc (oct-tree search, GAU_ANALYTIC likelihood, maximum-likelihood point).
     Metric: the hypocentre difference relative to the median distance from the hypocentre to the
     stations of the event (the length scale of the location problem); the 95th percentile of it
     over all events must be below 5 %.
  2. TomoSTAR against the published CAT1 hypocentres (whose picks include those of temporary
     stations and whose locations include station corrections, neither of which is public): the
     same metric, reported for information together with the same comparison for NonLinLoc run on
     TomoSTAR's input, which separates the effect of the input from that of the code.

Inputs (environment variables):
  TOMOSTAR_NORCIA_OUT      the out folder of examples/norcia2016/norcia.tomo (reloc/ holds the
                           INGV bulletin picks as TomoSTAR read them)
  TOMOSTAR_CHIARALUCE2022  the folder with the files CAT1 and CAT2 of Chiaraluce et al. (2022),
                           https://doi.org/10.5285/5afccfe5-142e-4e93-a6cc-55216fa1db06
  NLL_BIN                  the folder with the NonLinLoc programs Vel2Grid, Grid2Time and NLLoc
"""
import csv
import glob
import json
import math
import os
import subprocess
from datetime import datetime, timezone

import numpy as np

from common import record, tomostar, work

NORCIA = os.environ["TOMOSTAR_NORCIA_OUT"]
CATS = os.environ["TOMOSTAR_CHIARALUCE2022"]
NLL = os.environ["NLL_BIN"]
BOX = (12.75, 13.70, 42.35, 43.20)  # the grid of the Norcia example: lon min, lon max, lat min, lat max
LAT0, LON0 = 42.775, 13.225
R = 6371.0
# Carannante et al. (2013), 1-D P and S nodes of the authors' data set (linear between nodes), as in
# TomoSTAR's model library (id Carannante2013): depth km, Vp, Vs.
NODES = [(0, 5.63, 2.80), (4, 6.22, 3.34), (8, 6.23, 3.386), (12, 6.24, 3.39), (20, 6.26, 3.37),
         (30, 6.62, 3.64), (80, 7.92, 4.33)]


def read_cat(name):
    """A catalogue of Chiaraluce et al. (2022): INGV id -> (time, lat, lon, depth, errh, errv)."""
    rows = {}
    with open(os.path.join(CATS, name)) as f:
        next(f)
        for line in f:
            p = line.split()
            t = datetime.strptime(p[1] + " " + p[2], "%Y-%m-%d %H:%M:%S.%f").replace(tzinfo=timezone.utc)
            rows[p[0]] = (t, float(p[3]), float(p[4]), float(p[5]), float(p[6]), float(p[7]))
    return rows


def distance_km(lat1, lon1, lat2, lon2):
    """Great-circle distance on the sphere of radius R."""
    a, b = math.radians(lat1), math.radians(lat2)
    d = math.sin((b - a) / 2) ** 2 + math.cos(a) * math.cos(b) * math.sin(math.radians(lon2 - lon1) / 2) ** 2
    return 2 * R * math.asin(min(1.0, math.sqrt(d)))


def hypocentral_difference(a, b):
    """3-D distance between two hypocentres (lat, lon, depth), km."""
    return math.hypot(distance_km(a[0], a[1], b[0], b[1]), a[2] - b[2])


# ---- The common input: stations in the grid, one pick per station and phase ----------------------
cat1 = read_cat("CAT1")
stations = {}
with open(os.path.join(NORCIA, "reloc", "stations.csv")) as f:
    for s in csv.DictReader(f):
        lon, lat = float(s["lon"]), float(s["lat"])
        if BOX[0] <= lon <= BOX[1] and BOX[2] <= lat <= BOX[3]:
            stations[s["id"].upper()] = (lat, lon, float(s["elevation_m"]))
events = {}
with open(os.path.join(NORCIA, "reloc", "events.csv")) as f:
    for e in csv.DictReader(f):
        if e["id"] in cat1:
            events[e["id"]] = e
picks = {}
with open(os.path.join(NORCIA, "reloc", "picks.csv")) as f:
    for p in csv.DictReader(f):
        sid = p["station"].upper()
        if p["event"] not in events or sid not in stations or p["disabled"] == "1" or int(p["quality"]) >= 4:
            continue
        key = (p["event"], sid, p["phase"])
        # The rule of TomoSTAR's pick selection for picks of equal origin and uncertainty: the earliest.
        if key not in picks or p["time"] < picks[key]["time"]:
            picks[key] = p
by_event = {}
for (ev, sid, ph), p in picks.items():
    by_event.setdefault(ev, []).append(p)
keep = [ev for ev, ps in by_event.items() if sum(p["phase"] == "P" for p in ps) >= 4 and len(ps) >= 6]
keep.sort()
print(f"{len(keep)} events of CAT1, {sum(len(by_event[e]) for e in keep)} picks, {len(stations)} stations inside the grid.")

folder = work("location")
with open(os.path.join(folder, "stations.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["id", "lon", "lat", "elevation_m"])
    for sid, (lat, lon, elev) in sorted(stations.items()):
        w.writerow([sid, lon, lat, elev])
with open(os.path.join(folder, "events.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["id", "time", "lon", "lat", "depth_km", "magnitude", "fixed"])
    for ev in keep:
        e = events[ev]
        w.writerow([ev, e["time"], e["lon"], e["lat"], e["depth_km"], e["magnitude"], 0])
with open(os.path.join(folder, "picks.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["event", "station", "phase", "time", "sigma_s", "quality", "origin"])
    for ev in keep:
        for p in sorted(by_event[ev], key=lambda p: p["time"]):
            w.writerow([ev, p["station"].upper(), p["phase"], p["time"], p["sigma_s"], 0, "catalog"])

# ---- TomoSTAR: travel times by fast marching on a 0.01 deg x 0.5 km grid -------------------------
grid = {"MinLon": BOX[0], "MaxLon": BOX[1], "MinLat": BOX[2], "MaxLat": BOX[3], "MinDepthKm": -2.5, "MaxDepthKm": 25.0,
        "Nx": 96, "Ny": 86, "Nz": 56}
json.dump(grid, open(os.path.join(folder, "grid.json"), "w"), indent=1)
ts_out = os.path.join(folder, "tomostar")
if not os.path.exists(os.path.join(ts_out, "locations.csv")):
    tomostar("relocate", "--stations", os.path.join(folder, "stations.csv"), "--events", os.path.join(folder, "events.csv"),
             "--picks", os.path.join(folder, "picks.csv"), "--grid", os.path.join(folder, "grid.json"), "--model", "Carannante2013",
             "--set", "Locator.OutlierMads=1000000", "--set", "Locator.UseStationCorrections=false",
             "--set", "Locator.MinPhases=4", "--set", "Tomography.ForwardRefinement=1", "--out", ts_out, "--name", "tomostar", "--quiet")
ts = {}
with open(os.path.join(ts_out, "locations.csv")) as f:
    for r in csv.DictReader(f):
        if r["status"] == "ok":
            ts[r["id"]] = (float(r["lat"]), float(r["lon"]), float(r["depth_km"]), float(r["rms_s"]))

# ---- NonLinLoc: the same picks, the same model with its gradients, travel times by finite differences
# (Podvin and Lecomte, 1991) on a 0.1 km 2-D grid, a sphere of radius 6371 km -------------------------
nll = work("location", "nll")
for sub in ("model", "time", "loc"):
    os.makedirs(os.path.join(nll, sub), exist_ok=True)
labels = {sid: f"S{i:03d}" for i, sid in enumerate(sorted(stations))}
layers = ["LAYER -3.0 %.4f 0.0 %.4f 0.0 2.70 0.0" % (NODES[0][1], NODES[0][2])]
for (z0, p0, s0), (z1, p1, s1) in zip(NODES, NODES[1:]):
    layers.append("LAYER %.1f %.4f %.6f %.4f %.6f 2.70 0.0" % (z0, p0, (p1 - p0) / (z1 - z0), s0, (s1 - s0) / (z1 - z0)))
sources = [f"GTSRCE {labels[sid]} LATLON {lat:.6f} {lon:.6f} 0.0 {elev / 1000:.4f}" for sid, (lat, lon, elev) in sorted(stations.items())]
common = f"""CONTROL 1 54321
TRANS AZIMUTHAL_EQUIDIST Sphere {LAT0} {LON0} 0.0
VGOUT ./model/layer
VGTYPE P
VGTYPE S
VGGRID 2 1401 301 0.0 0.0 -3.0 0.1 0.1 0.1 SLOW_LEN
{chr(10).join(layers)}
GTMODE GRID2D ANGLES_NO
{chr(10).join(sources)}
GT_PLFD 1.0e-3 0
"""
obs = os.path.join(nll, "events.obs")
with open(obs, "w") as f:
    for ev in keep:
        f.write(f"PUBLIC_ID {ev}\n")
        for p in sorted(by_event[ev], key=lambda p: p["time"]):
            t = datetime.strptime(p["time"].rstrip("Z")[:26], "%Y-%m-%dT%H:%M:%S.%f")
            sec = t.second + t.microsecond / 1e6
            f.write(f"{labels[p['station'].upper()]:6s} ?    ?    ? {p['phase']:6s} ? {t:%Y%m%d} {t:%H%M} {sec:7.4f} GAU "
                    f"{float(p['sigma_s']):9.2e} -1.00e+00 -1.00e+00 -1.00e+00 1.0\n")
        f.write("\n")
loc = f"""LOCSIG TomoSTAR validation
LOCFILES ./events.obs NLLOC_OBS ./time/layer ./loc/norcia
LOCHYPOUT SAVE_NLLOC_SUM
LOCSEARCH OCT 20 20 6 0.005 30000 1000 0 1
LOCGRID 201 201 57 -50.0 -50.0 -3.0 0.5 0.5 0.5 PROB_DENSITY SAVE
LOCMETH GAU_ANALYTIC 9999.0 4 -1 -1 -1 -1 -1.0 1
LOCGAU 0.001 0.0
LOCPHASEID P P
LOCPHASEID S S
LOCQUAL2ERR 0.1 0.5 1.0 2.0 99999.9
LOCANGLES ANGLES_NO 5
"""
summary = os.path.join(nll, "loc", "norcia.sum.grid0.loc.hyp")
if not os.path.exists(summary):
    open(os.path.join(nll, "vel.in"), "w").write(common + "GTFILES ./model/layer ./time/layer P\n")
    subprocess.run([os.path.join(NLL, "Vel2Grid"), "vel.in"], cwd=nll, check=True, stdout=subprocess.DEVNULL)
    for phase in ("P", "S"):
        open(os.path.join(nll, f"time_{phase}.in"), "w").write(common + f"GTFILES ./model/layer ./time/layer {phase}\n")
        subprocess.run([os.path.join(NLL, "Grid2Time"), f"time_{phase}.in"], cwd=nll, check=True, stdout=subprocess.DEVNULL)
    open(os.path.join(nll, "loc.in"), "w").write(common + "GTFILES ./model/layer ./time/layer P\n" + loc)
    subprocess.run([os.path.join(NLL, "NLLoc"), "loc.in"], cwd=nll, check=True, stdout=subprocess.DEVNULL)
nl = {}
current = {}
for line in open(summary):
    p = line.split()
    if not p:
        continue
    if p[0] == "NLLOC":
        current = {"ok": "LOCATED" in line}
    elif p[0] == "PUBLIC_ID":
        current["id"] = p[1]
    elif p[0] == "GEOGRAPHIC":
        current["lat"], current["lon"], current["depth"] = float(p[9]), float(p[11]), float(p[13])
    elif p[0] == "QUALITY":
        current["rms"] = float(p[8])
    elif p[0] == "END_NLLOC" and current.get("ok") and "id" in current:
        nl[current["id"]] = (current["lat"], current["lon"], current["depth"], current["rms"])


def compare(reference):
    """Per-event differences of TomoSTAR from a reference: absolute (km) and relative to the scale."""
    rows = []
    for ev, a in ts.items():
        if ev not in reference:
            continue
        b = reference[ev]
        d = hypocentral_difference(a, b[:3])
        scale = float(np.median([math.hypot(distance_km(b[0], b[1], *stations[p["station"].upper()][:2]),
                                            b[2] + stations[p["station"].upper()][2] / 1000) for p in by_event[ev]]))
        rows.append((ev, d, d / scale, distance_km(a[0], a[1], b[0], b[1]), a[2] - b[2]))
    return rows


def stats(rows):
    d = np.array([r[1] for r in rows])
    rel = np.array([r[2] for r in rows])
    h = np.array([r[3] for r in rows])
    z = np.array([r[4] for r in rows])
    return {"events": len(rows), "median_km": float(np.median(d)), "p95_km": float(np.percentile(d, 95)),
            "median_horizontal_km": float(np.median(h)), "median_abs_depth_km": float(np.median(np.abs(z))),
            "mean_depth_difference_km": float(z.mean()), "median_relative": float(np.median(rel)),
            "p95_relative": float(np.percentile(rel, 95))}


# 1. Code against code.
rows = compare(nl)
s = stats(rows)
rms_ts = np.array([ts[r[0]][3] for r in rows])
rms_nl = np.array([nl[r[0]][3] for r in rows])
s["rms_tomostar_median_s"] = float(np.median(rms_ts))
s["rms_nonlinloc_median_s"] = float(np.median(rms_nl))
s["located_tomostar"] = len(ts)
s["located_nonlinloc"] = len(nl)
s["input_events"] = len(keep)
record("location", "absolute location (grid search + Levenberg-Marquardt), INGV picks, Carannante et al. 2013 model",
       "NonLinLoc (Lomax et al. 2000), same picks and model", "95th percentile of |hypocentre difference| / median hypocentral distance",
       s["p95_relative"], details=s)

# 2. Against the published catalogue CAT1 (NonLinLoc with station corrections and more picks).
cat1_pos = {ev: v[1:4] for ev, v in cat1.items()}
rows = compare(cat1_pos)
s = stats(rows)
# The same comparison for NonLinLoc run on the same input as TomoSTAR: if it differs from CAT1 as much,
# the difference comes from the input (station corrections, the temporary stations' picks), not the code.
ts_saved = ts
ts = {ev: v for ev, v in nl.items()}
s["nonlinloc_same_input_vs_cat1"] = {k: v for k, v in stats(compare(cat1_pos)).items() if k != "events"}
ts = ts_saved
errh = np.array([cat1[r[0]][4] for r in rows])
errv = np.array([cat1[r[0]][5] for r in rows])
s["published_median_errh_km"] = float(np.median(errh))
s["published_median_errv_km"] = float(np.median(errv))
record("location", "absolute location, INGV picks, Carannante et al. 2013 model, no station corrections",
       "published catalogue CAT1 of Chiaraluce et al. 2022 (NonLinLoc with station corrections)",
       "95th percentile of |hypocentre difference| / median hypocentral distance", s["p95_relative"], threshold=None, details=s)
