# Copyright 2026 Matteo Mangiagalli
# SPDX-License-Identifier: Apache-2.0
"""
Double-difference relocation on real, published data: the 2016-2017 central Italy sequence.

Reference case (peer reviewed): Michele et al. (2020), JGR Solid Earth 125, e2019JB018440, catalogue
CAT2 of Chiaraluce et al. (2022), relocated with HypoDD (Waldhauser and Ellsworth, 2000; Waldhauser,
2001) from the INGV bulletin picks and waveform cross-correlation delays in the 1-D model of
Carannante et al. (2013).

Two comparisons:
  1. code against code on identical input: HypoDD 2.1b and TomoSTAR relocate the events of the
     absolute-location test (v04_location.py) from the same starting hypocentres (TomoSTAR's absolute
     locations), with the catalogue differential times of the same INGV picks, the same pair
     selection (ph2dt parameters MAXSEP 10 km, MAXNGH 10, MINLNK 8, MINOBS 8, MAXOBS 50) and the same
     iteration sets (2, 4 and 6 iterations: no cutoff; 6 MADs and 10 km; 5 MADs and 6 km; P and S
     weights 1 and 0.5; TomoSTAR keeps the absolute times at a token weight of 0.01). HypoDD uses a layered
     version of the Carannante et al. (2013) gradient model (1 km layers, as Michele et al. 2020 did),
     TomoSTAR the gradient model itself. A double-difference solution fixes the relative positions of
     the events of a cluster better than its centroid, so the mean offset of each cluster is removed
     before comparing. Metric: the hypocentre difference relative to the median hypocentral distance
     of the event (as in v04); its 95th percentile must be below 5 %.
  2. TomoSTAR against the published CAT2 (which adds cross-correlation delays), the same metric,
     reported for information (the input differs), with the same comparison for HypoDD run on
     TomoSTAR's input.

Inputs: those of v04_location.py (which must have run), and HYPODD_BIN, the folder with the HypoDD
programs ph2dt and hypoDD.
"""
import csv
import json
import math
import os
import subprocess
from datetime import datetime, timezone

import numpy as np

from common import record, tomostar, work

CATS = os.environ["TOMOSTAR_CHIARALUCE2022"]
HYPODD = os.environ["HYPODD_BIN"]
R = 6371.0
NODES = [(0, 5.63, 2.80), (4, 6.22, 3.34), (8, 6.23, 3.386), (12, 6.24, 3.39), (20, 6.26, 3.37),
         (30, 6.62, 3.64), (80, 7.92, 4.33)]
base = work("location")
folder = work("dd")


def model(z):
    """Vp and Vs of the Carannante et al. (2013) nodes at depth z (constant above the first node)."""
    z = max(z, NODES[0][0])
    for (z0, p0, s0), (z1, p1, s1) in zip(NODES, NODES[1:]):
        if z <= z1:
            f = (z - z0) / (z1 - z0)
            return p0 + f * (p1 - p0), s0 + f * (s1 - s0)
    return NODES[-1][1], NODES[-1][2]


def parse_time(text):
    return datetime.strptime(text.rstrip("Z")[:26], "%Y-%m-%dT%H:%M:%S.%f").replace(tzinfo=timezone.utc)


def distance_km(lat1, lon1, lat2, lon2):
    a, b = math.radians(lat1), math.radians(lat2)
    d = math.sin((b - a) / 2) ** 2 + math.cos(a) * math.cos(b) * math.sin(math.radians(lon2 - lon1) / 2) ** 2
    return 2 * R * math.asin(min(1.0, math.sqrt(d)))


def local_xyz(lat, lon, depth, lat0, lon0):
    """East, north, down in km around (lat0, lon0) (small region: local tangent plane on the sphere)."""
    return np.array([math.radians(lon - lon0) * R * math.cos(math.radians(lat0)), math.radians(lat - lat0) * R, depth])


# ---- Common input: the stations, picks and absolute locations of v04 ------------------------------
stations = {}
with open(os.path.join(base, "stations.csv")) as f:
    for s in csv.DictReader(f):
        stations[s["id"]] = (float(s["lat"]), float(s["lon"]), float(s["elevation_m"]))
start = {}
with open(os.path.join(base, "tomostar", "locations.csv")) as f:
    for r in csv.DictReader(f):
        if r["status"] == "ok":
            start[r["id"]] = r
picks = {}
with open(os.path.join(base, "picks.csv")) as f:
    for p in csv.DictReader(f):
        if p["event"] in start:
            picks.setdefault(p["event"], []).append(p)
ids = sorted(start)
number = {ev: i + 1 for i, ev in enumerate(ids)}  # HypoDD needs integer event ids
print(f"{len(ids)} events from the absolute locations of v04.")

# ---- TomoSTAR ----------------------------------------------------------------------------------------
# HypoDD's iteration sets (below) in TomoSTAR's terms: differential P and S weights 1 and 0.5, the
# same residual (MADs) and separation (km) cutoffs, and the absolute times kept only at a token weight
# (HypoDD inverts differential times alone and fixes the cluster centroid instead).
config = {"Tomography": {"DoubleDifference": {"Sets": [
    {"Iterations": 2, "AbsoluteWeight": 0.01, "DifferentialWeightP": 1, "DifferentialWeightS": 0.5, "ResidualCutoff": 0, "MaxSeparationKm": 0},
    {"Iterations": 4, "AbsoluteWeight": 0.01, "DifferentialWeightP": 1, "DifferentialWeightS": 0.5, "ResidualCutoff": 6, "MaxSeparationKm": 10},
    {"Iterations": 6, "AbsoluteWeight": 0.01, "DifferentialWeightP": 1, "DifferentialWeightS": 0.5, "ResidualCutoff": 5, "MaxSeparationKm": 6}]}}}
json.dump(config, open(os.path.join(folder, "dd_config.json"), "w"), indent=1)
ts_out = os.path.join(folder, "tomostar")
if not os.path.exists(os.path.join(ts_out, "locations.csv")):
    tomostar("relocate", "--method", "dd", "--config", os.path.join(folder, "dd_config.json"), "--stations", os.path.join(base, "stations.csv"), "--events", os.path.join(base, "events.csv"),
             "--picks", os.path.join(base, "picks.csv"), "--hypocentres", os.path.join(base, "tomostar", "events_relocated.csv"),
             "--grid", os.path.join(base, "grid.json"), "--model", "Carannante2013", "--set", "Locator.UseStationCorrections=false",
             "--set", "Tomography.ForwardRefinement=1", "--set", "Tomography.StationTerms=false", "--out", ts_out, "--name", "tomostar_dd", "--quiet")
ts = {}
with open(os.path.join(ts_out, "locations.csv")) as f:
    for r in csv.DictReader(f):
        if r.get("status", "ok") == "ok":
            ts[r["id"]] = (float(r["lat"]), float(r["lon"]), float(r["depth_km"]))

# ---- HypoDD -------------------------------------------------------------------------------------------
hd = work("dd", "hypodd")
labels = {sid: f"S{i:03d}" for i, sid in enumerate(sorted(stations))}  # HypoDD labels have at most 7 characters
with open(os.path.join(hd, "station.dat"), "w") as f:
    for sid, (lat, lon, elev) in sorted(stations.items()):
        f.write(f"{labels[sid]:7s} {lat:10.5f} {lon:10.5f} {elev:7.1f}\n")
with open(os.path.join(hd, "phase.pha"), "w") as f:
    for ev in ids:
        r = start[ev]
        t0 = parse_time(r["time"])
        f.write(f"# {t0:%Y %m %d %H %M} {t0.second + t0.microsecond / 1e6:6.3f} {float(r['lat']):9.5f} {float(r['lon']):10.5f} "
                f"{float(r['depth_km']):7.3f} 0.0 0.0 0.0 {float(r['rms_s']):5.2f} {number[ev]:9d}\n")
        for p in picks[ev]:
            tt = (parse_time(p["time"]) - t0).total_seconds()
            f.write(f"{labels[p['station']]:7s} {tt:8.3f} 1.000 {p['phase']}\n")
open(os.path.join(hd, "ph2dt.inp"), "w").write(
    "station.dat\nphase.pha\n*MINWGHT MAXDIST MAXSEP MAXNGH MINLNK MINOBS MAXOBS\n0 200 10 10 8 8 50\n")
# 1 km layers where the events are, thicker below (HypoDD reads each input line into 220 characters).
tops = [-3.0] + [float(z) for z in range(0, 20)] + [20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 40.0, 60.0]
mid = [(a + b) / 2 if a >= 0 else 0.0 for a, b in zip(tops, tops[1:] + [tops[-1] + 20])]  # layer centres
vp = [model(z)[0] for z in mid]
ratio = [model(z)[0] / model(z)[1] for z in mid]
damp = 80
sets = [(2, -9, -9), (4, 6, 10), (6, 5, 6)]
inp = ["hypoDD_2", "dt.cc", "dt.ct", "event.sel", "station.dat", "hypoDD.loc", "hypoDD.reloc", "hypoDD.sta", "hypoDD.res", "hypoDD.src",
       "2 3 200", "0 8 -999 -999 -999", f"2 2 0 {len(sets)}"]
inp += [f"{n} -9 -9 -9 -9 1.0 0.5 {wr} {wd} {damp}" for n, wr, wd in sets]
inp += ["1", " ".join(f"{z:.1f}" for z in tops) + " -9", " ".join(f"{v:.3f}" for v in vp) + " -9",
        " ".join(f"{v:.3f}" for v in ratio) + " -9", "0", ""]
assert max(map(len, inp)) < 220
open(os.path.join(hd, "hypoDD.inp"), "w").write("\n".join(inp))
if not os.path.exists(os.path.join(hd, "hypoDD.reloc")):
    subprocess.run([os.path.join(HYPODD, "ph2dt"), "ph2dt.inp"], cwd=hd, check=True, stdout=subprocess.DEVNULL)
    subprocess.run([os.path.join(HYPODD, "hypoDD"), "hypoDD.inp"], cwd=hd, check=True, stdout=open(os.path.join(hd, "run.log"), "w"))
back = {v: k for k, v in number.items()}
hdd = {}
clusters = {}
for line in open(os.path.join(hd, "hypoDD.reloc")):
    p = line.split()
    ev = back[int(p[0])]
    hdd[ev] = (float(p[1]), float(p[2]), float(p[3]))
    clusters[ev] = int(p[23])


def compare(reference, cluster_of):
    """Per-event difference after removing each cluster's mean offset: absolute (km) and relative."""
    lat0, lon0 = np.mean([v[0] for v in reference.values()]), np.mean([v[1] for v in reference.values()])
    groups = {}
    for ev in ts:
        if ev in reference:
            groups.setdefault(cluster_of(ev), []).append(ev)
    rows = []
    for members in groups.values():
        if len(members) < 2:
            continue
        a = np.array([local_xyz(*ts[ev], lat0, lon0) for ev in members])
        b = np.array([local_xyz(*reference[ev], lat0, lon0) for ev in members])
        offset = (a - b).mean(axis=0)
        for ev, x, y in zip(members, a - offset, b):
            d = float(np.linalg.norm(x - y))
            lat, lon, dep = reference[ev]
            scale = float(np.median([math.hypot(distance_km(lat, lon, *stations[p["station"]][:2]), dep + stations[p["station"]][2] / 1000)
                                     for p in picks[ev]]))
            rows.append((d, d / scale))
    d = np.array([r[0] for r in rows])
    rel = np.array([r[1] for r in rows])
    return {"events": len(rows), "clusters": len(groups), "median_km": float(np.median(d)), "p95_km": float(np.percentile(d, 95)),
            "median_relative": float(np.median(rel)), "p95_relative": float(np.percentile(rel, 95))}


# 1. Code against code.
s = compare(hdd, lambda ev: clusters[ev])
shift = [math.hypot(distance_km(hdd[ev][0], hdd[ev][1], float(start[ev]["lat"]), float(start[ev]["lon"])),
                    hdd[ev][2] - float(start[ev]["depth_km"])) for ev in hdd]
s["hypodd_median_shift_from_start_km"] = float(np.median(shift))
s["relocated_hypodd"] = len(hdd)
s["relocated_tomostar"] = len(ts)
record("double_difference", "double-difference relocation, catalogue differential times of the INGV picks",
       "HypoDD 2.1b (Waldhauser 2001), same picks, pairs, weights and model",
       "95th percentile of |hypocentre difference| / median hypocentral distance (cluster offsets removed)", s["p95_relative"], details=s)

# 2. Against the published CAT2 (HypoDD with cross-correlation delays).
cat2 = {}
with open(os.path.join(CATS, "CAT2")) as f:
    next(f)
    for line in f:
        p = line.split()
        cat2[p[0]] = (float(p[3]), float(p[4]), float(p[5]))
s = compare(cat2, lambda ev: 0)
ts_saved = ts
ts = dict(hdd)
s["hypodd_same_input_vs_cat2"] = {k: v for k, v in compare(cat2, lambda ev: 0).items() if k not in ("events", "clusters")}
ts = ts_saved
s["published_formal_errors_km"] = {"east": 0.110, "north": 0.120, "vertical": 0.162}
record("double_difference", "double-difference relocation, catalogue differential times of the INGV picks",
       "published catalogue CAT2 (Michele et al. 2020; HypoDD with cross-correlation delays)",
       "95th percentile of |hypocentre difference| / median hypocentral distance (mean offset removed)", s["p95_relative"],
       threshold=None, details=s)
