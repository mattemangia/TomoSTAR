# Copyright 2026 Matteo Mangiagalli
# SPDX-License-Identifier: Apache-2.0
"""
Automatic P and S picking (STA/LTA trigger and AIC onset) on real data: the waveforms of the
M >= 3.5 earthquakes of the 2016-2017 central Italy sequence, against the picks of the INGV
analysts for the same records (the INGV bulletin, ISIDe Working Group 2007).

tomostar pick runs without the analyst picks, from the times the example's 3-D model predicts for
its relocated hypocentres, with a search window sized to that model's accuracy. Each automatic pick is matched with the analyst pick of the same event,
station and phase. Metric: |t_automatic - t_analyst| divided by the analyst's travel time (the
relative error the pick puts into the travel time the tomography fits); its 95th percentile must
be below 5 %, for P and for S, against the analyst picks of the bulletin's most precise class
(uncertainty 0.1 s): a reference whose own uncertainty is 0.3 to 1 s cannot certify 5 % of a travel
time of a few seconds. The comparison with all the analyst picks is reported for information.

Inputs: TOMOSTAR_NORCIA_OUT and TOMOSTAR_NORCIA_DATA (the data folder of examples/norcia2016).
"""
import csv
import os
from datetime import datetime

import numpy as np

from common import record, tomostar, work

NORCIA = os.environ["TOMOSTAR_NORCIA_OUT"]
DATA = os.environ["TOMOSTAR_NORCIA_DATA"]


def parse(t):
    return datetime.strptime(t.rstrip("Z")[:26], "%Y-%m-%dT%H:%M:%S.%f")


out = work("picker", "run")
if not os.path.exists(os.path.join(out, "picks.csv")):
    tomostar("pick", "--stations", os.path.join(NORCIA, "vel", "stations.csv"), "--events", os.path.join(NORCIA, "vel", "events_relocated.csv"),
             "--stationxml", os.path.join(DATA, "stations.xml"), "--waveforms", os.path.join(DATA, "waveforms"),
             "--grid", os.path.join(NORCIA, "grid", "grid.json"), "--model", os.path.join(NORCIA, "min1d", "model1d_minimum.txt"),
             "--vp", os.path.join(NORCIA, "vel", "volumes", "vel_Vp.qvol"), "--vs", os.path.join(NORCIA, "vel", "volumes", "vel_Vs.qvol"),
             # The search window sized to the predictions of the example's relocated 3-D model (RMS of its
             # residuals 0.10 s): +-1 s plus 3 % of the travel time, instead of the defaults (+-3 s plus 8 %)
             # meant for a starting 1-D model.
             "--set", "Picker.WindowSeconds=1", "--set", "Picker.WindowPerSecond=0.03",
             "--only-automatic", "--out", out, "--name", "pick", "--quiet")
origin = {}
with open(os.path.join(NORCIA, "vel", "events_relocated.csv")) as f:
    for e in csv.DictReader(f):
        origin[e["id"]] = parse(e["time"])
analyst = {}
with open(os.path.join(NORCIA, "reloc", "picks.csv")) as f:
    for p in csv.DictReader(f):
        if p["origin"] == "catalog" and p["disabled"] == "0":
            key = (p["event"], p["station"].upper(), p["phase"])
            t = parse(p["time"])
            if key not in analyst or t < analyst[key][0]:
                analyst[key] = (t, float(p["sigma_s"]))
results = {}
with open(os.path.join(out, "picks.csv")) as f:
    for p in csv.DictReader(f):
        if p["origin"] != "automatic":
            continue
        key = (p["event"], p["station"].upper(), p["phase"])
        if key not in analyst or p["event"] not in origin:
            continue
        t_ref, sigma = analyst[key]
        travel = (t_ref - origin[p["event"]]).total_seconds()
        if travel <= 0:
            continue
        d = (parse(p["time"]) - t_ref).total_seconds()
        results.setdefault(p["phase"], []).append((d, abs(d) / travel, sigma))


def summary(rows):
    d = np.array([x[0] for x in rows])
    rel = np.array([x[1] for x in rows])
    return float(np.percentile(rel, 95)), {
        "matched_picks": int(rel.size), "median_abs_s": float(np.median(np.abs(d))), "p95_abs_s": float(np.percentile(np.abs(d), 95)),
        "mean_s": float(d.mean()), "within_0.1_s": float(np.mean(np.abs(d) <= 0.1)), "median_relative": float(np.median(rel))}


for phase in ("P", "S"):
    rows = results.get(phase, [])
    value, details = summary([x for x in rows if x[2] <= 0.1])
    record("picker", f"automatic {phase} picks (STA/LTA trigger, AIC onset) on real waveforms",
           "INGV analyst picks of the same records, most precise class (0.1 s)", "95th percentile of |t_auto - t_analyst| / travel time",
           value, details=details)
    value, details = summary(rows)
    record("picker", f"automatic {phase} picks (STA/LTA trigger, AIC onset) on real waveforms",
           "INGV analyst picks of the same records, all classes (0.1 to 1 s)", "95th percentile of |t_auto - t_analyst| / travel time",
           value, threshold=None, details=details)
