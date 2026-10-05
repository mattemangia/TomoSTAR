# Copyright 2026 Matteo Mangiagalli
# SPDX-License-Identifier: Apache-2.0
"""
t* from P spectra inverted jointly per event, on real data: the M >= 3.5 earthquakes of the
2016-2017 central Italy sequence (waveforms of the INGV networks, picks of the INGV bulletin).

Reference software (peer reviewed method): AttenTIon, the t* inversion of Stachnik et al. (2004,
JGR 109, B10304) in the Python version of Wei and Wiens (2018, EPSL 502, 187-199), commit e34f1d8 of
https://github.com/swei-seismo/AttenTIon (environment variable ATTENTION_DIR). Its functions are
imported unchanged; only its multitaper module, which estimates spectra from SAC files with a
compiled library, is not loaded: both programs are given the same spectra, those TomoSTAR measures
on the real waveforms (the multitaper estimator itself is checked against SciPy in v03).

Two comparisons, for alpha = 0.27 (the default of AttenTIon) and non-negative t* in both:
  1. the joint inversion at the same corner frequency: AttenTIon's inversion() given TomoSTAR's
     corner frequency and the same fitted bands. Metric: for each event ||t*_TomoSTAR - t*_AttenTIon||
     / ||t*_AttenTIon||; its maximum over the events must be below 5 %.
  2. the whole event, each program with its own corner-frequency search (AttenTIon: bestfc(),
     minimum of ||r|| / sum(d); TomoSTAR: minimum of the mean squared residual), with the same search
     range (the corners of stress drops 0.1 to 100 MPa for the event's magnitude, TomoSTAR's range,
     given to AttenTIon through its stress-drop and magnitude parameters) and the same records for
     the search (TomoSTAR's high-quality subset). Reported for information: the criteria differ, so
     the difference measures the criterion. The default settings of AttenTIon (0.5 to 20 MPa, an mb
     to Mw conversion) put its corner at an end of its grid for many of these events; that run is
     reported in the details.
  3. the corner-frequency search itself: AttenTIon's corner grid, data vector and matrix (buildd,
     buildG) and nnls with the least-squares criterion, against TomoSTAR's search. Metric: the median
     over the events of the relative t* difference, below 5 %.
"""
import json
import os
import sys
import types

import numpy as np

from common import record, tomostar_validation, work

ATTENTION = os.environ["ATTENTION_DIR"]
NORCIA = os.environ["TOMOSTAR_NORCIA_OUT"]
NORCIA_DATA = os.environ["TOMOSTAR_NORCIA_DATA"]
ALPHA = 0.27
SCALE = 1e20  # amplitudes scaled so that ln(moment) is positive, as nnls requires; it shifts ln(moment) only

sys.modules["multitaper"] = types.ModuleType("multitaper")  # spectra are given, the compiled estimator is not needed
sys.path.insert(0, ATTENTION)
import tstar_inversion_function as tf  # noqa: E402
import tstar_parameters as tp  # noqa: E402
import tstarsub  # noqa: E402
from scipy.optimize import nnls  # noqa: E402

folder = work("tstar")
tp.resultdir = folder
tp.logfl = open(os.path.join(folder, "attention.log"), "w")
tp.fclist = open(os.path.join(folder, "attention_fc.lst"), "w")
param = tp.set_parameters()
param.update({"alpha": ALPHA, "doplotfcall": False, "doplotfcts": False, "source_para": 1})

out = os.path.join(folder, "tomostar.json")
if not os.path.exists(out):
    tomostar_validation("tstar-spectra", os.path.join(NORCIA, "vel", "stations.csv"), os.path.join(NORCIA, "vel", "events_relocated.csv"),
                        os.path.join(NORCIA, "reloc", "picks.csv"), os.path.join(NORCIA_DATA, "stations.xml"),
                        os.path.join(NORCIA_DATA, "waveforms"), ALPHA, out)
events = json.load(open(out))["Events"]


def corner_range(magnitude):
    """TomoSTAR's corner-frequency search range (TStarEstimator.CornerRange with the default settings)."""
    m0 = 10 ** (1.5 * magnitude + 9.1)
    fc = lambda mpa: 0.37 * 3500 * np.cbrt(16 * mpa * 1e6 / (7 * m0))
    return max(0.3, fc(0.1)), min(40.0, fc(100.0))


def attention_range(ev):
    """AttenTIon's magnitude and stress-drop parameters that give TomoSTAR's corner range in bestfc()."""
    lo, hi = corner_range(ev["Magnitude"])
    mb = (ev["Magnitude"] + 2.54) / 1.54  # bestfc converts mb to Mw = 1.54 mb - 2.54
    mo = 10 ** (1.5 * ev["Magnitude"] + 9.095)
    stress = lambda f: (f / (0.49 * param["beta"] * 100)) ** 3 * mo
    return mb, [stress(lo), stress(hi)]


def good(r):
    """TomoSTAR's records for the corner search: SNR at least twice the minimum, band at least 1.5 times."""
    f = np.array(r["Freq"])[np.array(r["Used"], bool)]
    return r.get("Snr", np.inf) >= 2 * 3 and f[-1] / f[0] >= 1.5 * 3


def saving_of(ev, search=None):
    """AttenTIon's per-station store, filled with TomoSTAR's spectra (fitted bins only) for cases 1 and 2."""
    saving, stations = {}, []
    for r in ev["Records"]:
        if search is not None and r["Station"] not in search:
            continue
        used = np.array(r["Used"], bool)
        f = np.array(r["Freq"])[used]
        a = np.array(r["Amplitude"])[used] * SCALE
        sta = r["Station"]
        saving[sta] = {"corr": [1.0 / r["DistanceKm"]], "Ptt": 1.0, 1: {"p": [f, a], "good": [True, False]},
                       2: {"p": [f, a], "good": [True, False]}, 3: {}}
        stations.append(sta)
    return saving, stations


same_fc, own_fc, fc_ratio, default_fc, default_ratio, ls_ratio, ls_fc, nrec = [], [], [], [], [], [], [], 0


def least_squares_corner(saving, stations, orig):
    """AttenTIon's corner grid and its own d and G, with the least-squares criterion (minimum ||r||)."""
    mw = 1.54 * orig["mb"] - 2.54
    mo = 10 ** (1.5 * mw + 9.095)
    lo = 0.49 * ((param["dstress"][0] / mo) ** (1 / 3)) * param["beta"] * 100
    hi = 0.49 * ((param["dstress"][1] / mo) ** (1 / 3)) * param["beta"] * 100
    grid = np.hstack((np.arange(lo, min(hi, 1.09), 0.02), np.arange(max(lo, 1.1), hi, 0.1), [hi]))
    g = tstarsub.buildG(saving, stations, param["alpha"], "P", 1, param)[:, :, 0]
    norms = [nnls(g, tstarsub.buildd(saving, stations, orig, "P", 1, param, fc)[:, 0])[1] for fc in grid]
    return float(grid[int(np.argmin(norms))])
for ev in events:
    if len(ev["Records"]) < 5:
        continue
    ours = np.array([r["TStar"] for r in ev["Records"]])
    saving, stations = saving_of(ev)
    orig = {"orid": ev["Id"], "mb": ev["Magnitude"], "mw": -1, "mo": -1, "fc": ev["CornerHz"]}
    # 1. AttenTIon's inversion at TomoSTAR's corner frequency.
    orig, saving = tf.inversion(ev["Id"], saving, stations, orig, "P", 2, param)
    theirs = np.array([saving[s][2]["tstar"][0] for s in stations])
    same_fc.append(np.linalg.norm(ours - theirs) / max(np.linalg.norm(theirs), 1e-9))
    nrec += len(stations)
    # 2. AttenTIon's own corner frequency search (on TomoSTAR's search records and range), then its inversion.
    for settings in ("same", "default"):
        if settings == "same":
            chosen = {r["Station"] for r in ev["Records"] if good(r)}
            if len(chosen) < 3:
                chosen = {r["Station"] for r in ev["Records"]}
            mb, param["dstress"] = attention_range(ev)
        else:
            chosen, mb, param["dstress"] = None, ev["Magnitude"], [0.5, 20.0]
        saving, stations = saving_of(ev, chosen)
        orig = {"orid": ev["Id"], "mb": mb, "mw": -1, "mo": -1, "fc": -1}
        orig, flag = tf.bestfc(ev["Id"], saving, stations, orig, "P", 1, param)
        if not flag:
            continue
        fc = orig["fc"]
        if settings == "same":
            # The same search with the least-squares criterion: the corner then differs only by the grids.
            fc_ls = least_squares_corner(saving, stations, orig)
            ls_ratio.append(ev["CornerHz"] / fc_ls)
            sv, st = saving_of(ev)
            o = {"orid": ev["Id"], "mb": mb, "mw": -1, "mo": -1, "fc": fc_ls}
            o, sv = tf.inversion(ev["Id"], sv, st, o, "P", 2, param)
            t_ls = np.array([sv[x][2]["tstar"][0] for x in st])
            ls_fc.append(np.linalg.norm(ours - t_ls) / max(np.linalg.norm(t_ls), 1e-9))
        saving, stations = saving_of(ev)
        orig = {"orid": ev["Id"], "mb": mb, "mw": -1, "mo": -1, "fc": fc}
        orig, saving = tf.inversion(ev["Id"], saving, stations, orig, "P", 2, param)
        theirs = np.array([saving[s][2]["tstar"][0] for s in stations])
        (own_fc if settings == "same" else default_fc).append(np.linalg.norm(ours - theirs) / max(np.linalg.norm(theirs), 1e-9))
        (fc_ratio if settings == "same" else default_ratio).append(ev["CornerHz"] / fc)

same_fc, own_fc, fc_ratio, default_fc = np.array(same_fc), np.array(own_fc), np.array(fc_ratio), np.array(default_fc)
ls_fc, ls_ratio = np.array(ls_fc), np.array(ls_ratio)
record("tstar", "joint multitaper t* inversion of an event (Brune source, shared ln Omega0, alpha 0.27, t* >= 0)",
       "AttenTIon inversion() (Stachnik et al. 2004; Wei and Wiens 2018), same spectra and corner frequency",
       "max over events of ||t*_TomoSTAR - t*_AttenTIon|| / ||t*_AttenTIon||", same_fc.max(),
       details={"events": int(len(same_fc)), "records": nrec, "median": float(np.median(same_fc))})
record("tstar", "joint multitaper t* inversion with its own corner-frequency search",
       "AttenTIon bestfc() + inversion(), same spectra, search records and corner range", "median over events of ||t*_TomoSTAR - t*_AttenTIon|| / ||t*_AttenTIon||",
       float(np.median(own_fc)), threshold=None,
       details={"events": int(len(own_fc)), "p90": float(np.percentile(own_fc, 90)),
                "corner_ratio_median": float(np.median(fc_ratio)), "corner_ratio_p10_p90": [float(np.percentile(fc_ratio, 10)), float(np.percentile(fc_ratio, 90))],
                "attention_default_settings": {"median": float(np.median(default_fc)), "p90": float(np.percentile(default_fc, 90)),
                                               "corner_ratio_median": float(np.median(default_ratio))},
                "note": "the two programs choose the corner frequency by different criteria (AttenTIon: ||r|| / sum(d); TomoSTAR: mean squared residual)"})
record("tstar", "joint multitaper t* inversion with its corner-frequency search (least squares)",
       "AttenTIon d, G and nnls on its corner grid with the least-squares criterion, same spectra, search records and range",
       "median over events of ||t*_TomoSTAR - t*_AttenTIon|| / ||t*_AttenTIon||", float(np.median(ls_fc)),
       details={"events": int(len(ls_fc)), "p90": float(np.percentile(ls_fc, 90)), "corner_ratio_median": float(np.median(ls_ratio)),
                "corner_ratio_p10_p90": [float(np.percentile(ls_ratio, 10)), float(np.percentile(ls_ratio, 90))]})
