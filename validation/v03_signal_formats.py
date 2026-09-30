# Copyright 2026 Matteo Mangiagalli
# SPDX-License-Identifier: Apache-2.0
"""
Readers and signal processing against ObsPy (Beyreuther et al. 2010; Krischer et al. 2015) and SciPy
(Virtanen et al. 2020), on real records of the 2016-2017 central Italy sequence:
  - miniSEED samples and start times            vs obspy.read
  - QuakeML origins and associated picks         vs obspy.read_events
  - StationXML coordinates and response amplitude vs obspy.read_inventory + evalresp
  - instrument response removal to velocity       vs obspy Trace.remove_response
  - Butterworth band-pass, causal and zero phase  vs scipy.signal.butter + sosfilt
  - Slepian (DPSS) tapers of the multitaper spectra vs scipy.signal.windows.dpss
  - recursive STA/LTA and AIC of the picker       vs obspy.signal.trigger
"""
import glob
import json
import os

import numpy as np
import obspy
from obspy.signal.trigger import aic_simple, recursive_sta_lta
from scipy import signal

from common import record, tomostar_validation, work
from sample_data import sample

qml, stationxml, mseeds = sample()

# ---- miniSEED -------------------------------------------------------------------------------
worst = 0.0
n_traces = n_samples = 0
for m in mseeds:
    out = work("formats", os.path.basename(m) + ".d")
    tomostar_validation("mseed", m, out)
    ours = json.load(open(os.path.join(out, "traces.json")))
    ref = obspy.read(m)
    ref.merge(method=-1) if False else None
    ref_by = {}
    for tr in ref:
        ref_by.setdefault((tr.id, str(tr.stats.starttime)), tr)
    for t in ours:
        # Segments are matched by channel and start time; ObsPy keeps every record group separately too.
        key = (f"{t['Network']}.{t['Station']}.{t['Location']}.{t['Channel']}", str(obspy.UTCDateTime(t["Start"])))
        tr = ref_by.get(key)
        if tr is None:
            worst = 1.0
            continue
        x = np.fromfile(os.path.join(out, t["File"]), "<f4")
        y = tr.data.astype(np.float32)
        n = min(len(x), len(y))
        diff = np.max(np.abs(x[:n] - y[:n])) / max(1e-30, np.max(np.abs(y[:n]))) if n else 1.0
        if len(x) != len(y):
            diff = max(diff, abs(len(x) - len(y)) / len(y))
        worst = max(worst, diff)
        n_traces += 1
        n_samples += n
record("formats", "miniSEED reader (samples, start times)", "ObsPy read", "max relative sample difference", worst,
       details={"traces": n_traces, "samples": n_samples})

# ---- QuakeML --------------------------------------------------------------------------------
out = work("formats", "quakeml.json")
files = sorted(glob.glob(os.path.join(qml, "*.xml")))[:40]
tmp = work("formats", "qml_subset")
for f in files:
    dst = os.path.join(tmp, os.path.basename(f))
    if not os.path.exists(dst):
        os.symlink(os.path.abspath(f), dst)
tomostar_validation("quakeml", tmp, out)
ours = {e["Id"]: e for e in json.load(open(out))}
worst_loc = worst_pick = 0.0
matched = 0
for f in files:
    cat = obspy.read_events(f)
    for ev in cat:
        eid = str(ev.resource_id).split("=")[-1].split("/")[-1]
        e = ours.get(eid)
        if e is None:
            continue
        o = ev.preferred_origin() or ev.origins[0]
        worst_loc = max(worst_loc, abs(e["Lat"] - o.latitude), abs(e["Lon"] - o.longitude), abs(e["DepthKm"] - o.depth / 1000) / 100,
                        abs(obspy.UTCDateTime(e["Time"]) - o.time))
        # The picks the preferred origin associates, with the phase of the arrival.
        picks = {str(p.resource_id): p for p in ev.picks}
        for a in o.arrivals:
            p = picks.get(str(a.pick_id))
            if p is None or not a.phase or a.phase[0] not in "PS" or a.phase not in ("P", "Pg", "Pn", "Pb", "S", "Sg", "Sn", "Sb"):
                continue
            sid = f"{p.waveform_id.network_code}.{p.waveform_id.station_code}"
            cands = [q for q in e["Picks"] if q["StationId"] == sid and q["Phase"] == a.phase[0]]
            if not cands:
                continue
            worst_pick = max(worst_pick, min(abs(obspy.UTCDateTime(q["Time"]) - p.time) for q in cands))
            matched += 1
record("formats", "QuakeML reader (origins)", "ObsPy read_events", "max difference (degrees, s, and depth / 100 km)", worst_loc,
       details={"events": len(ours)})
record("formats", "QuakeML reader (associated P and S picks)", "ObsPy read_events", "max time difference (s)", worst_pick,
       details={"picks": matched})

# ---- StationXML and responses --------------------------------------------------------------------
freqs = np.logspace(np.log10(0.2), np.log10(20), 25)
out = work("formats", "stationxml.json")
tomostar_validation("stationxml", stationxml, ",".join(f"{f:.6g}" for f in freqs), out)
ours = json.load(open(out))
inv = obspy.read_inventory(stationxml)
worst_coord = 0.0
coords = {}
for net in inv:
    for sta in net:
        coords.setdefault(f"{net.code}.{sta.code}", []).append((sta.latitude, sta.longitude, sta.elevation, sta.start_date))
for s in ours["Stations"]:
    # TomoSTAR keeps the coordinates of the latest station epoch.
    lat, lon, elev, _ = max(coords[s["Id"]], key=lambda c: c[3])
    worst_coord = max(worst_coord, abs(s["Lat"] - lat), abs(s["Lon"] - lon), abs(s["ElevationM"] - elev) / 1000)
record("formats", "StationXML reader (station coordinates)", "ObsPy read_inventory", "max difference (degrees, km)", worst_coord,
       details={"stations": len(ours["Stations"])})
worst_amp = 0.0
n_resp = 0
for r in ours["Responses"]:
    if not r["Start"]:
        continue
    t = obspy.UTCDateTime(r["Start"]) + 1
    try:
        resp = inv.get_response(r["Nslc"], t)
        sr = inv.select(*r["Nslc"].split("."), time=t)[0][0][0].sample_rate
        ref = resp.get_evalresp_response_for_frequencies(freqs, output="VEL" if "S**2" not in r["InputUnits"].upper() else "ACC")
    except Exception:
        continue
    # The band where the digital FIR stages are flat: from 0.5 Hz (0.05 of the sampling rate for
    # long-period channels) to 0.4 of the sampling rate.
    band = (freqs >= min(0.5, 0.05 * sr)) & (freqs <= 0.4 * sr)
    if not band.any():
        continue
    a = np.array(r["Amplitude"])[band]
    rel = np.abs(a - np.abs(ref[band])) / np.abs(ref[band])
    worst_amp = max(worst_amp, rel.max())
    n_resp += 1
record("formats", "instrument response amplitude (poles and zeros, sensitivity)", "ObsPy evalresp (all stages)", "max relative error in the pass band",
       worst_amp, details={"channels": n_resp, "frequencies_hz": [float(freqs[0]), float(freqs[-1])]})

# ---- Response removal ---------------------------------------------------------------------------
pre = (0.2, 0.5, 20.0, 30.0)
worst_rr = 0.0
worst_long = 0.0  # records long enough to resolve the low corner of the pre-filter (10 periods)
n_rr = 0
for m in mseeds:
    out = work("formats", os.path.basename(m) + ".vel")
    tomostar_validation("remove-response", m, stationxml, *pre, 60, out)
    for t in json.load(open(os.path.join(out, "traces.json"))):
        st = obspy.read(m)
        tr = [x for x in st if x.id == t["Nslc"] and abs(x.stats.starttime - obspy.UTCDateTime(t["Start"])) < 1e-3]
        if not tr:
            continue
        tr = tr[0].copy()
        # The same conventions as TomoSTAR: demean and linear detrend, then a cosine taper over 5% of the
        # record at each end (ObsPy's taper_fraction is the total of both ends).
        tr.detrend("demean")
        tr.detrend("linear")
        tr.remove_response(inventory=inv, output="VEL", pre_filt=pre, water_level=60, taper=True, taper_fraction=0.1, zero_mean=True)
        y = tr.data
        x = np.fromfile(os.path.join(out, t["File"]), "<f4")
        n = min(len(x), len(y))
        core = slice(int(0.1 * n), int(0.9 * n))  # away from the tapered ends
        rel = np.linalg.norm(x[core] - y[core]) / np.linalg.norm(y[core])
        worst_rr = max(worst_rr, rel)
        if len(y) / tr.stats.sampling_rate >= 10 / pre[0]:
            worst_long = max(worst_long, rel)
        n_rr += 1
record("signal", "instrument response removal to velocity", "ObsPy remove_response", "max relative L2 difference", worst_rr,
       details={"traces": n_rr, "pre_filter_hz": pre, "water_level_db": 60, "max_on_records_of_50_s_or_more": worst_long,
                "note": "short records are gap fragments; the residual difference is the FFT length (TomoSTAR pads to a power of 2) and the digital FIR stages TomoSTAR does not model"})

# ---- Filters, tapers, STA/LTA, AIC on a real trace -------------------------------------------------
tr = max((x for m in mseeds for x in obspy.read(m) if x.stats.channel.endswith("Z")), key=lambda x: x.stats.npts)
x = tr.data.astype(np.float32)
x -= x.mean()
fs = tr.stats.sampling_rate
inp = work("signal", "trace.f32")
x.astype("<f4").tofile(inp)
worst_f = 0.0
for order in (2, 4):
    for zero in (0, 1):
        out = work("signal", f"filtered_{order}_{zero}.f32")
        tomostar_validation("filter", inp, fs, 1.0, 15.0, order, zero, out)
        ours = np.fromfile(out, "<f4").astype(float)
        sos = np.vstack([signal.butter(order, 15.0, "lowpass", fs=fs, output="sos"), signal.butter(order, 1.0, "highpass", fs=fs, output="sos")])
        ref = signal.sosfilt(sos, x.astype(float))
        if zero:
            ref = signal.sosfilt(sos, ref[::-1])[::-1]
        worst_f = max(worst_f, np.linalg.norm(ours - ref) / np.linalg.norm(ref))
record("signal", "Butterworth band-pass 1-15 Hz, orders 2 and 4, causal and zero phase", "SciPy butter + sosfilt", "max relative L2 difference", worst_f)

worst_d = 0.0
for n, nw, k in ((256, 2.5, 4), (512, 4, 7), (1000, 3, 5)):
    out = work("signal", f"dpss_{n}.f64")
    tomostar_validation("dpss", n, nw, k, out)
    ours = np.fromfile(out).reshape(k, n)
    ref = signal.windows.dpss(n, nw, k)
    for i in range(k):
        a = ours[i] / np.linalg.norm(ours[i])
        b = ref[i] / np.linalg.norm(ref[i])
        worst_d = max(worst_d, min(np.linalg.norm(a - b), np.linalg.norm(a + b)))
record("signal", "Slepian (DPSS) tapers", "SciPy dpss", "max L2 difference of the unit-norm tapers", worst_d)

sta, lta = 0.5, 5.0
out = work("signal", "stalta.f32")
tomostar_validation("stalta", inp, fs, sta, lta, out)
ours = np.fromfile(out, "<f4").astype(float)
ref = recursive_sta_lta(x.astype(float), int(sta * fs), int(lta * fs))
# The two recursions are identical but start from different states (TomoSTAR seeds the LTA with the
# energy of the first second, ObsPy with zero); the difference decays as exp(-t / LTA), about 5% after
# 3 LTA and 0.03% after 8, so the comparison starts after 8 LTA.
start = int(8 * lta * fs)
assert len(x) > start + 10 * fs, "the trace is too short for the STA/LTA comparison"
rel = np.abs(ours[start:] - ref[start:]) / np.maximum(ref[start:], 1e-12)
record("signal", "recursive STA/LTA", "ObsPy recursive_sta_lta", "max relative difference after 8 LTA", rel.max(),
       details={"samples_compared": int(len(x) - start)})

seg = x[int(15 * fs):int(40 * fs)]
seg_file = work("signal", "segment.f32")
seg.astype("<f4").tofile(seg_file)
out = work("signal", "aic.f64")
tomostar_validation("aic", seg_file, out)
ours = np.fromfile(out)
ref = aic_simple(seg.astype(float))
# TomoSTAR's AIC(k) splits the window before sample k (k samples in the first part, as in Maeda 1985);
# ObsPy's aic_simple splits it after sample k, so its curve is the same one sample earlier.
inner = slice(2, len(seg) - 3)
shifted = slice(3, len(seg) - 2)
rel = np.abs(ours[shifted] - ref[inner]) / np.maximum(np.abs(ref[inner]), 1e-12)
record("signal", "AIC function of the onset picker (Maeda 1985)", "ObsPy aic_simple (aligned by one sample)", "max relative difference", rel.max(),
       details={"onset_tomostar": int(np.argmin(ours[2:-2]) + 2), "onset_obspy_plus_one": int(np.argmin(ref[2:-2]) + 3)})
