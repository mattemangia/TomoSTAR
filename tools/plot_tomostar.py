#!/usr/bin/env python3
# Copyright 2026 Matteo Mangiagalli
# SPDX-License-Identifier: Apache-2.0
"""
Figures of TomoSTAR results: depth slices and vertical sections of any .qvol volume, masked where
the rays do not sample the model, with stations and hypocentres. Needs numpy and matplotlib only.

    python3 plot_tomostar.py slices VOLUME.qvol --dws DWS.qvol --depths 4,8,12 --out fig.png
        [--events events.csv] [--stations stations.csv] [--cmap RdBu] [--limits -8,8]
        [--title TEXT] [--label TEXT] [--min-dws 20]
    python3 plot_tomostar.py section VOLUME.qvol --from LON,LAT --to LON,LAT --out fig.png
        [--dws DWS.qvol] [--events events.csv] [--width-km 5] ...
    python3 plot_tomostar.py pair TRUE.qvol RECOVERED.qvol --dws DWS.qvol --depths 4,8 --out fig.png
        (a resolution test: true pattern above, recovered below)

The .qvol format (little-endian): 8-byte magic, int32 header length, UTF-8 JSON header, padding to
64 bytes, then Nx*Ny*Nz float32 values, longitude fastest, depth slowest.
"""
import argparse
import csv
import json
import math
import struct

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

R_EARTH = 6371.0


def read_qvol(path):
    """Header (dict) and values as an array [depth, lat, lon]."""
    with open(path, "rb") as f:
        raw = f.read()
    if raw[:4] != b"QVOL":
        raise SystemExit(f"{path} is not a .qvol volume")
    n = struct.unpack("<i", raw[8:12])[0]
    header = json.loads(raw[12:12 + n].decode("utf-8"))
    offset = 12 + n
    offset += (64 - offset % 64) % 64
    g = header["Grid"]
    nx, ny, nz = g["Nx"], g["Ny"], g["Nz"]
    values = np.frombuffer(raw, dtype="<f4", count=nx * ny * nz, offset=offset).astype(float).reshape(nz, ny, nx)
    return header, values


def axes(g):
    lon = np.linspace(g["MinLon"], g["MaxLon"], g["Nx"])
    lat = np.linspace(g["MinLat"], g["MaxLat"], g["Ny"])
    dep = np.linspace(g["MinDepthKm"], g["MaxDepthKm"], g["Nz"])
    return lon, lat, dep


def trilinear(values, lon_ax, lat_ax, dep_ax, lon, lat, dep):
    """Trilinear sample of a volume at arrays of points (clamped to the grid)."""
    def frac(ax, x):
        f = (np.asarray(x) - ax[0]) / (ax[1] - ax[0])
        f = np.clip(f, 0, len(ax) - 1 - 1e-9)
        i = np.floor(f).astype(int)
        return i, f - i
    i, tx = frac(lon_ax, lon)
    j, ty = frac(lat_ax, lat)
    k, tz = frac(dep_ax, dep)
    out = np.zeros(np.shape(lon))
    for dk in (0, 1):
        for dj in (0, 1):
            for di in (0, 1):
                w = (tx if di else 1 - tx) * (ty if dj else 1 - ty) * (tz if dk else 1 - tz)
                out += w * values[k + dk, j + dj, i + di]
    return out


def read_points(path, lon_key="lon", lat_key="lat", dep_key="depth_km"):
    if not path:
        return None
    lon, lat, dep = [], [], []
    with open(path, newline="") as f:
        for row in csv.DictReader(line for line in f if not line.startswith("#")):
            try:
                lon.append(float(row[lon_key]))
                lat.append(float(row[lat_key]))
                dep.append(float(row.get(dep_key) or 0) if dep_key in row else -float(row.get("elevation_m") or 0) / 1000)
            except (KeyError, ValueError):
                continue
    return np.array(lon), np.array(lat), np.array(dep)


def fine_slice(values, dws, min_dws, lon, lat, dep, z, n=240):
    """
    A depth slice sampled trilinearly on a fine mesh, blank where the (also interpolated) DWS is
    below min_dws: the field is interpolated first and masked afterwards, so the edge of the sampled
    region is smooth instead of following the grid cells.
    """
    flon = np.linspace(lon[0], lon[-1], n)
    flat = np.linspace(lat[0], lat[-1], n)
    LON, LAT = np.meshgrid(flon, flat)
    Z = np.full_like(LON, z)
    field = trilinear(values, lon, lat, dep, LON, LAT, Z)
    if dws is not None:
        field = np.where(trilinear(dws, lon, lat, dep, LON, LAT, Z) >= min_dws, field, np.nan)
    return flon, flat, field


def limits(values, given, symmetric):
    if given:
        lo, hi = (float(x) for x in given.split(","))
        return lo, hi
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return -1, 1
    if symmetric:
        m = np.percentile(np.abs(finite), 98)
        return -m, m
    return np.percentile(finite, 2), np.percentile(finite, 98)


def map_panel(ax, lon, lat, field, vmin, vmax, cmap, title, ev=None, st=None, z=None, dz=None):
    aspect = 1 / math.cos(math.radians(0.5 * (lat[0] + lat[-1])))
    img = ax.pcolormesh(lon, lat, np.ma.masked_invalid(field), cmap=cmap, vmin=vmin, vmax=vmax, shading="auto", rasterized=True)
    if ev is not None and z is not None:
        sel = np.abs(ev[2] - z) <= dz
        ax.scatter(ev[0][sel], ev[1][sel], s=2, c="k", alpha=0.5, linewidths=0)
    if st is not None:
        ax.scatter(st[0], st[1], marker="^", s=18, c="white", edgecolors="k", linewidths=0.6)
    ax.set_xlim(lon[0], lon[-1])
    ax.set_ylim(lat[0], lat[-1])
    ax.set_aspect(aspect)
    ax.set_title(title, fontsize=10)
    ax.tick_params(labelsize=8)
    return img


def cmd_slices(a):
    h, v = read_qvol(a.volume)
    g = h["Grid"]
    lon, lat, dep = axes(g)
    dws = read_qvol(a.dws)[1] if a.dws else None
    ev = read_points(a.events)
    st = read_points(a.stations, dep_key="elevation_m")
    depths = [float(x) for x in a.depths.split(",")]
    cmap = a.cmap or ("RdBu" if h["Quantity"].startswith("d") else "viridis")
    fields = [fine_slice(v, dws, a.min_dws, lon, lat, dep, z) for z in depths]
    vmin, vmax = limits(np.concatenate([f[2].ravel() for f in fields]), a.limits, h["Quantity"].startswith("d"))
    n = len(depths)
    fig, axs = plt.subplots(1, n, figsize=(3.4 * n + 0.8, 4.0), constrained_layout=True, squeeze=False)
    dz = 0.5 * (dep[1] - dep[0])
    for ax, z, (flon, flat, field) in zip(axs[0], depths, fields):
        img = map_panel(ax, flon, flat, field, vmin, vmax, cmap, f"{z:g} km", ev, st, z, dz)
    cb = fig.colorbar(img, ax=axs[0].tolist(), shrink=0.8)
    cb.set_label(a.label or f"{h['Quantity']} ({h['Units']})" if h["Units"] else h["Quantity"])
    if a.title:
        fig.suptitle(a.title, fontsize=11)
    fig.savefig(a.out, dpi=a.dpi)


def cmd_pair(a):
    ht, t = read_qvol(a.true)
    hr, r = read_qvol(a.recovered)
    g = ht["Grid"]
    lon, lat, dep = axes(g)
    dws = read_qvol(a.dws)[1] if a.dws else None
    depths = [float(x) for x in a.depths.split(",")]
    vmin, vmax = limits(t, a.limits, True)
    n = len(depths)
    fig, axs = plt.subplots(2, n, figsize=(3.4 * n + 0.8, 7.4), constrained_layout=True, squeeze=False)
    st = read_points(a.stations, dep_key="elevation_m")
    for c, z in enumerate(depths):
        flon, flat, ft = fine_slice(t, None, 0, lon, lat, dep, z)
        _, _, fr = fine_slice(r, dws, a.min_dws, lon, lat, dep, z)
        img = map_panel(axs[0][c], flon, flat, ft, vmin, vmax, "RdBu", f"true, {z:g} km", st=st)
        map_panel(axs[1][c], flon, flat, fr, vmin, vmax, "RdBu", f"recovered, {z:g} km", st=st)
    cb = fig.colorbar(img, ax=axs.ravel().tolist(), shrink=0.6)
    cb.set_label(a.label or f"{ht['Quantity']} ({ht['Units']})")
    if a.title:
        fig.suptitle(a.title, fontsize=11)
    fig.savefig(a.out, dpi=a.dpi)


def cmd_section(a):
    h, v = read_qvol(a.volume)
    g = h["Grid"]
    lon_ax, lat_ax, dep_ax = axes(g)
    dws = read_qvol(a.dws)[1] if a.dws else None
    (lo1, la1), (lo2, la2) = ([float(x) for x in p.split(",")] for p in (a.from_, a.to))
    kx = R_EARTH * math.pi / 180 * math.cos(math.radians(0.5 * (la1 + la2)))
    ky = R_EARTH * math.pi / 180
    length = math.hypot((lo2 - lo1) * kx, (la2 - la1) * ky)
    s = np.linspace(0, 1, 300)
    z = np.linspace(dep_ax[0], dep_ax[-1], 150)
    S, Z = np.meshgrid(s, z)
    lon = lo1 + S * (lo2 - lo1)
    lat = la1 + S * (la2 - la1)
    field = trilinear(v, lon_ax, lat_ax, dep_ax, lon, lat, Z)
    if dws is not None:
        field = np.where(trilinear(dws, lon_ax, lat_ax, dep_ax, lon, lat, Z) >= a.min_dws, field, np.nan)
    cmap = a.cmap or ("RdBu" if h["Quantity"].startswith("d") else "viridis")
    vmin, vmax = limits(field, a.limits, h["Quantity"].startswith("d"))
    fig, ax = plt.subplots(figsize=(8, 3.6), constrained_layout=True)
    img = ax.pcolormesh(s * length, z, np.ma.masked_invalid(field), cmap=cmap, vmin=vmin, vmax=vmax, shading="auto", rasterized=True)
    ev = read_points(a.events)
    if ev is not None:
        # Distance along and across the section of every event; those within the width are drawn.
        ex = (ev[0] - lo1) * kx
        ey = (ev[1] - la1) * ky
        ux, uy = (lo2 - lo1) * kx / length, (la2 - la1) * ky / length
        along = ex * ux + ey * uy
        across = -ex * uy + ey * ux
        sel = (np.abs(across) <= a.width_km) & (along >= 0) & (along <= length)
        ax.scatter(along[sel], ev[2][sel], s=2, c="k", alpha=0.5, linewidths=0)
    ax.set_xlim(0, length)
    ax.set_ylim(dep_ax[-1], dep_ax[0])
    ax.set_xlabel("distance along the section (km)")
    ax.set_ylabel("depth (km)")
    ax.set_aspect("equal")
    cb = fig.colorbar(img, ax=ax, shrink=0.9)
    cb.set_label(a.label or (f"{h['Quantity']} ({h['Units']})" if h["Units"] else h["Quantity"]))
    ax.set_title(a.title or f"{lo1:.2f}E {la1:.2f}N to {lo2:.2f}E {la2:.2f}N", fontsize=10)
    fig.savefig(a.out, dpi=a.dpi)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--dws", help="DWS volume: nodes below --min-dws are left blank")
    common.add_argument("--min-dws", type=float, default=10)
    common.add_argument("--stations")
    common.add_argument("--events")
    common.add_argument("--cmap")
    common.add_argument("--limits")
    common.add_argument("--label")
    common.add_argument("--title")
    common.add_argument("--dpi", type=int, default=150)
    common.add_argument("--out", required=True)
    s = sub.add_parser("slices", parents=[common])
    s.add_argument("volume")
    s.add_argument("--depths", required=True)
    s.set_defaults(fn=cmd_slices)
    s = sub.add_parser("pair", parents=[common])
    s.add_argument("true")
    s.add_argument("recovered")
    s.add_argument("--depths", required=True)
    s.set_defaults(fn=cmd_pair)
    s = sub.add_parser("section", parents=[common])
    s.add_argument("volume")
    s.add_argument("--from", dest="from_", required=True)
    s.add_argument("--to", required=True)
    s.add_argument("--width-km", type=float, default=5)
    s.set_defaults(fn=cmd_section)
    a = p.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()
