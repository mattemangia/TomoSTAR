# Copyright 2026 Matteo Mangiagalli
# SPDX-License-Identifier: Apache-2.0
"""
Writes the table of results/*.json into docs/validation.md, between the lines
<!-- results begin --> and <!-- results end -->, so the document shows what the scripts measured.

    python3 report.py
"""
import glob
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DOC = os.path.join(HERE, "..", "docs", "validation.md")
ORDER = ["eikonal", "rays", "lsqr", "qtomography", "location", "double_difference", "tstar", "picker", "signal", "formats", "published_model"]


def fmt(v):
    if v is None:
        return ""
    if v == 0:
        return "0"
    if abs(v) < 1e-3:
        return f"{v:.1e}"
    return f"{100 * v:.2f} %" if abs(v) < 1 else f"{v:.3g}"


rows = []
for path in sorted(glob.glob(os.path.join(HERE, "results", "*.json")), key=lambda p: ORDER.index(os.path.basename(p)[:-5]) if os.path.basename(p)[:-5] in ORDER else 99):
    for r in json.load(open(path)):
        status = "information" if r["pass"] is None else ("pass" if r["pass"] else "FAIL")
        limit = "" if r["threshold"] is None else fmt(r["threshold"])
        cell = lambda t: str(t).replace("|", "\\|")  # a bar inside a cell would end it
        rows.append(f"| {cell(r['case'])} | {cell(r['method'])} | {cell(r['reference'])} | {cell(r['metric'])} | {fmt(r['value'])} | {limit} | {status} |")
table = "\n".join(["| Case | TomoSTAR method | Reference | Metric | Value | Limit | Result |", "|---|---|---|---|---|---|---|"] + rows)
text = open(DOC).read()
begin, end = "<!-- results begin -->", "<!-- results end -->"
head, rest = text.split(begin, 1)
_, tail = rest.split(end, 1)
open(DOC, "w").write(head + begin + "\n" + table + "\n" + end + tail)
print(f"{len(rows)} rows written to {os.path.relpath(DOC)}.")
