# Copyright 2026 Matteo Mangiagalli
# SPDX-License-Identifier: Apache-2.0
"""Shared helpers of the validation scripts: running the TomoSTAR exporter and recording results."""
import json
import os
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.environ.get("TOMOSTAR_VALIDATION_WORK", os.path.join(HERE, "work"))
RESULTS = os.path.join(HERE, "results")
EXPORTER = os.path.join(HERE, "TomoStar.Validation", "bin", "Release", "net10.0", "tomostar-validation")
CLI = os.path.join(HERE, "..", "src", "TomoStar.Cli", "bin", "Release", "net10.0", "tomostar")
THRESHOLD = 0.05  # the acceptance limit of every comparison: 5 % relative error


def tomostar_validation(*args):
    """Runs one command of the TomoSTAR validation exporter."""
    subprocess.run([EXPORTER, *map(str, args)], check=True)


def tomostar(*args):
    """Runs the tomostar program itself."""
    subprocess.run([CLI, *map(str, args)], check=True)


def work(*parts):
    path = os.path.join(WORK, *parts)
    os.makedirs(os.path.dirname(path) if os.path.splitext(path)[1] else path, exist_ok=True)
    return path


def record(case, method, reference, metric, value, threshold=THRESHOLD, details=None):
    """
    Stores the outcome of one comparison in results/<case>.json: the method of TomoSTAR, the
    reference it is compared with, the error metric, its value, the acceptance threshold and
    whether it passes. Every value is a fraction (0.01 = 1 %) unless the metric says otherwise.
    """
    os.makedirs(RESULTS, exist_ok=True)
    path = os.path.join(RESULTS, f"{case}.json")
    rows = json.load(open(path)) if os.path.exists(path) else []
    rows = [r for r in rows if not (r["method"] == method and r["reference"] == reference and r["metric"] == metric)]
    # threshold None: a comparison reported for information (the reference used other input, so the
    # difference measures the input rather than the code); it has no pass or fail.
    rows.append({"case": case, "method": method, "reference": reference, "metric": metric,
                 "value": float(value), "threshold": None if threshold is None else float(threshold),
                 "pass": None if threshold is None else bool(value < threshold), "details": details or {}})
    json.dump(rows, open(path, "w"), indent=2)
    status = "INFO" if threshold is None else ("PASS" if value < threshold else "FAIL")
    limit = "" if threshold is None else f" (limit {threshold:g})"
    print(f"[{status}] {case}: {method} vs {reference}: {metric} = {value:.4g}{limit}")
