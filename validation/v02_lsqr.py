# Copyright 2026 Matteo Mangiagalli
# SPDX-License-Identifier: Apache-2.0
"""
LSQR (Paige and Saunders 1982) against SciPy's scipy.sparse.linalg.lsqr (Virtanen et al. 2020) on
the real tomography problem of the 2016-2017 central Italy example: the matrix G of the last
iteration of its Vp and Vs inversion (velocity, hypocentre, station-term and regularisation rows;
252,053 rows, 48,256 columns) and, as the right-hand side, its real travel-time residuals divided by
their uncertainty (zero on the regularisation rows), for three damping values.

Inputs: TOMOSTAR_NORCIA_OUT, the out folder of examples/norcia2016/norcia.tomo.
"""
import csv
import json
import os

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import lsqr

from common import record, tomostar_validation, work

RUN = os.path.join(os.environ["TOMOSTAR_NORCIA_OUT"], "vel")


def read_qcsr(path):
    with open(path, "rb") as f:
        assert f.read(4) == b"QCSR"
        rows, cols, nnz = np.frombuffer(f.read(12), "<i4")
        start = np.frombuffer(f.read(4 * (rows + 1)), "<i4")
        ci = np.frombuffer(f.read(4 * nnz), "<i4")
        v = np.frombuffer(f.read(4 * nnz), "<f4").astype(float)
    return csr_matrix((v, ci, start), shape=(rows, cols))


path = os.path.join(RUN, "G.qcsr")
a = read_qcsr(path)
layout = json.load(open(os.path.join(RUN, "G.layout.json")))
with open(os.path.join(RUN, "residuals.csv")) as f:
    residuals = [(float(r["final_residual_s"]), float(r["sigma_s"])) for r in csv.DictReader(f)]
b = np.zeros(a.shape[0])
for row, arrival in enumerate(layout["RowArrival"]):
    res, sigma = residuals[arrival]
    b[row] = res / sigma
rhs = work("lsqr", "b.f64")
b.astype("<f8").tofile(rhs)
worst = 0.0
details = {}
for damp in (0.0, 1.0, 10.0):
    ref = lsqr(a, b, damp=damp, atol=1e-12, btol=1e-12, conlim=1e12, iter_lim=20000)[0]
    out = work("lsqr", "x.f64")
    tomostar_validation("lsqr", path, rhs, damp, 20000, out)
    x = np.fromfile(out)
    rel = np.linalg.norm(x - ref) / np.linalg.norm(ref)
    worst = max(worst, rel)
    details[f"damping {damp:g}"] = {"rows": a.shape[0], "columns": a.shape[1], "nonzeros": int(a.nnz), "data_rows": layout["DataRows"],
                                    "relative_difference": float(rel)}
record("lsqr", "damped LSQR on the real tomography problem (central Italy 2016-2017)", "SciPy lsqr", "max ||x - x_ref|| / ||x_ref||",
       worst, details=details)
