# Copyright 2026 Matteo Mangiagalli
# SPDX-License-Identifier: Apache-2.0
"""
LSQR (Paige and Saunders 1982) against SciPy's scipy.sparse.linalg.lsqr (Virtanen et al. 2020) on a
real tomography matrix: the sensitivity matrix G of a TomoSTAR run (P and S velocity, hypocentre
and regularisation rows) on a synthetic data set, with right-hand sides b = G x + noise, for three
damping values. The matrix of the Norcia example is used too when TOMOSTAR_NORCIA_G points to it.
"""
import os

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import lsqr

from common import record, tomostar, tomostar_validation, work


def read_qcsr(path):
    with open(path, "rb") as f:
        assert f.read(4) == b"QCSR"
        rows, cols, nnz = np.frombuffer(f.read(12), "<i4")
        start = np.frombuffer(f.read(4 * (rows + 1)), "<i4")
        ci = np.frombuffer(f.read(4 * nnz), "<i4")
        v = np.frombuffer(f.read(4 * nnz), "<f4").astype(float)
    return csr_matrix((v, ci, start), shape=(rows, cols))


data = work("lsqr", "data")
run = work("lsqr", "run")
if not os.path.exists(os.path.join(run, "G.qcsr")):
    tomostar("synth", "--out", data, "--events", 80, "--stations", 20, "--quiet")
    tomostar("invert", data, "--grid", os.path.join(data, "grid.json"), "--model", os.path.join(data, "model1d.txt"),
             "--iterations", 1, "--out", run, "--no-opencl", "--quiet")
matrices = [("synthetic", os.path.join(run, "G.qcsr"))]
if os.environ.get("TOMOSTAR_NORCIA_G"):
    matrices.append(("Norcia 2016", os.environ["TOMOSTAR_NORCIA_G"]))

worst = 0.0
details = {}
rng = np.random.default_rng(1)
for name, path in matrices:
    a = read_qcsr(path)
    x_true = rng.normal(size=a.shape[1])
    b = a @ x_true
    b += 0.05 * np.linalg.norm(b) / np.sqrt(len(b)) * rng.normal(size=len(b))
    rhs = work("lsqr", "b.f64")
    b.astype("<f8").tofile(rhs)
    for damp in (0.0, 1.0, 10.0):
        ref = lsqr(a, b, damp=damp, atol=1e-12, btol=1e-12, conlim=1e12, iter_lim=20000)[0]
        out = work("lsqr", "x.f64")
        tomostar_validation("lsqr", path, rhs, damp, 20000, out)
        x = np.fromfile(out)
        rel = np.linalg.norm(x - ref) / np.linalg.norm(ref)
        worst = max(worst, rel)
        details[f"{name}, damping {damp:g}"] = {"rows": a.shape[0], "columns": a.shape[1], "nonzeros": int(a.nnz), "relative_difference": float(rel)}
record("lsqr", "damped LSQR on tomography matrices", "SciPy lsqr", "max ||x - x_ref|| / ||x_ref||", worst, details=details)
