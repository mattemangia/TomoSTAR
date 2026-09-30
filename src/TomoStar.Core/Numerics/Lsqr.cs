// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

namespace TomoStar.Core.Numerics;

/// <summary>How an LSQR solve ended.</summary>
public sealed record LsqrResult(
    double[] Solution,
    int Iterations,
    double InitialResidualNorm,
    double ResidualNorm,
    double EstimatedConditionNumber,
    string StopReason);

/// <summary>
/// Damped LSQR (Paige &amp; Saunders 1982, ACM TOMS 8(1), 43-71; Algorithm 583, ACM TOMS 8(2),
/// 195-209). Solves
///
///     min ‖A x − b‖² + damp² ‖x‖²
///
/// by Golub-Kahan bidiagonalisation with plane rotations, following the SOL reference kernel that
/// SciPy's <c>scipy.sparse.linalg.lsqr</c> and MATLAB's <c>lsqr</c> implement. It is the solver of
/// SIMULPS-family local-earthquake tomography.
///
/// All three of the paper's stopping rules are applied: <c>atol</c> on ‖Aᵀr‖/(‖A‖‖r‖) (the test that
/// decides an inconsistent, i.e. every real, tomography system), <c>btol</c> on ‖r‖/‖b‖, and
/// <c>conlim</c> on the estimated condition number. A solve that stops on its iteration budget says
/// so, because a truncated LSQR is itself a regulariser nobody chose.
///
/// The dense vector updates run on <see cref="SimdVector"/>; the sparse products are row-parallel.
/// </summary>
public static class Lsqr
{
    public static LsqrResult Solve(
        CsrMatrix a,
        double[] b,
        double damp = 0,
        int maxIterations = 500,
        double atol = 1e-8,
        double btol = 1e-8,
        double conlim = 1e8,
        CancellationToken ct = default,
        Action<int, double>? onIteration = null)
    {
        var m = a.RowCount;
        var n = a.ColumnCount;
        if (b.Length != m) throw new ArgumentException("Right-hand side length differs from the row count.");
        var x = new double[n];
        if (m == 0 || n == 0) return new LsqrResult(x, 0, 0, 0, 0, "empty system");

        var u = (double[])b.Clone();
        var v = new double[n];
        var w = new double[n];
        var tmpM = new double[m];
        var tmpN = new double[n];

        var beta = SimdVector.Norm(u);
        if (beta == 0) return new LsqrResult(x, 0, 0, 0, 0, "zero right-hand side");
        var bnorm = beta;
        SimdVector.Scale(u, 1.0 / beta);
        a.TransposeMultiplyInto(u, v);
        var alpha = SimdVector.Norm(v);
        if (alpha == 0) return new LsqrResult(x, 0, bnorm, bnorm, 0, "Aᵀb = 0: the data do not see the model");
        SimdVector.Scale(v, 1.0 / alpha);
        Array.Copy(v, w, n);

        double rhobar = alpha, phibar = beta;
        double ddnorm = 0, xxnorm = 0, z = 0, cs2 = -1, sn2 = 0, anorm2 = 0, dampedRes2 = 0;
        double rnorm = beta, acond = 0;
        var ctol = conlim > 0 ? 1.0 / conlim : 0.0;
        var stop = "iteration budget exhausted";
        var iterations = 0;

        for (var it = 0; it < maxIterations; it++)
        {
            ct.ThrowIfCancellationRequested();
            iterations = it + 1;

            // β_{k+1} u_{k+1} = A v_k − α_k u_k
            a.MultiplyInto(v, tmpM);
            SimdVector.Xmby(tmpM, alpha, u);
            beta = SimdVector.Norm(u);
            if (beta > 0)
            {
                SimdVector.Scale(u, 1.0 / beta);
                anorm2 += alpha * alpha + beta * beta + damp * damp;
                // α_{k+1} v_{k+1} = Aᵀ u_{k+1} − β_{k+1} v_k
                a.TransposeMultiplyInto(u, tmpN);
                SimdVector.Xmby(tmpN, beta, v);
                alpha = SimdVector.Norm(v);
                if (alpha > 0) SimdVector.Scale(v, 1.0 / alpha);
            }

            // Fold the damping into the bidiagonal.
            var rhobar1 = Math.Sqrt(rhobar * rhobar + damp * damp);
            var cs1 = rhobar / rhobar1;
            var sn1 = damp / rhobar1;
            var psi = sn1 * phibar;
            phibar = cs1 * phibar;

            // Eliminate the subdiagonal β.
            var rho = Math.Sqrt(rhobar1 * rhobar1 + beta * beta);
            var cs = rhobar1 / rho;
            var sn = beta / rho;
            var theta = sn * alpha;
            rhobar = -cs * alpha;
            var phi = cs * phibar;
            phibar = sn * phibar;
            var tau = sn * phi;

            // x ← x + (φ/ρ) w ;  w ← v − (θ/ρ) w
            ddnorm += SimdVector.Dot(w, w) / (rho * rho);
            SimdVector.Axpy(phi / rho, w, x);
            SimdVector.Xpby(v, -theta / rho, w);

            // Norm estimates.
            var delta = sn2 * rho;
            var gambar = -cs2 * rho;
            var rhs = phi - delta * z;
            var gamma = Math.Sqrt(gambar * gambar + theta * theta);
            cs2 = gambar / gamma;
            sn2 = theta / gamma;
            z = rhs / gamma;
            xxnorm += z * z;

            // ‖r‖ of the augmented system includes what the damping absorbed (Σψ²).
            dampedRes2 += psi * psi;
            var anorm = Math.Sqrt(anorm2);
            acond = anorm * Math.Sqrt(ddnorm);
            var arnorm = alpha * Math.Abs(tau);
            rnorm = Math.Sqrt(phibar * phibar + dampedRes2);
            var xnorm = Math.Sqrt(xxnorm);

            var test1 = rnorm / bnorm;
            var test2 = arnorm / Math.Max(1e-300, anorm * rnorm);
            var test3 = 1.0 / Math.Max(1e-300, acond);
            var t1 = test1 / (1 + anorm * xnorm / bnorm);
            var rtol = btol + atol * anorm * xnorm / bnorm;

            onIteration?.Invoke(iterations, rnorm);

            if (test1 <= rtol) { stop = "residual reduced to tolerance (btol)"; break; }
            if (test2 <= atol) { stop = "normal equations satisfied (atol)"; break; }
            if (ctol > 0 && test3 <= ctol) { stop = "condition-number limit reached (conlim)"; break; }
            if (1 + t1 <= 1) { stop = "residual at machine precision"; break; }
            if (1 + test2 <= 1) { stop = "normal equations at machine precision"; break; }
            if (1 + test3 <= 1) { stop = "condition number at machine precision"; break; }
            if (beta == 0 || alpha == 0) { stop = "bidiagonalisation terminated"; break; }
        }

        return new LsqrResult(x, iterations, bnorm, rnorm, acond, stop);
    }
}
