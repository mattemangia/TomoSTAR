// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using TomoStar.Core.Forward;
using TomoStar.Core.Geo;
using TomoStar.Core.Numerics;
using TomoStar.Core.Tomography;

namespace TomoStar.Tests;

/// <summary>The numerical kernels against answers known in closed form.</summary>
public class NumericsTests
{
    [Fact]
    public void SimdKernelsAgreeWithScalarLoops()
    {
        var rnd = new Random(3);
        // A length that is not a multiple of any vector width, so the scalar tail is exercised too.
        var a = Enumerable.Range(0, 1037).Select(_ => rnd.NextDouble() - 0.5).ToArray();
        var b = Enumerable.Range(0, 1037).Select(_ => rnd.NextDouble() - 0.5).ToArray();
        var dot = a.Zip(b, (x, y) => x * y).Sum();
        Assert.Equal(dot, SimdVector.Dot(a, b), 10);
        Assert.Equal(Math.Sqrt(a.Sum(x => x * x)), SimdVector.Norm(a), 10);
        var y2 = (double[])b.Clone();
        SimdVector.Axpy(0.7, a, y2);
        for (var i = 0; i < a.Length; i++) Assert.Equal(b[i] + 0.7 * a[i], y2[i], 12);
    }

    [Fact]
    public void LsqrSolvesAnOverdeterminedSystemAsTheNormalEquations()
    {
        // A random 60 x 8 system: LSQR must reach the least-squares solution of the normal equations.
        var rnd = new Random(7);
        const int m = 60, n = 8;
        var dense = new double[m, n];
        var builder = new CsrBuilder(n);
        var rhs = new double[m];
        for (var i = 0; i < m; i++)
        {
            var cols = new List<int>();
            var vals = new List<double>();
            for (var j = 0; j < n; j++)
            {
                dense[i, j] = rnd.NextDouble() - 0.5;
                cols.Add(j);
                vals.Add(dense[i, j]);
            }
            rhs[i] = rnd.NextDouble();
            builder.AddRow(cols.ToArray(), vals.ToArray(), rhs[i]);
        }
        var (g, b) = builder.Build();
        var x = Lsqr.Solve(g, b, 0, 200, 1e-12, 1e-12).Solution;

        // Normal equations by Gaussian elimination.
        var ata = new double[n, n + 1];
        for (var r = 0; r < n; r++)
        {
            for (var c = 0; c < n; c++)
                for (var i = 0; i < m; i++) ata[r, c] += dense[i, r] * dense[i, c];
            for (var i = 0; i < m; i++) ata[r, n] += dense[i, r] * rhs[i];
        }
        for (var c = 0; c < n; c++)
        for (var r = 0; r < n; r++)
        {
            if (r == c) continue;
            var f = ata[r, c] / ata[c, c];
            for (var k = c; k <= n; k++) ata[r, k] -= f * ata[c, k];
        }
        for (var j = 0; j < n; j++) Assert.Equal(ata[j, n] / ata[j, j], x[j], 6);
    }

    [Fact]
    public void DampedLsqrShrinksTheSolution()
    {
        var b = new CsrBuilder(2);
        b.AddRow([0, 1], [1.0, 1.0], 2);
        b.AddRow([0, 1], [1.0, -1.0], 0);
        var (g, rhs) = b.Build();
        var free = Lsqr.Solve(g, rhs, 0, 50, 1e-12, 1e-12).Solution;
        var damped = Lsqr.Solve(g, rhs, 1, 50, 1e-12, 1e-12).Solution;
        Assert.Equal(1, free[0], 8);
        // With damping 1: (A'A + I) x = A'b, A'A = 2I, so x = A'b / 3 = (2/3, 2/3).
        Assert.Equal(2.0 / 3, damped[0], 6);
        Assert.Equal(2.0 / 3, damped[1], 6);
    }

    [Fact]
    public void FastMarchingMatchesStraightRaysInAHomogeneousModel()
    {
        // 60 x 60 x 30 km at 1 km spacing, 6 km/s: first-arrival times are distance / velocity.
        var def = GridDefinition.FromSpacing(13, 13.74, 43, 43.54, 0, 30, 1, 1);
        var g = new SphericalGrid(def);
        var s = Enumerable.Repeat(1 / 6.0, g.Count).ToArray();
        var seeds = SphericalEikonal.PointSource(g, s, 13.37, 43.27, 10);
        var t = SphericalEikonal.Solve(new GridMetric(g), s, seeds);
        var src = GeoMath.ToCartesian(13.37, 43.27, 10);
        double worst = 0;
        foreach (var (i, j, k) in new[] { (0, 0, 0), (g.Nx - 1, g.Ny - 1, 0), (g.Nx - 1, 0, g.Nz - 1), (g.Nx / 2, g.Ny - 1, g.Nz / 2) })
        {
            var exact = (g.Cartesian(i, j, k) - src).Length / 6;
            worst = Math.Max(worst, Math.Abs(t[g.Index(i, j, k)] - exact) / exact);
        }
        Assert.True(worst < 0.01, $"relative error {worst:P2}");
    }

    [Fact]
    public void TrilinearWeightsSumToOneAndInterpolateLinearFieldsExactly()
    {
        var g = new SphericalGrid(new GridDefinition { MinLon = 10, MaxLon = 11, MinLat = 40, MaxLat = 41, MinDepthKm = 0, MaxDepthKm = 20, Nx = 5, Ny = 6, Nz = 7 });
        var field = new double[g.Count];
        for (var n = 0; n < g.Count; n++)
        {
            var (i, j, k) = g.Decompose(n);
            field[n] = 2 * g.LonDeg[i] - g.LatDeg[j] + 0.1 * g.DepthKm[k];
        }
        Span<int> nodes = stackalloc int[8];
        Span<double> w = stackalloc double[8];
        g.TrilinearWeights(10.37, 40.61, 7.3, nodes, w);
        double sum = 0;
        foreach (var x in w) sum += x;
        Assert.Equal(1, sum, 12);
        Assert.Equal(2 * 10.37 - 40.61 + 0.73, g.Interpolate(field, 10.37, 40.61, 7.3), 9);
    }

    /// <summary>
    /// The L-curve corner is the point farthest from the chord of the scaled log-log curve, and does not move when the
    /// residual norms of the flat branch change in their third decimal (the finite-difference curvature it replaced did).
    /// </summary>
    [Fact]
    public void LCurveCornerIsStableOnTheBranches()
    {
        // Damping ascending: model norm falls steeply, then the residual grows with little change of the model.
        (double, double)[] curve = [(10.00, 100), (10.01, 30), (10.03, 9), (10.1, 3), (15, 2.5), (30, 2.2), (60, 2)];
        Assert.Equal(3, TravelTimeTomography.Corner(curve));
        var perturbed = curve.Select((p, i) => i < 3 ? (p.Item1 * (1 + 0.0005 * (i % 2 == 0 ? 1 : -1)), p.Item2) : p).ToArray();
        Assert.Equal(3, TravelTimeTomography.Corner(perturbed));
    }
}
