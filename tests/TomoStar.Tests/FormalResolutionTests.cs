// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using TomoStar.Core.Geo;
using TomoStar.Core.Numerics;
using TomoStar.Core.Tomography;

namespace TomoStar.Tests;

/// <summary>
/// The exact resolution diagonal and standard deviations, with the hypocentres separated by Schur
/// complements, against the kept block of the full R = N⁻¹D and C = N⁻¹ computed by brute force.
/// </summary>
public class FormalResolutionTests
{
    private const int Model = 230, Events = 12, Stations = 9;

    /// <summary>
    /// A system shaped like the joint inversion: data rows over a random run of model columns, the 4
    /// hypocentre columns of one event and one station term; damping and first-difference smoothing on
    /// the model, damping on the hypocentres and the station terms, and the zero-sum row of the terms.
    /// </summary>
    private static (CsrMatrix A, int DataRows, List<int[]> Groups) System(int seed, bool coupleEvents = false)
    {
        var rnd = new Random(seed);
        var offH = Model;
        var offC = Model + 4 * Events;
        var b = new CsrBuilder(offC + Stations);
        for (var r = 0; r < 600; r++)
        {
            var e = rnd.Next(Events);
            var start = rnd.Next(Model - 40);
            var len = 5 + rnd.Next(35);
            var row = new SortedDictionary<int, double>();
            for (var k = start; k < start + len; k++) row[k] = 0.2 + rnd.NextDouble();
            for (var h = 0; h < 4; h++) row[offH + 4 * e + h] = h == 0 ? 1 : rnd.NextDouble() - 0.5;
            if (coupleEvents && r == 0) row[offH + 4 * ((e + 1) % Events)] = -1;
            row[offC + rnd.Next(Stations)] = 1;
            b.AddRow(row, rnd.NextDouble() - 0.5, 1 / (0.05 + 0.1 * rnd.NextDouble()));
        }
        var dataRows = b.RowCount;
        Regularization.AddDamping(b, 0, Model, 3);
        for (var k = 0; k + 1 < Model; k++) b.AddRow([k, k + 1], [1.0, -1.0], 0, 5);
        Regularization.AddDamping(b, offH, 4 * Events, 0.5);
        Regularization.AddDamping(b, offC, Stations, 2);
        b.AddRow(Enumerable.Range(offC, Stations).ToArray(), Enumerable.Repeat(1.0, Stations).ToArray(), 0, 10);
        var groups = Enumerable.Range(0, Events).Select(e => Enumerable.Range(offH + 4 * e, 4).ToArray()).ToList();
        return (b.Build().Matrix, dataRows, groups);
    }

    private static double[,] Normal(CsrMatrix a, Func<int, bool> rows)
    {
        var n = a.ColumnCount;
        var m = new double[n, n];
        for (var r = 0; r < a.RowCount; r++)
        {
            if (!rows(r)) continue;
            for (var p = a.RowStart[r]; p < a.RowStart[r + 1]; p++)
            for (var q = a.RowStart[r]; q < a.RowStart[r + 1]; q++)
                m[a.ColumnIndex[p], a.ColumnIndex[q]] += a.Values[p] * a.Values[q];
        }
        return m;
    }

    private static double[,] Inverse(double[,] m)
    {
        var n = m.GetLength(0);
        var a = (double[,])m.Clone();
        var inv = new double[n, n];
        for (var i = 0; i < n; i++) inv[i, i] = 1;
        for (var c = 0; c < n; c++)
        {
            var piv = Enumerable.Range(c, n - c).MaxBy(r => Math.Abs(a[r, c]));
            for (var k = 0; k < n; k++)
            {
                (a[c, k], a[piv, k]) = (a[piv, k], a[c, k]);
                (inv[c, k], inv[piv, k]) = (inv[piv, k], inv[c, k]);
            }
            var d = a[c, c];
            for (var k = 0; k < n; k++) { a[c, k] /= d; inv[c, k] /= d; }
            for (var r = 0; r < n; r++)
            {
                if (r == c || a[r, c] == 0) continue;
                var f = a[r, c];
                for (var k = 0; k < n; k++) { a[r, k] -= f * a[c, k]; inv[r, k] -= f * inv[c, k]; }
            }
        }
        return inv;
    }

    [Fact]
    public void MatchesTheKeptBlockOfTheFullResolutionAndCovariance()
    {
        var (a, dataRows, groups) = System(7);
        var n = a.ColumnCount;
        var nInv = Inverse(Normal(a, _ => true));
        var d = Normal(a, r => r < dataRows);
        var res = FormalResolution.Compute(a, r => r < dataRows, groups, 10_000);
        Assert.Equal(Model + Stations, res.Columns.Length);
        double worstR = 0, worstS = 0;
        for (var k = 0; k < res.Columns.Length; k++)
        {
            var c = res.Columns[k];
            double rcc = 0;
            for (var j = 0; j < n; j++) rcc += nInv[c, j] * d[j, c];
            worstR = Math.Max(worstR, Math.Abs(res.Diagonal[k] - rcc));
            worstS = Math.Max(worstS, Math.Abs(res.StandardDeviation[k] / Math.Sqrt(nInv[c, c]) - 1));
        }
        Assert.True(worstR < 1e-9, $"largest error of R_ii {worstR:G3}");
        Assert.True(worstS < 1e-9, $"largest relative error of the standard deviation {worstS:G3}");
        // The data resolve the model columns partly: the diagonal lies in [0, 1) (0 where no row reaches) and is not trivial.
        var model = res.Diagonal.Take(Model).ToArray();
        Assert.True(model.Min() > -1e-9 && model.Max() < 1 && model.Max() - model.Min() > 0.1, $"R_ii from {model.Min():0.000} to {model.Max():0.000}");
    }

    /// <summary>The resolution length from the full rows of R, over the placed columns of the same block only.</summary>
    [Fact]
    public void ResolutionLengthsMatchTheFullRows()
    {
        var (a, dataRows, groups) = System(13);
        var n = a.ColumnCount;
        var nInv = Inverse(Normal(a, _ => true));
        var d = Normal(a, r => r < dataRows);
        // Model columns along a line, 1 km apart, in two blocks; station terms without a place.
        var places = new ColumnPlace?[n];
        for (var c = 0; c < Model; c++) places[c] = new ColumnPlace(new Vec3(c % (Model / 2), 0, 0), c < Model / 2 ? 0 : 1);
        var res = FormalResolution.Compute(a, r => r < dataRows, groups, 10_000, places: places);
        var kept = res.Columns;
        double worst = 0;
        for (var k = 0; k < kept.Length; k++)
        {
            var i = kept[k];
            if (places[i] is not { } pi) { Assert.True(double.IsNaN(res.LengthKm[k])); continue; }
            double num = 0, den = 0;
            foreach (var j in kept)
            {
                if (places[j] is not { } pj || pj.Block != pi.Block) continue;
                double rij = 0;
                for (var q = 0; q < n; q++) rij += nInv[i, q] * d[q, j];
                num += rij * rij * Math.Pow(pj.Position.X - pi.Position.X, 2);
                den += rij * rij;
            }
            var expected = den > 0 ? Math.Sqrt(num / den) : double.NaN;
            if (double.IsNaN(expected)) continue;
            worst = Math.Max(worst, Math.Abs(res.LengthKm[k] - expected) / Math.Max(1, expected));
        }
        Assert.True(worst < 1e-8, $"largest relative error of the resolution length {worst:G3}");
    }

    [Fact]
    public void WithoutEliminatedColumnsItIsTheFullResolution()
    {
        var (a, dataRows, _) = System(11);
        var nInv = Inverse(Normal(a, _ => true));
        var d = Normal(a, r => r < dataRows);
        var res = FormalResolution.Compute(a, r => r < dataRows, [], 10_000);
        for (var c = 0; c < a.ColumnCount; c++)
        {
            double rcc = 0;
            for (var j = 0; j < a.ColumnCount; j++) rcc += nInv[c, j] * d[j, c];
            Assert.InRange(res.Diagonal[c] - rcc, -1e-9, 1e-9);
        }
    }

    [Fact]
    public void RowsCouplingTwoEventsAreRefused()
    {
        var (a, dataRows, groups) = System(3, coupleEvents: true);
        Assert.Throws<NotSupportedException>(() => FormalResolution.Compute(a, r => r < dataRows, groups, 10_000));
    }

    [Fact]
    public void TooManyUnknownsAreRefusedBeforeAllocating()
    {
        var (a, dataRows, groups) = System(5);
        var e = Assert.Throws<InvalidOperationException>(() => FormalResolution.Compute(a, r => r < dataRows, groups, 100));
        Assert.Contains("adaptive grid", e.Message);
    }
}
