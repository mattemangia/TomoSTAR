// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using System.Numerics;
using TomoStar.Core.Geo;
using TomoStar.Core.Numerics;

namespace TomoStar.Core.Tomography;

/// <summary>
/// Diagonal of the model resolution matrix and model standard deviations, computed exactly, for the
/// kept columns of a system (the model and the station terms), with the eliminated columns (the
/// hypocentres) separated out.
/// </summary>
/// <param name="Columns">The kept columns of the system, in order.</param>
/// <param name="Diagonal">Rᵢᵢ of each kept column.</param>
/// <param name="StandardDeviation">√(Cᵢᵢ), C = S⁻¹: the posterior standard deviation of the unknown, in its own units
/// (fractional change), from the data uncertainties the rows are weighted by, with the regularisation as the prior.</param>
/// <param name="LengthKm">Resolution length of each kept column with a place, km (NaN without one): the spread of its row of R
/// over the columns of the same block, √(Σⱼ Rᵢⱼ² dᵢⱼ² / Σⱼ Rᵢⱼ²).</param>
public sealed record FormalResolutionResult(int[] Columns, double[] Diagonal, double[] StandardDeviation, double[] LengthKm, double Seconds);

/// <summary>Where a column of the system stands (the centre of its node or cell) and which block (P, S, Q) it belongs to.</summary>
public readonly record struct ColumnPlace(Vec3 Position, int Block);

/// <summary>
/// The model resolution matrix and the posterior covariance of a damped, smoothed least-squares
/// inversion, computed exactly from the weighted system A = [W·G; regularisation] of the last
/// linearised step.
///
/// With N = AᵀA, D its data part and Γ = N − D its regularisation part, R = N⁻¹D and C = N⁻¹
/// (Menke 2012, Geophysical Data Analysis, ch. 4; Tarantola 2005, ch. 3, for C as the posterior covariance with
/// the regularisation as a Gaussian prior). Hypocentre columns are separated as in parameter separation (Pavlis
/// &amp; Booker 1980, J. Geophys. Res. 85(B9), 4801-4810; Spencer &amp; Gubbins 1980, Geophys. J. R. Astr. Soc. 63,
/// 95-116): N is reduced to the kept columns by the Schur complement S = N_kk − N_kh N_hh⁻¹ N_hk, one event (4×4
/// block) at a time. The regularisation never couples a kept with an eliminated column, so the data part reduces
/// to S − Γ_kk and the kept block of the resolution matrix is R = I − S⁻¹Γ_kk. Γ is sparse (damping, smoothing,
/// the zero-sum row of the station terms), so with S = LLᵀ and Y = L⁻¹:
/// Cᵢᵢ = ‖Y_{:,i}‖², Rᵢᵢ = 1 − Σⱼ Γᵢⱼ (Y_{:,i}·Y_{:,j}).
///
/// The cost is that of the dense factor: n² memory and about n³/3 + n³/6 operations for n kept columns, which is
/// why it is offered for the adaptive cells (a few thousand unknowns) rather than for fine regular grids. Double
/// differences, which couple the hypocentres of two events, are not supported.
/// </summary>
/// <summary>Whether to compute the exact resolution diagonal and standard deviations on the last system (<see cref="FormalResolution"/>).</summary>
public sealed class FormalResolutionSettings
{
    public bool Enabled { get; set; }

    /// <summary>Largest number of unknowns (model and station terms, without the hypocentres): memory is 16 n² bytes.</summary>
    public int MaxUnknowns { get; set; } = 20000;

    public FormalResolutionSettings Clone() => (FormalResolutionSettings)MemberwiseClone();
}

/// <summary>The resolution diagonal and the standard deviation of one block of unknowns, on the grid nodes.</summary>
/// <param name="Resolution">Rᵢᵢ of the unknown each node belongs to (its cell on the adaptive grid).</param>
/// <param name="Sigma">Its posterior standard deviation, as a fractional change.</param>
/// <param name="LengthKm">Its resolution length, km.</param>
public sealed record FormalResolutionMaps(double[] Resolution, double[] Sigma, double[] LengthKm);

public static class FormalResolution
{
    /// <summary>The values of a block on the nodes: node i takes those of column offset + columnOfNode(i).</summary>
    public static FormalResolutionMaps OnNodes(FormalResolutionResult r, int columnCount, int nodes, int offset, Func<int, int> columnOfNode)
    {
        var at = new int[columnCount];
        Array.Fill(at, -1);
        for (var k = 0; k < r.Columns.Length; k++) at[r.Columns[k]] = k;
        var res = new double[nodes];
        var sd = new double[nodes];
        var len = new double[nodes];
        for (var i = 0; i < nodes; i++)
        {
            var k = at[offset + columnOfNode(i)];
            res[i] = k >= 0 ? r.Diagonal[k] : double.NaN;
            sd[i] = k >= 0 ? r.StandardDeviation[k] : double.NaN;
            len[i] = k >= 0 ? r.LengthKm[k] : double.NaN;
        }
        return new FormalResolutionMaps(res, sd, len);
    }

    /// <summary>
    /// Places of the columns of a block: each column at the centroid of the grid nodes it stands for (a node, an
    /// adaptive cell, a layer), in geocentric km.
    /// </summary>
    public static void AddPlaces(ColumnPlace?[] places, SphericalGrid grid, int offset, int block, Func<int, int> columnOfNode)
    {
        var sum = new Dictionary<int, (Vec3 S, int N)>();
        for (var i = 0; i < grid.Count; i++)
        {
            var (x, y, z) = grid.Decompose(i);
            var c = offset + columnOfNode(i);
            var p = GeoMath.ToCartesian(grid.LonDeg[x], grid.LatDeg[y], grid.DepthKm[z]);
            var e = sum.GetValueOrDefault(c);
            sum[c] = (e.S + p, e.N + 1);
        }
        foreach (var (c, (v, n)) in sum) places[c] = new ColumnPlace(v / n, block);
    }

    /// <summary>Memory of the dense matrices, GB, for n kept columns.</summary>
    public static double MemoryGb(int n) => 2.0 * n * (double)n * 8 / 1e9;

    /// <param name="a">The weighted system.</param>
    /// <param name="isDataRow">Whether a row is data (the others are regularisation).</param>
    /// <param name="eliminate">Groups of columns to separate out (the 4 hypocentre columns of each event).</param>
    /// <param name="maxUnknowns">Largest number of kept columns accepted (memory is 16 n² bytes).</param>
    /// <param name="places">Optional place of each column of the system: with it, the resolution length of every placed column,
    /// from its full row of R (about n³/3 more operations).</param>
    public static FormalResolutionResult Compute(CsrMatrix a, Func<int, bool> isDataRow, IReadOnlyList<int[]> eliminate,
        int maxUnknowns, IProgress<(double, string)>? progress = null, CancellationToken ct = default, ColumnPlace?[]? places = null)
    {
        var clock = System.Diagnostics.Stopwatch.StartNew();
        var group = new int[a.ColumnCount];
        Array.Fill(group, -1);
        for (var g = 0; g < eliminate.Count; g++)
            foreach (var c in eliminate[g]) group[c] = g;
        var kept = Enumerable.Range(0, a.ColumnCount).Where(c => group[c] < 0).ToArray();
        var n = kept.Length;
        if (n > maxUnknowns)
            throw new InvalidOperationException(
                $"Formal resolution: {n} unknowns, above the limit of {maxUnknowns} ({MemoryGb(n):0.#} GB of dense matrices); " +
                "use the adaptive grid, or the checkerboard and spike tests.");
        var index = new int[a.ColumnCount];
        Array.Fill(index, -1);
        for (var k = 0; k < n; k++) index[kept[k]] = k;

        // Γ on the kept columns, and the checks that the separation is exact.
        var gamma = new Dictionary<int, double>[n];
        for (var k = 0; k < n; k++) gamma[k] = [];
        for (var r = 0; r < a.RowCount; r++)
        {
            int s0 = a.RowStart[r], s1 = a.RowStart[r + 1];
            var rowGroup = -1;
            bool touchesKept = false, touchesElim = false;
            for (var p = s0; p < s1; p++)
            {
                var g = group[a.ColumnIndex[p]];
                if (g < 0) { touchesKept = true; continue; }
                touchesElim = true;
                if (rowGroup >= 0 && g != rowGroup)
                    throw new NotSupportedException("Formal resolution: a row couples the hypocentres of two events (double differences); not supported.");
                rowGroup = g;
            }
            if (isDataRow(r)) continue;
            if (touchesKept && touchesElim)
                throw new InvalidOperationException("Formal resolution: a regularisation row couples model and hypocentre columns.");
            if (!touchesKept) continue;
            for (var p = s0; p < s1; p++)
            {
                var i = index[a.ColumnIndex[p]];
                for (var q = s0; q < s1; q++)
                {
                    var j = index[a.ColumnIndex[q]];
                    gamma[i][j] = gamma[i].GetValueOrDefault(j) + a.Values[p] * a.Values[q];
                }
            }
        }

        // S = N_kk (lower triangle), row by row from the transpose: no two threads write the same row.
        progress?.Report((0.05, "Formal resolution: normal matrix"));
        var t = a.Transpose();
        var S = new double[(long)n * n];
        Parallel.For(0, n, new ParallelOptions { CancellationToken = ct }, i =>
        {
            var c = kept[i];
            var row = (long)i * n;
            for (var p = t.RowStart[c]; p < t.RowStart[c + 1]; p++)
            {
                var r = t.ColumnIndex[p];
                var v = t.Values[p];
                for (var q = a.RowStart[r]; q < a.RowStart[r + 1]; q++)
                {
                    var j = index[a.ColumnIndex[q]];
                    if (j >= 0 && j <= i) S[row + j] += v * a.Values[q];
                }
            }
        });

        // Schur complement of each event: S −= U Uᵀ with U = N_kg L_g⁻ᵀ, L_g Lᵀ_g = N_gg.
        if (eliminate.Count > 0)
        {
            progress?.Report((0.15, "Formal resolution: separating the hypocentres"));
            var blocks = new (int[] Rows, double[] U, int M)[eliminate.Count];
            Parallel.For(0, eliminate.Count, new ParallelOptions { CancellationToken = ct }, g =>
            {
                var cols = eliminate[g];
                var m = cols.Length;
                var ngg = new double[m * m];
                var nkg = new Dictionary<int, double[]>();
                for (var l = 0; l < m; l++)
                for (var p = t.RowStart[cols[l]]; p < t.RowStart[cols[l] + 1]; p++)
                {
                    var r = t.ColumnIndex[p];
                    var v = t.Values[p];
                    for (var q = a.RowStart[r]; q < a.RowStart[r + 1]; q++)
                    {
                        var c = a.ColumnIndex[q];
                        if (group[c] == g) { ngg[l * m + Array.IndexOf(cols, c)] += v * a.Values[q]; continue; }
                        var k = index[c];
                        if (k < 0) continue;
                        if (!nkg.TryGetValue(k, out var e)) nkg[k] = e = new double[m];
                        e[l] += v * a.Values[q];
                    }
                }
                // Columns with no entry (a fixed event without damping) carry no information: left out.
                var live = Enumerable.Range(0, m).Where(l => ngg[l * m + l] > 0).ToArray();
                var ml = live.Length;
                var lg = new double[ml * ml];
                for (var x = 0; x < ml; x++)
                for (var y = 0; y < ml; y++) lg[x * ml + y] = ngg[live[x] * m + live[y]];
                if (ml == 0 || !SmallCholesky(lg, ml)) { blocks[g] = ([], [], 0); return; }
                var rows = nkg.Keys.Order().ToArray();
                var u = new double[rows.Length * ml];
                for (var z = 0; z < rows.Length; z++)
                {
                    var e = nkg[rows[z]];
                    for (var x = 0; x < ml; x++)
                    {
                        var sum = e[live[x]];
                        for (var y = 0; y < x; y++) sum -= lg[x * ml + y] * u[z * ml + y];
                        u[z * ml + x] = sum / lg[x * ml + x];
                    }
                }
                blocks[g] = (rows, u, ml);
            });
            // U by kept row, to subtract row by row without contention.
            var count = new int[n + 1];
            foreach (var (rows, _, m) in blocks)
                foreach (var k in rows) count[k + 1] += m;
            for (var k = 0; k < n; k++) count[k + 1] += count[k];
            var next = (int[])count.Clone();
            var uBlock = new int[count[n]];
            var uPos = new int[count[n]];
            for (var g = 0; g < blocks.Length; g++)
            {
                var (rows, _, m) = blocks[g];
                for (var z = 0; z < rows.Length; z++)
                for (var x = 0; x < m; x++)
                {
                    var p = next[rows[z]]++;
                    uBlock[p] = g;
                    uPos[p] = z * m + x;
                }
            }
            Parallel.For(0, n, new ParallelOptions { CancellationToken = ct }, i =>
            {
                var row = (long)i * n;
                for (var p = count[i]; p < count[i + 1]; p++)
                {
                    var (rows, u, m) = blocks[uBlock[p]];
                    var x = uPos[p] % m;
                    var v = u[uPos[p]];
                    for (var z = 0; z < rows.Length; z++)
                    {
                        var j = rows[z];
                        if (j > i) break;
                        S[row + j] -= v * u[z * m + x];
                    }
                }
            });
        }

        progress?.Report((0.3, "Formal resolution: Cholesky factor"));
        if (!Cholesky(S, n, ct))
            throw new InvalidOperationException("Formal resolution: the reduced normal matrix is not positive definite (an unknown without data and without damping).");
        progress?.Report((0.6, "Formal resolution: inverse of the factor"));
        // Y = Lᵀ⁻¹ stored by rows: row i of Y is column i of L⁻¹, entries i..n−1.
        var Y = new double[(long)n * n];
        Parallel.For(0, n, new ParallelOptions { CancellationToken = ct }, () => new double[n], (c, _, x) =>
        {
            Array.Clear(x);
            for (var k = c; k < n; k++)
            {
                var sum = k == c ? 1.0 : 0.0;
                sum -= DotLong(S, (long)k * n + c, x, c, k - c);
                x[k] = sum / S[(long)k * n + k];
            }
            Array.Copy(x, c, Y, (long)c * n + c, n - c);
            return x;
        }, _ => { });
        S = [];
        progress?.Report((0.9, "Formal resolution: diagonals"));
        var diag = new double[n];
        var sd = new double[n];
        Parallel.For(0, n, i =>
        {
            var yi = (long)i * n;
            sd[i] = Math.Sqrt(DotLong(Y, yi + i, Y, yi + i, n - i));
            double sum = 0;
            foreach (var (j, gij) in gamma[i])
            {
                var m = Math.Max(i, j);
                sum += gij * DotLong(Y, yi + m, Y, (long)j * n + m, n - m);
            }
            diag[i] = 1 - sum;
        });
        var length = new double[n];
        Array.Fill(length, double.NaN);
        if (places != null)
        {
            progress?.Report((0.93, "Formal resolution: resolution lengths"));
            var place = kept.Select(c => places[c]).ToArray();
            // Row i of R = eᵢ − (S⁻¹)ᵢ Γ, with (S⁻¹)ᵢⱼ = Yᵢ·Yⱼ; Γ is symmetric, so its rows serve as its columns.
            Parallel.For(0, n, new ParallelOptions { CancellationToken = ct }, () => (new double[n], new double[n]), (i, _, buf) =>
            {
                if (place[i] is not { } pi) return buf;
                var (sinv, row) = buf;
                var yi = (long)i * n;
                for (var j = 0; j < n; j++)
                {
                    var m = Math.Max(i, j);
                    sinv[j] = DotLong(Y, yi + m, Y, (long)j * n + m, n - m);
                }
                Array.Clear(row);
                row[i] = 1;
                for (var k = 0; k < n; k++)
                {
                    var v = sinv[k];
                    if (v == 0) continue;
                    foreach (var (j, g) in gamma[k]) row[j] -= v * g;
                }
                double num = 0, den = 0;
                for (var j = 0; j < n; j++)
                {
                    if (place[j] is not { } pj || pj.Block != pi.Block) continue;
                    var r2 = row[j] * row[j];
                    var d = pj.Position - pi.Position;
                    num += r2 * (d.X * d.X + d.Y * d.Y + d.Z * d.Z);
                    den += r2;
                }
                length[i] = den > 0 ? Math.Sqrt(num / den) : double.NaN;
                return buf;
            }, _ => { });
        }
        progress?.Report((1, "Formal resolution: done"));
        return new FormalResolutionResult(kept, diag, sd, length, clock.Elapsed.TotalSeconds);
    }

    private static double DotLong(double[] a, long ia, double[] b, long ib, int len)
    {
        var x = new ReadOnlySpan<double>(a, (int)ia, len);
        var y = new ReadOnlySpan<double>(b, (int)ib, len);
        var w = Vector<double>.Count;
        var acc = Vector<double>.Zero;
        var k = 0;
        for (; k <= len - w; k += w) acc += new Vector<double>(x.Slice(k, w)) * new Vector<double>(y.Slice(k, w));
        var sum = Vector.Dot(acc, Vector<double>.One);
        for (; k < len; k++) sum += x[k] * y[k];
        return sum;
    }

    /// <summary>In-place Cholesky of a small dense matrix (lower triangle); false when not positive definite.</summary>
    private static bool SmallCholesky(double[] m, int n)
    {
        for (var j = 0; j < n; j++)
        {
            var d = m[j * n + j];
            for (var k = 0; k < j; k++) d -= m[j * n + k] * m[j * n + k];
            if (d <= 0) return false;
            m[j * n + j] = Math.Sqrt(d);
            for (var i = j + 1; i < n; i++)
            {
                var s = m[i * n + j];
                for (var k = 0; k < j; k++) s -= m[i * n + k] * m[j * n + k];
                m[i * n + j] = s / m[j * n + j];
            }
        }
        return true;
    }

    /// <summary>
    /// Blocked, right-looking Cholesky factorisation in place (lower triangle of a row-major n × n matrix): factor
    /// a diagonal block, solve the panel below it, subtract the panel's outer product from the trailing matrix
    /// (row-parallel, so the result does not depend on the number of threads).
    /// </summary>
    internal static bool Cholesky(double[] A, int n, CancellationToken ct = default)
    {
        const int B = 96;
        for (var k0 = 0; k0 < n; k0 += B)
        {
            ct.ThrowIfCancellationRequested();
            var k1 = Math.Min(n, k0 + B);
            // Diagonal block.
            for (var j = k0; j < k1; j++)
            {
                var rj = (long)j * n;
                var d = A[rj + j] - DotLong(A, rj + k0, A, rj + k0, j - k0);
                if (!(d > 0)) return false;
                d = Math.Sqrt(d);
                A[rj + j] = d;
                for (var i = j + 1; i < k1; i++)
                {
                    var ri = (long)i * n;
                    A[ri + j] = (A[ri + j] - DotLong(A, ri + k0, A, rj + k0, j - k0)) / d;
                }
            }
            if (k1 == n) break;
            // Panel: rows below the block.
            Parallel.For(k1, n, i =>
            {
                var ri = (long)i * n;
                for (var j = k0; j < k1; j++)
                {
                    var rj = (long)j * n;
                    A[ri + j] = (A[ri + j] - DotLong(A, ri + k0, A, rj + k0, j - k0)) / A[rj + j];
                }
            });
            // Trailing update, lower triangle.
            var w = k1 - k0;
            Parallel.For(k1, n, i =>
            {
                var ri = (long)i * n;
                for (var j = k1; j <= i; j++)
                    A[ri + j] -= DotLong(A, ri + k0, A, (long)j * n + k0, w);
            });
        }
        return true;
    }
}
