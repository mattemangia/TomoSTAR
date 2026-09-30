using TomoStar.Core.Numerics;

namespace TomoStar.Core.Tomography;

/// <summary>Diagonal of the model resolution matrix, estimated, with the standard error of each entry.</summary>
public sealed record ResolutionDiagonalResult(double[] Diagonal, double[] StandardError, int Probes, int MeanLsqrIterations);

/// <summary>
/// Stochastic estimate of the diagonal of the model resolution matrix of a damped, smoothed
/// least-squares inversion, without forming the matrix.
///
/// The system solved by the inversion is the stacked, weighted matrix A = [W·G; λD; μL]. Its
/// resolution matrix is R = (AᵀA)⁻¹ (WG)ᵀ(WG). For any vector z, R·z is the LSQR solution of
/// A·m ≈ [W·G·z; 0]: the model recovered from the data that the model z would produce, with the same
/// regularisation. With random probes z of independent ±1 entries, diag(R) is estimated by the
/// average of z ⊙ R·z (Bekas, Kokiopoulou &amp; Saad 2007, Appl. Numer. Math. 57(11-12), 1214-1229,
/// doi:10.1016/j.apnum.2007.01.003), the approach MacCarthy, Borchers &amp; Aster (2011, J. Geophys.
/// Res. 116(B10), doi:10.1029/2011JB008234) applied to large geophysical inverse problems.
///
/// Each probe costs one LSQR solve. The error of an entry falls as 1/√(probes) and grows with the
/// off-diagonal mass of its row of R (smearing), so the standard error across probes is returned
/// with the estimate: values within about two standard errors of zero are not resolved.
/// </summary>
public static class ResolutionDiagonal
{
    /// <param name="a">The weighted system of the last linearised step, data rows first.</param>
    /// <param name="dataRows">How many leading rows of <paramref name="a"/> are data (the rest is regularisation).</param>
    public static ResolutionDiagonalResult Estimate(CsrMatrix a, int dataRows, int probes, int lsqrIterations, double tolerance, int seed,
        IProgress<(double, string)>? progress = null, CancellationToken ct = default)
    {
        if (dataRows <= 0 || dataRows > a.RowCount) throw new ArgumentOutOfRangeException(nameof(dataRows));
        probes = Math.Max(2, probes);
        var n = a.ColumnCount;
        var sum = new double[n];
        var sum2 = new double[n];
        var rng = new Random(seed);
        var z = new double[n];
        var az = new double[a.RowCount];
        var rhs = new double[a.RowCount];
        long iterations = 0;
        for (var k = 0; k < probes; k++)
        {
            ct.ThrowIfCancellationRequested();
            progress?.Report(((double)k / probes, $"Resolution diagonal: probe {k + 1}/{probes}"));
            for (var i = 0; i < n; i++) z[i] = rng.Next(2) == 0 ? -1 : 1;
            a.MultiplyInto(z, az);
            // Synthetic data of the probe model on the data rows only; regularisation rows ask for zero.
            Array.Clear(rhs);
            Array.Copy(az, rhs, dataRows);
            var sol = Lsqr.Solve(a, rhs, 0, lsqrIterations, tolerance, tolerance, 1e8, ct);
            iterations += sol.Iterations;
            for (var i = 0; i < n; i++)
            {
                var e = z[i] * sol.Solution[i];
                sum[i] += e;
                sum2[i] += e * e;
            }
        }
        var diag = new double[n];
        var err = new double[n];
        for (var i = 0; i < n; i++)
        {
            var mean = sum[i] / probes;
            diag[i] = mean;
            var variance = Math.Max(0, (sum2[i] - probes * mean * mean) / (probes - 1));
            err[i] = Math.Sqrt(variance / probes);
        }
        progress?.Report((1, "Resolution diagonal: done"));
        return new ResolutionDiagonalResult(diag, err, probes, (int)(iterations / probes));
    }
}
