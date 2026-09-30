using System.Numerics;

namespace TomoStar.Core.Signal;

/// <summary>
/// Multitaper spectral estimation (Thomson, 1982, Proc. IEEE 70(9), 1055-1096).
///
/// The discrete prolate spheroidal sequences (Slepian, 1978, Bell Syst. Tech. J. 57(5), 1371-1430)
/// of length N and half-bandwidth W = NW/N are the eigenvectors of the symmetric tridiagonal matrix
///
///     diagonal     ((N − 1)/2 − n)² cos 2πW,        n = 0 … N−1
///     off-diagonal n (N − n) / 2,                   n = 1 … N−1
///
/// (Percival &amp; Walden, 1993, Spectral Analysis for Physical Applications, §8.3), ordered by
/// decreasing eigenvalue. The eigenvalues are found by bisection on the Sturm sequence and the
/// vectors by inverse iteration, O(N) per step; the concentration ratios λₖ (energy inside ±W) come
/// from the tapers themselves. The K = 2NW − 1 best-concentrated tapers give nearly independent
/// spectra of the same record, whose weighted mean has about 1/K the variance of a single
/// tapered periodogram, at the cost of a resolution of ±W.
/// </summary>
public static class Multitaper
{
    /// <summary>The first <paramref name="k"/> Slepian tapers of length <paramref name="n"/> (unit energy) and their concentrations.</summary>
    public static (double[][] Tapers, double[] Concentration) Tapers(int n, double nw, int k)
    {
        if (n < 8) throw new ArgumentOutOfRangeException(nameof(n), "At least 8 samples are needed.");
        if (nw <= 0 || nw >= n / 2.0) throw new ArgumentOutOfRangeException(nameof(nw));
        k = Math.Clamp(k, 1, n);
        var w = nw / n;
        var diag = new double[n];
        var off = new double[n]; // off[i] couples i−1 and i
        var cos = Math.Cos(2 * Math.PI * w);
        for (var i = 0; i < n; i++)
        {
            var c = (n - 1) / 2.0 - i;
            diag[i] = c * c * cos;
            if (i > 0) off[i] = i * (double)(n - i) / 2;
        }
        // Gershgorin bounds for the bisection.
        double lo = double.MaxValue, hi = double.MinValue;
        for (var i = 0; i < n; i++)
        {
            var r = (i > 0 ? Math.Abs(off[i]) : 0) + (i + 1 < n ? Math.Abs(off[i + 1]) : 0);
            lo = Math.Min(lo, diag[i] - r);
            hi = Math.Max(hi, diag[i] + r);
        }
        var tapers = new double[k][];
        var lambda = new double[k];
        for (var m = 0; m < k; m++)
        {
            // The (m+1)-th largest eigenvalue: the number of eigenvalues above it is m.
            var ev = Eigenvalue(diag, off, n - 1 - m, lo, hi);
            var v = InverseIteration(diag, off, ev, m);
            // Sign convention (Percival & Walden): symmetric tapers sum positive, antisymmetric
            // ones start positive.
            var sum = v.Sum();
            if ((m % 2 == 0 && sum < 0) || (m % 2 == 1 && FirstLobe(v) < 0))
                for (var i = 0; i < n; i++) v[i] = -v[i];
            tapers[m] = v;
            lambda[m] = Concentration(v, w);
        }
        return (tapers, lambda);
    }

    /// <summary>
    /// Multitaper amplitude spectrum √S(f) of a record, with S the eigenvalue-weighted mean of the K
    /// eigenspectra, scaled like <see cref="Fft.AmplitudeSpectrum"/> (≈ continuous Fourier amplitude).
    /// </summary>
    public static (double[] Freq, double[] Amp) AmplitudeSpectrum(ReadOnlySpan<float> x, double sampleRate, double nw = 4, int k = 0, int minLength = 0)
    {
        var n = x.Length;
        if (k <= 0) k = Math.Max(1, (int)Math.Floor(2 * nw) - 1);
        var (tapers, lambda) = Tapers(n, nw, k);
        var nfft = Fft.NextPow2(Math.Max(n, minLength));
        var half = nfft / 2;
        var power = new double[half];
        double lambdaSum = 0;
        var buf = new Complex[nfft];
        for (var m = 0; m < tapers.Length; m++)
        {
            Array.Clear(buf);
            // Tapers have unit energy; √N makes a flat taper of the same energy equal 1, so the
            // result compares with an untapered transform of the record.
            var scale = Math.Sqrt(n);
            for (var i = 0; i < n; i++) buf[i] = new Complex(x[i] * tapers[m][i] * scale, 0);
            Fft.Transform(buf);
            for (var f = 0; f < half; f++) power[f] += lambda[m] * (buf[f].Real * buf[f].Real + buf[f].Imaginary * buf[f].Imaginary);
            lambdaSum += lambda[m];
        }
        var freq = new double[half];
        var amp = new double[half];
        for (var f = 0; f < half; f++)
        {
            freq[f] = f * sampleRate / nfft;
            amp[f] = Math.Sqrt(power[f] / lambdaSum) / sampleRate;
        }
        return (freq, amp);
    }

    /// <summary>Number of eigenvalues of the tridiagonal matrix below <paramref name="x"/> (Sturm count).</summary>
    private static int CountBelow(double[] d, double[] e, double x)
    {
        var count = 0;
        var q = d[0] - x;
        if (q < 0) count++;
        for (var i = 1; i < d.Length; i++)
        {
            if (Math.Abs(q) < 1e-300) q = 1e-300;
            q = d[i] - x - e[i] * e[i] / q;
            if (q < 0) count++;
        }
        return count;
    }

    /// <summary>The eigenvalue with <paramref name="index"/> eigenvalues below it (0 = smallest), by bisection.</summary>
    private static double Eigenvalue(double[] d, double[] e, int index, double lo, double hi)
    {
        for (var it = 0; it < 200; it++)
        {
            var mid = 0.5 * (lo + hi);
            if (CountBelow(d, e, mid) > index) hi = mid; else lo = mid;
            if (hi - lo <= 1e-13 * Math.Max(1, Math.Abs(mid))) break;
        }
        return 0.5 * (lo + hi);
    }

    /// <summary>Eigenvector of a tridiagonal matrix for a known eigenvalue, by inverse iteration.</summary>
    private static double[] InverseIteration(double[] d, double[] e, double lambda, int seed)
    {
        var n = d.Length;
        // Shift a hair off the eigenvalue so the system stays solvable.
        var shift = lambda + 1e-10 * Math.Max(1, Math.Abs(lambda));
        var rnd = new Random(1234 + seed);
        var v = Enumerable.Range(0, n).Select(_ => rnd.NextDouble() - 0.5).ToArray();
        var a = new double[n];
        var b = new double[n];
        var c = new double[n];
        for (var it = 0; it < 6; it++)
        {
            for (var i = 0; i < n; i++)
            {
                a[i] = i > 0 ? e[i] : 0;
                b[i] = d[i] - shift;
                c[i] = i + 1 < n ? e[i + 1] : 0;
            }
            v = SolveTridiagonal(a, b, c, v);
            var norm = Math.Sqrt(v.Sum(x => x * x));
            for (var i = 0; i < n; i++) v[i] /= norm;
        }
        return v;
    }

    /// <summary>Tridiagonal solve with partial pivoting (the shifted matrix is nearly singular by design).</summary>
    private static double[] SolveTridiagonal(double[] a, double[] b, double[] c, double[] rhs)
    {
        // Banded Gaussian elimination with row swaps: after pivoting a row can reach two
        // super-diagonals, so a (n × 3) band is kept.
        var n = b.Length;
        var u0 = new double[n];
        var u1 = new double[n];
        var u2 = new double[n];
        var r = (double[])rhs.Clone();
        for (var i = 0; i < n; i++) { u0[i] = b[i]; u1[i] = c[i]; u2[i] = 0; }
        var sub = (double[])a.Clone();
        for (var i = 0; i < n - 1; i++)
        {
            if (Math.Abs(sub[i + 1]) > Math.Abs(u0[i]))
            {
                // Swap rows i and i+1.
                (u0[i], sub[i + 1]) = (sub[i + 1], u0[i]);
                (u1[i], u0[i + 1]) = (u0[i + 1], u1[i]);
                (u2[i], u1[i + 1]) = (u1[i + 1], u2[i]);
                (r[i], r[i + 1]) = (r[i + 1], r[i]);
            }
            if (Math.Abs(u0[i]) < 1e-300) u0[i] = 1e-300;
            var m = sub[i + 1] / u0[i];
            u0[i + 1] -= m * u1[i];
            u1[i + 1] -= m * u2[i];
            r[i + 1] -= m * r[i];
        }
        if (Math.Abs(u0[n - 1]) < 1e-300) u0[n - 1] = 1e-300;
        var x = new double[n];
        for (var i = n - 1; i >= 0; i--)
        {
            var s = r[i];
            if (i + 1 < n) s -= u1[i] * x[i + 1];
            if (i + 2 < n) s -= u2[i] * x[i + 2];
            x[i] = s / u0[i];
        }
        return x;
    }

    private static double FirstLobe(double[] v)
    {
        foreach (var x in v) if (Math.Abs(x) > 1e-12) return x;
        return 0;
    }

    /// <summary>Fraction of a taper's energy inside the band ±W (its concentration λ).</summary>
    private static double Concentration(double[] v, double w)
    {
        // λ = Σᵢ Σⱼ vᵢ vⱼ sin(2πW(i−j)) / (π(i−j)), with 2W on the diagonal; the sum over lags uses the
        // autocorrelation of the taper, O(N²) but N is a few hundred samples.
        var n = v.Length;
        double total = 2 * w; // lag 0 (Σ vᵢ² = 1)
        for (var lag = 1; lag < n; lag++)
        {
            double acf = 0;
            for (var i = 0; i + lag < n; i++) acf += v[i] * v[i + lag];
            total += 2 * acf * Math.Sin(2 * Math.PI * w * lag) / (Math.PI * lag);
        }
        return total;
    }
}
