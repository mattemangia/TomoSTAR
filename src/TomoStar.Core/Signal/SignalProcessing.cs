// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using System.Numerics;

namespace TomoStar.Core.Signal;

/// <summary>Elementary operations on sampled signals.</summary>
public static class SignalProcessing
{
    /// <summary>Removes the least-squares straight line.</summary>
    public static void Detrend(Span<float> x)
    {
        var n = x.Length;
        if (n < 2) return;
        double sx = 0, sy = 0, sxx = 0, sxy = 0;
        for (var i = 0; i < n; i++)
        {
            sx += i; sy += x[i]; sxx += (double)i * i; sxy += i * (double)x[i];
        }
        var den = n * sxx - sx * sx;
        var slope = den != 0 ? (n * sxy - sx * sy) / den : 0;
        var intercept = (sy - slope * sx) / n;
        for (var i = 0; i < n; i++) x[i] = (float)(x[i] - (intercept + slope * i));
    }

    public static void Demean(Span<float> x)
    {
        if (x.Length == 0) return;
        double s = 0;
        foreach (var v in x) s += v;
        var m = (float)(s / x.Length);
        for (var i = 0; i < x.Length; i++) x[i] -= m;
    }

    /// <summary>Cosine (Tukey) taper of the given fraction at each end.</summary>
    public static void Taper(Span<float> x, double fraction = 0.05)
    {
        var n = x.Length;
        var m = (int)(Math.Clamp(fraction, 0, 0.5) * n);
        for (var i = 0; i < m; i++)
        {
            var w = (float)(0.5 * (1 - Math.Cos(Math.PI * i / m)));
            x[i] *= w;
            x[n - 1 - i] *= w;
        }
    }

    public static float MaxAbs(ReadOnlySpan<float> x)
    {
        float m = 0;
        foreach (var v in x) m = Math.Max(m, Math.Abs(v));
        return m;
    }

    public static double Rms(ReadOnlySpan<float> x)
    {
        if (x.Length == 0) return 0;
        double s = 0;
        foreach (var v in x) s += (double)v * v;
        return Math.Sqrt(s / x.Length);
    }

    /// <summary>
    /// The standard processing before picking or display: demean, detrend, taper, band-pass
    /// (zero-phase by default; causal for onset picking, see <see cref="Filters"/>).
    /// </summary>
    public static float[] Prepare(ReadOnlySpan<float> data, double sampleRate, double lowHz, double highHz, int order = 4, bool zeroPhase = true)
    {
        var x = data.ToArray();
        Demean(x);
        Detrend(x);
        Taper(x, 0.02);
        if (lowHz > 0 || highHz > 0) Butterworth.Filter(x, sampleRate, lowHz, highHz, order, zeroPhase);
        return x;
    }
}

/// <summary>
/// Butterworth filters as cascades of second-order sections designed by the bilinear transform
/// with frequency pre-warping, applied forward and backward for zero phase (the effective order
/// doubles and the −3 dB points become −6 dB). Standard design: Oppenheim &amp; Schafer,
/// Discrete-Time Signal Processing (3rd ed., 2010), §7.1.
/// </summary>
public static class Butterworth
{
    /// <summary>A second-order section, direct form II transposed, a0 normalised to 1.</summary>
    internal readonly record struct Section(double B0, double B1, double B2, double A1, double A2);

    /// <param name="lowHz">High-pass corner, 0 to skip.</param>
    /// <param name="highHz">Low-pass corner, 0 (or ≥ Nyquist) to skip.</param>
    public static void BandpassZeroPhase(Span<float> x, double sampleRate, double lowHz, double highHz, int order = 4) =>
        Filter(x, sampleRate, lowHz, highHz, order, zeroPhase: true);

    /// <summary>
    /// High-pass at <paramref name="lowHz"/> and/or low-pass at <paramref name="highHz"/> (0 skips
    /// either), causal (forward only, the order given) or zero-phase (forward and backward).
    /// </summary>
    public static void Filter(Span<float> x, double sampleRate, double lowHz, double highHz, int order, bool zeroPhase)
    {
        var nyq = 0.5 * sampleRate;
        var sections = new List<Section>();
        if (highHz > 0 && highHz < 0.99 * nyq) sections.AddRange(Design(order, highHz, sampleRate, lowPass: true));
        if (lowHz > 0 && lowHz < 0.99 * nyq) sections.AddRange(Design(order, lowHz, sampleRate, lowPass: false));
        if (sections.Count == 0) return;
        foreach (var s in sections) Run(x, s, reverse: false);
        if (zeroPhase) foreach (var s in sections) Run(x, s, reverse: true);
    }

    private static IEnumerable<Section> Design(int order, double fc, double fs, bool lowPass)
    {
        order = Math.Max(2, order + (order & 1)); // even order: whole biquads
        var k = Math.Tan(Math.PI * fc / fs);       // pre-warped analogue corner
        var k2 = k * k;
        for (var i = 0; i < order / 2; i++)
        {
            // Pole pair angle of the analogue prototype.
            var theta = Math.PI * (2.0 * i + 1) / (2.0 * order);
            var q = 1.0 / (2.0 * Math.Sin(theta));
            var norm = 1.0 / (1 + k / q + k2);
            if (lowPass)
            {
                var b0 = k2 * norm;
                yield return new Section(b0, 2 * b0, b0, 2 * (k2 - 1) * norm, (1 - k / q + k2) * norm);
            }
            else
            {
                yield return new Section(norm, -2 * norm, norm, 2 * (k2 - 1) * norm, (1 - k / q + k2) * norm);
            }
        }
    }

    internal static void Run(Span<float> x, Section f, bool reverse)
    {
        double z1 = 0, z2 = 0;
        var n = x.Length;
        for (var t = 0; t < n; t++)
        {
            var i = reverse ? n - 1 - t : t;
            double input = x[i];
            // Direct form II transposed.
            var y = f.B0 * input + z1;
            z1 = f.B1 * input - f.A1 * y + z2;
            z2 = f.B2 * input - f.A2 * y;
            x[i] = (float)y;
        }
    }
}

/// <summary>Radix-2 complex FFT (Cooley-Tukey), in place.</summary>
public static class Fft
{
    public static int NextPow2(int n)
    {
        var p = 1;
        while (p < n) p <<= 1;
        return p;
    }

    public static void Transform(Complex[] a, bool inverse = false)
    {
        var n = a.Length;
        if ((n & (n - 1)) != 0) throw new ArgumentException("FFT length must be a power of two.");
        for (int i = 1, j = 0; i < n; i++)
        {
            var bit = n >> 1;
            for (; (j & bit) != 0; bit >>= 1) j ^= bit;
            j ^= bit;
            if (i < j) (a[i], a[j]) = (a[j], a[i]);
        }
        for (var len = 2; len <= n; len <<= 1)
        {
            var ang = 2 * Math.PI / len * (inverse ? 1 : -1);
            var wl = new Complex(Math.Cos(ang), Math.Sin(ang));
            for (var i = 0; i < n; i += len)
            {
                var w = Complex.One;
                for (var j = 0; j < len / 2; j++)
                {
                    var u = a[i + j];
                    var v = a[i + j + len / 2] * w;
                    a[i + j] = u + v;
                    a[i + j + len / 2] = u - v;
                    w *= wl;
                }
            }
        }
        if (inverse)
            for (var i = 0; i < n; i++) a[i] /= n;
    }

    /// <summary>
    /// One-sided amplitude spectrum of a real window (tapered, zero-padded), in units of the
    /// continuous Fourier amplitude (samples·s). <paramref name="taperFraction"/> is the Tukey
    /// fraction of the window that is tapered, half at each end: 1 is a Hann window (spectral
    /// analysis of a stationary signal), 0.1 keeps the window flat except at its edges, which is
    /// what a window starting at a P onset needs, since a Hann window would all but erase the onset.
    /// Returns frequencies (Hz) and amplitudes.
    /// </summary>
    public static (double[] Freq, double[] Amp) AmplitudeSpectrum(ReadOnlySpan<float> x, double sampleRate, int minLength = 0, double taperFraction = 1)
    {
        var n = NextPow2(Math.Max(x.Length, minLength));
        var a = new Complex[n];
        var ramp = Math.Clamp(taperFraction, 0, 1) * 0.5 * (x.Length - 1); // samples tapered at each end
        for (var i = 0; i < x.Length; i++)
        {
            var edge = Math.Min(i, x.Length - 1 - i);
            var w = ramp <= 0 || edge >= ramp ? 1 : 0.5 * (1 - Math.Cos(Math.PI * edge / ramp));
            a[i] = new Complex(x[i] * w, 0);
        }
        Transform(a);
        var half = n / 2;
        var f = new double[half];
        var amp = new double[half];
        for (var k = 0; k < half; k++)
        {
            f[k] = k * sampleRate / n;
            amp[k] = a[k].Magnitude / sampleRate; // ≈ continuous Fourier amplitude, s·units
        }
        return (f, amp);
    }
}
