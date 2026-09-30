// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using System.Numerics;
using System.Text.Json.Serialization;

namespace TomoStar.Core.Signal;

public enum GroundMotion
{
    Displacement,
    Velocity,
    Acceleration
}

/// <summary>A poles-and-zeros stage of a StationXML response (analogue, Laplace transform).</summary>
public sealed class PolesZerosStage
{
    /// <summary>"LAPLACE (RADIANS/SECOND)" or "LAPLACE (HERTZ)" (digital stages are not stored).</summary>
    public string TransferFunction { get; set; } = "LAPLACE (RADIANS/SECOND)";

    public double NormalizationFactor { get; set; } = 1;
    public double NormalizationFrequency { get; set; } = 1;

    /// <summary>Zeros and poles as [real, imaginary] pairs.</summary>
    public List<double[]> Zeros { get; set; } = [];

    public List<double[]> Poles { get; set; } = [];

    [JsonIgnore] public bool InHertz => TransferFunction.Contains("HERTZ", StringComparison.OrdinalIgnoreCase);

    /// <summary>A0 · Π(s − zᵢ) / Π(s − pⱼ) at frequency f, with s = 2πif (or if for a Hertz stage).</summary>
    public Complex Evaluate(double f)
    {
        var s = InHertz ? new Complex(0, f) : new Complex(0, 2 * Math.PI * f);
        var h = new Complex(NormalizationFactor, 0);
        foreach (var z in Zeros) h *= s - new Complex(z[0], z[1]);
        foreach (var p in Poles) h /= s - new Complex(p[0], p[1]);
        return h;
    }
}

/// <summary>
/// The response of a channel from StationXML: the overall sensitivity (counts per input unit at a
/// frequency) and the shape of the analogue poles-and-zeros stages. The complex response is
/// R(f) = S · H(f) / |H(f_S)|, H being the product of the poles-and-zeros stages: the amplitude is
/// exactly the sensitivity at its frequency and the shape and phase are those of the analogue
/// stages. Digital FIR stages are not modelled; their passband is flat to within a few percent up
/// to about 80% of the Nyquist frequency, where the pre-filter of the deconvolution cuts anyway.
/// </summary>
public sealed class InstrumentResponse
{
    public string InputUnits { get; set; } = "";
    public double Sensitivity { get; set; } = 1;
    public double SensitivityFrequency { get; set; } = 1;
    public List<PolesZerosStage> Stages { get; set; } = [];

    [JsonIgnore] public bool HasShape => Stages.Count > 0;

    /// <summary>0 for displacement (m), 1 for velocity (m/s), 2 for acceleration (m/s²), −1 when unknown.</summary>
    [JsonIgnore]
    public int InputPower => UnitPower(InputUnits);

    public static int UnitPower(string units)
    {
        var u = units.Trim().ToUpperInvariant().Replace(" ", "");
        return u switch
        {
            "M" or "NM" or "MM" or "CM" => 0,
            "M/S" or "NM/S" or "MM/S" or "CM/S" => 1,
            "M/S**2" or "M/S/S" or "M/S2" or "M/S^2" or "NM/S**2" or "NM/S/S" or "MM/S**2" or "CM/S**2" or "G" => 2,
            _ => -1
        };
    }

    /// <summary>Counts per input unit, complex, at frequency f.</summary>
    public Complex Evaluate(double f)
    {
        if (!HasShape) return new Complex(Sensitivity, 0);
        var h = Complex.One;
        foreach (var s in Stages) h *= s.Evaluate(f);
        var h0 = Complex.One;
        foreach (var s in Stages) h0 *= s.Evaluate(SensitivityFrequency);
        var m = h0.Magnitude;
        return m > 0 ? Sensitivity * h / m : new Complex(Sensitivity, 0);
    }
}

/// <summary>
/// Removal of the instrument response by spectral division, as usual in seismology (e.g. ObsPy's
/// remove_response): the record is demeaned and tapered, transformed, divided by the complex
/// response with a water level (|R| is not allowed below max|R|·10^(−level/20), which bounds the
/// amplification of noise where the instrument is insensitive), converted to the requested ground
/// motion by powers of iω, band-limited by a cosine pre-filter f1 &lt; f2 … f3 &lt; f4, and transformed back.
/// </summary>
public static class ResponseRemoval
{
    public static float[] Remove(ReadOnlySpan<float> data, double sampleRate, InstrumentResponse response, GroundMotion output,
        double f1, double f2, double f3, double f4, double waterLevelDb = 60)
    {
        var n = data.Length;
        if (n < 8) return data.ToArray();
        var x = data.ToArray();
        SignalProcessing.Demean(x);
        SignalProcessing.Detrend(x);
        SignalProcessing.Taper(x, 0.05);
        var nfft = Fft.NextPow2(2 * n);
        var a = new Complex[nfft];
        for (var i = 0; i < n; i++) a[i] = new Complex(x[i], 0);
        Fft.Transform(a);

        var half = nfft / 2;
        var r = new Complex[half + 1];
        double rmax = 0;
        for (var k = 1; k <= half; k++)
        {
            r[k] = response.Evaluate(k * sampleRate / nfft);
            rmax = Math.Max(rmax, r[k].Magnitude);
        }
        var floor = rmax * Math.Pow(10, -waterLevelDb / 20);
        var inPow = response.InputPower < 0 ? 1 : response.InputPower;
        var outPow = output switch { GroundMotion.Displacement => 0, GroundMotion.Velocity => 1, _ => 2 };
        a[0] = Complex.Zero;
        for (var k = 1; k <= half; k++)
        {
            var f = k * sampleRate / nfft;
            var w = Taper(f, f1, f2, f3, f4);
            Complex v;
            if (w == 0) v = Complex.Zero;
            else
            {
                var rk = r[k];
                if (rk.Magnitude < floor) rk = rk.Magnitude > 0 ? rk * (floor / rk.Magnitude) : new Complex(floor, 0);
                v = a[k] / rk * w;
                // From the input quantity to the requested one: ×(iω) per derivative.
                var iw = new Complex(0, 2 * Math.PI * f);
                for (var p = inPow; p < outPow; p++) v *= iw;
                for (var p = outPow; p < inPow; p++) v /= iw;
            }
            a[k] = v;
            if (k < half) a[nfft - k] = Complex.Conjugate(v);
        }
        a[half] = new Complex(a[half].Real, 0);
        Fft.Transform(a, inverse: true);
        var y = new float[n];
        for (var i = 0; i < n; i++) y[i] = (float)a[i].Real;
        return y;
    }

    /// <summary>Cosine taper: 0 below f1 and above f4, 1 between f2 and f3.</summary>
    public static double Taper(double f, double f1, double f2, double f3, double f4)
    {
        if (f <= f1 || f >= f4) return 0;
        if (f < f2) return 0.5 * (1 - Math.Cos(Math.PI * (f - f1) / (f2 - f1)));
        if (f > f3) return 0.5 * (1 + Math.Cos(Math.PI * (f - f3) / (f4 - f3)));
        return 1;
    }

    /// <summary>Default pre-filter for a record: from twice its lowest resolvable frequency to 80% of Nyquist.</summary>
    public static (double F1, double F2, double F3, double F4) DefaultPreFilter(double sampleRate, double durationSeconds)
    {
        var lowest = 1 / Math.Max(1, durationSeconds);
        var nyq = 0.5 * sampleRate;
        return (Math.Max(lowest, 0.005), Math.Max(2 * lowest, 0.01), 0.8 * nyq, 0.9 * nyq);
    }

    public static string Units(GroundMotion g) => g switch { GroundMotion.Displacement => "m", GroundMotion.Velocity => "m/s", _ => "m/s^2" };
}
