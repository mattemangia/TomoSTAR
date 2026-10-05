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

    /// <summary>Frequency at which the stage gain is declared, Hz (NaN: the normalisation frequency).</summary>
    public double GainFrequency { get; set; } = double.NaN;

    /// <summary>
    /// The stage scaled so that the stage gain holds at its own frequency when that differs from the
    /// normalisation frequency: A0 H(f) · |A0 H(f_norm)| / |A0 H(f_gain)| (the convention of evalresp).
    /// </summary>
    public Complex EvaluateWithGainFrequency(double f)
    {
        var h = Evaluate(f);
        if (!double.IsFinite(GainFrequency) || Math.Abs(GainFrequency - NormalizationFrequency) < 1e-9) return h;
        var atGain = Evaluate(GainFrequency).Magnitude;
        return atGain > 0 ? h * (Evaluate(NormalizationFrequency).Magnitude / atGain) : h;
    }

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
/// A digital stage of a StationXML response (FIR or coefficients, numerator only), normalised to 1 at
/// zero frequency as evalresp normalises FIR filters. Symmetric (linear-phase) filters are a pure
/// amplitude; asymmetric ones keep their phase less the digitiser's time-stamp correction (IV.GIGS in
/// 2016 has an asymmetric stage without correction: 0.03 s of delay, three samples at 100 Hz).
/// </summary>
public sealed class DigitalStage
{
    /// <summary>Sampling rate at the input of the stage, Hz.</summary>
    public double InputSampleRate { get; set; }

    /// <summary>The full set of numerator coefficients (symmetric halves already expanded).</summary>
    public double[] Coefficients { get; set; } = [];

    /// <summary>Delay of the stage and the correction applied to the time stamps for it, s.</summary>
    public double Delay { get; set; }

    public double Correction { get; set; }

    /// <summary>
    /// The stage's complex factor at f, by the convention of evalresp: a symmetric filter is taken as
    /// zero phase, Σ bₖ cos(2πf (k − (N−1)/2) / fs) / Σ bₖ (its delay is the one digitisers correct);
    /// an asymmetric one keeps its phase, Σ bₖ e^(−2πi f k / fs) / Σ bₖ, times e^(2πi f c) for the
    /// correction c applied to the time stamps.
    /// </summary>
    public Complex Evaluate(double f)
    {
        var n = Coefficients.Length;
        if (n < 2 || InputSampleRate <= 0) return Complex.One;
        var sum = Coefficients.Sum();
        if (Math.Abs(sum) == 0) return Complex.One;
        var w = 2 * Math.PI * f / InputSampleRate;
        if (Symmetric)
        {
            double re = 0;
            for (var k = 0; k < n; k++) re += Coefficients[k] * Math.Cos(w * (k - 0.5 * (n - 1)));
            return new Complex(re / Math.Abs(sum), 0);
        }
        var h = Complex.Zero;
        for (var k = 0; k < n; k++) h += Coefficients[k] * Complex.FromPolarCoordinates(1, -w * k);
        return h / Math.Abs(sum) * Complex.FromPolarCoordinates(1, 2 * Math.PI * f * Correction);
    }

    /// <summary>The coefficients read the same backwards (to a part in 10⁶ of the largest).</summary>
    public bool Symmetric
    {
        get
        {
            var n = Coefficients.Length;
            var tol = 1e-6 * Coefficients.Max(Math.Abs);
            for (var k = 0; k < n / 2; k++)
                if (Math.Abs(Coefficients[k] - Coefficients[n - 1 - k]) > tol) return false;
            return true;
        }
    }
}

/// <summary>
/// The response of a channel from StationXML: the overall sensitivity (counts per input unit at a
/// frequency), the product of the stage gains and the analogue poles-and-zeros stages. When every
/// stage declares its gain, the complex response is R(f) = G · Π A0ₖ Hₖ(f), G being the product of
/// the stage gains and A0ₖ Hₖ the normalised poles-and-zeros stages (each rescaled so that its gain
/// holds at its declared gain frequency), as evalresp (and so ObsPy) computes it; otherwise R(f) = S · H(f) / |H(f_S)|, the shape of the analogue stages scaled to the
/// sensitivity at its frequency. The two agree when the metadata are consistent; where they are
/// not, the stage description is the one to trust: on the INGV metadata of 2012 for IV.SAP2 the
/// sensitivity of a 1 Hz geophone is declared at 0.2 Hz but holds at 1 Hz, which the second form
/// turned into a response 25 times too high. The amplitude of the digital stages multiplies both
/// forms: most decimation filters are flat in their passband, but some are not (the binomial
/// filters of the IV digitisers of 2008-2009 halve the amplitude at 20 Hz on a 100 Hz channel).
/// </summary>
public sealed class InstrumentResponse
{
    public string InputUnits { get; set; } = "";
    public double Sensitivity { get; set; } = 1;
    public double SensitivityFrequency { get; set; } = 1;

    /// <summary>Product of the gains of all the stages, NaN when a stage declares none.</summary>
    public double StageGain { get; set; } = double.NaN;

    public List<PolesZerosStage> Stages { get; set; } = [];

    /// <summary>Digital stages with coefficients (FIR and numerator-only coefficient stages).</summary>
    public List<DigitalStage> DigitalStages { get; set; } = [];

    private Complex Digital(double f)
    {
        var a = Complex.One;
        foreach (var d in DigitalStages) a *= d.Evaluate(f);
        return a;
    }

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
        if (!HasShape) return Sensitivity * Digital(f) / Math.Max(1e-12, Digital(SensitivityFrequency).Magnitude);
        if (double.IsFinite(StageGain) && StageGain > 0)
        {
            var g = Complex.One;
            foreach (var s in Stages) g *= s.EvaluateWithGainFrequency(f);
            return StageGain * g * Digital(f);
        }
        var h = Complex.One;
        foreach (var s in Stages) h *= s.Evaluate(f);
        h *= Digital(f);
        var h0 = Complex.One;
        foreach (var s in Stages) h0 *= s.Evaluate(SensitivityFrequency);
        var m = h0.Magnitude * Digital(SensitivityFrequency).Magnitude;
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
