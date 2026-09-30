using TomoStar.Core.IO;
using TomoStar.Core.Model;

namespace TomoStar.Core.Signal;

public sealed class AutoPickerSettings
{
    public double FilterLowHz { get; set; } = 2;
    public double FilterHighHz { get; set; } = 15;
    public double StaSeconds { get; set; } = 0.5;
    public double LtaSeconds { get; set; } = 5;
    public double TriggerOn { get; set; } = 3.5;

    /// <summary>Half-width of the search window around a predicted arrival, s (scaled up with travel time).</summary>
    public double WindowSeconds { get; set; } = 3;

    /// <summary>Extra window per second of predicted travel time (model error grows with distance).</summary>
    public double WindowPerSecond { get; set; } = 0.08;

    public bool PickS { get; set; } = true;

    /// <summary>Picks with a signal-to-noise ratio below this are not kept.</summary>
    public double MinSnr { get; set; } = 3;

    /// <summary>
    /// Filter forward and backward (zero phase). Off by default: a causal filter puts no energy
    /// before the onset, so the AIC onset is not pulled early by the filter's symmetric response.
    /// </summary>
    public bool ZeroPhaseFilter { get; set; }
}

/// <summary>An automatic pick with the evidence behind it, for the diagnostics view.</summary>
public sealed record AutoPickResult(
    Phase Phase, DateTime Time, double Snr, int Quality, double Uncertainty,
    string Channel, int TriggerIndex, int OnsetIndex);

/// <summary>
/// Automatic P and S picker.
///
/// 1. A recursive STA/LTA (Allen 1978, BSSA 68(5), 1521-1532; recursive form after Withers et al.
///    1998, BSSA 88(1), 95-106) finds the first trigger inside a window around the time the current
///    model predicts; the window keeps the picker from locking onto noise bursts or later phases.
/// 2. The onset is refined with the AIC picker of Maeda (1985, Zisin 38, 365-379): on the segment
///    around the trigger, AIC(k) = k·log var(x[0..k]) + (N−k−1)·log var(x[k+1..N]) is minimum where
///    the series splits best into "noise" and "signal".
/// 3. The quality class follows the signal-to-noise ratio (HYPO71-style weights 0-4) and gives the
///    pick uncertainty used as its inversion weight.
///
/// P is picked on the vertical; S on the horizontals (their summed energy) when present, otherwise
/// on the vertical.
/// </summary>
public static class AutoPicker
{
    public static List<AutoPickResult> PickStation(
        IReadOnlyList<Trace> traces, DateTime? predictedP, DateTime? predictedS, AutoPickerSettings s)
    {
        var results = new List<AutoPickResult>();
        var z = traces.FirstOrDefault(t => t.IsVertical) ?? traces.FirstOrDefault();
        if (z == null || z.SampleRate <= 0) return results;
        var horizontals = traces.Where(t => !t.IsVertical && Math.Abs(t.SampleRate - z.SampleRate) < 1e-6).ToList();

        var zf = SignalProcessing.Prepare(z.Data, z.SampleRate, s.FilterLowHz, s.FilterHighHz, 4, s.ZeroPhaseFilter);
        var p = PickOne(zf, z, predictedP, Phase.P, s);
        if (p != null) results.Add(p);

        if (s.PickS)
        {
            float[] energy;
            Trace reference;
            if (horizontals.Count > 0)
            {
                reference = horizontals[0];
                var start = horizontals.Max(h => h.StartTime);
                var len = horizontals.Min(h => h.Data.Length - h.IndexOf(start));
                if (len < 10) return results;
                energy = new float[len];
                foreach (var h in horizontals)
                {
                    var f = SignalProcessing.Prepare(h.Data, h.SampleRate, s.FilterLowHz, s.FilterHighHz, 4, s.ZeroPhaseFilter);
                    var o = h.IndexOf(start);
                    for (var i = 0; i < len; i++) energy[i] += f[o + i] * f[o + i];
                }
                for (var i = 0; i < len; i++) energy[i] = MathF.Sqrt(energy[i]);
                reference = reference.Clone(energy);
                reference.StartTime = start;
            }
            else
            {
                reference = z;
                energy = zf;
            }
            // S must come after P, and the P wave and its coda also reach the horizontals: the
            // search starts halfway between the P pick and the predicted S (at least 0.5 s after P),
            // otherwise the STA/LTA still raised by the P wave triggers at once and the "S" lands
            // on the P coda (5 s early on the Amatrice 2016 data).
            DateTime? notBefore = null;
            if (p != null)
            {
                var gap = predictedS is { } ps && predictedP is { } pp ? 0.5 * (ps - pp).TotalSeconds : 0;
                notBefore = p.Time.AddSeconds(Math.Max(0.5, gap));
            }
            var sPick = PickOne(energy, reference, predictedS, Phase.S, s, notBefore);
            if (sPick != null) results.Add(sPick);
        }
        return results;
    }

    private static AutoPickResult? PickOne(float[] x, Trace trace, DateTime? predicted, Phase phase, AutoPickerSettings s, DateTime? notBefore = null)
    {
        var fs = trace.SampleRate;
        int lo, hi;
        if (predicted is { } tp)
        {
            var travel = Math.Max(0, (tp - trace.StartTime).TotalSeconds);
            var half = s.WindowSeconds + s.WindowPerSecond * travel;
            lo = trace.IndexOf(tp.AddSeconds(-half));
            hi = trace.IndexOf(tp.AddSeconds(half));
        }
        else
        {
            lo = 0;
            hi = x.Length - 1;
        }
        if (notBefore is { } nb) lo = Math.Max(lo, trace.IndexOf(nb));
        lo = Math.Clamp(lo, 0, x.Length - 1);
        hi = Math.Clamp(hi, 0, x.Length - 1);
        if (hi - lo < (int)(2 * fs)) return null;

        var cf = StaLta(x, fs, s.StaSeconds, s.LtaSeconds);
        var warm = (int)(s.LtaSeconds * fs);
        // The first upward crossing of the threshold: a window that opens while an earlier arrival
        // still holds the ratio up must not trigger on that arrival.
        var trigger = -1;
        for (var i = Math.Max(Math.Max(lo, warm), 1); i <= hi; i++)
            if (cf[i] >= s.TriggerOn && cf[i - 1] < s.TriggerOn) { trigger = i; break; }
        if (trigger < 0) return null;

        // AIC on a segment around the trigger: onsets lie before the STA/LTA trigger.
        var a = Math.Max(0, trigger - (int)(Math.Max(1.0, 2 * s.StaSeconds + 1) * fs));
        var b = Math.Min(x.Length - 1, trigger + (int)(0.5 * fs));
        var onset = AicOnset(x, a, b);
        if (onset < 0) onset = trigger;

        var snr = Snr(x, onset, fs);
        if (snr < s.MinSnr) return null;
        var quality = snr >= 20 ? 0 : snr >= 10 ? 1 : snr >= 5 ? 2 : 3;
        // Uncertainty grows as the onset gets emergent; a floor of two samples.
        var uncertainty = Math.Max(2 / fs, quality switch { 0 => 0.03, 1 => 0.06, 2 => 0.12, _ => 0.25 } * (phase == Phase.S ? 1.5 : 1));
        return new AutoPickResult(phase, trace.TimeOf(onset), snr, quality, uncertainty, trace.Channel, trigger, onset);
    }

    /// <summary>Recursive STA/LTA characteristic function of the squared signal.</summary>
    public static float[] StaLta(ReadOnlySpan<float> x, double fs, double sta, double lta)
    {
        var cf = new float[x.Length];
        var a = Math.Min(1, 1 / (sta * fs));
        var b = Math.Min(1, 1 / (lta * fs));
        double s = 0, l = 1e-30;
        // Seed the LTA with the first second's mean energy so the ratio is not huge at the start.
        var seed = Math.Min(x.Length, (int)Math.Max(1, fs));
        for (var i = 0; i < seed; i++) l += (double)x[i] * x[i] / seed;
        s = l;
        for (var i = 0; i < x.Length; i++)
        {
            var e = (double)x[i] * x[i];
            s = a * e + (1 - a) * s;
            l = b * e + (1 - b) * l;
            cf[i] = (float)(l > 1e-30 ? s / l : 0);
        }
        return cf;
    }

    /// <summary>Maeda (1985) AIC function on x[a..b]; returns it for display.</summary>
    public static double[] Aic(ReadOnlySpan<float> x, int a, int b)
    {
        var n = b - a + 1;
        var aic = new double[n];
        if (n < 4) return aic;
        // Prefix sums of x and x² make every variance O(1).
        var s1 = new double[n + 1];
        var s2 = new double[n + 1];
        for (var i = 0; i < n; i++)
        {
            s1[i + 1] = s1[i] + x[a + i];
            s2[i + 1] = s2[i] + (double)x[a + i] * x[a + i];
        }
        for (var k = 1; k < n - 1; k++)
        {
            double n1 = k, n2 = n - k;
            var v1 = Math.Max(1e-30, s2[k] / n1 - Math.Pow(s1[k] / n1, 2));
            var v2 = Math.Max(1e-30, (s2[n] - s2[k]) / n2 - Math.Pow((s1[n] - s1[k]) / n2, 2));
            aic[k] = k * Math.Log(v1) + (n - k - 1) * Math.Log(v2);
        }
        aic[0] = aic[1];
        aic[n - 1] = aic[n - 2];
        return aic;
    }

    public static int AicOnset(ReadOnlySpan<float> x, int a, int b)
    {
        var aic = Aic(x, a, b);
        if (aic.Length < 4) return -1;
        var best = 1;
        for (var k = 2; k < aic.Length - 1; k++)
            if (aic[k] < aic[best]) best = k;
        return a + best;
    }

    /// <summary>Ratio of the maximum amplitude in 1 s after the onset to the RMS of 2 s before it.</summary>
    public static double Snr(ReadOnlySpan<float> x, int onset, double fs)
    {
        var n0 = Math.Max(0, onset - (int)(2 * fs));
        var n1 = Math.Min(x.Length, onset + (int)fs);
        if (onset - n0 < 5 || n1 - onset < 5) return 0;
        var noise = SignalProcessing.Rms(x[n0..onset]);
        var signal = SignalProcessing.MaxAbs(x[onset..n1]);
        return noise > 0 ? signal / noise : 0;
    }

    /// <summary>Converts a result to an automatic pick of the catalogue.</summary>
    public static PickRecord ToPick(AutoPickResult r, Trace trace) => new()
    {
        StationId = trace.StationId, Channel = r.Channel, Phase = r.Phase, Time = r.Time, Sigma = r.Uncertainty,
        Quality = r.Quality, Origin = PickOrigin.Automatic, Snr = r.Snr
    };
}
