// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using TomoStar.Core.IO;

namespace TomoStar.Core.Signal;

/// <summary>The quality of one channel of one event, with the reasons it may be unusable.</summary>
public sealed record TraceQualityReport(
    string Nslc, int Samples, double DurationSeconds, int Segments, double GapSeconds, double OverlapSeconds,
    double MaxAbs, double Rms, bool Flat, bool Clipped, int Spikes, double SnrAtP, string Issues)
{
    public bool Usable => !Flat && !Clipped;
}

/// <summary>
/// Automatic checks of a waveform before it is picked or inverted:
/// <list type="bullet">
/// <item>segments, gaps and overlaps between the records of the channel;</item>
/// <item>flat (dead) channel: no variation, or one value for more than half the record;</item>
/// <item>clipping: the largest absolute value held for three or more consecutive samples, at least
///   twice, as a saturated digitiser or sensor produces flat-topped peaks;</item>
/// <item>spikes: isolated samples far from both neighbours (more than 25 median absolute deviations
///   of the sample-to-sample differences), typical of telemetry errors;</item>
/// <item>signal-to-noise ratio at the P arrival: RMS in 2 s after P over RMS in 5 s before it.</item>
/// </list>
/// </summary>
public static class TraceQuality
{
    public static TraceQualityReport Assess(IReadOnlyList<Trace> segments, DateTime? pArrival = null)
    {
        var main = segments.OrderByDescending(s => s.Data.Length).First();
        var x = main.Data;
        var n = x.Length;
        var issues = new List<string>();

        // Gaps and overlaps between consecutive segments.
        double gap = 0, overlap = 0;
        var ordered = segments.OrderBy(s => s.StartTime).ToList();
        for (var i = 1; i < ordered.Count; i++)
        {
            var d = (ordered[i].StartTime - ordered[i - 1].EndTime).TotalSeconds - 1 / ordered[i - 1].SampleRate;
            if (d > 0.5 / ordered[i].SampleRate) gap += d;
            else if (d < -0.5 / ordered[i].SampleRate) overlap += -d;
        }
        if (gap > 0) issues.Add($"gaps {gap:0.##} s");
        if (overlap > 0) issues.Add($"overlaps {overlap:0.##} s");

        double mean = 0;
        foreach (var v in x) mean += v;
        mean = n > 0 ? mean / n : 0;
        double ss = 0;
        float maxAbs = 0;
        foreach (var v in x) { ss += (v - mean) * (v - mean); maxAbs = Math.Max(maxAbs, Math.Abs(v)); }
        var rms = n > 0 ? Math.Sqrt(ss / n) : 0;

        // Longest run of one value.
        var longest = 0;
        for (int i = 0, run = 1; i < n; i++)
        {
            run = i > 0 && x[i] == x[i - 1] ? run + 1 : 1;
            longest = Math.Max(longest, run);
        }
        var flat = n < 2 || rms <= 1e-9 * Math.Max(1, Math.Abs(mean)) || longest > n / 2;
        if (flat) issues.Add("flat");

        // Flat tops at the maximum.
        var plateaus = 0;
        if (!flat && maxAbs > 0)
        {
            for (var i = 0; i < n;)
            {
                if (Math.Abs(x[i]) < maxAbs) { i++; continue; }
                var j = i;
                while (j + 1 < n && x[j + 1] == x[i]) j++;
                if (j - i + 1 >= 3) plateaus++;
                i = j + 1;
            }
        }
        var clipped = plateaus >= 2;
        if (clipped) issues.Add("clipped");

        // Isolated spikes.
        var spikes = 0;
        if (!flat && n > 10)
        {
            var diffs = new float[n - 1];
            for (var i = 1; i < n; i++) diffs[i - 1] = Math.Abs(x[i] - x[i - 1]);
            var sorted = (float[])diffs.Clone();
            Array.Sort(sorted);
            var mad = Math.Max(1e-12, sorted[sorted.Length / 2]);
            for (var i = 1; i < n - 1; i++)
            {
                double a = x[i] - x[i - 1], b = x[i] - x[i + 1];
                if (Math.Sign(a) == Math.Sign(b) && Math.Min(Math.Abs(a), Math.Abs(b)) > 25 * mad
                    && Math.Abs(x[i + 1] - x[i - 1]) < 5 * mad) spikes++;
            }
        }
        if (spikes > 0) issues.Add($"{spikes} spike{(spikes == 1 ? "" : "s")}");

        var snr = double.NaN;
        if (pArrival is { } tp && !flat)
        {
            var ip = main.IndexOf(tp);
            var fs = main.SampleRate;
            int a0 = ip - (int)(5 * fs), a1 = ip - (int)(0.2 * fs), b1 = ip + (int)(2 * fs);
            if (a0 >= 0 && b1 < n && a1 > a0)
            {
                snr = Rms(x, ip, b1) / Math.Max(1e-30, Rms(x, a0, a1));
                if (snr < 2) issues.Add($"low SNR at P ({snr:0.0})");
            }
        }
        return new TraceQualityReport(main.Nslc, n, n / main.SampleRate, segments.Count, gap, overlap, maxAbs, rms, flat, clipped, spikes, snr,
            issues.Count == 0 ? "ok" : string.Join("; ", issues));
    }

    private static double Rms(float[] x, int from, int to)
    {
        double m = 0, s = 0;
        for (var i = from; i < to; i++) m += x[i];
        m /= Math.Max(1, to - from);
        for (var i = from; i < to; i++) s += (x[i] - m) * (x[i] - m);
        return Math.Sqrt(s / Math.Max(1, to - from));
    }
}
