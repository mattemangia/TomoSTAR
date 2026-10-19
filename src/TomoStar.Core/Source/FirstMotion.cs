// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using TomoStar.Core.IO;

namespace TomoStar.Core.Source;

/// <summary>One P first-motion reading at a station, with the ray's direction at the source.</summary>
public sealed record PolarityObservation(string StationId, double AzimuthDeg, double TakeoffDeg, int Polarity, double Weight = 1)
{
    /// <summary>Unit direction of the ray leaving the source (north, east, down); takeoff from the downward vertical.</summary>
    public (double N, double E, double D) Direction
    {
        get
        {
            double i = TakeoffDeg * Math.PI / 180, a = AzimuthDeg * Math.PI / 180;
            return (Math.Sin(i) * Math.Cos(a), Math.Sin(i) * Math.Sin(a), Math.Cos(i));
        }
    }
}

public sealed class FirstMotionSettings
{
    /// <summary>Grid step in strike, dip and rake, degrees.</summary>
    public double GridStepDeg { get; set; } = 5;

    /// <summary>Fraction of polarities expected to be wrong; sets how many misfits an acceptable mechanism may have.</summary>
    public double ExpectedErrorFraction { get; set; } = 0.1;

    public int MinPolarities { get; set; } = 8;
}

/// <summary>A first-motion focal mechanism with its uncertainty and quality.</summary>
public sealed class FirstMotionResult
{
    public required NodalPlane Plane1 { get; init; }
    public required NodalPlane Plane2 { get; init; }
    public required int Polarities { get; init; }

    /// <summary>Weighted fraction of polarities the preferred mechanism does not fit.</summary>
    public required double MisfitFraction { get; init; }

    /// <summary>RMS angular difference (Kagan) of the acceptable mechanisms from the preferred one, degrees.</summary>
    public required double UncertaintyDeg { get; init; }

    /// <summary>Fraction of acceptable mechanisms within 30° of the preferred one.</summary>
    public required double Probability { get; init; }

    /// <summary>Station distribution ratio: mean |P radiation| at the stations, 0 (all nodal) to 1.</summary>
    public required double Stdr { get; init; }

    public required double AzimuthalGapDeg { get; init; }
    public required double TakeoffGapDeg { get; init; }
    public required string Quality { get; init; }

    /// <summary>The acceptable set (one plane of each), for plotting the uncertainty.</summary>
    public required List<NodalPlane> Acceptable { get; init; }
}

/// <summary>
/// Double-couple focal mechanism from P first-motion polarities by grid search over strike, dip and
/// rake. The misfit of a mechanism is the weighted fraction of polarities it predicts wrongly. The
/// acceptable set holds every mechanism within the misfit a given rate of wrong polarities allows of
/// the best one; the preferred mechanism is the member of that set closest to its average, and its
/// uncertainty the RMS rotation angle of the set from it, with quality classes A–D from misfit,
/// uncertainty, station distribution ratio and the azimuth and takeoff gaps (Hardebeck &amp; Shearer
/// 2002, Bull. Seismol. Soc. Am. 92(6), 2264–2276, doi:10.1785/0120010200).
/// </summary>
public static class FirstMotionInversion
{
    public static FirstMotionResult? Invert(IReadOnlyList<PolarityObservation> obs, FirstMotionSettings s)
    {
        var data = obs.Where(o => o.Polarity != 0 && o.Weight > 0).ToList();
        if (data.Count < Math.Max(4, s.MinPolarities)) return null;
        var dirs = data.Select(o => o.Direction).ToArray();
        var pol = data.Select(o => (double)Math.Sign(o.Polarity)).ToArray();
        var w = data.Select(o => o.Weight).ToArray();
        var wsum = w.Sum();

        var step = Math.Clamp(s.GridStepDeg, 1, 20);
        var strikes = Enumerable.Range(0, (int)Math.Round(360 / step)).Select(i => i * step).ToArray();
        var dips = Enumerable.Range(1, (int)Math.Round(90 / step)).Select(i => i * step).ToArray();
        var rakes = Enumerable.Range(0, (int)Math.Round(360 / step)).Select(i => -180 + i * step).ToArray();
        var misfits = new double[strikes.Length, dips.Length, rakes.Length];
        Parallel.For(0, strikes.Length, a =>
        {
            for (var b = 0; b < dips.Length; b++)
            for (var c = 0; c < rakes.Length; c++)
            {
                var m = MechanismGeometry.Tensor(new NodalPlane(strikes[a], dips[b], rakes[c]));
                double bad = 0;
                for (var k = 0; k < dirs.Length; k++)
                {
                    var r = m.Radiation(dirs[k].N, dirs[k].E, dirs[k].D);
                    if (Math.Sign(r) != pol[k]) bad += w[k];
                }
                misfits[a, b, c] = bad / wsum;
            }
        });

        var min = double.MaxValue;
        foreach (var x in misfits) min = Math.Min(min, x);
        // Allowed extra misfit: the binomial spread of the number of wrong polarities, at least one.
        var n = data.Count;
        var f = Math.Clamp(s.ExpectedErrorFraction, 0.01, 0.5);
        var extra = Math.Max(1, Math.Sqrt(n * f * (1 - f))) / n;
        var acceptable = new List<NodalPlane>();
        for (var a = 0; a < strikes.Length; a++)
        for (var b = 0; b < dips.Length; b++)
        for (var c = 0; c < rakes.Length; c++)
            if (misfits[a, b, c] <= min + extra + 1e-12) acceptable.Add(new NodalPlane(strikes[a], dips[b], rakes[c]));

        // Average of the set: tensors aligned to the first, summed; its best double couple.
        var reference = acceptable[0];
        var sum = new MomentTensor(0, 0, 0, 0, 0, 0);
        foreach (var p in acceptable)
        {
            var t = MechanismGeometry.Tensor(p);
            sum = new MomentTensor(sum.Nn + t.Nn, sum.Ee + t.Ee, sum.Dd + t.Dd, sum.Ne + t.Ne, sum.Nd + t.Nd, sum.Ed + t.Ed);
        }
        var (avg, _) = MechanismGeometry.BestDoubleCouple(sum);
        var preferred = acceptable.MinBy(p => MechanismGeometry.KaganAngle(p, avg));
        _ = reference;
        var angles = acceptable.Select(p => MechanismGeometry.KaganAngle(p, preferred)).ToList();
        var rms = Math.Sqrt(angles.Average(x => x * x));
        var prob = angles.Count(x => x <= 30) / (double)angles.Count;

        var tensor = MechanismGeometry.Tensor(preferred);
        double bad2 = 0, amp = 0;
        for (var k = 0; k < dirs.Length; k++)
        {
            var r = tensor.Radiation(dirs[k].N, dirs[k].E, dirs[k].D);
            if (Math.Sign(r) != pol[k]) bad2 += w[k];
            amp += w[k] * Math.Abs(r);
        }
        var misfit = bad2 / wsum;
        var stdr = amp / wsum;
        // Gaps on the lower focal hemisphere: an upgoing ray is folded to its antipode (takeoff
        // 180 − i, azimuth + 180), where the P radiation is the same.
        var folded = data.Select(o => o.TakeoffDeg > 90 ? (Az: o.AzimuthDeg + 180, To: 180 - o.TakeoffDeg) : (Az: o.AzimuthDeg, To: o.TakeoffDeg)).ToList();
        var azGap = Gap(folded.Select(x => x.Az));
        var toGap = TakeoffGap(folded.Select(x => x.To));
        var quality =
            prob > 0.8 && rms <= 25 && misfit <= 0.15 && stdr >= 0.5 && azGap <= 90 && toGap <= 60 ? "A" :
            prob > 0.6 && rms <= 35 && misfit <= 0.20 && stdr >= 0.4 && azGap <= 90 && toGap <= 60 ? "B" :
            prob > 0.5 && rms <= 45 && misfit <= 0.30 && stdr >= 0.3 && azGap <= 90 && toGap <= 60 ? "C" : "D";
        return new FirstMotionResult
        {
            Plane1 = preferred, Plane2 = MechanismGeometry.Auxiliary(preferred), Polarities = n, MisfitFraction = misfit,
            UncertaintyDeg = rms, Probability = prob, Stdr = stdr, AzimuthalGapDeg = azGap, TakeoffGapDeg = toGap, Quality = quality,
            Acceptable = acceptable.Count > 400 ? acceptable.Where((_, i) => i % (acceptable.Count / 400 + 1) == 0).ToList() : acceptable
        };
    }

    private static double Gap(IEnumerable<double> azimuths)
    {
        var a = azimuths.Select(MechanismGeometry.Wrap360).OrderBy(x => x).ToArray();
        if (a.Length < 2) return 360;
        var gap = 360 - a[^1] + a[0];
        for (var i = 1; i < a.Length; i++) gap = Math.Max(gap, a[i] - a[i - 1]);
        return gap;
    }

    /// <summary>Largest gap in (folded) takeoff angle over 0–90°, the ends included.</summary>
    private static double TakeoffGap(IEnumerable<double> takeoffs)
    {
        var t = takeoffs.OrderBy(x => x).ToArray();
        if (t.Length == 0) return 90;
        var gap = Math.Max(t[0], 90 - t[^1]);
        for (var i = 1; i < t.Length; i++) gap = Math.Max(gap, t[i] - t[i - 1]);
        return gap;
    }
}

/// <summary>
/// Automatic P first-motion polarity at a pick: the sign of the first swing that rises clearly
/// above the pre-event noise, read on the vertical with the mean of the noise window removed and
/// the channel's orientation applied (SEED: a vertical with dip −90° is positive up).
/// </summary>
public static class PolarityReader
{
    public sealed record Reading(int Polarity, double Ratio);

    /// <param name="dipDeg">Channel dip from the metadata: −90 positive up, +90 positive down.</param>
    /// <param name="minRatio">The first swing must exceed this many noise standard deviations.</param>
    public static Reading? Read(Trace t, DateTime pick, double dipDeg = -90, double pickUncertainty = 0.05, double minRatio = 4)
    {
        var fs = t.SampleRate;
        var ip = t.IndexOf(pick);
        var n0 = ip - (int)(2 * fs);
        var n1 = ip - (int)(Math.Max(0.1, 2 * pickUncertainty) * fs);
        if (n0 < 0 || n1 <= n0 + 10 || ip + (int)(0.5 * fs) >= t.Data.Length) return null;
        double mean = 0;
        for (var i = n0; i < n1; i++) mean += t.Data[i];
        mean /= n1 - n0;
        double var = 0;
        for (var i = n0; i < n1; i++) var += (t.Data[i] - mean) * (t.Data[i] - mean);
        var sd = Math.Sqrt(var / (n1 - n0));
        if (sd <= 0) sd = 1e-12;
        // Search from just before the pick (its uncertainty) to 0.3 s after it.
        var a = Math.Max(n1, ip - (int)(pickUncertainty * fs));
        var b = Math.Min(t.Data.Length - 1, ip + (int)(0.3 * fs));
        for (var i = a; i <= b; i++)
        {
            var x = t.Data[i] - mean;
            if (Math.Abs(x) < minRatio * sd) continue;
            // The first sample above the threshold: its sign is the first motion. Its size at the
            // following extremum gives the confidence.
            var sign = Math.Sign(x);
            var peak = Math.Abs(x);
            for (var j = i + 1; j <= Math.Min(b, i + (int)(0.2 * fs)); j++)
            {
                var y = (t.Data[j] - mean) * sign;
                if (y <= 0) break;
                peak = Math.Max(peak, y);
            }
            var orientation = dipDeg > 0 ? -1 : 1;
            return new Reading(sign * orientation, peak / sd);
        }
        return null;
    }
}
