// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using System.Globalization;
using TomoStar.Core.Model;

namespace TomoStar.Core.Attenuation;

/// <summary>
/// A minimum magnitude for the t* measurement, from the data. Small earthquakes have their corner
/// frequency inside the fit band, where it trades off with t*; how small is too small depends on the
/// band the records allow (noise, instruments, distance, attenuation) and on the stress drop.
///
/// Per event (preferred, with enough events): the lowest magnitude from which, in every window of
/// consecutive magnitudes, most events have their corner resolved
/// (<see cref="TStarMeasurement.CornerResolved"/>); the unresolved events above it are flagged one by
/// one. Otherwise from Brune scaling: the stress drop of
/// the resolved events (f_c = 0.37 β (16 Δσ / 7 M₀)^⅓, Brune 1970, β = 3.5 km/s,
/// M₀ = 10^(1.5 M + 9.1) N·m) and the magnitude whose corner reaches a third of the usual upper end
/// of the fit band.
/// </summary>
public static class TStarQuality
{
    public const double ResolvedFraction = 0.5;
    public const int MinEventsForPerEvent = 8;

    public sealed record Recommendation(double MinMagnitude, string Method, string Explanation, int Events, int Unresolved, int EventsKept);

    private sealed record EventRow(string Id, double Magnitude, double Fc, bool Resolved, double FitMaxHz);

    public static Recommendation? Recommend(IEnumerable<TStarMeasurement> measurements, IReadOnlyDictionary<string, EventRecord> events)
    {
        var rows = measurements.Where(m => m.Phase == Phase.P && double.IsFinite(m.CornerTradeOff))
            .GroupBy(m => m.EventId)
            .Select(g => new EventRow(g.Key, events.TryGetValue(g.Key, out var e) ? e.Magnitude : double.NaN, g.First().CornerFrequencyHz,
                g.First().CornerResolved, Median(g.Select(m => m.FitMaxHz))))
            .ToList();
        if (rows.Count == 0) return null;
        var withM = rows.Where(r => double.IsFinite(r.Magnitude)).OrderBy(r => r.Magnitude).ToList();
        var unresolved = rows.Count(r => !r.Resolved);
        var inv = CultureInfo.InvariantCulture;

        if (withM.Count >= MinEventsForPerEvent)
        {
            // Resolved fraction in windows of consecutive magnitudes: the recommendation is the lowest
            // magnitude from which every window has most of its events resolved. Unresolved events
            // above it are flagged one by one.
            var k = Math.Max(5, withM.Count / 6);
            var fraction = new double[withM.Count - k + 1];
            for (var i = 0; i < fraction.Length; i++) fraction[i] = withM.Skip(i).Take(k).Count(r => r.Resolved) / (double)k;
            var first = fraction.Length;
            while (first > 0 && fraction[first - 1] >= ResolvedFraction) first--;
            if (first == 0)
                return new Recommendation(0, "per event", string.Create(inv,
                    $"Most events are resolved at every magnitude (M {withM[0].Magnitude:0.0}-{withM[^1].Magnitude:0.0}, {withM.Count(r => r.Resolved)} of {withM.Count}): " +
                    $"no minimum magnitude is needed; the {unresolved} unresolved events are flagged one by one."),
                    rows.Count, unresolved, rows.Count);
            if (first < fraction.Length)
            {
                var mMin = Math.Floor(withM[first].Magnitude * 10 + 1e-6) / 10;
                var below = withM.Where(r => r.Magnitude < mMin).ToList();
                var above = withM.Where(r => r.Magnitude >= mMin).ToList();
                return new Recommendation(mMin, "per event", string.Create(inv,
                    $"Below M {mMin:0.0} most events have an unresolved corner frequency ({below.Count(r => !r.Resolved)} of {below.Count}); from M {mMin:0.0} up, " +
                    $"{above.Count(r => r.Resolved)} of {above.Count} are resolved. Unresolved: corner and t* trade off (corner inside the fit band, which ends near {Median(rows.Select(x => x.FitMaxHz)):0} Hz)."),
                    rows.Count, unresolved, rows.Count - below.Count);
            }
        }

        // Brune scaling from the resolved events.
        var resolved = withM.Where(r => r.Resolved).ToList();
        if (resolved.Count == 0) return null;
        const double beta = 3500;
        var logStress = resolved.Select(r =>
        {
            var m0 = Math.Pow(10, 1.5 * r.Magnitude + 9.1);
            return Math.Log10(7.0 / 16 * m0 * Math.Pow(r.Fc / (0.37 * beta), 3));
        }).ToList();
        var stress = Math.Pow(10, Median(logStress)); // Pa
        var band = Median(rows.Select(r => r.FitMaxHz));
        var fcLimit = band / 3;
        var m0Limit = 16 * stress / 7 * Math.Pow(0.37 * beta / fcLimit, 3);
        var mBrune = Math.Ceiling(((Math.Log10(m0Limit) - 9.1) / 1.5) * 10 - 1e-6) / 10;
        return new Recommendation(mBrune, "Brune scaling", string.Create(inv,
            $"Stress drop {stress / 1e6:0.0} MPa (median of {resolved.Count} events with a resolved corner), fit band up to {band:0} Hz: " +
            $"the corner reaches {fcLimit:0.0} Hz, a third of the band, at M {mBrune:0.0}."),
            rows.Count, unresolved, rows.Count(r => !double.IsFinite(r.Magnitude) || r.Magnitude >= mBrune));
    }

    private static double Median(IEnumerable<double> values)
    {
        var v = values.Where(double.IsFinite).OrderBy(x => x).ToArray();
        return v.Length == 0 ? double.NaN : v[v.Length / 2];
    }
}
