using TomoStar.Core.Model;

namespace TomoStar.Core.Attenuation;

/// <summary>
/// How well an event's corner frequency separates from t*. The corner frequency is found by a grid
/// search on the misfit; its interval is the range of grid values whose misfit stays within the
/// Δχ² = 3.84 bound of one parameter (95 %) with the variance estimated from the best fit, and the
/// change of each t* across that interval is added to its uncertainty (for a large event whose corner
/// lies below the band the corner is poorly resolved yet t* hardly changes). Whether the corner can be
/// trusted at all follows the bandwidth rule of Abercrombie (2015): a corner above a third of the
/// highest usable frequency is underestimated, and t* absorbs the error (the M 2.5 event of the coda
/// synthetic: true corner 12 Hz in a band ending near 28 Hz, found at 40 Hz, t* 20 ms too high).
/// </summary>
public static class CornerTradeOff
{
    /// <summary>Indices of the contiguous grid range around <paramref name="best"/> within the misfit bound.</summary>
    public static (int Lo, int Hi) Interval(IReadOnlyList<double> rss, int best, int nData, int nParams)
    {
        var nu = Math.Max(1, nData - nParams);
        var limit = rss[best] * (1 + 3.84 / nu);
        int lo = best, hi = best;
        while (lo > 0 && rss[lo - 1] <= limit) lo--;
        while (hi < rss.Count - 1 && rss[hi + 1] <= limit) hi++;
        return (lo, hi);
    }

    /// <summary>
    /// Fills the interval, the trade-off and the resolved flag of an event's measurements (one phase);
    /// <paramref name="spread"/> holds each record's largest t* change across the interval, which is
    /// added to its uncertainty. The corner is resolved when it lies <see cref="TStarSettings.CornerBandFactor"/>
    /// times below the highest fitted frequency (median of the records); unresolved events are
    /// disabled when the settings say so.
    /// </summary>
    public static void Apply(IReadOnlyList<TStarMeasurement> records, IReadOnlyList<double> spread, double lowHz, double highHz, double fc, TStarSettings s,
        bool forceUnresolved = false)
    {
        if (records.Count == 0) return;
        var fmax = records.Select(m => m.FitMaxHz).OrderBy(x => x).ElementAt(records.Count / 2);
        var resolved = !forceUnresolved && fc <= fmax / s.CornerBandFactor;
        for (var i = 0; i < records.Count; i++)
        {
            var m = records[i];
            m.CornerLowHz = lowHz;
            m.CornerHighHz = highHz;
            m.CornerTradeOff = spread[i];
            m.CornerResolved = resolved;
            m.Uncertainty = Math.Sqrt(m.Uncertainty * m.Uncertainty + spread[i] * spread[i]);
            if (!resolved && s.ExcludeUnresolvedCorners) m.Disabled = true;
        }
    }
}
