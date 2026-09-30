// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

namespace TomoStar.Core.Model;

/// <summary>A seismic station: position and the static corrections of an earlier run, if any.</summary>
public sealed class StationRecord
{
    /// <summary>Identifier, usually NET.STA; matched without regard to case.</summary>
    public required string Id { get; init; }

    public double Lon { get; set; }
    public double Lat { get; set; }

    /// <summary>Elevation above sea level, m.</summary>
    public double ElevationM { get; set; }

    /// <summary>Depth below sea level, km (negative above sea level), the vertical coordinate of the grid.</summary>
    public double DepthKm => -ElevationM / 1000.0;

    /// <summary>Static P correction, s (only used when station terms are inverted).</summary>
    public double CorrectionP { get; set; }

    /// <summary>Static S correction, s.</summary>
    public double CorrectionS { get; set; }
}

/// <summary>Who made a pick, which decides the pick used when a station has several for one phase.</summary>
public enum PickOrigin
{
    /// <summary>From a bulletin or catalogue (analyst reviewed).</summary>
    Catalog,

    /// <summary>From an automatic picker.</summary>
    Automatic,

    /// <summary>Picked by hand.</summary>
    Manual
}

/// <summary>When automatic picks enter an inversion.</summary>
public enum AutomaticPicks
{
    /// <summary>Only for events the analysts picked little (default).</summary>
    IfFewAnalystPicks,

    /// <summary>Wherever there is no analyst pick for the station and phase.</summary>
    Always,

    /// <summary>Never.</summary>
    Never
}

/// <summary>An arrival-time reading.</summary>
public sealed class PickRecord
{
    public required string StationId { get; init; }
    public Phase Phase { get; init; }

    /// <summary>Absolute arrival time, UTC.</summary>
    public DateTime Time { get; init; }

    /// <summary>One-sigma pick uncertainty, s.</summary>
    public double Sigma { get; init; } = 0.1;

    /// <summary>HYPO71-style quality class 0 (best) to 4 (rejected).</summary>
    public int Quality { get; init; }

    public PickOrigin Origin { get; init; } = PickOrigin.Catalog;

    /// <summary>Excluded from inversions.</summary>
    public bool Disabled { get; init; }

    /// <summary>Usable at all: not disabled and of quality below 4.</summary>
    public bool Usable => !Disabled && Quality < 4 && Sigma > 0 && double.IsFinite(Sigma);

    /// <summary>Weight of the pick in a location: 1/sigma (sigma at least 10 ms), zero when not usable.</summary>
    public double Weight => Usable ? 1.0 / Math.Max(0.01, Sigma) : 0;

    /// <summary>Last computed residual (observed minus calculated), s; NaN when never computed.</summary>
    public double Residual { get; set; } = double.NaN;

    /// <summary>Channel the pick was made on (automatic picks), for the output tables.</summary>
    public string Channel { get; init; } = "";

    /// <summary>Signal-to-noise ratio of an automatic pick (NaN for other picks).</summary>
    public double Snr { get; init; } = double.NaN;
}

/// <summary>A hypocentre with its uncertainties, as produced by the event locator.</summary>
public sealed class Hypocentre
{
    public DateTime OriginTime { get; set; }
    public double Lat { get; set; }
    public double Lon { get; set; }
    public double DepthKm { get; set; }

    /// <summary>RMS residual of the solution, s; NaN when unknown.</summary>
    public double Rms { get; set; } = double.NaN;

    /// <summary>1-sigma semi-major axis of the horizontal error ellipse, km.</summary>
    public double ErrorHorizontalKm { get; set; } = double.NaN;

    /// <summary>1-sigma semi-minor axis of the horizontal error ellipse, km.</summary>
    public double ErrorMinorKm { get; set; } = double.NaN;

    /// <summary>Azimuth of the major axis, degrees clockwise from north in [0, 180).</summary>
    public double ErrorAzimuthDeg { get; set; } = double.NaN;

    public double ErrorDepthKm { get; set; } = double.NaN;
    public double ErrorTimeS { get; set; } = double.NaN;
    public double AzimuthalGapDeg { get; set; } = double.NaN;
    public int PhaseCount { get; set; }

    public Hypocentre Clone() => (Hypocentre)MemberwiseClone();
}

/// <summary>A t* measurement with the details of its spectral fit, as produced by the t* estimators.</summary>
public sealed class TStarMeasurement
{
    public string EventId { get; set; } = "";
    public string StationId { get; set; } = "";
    public Phase Phase { get; set; } = Phase.P;

    /// <summary>Whole-path attenuation operator t* = integral of ds / (v Q), s.</summary>
    public double TStar { get; set; }

    public double Uncertainty { get; set; }
    public double CornerFrequencyHz { get; set; }
    public double LogOmega0 { get; set; }

    /// <summary>Frequency dependence used: t* is the value at 1 Hz, t*(f) = t* f^(-alpha) (0 = constant Q).</summary>
    public double Alpha { get; set; }

    public double FitMinHz { get; set; }
    public double FitMaxHz { get; set; }
    public double Snr { get; set; }

    /// <summary>Confidence interval of the event's corner frequency from its misfit curve (NaN: not computed).</summary>
    public double CornerLowHz { get; set; } = double.NaN;

    public double CornerHighHz { get; set; } = double.NaN;

    /// <summary>
    /// Largest change of this t* when the corner frequency moves within its interval, s (the
    /// corner and t* trade-off; NaN: not computed). Included in <see cref="Uncertainty"/>.
    /// </summary>
    public double CornerTradeOff { get; set; } = double.NaN;

    /// <summary>
    /// False when the event's corner frequency does not separate from t* (median trade-off above the
    /// limit, or corner at the top of its search range): typical of events too small for the band.
    /// </summary>
    public bool CornerResolved { get; set; } = true;

    public bool Disabled { get; set; }
    public double Residual { get; set; } = double.NaN;

    /// <summary>The measurement as an input record of the Q tomography.</summary>
    public TStarRecord ToRecord() => new()
    {
        EventId = EventId, StationId = StationId, Phase = Phase, TStar = TStar, Sigma = Math.Max(0.002, Uncertainty), Disabled = Disabled
    };
}

/// <summary>An earthquake (or a shot): its starting hypocentre and its picks.</summary>
public sealed class EventRecord
{
    public required string Id { get; init; }

    /// <summary>Origin time, UTC; the reference of every travel time of the event.</summary>
    public DateTime OriginTime { get; set; }

    public double Lon { get; set; }
    public double Lat { get; set; }
    public double DepthKm { get; set; }
    public double Magnitude { get; set; } = double.NaN;

    /// <summary>Location held fixed during the inversion (a shot, a mine blast of known position).</summary>
    public bool Fixed { get; set; }

    public List<PickRecord> Picks { get; } = [];
}

/// <summary>A t* measurement for one event and station (whole-path attenuation operator, s).</summary>
public sealed class TStarRecord
{
    public required string EventId { get; init; }
    public required string StationId { get; init; }
    public Phase Phase { get; init; } = Phase.P;

    /// <summary>t* = integral of ds / (v Q) along the ray, s.</summary>
    public double TStar { get; init; }

    /// <summary>One-sigma uncertainty, s.</summary>
    public double Sigma { get; init; } = 0.005;

    public bool Disabled { get; init; }
}

/// <summary>
/// Everything an inversion reads: stations, events with their picks, and t* measurements. It is
/// filled by the readers of <see cref="IO.CatalogueReader"/> (plain CSV files or a QUIVER project)
/// and turned into the index-based <see cref="Tomography.ObservationSet"/> the solvers work on.
/// </summary>
public sealed class Catalogue
{
    public List<StationRecord> Stations { get; } = [];
    public List<EventRecord> Events { get; } = [];
    public List<TStarRecord> TStar { get; } = [];

    /// <summary>Problems met while reading (unknown stations, unparsable lines), reported to the user.</summary>
    public List<string> Warnings { get; } = [];

    /// <summary>Stations by identifier, case-insensitive; the first of duplicated identifiers wins.</summary>
    public Dictionary<string, StationRecord> StationIndex()
    {
        var index = new Dictionary<string, StationRecord>(StringComparer.OrdinalIgnoreCase);
        foreach (var s in Stations) index.TryAdd(s.Id, s);
        return index;
    }

    /// <summary>Events by identifier.</summary>
    public Dictionary<string, EventRecord> EventIndex()
    {
        var index = new Dictionary<string, EventRecord>(StringComparer.Ordinal);
        foreach (var e in Events) index.TryAdd(e.Id, e);
        return index;
    }

    /// <summary>Number of usable P and S picks.</summary>
    public (int P, int S) PickCounts()
    {
        int p = 0, s = 0;
        foreach (var e in Events)
        foreach (var k in e.Picks)
        {
            if (!k.Usable) continue;
            if (k.Phase == Phase.P) p++;
            else s++;
        }
        return (p, s);
    }
}

/// <summary>How several picks of one station and phase reduce to the one an inversion uses.</summary>
public static class PickSelection
{
    /// <summary>
    /// The usable picks, one per station and phase: a manual pick first, then a catalogue one, then an
    /// automatic one; within the same origin the smallest uncertainty, and among equals the earliest
    /// (the first arrival). With <see cref="AutomaticPicks.IfFewAnalystPicks"/> automatic picks are used
    /// only for events with fewer than <paramref name="enoughAnalyst"/> analyst picks: where analysts
    /// picked an event well, the stations they left out are usually those with doubtful records.
    /// </summary>
    public static IEnumerable<PickRecord> Best(IEnumerable<PickRecord> picks, AutomaticPicks automatic = AutomaticPicks.IfFewAnalystPicks, int enoughAnalyst = 8)
    {
        var usable = picks.Where(p => p.Usable).ToList();
        if (automatic == AutomaticPicks.Never
            || (automatic == AutomaticPicks.IfFewAnalystPicks && usable.Count(p => Rank(p) < 2) >= enoughAnalyst))
            usable.RemoveAll(p => Rank(p) == 2);
        return usable.GroupBy(p => (Station: p.StationId.ToUpperInvariant(), p.Phase))
            .Select(g => g.OrderBy(Rank).ThenBy(p => p.Sigma).ThenBy(p => p.Time).First());
    }

    /// <summary>Preference of an origin: 0 manual, 1 catalogue, 2 automatic.</summary>
    public static int Rank(PickRecord p) => p.Origin switch
    {
        PickOrigin.Manual => 0,
        PickOrigin.Catalog => 1,
        _ => 2
    };
}
