using TomoStar.Core.Geo;
using TomoStar.Core.Model;

namespace TomoStar.Core.Tomography;

/// <summary>
/// One set of iterations of a double-difference inversion, with its weights (the weighting scheme of
/// Zhang &amp; Thurber 2003, 2006, and of Waldhauser &amp; Ellsworth 2000). The usual schedule starts with the absolute
/// times weighted up, so they set the large-scale model and the absolute locations, and ends with
/// the differential times weighted up, so they sharpen the source region and the relative locations.
/// </summary>
public sealed class DoubleDifferenceIterationSet
{
    /// <summary>Iterations with these weights.</summary>
    public int Iterations { get; set; } = 2;

    /// <summary>A priori weight of the absolute times relative to the differential ones (applied to the absolute rows).</summary>
    public double AbsoluteWeight { get; set; } = 1;

    /// <summary>A priori weight of the catalogue differential P and S times; 0 leaves them out.</summary>
    public double DifferentialWeightP { get; set; } = 1;

    public double DifferentialWeightS { get; set; } = 0.5;

    /// <summary>
    /// Residual cutoff of the differential times: at or above 1 a factor of the scaled median
    /// absolute deviation (MAD / 0.6745) of their residuals, below 1 a fixed cutoff in seconds, 0 or
    /// less none. Inside the cutoff the weight falls off as (1 − (|r| / c)³)³.
    /// </summary>
    public double ResidualCutoff { get; set; }

    /// <summary>Largest event separation, km, for a differential time; its weight falls off as (1 − (s / d)³)³. 0 or less: no limit.</summary>
    public double MaxSeparationKm { get; set; }

    public DoubleDifferenceIterationSet Clone() => (DoubleDifferenceIterationSet)MemberwiseClone();
}

/// <summary>
/// Double-difference tomography (Zhang &amp; Thurber 2003, BSSA 93(5), 1875-1889): differential
/// travel times of event pairs recorded at a common station are inverted together with the absolute
/// times for velocity and hypocentres. Most of the path is shared by the two rays, so the differential
/// time is sensitive to the structure near the sources and to the relative position of the events,
/// and insensitive to the structure near the station and to the station term.
/// </summary>
public sealed class DoubleDifferenceSettings
{
    public bool Enabled { get; set; }

    // Pair selection.

    /// <summary>Largest distance between an event pair and the station, km.</summary>
    public double MaxStationDistanceKm { get; set; } = 200;

    /// <summary>Largest hypocentral separation of a pair, km.</summary>
    public double MaxPairSeparationKm { get; set; } = 10;

    /// <summary>Neighbours sought for each event.</summary>
    public int MaxNeighbours { get; set; } = 10;

    /// <summary>Phase pairs needed for a neighbour to count as strongly linked.</summary>
    public int MinLinks { get; set; } = 8;

    /// <summary>Fewest phase pairs kept for a pair; pairs with fewer are dropped.</summary>
    public int MinObservations { get; set; } = 8;

    /// <summary>Most phase pairs kept for a pair, taken from the stations closest to the pair.</summary>
    public int MaxObservations { get; set; } = 50;

    /// <summary>Rows whose relative weight (a priori × separation × residual) falls below this are left out.</summary>
    public double MinWeight { get; set; } = 1e-3;

    /// <summary>
    /// Iteration sets, run in order; they replace <see cref="TomographySettings.Iterations"/>.
    /// The default follows the usual progression: absolute data first, differential data
    /// weighted up and reweighted by residual and separation afterwards.
    /// </summary>
    public List<DoubleDifferenceIterationSet> Sets { get; set; } =
    [
        new() { Iterations = 2, AbsoluteWeight = 1, DifferentialWeightP = 0.1, DifferentialWeightS = 0.05 },
        new() { Iterations = 2, AbsoluteWeight = 0.3, DifferentialWeightP = 1, DifferentialWeightS = 0.5, ResidualCutoff = 6, MaxSeparationKm = 10 },
        new() { Iterations = 2, AbsoluteWeight = 0.1, DifferentialWeightP = 1, DifferentialWeightS = 0.5, ResidualCutoff = 5, MaxSeparationKm = 6 }
    ];

    public int TotalIterations => Sets.Sum(x => Math.Max(0, x.Iterations));

    /// <summary>The iteration set of nonlinear iteration <paramref name="iteration"/> (0-based); the last set beyond the end.</summary>
    public DoubleDifferenceIterationSet SetOf(int iteration)
    {
        var k = 0;
        foreach (var set in Sets)
        {
            k += Math.Max(0, set.Iterations);
            if (iteration < k) return set;
        }
        return Sets.Count > 0 ? Sets[^1] : new DoubleDifferenceIterationSet();
    }

    public DoubleDifferenceSettings Clone()
    {
        var c = (DoubleDifferenceSettings)MemberwiseClone();
        c.Sets = Sets.Select(x => x.Clone()).ToList();
        return c;
    }

    public string Describe() => Enabled
        ? $"double difference: pairs within {MaxPairSeparationKm:0.#} km, {MaxNeighbours} neighbours, {MinObservations}-{MaxObservations} phase pairs; {Sets.Count} iteration sets ({TotalIterations} iterations)"
        : "double difference off";
}

/// <summary>A differential time: the same phase at the same station for two events (<see cref="ArrivalA"/> − <see cref="ArrivalB"/>).</summary>
public readonly record struct DifferentialTime(int EventA, int EventB, int ArrivalA, int ArrivalB);

/// <summary>Catalogue differential times of event pairs (nearest-neighbour selection after Waldhauser &amp; Ellsworth 2000).</summary>
public static class DoubleDifferencePairs
{
    /// <summary>
    /// For each event, its nearest events within <see cref="DoubleDifferenceSettings.MaxPairSeparationKm"/>
    /// are linked until <see cref="DoubleDifferenceSettings.MaxNeighbours"/> of them share at least
    /// <see cref="DoubleDifferenceSettings.MinLinks"/> phase readings (weaker neighbours met on the
    /// way are linked too, but do not count). For each pair, the common readings within
    /// <see cref="DoubleDifferenceSettings.MaxStationDistanceKm"/> are taken from the station closest
    /// to the pair outwards, up to <see cref="DoubleDifferenceSettings.MaxObservations"/>; a pair with
    /// fewer than <see cref="DoubleDifferenceSettings.MinObservations"/> is dropped. A delay longer
    /// than the time a wave needs to travel between the two sources (4 km/s for P, 2.3 km/s for S,
    /// plus 0.5 s for the location errors) is an outlier and is left out.
    /// </summary>
    public static List<DifferentialTime> Build(ObservationSet data, DoubleDifferenceSettings s, Action<string>? log = null)
    {
        var n = data.Events.Count;
        var pos = data.Events.Select(e => GeoMath.ToCartesian(e.Lon, e.Lat, e.DepthKm)).ToArray();
        var staPos = data.Stations.Select(x => GeoMath.ToCartesian(x.Lon, x.Lat, x.DepthKm)).ToArray();
        // Readings of each event by (station, phase).
        var readings = new Dictionary<(int Station, Phase Phase), int>[n];
        for (var i = 0; i < n; i++) readings[i] = [];
        for (var k = 0; k < data.Arrivals.Count; k++)
        {
            var a = data.Arrivals[k];
            if (a.Sigma <= 0 || double.IsNaN(a.Time)) continue;
            readings[a.Event].TryAdd((a.Station, a.Phase), k);
        }
        // Neighbour search on a uniform grid of cells as wide as the search radius.
        var cellKm = Math.Max(0.1, s.MaxPairSeparationKm);
        var cells = new Dictionary<(int, int, int), List<int>>();
        (int, int, int) Cell(Vec3 p) => ((int)Math.Floor(p.X / cellKm), (int)Math.Floor(p.Y / cellKm), (int)Math.Floor(p.Z / cellKm));
        for (var i = 0; i < n; i++)
        {
            if (data.Events[i].Outside) continue;
            var c = Cell(pos[i]);
            if (!cells.TryGetValue(c, out var list)) cells[c] = list = [];
            list.Add(i);
        }

        // The phase pairs of one event pair (a, b with a < b): the same whichever event found the other.
        (List<DifferentialTime> Times, int Outliers, bool Strong) PairTimes(int a, int b, double d, List<(int Station, Phase Phase)> common)
        {
            var mid = 0.5 * (pos[a] + pos[b]);
            var ea = data.Events[a];
            var chosen = new List<DifferentialTime>();
            var outliers = 0;
            foreach (var (st, ph) in common.OrderBy(c => (staPos[c.Station] - mid).Length).ThenBy(c => c.Station).ThenBy(c => c.Phase))
            {
                if (chosen.Count >= s.MaxObservations) break;
                var sta = data.Stations[st];
                if (GeoMath.SurfaceDistanceKm(ea.Lon, ea.Lat, sta.Lon, sta.Lat) > s.MaxStationDistanceKm) continue;
                var ka = readings[a][(st, ph)];
                var kb = readings[b][(st, ph)];
                var delay = (data.Arrivals[ka].Time - ea.T0) - (data.Arrivals[kb].Time - data.Events[b].T0);
                if (Math.Abs(delay) > d / (ph == Phase.P ? 4.0 : 2.3) + 0.5) { outliers++; continue; }
                chosen.Add(new DifferentialTime(a, b, ka, kb));
            }
            return (chosen, outliers, common.Count >= s.MinLinks);
        }

        // Each event's neighbours in parallel; a pair found from both ends is computed identically, so
        // the merge below keeps one copy and the result does not depend on the thread schedule.
        var found = new List<(int A, int B, double D, List<DifferentialTime> Times, int Outliers)>[n];
        var weak = new bool[n];
        Parallel.For(0, n, i =>
        {
            found[i] = [];
            if (data.Events[i].Outside || readings[i].Count == 0) return;
            var (cx, cy, cz) = Cell(pos[i]);
            var candidates = new List<(int J, double D)>();
            for (var dx = -1; dx <= 1; dx++)
            for (var dy = -1; dy <= 1; dy++)
            for (var dz = -1; dz <= 1; dz++)
            {
                if (!cells.TryGetValue((cx + dx, cy + dy, cz + dz), out var list)) continue;
                foreach (var j in list)
                {
                    if (j == i) continue;
                    var d = (pos[j] - pos[i]).Length;
                    if (d <= s.MaxPairSeparationKm) candidates.Add((j, d));
                }
            }
            var strong = 0;
            foreach (var (j, d) in candidates.OrderBy(c => c.D).ThenBy(c => c.J))
            {
                if (strong >= s.MaxNeighbours) break;
                var common = readings[i].Keys.Where(readings[j].ContainsKey).ToList();
                if (common.Count >= s.MinLinks) strong++;
                var (a, b) = i < j ? (i, j) : (j, i);
                var (times, outliers, _) = PairTimes(a, b, d, common);
                found[i].Add((a, b, d, times, outliers));
            }
            weak[i] = strong == 0;
        });

        var pairs = new HashSet<(int, int)>();
        var result = new List<DifferentialTime>();
        int strongPairs = 0, outlierCount = 0;
        double sepSum = 0;
        for (var i = 0; i < n; i++)
            foreach (var (a, b, d, times, outliers) in found[i])
            {
                if (!pairs.Add((a, b))) continue;
                outlierCount += outliers;
                if (times.Count < s.MinObservations) continue;
                result.AddRange(times);
                if (times.Count >= s.MinLinks) { strongPairs++; sepSum += d; }
            }
        log?.Invoke(string.Create(System.Globalization.CultureInfo.InvariantCulture,
            $"Double difference: {result.Count} differential times for {pairs.Count} event pairs ({strongPairs} strongly linked, mean separation {(strongPairs > 0 ? sepSum / strongPairs : 0):0.0} km); {weak.Count(w => w)} events without a strong neighbour; {outlierCount} delays longer than the travel time between the sources left out."));
        return result;
    }

    /// <summary>Median and scaled median absolute deviation (MAD / 0.6745, σ for Gaussian noise).</summary>
    public static (double Median, double Sigma) RobustSpread(IReadOnlyList<double> values)
    {
        if (values.Count == 0) return (0, 0);
        var sorted = values.OrderBy(v => v).ToArray();
        var med = sorted.Length % 2 == 1 ? sorted[sorted.Length / 2] : 0.5 * (sorted[sorted.Length / 2 - 1] + sorted[sorted.Length / 2]);
        var dev = values.Select(v => Math.Abs(v - med)).OrderBy(v => v).ToArray();
        var mad = dev.Length % 2 == 1 ? dev[dev.Length / 2] : 0.5 * (dev[dev.Length / 2 - 1] + dev[dev.Length / 2]);
        return (med, mad / 0.67449);
    }

    /// <summary>Tricube taper (1 − (x / limit)³)³ for 0 ≤ x ≤ limit, 0 beyond: the distance and residual reweighting.</summary>
    public static double Tricube(double x, double limit)
    {
        if (limit <= 0) return 1;
        var u = Math.Abs(x) / limit;
        if (u >= 1) return 0;
        var t = 1 - u * u * u;
        return t * t * t;
    }
}
