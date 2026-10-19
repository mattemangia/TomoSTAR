// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using System.Globalization;
using TomoStar.Core.Geo;
using TomoStar.Core.IO;
using TomoStar.Core.Model;
using TomoStar.Core.Signal;

namespace TomoStar.Core.Location;

/// <summary>Settings of the waveform cross-correlation of event pairs.</summary>
public sealed class CrossCorrelationSettings
{
    /// <summary>Largest hypocentral separation of a pair, km.</summary>
    public double MaxSeparationKm { get; set; } = 5;

    /// <summary>Nearest neighbours correlated with each event (within <see cref="MaxSeparationKm"/>).</summary>
    public int MaxNeighbours { get; set; } = 30;

    /// <summary>Largest distance between an event and the station, km.</summary>
    public double MaxStationDistanceKm { get; set; } = 150;

    /// <summary>Band-pass of the records before correlation (zero phase), Hz.</summary>
    public double LowHz { get; set; } = 1.5;

    public double HighHz { get; set; } = 15;

    /// <summary>P window: from this long before the pick to <see cref="PAfterSeconds"/> after it, on the vertical.</summary>
    public double PBeforeSeconds { get; set; } = 0.2;

    public double PAfterSeconds { get; set; } = 1.0;

    /// <summary>S window, on the horizontals (their correlations summed); cut before the S pick at the P pick if the two are closer.</summary>
    public double SBeforeSeconds { get; set; } = 0.3;

    public double SAfterSeconds { get; set; } = 1.5;

    /// <summary>Correlate S on the horizontals too.</summary>
    public bool IncludeS { get; set; } = true;

    /// <summary>Largest shift searched around the alignment on the picks, s; a maximum at the edge of the range is rejected.</summary>
    public double MaxLagSeconds { get; set; } = 0.3;

    /// <summary>Lowest normalised correlation coefficient kept.</summary>
    public double MinCoefficient { get; set; } = 0.7;

    /// <summary>When automatic picks may align the windows.</summary>
    public AutomaticPicks AutomaticPicks { get; set; } = AutomaticPicks.Always;

    public CrossCorrelationSettings Clone() => (CrossCorrelationSettings)MemberwiseClone();

    public string Describe() => string.Create(CultureInfo.InvariantCulture,
        $"cross-correlation: pairs within {MaxSeparationKm:0.#} km ({MaxNeighbours} neighbours), stations within {MaxStationDistanceKm:0} km, {LowHz:0.##}–{HighHz:0.##} Hz, " +
        $"P {-PBeforeSeconds:0.##}…{PAfterSeconds:0.##} s{(IncludeS ? string.Create(CultureInfo.InvariantCulture, $", S {-SBeforeSeconds:0.##}…{SAfterSeconds:0.##} s") : "")}, lag ±{MaxLagSeconds:0.##} s, CC ≥ {MinCoefficient:0.##}");
}

/// <summary>
/// A differential arrival time measured by cross-correlation: the same phase of events A and B at one
/// station. <see cref="ArrivalDifference"/> is t_A − t_B of the absolute arrival times, s (origin times
/// are not in it, so it stays valid when the events are relocated); <see cref="Lag"/> is how much it
/// differs from the difference of the two picks.
/// </summary>
public sealed record CcMeasurement(string EventA, string EventB, string StationId, Phase Phase, double ArrivalDifference, double Lag, double Coefficient);

/// <summary>
/// Differential arrival times of neighbouring events from the cross-correlation of their waveforms
/// (Poupinet et al. 1984; as used for double-difference location by Waldhauser &amp; Ellsworth 2000 and
/// Schaff et al. 2004). Records of the two events at a station are band-passed, cut around the picks
/// (P on the vertical, S on the horizontals, whose correlation functions are summed) and correlated in
/// the time domain; the lag of the largest normalised coefficient, refined by a parabola through the
/// three highest samples, corrects the difference of the picks. Only positive maxima inside the search
/// range count (a maximum at its edge is a cycle skip or no similarity at all).
/// </summary>
public static class WaveformCrossCorrelation
{
    /// <summary>A cut, filtered window of one event's phase at one station.</summary>
    internal sealed class Window
    {
        public required string Instrument { get; init; }
        public required double SampleRate { get; init; }

        /// <summary>Components (by their last channel letter) cut from pick − before − lag to pick + after + lag.</summary>
        public required Dictionary<char, float[]> Components { get; init; }
        /// <summary>Time of the sample the window is aligned on (the one nearest the pick).</summary>
        public required DateTime Pick { get; init; }

        /// <summary>Samples before the pick in the window (before + lag).</summary>
        public required int Lead { get; init; }

        /// <summary>Samples of the template (before + after).</summary>
        public required int Length { get; init; }
    }

    /// <summary>
    /// Event pairs to correlate: for each event its nearest events within
    /// <see cref="CrossCorrelationSettings.MaxSeparationKm"/>, up to <see cref="CrossCorrelationSettings.MaxNeighbours"/>,
    /// as (i, j) with i &lt; j into <paramref name="events"/>.
    /// </summary>
    public static List<(int A, int B)> Pairs(IReadOnlyList<EventRecord> events, CrossCorrelationSettings s)
    {
        var pos = events.Select(e => GeoMath.ToCartesian(e.Lon, e.Lat, e.DepthKm)).ToArray();
        var cellKm = Math.Max(0.1, s.MaxSeparationKm);
        var cells = new Dictionary<(int, int, int), List<int>>();
        (int, int, int) Cell(Vec3 p) => ((int)Math.Floor(p.X / cellKm), (int)Math.Floor(p.Y / cellKm), (int)Math.Floor(p.Z / cellKm));
        for (var i = 0; i < pos.Length; i++)
        {
            var c = Cell(pos[i]);
            if (!cells.TryGetValue(c, out var list)) cells[c] = list = [];
            list.Add(i);
        }
        var set = new HashSet<(int, int)>();
        for (var i = 0; i < pos.Length; i++)
        {
            var (cx, cy, cz) = Cell(pos[i]);
            var near = new List<(int J, double D)>();
            for (var dx = -1; dx <= 1; dx++)
            for (var dy = -1; dy <= 1; dy++)
            for (var dz = -1; dz <= 1; dz++)
            {
                if (!cells.TryGetValue((cx + dx, cy + dy, cz + dz), out var list)) continue;
                foreach (var j in list)
                {
                    if (j == i) continue;
                    var d = (pos[j] - pos[i]).Length;
                    if (d <= s.MaxSeparationKm) near.Add((j, d));
                }
            }
            foreach (var (j, _) in near.OrderBy(x => x.D).ThenBy(x => x.J).Take(Math.Max(0, s.MaxNeighbours)))
                set.Add(i < j ? (i, j) : (j, i));
        }
        return set.OrderBy(p => p.Item1).ThenBy(p => p.Item2).ToList();
    }

    /// <summary>
    /// Measures every pair at every station, one station at a time (so only one station's records are in
    /// memory): <paramref name="loadStation"/> gives the records of a station as (event id, trace).
    /// </summary>
    public static List<CcMeasurement> Measure(IReadOnlyList<EventRecord> events, IReadOnlyDictionary<string, StationRecord> stations,
        Func<string, IEnumerable<(string EventId, Trace Trace)>> loadStation, CrossCorrelationSettings s,
        IProgress<(double, string)>? progress = null, Action<string>? log = null, CancellationToken ct = default)
    {
        var pairs = Pairs(events, s);
        log?.Invoke($"Cross-correlation: {pairs.Count} event pairs within {s.MaxSeparationKm.ToString("0.#", CultureInfo.InvariantCulture)} km.");
        // Picks of every event by station and phase.
        var picks = events.Select(e => PickSelection.Best(e.Picks, s.AutomaticPicks)
            .Where(p => p.Phase == Phase.P || s.IncludeS)
            .GroupBy(p => (p.StationId.ToUpperInvariant(), p.Phase)).ToDictionary(g => g.Key, g => g.First())).ToArray();
        var pairsAt = new Dictionary<string, List<(int A, int B)>>(StringComparer.OrdinalIgnoreCase);
        foreach (var (a, b) in pairs)
            foreach (var key in picks[a].Keys.Where(picks[b].ContainsKey).Select(k => k.Item1).Distinct())
            {
                if (!pairsAt.TryGetValue(key, out var list)) pairsAt[key] = list = [];
                list.Add((a, b));
            }
        var indexOf = new Dictionary<string, int>(StringComparer.Ordinal);
        for (var i = 0; i < events.Count; i++) indexOf.TryAdd(events[i].Id, i);

        var result = new List<CcMeasurement>();
        var stationIds = pairsAt.Keys.OrderBy(k => k, StringComparer.Ordinal).ToList();
        var done = 0;
        long tried = 0;
        var gate = new object();
        Parallel.ForEach(stationIds, Compute.ComputeSettings.Options(ct), sid =>
        {
            var station = stations.TryGetValue(sid, out var st) ? st : null;
            var near = pairsAt[sid].Where(p => station == null || Distance(events[p.A], station) <= s.MaxStationDistanceKm).ToList();
            var local = new List<CcMeasurement>();
            long localTried = 0;
            if (near.Count > 0)
            {
                // Cut windows once per event and phase.
                var needed = near.SelectMany(p => new[] { p.A, p.B }).ToHashSet();
                var traces = new Dictionary<int, List<Trace>>();
                foreach (var (id, tr) in loadStation(sid))
                    if (indexOf.TryGetValue(id, out var i) && needed.Contains(i))
                    {
                        if (!traces.TryGetValue(i, out var l)) traces[i] = l = [];
                        l.Add(tr);
                    }
                var windows = new Dictionary<(int, Phase), Window?>();
                Window? W(int e, Phase ph)
                {
                    if (windows.TryGetValue((e, ph), out var w)) return w;
                    w = null;
                    if (traces.TryGetValue(e, out var list) && picks[e].TryGetValue((sid.ToUpperInvariant(), ph), out var pk))
                    {
                        DateTime? other = null;
                        if (ph == Phase.S && picks[e].TryGetValue((sid.ToUpperInvariant(), Phase.P), out var pp)) other = pp.Time;
                        w = Cut(list, pk.Time, ph, other, s);
                    }
                    windows[(e, ph)] = w;
                    return w;
                }
                foreach (var (a, b) in near)
                {
                    ct.ThrowIfCancellationRequested();
                    foreach (var ph in s.IncludeS ? new[] { Phase.P, Phase.S } : [Phase.P])
                    {
                        if (!picks[a].ContainsKey((sid.ToUpperInvariant(), ph)) || !picks[b].ContainsKey((sid.ToUpperInvariant(), ph))) continue;
                        var wa = W(a, ph);
                        var wb = W(b, ph);
                        if (wa == null || wb == null) continue;
                        localTried++;
                        if (Correlate(wa, wb, s) is not { } m) continue;
                        var pickDiff = (wa.Pick - wb.Pick).TotalSeconds;
                        local.Add(new CcMeasurement(events[a].Id, events[b].Id, sid, ph, pickDiff - m.Lag, -m.Lag, m.Coefficient));
                    }
                }
            }
            lock (gate)
            {
                result.AddRange(local);
                tried += localTried;
                done++;
                progress?.Report(((double)done / stationIds.Count, $"Cross-correlation: {done}/{stationIds.Count} stations"));
            }
        });
        result.Sort((x, y) =>
        {
            var c = indexOf[x.EventA].CompareTo(indexOf[y.EventA]);
            if (c == 0) c = indexOf[x.EventB].CompareTo(indexOf[y.EventB]);
            if (c == 0) c = string.CompareOrdinal(x.StationId, y.StationId);
            return c != 0 ? c : x.Phase.CompareTo(y.Phase);
        });
        log?.Invoke(string.Create(CultureInfo.InvariantCulture,
            $"Cross-correlation: {result.Count} differential times kept of {tried} correlated (P {result.Count(m => m.Phase == Phase.P)}, S {result.Count(m => m.Phase == Phase.S)}); " +
            $"median coefficient {Median(result.Select(m => m.Coefficient)):0.00}, median |lag| {Median(result.Select(m => Math.Abs(m.Lag))):0.000} s."));
        return result;
    }

    private static double Distance(EventRecord e, StationRecord st) =>
        GeoMath.SurfaceDistanceKm(e.Lon, e.Lat, st.Lon, st.Lat);

    private static double Median(IEnumerable<double> v)
    {
        var a = v.OrderBy(x => x).ToArray();
        return a.Length == 0 ? double.NaN : a.Length % 2 == 1 ? a[a.Length / 2] : 0.5 * (a[a.Length / 2 - 1] + a[a.Length / 2]);
    }

    /// <summary>
    /// The window of one phase: the components of one instrument (location and the first two letters of
    /// the channel, the highest sample rate first), P on the vertical, S on the horizontals, filtered on
    /// a longer stretch and cut to [pick − before − lag, pick + after + lag].
    /// </summary>
    internal static Window? Cut(List<Trace> traces, DateTime pick, Phase phase, DateTime? pBeforeS, CrossCorrelationSettings s)
    {
        var before = phase == Phase.P ? s.PBeforeSeconds : s.SBeforeSeconds;
        var after = phase == Phase.P ? s.PAfterSeconds : s.SAfterSeconds;
        // An S window must not reach back into the P onset.
        if (phase == Phase.S && pBeforeS is { } tp) before = Math.Min(before, Math.Max(0.05, 0.5 * (pick - tp).TotalSeconds));
        foreach (var group in traces.GroupBy(t => $"{t.Location}.{(t.Channel.Length >= 2 ? t.Channel[..2] : t.Channel)}")
                     .OrderByDescending(g => g.Max(t => t.SampleRate)).ThenBy(g => g.Key, StringComparer.Ordinal))
        {
            var comps = new Dictionary<char, float[]>();
            var aligned = pick;
            double fs = 0;
            int lead = 0, length = 0;
            foreach (var tr in group.OrderBy(t => t.Channel, StringComparer.Ordinal))
            {
                if (phase == Phase.P ? !tr.IsVertical : tr.IsVertical) continue;
                if (comps.ContainsKey(tr.Component)) continue;
                if (fs > 0 && Math.Abs(tr.SampleRate - fs) > 1e-6) continue;
                var pad = 2.0 + 1.0 / Math.Max(0.05, s.LowHz);
                var from = tr.IndexOf(pick.AddSeconds(-before - s.MaxLagSeconds - pad));
                var to = tr.IndexOf(pick.AddSeconds(after + s.MaxLagSeconds + pad));
                if (from < 0 || to >= tr.Data.Length) continue;
                var seg = SignalProcessing.Prepare(tr.Data.AsSpan(from, to - from + 1), tr.SampleRate, s.LowHz, s.HighHz, 4, zeroPhase: true);
                var lagN = (int)Math.Round(s.MaxLagSeconds * tr.SampleRate);
                var nb = (int)Math.Round(before * tr.SampleRate);
                var na = (int)Math.Round(after * tr.SampleRate);
                var i0 = tr.IndexOf(pick) - from - nb - lagN;
                var n = nb + na + 2 * lagN + 1;
                if (i0 < 0 || i0 + n > seg.Length) continue;
                comps[tr.Component] = seg.AsSpan(i0, n).ToArray();
                // The window is cut on the sample nearest the pick: the lag is measured from that sample's time.
                aligned = tr.TimeOf(tr.IndexOf(pick));
                fs = tr.SampleRate;
                lead = nb + lagN;
                length = nb + na + 1;
            }
            if (comps.Count == 0) continue;
            return new Window { Instrument = group.Key, SampleRate = fs, Components = comps, Pick = aligned, Lead = lead, Length = length };
        }
        return null;
    }

    /// <summary>
    /// Normalised cross-correlation of A's template with B's window, summed over the components both have:
    /// c(τ) = Σ a(t) b(t + τ) / √(Σ a² Σ b²(· + τ)), t relative to each pick. The lag returned is τ at the
    /// maximum (s, refined by a parabola): B's arrival is at its pick + τ in A's convention, so
    /// t_A − t_B = (pick_A − pick_B) − τ.
    /// </summary>
    internal static (double Lag, double Coefficient)? Correlate(Window a, Window b, CrossCorrelationSettings s)
    {
        if (Math.Abs(a.SampleRate - b.SampleRate) > 1e-6 || a.Length != b.Length || a.Lead != b.Lead) return null;
        var comps = a.Components.Keys.Where(b.Components.ContainsKey).OrderBy(c => c).ToArray();
        if (comps.Length == 0) return null;
        var lagN = (a.Components[comps[0]].Length - a.Length) / 2;
        var start = lagN; // template: samples [lagN, lagN + Length) of A's window
        double ea = 0;
        foreach (var c in comps)
        {
            var x = a.Components[c];
            for (var i = 0; i < a.Length; i++) ea += x[start + i] * (double)x[start + i];
        }
        if (ea <= 0) return null;
        var cc = new double[2 * lagN + 1];
        for (var k = -lagN; k <= lagN; k++)
        {
            double num = 0, eb = 0;
            foreach (var c in comps)
            {
                var x = a.Components[c];
                var y = b.Components[c];
                for (var i = 0; i < a.Length; i++)
                {
                    var v = y[start + k + i];
                    num += x[start + i] * (double)v;
                    eb += v * (double)v;
                }
            }
            cc[k + lagN] = eb > 0 ? num / Math.Sqrt(ea * eb) : 0;
        }
        var best = 0;
        for (var k = 1; k < cc.Length; k++) if (cc[k] > cc[best]) best = k;
        if (best == 0 || best == cc.Length - 1) return null; // at the edge of the search range
        var coefficient = cc[best];
        if (coefficient < s.MinCoefficient) return null;
        // Parabola through the maximum and its neighbours.
        double y0 = cc[best - 1], y1 = cc[best], y2 = cc[best + 1];
        var denom = y0 - 2 * y1 + y2;
        var frac = denom < 0 ? 0.5 * (y0 - y2) / denom : 0;
        frac = Math.Clamp(frac, -0.5, 0.5);
        var peak = y1 - 0.25 * (y0 - y2) * frac;
        return ((best - lagN + frac) / a.SampleRate, Math.Min(1, peak));
    }

    // ─── Storage ─────────────────────────────────────────────────────────

    /// <summary>Default file name of the cross-correlation times.</summary>
    public const string FileName = "dt_cc.csv";

    private const string Header = "event_a,event_b,station,phase,arrival_difference_s,lag_s,cc";

    public static void Write(string path, IEnumerable<CcMeasurement> list)
    {
        var tmp = path + ".tmp";
        using (var w = new StreamWriter(tmp))
        {
            w.WriteLine(Header);
            foreach (var m in list)
                w.WriteLine(string.Create(CultureInfo.InvariantCulture,
                    $"{m.EventA},{m.EventB},{m.StationId},{m.Phase},{m.ArrivalDifference:R},{m.Lag:0.######},{m.Coefficient:0.####}"));
        }
        File.Move(tmp, path, overwrite: true);
    }

    public static List<CcMeasurement> Read(string path)
    {
        var list = new List<CcMeasurement>();
        if (!File.Exists(path)) return list;
        foreach (var line in File.ReadLines(path).Skip(1))
        {
            var f = line.Split(',');
            if (f.Length < 7 || !Enum.TryParse<Phase>(f[3], out var ph)) continue;
            list.Add(new CcMeasurement(f[0], f[1], f[2], ph, double.Parse(f[4], CultureInfo.InvariantCulture),
                double.Parse(f[5], CultureInfo.InvariantCulture), double.Parse(f[6], CultureInfo.InvariantCulture)));
        }
        return list;
    }
}
