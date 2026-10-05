// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using TomoStar.Core.Geo;
using TomoStar.Core.Model;

namespace TomoStar.Core.Tomography;

/// <summary>The state of one event during an inversion: position and origin-time shift.</summary>
public sealed class EventState
{
    public required string Id { get; init; }
    public double Lon { get; set; }
    public double Lat { get; set; }
    public double DepthKm { get; set; }

    /// <summary>Origin time, s relative to <see cref="Reference"/>.</summary>
    public double T0 { get; set; }

    public DateTime Reference { get; init; }

    /// <summary>Held fixed (e.g. a shot, or a synthetic test that keeps hypocentres).</summary>
    public bool Fixed { get; set; }

    /// <summary>
    /// Outside the inversion grid: the position stays at the catalogue location (the network seldom
    /// surrounds such an event), the origin time is still solved for.
    /// </summary>
    public bool Outside { get; set; }

    public EventState Clone() => (EventState)MemberwiseClone();
}

public sealed class StationState
{
    public required string Id { get; init; }
    public double Lon { get; init; }
    public double Lat { get; init; }
    public double DepthKm { get; init; }
    public double CorrectionP { get; set; }
    public double CorrectionS { get; set; }

    /// <summary>Outside the inversion grid (used because its rays cross it).</summary>
    public bool Outside { get; init; }

    public double Correction(Phase p) => p == Phase.P ? CorrectionP : CorrectionS;
}

/// <summary>One arrival: seconds after the event reference time, with its one-sigma error.</summary>
public sealed class Observation
{
    public required int Event { get; init; }
    public required int Station { get; init; }
    public required Phase Phase { get; init; }
    public double Time { get; set; }
    public double Sigma { get; set; } = 0.1;

    public bool Rejected { get; set; }

    /// <summary>Why the arrival was left out of the last iteration (empty when used).</summary>
    public string RejectReason { get; set; } = "";

    public double Residual { get; set; } = double.NaN;
    public double InitialResidual { get; set; } = double.NaN;
}

/// <summary>
/// The data of an inversion, detached from the project: events, stations and arrivals with index
/// references. Real data come from the project picks; synthetic tests build the same structure
/// with computed times, so the engine never knows which it is inverting.
/// </summary>
public sealed class ObservationSet
{
    public List<EventState> Events { get; } = [];
    public List<StationState> Stations { get; } = [];
    public List<Observation> Arrivals { get; } = [];
    public List<string> Warnings { get; } = [];

    public ObservationSet CloneGeometry()
    {
        var o = new ObservationSet();
        o.Events.AddRange(Events.Select(e => e.Clone()));
        o.Stations.AddRange(Stations.Select(s => new StationState { Id = s.Id, Lon = s.Lon, Lat = s.Lat, DepthKm = s.DepthKm, Outside = s.Outside }));
        o.Arrivals.AddRange(Arrivals.Select(a => new Observation { Event = a.Event, Station = a.Station, Phase = a.Phase, Time = a.Time, Sigma = a.Sigma }));
        return o;
    }

    /// <summary>
    /// The <paramref name="maxEvents"/> events with the most arrivals (ties: in their order), with
    /// their arrivals re-indexed; all stations are kept. Used to preview a regularisation quickly.
    /// </summary>
    public ObservationSet Subset(int maxEvents)
    {
        var counts = new int[Events.Count];
        foreach (var a in Arrivals) counts[a.Event]++;
        var keep = Enumerable.Range(0, Events.Count).Where(i => counts[i] > 0)
            .OrderByDescending(i => counts[i]).ThenBy(i => i).Take(Math.Max(1, maxEvents)).OrderBy(i => i).ToList();
        var map = new Dictionary<int, int>();
        var o = new ObservationSet();
        foreach (var i in keep)
        {
            map[i] = o.Events.Count;
            o.Events.Add(Events[i].Clone());
        }
        o.Stations.AddRange(Stations.Select(s => new StationState { Id = s.Id, Lon = s.Lon, Lat = s.Lat, DepthKm = s.DepthKm, Outside = s.Outside }));
        foreach (var a in Arrivals)
            if (map.TryGetValue(a.Event, out var e))
                o.Arrivals.Add(new Observation { Event = e, Station = a.Station, Phase = a.Phase, Time = a.Time, Sigma = a.Sigma });
        return o;
    }

    /// <summary>Every event and station position, for sizing the forward domain.</summary>
    public IEnumerable<(double Lon, double Lat, double DepthKm)> Points() =>
        Events.Select(e => (e.Lon, e.Lat, e.DepthKm)).Concat(Stations.Select(s => (s.Lon, s.Lat, s.DepthKm)));

    public bool HasOutsidePoints => Events.Any(e => e.Outside) || Stations.Any(s => s.Outside);

    /// <summary>
    /// Gathers every usable pick of a catalogue. Events and stations inside the grid are always used.
    /// Outside it, a pick is used only when <paramref name="outside"/> allows it and the straight
    /// path from event to station crosses the grid; the rest are left out with a warning, never
    /// silently moved inside the grid. Each event keeps one arrival per station and phase
    /// (<see cref="PickSelection.Best"/>), and events with fewer than
    /// <paramref name="minPhasesPerEvent"/> arrivals are left out.
    /// </summary>
    public static ObservationSet FromCatalogue(Catalogue catalogue, SphericalGrid grid, bool includeP, bool includeS, int minPhasesPerEvent = 4,
        OutsideData outside = OutsideData.Exclude, AutomaticPicks automatic = AutomaticPicks.IfFewAnalystPicks)
    {
        var set = new ObservationSet();
        set.Warnings.AddRange(catalogue.Warnings);
        var stationIndex = new Dictionary<string, int>(StringComparer.OrdinalIgnoreCase);
        var outsideStations = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
        var unknownStations = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
        var stations = catalogue.StationIndex();
        int outsideEvents = 0, fewPhases = 0, missed = 0;
        var top = grid.DepthKm[0];

        foreach (var ev in catalogue.Events)
        {
            var evInside = grid.Contains(ev.Lon, ev.Lat, ev.DepthKm);
            if (!evInside && outside == OutsideData.Exclude) { outsideEvents++; continue; }
            var picks = PickSelection.Best(ev.Picks, automatic).Where(p => p.Phase == Phase.P ? includeP : includeS).ToList();
            var usable = new List<(PickRecord Pick, int Station)>();
            foreach (var p in picks)
            {
                if (!stations.TryGetValue(p.StationId, out var st)) { unknownStations.Add(p.StationId); continue; }
                var depth = Math.Max(top, st.DepthKm);
                var stInside = grid.Contains(st.Lon, st.Lat, depth);
                if (!stInside && outside == OutsideData.Exclude) { outsideStations.Add(st.Id); continue; }
                if ((!stInside || !evInside) && !PathCrosses(grid, ev.Lon, ev.Lat, ev.DepthKm, st.Lon, st.Lat, depth))
                {
                    missed++;
                    continue;
                }
                if (!stationIndex.TryGetValue(st.Id, out var si))
                {
                    si = set.Stations.Count;
                    stationIndex[st.Id] = si;
                    set.Stations.Add(new StationState
                    {
                        Id = st.Id, Lon = st.Lon, Lat = st.Lat,
                        // A station above the top of the grid sits on it (warned below).
                        DepthKm = depth, Outside = !stInside,
                        CorrectionP = st.CorrectionP, CorrectionS = st.CorrectionS
                    });
                    if (st.DepthKm < top) set.Warnings.Add($"Station {st.Id} (elevation {st.ElevationM:0} m) is above the grid top ({-top * 1000:0} m); placed on the top surface.");
                }
                usable.Add((p, si));
            }
            if (usable.Count < minPhasesPerEvent) { fewPhases++; continue; }
            var ei = set.Events.Count;
            set.Events.Add(new EventState { Id = ev.Id, Lon = ev.Lon, Lat = ev.Lat, DepthKm = ev.DepthKm, Reference = ev.OriginTime, Outside = !evInside, Fixed = ev.Fixed });
            // One arrival per station and phase (already reduced by PickSelection.Best).
            foreach (var (pick, station) in usable)
                set.Arrivals.Add(new Observation
                {
                    Event = ei, Station = station, Phase = pick.Phase,
                    Time = (pick.Time - ev.OriginTime).TotalSeconds,
                    Sigma = Math.Max(0.01, pick.Sigma)
                });
        }
        if (unknownStations.Count > 0) set.Warnings.Add($"{unknownStations.Count} stations with picks are missing from the station list and their picks were left out: {string.Join(", ", unknownStations.Take(10))}{(unknownStations.Count > 10 ? ", ..." : "")}");
        if (outsideEvents > 0) set.Warnings.Add($"{outsideEvents} events lie outside the grid and were left out.");
        if (outsideStations.Count > 0) set.Warnings.Add($"{outsideStations.Count} stations lie outside the grid and were left out: {string.Join(", ", outsideStations.Take(10))}{(outsideStations.Count > 10 ? ", ..." : "")}");
        if (missed > 0) set.Warnings.Add($"{missed} picks between an event and a station outside the grid have a path that does not cross it and were left out.");
        var oe = set.Events.Count(e => e.Outside);
        var os = set.Stations.Count(s => s.Outside);
        if (oe + os > 0) set.Warnings.Add($"Used from outside the grid because their rays cross it: {oe} events (position fixed, origin time solved) and {os} stations.");
        if (fewPhases > 0) set.Warnings.Add($"{fewPhases} events have fewer than {minPhasesPerEvent} usable phases and were left out.");
        return set;
    }

    /// <summary>
    /// Does the straight chord between two points pass through the grid volume? Sampled every
    /// half node spacing. Curved rays bend, but a chord that misses the box by more than the
    /// bending carries no information on the model and is not worth tracing.
    /// </summary>
    public static bool PathCrosses(SphericalGrid grid, double lon1, double lat1, double dep1, double lon2, double lat2, double dep2)
    {
        if (grid.Contains(lon1, lat1, dep1) || grid.Contains(lon2, lat2, dep2)) return true;
        var a = GeoMath.ToCartesian(lon1, lat1, dep1);
        var b = GeoMath.ToCartesian(lon2, lat2, dep2);
        var len = (b - a).Length;
        var n = Math.Max(2, (int)Math.Ceiling(len / Math.Max(0.1, 0.5 * grid.MinSpacingKm())));
        for (var i = 1; i < n; i++)
        {
            var p = GeoMath.FromCartesian(a + (b - a) * ((double)i / n));
            // The chord of a long path dips below the surface; a point just above the top counts.
            if (grid.Contains(p.Lon, p.Lat, Math.Max(p.DepthKm, grid.DepthKm[0]))) return true;
        }
        return false;
    }
}

/// <summary>How the smoothing (Laplacian) weighs the grid axes.</summary>
public enum SmoothingScale
{
    /// <summary>
    /// Second differences in physical units (each axis weighted by 1/Δ², Δ in km, normalised to the
    /// horizontal spacing): the same structure costs the same whatever the node spacing. On a grid
    /// of 10 km × 2 km cells this is what keeps the smoothing from preferring horizontal layers.
    /// </summary>
    Kilometres,

    /// <summary>
    /// Second differences between neighbouring nodes, whatever their distance (the usual choice of
    /// node-based tomography codes). With cells longer than they are high, a lateral gradient
    /// costs (Δh/Δz)² times more than a vertical one of the same length: models come out layered.
    /// </summary>
    Nodes
}

/// <summary>How the roughness of the model is penalised.</summary>
public enum SmoothingMethod
{
    /// <summary>Second differences, least squares: smooth models (the default).</summary>
    Laplacian,

    /// <summary>First differences, least squares: the flattest model; gradients cost, curvature does not.</summary>
    Gradient,

    /// <summary>
    /// First differences in the L1 sense (total variation, Rudin et al. 1992), by iteratively
    /// reweighted least squares: blocky models, sharp contacts between nearly uniform units.
    /// </summary>
    TotalVariation,

    /// <summary>
    /// Second differences weighted by 1/(1 + |∇m|²/ε²) from the current model (anisotropic
    /// diffusion, Perona &amp; Malik 1990; edge-preserving regularisation, Charbonnier et al. 1997):
    /// smooth within units, little smoothing across the contacts the data ask for.
    /// </summary>
    EdgePreserving
}

/// <summary>Regularisation rows shared by the velocity and the attenuation tomography.</summary>
public static class Regularization
{
    /// <summary>
    /// Roughness rows of the chosen <paramref name="method"/> on a node grid block (see
    /// <see cref="AddLaplacian"/> for the conventions). For the reweighted methods the weights
    /// come from <paramref name="weightModel"/> (the current deviation; null or all zero: plain
    /// least squares, as at the first iteration from a 1-D start): a row whose roughness r exceeds
    /// the scale ε gets the weight ε/√(r² + ε²) (total variation) or 1/(1 + r²/ε²) (edge-preserving), with
    /// ε = <paramref name="edgeScale"/> × the robust spread (1.4826 × median |r|) of the rows; r is
    /// the gradient magnitude at the node for both.
    /// </summary>
    public static void AddSmoothing(Numerics.CsrBuilder b, SphericalGrid g, int offset, double lambda, double verticalWeight,
        double[]? current, SmoothingMethod method, double edgeScale = 1, double[]? weightModel = null, SmoothingScale scale = SmoothingScale.Kilometres)
    {
        if (lambda <= 0) return;
        switch (method)
        {
            case SmoothingMethod.Laplacian:
                AddLaplacian(b, g, offset, lambda, verticalWeight, current, scale: scale);
                return;
            case SmoothingMethod.EdgePreserving:
            {
                // Perona-Malik: the Laplacian of each node weighted down where the gradient is large,
                // so the model is smoothed within units and not across their contacts.
                var r = weightModel != null ? GradientMagnitude(g, weightModel, verticalWeight, scale) : null;
                var w = Reweight(r, edgeScale, cauchy: true);
                AddLaplacian(b, g, offset, lambda, verticalWeight, current, scale: scale, rowWeights: w);
                return;
            }
            default:
            {
                double[]? w = null;
                if (method == SmoothingMethod.TotalVariation && weightModel != null)
                    w = Reweight(GradientMagnitude(g, weightModel, verticalWeight, scale), edgeScale, cauchy: false);
                AddGradient(b, g, offset, lambda, verticalWeight, current, scale, w);
                return;
            }
        }
    }

    /// <summary>√(IRLS weight) per node from its roughness (null: none, least squares).</summary>
    private static double[]? Reweight(double[]? r, double edgeScale, bool cauchy)
    {
        if (r == null) return null;
        var abs = r.Select(Math.Abs).Where(x => x > 0).OrderBy(x => x).ToArray();
        if (abs.Length < 8) return null;
        var eps = Math.Max(1e-12, Math.Max(0.01, edgeScale) * 1.4826 * abs[abs.Length / 2]);
        return r.Select(x => cauchy ? 1 / Math.Sqrt(1 + x * x / (eps * eps)) : Math.Sqrt(eps / Math.Sqrt(x * x + eps * eps))).ToArray();
    }

    /// <summary>The weighted Laplacian of a node field (the residual of each smoothing row on it).</summary>
    public static double[] LaplacianOf(SphericalGrid g, double[] x, double verticalWeight, SmoothingScale scale = SmoothingScale.Kilometres)
    {
        var (wx, wy, wz) = AxisWeights(g, verticalWeight, scale);
        var lap = new double[g.Count];
        Parallel.For(0, g.Nz, k =>
        {
            for (var j = 0; j < g.Ny; j++)
            for (var i = 0; i < g.Nx; i++)
            {
                double v = 0, diag = 0;
                void Add(int di, int dj, int dk, double w)
                {
                    var ii = i + di; var jj = j + dj; var kk = k + dk;
                    if (!g.InRange(ii, jj, kk)) return;
                    v -= w * x[g.Index(ii, jj, kk)];
                    diag += w;
                }
                Add(-1, 0, 0, wx); Add(1, 0, 0, wx); Add(0, -1, 0, wy); Add(0, 1, 0, wy);
                Add(0, 0, -1, wz); Add(0, 0, 1, wz);
                var n = g.Index(i, j, k);
                lap[n] = v + diag * x[n];
            }
        });
        return lap;
    }

    /// <summary>|∇x| at each node from forward differences, the axes weighted as in <see cref="AxisWeights"/> (square roots).</summary>
    public static double[] GradientMagnitude(SphericalGrid g, double[] x, double verticalWeight, SmoothingScale scale = SmoothingScale.Kilometres)
    {
        var (wx, wy, wz) = AxisWeights(g, verticalWeight, scale);
        var m = new double[g.Count];
        Parallel.For(0, g.Nz, k =>
        {
            for (var j = 0; j < g.Ny; j++)
            for (var i = 0; i < g.Nx; i++)
            {
                var n = g.Index(i, j, k);
                double s = 0;
                if (i + 1 < g.Nx) s += wx * Math.Pow(x[g.Index(i + 1, j, k)] - x[n], 2);
                if (j + 1 < g.Ny) s += wy * Math.Pow(x[g.Index(i, j + 1, k)] - x[n], 2);
                if (k + 1 < g.Nz) s += wz * Math.Pow(x[g.Index(i, j, k + 1)] - x[n], 2);
                m[n] = Math.Sqrt(s);
            }
        });
        return m;
    }

    /// <summary>
    /// First-difference rows: for every node and axis with a forward neighbour,
    /// λ·√wₐ·(xₙ₊ₐ − xₙ) = −λ·√wₐ·(same on <paramref name="current"/>), each row times the node's
    /// weight when <paramref name="rowWeights"/> is given.
    /// </summary>
    public static void AddGradient(Numerics.CsrBuilder b, SphericalGrid g, int offset, double lambda, double verticalWeight,
        double[]? current, SmoothingScale scale = SmoothingScale.Kilometres, double[]? rowWeights = null)
    {
        if (lambda <= 0) return;
        var (wx, wy, wz) = AxisWeights(g, verticalWeight, scale);
        Span<int> cols = stackalloc int[2];
        Span<double> vals = stackalloc double[2];
        for (var k = 0; k < g.Nz; k++)
        for (var j = 0; j < g.Ny; j++)
        for (var i = 0; i < g.Nx; i++)
        {
            var n = g.Index(i, j, k);
            var rw = rowWeights?[n] ?? 1;
            foreach (var (di, dj, dk, w) in new[] { (1, 0, 0, wx), (0, 1, 0, wy), (0, 0, 1, wz) })
            {
                if (!g.InRange(i + di, j + dj, k + dk) || w <= 0) continue;
                var m = g.Index(i + di, j + dj, k + dk);
                var a = Math.Sqrt(w);
                cols[0] = offset + n; vals[0] = -a;
                cols[1] = offset + m; vals[1] = a;
                var rhs = current != null ? -a * (current[m] - current[n]) : 0;
                b.AddRow(cols, vals, rhs, lambda * rw);
            }
        }
    }

    /// <summary>
    /// Weights of the three axes in the Laplacian. The vertical weight states how much less (below 1)
    /// or more a vertical second difference costs than a horizontal one of the same length.
    /// </summary>
    public static (double Lon, double Lat, double Depth) AxisWeights(SphericalGrid g, double verticalWeight, SmoothingScale scale)
    {
        if (scale == SmoothingScale.Nodes) return (1, 1, verticalWeight);
        // Spacings at the centre of the grid near the surface.
        var (hx, hy, hz) = g.Spacing(g.Nx / 2, g.Ny / 2, 0);
        var h = Math.Sqrt(hx * hy);
        return (Math.Pow(h / hx, 2), Math.Pow(h / hy, 2), verticalWeight * Math.Pow(h / hz, 2));
    }

    /// <summary>
    /// Second-difference (Laplacian) rows on a node grid block starting at column
    /// <paramref name="offset"/>: for every node, λ·Σ_axes wₐ(2xₙ − xₙ₋ₐ − xₙ₊ₐ) = −λ·(same on
    /// <paramref name="current"/>), i.e. the roughness of the TOTAL deviation from the reference is
    /// penalised, not only that of the update, so roughness cannot accumulate over iterations.
    /// The axis weights come from <see cref="AxisWeights"/>. At the edges the missing neighbour is
    /// dropped (one-sided).
    /// </summary>
    public static void AddLaplacian(Numerics.CsrBuilder b, SphericalGrid g, int offset, double lambda, double verticalWeight,
        double[]? current, bool[]? active = null, SmoothingScale scale = SmoothingScale.Kilometres, double[]? rowWeights = null)
    {
        if (lambda <= 0) return;
        var (wx, wy, wz) = AxisWeights(g, verticalWeight, scale);
        var cols = new int[7];
        var vals = new double[7];
        for (var k = 0; k < g.Nz; k++)
        for (var j = 0; j < g.Ny; j++)
        for (var i = 0; i < g.Nx; i++)
        {
            var n = g.Index(i, j, k);
            if (active != null && !active[n]) continue;
            var c = 0;
            double diag = 0;
            void Add(int di, int dj, int dk, double w)
            {
                var ii = i + di; var jj = j + dj; var kk = k + dk;
                if (!g.InRange(ii, jj, kk)) return;
                cols[c] = offset + g.Index(ii, jj, kk);
                vals[c] = -w;
                c++;
                diag += w;
            }
            Add(-1, 0, 0, wx); Add(1, 0, 0, wx); Add(0, -1, 0, wy); Add(0, 1, 0, wy);
            Add(0, 0, -1, wz); Add(0, 0, 1, wz);
            cols[c] = offset + n;
            vals[c] = diag;
            c++;
            double rhs = 0;
            if (current != null)
            {
                for (var q = 0; q < c; q++) rhs -= vals[q] * current[cols[q] - offset];
            }
            b.AddRow(cols.AsSpan(0, c), vals.AsSpan(0, c), rhs, lambda * (rowWeights?[n] ?? 1));
        }
    }

    /// <summary>
    /// Second differences in depth on a block of <paramref name="layers"/> columns, one per layer:
    /// the smoothing of a layered (1-D) model. Same convention as <see cref="AddLaplacian"/>.
    /// </summary>
    public static void AddLaplacian1D(Numerics.CsrBuilder b, int layers, int offset, double lambda, double[]? current)
    {
        if (lambda <= 0) return;
        for (var k = 0; k < layers; k++)
        {
            var cols = new List<int>();
            var vals = new List<double>();
            double diag = 0;
            foreach (var d in new[] { -1, 1 })
            {
                if (k + d < 0 || k + d >= layers) continue;
                cols.Add(offset + k + d);
                vals.Add(-1);
                diag += 1;
            }
            if (diag == 0) continue;
            cols.Add(offset + k);
            vals.Add(diag);
            double rhs = 0;
            if (current != null) for (var q = 0; q < cols.Count; q++) rhs -= vals[q] * current[cols[q] - offset];
            var order = cols.Select((c, i) => (c, v: vals[i])).OrderBy(x => x.c).ToArray();
            b.AddRow(order.Select(x => x.c).ToArray(), order.Select(x => x.v).ToArray(), rhs, lambda);
        }
    }

    /// <summary>λ·xₙ = 0 for every column in [offset, offset+count): damping of the update.</summary>
    public static void AddDamping(Numerics.CsrBuilder b, int offset, int count, double lambda)
    {
        if (lambda <= 0) return;
        Span<int> c = stackalloc int[1];
        Span<double> v = stackalloc double[1];
        for (var n = 0; n < count; n++)
        {
            c[0] = offset + n;
            v[0] = 1;
            b.AddRow(c, v, 0, lambda);
        }
    }

    /// <summary>Roughness ‖L x‖ of a node field (for L-curves).</summary>
    public static double Roughness(SphericalGrid g, double[] x, double verticalWeight, SmoothingScale scale = SmoothingScale.Kilometres)
    {
        var (wx, wy, wz) = AxisWeights(g, verticalWeight, scale);
        double s = 0;
        for (var k = 0; k < g.Nz; k++)
        for (var j = 0; j < g.Ny; j++)
        for (var i = 0; i < g.Nx; i++)
        {
            double lap = 0, diag = 0;
            void Add(int di, int dj, int dk, double w)
            {
                var ii = i + di; var jj = j + dj; var kk = k + dk;
                if (!g.InRange(ii, jj, kk)) return;
                lap -= w * x[g.Index(ii, jj, kk)];
                diag += w;
            }
            Add(-1, 0, 0, wx); Add(1, 0, 0, wx); Add(0, -1, 0, wy); Add(0, 1, 0, wy);
            Add(0, 0, -1, wz); Add(0, 0, 1, wz);
            lap += diag * x[g.Index(i, j, k)];
            s += lap * lap;
        }
        return Math.Sqrt(s);
    }
}
