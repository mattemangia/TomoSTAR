// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using TomoStar.Core.Forward;
using TomoStar.Core.Geo;
using TomoStar.Core.Model;
using TomoStar.Core.Tomography;

namespace TomoStar.Core.Location;

/// <summary>The relocation of one event: where it started, where it ended, and how well.</summary>
public sealed record RelocatedEvent(EventRecord Start, Hypocentre Result, int Used, int Rejected, int Iterations, string Message)
{
    /// <summary>The relocation succeeded (enough phases, a finite solution).</summary>
    public bool Ok => Message == "ok";
}

/// <summary>
/// Hypocentre relocation in a fixed velocity model, 1-D or 3-D, in two ways.
///
/// <b>Absolute</b> (<see cref="Absolute"/>): every event is located on its own by
/// <see cref="EventLocator"/>, a grid search on the travel-time tables of the stations (misfit L1
/// with the weighted-median origin time, robust to single wrong picks) refined by Geiger's
/// linearised least squares with Levenberg-Marquardt damping, with the error ellipse, depth and time
/// errors from the scaled covariance and outlier rejection. The tables are computed once per station
/// and phase with fast marching (CPU) or fast sweeping (OpenCL) on the forward grid, so the cost per
/// event is small and events are located in parallel.
///
/// <b>Double difference</b> (<see cref="DoubleDifference"/>): the method of Waldhauser &amp; Ellsworth
/// (2000, BSSA 90(6), 1353-1368) in a 3-D model as in tomoDD (Zhang &amp; Thurber 2003): differential
/// times of neighbouring events at common stations are inverted together with the absolute times for
/// the four hypocentral parameters of all events at once, velocities held fixed. The travel-time
/// tomography engine is used with its velocity unknowns switched off, so the rays are the curved
/// fast-marching rays of the model, recomputed at each iteration, and the weighting schedule of the
/// iteration sets (absolute times first, then differential times reweighted by residual and by
/// separation) is the one of the double-difference tomography.
/// </summary>
public static class Relocation
{
    /// <summary>
    /// Locates every event of the catalogue independently. The model is given on
    /// <paramref name="grid"/> (a 1-D model sampled on it, or the result of a tomography); the tables
    /// are computed on the grid refined <paramref name="refinement"/> times. Events keep their
    /// identifiers; the catalogue itself is not changed.
    /// </summary>
    public static List<RelocatedEvent> Absolute(Catalogue catalogue, SphericalGrid grid, double[] vp, double[] vs, int refinement,
        LocatorSettings settings, string workFolder, bool useOpenCl, Action<string>? log = null, IProgress<(double, string)>? progress = null,
        CancellationToken ct = default)
    {
        var stations = catalogue.StationIndex();
        var top = grid.DepthKm[0];
        // Tables only for the stations and phases that have picks.
        var needed = new HashSet<(string Station, Phase Phase)>();
        foreach (var e in catalogue.Events)
        foreach (var p in e.Picks)
            if (p.Usable && stations.ContainsKey(p.StationId)) needed.Add((stations[p.StationId].Id, p.Phase));
        var sources = needed.Select(x => x.Station).Distinct(StringComparer.OrdinalIgnoreCase)
            .Select(id => stations[id])
            .Where(s => grid.Contains(s.Lon, s.Lat, Math.Max(top, s.DepthKm)))
            .Select(s => new TableSource(s.Id, s.Lon, s.Lat, Math.Max(top, s.DepthKm))).ToList();
        var outside = needed.Select(x => x.Station).Distinct(StringComparer.OrdinalIgnoreCase).Count() - sources.Count;
        if (outside > 0) log?.Invoke($"Relocation: {outside} stations outside the grid are not used (the grid holds the tables).");
        var phases = needed.Select(x => x.Phase).Distinct().ToList();
        var fm = ForwardModel.FromInversionGrid(grid, vp, vs, grid.Definition.Refined(refinement));
        using var gpu = useOpenCl ? EikonalOpenCl.TryCreate(log) : null;
        using var tables = TravelTimeTableSet.Build(fm, sources, phases, Path.Combine(workFolder, "tables"), gpu,
            progress == null ? null : new Progress<(double, string)>(p => progress.Report((0.5 * p.Item1, p.Item2))), log, ct);
        var locator = new EventLocator(tables, stations, settings);
        var results = new RelocatedEvent?[catalogue.Events.Count];
        var done = 0;
        Parallel.For(0, catalogue.Events.Count, Compute.ComputeSettings.Options(ct), i =>
        {
            var ev = catalogue.Events[i];
            if (ev.Fixed)
            {
                results[i] = new RelocatedEvent(ev, new Hypocentre { OriginTime = ev.OriginTime, Lon = ev.Lon, Lat = ev.Lat, DepthKm = ev.DepthKm }, 0, 0, 0, "fixed");
                return;
            }
            LocationResult r;
            try { r = locator.Locate(ev, ct); }
            catch (Exception ex) when (ex is not OperationCanceledException)
            {
                r = new LocationResult(new Hypocentre { OriginTime = ev.OriginTime, Lon = ev.Lon, Lat = ev.Lat, DepthKm = ev.DepthKm }, 0, 0, 0, ex.Message);
            }
            results[i] = new RelocatedEvent(ev, r.Hypocentre, r.Used, r.Rejected, r.Iterations, r.Message);
            var d = Interlocked.Increment(ref done);
            if (d % 50 == 0) progress?.Report((0.5 + 0.5 * d / catalogue.Events.Count, $"Located {d}/{catalogue.Events.Count}"));
        });
        var list = results.Select(x => x!).ToList();
        var ok = list.Count(x => x.Ok);
        log?.Invoke($"Relocation: {ok} of {list.Count} events located" +
                    (ok > 0 ? string.Create(System.Globalization.CultureInfo.InvariantCulture,
                        $", median RMS {Median(list.Where(x => x.Ok).Select(x => x.Result.Rms)):0.000} s, median shift {Median(list.Where(x => x.Ok).Select(Shift)):0.00} km.") : "."));
        return list;
    }

    /// <summary>
    /// Double-difference relocation of all events at once in a fixed model (see the class summary).
    /// <paramref name="settings"/> gives the forward method, the outlier rules, the hypocentre damping
    /// and the double-difference pairs and iteration sets; its velocity and station-term options are
    /// ignored. Returns the relocated events (in the order of the observation set) and the run of the
    /// engine, whose residuals and statistics can be written like those of a tomography.
    /// </summary>
    public static (List<RelocatedEvent> Events, TomographyResult Run) DoubleDifference(Catalogue catalogue, SphericalGrid grid, double[] vp, double[] vs,
        TomographySettings settings, string workFolder, int minPhasesPerEvent = 4, OutsideData outside = OutsideData.Exclude,
        VelocityModel1D? background = null, Action<string>? log = null, IProgress<(double, string)>? progress = null, CancellationToken ct = default)
    {
        var s = settings.Clone();
        s.InvertP = false;
        s.InvertS = false;
        s.StationCorrections = false;
        s.JointHypocentres = true;
        s.Layered = false;
        s.Adaptive.Enabled = false;
        s.Lattice.Enabled = false;
        s.DoubleDifference.Enabled = true;
        var data = ObservationSet.FromCatalogue(catalogue, grid, true, true, minPhasesPerEvent, outside);
        if (data.Events.Count < 2) throw new InvalidOperationException("Double-difference relocation needs at least two events with enough phases.");
        log?.Invoke("Double-difference relocation: " + s.DoubleDifference.Describe() + "; velocities fixed.");
        var run = new TravelTimeTomography(grid, s, workFolder, log) { Background = background }.Run(data, vp, vs, progress, ct);
        var byId = catalogue.EventIndex();
        var count = new int[data.Events.Count];
        var rejected = new int[data.Events.Count];
        var sum = new double[data.Events.Count];
        foreach (var a in data.Arrivals)
        {
            if (a.Rejected) { rejected[a.Event]++; continue; }
            if (double.IsNaN(a.Residual)) continue;
            count[a.Event]++;
            sum[a.Event] += a.Residual * a.Residual;
        }
        var events = data.Events.Select((e, i) => new RelocatedEvent(byId[e.Id],
            new Hypocentre
            {
                OriginTime = e.Reference.AddSeconds(e.T0), Lon = e.Lon, Lat = e.Lat, DepthKm = e.DepthKm,
                Rms = count[i] > 0 ? Math.Sqrt(sum[i] / count[i]) : double.NaN, PhaseCount = count[i]
            }, count[i], rejected[i], run.Iterations.Count - 1, "ok")).ToList();
        return (events, run);
    }

    /// <summary>Horizontal and vertical distance moved, km.</summary>
    public static double Shift(RelocatedEvent r)
    {
        var h = GeoMath.SurfaceDistanceKm(r.Start.Lon, r.Start.Lat, r.Result.Lon, r.Result.Lat);
        var z = r.Result.DepthKm - r.Start.DepthKm;
        return Math.Sqrt(h * h + z * z);
    }

    /// <summary>
    /// A copy of the catalogue with the relocated hypocentres (events that failed keep their start),
    /// ready for the next step: the tomography starts from the relocated events.
    /// </summary>
    public static Catalogue Apply(Catalogue catalogue, IEnumerable<RelocatedEvent> relocated)
    {
        var map = relocated.Where(r => r.Ok).ToDictionary(r => r.Start.Id);
        var c = new Catalogue();
        c.Stations.AddRange(catalogue.Stations);
        c.TStar.AddRange(catalogue.TStar);
        c.Warnings.AddRange(catalogue.Warnings);
        foreach (var e in catalogue.Events)
        {
            var copy = new EventRecord { Id = e.Id, OriginTime = e.OriginTime, Lon = e.Lon, Lat = e.Lat, DepthKm = e.DepthKm, Magnitude = e.Magnitude, Fixed = e.Fixed };
            if (map.TryGetValue(e.Id, out var r))
                (copy.OriginTime, copy.Lon, copy.Lat, copy.DepthKm) = (r.Result.OriginTime, r.Result.Lon, r.Result.Lat, r.Result.DepthKm);
            copy.Picks.AddRange(e.Picks);
            c.Events.Add(copy);
        }
        return c;
    }

    private static double Median(IEnumerable<double> values)
    {
        var a = values.Where(double.IsFinite).OrderBy(x => x).ToArray();
        return a.Length == 0 ? double.NaN : a.Length % 2 == 1 ? a[a.Length / 2] : 0.5 * (a[a.Length / 2 - 1] + a[a.Length / 2]);
    }
}
