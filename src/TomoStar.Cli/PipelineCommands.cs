// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using System.Collections.Concurrent;
using System.Globalization;
using TomoStar.Core.Attenuation;
using TomoStar.Core.Forward;
using TomoStar.Core.Geo;
using TomoStar.Core.IO;
using TomoStar.Core.Location;
using TomoStar.Core.Model;
using TomoStar.Core.Signal;
using TomoStar.Core.Tomography;

namespace TomoStar.Cli;

/// <summary>The steps before the inversion: automatic picking, relocation and t* measurement.</summary>
public static class PipelineCommands
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>The waveform archive: --waveforms, else the waveforms of the QUIVER project given as data.</summary>
    private static WaveformArchive Archive(CommandContext c)
    {
        WaveformArchive a;
        if (c.Args.Get("waveforms") is { } dir) a = WaveformArchive.FromFolder(dir);
        else if (c.Quiver != null) a = WaveformArchive.FromQuiverProject(c.Quiver.Folder);
        else throw new UsageException("Give the waveforms with --waveforms DIR (a folder per event, or any miniSEED files).");
        a.BeforeSeconds = c.Config.Waveforms.BeforeSeconds;
        a.AfterSeconds = c.Config.Waveforms.AfterSeconds;
        return a;
    }

    // ---- pick ------------------------------------------------------------------------------------

    /// <summary>
    /// STA/LTA and AIC picks on the waveforms of every event, searched around the times the starting
    /// model predicts (travel-time tables on the grid; the whole trace for stations or events
    /// outside it). The existing picks are kept (analyst picks are preferred to automatic ones by
    /// the inversions) unless --only-automatic.
    /// </summary>
    public static int Pick(CommandContext c)
    {
        var cat = c.Catalogue;
        var archive = Archive(c);
        var s = c.Config.Picker;
        var g = new SphericalGrid(c.Grid);
        var (vp, vs, _) = c.Velocities();
        var top = g.DepthKm[0];
        var stations = cat.StationIndex();
        var sources = cat.Stations.Where(st => g.Contains(st.Lon, st.Lat, Math.Max(top, st.DepthKm)))
            .Select(st => new TableSource(st.Id, st.Lon, st.Lat, Math.Max(top, st.DepthKm))).ToList();
        var phases = s.PickS ? new List<Phase> { Phase.P, Phase.S } : [Phase.P];
        c.Log($"Picking: travel-time tables for {sources.Count} stations to predict the arrivals.");
        using var gpu = c.Config.UseOpenCl ? EikonalOpenCl.TryCreate(c.Log) : null;
        using var tables = TravelTimeTableSet.Build(ForwardModel.FromInversionGrid(g, vp, vs, g.Definition.Refined(c.Config.Tomography.ForwardRefinement)),
            sources, phases, Path.Combine(c.WorkFolder, "tables"), gpu, c.Progress(), c.Log, c.Cancel);
        var onlyAutomatic = c.Args.Flag("only-automatic");
        var automatic = new ConcurrentBag<(string Event, PickRecord Pick)>();
        var noWaveforms = 0;
        var done = 0;
        var progress = c.Progress();
        Parallel.ForEach(cat.Events, Core.Compute.ComputeSettings.Options(c.Cancel), ev =>
        {
            var traces = archive.ForEvent(ev);
            if (traces.Count == 0) { Interlocked.Increment(ref noWaveforms); return; }
            var inside = g.Contains(ev.Lon, ev.Lat, ev.DepthKm);
            foreach (var (stationId, list) in traces)
            {
                if (!stations.ContainsKey(stationId)) continue;
                DateTime? Predicted(Phase ph)
                {
                    if (!inside || tables.Get(stationId, ph) is not { } t) return null;
                    var tt = t.Time(ev.Lon, ev.Lat, ev.DepthKm);
                    return double.IsFinite(tt) ? ev.OriginTime.AddSeconds(tt) : null;
                }
                foreach (var r in AutoPicker.PickStation(list, Predicted(Phase.P), s.PickS ? Predicted(Phase.S) : null, s, ev.OriginTime))
                {
                    var trace = list.FirstOrDefault(t => t.Channel == r.Channel) ?? list[0];
                    automatic.Add((ev.Id, AutoPicker.ToPick(r, trace)));
                }
            }
            var d = Interlocked.Increment(ref done);
            progress.Report(((double)d / cat.Events.Count, $"Picking: {d}/{cat.Events.Count} events"));
        });
        var byEvent = automatic.GroupBy(x => x.Event).ToDictionary(x => x.Key, x => x.Select(y => y.Pick).ToList());
        foreach (var ev in cat.Events)
        {
            if (onlyAutomatic) ev.Picks.Clear();
            else ev.Picks.RemoveAll(p => p.Origin == PickOrigin.Automatic); // replaced by the new ones
            if (byEvent.TryGetValue(ev.Id, out var list))
                ev.Picks.AddRange(list.OrderBy(p => p.StationId, StringComparer.Ordinal).ThenBy(p => p.Phase));
        }
        var nP = automatic.Count(x => x.Pick.Phase == Phase.P);
        var nS = automatic.Count - nP;
        CatalogueWriter.WritePicks(Path.Combine(c.RunFolder, "picks.csv"), cat.Events);
        CatalogueWriter.WriteEvents(Path.Combine(c.RunFolder, "events.csv"), cat.Events);
        CatalogueWriter.WriteStations(Path.Combine(c.RunFolder, "stations.csv"), cat.Stations);
        RunWriter.WriteSummary(c.RunFolder, new
        {
            Tool = "Automatic picking", Software = "TomoSTAR", Picker = s, AutomaticP = nP, AutomaticS = nS,
            EventsWithoutWaveforms = noWaveforms, Grid = c.Grid, StartingModel = c.Model.Name
        });
        c.Log($"Picking: {nP} P and {nS} S automatic picks{(noWaveforms > 0 ? $"; {noWaveforms} events without waveforms" : "")}.");
        c.Log($"Tables stations.csv, events.csv and picks.csv in {c.RunFolder} (use it as --data of the next step).");
        RunWriter.WriteLog(c.RunFolder, c.LogLines);
        return 0;
    }

    // ---- relocate --------------------------------------------------------------------------------

    public static int Relocate(CommandContext c)
    {
        var method = (c.Args.Get("method") ?? c.Config.Relocation.Method).ToLowerInvariant();
        var cat = c.Catalogue;
        var g = new SphericalGrid(c.Grid);
        var (vp, vs, threeD) = c.Velocities();
        c.Log($"Relocation ({method}) in the {(threeD ? "3-D" : "1-D")} model.");
        List<RelocatedEvent> located;
        if (method is "absolute" or "abs")
        {
            var ls = c.Config.Locator;
            if (c.Args.Flag("fix-depth")) ls.FixDepth = true;
            ls.AutomaticPicks = c.Config.AutomaticPicks;
            located = Relocation.Absolute(cat, g, vp, vs, c.Config.Tomography.ForwardRefinement, ls, c.WorkFolder, c.Config.UseOpenCl, c.Log, c.Progress(), c.Cancel);
            CatalogueWriter.WriteLocations(Path.Combine(c.RunFolder, "locations.csv"), located);
            RunWriter.WriteSummary(c.RunFolder, new
            {
                Tool = "Relocation", Software = "TomoSTAR", Method = "absolute", Locator = ls, Grid = c.Grid, Model = threeD ? c.Args.Get("vp") : c.Model.Name,
                Located = located.Count(x => x.Ok), Failed = located.Count(x => !x.Ok),
                MedianRmsS = Median(located.Where(x => x.Ok).Select(x => x.Result.Rms)),
                MedianShiftKm = Median(located.Where(x => x.Ok).Select(Relocation.Shift))
            });
        }
        else if (method is "dd" or "double-difference")
        {
            var s = c.Config.Tomography.Clone();
            if (c.Args.Int("iterations") is { } it) s.DoubleDifference.Sets[^1].Iterations = Math.Max(1, it);
            var (list, run) = Relocation.DoubleDifference(cat, g, vp, vs, s, c.WorkFolder, c.Config.MinPhasesPerEvent, c.Config.OutsideData, c.Model, c.Log, c.Progress(), c.Cancel);
            located = list;
            CatalogueWriter.WriteLocations(Path.Combine(c.RunFolder, "locations.csv"), located);
            RunWriter.WriteIterations(Path.Combine(c.RunFolder, "iterations.csv"), run.Iterations);
            RunWriter.WriteResiduals(Path.Combine(c.RunFolder, "residuals.csv"), run.Data);
            RunWriter.WriteHypocentres(Path.Combine(c.RunFolder, "hypocentres.csv"), cat, run.Data);
            RunWriter.WriteSummary(c.RunFolder, new
            {
                Tool = "Relocation", Software = "TomoSTAR", Method = "double difference", s.DoubleDifference, Grid = c.Grid,
                Model = threeD ? c.Args.Get("vp") : c.Model.Name, run.Iterations, MedianShiftKm = Median(located.Select(Relocation.Shift))
            });
        }
        else throw new UsageException($"--method {method}: use absolute or dd.");
        var relocated = Relocation.Apply(cat, located);
        CatalogueWriter.WriteEvents(Path.Combine(c.RunFolder, "events_relocated.csv"), relocated.Events);
        // A complete data folder for the next step: the same stations and picks with the new events.
        CatalogueWriter.WriteEvents(Path.Combine(c.RunFolder, "events.csv"), relocated.Events);
        CatalogueWriter.WriteStations(Path.Combine(c.RunFolder, "stations.csv"), relocated.Stations);
        CatalogueWriter.WritePicks(Path.Combine(c.RunFolder, "picks.csv"), relocated.Events);
        if (cat.TStar.Count > 0) WriteTStarRecords(Path.Combine(c.RunFolder, "tstar.csv"), cat.TStar);
        c.Log(string.Create(Inv, $"Relocated {located.Count(x => x.Ok)} of {located.Count} events, median shift {Median(located.Where(x => x.Ok).Select(Relocation.Shift)):0.00} km."));
        c.Log($"locations.csv and a complete data folder (stations, events, picks) in {c.RunFolder}.");
        RunWriter.WriteLog(c.RunFolder, c.LogLines);
        return 0;
    }

    private static void WriteTStarRecords(string path, IEnumerable<TStarRecord> records) =>
        CatalogueWriter.WriteTStar(path, records.Select(t => new TStarMeasurement
        {
            EventId = t.EventId, StationId = t.StationId, Phase = t.Phase, TStar = t.TStar, Uncertainty = t.Sigma, Disabled = t.Disabled
        }));

    // ---- tstar ------------------------------------------------------------------------------------

    /// <summary>t* of every event from its waveforms and P (and S) picks, by the method of the TStar settings.</summary>
    public static int TStar(CommandContext c)
    {
        var cat = c.Catalogue;
        var archive = Archive(c);
        var s = c.Config.TStar;
        var stations = cat.StationIndex();
        var responses = c.Responses;
        if (responses.Count == 0) c.Log("No StationXML given (--stationxml): the responses are taken as flat over the fit band (velocity or acceleration channels).");
        InstrumentResponse? ResponseOf(Trace t) => StationXml.ResponseAt(responses, t.Nslc, t.StartTime);
        var all = new ConcurrentBag<TStarMeasurement>();
        var messages = new ConcurrentBag<string>();
        var done = 0;
        var progress = c.Progress();
        Parallel.ForEach(cat.Events, Core.Compute.ComputeSettings.Options(c.Cancel), ev =>
        {
            var traces = archive.ForEvent(ev);
            if (traces.Count > 0)
            {
                var (results, _, message) = TStarEstimator.MeasureEvent(ev, traces, s, responses.Count > 0 ? ResponseOf : null,
                    id => stations.GetValueOrDefault(id));
                foreach (var m in results) all.Add(m);
                messages.Add(message);
            }
            var d = Interlocked.Increment(ref done);
            progress.Report(((double)d / cat.Events.Count, $"t*: {d}/{cat.Events.Count} events"));
        });
        var list = all.OrderBy(m => m.EventId, StringComparer.Ordinal).ThenBy(m => m.StationId, StringComparer.Ordinal).ThenBy(m => m.Phase).ToList();
        CatalogueWriter.WriteTStar(Path.Combine(c.RunFolder, "tstar.csv"), list);
        File.WriteAllLines(Path.Combine(c.RunFolder, "tstar_events.txt"), messages.OrderBy(x => x, StringComparer.Ordinal));
        var rec = TStarQuality.Recommend(list, cat.EventIndex());
        RunWriter.WriteSummary(c.RunFolder, new
        {
            Tool = "t* measurement", Software = "TomoSTAR", Settings = s, Measurements = list.Count, Disabled = list.Count(m => m.Disabled),
            Events = list.Select(m => m.EventId).Distinct().Count(), MinimumMagnitude = rec
        });
        c.Log($"t*: {list.Count} measurements ({list.Count(m => m.Disabled)} flagged) from {list.Select(m => m.EventId).Distinct().Count()} events.");
        if (rec != null) c.Log("t* quality: " + rec.Explanation);
        c.Log($"tstar.csv in {c.RunFolder} (the input of 'tomostar qtomo --tstar').");
        RunWriter.WriteLog(c.RunFolder, c.LogLines);
        return 0;
    }

    private static double Median(IEnumerable<double> values)
    {
        var a = values.Where(double.IsFinite).OrderBy(x => x).ToArray();
        return a.Length == 0 ? double.NaN : a.Length % 2 == 1 ? a[a.Length / 2] : 0.5 * (a[a.Length / 2 - 1] + a[a.Length / 2]);
    }
}
