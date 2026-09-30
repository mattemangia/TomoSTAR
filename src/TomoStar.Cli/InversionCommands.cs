using System.Globalization;
using TomoStar.Core.Attenuation;
using TomoStar.Core.Geo;
using TomoStar.Core.IO;
using TomoStar.Core.Model;
using TomoStar.Core.Tomography;

namespace TomoStar.Cli;

/// <summary>The inversions: velocity tomography, minimum 1-D model, Q tomography, trade-off curves and resolution tests.</summary>
public static class InversionCommands
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>The arrivals of the catalogue on the grid, with the configuration's selection rules.</summary>
    private static ObservationSet Arrivals(CommandContext c, TomographySettings s, SphericalGrid g)
    {
        var data = ObservationSet.FromCatalogue(c.Catalogue, g, s.InvertP, s.InvertS, c.Config.MinPhasesPerEvent, c.Config.OutsideData, c.Config.AutomaticPicks);
        foreach (var w in data.Warnings) c.Log("Warning: " + w);
        if (data.Arrivals.Count == 0) throw new InvalidOperationException("No arrivals inside the grid: check the grid, the picks and MinPhasesPerEvent.");
        c.Log($"Inversion data: {data.Events.Count} events, {data.Stations.Count} stations, {data.Arrivals.Count(a => a.Phase == Phase.P)} P and {data.Arrivals.Count(a => a.Phase == Phase.S)} S arrivals.");
        return data;
    }

    /// <summary>The t* of the catalogue on the grid (phase from the Attenuation settings).</summary>
    private static ObservationSet TStarData(CommandContext c, SphericalGrid g, QTomographySettings s)
    {
        var data = QTomography.FromCatalogue(c.Catalogue, g, c.Config.OutsideData, s.Phase);
        foreach (var w in data.Warnings) c.Log("Warning: " + w);
        if (data.Arrivals.Count == 0) throw new InvalidOperationException($"No {s.Phase} t* inside the grid: give --tstar (or run 'tomostar tstar' first).");
        c.Log($"Q data: {data.Arrivals.Count} {s.Phase} t* from {data.Events.Count} events at {data.Stations.Count} stations.");
        return data;
    }

    // ---- invert / minimum1d ------------------------------------------------------------------------

    public static int Invert(CommandContext c) => Velocity(c, c.Args.Flag("layered"));

    public static int Minimum1D(CommandContext c) => Velocity(c, true);

    private static int Velocity(CommandContext c, bool layered)
    {
        var s = c.Config.Tomography.Clone();
        s.Layered = layered;
        if (c.Args.Int("iterations") is { } it) s.Iterations = it;
        if (c.Args.Double("damping") is { } d) s.DampingVelocity = d;
        if (c.Args.Double("smoothing") is { } sm) s.Smoothing = sm;
        var g = new SphericalGrid(c.Grid);
        var (vp, vs, _) = c.Velocities();
        var data = Arrivals(c, s, g);
        var background = Wadati.ApplyStart(s, vp, vs, c.Model, data, c.Log);
        var tomo = new TravelTimeTomography(g, s, c.WorkFolder, c.Log) { Background = background };
        var clock = System.Diagnostics.Stopwatch.StartNew();
        var r = tomo.Run(data, vp, vs, c.Progress(), c.Cancel);
        var volumes = RunWriter.SaveVelocity(c.RunFolder, c.RunName, r, s, c.Catalogue, c.Model, c.LogLines,
            new { ElapsedSeconds = clock.Elapsed.TotalSeconds });
        if (layered)
        {
            var m = TravelTimeTomography.LayeredModel(r, c.Model, $"minimum 1-D model ({c.RunName})");
            var file = Path.Combine(c.RunFolder, "model1d_minimum.txt");
            File.WriteAllText(file, $"# {m.Name}\n# {m.Reference}\n# depth_km vp_km_s vs_km_s\n{m.Format()}\n");
            c.Log($"Minimum 1-D model: {file} (use it with --model {file}).");
        }
        var last = r.Iterations[^1];
        c.Log(string.Create(Inv, $"Done in {clock.Elapsed.TotalSeconds:0} s: RMS {r.Iterations[0].Rms:0.000} to {last.Rms:0.000} s, variance reduction {last.VarianceReductionPercent:0.0} %."));
        c.Log($"Results in {c.RunFolder}.");
        RunWriter.WriteLog(c.RunFolder, c.LogLines);
        c.RegisterIfAsked(volumes);
        return 0;
    }

    // ---- qtomo ----------------------------------------------------------------------------------

    public static int QTomo(CommandContext c)
    {
        var s = c.Config.Attenuation.Clone();
        if (c.Args.Double("damping") is { } d) s.Damping = d;
        if (c.Args.Double("smoothing") is { } sm) s.Smoothing = sm;
        if (c.Args.Get("phase") is { } ph) s.Phase = Enum.Parse<Phase>(ph, true);
        var g = new SphericalGrid(c.Grid);
        var (vp, vs, threeD) = c.Velocities();
        if (!threeD) c.Log("Note: the rays are traced in the 1-D model; give the velocity model of the tomography with --vp and --vs.");
        var data = TStarData(c, g, s);
        var q = new QTomography(g, s, c.WorkFolder, c.Log) { Background = c.Model };
        var r = q.Run(data, vp, vs, c.Progress(), c.Cancel);
        var volumes = RunWriter.SaveQ(c.RunFolder, c.RunName, r, s, c.LogLines);
        c.Log(string.Create(Inv, $"Q tomography: reference Q {r.ReferenceQ:0}, t* RMS {r.RmsBefore * 1000:0.0} to {r.RmsAfter * 1000:0.0} ms. Results in {c.RunFolder}."));
        RunWriter.WriteLog(c.RunFolder, c.LogLines);
        c.RegisterIfAsked(volumes);
        return 0;
    }

    // ---- lcurve ---------------------------------------------------------------------------------

    public static int LCurve(CommandContext c)
    {
        var kind = (c.Args.Get("kind") ?? "velocity").ToLowerInvariant();
        var dampings = c.Args.Numbers("dampings") ?? c.Config.LCurve.Dampings;
        var smoothings = c.Args.Numbers("smoothings") ?? c.Config.LCurve.Smoothings;
        var saveModels = c.Args.Flag("save-models") || c.Config.LCurve.SaveModels;
        var pairs = smoothings.SelectMany(sm => dampings.Select(d => (d, sm))).ToList();
        if (pairs.Count < 3) throw new UsageException("An L-curve needs at least three damping values.");
        var g = new SphericalGrid(c.Grid);
        var (vp, vs, _) = c.Velocities();
        List<TravelTimeTomography.TradeOffSummary> summaries;
        List<(double, double)> variances;
        var volumes = new List<WrittenVolume>();
        object settings;
        string tool;
        if (kind is "velocity" or "v" or "vp")
        {
            var s = c.Config.Tomography.Clone();
            var data = Arrivals(c, s, g);
            var background = Wadati.ApplyStart(s, vp, vs, c.Model, data, c.Log);
            c.Log($"L-curve (velocity): first linearised step for {pairs.Count} pairs of damping and smoothing.");
            var points = new TravelTimeTomography(g, s, c.WorkFolder, c.Log) { Background = background }.TradeOff(data, vp, vs, pairs, c.Progress(), c.Cancel);
            summaries = points.Select(p => new TravelTimeTomography.TradeOffSummary(p.Damping, p.Smoothing, p.ResidualNorm, p.ModelNorm, p.Roughness)).ToList();
            variances = points.Select(p => (p.DataVariance, p.ModelVariance)).ToList();
            if (saveModels)
                foreach (var p in points)
                {
                    var tag = string.Create(Inv, $"{c.RunName}_d{p.Damping:G4}_s{p.Smoothing:G4}");
                    volumes.Add(RunWriter.SaveVolume(c.RunFolder, tag, "Vp", "km/s", g, p.Vp));
                    volumes.Add(RunWriter.SaveVolume(c.RunFolder, tag, "dVp", "%", g, p.Vp.Select((v, i) => 100 * (v / vp[i] - 1)).ToArray()));
                    if (s.InvertS) volumes.Add(RunWriter.SaveVolume(c.RunFolder, tag, "Vs", "km/s", g, p.Vs));
                }
            settings = s;
            tool = "L-curve (velocity)";
        }
        else if (kind is "q" or "attenuation")
        {
            var s = c.Config.Attenuation.Clone();
            var data = TStarData(c, g, s);
            c.Log($"L-curve (Q): complete linear inversion for {pairs.Count} pairs of damping and smoothing.");
            var points = new QTomography(g, s, c.WorkFolder, c.Log) { Background = c.Model }.TradeOff(data, vp, vs, pairs, c.Progress(), c.Cancel);
            summaries = points.Select(p => new TravelTimeTomography.TradeOffSummary(p.Damping, p.Smoothing, p.ResidualNorm, p.ModelNorm, p.Roughness)).ToList();
            variances = points.Select(p => (p.DataVariance, p.ModelVariance)).ToList();
            if (saveModels)
                foreach (var p in points)
                    volumes.Add(RunWriter.SaveVolume(c.RunFolder, string.Create(Inv, $"{c.RunName}_d{p.Damping:G4}_s{p.Smoothing:G4}"), s.Phase == Phase.S ? "Qs" : "Qp", "", g, p.Q));
            settings = s;
            tool = "L-curve (Q)";
        }
        else throw new UsageException($"--kind {kind}: use velocity or q.");

        var chosen = TravelTimeTomography.ChooseRegularisation(summaries);
        var corners = summaries.GroupBy(p => p.Smoothing).OrderBy(x => x.Key).Select(gr =>
        {
            var list = gr.OrderBy(p => p.Damping).ToList();
            var k = TravelTimeTomography.Corner(list.Select(p => (p.ResidualNorm, p.ModelNorm)).ToList());
            return (gr.Key, list[Math.Clamp(k, 0, list.Count - 1)].Damping);
        }).ToList();
        foreach (var (sm, d) in corners) c.Log(string.Create(Inv, $"  smoothing {sm:G6}: L-curve corner at damping {d:G6}"));
        c.Log(string.Create(Inv, $"Recommended: damping {chosen.Damping:G6}, smoothing {chosen.Smoothing:G6}."));
        RunWriter.SaveTradeOff(c.RunFolder, tool, summaries, variances, chosen, corners, settings, c.LogLines);
        Plots.LCurveSvg(Path.Combine(c.RunFolder, "lcurve.svg"), tool, summaries, chosen);
        if (c.Args.Get("update-config") is { } cfg)
        {
            var (dPath, sPath) = kind.StartsWith('q') || kind == "attenuation" ? ("Attenuation.Damping", "Attenuation.Smoothing") : ("Tomography.DampingVelocity", "Tomography.Smoothing");
            TomoConfig.UpdateFile(cfg, [(dPath, chosen.Damping), (sPath, chosen.Smoothing)]);
            c.Log($"Wrote the recommended damping and smoothing into {cfg}.");
        }
        c.Log($"Table lcurve.csv, figure lcurve.svg and summary.json in {c.RunFolder}.");
        RunWriter.WriteLog(c.RunFolder, c.LogLines);
        if (volumes.Count > 0) c.RegisterIfAsked(volumes, "L-curve");
        return 0;
    }

    // ---- resolution -------------------------------------------------------------------------------

    public static int Resolution(CommandContext c)
    {
        var rc = c.Config.Resolution;
        var patternName = (c.Args.Get("pattern") ?? rc.Pattern).ToLowerInvariant();
        var target = (c.Args.Get("target") ?? rc.Target).ToLowerInvariant();
        var cell = c.Args.Double("cell") ?? rc.CellKm;
        var cellZ = c.Args.Double("cell-depth") ?? rc.CellDepthKm;
        var amp = c.Args.Double("amplitude") ?? rc.AmplitudePercent;
        var spikes = rc.Spikes.ToList();
        foreach (var s in c.Args.All("spike"))
        {
            var v = Arguments.ParseNumbers(s, "--spike");
            if (v.Length != 5) throw new UsageException("--spike lon,lat,depth_km,radius_km,amplitude_percent");
            spikes.Add(new Spike(v[0], v[1], v[2], v[3], v[4]));
        }
        var bodies = rc.Bodies.ToList();
        foreach (var s in c.Args.All("body"))
        {
            var v = Arguments.ParseNumbers(s, "--body");
            if (v.Length != 9) throw new UsageException("--body lon,lat,depth_km,strike_deg,dip_deg,length_km,width_km,thickness_km,amplitude_percent");
            bodies.Add(new TabularBody(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8]));
        }
        var pattern = patternName switch
        {
            "checkerboard" or "checker" => new ResolutionPattern { IsCheckerboard = true, CellHorizontalKm = cell, CellVerticalKm = cellZ, AmplitudePercent = amp },
            "spike" or "spikes" or "body" or "bodies" => new ResolutionPattern { IsCheckerboard = false, Spikes = spikes, Bodies = bodies },
            _ => throw new UsageException($"--pattern {patternName}: use checkerboard, spike or body.")
        };
        if (!pattern.IsCheckerboard && spikes.Count == 0 && bodies.Count == 0)
            throw new UsageException("A spike or body test needs --spike or --body (or Resolution.Spikes / Resolution.Bodies in the configuration).");
        var g = new SphericalGrid(c.Grid);
        var (vp, vs, _) = c.Velocities();
        var seed = c.Args.Int("seed") ?? rc.Seed;
        ResolutionTestResult result;
        if (target is "q" or "qleakage" or "q-leakage")
        {
            var s = c.Config.Attenuation.Clone();
            var geometry = TStarData(c, g, s);
            result = ResolutionTests.Attenuation(g, s, geometry, vp, vs, pattern, rc.NoiseTStar, seed, c.WorkFolder, c.Progress(), c.Log, c.Cancel,
                c.Model, velocityLeakage: target != "q");
        }
        else
        {
            var t = target switch
            {
                "vp" => ResolutionTarget.Vp,
                "vs" => ResolutionTarget.Vs,
                "vpvs" => ResolutionTarget.VpVs,
                "leakage" or "vpvsleakage" => ResolutionTarget.VpVsLeakage,
                _ => throw new UsageException($"--target {target}: use Vp, Vs, VpVs, Leakage, Q or QLeakage.")
            };
            var s = c.Config.Tomography.Clone();
            var geometry = Arrivals(c, s, g);
            result = ResolutionTests.Velocity(g, s, geometry, vp, vs, pattern, t, rc.NoiseP, rc.NoiseS, seed, c.WorkFolder, c.Progress(), c.Log, c.Cancel, c.Model);
        }
        c.Log(string.Create(Inv, $"Resolution test ({result.Description}): correlation {result.Correlation:0.000} over the sampled nodes, {result.CorrelationWellSampled:0.000} over the well-sampled ones, error RMS {100 * result.ErrorRms:0.00} %."));
        var volumes = RunWriter.SaveResolutionTest(c.RunFolder, c.RunName, result, g, c.LogLines);
        c.Log($"Results in {c.RunFolder}.");
        RunWriter.WriteLog(c.RunFolder, c.LogLines);
        c.RegisterIfAsked(volumes);
        return 0;
    }
}
