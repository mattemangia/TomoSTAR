using System.Globalization;
using System.Numerics;
using System.Text.Json;
using TomoStar.Core.Compute;
using TomoStar.Core.Forward;
using TomoStar.Core.Geo;
using TomoStar.Core.IO;
using TomoStar.Core.Model;
using TomoStar.Core.Numerics;
using TomoStar.Core.Tomography;

namespace TomoStar.Cli;

/// <summary>Setting up a run: synthetic data, grid and 1-D profile, the model library, configuration, devices, conversions.</summary>
public static class SetupCommands
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    // ---- synth -------------------------------------------------------------------------------------

    public static int Synth(CommandContext c)
    {
        var o = c.Config.Synthetic;
        if (c.Args.Flag("with-waveforms")) o.Waveforms = true;
        if (c.Args.Int("events") is { } ne) o.Events = ne;
        if (c.Args.Int("stations") is { } ns) o.Stations = ns;
        if (c.Args.Int("seed") is { } seed) o.Seed = seed;
        var r = SyntheticData.Generate(o, c.RunFolder, c.Log, c.Cancel);
        c.Log($"Synthetic data set: {r.Catalogue.Stations.Count} stations, {r.Catalogue.Events.Count} events, {r.Arrivals} picks, {r.TStar} t*{(o.Waveforms ? ", waveforms" : "")}.");
        c.Log($"Written to {c.RunFolder}: stations.csv, events.csv (perturbed starting locations), events_true.csv, picks.csv, tstar.csv, grid.json, model1d.txt, truth/*.qvol{(o.Waveforms ? ", waveforms/<event>/<event>.mseed" : "")}.");
        RunWriter.WriteLog(c.RunFolder, c.LogLines);
        return 0;
    }

    // ---- grid ----------------------------------------------------------------------------------------

    /// <summary>
    /// The inversion grid and the 1-D profile on it. From the data (the advisor: spacing from the
    /// stations, extent from stations and events, depth from the events and the turning depth of the
    /// rays), or from --bounds. Writes grid.json, advice.txt, profile1d.csv, model1d.txt, the starting
    /// volumes and a configuration with the grid and the proposed settings.
    /// </summary>
    public static int Grid(CommandContext c)
    {
        var adv = c.Config.GridAdvice;
        adv.HorizontalSpacingKm = c.Args.Double("h") ?? adv.HorizontalSpacingKm;
        adv.VerticalSpacingKm = c.Args.Double("dz") ?? adv.VerticalSpacingKm;
        adv.MarginKm = c.Args.Double("margin") ?? adv.MarginKm;
        adv.TopKm = c.Args.Double("top") ?? adv.TopKm;
        adv.BottomKm = c.Args.Double("bottom") ?? adv.BottomKm;
        GridDefinition def;
        TomographySettings settings;
        var advice = new List<string>();
        if (c.Args.Numbers("bounds") is { } b)
        {
            if (b.Length != 4) throw new UsageException("--bounds min_lon,max_lon,min_lat,max_lat");
            def = GridDefinition.FromSpacing(b[0], b[1], b[2], b[3], adv.TopKm ?? -2, adv.BottomKm ?? 30, adv.HorizontalSpacingKm ?? 5, adv.VerticalSpacingKm ?? 2);
            settings = c.Config.Tomography.Clone();
            advice.Add("Grid from --bounds and the given (or default) spacing and depth range.");
        }
        else
        {
            var model = c.Args.Get("model") is { } m ? CatalogueReader.ReadModel1D(m) : null;
            // The model used to estimate the turning depths: the given one, else the one for the station area.
            if (model == null)
            {
                var cat = c.Catalogue;
                var box = new GeoBox(cat.Stations.Min(s => s.Lon), cat.Stations.Max(s => s.Lon), cat.Stations.Min(s => s.Lat), cat.Stations.Max(s => s.Lat));
                model = c.Config.StartingModel is "auto" ? ReferenceModelLibrary.Select(box).Model.Model()
                    : c.Config.StartingModel is "quiver" ? c.Quiver?.StartingModel ?? ReferenceModelLibrary.Select(box).Model.Model()
                    : CatalogueReader.ReadModel1D(c.Config.StartingModel);
            }
            var report = ParameterAdvisor.Advise(c.Catalogue, model, adv);
            def = report.Grid;
            settings = report.Settings;
            advice.AddRange(report.Summary().Split(Environment.NewLine));
        }
        def.Validate();
        File.WriteAllText(Path.Combine(c.RunFolder, "grid.json"), JsonSerializer.Serialize(def, TomoJson.Options));
        // From here on the context resolves the grid we just made.
        var g = new SphericalGrid(def);
        var chosen = c.Args.Get("model") is { } spec ? CatalogueReader.ReadModel1D(spec)
            : c.Config.StartingModel is "auto" ? ReferenceModelLibrary.Select(new GeoBox(def.MinLon, def.MaxLon, def.MinLat, def.MaxLat)).Model.Model()
            : c.Config.StartingModel is "quiver" && c.Quiver?.StartingModel is { } qm ? qm
            : CatalogueReader.ReadModel1D(c.Config.StartingModel);
        advice.Add($"Starting model = {chosen.Name}{(chosen.Reference.Length > 0 ? ": " + chosen.Reference : "")}");
        File.WriteAllLines(Path.Combine(c.RunFolder, "advice.txt"), advice);
        WriteProfile(Path.Combine(c.RunFolder, "profile1d.csv"), chosen, g.DepthKm);
        File.WriteAllText(Path.Combine(c.RunFolder, "model1d.txt"), $"# {chosen.Name}\n# {chosen.Reference}\n# depth_km vp_km_s vs_km_s\n{chosen.Format()}\n");
        var (vp, vs) = TravelTimeTomography.StartingModel(g, chosen);
        var meta = new Dictionary<string, string> { ["model"] = chosen.Name };
        RunWriter.SaveVolume(c.RunFolder, "start", "Vp", "km/s", g, vp, meta);
        RunWriter.SaveVolume(c.RunFolder, "start", "Vs", "km/s", g, vs, meta);
        // A configuration to start from: the grid, the model and the proposed settings.
        var cfg = new TomoConfig { Grid = def, StartingModel = "model1d.txt", Tomography = settings };
        File.WriteAllText(Path.Combine(c.RunFolder, "config.json"), cfg.ToJson());
        foreach (var line in advice) c.Log(line);
        c.Log(string.Create(Inv, $"Grid {def.Nx} x {def.Ny} x {def.Nz} = {def.Count} nodes, spacing {g.Spacing(def.Nx / 2, def.Ny / 2, 0).HLon:0.##} x {g.Spacing(def.Nx / 2, def.Ny / 2, 0).HLat:0.##} x {def.DDepth:0.##} km."));
        c.Log($"grid.json, profile1d.csv, model1d.txt, volumes/start_Vp.qvol, volumes/start_Vs.qvol and config.json in {c.RunFolder}.");
        RunWriter.WriteLog(c.RunFolder, c.LogLines);
        return 0;
    }

    /// <summary>The 1-D model at the given depths: depth, Vp, Vs and Vp/Vs.</summary>
    private static void WriteProfile(string path, VelocityModel1D m, IEnumerable<double> depths) =>
        File.WriteAllLines(path, ["depth_km,vp_km_s,vs_km_s,vp_vs", .. depths.Select(z =>
            string.Create(Inv, $"{z:0.###},{m.Vp(z):0.####},{m.Vs(z):0.####},{m.Vp(z) / m.Vs(z):0.####}"))]);

    // ---- model1d ------------------------------------------------------------------------------------

    public static int Model1D(CommandContext c)
    {
        if (c.Args.Flag("list") || (c.Args.Positional.Count == 0 && c.Args.Get("area") == null))
        {
            Console.WriteLine($"{"id",-26} {"scope",-9} {"valid to",9}  name / DOI");
            foreach (var m in ReferenceModelLibrary.All)
                Console.WriteLine(string.Create(Inv, $"{m.Id,-26} {m.Scope,-9} {m.ValidToKm,6:0} km  {m.Name}  doi:{m.Doi}"));
            return 0;
        }
        if (c.Args.Numbers("area") is { } a)
        {
            if (a.Length != 4) throw new UsageException("--area min_lon,max_lon,min_lat,max_lat");
            var (best, reason) = ReferenceModelLibrary.Select(new GeoBox(a[0], a[1], a[2], a[3]));
            Console.WriteLine(reason);
            if (c.Args.Positional.Count == 0) c.Args.Positional.Add(best.Id);
        }
        var model = CatalogueReader.ReadModel1D(c.Args.Positional[0]);
        var text = $"# {model.Name}\n# {model.Reference}\n# depth_km vp_km_s vs_km_s\n{model.Format()}\n";
        if (c.Args.Get("out") is { } outFile)
        {
            File.WriteAllText(outFile, text);
            Console.WriteLine($"Written {outFile}.");
        }
        else Console.Write(text);
        if (c.Args.Get("grid") is { } gridFile)
        {
            var g = new SphericalGrid(CatalogueReader.ReadGrid(gridFile));
            var profile = c.Args.Get("profile") ?? "profile1d.csv";
            WriteProfile(profile, model, g.DepthKm);
            Console.WriteLine($"Profile at the {g.Nz} grid depths written to {profile}.");
        }
        return 0;
    }

    // ---- config -------------------------------------------------------------------------------------

    public static int Config(CommandContext c)
    {
        if (c.Args.Flag("template"))
        {
            var path = c.Args.Get("out") ?? "tomostar.json";
            File.WriteAllText(path, new TomoConfig().ToJson());
            Console.WriteLine($"Configuration template with every setting and its default written to {path}.");
            return 0;
        }
        // Default: show the configuration a command would run with (defaults, project, file, --set).
        Console.WriteLine(c.Config.ToJson());
        return 0;
    }

    // ---- info -----------------------------------------------------------------------------------------

    public static int Info(CommandContext c)
    {
        Console.WriteLine($"TomoSTAR {typeof(SetupCommands).Assembly.GetName().Version?.ToString(3)} on .NET {Environment.Version}, {Environment.OSVersion}");
        Console.WriteLine($"CPU: {Environment.ProcessorCount} logical processors; threads used: {ComputeSettings.Threads}");
        Console.WriteLine($"SIMD: Vector<double> holds {SimdVector.Width} values, hardware accelerated: {Vector.IsHardwareAccelerated}");
        var devices = OpenClContext.Enumerate(out var error);
        if (error != null) Console.WriteLine($"OpenCL: {error}");
        foreach (var d in devices) Console.WriteLine($"OpenCL device: {d}, {(d.IsGpu ? "GPU" : "CPU")} on {d.PlatformName}");
        if (devices.Count > 0)
        {
            Console.WriteLine("Self-tests (a device is used only when its results match the CPU solvers):");
            using (var e = EikonalOpenCl.TryCreate(m => Console.WriteLine("  " + m))) Console.WriteLine(e != null ? $"  eikonal: OK on {e.Device}" : "  eikonal: CPU");
            using (var l = LsqrOpenCl.TryCreate(m => Console.WriteLine("  " + m))) Console.WriteLine(l != null ? $"  LSQR: OK on {l.Device}" : "  LSQR: CPU");
        }
        return 0;
    }

    // ---- export -------------------------------------------------------------------------------------

    public static int Export(CommandContext c)
    {
        if (c.Args.Positional.Count == 0) throw new UsageException("export FILE.qvol [FILE.qvol ...] [--csv] [--vtk] [--out DIR]");
        var csv = c.Args.Flag("csv");
        var vtk = c.Args.Flag("vtk");
        if (!csv && !vtk) csv = true;
        var outDir = c.Args.Get("out");
        foreach (var file in c.Args.Positional)
        {
            using var r = new VolumeReader(file);
            var values = r.ReadAll();
            var target = Path.Combine(outDir ?? Path.GetDirectoryName(Path.GetFullPath(file))!, Path.GetFileNameWithoutExtension(file));
            if (outDir != null) Directory.CreateDirectory(outDir);
            if (csv) VolumeExport.WriteCsv(target + ".csv", r.Header, values);
            if (vtk) VolumeExport.WriteVtk(target + ".vtk", r.Header, values);
            Console.WriteLine($"{file}: {r.Header.Quantity} {r.Header.Grid.Nx}x{r.Header.Grid.Ny}x{r.Header.Grid.Nz} to {target}{(csv ? ".csv " : " ")}{(vtk ? target + ".vtk" : "")}");
        }
        return 0;
    }

    // ---- register -----------------------------------------------------------------------------------

    /// <summary>Adds an existing run folder (its volumes/*.qvol) to a QUIVER project.</summary>
    public static int Register(CommandContext c)
    {
        if (c.Args.Positional.Count == 0) throw new UsageException("register RUN_FOLDER --project QUIVER_PROJECT");
        var run = Path.GetFullPath(c.Args.Positional[0]);
        var project = c.Args.Require("project", "the QUIVER project folder");
        var volumes = Directory.Exists(Path.Combine(run, "volumes"))
            ? Directory.GetFiles(Path.Combine(run, "volumes"), "*.qvol").Select(f =>
            {
                var h = VolumeFile.ReadHeader(f);
                return new WrittenVolume(f, h.Name, h.Quantity, h.Units, h.Source);
            }).ToList()
            : [];
        var name = c.Args.Get("name") ?? Path.GetFileName(run.TrimEnd(Path.DirectorySeparatorChar));
        var n = RunWriter.RegisterInQuiverProject(project, run, name, volumes, c.Args.Get("group") ?? "");
        Console.WriteLine($"Registered {run} in {project}: {n} volumes. Close the project in QUIVER before registering, and reopen it to see them.");
        return 0;
    }
}
