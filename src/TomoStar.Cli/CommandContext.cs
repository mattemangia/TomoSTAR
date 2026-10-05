// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using System.Globalization;
using TomoStar.Core.Compute;
using TomoStar.Core.Geo;
using TomoStar.Core.IO;
using TomoStar.Core.Model;
using TomoStar.Core.Tomography;

namespace TomoStar.Cli;

/// <summary>
/// What every command shares: its arguments and configuration, the log, the run folder, and the
/// loading of the data, the grid and the starting model, each done once and in the same way by
/// every command, so that one set of options means the same thing everywhere.
/// </summary>
public sealed class CommandContext
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    private readonly List<string> _log = [];
    private readonly object _logGate = new();
    private Catalogue? _catalogue;
    private GridDefinition? _grid;
    private VelocityModel1D? _model;
    private string? _runFolder;

    public CommandContext(string command, Arguments args, CancellationToken ct)
    {
        Command = command;
        Args = args;
        Cancel = ct;
        Quiet = args.Flag("quiet");
        var data = args.Get("data");
        if (data != null && CatalogueReader.IsQuiverProject(data)) Quiver = CatalogueReader.ReadQuiverProjectInfo(data);
        Config = TomoConfig.Build(args.Get("config"), args.All("set"), Quiver);
        if (args.Flag("no-opencl")) Config.UseOpenCl = false;
        if (args.Int("threads") is { } threads) Config.Threads = threads;
        if (args.Flag("adaptive")) { Config.Tomography.Adaptive.Enabled = true; Config.Attenuation.Adaptive.Enabled = true; }
        if (args.Flag("straight")) { Config.Tomography.RayMethod = Core.Forward.RayMethod.Straight; Config.Attenuation.RayMethod = Core.Forward.RayMethod.Straight; }
        if (args.Flag("csv") && !Config.OutputFormats.Contains("csv")) Config.OutputFormats.Add("csv");
        if (args.Flag("vtk") && !Config.OutputFormats.Contains("vtk")) Config.OutputFormats.Add("vtk");
        if (args.Get("formats") is { } formats) Config.OutputFormats = formats.Split(',', StringSplitOptions.RemoveEmptyEntries).ToList();
        Config.Tomography.UseOpenCl = Config.UseOpenCl;
        Config.Attenuation.UseOpenCl = Config.UseOpenCl;
        ComputeSettings.MaxThreads = Config.Threads;
        RunWriter.Formats = Config.Formats();
    }

    public string Command { get; }
    public Arguments Args { get; }
    public TomoConfig Config { get; }
    public QuiverProjectInfo? Quiver { get; }
    public CancellationToken Cancel { get; }
    public bool Quiet { get; }

    /// <summary>Every message of the command, written to log.txt in the run folder.</summary>
    public IReadOnlyList<string> LogLines
    {
        get { lock (_logGate) return _log.ToList(); }
    }

    public void Log(string message)
    {
        lock (_logGate)
        {
            _log.Add(message);
            if (!Quiet) Console.WriteLine(message);
        }
    }

    /// <summary>A progress reporter that prints at most every few seconds, or when the stage changes.</summary>
    public IProgress<(double, string)> Progress() => new ThrottledProgress(this);

    private sealed class ThrottledProgress(CommandContext c) : IProgress<(double, string)>
    {
        private readonly object _gate = new();
        private DateTime _last = DateTime.MinValue;
        private string _stage = "";

        public void Report((double, string) value)
        {
            if (c.Quiet) return;
            var (f, text) = value;
            var colon = text.IndexOf(':');
            var stage = colon > 0 ? text[..colon] : new string(text.TakeWhile(ch => !char.IsDigit(ch)).ToArray());
            lock (_gate)
            {
                var now = DateTime.UtcNow;
                if (stage == _stage && (now - _last).TotalSeconds < 3) return;
                _stage = stage;
                _last = now;
                Console.WriteLine(string.Create(Inv, $"  [{Math.Clamp(f, 0, 1) * 100,3:0}%] {text}"));
            }
        }
    }

    // ---- Run folder ------------------------------------------------------------------------------

    /// <summary>The folder of this run: --out, else runs/&lt;date_time&gt;_&lt;command&gt; in the current folder.</summary>
    public string RunFolder
    {
        get
        {
            if (_runFolder != null) return _runFolder;
            _runFolder = Path.GetFullPath(Args.Get("out") ?? Path.Combine("runs", $"{DateTime.Now:yyyyMMdd_HHmmss}_{Command}"));
            Directory.CreateDirectory(_runFolder);
            return _runFolder;
        }
    }

    /// <summary>The run name used in volume names: --name, else the folder name.</summary>
    public string RunName => Args.Get("name") ?? Path.GetFileName(RunFolder.TrimEnd(Path.DirectorySeparatorChar, Path.AltDirectorySeparatorChar));

    /// <summary>Scratch space for travel-time tables; removed at the end unless --keep-work.</summary>
    public string WorkFolder => Path.Combine(RunFolder, "work");

    public void CleanUp()
    {
        if (_runFolder == null || Args.Flag("keep-work")) return;
        try
        {
            if (Directory.Exists(WorkFolder)) Directory.Delete(WorkFolder, true);
        }
        catch (IOException) { }
        catch (UnauthorizedAccessException) { }
    }

    // ---- Data ----------------------------------------------------------------------------------

    /// <summary>StationXML responses by channel, when StationXML files were given.</summary>
    public Dictionary<string, List<StationXml.ChannelResponse>> Responses { get; } = new(StringComparer.OrdinalIgnoreCase);

    /// <summary>
    /// The data: --data (a folder with stations.csv, events.csv, picks.csv, tstar.csv, or a QUIVER
    /// project), the single tables (--stations, --events, --picks, --tstar, each replacing the one of
    /// --data), QuakeML and StationXML files (--quakeml, --stationxml, added), and --hypocentres (an
    /// event table whose positions and origin times replace those of the same events).
    /// </summary>
    public Catalogue Catalogue => _catalogue ??= LoadCatalogue();

    private Catalogue LoadCatalogue()
    {
        var data = Args.Get("data");
        Catalogue c;
        string? Table(string option, string file)
        {
            var explicitPath = Args.Get(option);
            if (explicitPath != null) return explicitPath;
            if (data != null && Directory.Exists(data) && Quiver == null)
            {
                var p = Path.Combine(data, file);
                return File.Exists(p) ? p : null;
            }
            return null;
        }
        if (Quiver != null)
        {
            c = CatalogueReader.ReadQuiverProject(data!, Args.Flag("catalog-locations"));
            Log($"Data: QUIVER project {Quiver.Folder} ({c.Stations.Count} stations, {c.Events.Count} events, {c.TStar.Count} t*).");
            if (Args.Get("stations") is { } s) CatalogueReader.ReadStations(s, c);
            if (Args.Get("events") is { } e) CatalogueReader.ReadEvents(e, c);
            if (Args.Get("picks") is { } p) CatalogueReader.ReadPicks(p, c);
            if (Args.Get("tstar") is { } t) CatalogueReader.ReadTStar(t, c);
        }
        else
        {
            if (data != null && !Directory.Exists(data)) throw new UsageException($"--data {data}: neither a folder nor a QUIVER project.");
            c = new Catalogue();
            if (Table("stations", "stations.csv") is { } s) CatalogueReader.ReadStations(s, c);
            if (Table("events", "events.csv") is { } e) CatalogueReader.ReadEvents(e, c);
            if (Table("picks", "picks.csv") is { } p) CatalogueReader.ReadPicks(p, c);
            if (Table("tstar", "tstar.csv") is { } t) CatalogueReader.ReadTStar(t, c);
        }
        // A --quakeml folder stands for every QuakeML file in it (one file per event, as some event
        // services return the picks of one event at a time).
        var quakeml = Args.All("quakeml").SelectMany<string, string>(q => Directory.Exists(q)
            ? Directory.EnumerateFiles(q).Where(f => f.EndsWith(".xml", StringComparison.OrdinalIgnoreCase) || f.EndsWith(".quakeml", StringComparison.OrdinalIgnoreCase)).Order(StringComparer.Ordinal)
            : [q]).ToList();
        var fromQuakeMl = (Events: 0, Picks: 0);
        // Parsed in parallel, added in file order (the same catalogue on any number of threads).
        var parsed = new List<EventRecord>[quakeml.Count];
        Parallel.For(0, quakeml.Count, ComputeSettings.Options(Cancel), i => parsed[i] = QuakeMl.ReadFile(quakeml[i]));
        var knownEvents = c.EventIndex();
        for (var i = 0; i < quakeml.Count; i++)
        {
            var added = 0;
            foreach (var ev in parsed[i].Where(ev => knownEvents.TryAdd(ev.Id, ev))) { c.Events.Add(ev); added++; }
            fromQuakeMl = (fromQuakeMl.Events + added, fromQuakeMl.Picks + parsed[i].Sum(ev => ev.Picks.Count));
            if (quakeml.Count <= 10) Log($"QuakeML {quakeml[i]}: {added} events, {parsed[i].Sum(ev => ev.Picks.Count)} picks.");
        }
        if (quakeml.Count > 10) Log($"QuakeML: {quakeml.Count} files, {fromQuakeMl.Events} events, {fromQuakeMl.Picks} picks.");
        foreach (var file in Args.All("stationxml"))
        {
            var r = StationXml.ReadFile(file);
            var known = c.StationIndex();
            var added = 0;
            foreach (var s in r.Stations.Where(s => !known.ContainsKey(s.Id))) { c.Stations.Add(s); added++; }
            foreach (var (k, v) in r.Responses)
            {
                if (!Responses.TryGetValue(k, out var list)) Responses[k] = list = [];
                list.AddRange(v);
            }
            Log($"StationXML {file}: {added} stations added, {r.Responses.Count} channel responses.");
        }
        if (Args.Get("hypocentres") is { } hyp)
        {
            var h = new Catalogue();
            CatalogueReader.ReadEvents(hyp, h);
            var byId = c.EventIndex();
            var moved = 0;
            foreach (var e in h.Events)
            {
                if (!byId.TryGetValue(e.Id, out var ev)) continue;
                (ev.OriginTime, ev.Lon, ev.Lat, ev.DepthKm) = (e.OriginTime, e.Lon, e.Lat, e.DepthKm);
                moved++;
            }
            Log($"Hypocentres from {hyp}: {moved} of {c.Events.Count} events updated.");
        }
        var relabelled = c.ResolvePickNetworks();
        if (relabelled > 0) Log($"{relabelled} picks took the network of the only station with their station code (bulletin and metadata label the station differently).");
        if (c.Stations.Count == 0) throw new UsageException("No stations: give --data (a folder or a QUIVER project), --stations or --stationxml.");
        var (np, ns) = c.PickCounts();
        if (Quiver == null) Log($"Data: {c.Stations.Count} stations, {c.Events.Count} events, {np} P and {ns} S picks, {c.TStar.Count} t*.");
        foreach (var w in c.Warnings) Log("Warning: " + w);
        return c;
    }

    // ---- Grid and models ------------------------------------------------------------------------

    /// <summary>
    /// The inversion grid: --grid (a grid.json or a QUIVER project), else the configuration's Grid
    /// (or the QUIVER project's), else the grid of --vp, else one proposed from the data.
    /// </summary>
    public GridDefinition Grid => _grid ??= ResolveGrid();

    private GridDefinition ResolveGrid()
    {
        GridDefinition g;
        if (Args.Get("grid") is { } path)
        {
            g = CatalogueReader.ReadGrid(path);
            Log($"Grid from {path}.");
        }
        else if (Config.Grid != null) g = Config.Grid;
        else if (Args.Get("vp") is { } vp)
        {
            g = VolumeFile.ReadHeader(vp).Grid;
            Log($"Grid of {vp}.");
        }
        else
        {
            var model = ModelFor(StationBox());
            var report = ParameterAdvisor.Advise(Catalogue, model, Config.GridAdvice);
            g = report.Grid;
            Log("No grid given: proposed from the data (run 'tomostar grid' to see the reasons and keep it).");
        }
        g.Validate();
        Log(string.Create(Inv, $"Grid: {g.Nx} x {g.Ny} x {g.Nz} nodes, {g.MinLon:0.###} to {g.MaxLon:0.###} E, {g.MinLat:0.###} to {g.MaxLat:0.###} N, {g.MinDepthKm:0.##} to {g.MaxDepthKm:0.##} km."));
        return g;
    }

    private GeoBox StationBox()
    {
        var c = Catalogue;
        var lons = c.Stations.Select(s => s.Lon).Concat(c.Events.Select(e => e.Lon)).ToArray();
        var lats = c.Stations.Select(s => s.Lat).Concat(c.Events.Select(e => e.Lat)).ToArray();
        return new GeoBox(lons.Min(), lons.Max(), lats.Min(), lats.Max());
    }

    /// <summary>The 1-D starting and background model (--model, else the configuration's StartingModel).</summary>
    public VelocityModel1D Model => _model ??= ModelFor(new GeoBox(Grid.MinLon, Grid.MaxLon, Grid.MinLat, Grid.MaxLat));

    private VelocityModel1D ModelFor(GeoBox area)
    {
        if (_model != null) return _model;
        var spec = Args.Get("model") ?? Config.StartingModel;
        VelocityModel1D m;
        if (spec.Equals("quiver", StringComparison.OrdinalIgnoreCase))
        {
            m = Quiver?.StartingModel ?? throw new UsageException("StartingModel = quiver needs a QUIVER project with a starting model as --data.");
            Log($"Starting model: '{m.Name}' (from the QUIVER project).");
        }
        else if (spec.Equals("auto", StringComparison.OrdinalIgnoreCase))
        {
            var (chosen, reason) = ReferenceModelLibrary.Select(area);
            m = chosen.Model();
            Log($"Starting model: {reason}");
        }
        else
        {
            m = CatalogueReader.ReadModel1D(spec);
            Log($"Starting model: '{m.Name}'.");
        }
        return _model = m;
    }

    /// <summary>
    /// Velocities on the grid: --vp and --vs volumes (a previous result; resampled when on another
    /// grid; Vs from Vp and the Vp/Vs of the 1-D model when only --vp is given), else the 1-D model.
    /// </summary>
    public (double[] Vp, double[] Vs, bool ThreeD) Velocities()
    {
        var g = new SphericalGrid(Grid);
        var (vp, vs) = TravelTimeTomography.StartingModel(g, Model);
        if (Args.Get("vp") is not { } vpFile) return (vp, vs, false);
        var vp3 = CatalogueReader.ReadVolumeOnGrid(vpFile, Grid);
        double[] vs3;
        if (Args.Get("vs") is { } vsFile) vs3 = CatalogueReader.ReadVolumeOnGrid(vsFile, Grid);
        else vs3 = vp3.Select((v, i) => v * vs[i] / vp[i]).ToArray();
        for (var i = 0; i < vp3.Length; i++)
        {
            if (!double.IsFinite(vp3[i])) vp3[i] = vp[i];
            if (!double.IsFinite(vs3[i])) vs3[i] = vs[i];
        }
        Log($"Velocity model from {vpFile}{(Args.Get("vs") is { } f ? " and " + f : " (Vs from the Vp/Vs of the 1-D model)")}.");
        return (vp3, vs3, true);
    }

    // ---- Registration ---------------------------------------------------------------------------

    /// <summary>Copies the run into a QUIVER project when --register (or --data as a project with --register=data) asks.</summary>
    public void RegisterIfAsked(IReadOnlyList<WrittenVolume> volumes, string group = "")
    {
        var target = Args.Get("register");
        if (target == null) return;
        if (target.Equals("data", StringComparison.OrdinalIgnoreCase)) target = Quiver?.Folder ?? throw new UsageException("--register data needs a QUIVER project as --data.");
        var n = RunWriter.RegisterInQuiverProject(target, RunFolder, RunName, volumes, group);
        Log($"Registered in the QUIVER project {target}: run folder runs/{RunName}, {n} volumes (open the project in QUIVER to see them).");
    }
}
