using System.Text.Json;
using System.Text.Json.Nodes;
using Microsoft.Data.Sqlite;
using TomoStar.Core.Attenuation;
using TomoStar.Core.Geo;
using TomoStar.Core.Model;
using TomoStar.Core.Tomography;

namespace TomoStar.Core.IO;

/// <summary>
/// What a QUIVER project holds besides its catalogue: the grid, the starting 1-D model and the
/// inversion settings last chosen in the application. Every field is optional; the command line
/// fills in what the project does not give.
/// </summary>
public sealed class QuiverProjectInfo
{
    public required string Folder { get; init; }
    public string Name { get; init; } = "";
    public GridDefinition? Grid { get; init; }
    public VelocityModel1D? StartingModel { get; init; }
    public TomographySettings? Tomography { get; init; }
    public QTomographySettings? Attenuation { get; init; }
    public OutsideData? OutsideData { get; init; }
    public int? ForwardRefinement { get; init; }
    public bool? UseOpenCl { get; init; }
}

/// <summary>
/// Readers of the input data: plain CSV tables (see the README for the columns) or a QUIVER project
/// folder (<c>*.quiver</c>) read in place, without QUIVER itself: <c>project.json</c> for the grid, the
/// starting model and the settings, <c>stations.json</c> for the stations, and <c>catalog.sqlite</c>
/// (opened read-only) for events, picks and t* measurements.
/// </summary>
public static class CatalogueReader
{
    // ---- CSV ---------------------------------------------------------------------------------

    /// <summary>
    /// Reads stations, events and picks (and optionally t*) from CSV tables. Any of the paths may be
    /// null except the stations: a Q inversion needs no picks, a velocity inversion no t*.
    /// </summary>
    public static Catalogue ReadCsv(string stationsPath, string? eventsPath, string? picksPath, string? tstarPath = null)
    {
        var c = new Catalogue();
        ReadStations(stationsPath, c);
        if (eventsPath != null) ReadEvents(eventsPath, c);
        if (picksPath != null) ReadPicks(picksPath, c);
        if (tstarPath != null) ReadTStar(tstarPath, c);
        return c;
    }

    /// <summary>
    /// Station table. Columns: <c>id</c> (or <c>network</c> and <c>station</c>, joined as NET.STA),
    /// <c>lon</c>, <c>lat</c>, and <c>elevation_m</c> (or <c>depth_km</c>, positive down); optional
    /// <c>correction_p_s</c>, <c>correction_s_s</c>.
    /// </summary>
    public static void ReadStations(string path, Catalogue c)
    {
        var t = CsvTable.Read(path);
        var cId = t.Find("id", "station_id", "code");
        var cNet = t.Find("network", "net");
        var cSta = t.Find("station", "sta");
        if (cId < 0 && cSta < 0) throw new InvalidDataException($"{path}: no station id column (id, or network and station).");
        var cLon = t.Require("longitude", "lon", "longitude", "x");
        var cLat = t.Require("latitude", "lat", "latitude", "y");
        var cElev = t.Find("elevation_m", "elevation", "elev_m", "elev");
        var cDepth = t.Find("depth_km", "depth");
        var cCp = t.Find("correction_p_s", "correction_p");
        var cCs = t.Find("correction_s_s", "correction_s");
        foreach (var (line, cells) in t.Rows)
        {
            var id = cId >= 0 ? CsvTable.Cell(cells, cId)
                : cNet >= 0 && CsvTable.Cell(cells, cNet).Length > 0 ? $"{CsvTable.Cell(cells, cNet)}.{CsvTable.Cell(cells, cSta)}" : CsvTable.Cell(cells, cSta);
            if (id.Length == 0) throw new InvalidDataException($"{path}, line {line}: empty station id.");
            var elevation = cElev >= 0 ? t.Number(cells, cElev, line, 0)
                : cDepth >= 0 ? -1000 * t.Number(cells, cDepth, line, 0) : 0;
            c.Stations.Add(new StationRecord
            {
                Id = id, Lon = t.Number(cells, cLon, line), Lat = t.Number(cells, cLat, line), ElevationM = elevation,
                CorrectionP = t.Number(cells, cCp, line, 0), CorrectionS = t.Number(cells, cCs, line, 0)
            });
        }
    }

    /// <summary>
    /// Event table. Columns: <c>id</c>, <c>time</c> (origin time, UTC), <c>lon</c>, <c>lat</c>,
    /// <c>depth_km</c>; optional <c>magnitude</c> and <c>fixed</c> (1 to hold the hypocentre).
    /// </summary>
    public static void ReadEvents(string path, Catalogue c)
    {
        var t = CsvTable.Read(path);
        var cId = t.Require("event id", "id", "event", "event_id");
        var cTime = t.Require("origin time", "time", "origin_time", "origin_time_utc", "time_utc");
        var cLon = t.Require("longitude", "lon", "longitude");
        var cLat = t.Require("latitude", "lat", "latitude");
        var cDep = t.Require("depth", "depth_km", "depth");
        var cMag = t.Find("magnitude", "mag");
        var cFix = t.Find("fixed");
        var seen = new HashSet<string>();
        foreach (var (line, cells) in t.Rows)
        {
            var id = CsvTable.Cell(cells, cId);
            if (!seen.Add(id)) throw new InvalidDataException($"{path}, line {line}: event id '{id}' appears twice.");
            c.Events.Add(new EventRecord
            {
                Id = id, OriginTime = t.Time(cells, cTime, line), Lon = t.Number(cells, cLon, line), Lat = t.Number(cells, cLat, line),
                DepthKm = t.Number(cells, cDep, line), Magnitude = t.Number(cells, cMag, line), Fixed = CsvTable.Flag(cells, cFix)
            });
        }
    }

    /// <summary>
    /// Pick table. Columns: <c>event</c>, <c>station</c>, <c>phase</c> (P or S; Pg, Pn, Sg, Sn... are
    /// read by their first letter), and either <c>time</c> (absolute UTC) or <c>travel_time_s</c>
    /// (seconds after the origin time); optional <c>sigma_s</c> (default 0.1 s for P, 0.2 s for S),
    /// <c>quality</c> (0 to 4, 4 = unusable), <c>origin</c> (manual, catalog, automatic) and
    /// <c>disabled</c>.
    /// </summary>
    public static void ReadPicks(string path, Catalogue c)
    {
        var t = CsvTable.Read(path);
        var cEv = t.Require("event id", "event", "event_id", "id");
        var cSt = t.Require("station", "station", "station_id", "sta");
        var cPh = t.Require("phase", "phase");
        var cTime = t.Find("time", "time_utc", "arrival_time");
        var cTt = t.Find("travel_time_s", "travel_time", "tt_s", "tt");
        if (cTime < 0 && cTt < 0) throw new InvalidDataException($"{path}: no time column (time for absolute UTC, or travel_time_s).");
        var cSig = t.Find("sigma_s", "sigma", "uncertainty", "uncertainty_s");
        var cQ = t.Find("quality", "weight_class");
        var cOrigin = t.Find("origin", "source_type");
        var cDis = t.Find("disabled");
        var events = c.EventIndex();
        var unknown = 0;
        var otherPhase = 0;
        foreach (var (line, cells) in t.Rows)
        {
            if (!events.TryGetValue(CsvTable.Cell(cells, cEv), out var ev)) { unknown++; continue; }
            var phase = ParsePhase(CsvTable.Cell(cells, cPh));
            if (phase == null) { otherPhase++; continue; }
            var time = cTime >= 0 && CsvTable.Cell(cells, cTime).Length > 0
                ? t.Time(cells, cTime, line)
                : ev.OriginTime.AddTicks((long)Math.Round(t.Number(cells, cTt, line) * TimeSpan.TicksPerSecond));
            var origin = CsvTable.Cell(cells, cOrigin).ToLowerInvariant() switch
            {
                "manual" or "m" => PickOrigin.Manual,
                "automatic" or "auto" or "a" => PickOrigin.Automatic,
                _ => PickOrigin.Catalog
            };
            ev.Picks.Add(new PickRecord
            {
                StationId = CsvTable.Cell(cells, cSt), Phase = phase.Value, Time = time,
                Sigma = t.Number(cells, cSig, line, phase == Phase.P ? 0.1 : 0.2),
                Quality = (int)t.Number(cells, cQ, line, 0), Origin = origin, Disabled = CsvTable.Flag(cells, cDis)
            });
        }
        if (unknown > 0) c.Warnings.Add($"{path}: {unknown} picks refer to events not in the event table and were skipped.");
        if (otherPhase > 0) c.Warnings.Add($"{path}: {otherPhase} picks of phases other than P and S were skipped.");
    }

    /// <summary>
    /// t* table. Columns: <c>event</c>, <c>station</c>, <c>tstar_s</c>; optional <c>phase</c>
    /// (default P), <c>sigma_s</c> (default 0.005 s) and <c>disabled</c>.
    /// </summary>
    public static void ReadTStar(string path, Catalogue c)
    {
        var t = CsvTable.Read(path);
        var cEv = t.Require("event id", "event", "event_id", "id");
        var cSt = t.Require("station", "station", "station_id", "sta");
        var cT = t.Require("t*", "tstar_s", "tstar", "t_star");
        var cPh = t.Find("phase");
        var cSig = t.Find("sigma_s", "sigma", "uncertainty", "uncertainty_s");
        var cDis = t.Find("disabled");
        foreach (var (line, cells) in t.Rows)
        {
            var phase = cPh >= 0 ? ParsePhase(CsvTable.Cell(cells, cPh)) ?? Phase.P : Phase.P;
            c.TStar.Add(new TStarRecord
            {
                EventId = CsvTable.Cell(cells, cEv), StationId = CsvTable.Cell(cells, cSt), Phase = phase,
                TStar = t.Number(cells, cT, line), Sigma = t.Number(cells, cSig, line, 0.005), Disabled = CsvTable.Flag(cells, cDis)
            });
        }
    }

    /// <summary>P or S from a phase label (P, Pg, Pn, Pb, S, Sg, Sn...); null for anything else.</summary>
    public static Phase? ParsePhase(string label)
    {
        if (label.Length == 0) return null;
        return char.ToUpperInvariant(label[0]) switch
        {
            'P' => Phase.P,
            'S' => Phase.S,
            _ => null
        };
    }

    // ---- QUIVER project ------------------------------------------------------------------------

    /// <summary>Is this path a QUIVER project folder (or its project.json)?</summary>
    public static bool IsQuiverProject(string path) =>
        Directory.Exists(path) ? File.Exists(System.IO.Path.Combine(path, "project.json"))
            : File.Exists(path) && System.IO.Path.GetFileName(path).Equals("project.json", StringComparison.OrdinalIgnoreCase);

    /// <summary>The project folder of a path that may name the folder or its project.json.</summary>
    public static string ProjectFolder(string path) =>
        File.Exists(path) ? System.IO.Path.GetDirectoryName(System.IO.Path.GetFullPath(path))! : System.IO.Path.GetFullPath(path);

    /// <summary>Grid, starting model and settings stored in a QUIVER project.</summary>
    public static QuiverProjectInfo ReadQuiverProjectInfo(string path)
    {
        var folder = ProjectFolder(path);
        var file = System.IO.Path.Combine(folder, "project.json");
        if (!File.Exists(file)) throw new FileNotFoundException($"{folder} is not a QUIVER project (no project.json).", file);
        var root = JsonNode.Parse(File.ReadAllText(file), documentOptions: new JsonDocumentOptions { CommentHandling = JsonCommentHandling.Skip })?.AsObject()
                   ?? throw new InvalidDataException($"{file} is empty.");
        T? Get<T>(string name) where T : class
        {
            var node = Property(root, name);
            if (node == null) return null;
            try { return node.Deserialize<T>(TomoJson.Options); }
            catch (JsonException) { return null; } // a field of another schema version: leave it to the command line
        }
        OutsideData? outside = null;
        if (Property(root, "OutsideData")?.GetValue<string>() is { } o && Enum.TryParse<OutsideData>(o, true, out var od)) outside = od;
        return new QuiverProjectInfo
        {
            Folder = folder,
            Name = Property(root, "Name")?.GetValue<string>() ?? System.IO.Path.GetFileNameWithoutExtension(folder),
            Grid = Get<GridDefinition>("Grid"),
            StartingModel = Get<VelocityModel1D>("StartingModel"),
            Tomography = Get<TomographySettings>("Tomography"),
            Attenuation = Get<QTomographySettings>("Attenuation"),
            OutsideData = outside,
            ForwardRefinement = Property(root, "ForwardRefinement")?.GetValue<int>(),
            UseOpenCl = Property(root, "UseOpenCl")?.GetValue<bool>()
        };
    }

    private static JsonNode? Property(JsonObject o, string name)
    {
        foreach (var (k, v) in o)
            if (string.Equals(k, name, StringComparison.OrdinalIgnoreCase)) return v;
        return null;
    }

    // The parts of QUIVER's records that an inversion needs. Unknown fields are ignored, so the
    // reader follows additions to QUIVER's schema without change.
    private sealed class QStation
    {
        public string Network { get; set; } = "";
        public string Code { get; set; } = "";
        public double Latitude { get; set; }
        public double Longitude { get; set; }
        public double ElevationM { get; set; }
        public double CorrectionP { get; set; }
        public double CorrectionS { get; set; }
    }

    private sealed class QHypocentre
    {
        public DateTime OriginTime { get; set; }
        public double Latitude { get; set; }
        public double Longitude { get; set; }
        public double DepthKm { get; set; }
    }

    private sealed class QEvent
    {
        public string Id { get; set; } = "";
        public QHypocentre Catalog { get; set; } = new();
        public QHypocentre? Relocated { get; set; }
        public double Magnitude { get; set; } = double.NaN;
    }

    private sealed class QTStar
    {
        public string EventId { get; set; } = "";
        public string StationId { get; set; } = "";
        public Phase Phase { get; set; } = Phase.P;
        public double TStar { get; set; }
        public double Uncertainty { get; set; }
        public bool Disabled { get; set; }
    }

    /// <summary>
    /// Reads the catalogue of a QUIVER project. Each event starts from its latest QUIVER relocation
    /// when it has one (unless <paramref name="catalogLocations"/>), otherwise from the catalogue
    /// location. Picks keep their uncertainty, quality, origin (manual, catalogue, automatic) and
    /// disabled flag, so the pick selection is the one QUIVER itself makes.
    /// </summary>
    public static Catalogue ReadQuiverProject(string path, bool catalogLocations = false)
    {
        var folder = ProjectFolder(path);
        var c = new Catalogue();
        var stationsFile = System.IO.Path.Combine(folder, "stations.json");
        if (File.Exists(stationsFile))
        {
            var list = JsonSerializer.Deserialize<List<QStation>>(File.ReadAllText(stationsFile), TomoJson.Options) ?? [];
            foreach (var s in list)
                c.Stations.Add(new StationRecord
                {
                    Id = $"{s.Network}.{s.Code}", Lon = s.Longitude, Lat = s.Latitude, ElevationM = s.ElevationM,
                    CorrectionP = s.CorrectionP, CorrectionS = s.CorrectionS
                });
        }
        else c.Warnings.Add($"{folder}: no stations.json.");

        var db = System.IO.Path.Combine(folder, "catalog.sqlite");
        if (!File.Exists(db))
        {
            c.Warnings.Add($"{folder}: no catalog.sqlite (the project has no events yet).");
            return c;
        }
        var cs = new SqliteConnectionStringBuilder { DataSource = db, Mode = SqliteOpenMode.ReadOnly, Pooling = false }.ToString();
        using var conn = new SqliteConnection(cs);
        conn.Open();

        var byRid = new Dictionary<long, EventRecord>();
        using (var cmd = conn.CreateCommand())
        {
            cmd.CommandText = "SELECT rid, data FROM events ORDER BY ord, rid";
            using var r = cmd.ExecuteReader();
            while (r.Read())
            {
                var e = JsonSerializer.Deserialize<QEvent>(r.GetFieldValue<byte[]>(1), TomoJson.Options);
                if (e == null) continue;
                var h = !catalogLocations && e.Relocated != null ? e.Relocated : e.Catalog;
                var ev = new EventRecord
                {
                    Id = e.Id, OriginTime = DateTime.SpecifyKind(h.OriginTime, DateTimeKind.Utc),
                    Lon = h.Longitude, Lat = h.Latitude, DepthKm = h.DepthKm, Magnitude = e.Magnitude
                };
                byRid[r.GetInt64(0)] = ev;
                c.Events.Add(ev);
            }
        }
        var channels = new Dictionary<long, string>();
        using (var cmd = conn.CreateCommand())
        {
            cmd.CommandText = "SELECT cid, network, station FROM channels";
            using var r = cmd.ExecuteReader();
            while (r.Read()) channels[r.GetInt64(0)] = $"{r.GetString(1)}.{r.GetString(2)}";
        }
        using (var cmd = conn.CreateCommand())
        {
            // Pick flags (as QUIVER stores them): phase in bits 0-3, quality 4-7, origin 8-11
            // (0 catalogue, 1 automatic, 2 manual), polarity 12-14, disabled 15.
            cmd.CommandText = "SELECT event_rid, cid, time, uncertainty, flags FROM picks ORDER BY event_rid, seq";
            using var r = cmd.ExecuteReader();
            while (r.Read())
            {
                if (!byRid.TryGetValue(r.GetInt64(0), out var ev) || !channels.TryGetValue(r.GetInt64(1), out var station)) continue;
                var f = r.GetInt64(4);
                var phase = (int)(f & 15);
                if (phase > 1) continue;
                ev.Picks.Add(new PickRecord
                {
                    StationId = station, Phase = (Phase)phase, Time = new DateTime(r.GetInt64(2), DateTimeKind.Utc),
                    Sigma = r.GetDouble(3), Quality = (int)((f >> 4) & 15),
                    Origin = ((f >> 8) & 15) switch { 1 => PickOrigin.Automatic, 2 => PickOrigin.Manual, _ => PickOrigin.Catalog },
                    Disabled = ((f >> 15) & 1) != 0
                });
            }
        }
        using (var cmd = conn.CreateCommand())
        {
            cmd.CommandText = "SELECT data FROM tstar ORDER BY rid";
            using var r = cmd.ExecuteReader();
            while (r.Read())
            {
                var m = JsonSerializer.Deserialize<QTStar>(r.GetFieldValue<byte[]>(0), TomoJson.Options);
                if (m == null) continue;
                c.TStar.Add(new TStarRecord
                {
                    EventId = m.EventId, StationId = m.StationId, Phase = m.Phase, TStar = m.TStar,
                    Sigma = m.Uncertainty > 0 ? m.Uncertainty : 0.005, Disabled = m.Disabled
                });
            }
        }
        return c;
    }

    // ---- Models ----------------------------------------------------------------------------------

    /// <summary>
    /// A grid definition from JSON: either a bare <see cref="GridDefinition"/> (as written by
    /// <c>tomostar grid</c>) or a QUIVER project.json, whose "Grid" is taken.
    /// </summary>
    public static GridDefinition ReadGrid(string path)
    {
        if (IsQuiverProject(path))
            return ReadQuiverProjectInfo(path).Grid ?? throw new InvalidDataException($"{path}: the project has no grid yet.");
        var node = JsonNode.Parse(File.ReadAllText(path), documentOptions: new JsonDocumentOptions { CommentHandling = JsonCommentHandling.Skip, AllowTrailingCommas = true })!.AsObject();
        var gridNode = Property(node, "Grid") ?? node;
        var g = gridNode.Deserialize<GridDefinition>(TomoJson.Options) ?? throw new InvalidDataException($"{path}: no grid.");
        g.Validate();
        return g;
    }

    /// <summary>
    /// A 1-D model from a text file of "depth_km vp vs" lines (whitespace or commas, '#' comments;
    /// a repeated depth is a discontinuity), from a JSON <see cref="VelocityModel1D"/>, or by name
    /// from the library of published models (ak135, iasp91, PREM and the regional ones).
    /// </summary>
    public static VelocityModel1D ReadModel1D(string pathOrName)
    {
        if (!File.Exists(pathOrName))
        {
            if (ReferenceModelLibrary.TryGet(pathOrName, out var published)) return published.Model();
            throw new FileNotFoundException($"'{pathOrName}' is neither a file nor a known 1-D model (see 'tomostar model1d --list').");
        }
        var text = File.ReadAllText(pathOrName);
        if (text.TrimStart().StartsWith('{'))
        {
            var m = JsonSerializer.Deserialize<VelocityModel1D>(text, TomoJson.Options) ?? throw new InvalidDataException($"{pathOrName}: empty model.");
            m.Validate();
            return m;
        }
        return VelocityModel1D.Parse(System.IO.Path.GetFileNameWithoutExtension(pathOrName), text);
    }

    /// <summary>Reads a <c>.qvol</c> onto a grid; the volume must have been written on the same grid.</summary>
    public static double[] ReadVolumeOnGrid(string path, GridDefinition grid)
    {
        using var r = new VolumeReader(path);
        var g = r.Header.Grid;
        if (g.Nx != grid.Nx || g.Ny != grid.Ny || g.Nz != grid.Nz
            || Math.Abs(g.MinLon - grid.MinLon) > 1e-6 || Math.Abs(g.MaxLon - grid.MaxLon) > 1e-6
            || Math.Abs(g.MinLat - grid.MinLat) > 1e-6 || Math.Abs(g.MaxLat - grid.MaxLat) > 1e-6
            || Math.Abs(g.MinDepthKm - grid.MinDepthKm) > 1e-6 || Math.Abs(g.MaxDepthKm - grid.MaxDepthKm) > 1e-6)
        {
            // Another grid: resample by trilinear interpolation (clamped at the edges).
            var src = new SphericalGrid(g);
            var values = r.ReadAll();
            return new SphericalGrid(grid).ResampleFrom(src, values);
        }
        return r.ReadAll();
    }
}
