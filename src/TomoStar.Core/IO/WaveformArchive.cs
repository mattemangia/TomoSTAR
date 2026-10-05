// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using Microsoft.Data.Sqlite;
using TomoStar.Core.Model;

namespace TomoStar.Core.IO;

/// <summary>
/// Access to the waveforms of the events, for the automatic picker and the t* measurement. Three
/// layouts are understood, tried in this order for each event:
///
/// 1. A folder per event: <c>&lt;root&gt;/&lt;event id&gt;/*.mseed</c> (any file name; every file in the
///    folder is read). This is the layout <c>tomostar</c> documents and the quickest to read.
/// 2. The waveform references of a QUIVER project (table <c>waveforms</c> of its catalogue), paths
///    relative to the project folder; miniSEED only.
/// 3. Any other miniSEED files under the root (continuous data, or files named freely): they are
///    indexed once (channel and time span of every segment) and the segments overlapping the time
///    window of an event are read.
///
/// Every channel keeps its longest segment within the window. Traces are returned grouped by station
/// (NET.STA), which is what the picker and the t* estimators take.
/// </summary>
public sealed class WaveformArchive
{
    private readonly string? _root;
    private readonly string? _project;
    private List<(string Path, string StationId, string Channel, DateTime Start, DateTime End)>? _index;
    private Dictionary<string, List<(string Nslc, string Path)>>? _projectRefs;
    private readonly object _gate = new();

    /// <summary>An archive over a folder of miniSEED files.</summary>
    public static WaveformArchive FromFolder(string root)
    {
        if (!Directory.Exists(root)) throw new DirectoryNotFoundException($"Waveform folder not found: {root}");
        return new WaveformArchive(root, null);
    }

    /// <summary>An archive over the waveforms referenced by a QUIVER project.</summary>
    public static WaveformArchive FromQuiverProject(string projectFolder) => new(null, CatalogueReader.ProjectFolder(projectFolder));

    private WaveformArchive(string? root, string? project)
    {
        _root = root;
        _project = project;
    }

    /// <summary>Seconds of data read before the origin time (the noise window of the picker and of t*).</summary>
    public double BeforeSeconds { get; set; } = 30;

    /// <summary>Seconds of data read after the origin time (must hold the S wave at the farthest station).</summary>
    public double AfterSeconds { get; set; } = 120;

    /// <summary>The traces of one event, grouped by station.</summary>
    public Dictionary<string, List<Trace>> ForEvent(EventRecord ev, Action<string>? log = null)
    {
        var from = ev.OriginTime.AddSeconds(-BeforeSeconds);
        var to = ev.OriginTime.AddSeconds(AfterSeconds);
        IEnumerable<string> files;
        if (_root != null && Directory.Exists(Path.Combine(_root, Safe(ev.Id))))
            files = Directory.EnumerateFiles(Path.Combine(_root, Safe(ev.Id)), "*", SearchOption.AllDirectories);
        else if (_project != null)
            files = ProjectFiles(ev.Id);
        else
            files = Index(log).Where(x => x.End >= from && x.Start <= to).Select(x => x.Path).Distinct();

        // Longest segment of every channel inside the window.
        var best = new Dictionary<string, Trace>(StringComparer.OrdinalIgnoreCase);
        foreach (var file in files)
        {
            List<Trace> traces;
            try { traces = MiniSeed.ReadFile(file); }
            catch (Exception ex) when (ex is IOException or InvalidDataException or ArgumentException)
            {
                log?.Invoke($"Skipped {file}: {ex.Message}");
                continue;
            }
            foreach (var t in traces)
            {
                if (t.SampleRate <= 0 || t.Data.Length == 0 || t.EndTime < from || t.StartTime > to) continue;
                var cut = t.StartTime < from || t.EndTime > to ? t.Slice(from, to) : t;
                if (cut.Data.Length < 10) continue;
                if (!best.TryGetValue(cut.Nslc, out var old) || old.Data.Length < cut.Data.Length) best[cut.Nslc] = cut;
            }
        }
        var byStation = new Dictionary<string, List<Trace>>(StringComparer.OrdinalIgnoreCase);
        foreach (var t in best.Values)
        {
            if (!byStation.TryGetValue(t.StationId, out var list)) byStation[t.StationId] = list = [];
            list.Add(t);
        }
        // Several location codes or instruments at one station: keep the band with the highest sample
        // rate, velocity (H) before acceleration (N, G), so the picker reads one consistent set.
        foreach (var (id, list) in byStation.ToList())
        {
            var groups = list.GroupBy(t => t.Location + "." + (t.Channel.Length >= 2 ? t.Channel[..2] : t.Channel)).ToList();
            if (groups.Count <= 1) continue;
            var keep = groups.OrderByDescending(g => g.Any(t => t.IsVertical))
                .ThenBy(g => g.First().IsAccelerometer)
                .ThenByDescending(g => g.First().SampleRate).First();
            byStation[id] = keep.ToList();
        }
        return byStation;
    }

    private IEnumerable<string> ProjectFiles(string eventId)
    {
        lock (_gate)
        {
            if (_projectRefs == null)
            {
                _projectRefs = new Dictionary<string, List<(string, string)>>(StringComparer.Ordinal);
                var db = Path.Combine(_project!, "catalog.sqlite");
                if (File.Exists(db))
                {
                    var cs = new SqliteConnectionStringBuilder { DataSource = db, Mode = SqliteOpenMode.ReadOnly, Pooling = false }.ToString();
                    using var c = new SqliteConnection(cs);
                    c.Open();
                    using var cmd = c.CreateCommand();
                    cmd.CommandText = "SELECT event_id, network, station, location, channel, path, format FROM waveforms";
                    using var r = cmd.ExecuteReader();
                    while (r.Read())
                    {
                        if (!r.GetString(6).Contains("seed", StringComparison.OrdinalIgnoreCase)) continue;
                        var path = r.GetString(5);
                        if (!Path.IsPathRooted(path)) path = Path.Combine(_project!, path);
                        if (!_projectRefs.TryGetValue(r.GetString(0), out var list)) _projectRefs[r.GetString(0)] = list = [];
                        list.Add(($"{r.GetString(1)}.{r.GetString(2)}.{r.GetString(3)}.{r.GetString(4)}", path));
                    }
                }
            }
        }
        return _projectRefs.TryGetValue(eventId, out var refs) ? refs.Select(x => x.Item2).Where(File.Exists).Distinct() : [];
    }

    /// <summary>Channel and time span of every segment of every file under the root, read once.</summary>
    private List<(string Path, string StationId, string Channel, DateTime Start, DateTime End)> Index(Action<string>? log)
    {
        lock (_gate)
        {
            if (_index != null) return _index;
            _index = [];
            var files = Directory.EnumerateFiles(_root!, "*", SearchOption.AllDirectories).ToList();
            log?.Invoke($"Indexing {files.Count} waveform files under {_root}.");
            var bag = new System.Collections.Concurrent.ConcurrentBag<(string, string, string, DateTime, DateTime)>();
            Parallel.ForEach(files, Compute.ComputeSettings.Options(), f =>
            {
                try
                {
                    foreach (var t in MiniSeed.ReadFile(f))
                        if (t.SampleRate > 0 && t.Data.Length > 0) bag.Add((f, t.StationId, t.Channel, t.StartTime, t.EndTime));
                }
                catch (Exception ex) when (ex is IOException or InvalidDataException or ArgumentException)
                {
                    // Not miniSEED (a StationXML, a README...): ignored.
                }
            });
            _index = bag.ToList();
            return _index;
        }
    }

    /// <summary>An event id as a folder name (characters not allowed in file names replaced by '_').</summary>
    public static string Safe(string id) =>
        new(id.Select(c => char.IsLetterOrDigit(c) || c is '.' or '_' or '-' ? c : '_').ToArray());
}
