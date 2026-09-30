// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using System.Globalization;
using System.Text.Json;
using System.Text.Json.Nodes;
using TomoStar.Core.Attenuation;
using TomoStar.Core.Geo;
using TomoStar.Core.IO;
using TomoStar.Core.Location;
using TomoStar.Core.Model;
using TomoStar.Core.Signal;
using TomoStar.Core.Tomography;

namespace TomoStar.Cli;

/// <summary>The damping and smoothing values a trade-off run tries.</summary>
public sealed class LCurveConfig
{
    public double[] Dampings { get; set; } = [1, 2, 5, 10, 20, 50, 100, 200, 500];
    public double[] Smoothings { get; set; } = [5, 20, 80];

    /// <summary>Keep the model of every point as volumes (grouped as "L-curve" in QUIVER).</summary>
    public bool SaveModels { get; set; }
}

/// <summary>The synthetic pattern and noise of a resolution test.</summary>
public sealed class ResolutionConfig
{
    /// <summary>checkerboard, spike or body.</summary>
    public string Pattern { get; set; } = "checkerboard";

    /// <summary>Vp, Vs, VpVs, Leakage (same pattern in Vp and Vs), Q or QLeakage (pattern in velocity, measured in 1/Q).</summary>
    public string Target { get; set; } = "Vp";

    public double CellKm { get; set; } = 10;
    public double CellDepthKm { get; set; } = 6;
    public double AmplitudePercent { get; set; } = 5;
    public List<Spike> Spikes { get; set; } = [];
    public List<TabularBody> Bodies { get; set; } = [];

    /// <summary>Gaussian noise added to the synthetic P and S times and to the t*, s.</summary>
    public double NoiseP { get; set; } = 0.03;

    public double NoiseS { get; set; } = 0.06;
    public double NoiseTStar { get; set; } = 0.002;
    public int Seed { get; set; } = 1;
}

/// <summary>How much waveform is read around each event.</summary>
public sealed class WaveformConfig
{
    public double BeforeSeconds { get; set; } = 30;
    public double AfterSeconds { get; set; } = 120;
}

/// <summary>Hypocentre relocation.</summary>
public sealed class RelocationConfig
{
    /// <summary>absolute (event by event) or dd (double difference, all events at once).</summary>
    public string Method { get; set; } = "absolute";
}

/// <summary>
/// Every setting of TomoSTAR in one JSON document. A command starts from the defaults, then applies
/// in turn the settings of a QUIVER project given as data (grid, starting model, tomography and
/// attenuation settings), the configuration file (<c>--config</c>) and the <c>--set Path=value</c>
/// overrides. <c>tomostar config --template</c> writes the whole document with its defaults.
/// </summary>
public sealed class TomoConfig
{
    /// <summary>The inversion grid; null to propose one from the data (see GridAdvice).</summary>
    public GridDefinition? Grid { get; set; }

    /// <summary>What the automatic grid proposal may not choose freely.</summary>
    public GridAdviceOptions GridAdvice { get; set; } = new();

    /// <summary>
    /// The 1-D starting (and background) model: a name from the library (ak135, iasp91, PREM or a
    /// regional model, see 'tomostar model1d --list'), a file of depth vp vs lines or a JSON model,
    /// "auto" for the published model that suits the grid area best, or "quiver" for the model of
    /// the QUIVER project given as data.
    /// </summary>
    public string StartingModel { get; set; } = "auto";

    public OutsideData OutsideData { get; set; } = OutsideData.Exclude;
    public int MinPhasesPerEvent { get; set; } = 4;
    public AutomaticPicks AutomaticPicks { get; set; } = AutomaticPicks.IfFewAnalystPicks;

    /// <summary>Use OpenCL devices (eikonal and LSQR) when one passes its self-test.</summary>
    public bool UseOpenCl { get; set; } = true;

    /// <summary>CPU threads (0 = all).</summary>
    public int Threads { get; set; }

    /// <summary>Formats written for every volume besides .qvol: csv, vtk.</summary>
    public List<string> OutputFormats { get; set; } = ["csv"];

    public TomographySettings Tomography { get; set; } = new();
    public QTomographySettings Attenuation { get; set; } = new();
    public LocatorSettings Locator { get; set; } = new();
    public RelocationConfig Relocation { get; set; } = new();
    public AutoPickerSettings Picker { get; set; } = new();
    public TStarSettings TStar { get; set; } = new();
    public LCurveConfig LCurve { get; set; } = new();
    public ResolutionConfig Resolution { get; set; } = new();
    public WaveformConfig Waveforms { get; set; } = new();
    public SyntheticOptions Synthetic { get; set; } = new();

    /// <summary>The configuration a command runs with (see the class summary for the order).</summary>
    public static TomoConfig Build(string? configFile, IEnumerable<string> sets, QuiverProjectInfo? quiver)
    {
        var node = JsonSerializer.SerializeToNode(new TomoConfig(), TomoJson.Options)!.AsObject();
        if (quiver != null)
        {
            var q = new JsonObject();
            if (quiver.Grid != null) q["Grid"] = JsonSerializer.SerializeToNode(quiver.Grid, TomoJson.Options);
            if (quiver.StartingModel != null) q["StartingModel"] = "quiver";
            if (quiver.Tomography != null) q["Tomography"] = JsonSerializer.SerializeToNode(quiver.Tomography, TomoJson.Options);
            if (quiver.Attenuation != null) q["Attenuation"] = JsonSerializer.SerializeToNode(quiver.Attenuation, TomoJson.Options);
            if (quiver.OutsideData != null) q["OutsideData"] = quiver.OutsideData.ToString();
            if (quiver.UseOpenCl != null) q["UseOpenCl"] = quiver.UseOpenCl.Value;
            Merge(node, q);
        }
        if (configFile != null)
        {
            if (!File.Exists(configFile)) throw new UsageException($"Configuration file not found: {configFile}");
            var file = JsonNode.Parse(File.ReadAllText(configFile), documentOptions: new JsonDocumentOptions { CommentHandling = JsonCommentHandling.Skip, AllowTrailingCommas = true });
            if (file is not JsonObject fo) throw new UsageException($"{configFile}: the configuration must be a JSON object.");
            // A relative model file is relative to the configuration file.
            if (Property(fo, "StartingModel") is JsonValue sm && sm.TryGetValue<string>(out var name) && !Path.IsPathRooted(name))
            {
                var candidate = Path.Combine(Path.GetDirectoryName(Path.GetFullPath(configFile))!, name);
                if (File.Exists(candidate)) SetProperty(fo, "StartingModel", candidate);
            }
            Merge(node, fo);
        }
        foreach (var s in sets) Apply(node, s);
        try
        {
            return node.Deserialize<TomoConfig>(TomoJson.Options) ?? new TomoConfig();
        }
        catch (JsonException ex)
        {
            throw new UsageException($"Invalid configuration: {ex.Message}");
        }
    }

    /// <summary>
    /// Applies one "Path.To.Property=value" override. The value is read as JSON when it is valid
    /// JSON (numbers, true/false, [arrays], {objects}, "strings"), otherwise as a plain string.
    /// </summary>
    public static void Apply(JsonObject root, string assignment)
    {
        var eq = assignment.IndexOf('=');
        if (eq <= 0) throw new UsageException($"--set needs Path=value, got '{assignment}'.");
        var path = assignment[..eq].Trim().Split('.', StringSplitOptions.RemoveEmptyEntries);
        var text = assignment[(eq + 1)..].Trim();
        JsonNode? value;
        try { value = JsonNode.Parse(text); }
        catch (JsonException) { value = JsonValue.Create(text); }
        var o = root;
        for (var i = 0; i < path.Length - 1; i++)
        {
            var child = Property(o, path[i]);
            if (child is not JsonObject co)
            {
                co = new JsonObject();
                SetProperty(o, path[i], co);
            }
            o = co;
        }
        SetProperty(o, path[^1], value);
    }

    /// <summary>Deep merge: objects are merged key by key, anything else replaces.</summary>
    public static void Merge(JsonObject target, JsonObject source)
    {
        foreach (var (key, value) in source.ToList())
        {
            if (value is JsonObject so && Property(target, key) is JsonObject to) Merge(to, so);
            else SetProperty(target, key, value?.DeepClone());
        }
    }

    private static JsonNode? Property(JsonObject o, string name)
    {
        foreach (var (k, v) in o)
            if (string.Equals(k, name, StringComparison.OrdinalIgnoreCase)) return v;
        return null;
    }

    private static void SetProperty(JsonObject o, string name, JsonNode? value)
    {
        var key = o.Select(kv => kv.Key).FirstOrDefault(k => string.Equals(k, name, StringComparison.OrdinalIgnoreCase)) ?? name;
        o[key] = value;
    }

    /// <summary>The volume formats asked for.</summary>
    public VolumeFormats Formats()
    {
        var f = VolumeFormats.Qvol;
        foreach (var s in OutputFormats)
            f |= s.Trim().ToLowerInvariant() switch
            {
                "csv" => VolumeFormats.Csv,
                "vtk" => VolumeFormats.Vtk,
                "qvol" or "" => VolumeFormats.Qvol,
                _ => throw new UsageException($"Unknown output format '{s}' (csv, vtk).")
            };
        return f;
    }

    public string ToJson() => JsonSerializer.Serialize(this, TomoJson.Options);

    /// <summary>Sets one value in a configuration file on disk, creating the file if needed (used by 'lcurve --update-config').</summary>
    public static void UpdateFile(string path, IEnumerable<(string Path, double Value)> values)
    {
        JsonObject root = File.Exists(path)
            ? JsonNode.Parse(File.ReadAllText(path), documentOptions: new JsonDocumentOptions { CommentHandling = JsonCommentHandling.Skip, AllowTrailingCommas = true }) as JsonObject ?? new JsonObject()
            : new JsonObject();
        foreach (var (p, v) in values) Apply(root, $"{p}={v.ToString("R", CultureInfo.InvariantCulture)}");
        File.WriteAllText(path, root.ToJsonString(new JsonSerializerOptions { WriteIndented = true }));
    }
}
