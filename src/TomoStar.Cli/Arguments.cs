using System.Globalization;

namespace TomoStar.Cli;

/// <summary>
/// The arguments of one command: positional values, options (<c>--name value</c> or
/// <c>--name=value</c>) and flags (<c>--name</c> alone). Which names are flags is fixed in
/// <see cref="Flags"/>, so an option value may itself start with a minus sign (a longitude).
/// Options listed in <see cref="Repeatable"/> may be given several times.
/// </summary>
public sealed class Arguments
{
    /// <summary>Options that take no value.</summary>
    public static readonly HashSet<string> Flags = new(StringComparer.OrdinalIgnoreCase)
    {
        "help", "quiet", "no-opencl", "list", "layered", "only-automatic", "catalog-locations", "save-models", "keep-going",
        "with-waveforms", "csv", "vtk", "keep-work", "adaptive", "show", "template", "dry-run", "fix-depth", "straight"
    };

    /// <summary>Options that may be repeated; their values are kept in order.</summary>
    public static readonly HashSet<string> Repeatable = new(StringComparer.OrdinalIgnoreCase) { "set", "var", "spike", "body", "quakeml", "stationxml" };

    private readonly Dictionary<string, List<string>> _options = new(StringComparer.OrdinalIgnoreCase);
    private readonly HashSet<string> _flags = new(StringComparer.OrdinalIgnoreCase);
    private readonly HashSet<string> _used = new(StringComparer.OrdinalIgnoreCase);

    public List<string> Positional { get; } = [];

    public static Arguments Parse(IReadOnlyList<string> tokens)
    {
        var a = new Arguments();
        for (var i = 0; i < tokens.Count; i++)
        {
            var t = tokens[i];
            if (t.StartsWith("--", StringComparison.Ordinal) && t.Length > 2)
            {
                var name = t[2..];
                string? value = null;
                var eq = name.IndexOf('=');
                if (eq >= 0) { value = name[(eq + 1)..]; name = name[..eq]; }
                if (Flags.Contains(name) && value == null) { a._flags.Add(name); continue; }
                if (value == null)
                {
                    if (i + 1 >= tokens.Count) throw new UsageException($"Option --{name} needs a value.");
                    value = tokens[++i];
                }
                if (!a._options.TryGetValue(name, out var list)) a._options[name] = list = [];
                if (list.Count > 0 && !Repeatable.Contains(name)) list.Clear(); // the last one wins
                list.Add(value);
            }
            else a.Positional.Add(t);
        }
        return a;
    }

    public bool Flag(string name)
    {
        _used.Add(name);
        return _flags.Contains(name);
    }

    public bool Has(string name) => _options.ContainsKey(name) || _flags.Contains(name);

    public string? Get(string name)
    {
        _used.Add(name);
        return _options.TryGetValue(name, out var v) ? v[^1] : null;
    }

    public IReadOnlyList<string> All(string name)
    {
        _used.Add(name);
        return _options.TryGetValue(name, out var v) ? v : [];
    }

    public string Require(string name, string what) => Get(name) ?? throw new UsageException($"Missing --{name} ({what}).");

    public double? Double(string name)
    {
        var s = Get(name);
        if (s == null) return null;
        return double.TryParse(s, NumberStyles.Float, CultureInfo.InvariantCulture, out var v) ? v : throw new UsageException($"--{name}: '{s}' is not a number.");
    }

    public int? Int(string name)
    {
        var s = Get(name);
        if (s == null) return null;
        return int.TryParse(s, NumberStyles.Integer, CultureInfo.InvariantCulture, out var v) ? v : throw new UsageException($"--{name}: '{s}' is not an integer.");
    }

    /// <summary>A comma-separated list of numbers.</summary>
    public double[]? Numbers(string name)
    {
        var s = Get(name);
        if (s == null) return null;
        return ParseNumbers(s, name);
    }

    public static double[] ParseNumbers(string s, string what) =>
        s.Split([',', ' ', ';'], StringSplitOptions.RemoveEmptyEntries).Select(x =>
            double.TryParse(x, NumberStyles.Float, CultureInfo.InvariantCulture, out var v) ? v : throw new UsageException($"{what}: '{x}' is not a number.")).ToArray();

    /// <summary>
    /// Options every command accepts, whether it needs them or not, so that a script can set them
    /// once for all its commands (the 'options' statement) without a warning from each command.
    /// </summary>
    public static readonly HashSet<string> Common = new(StringComparer.OrdinalIgnoreCase)
    {
        "data", "stations", "events", "picks", "tstar", "quakeml", "stationxml", "hypocentres", "catalog-locations", "grid", "model",
        "vp", "vs", "config", "set", "out", "name", "formats", "csv", "vtk", "register", "adaptive", "straight", "no-opencl", "threads",
        "quiet", "keep-work", "waveforms"
    };

    /// <summary>Options given but never read by the command (other than the common ones): most likely a typing error.</summary>
    public IEnumerable<string> Unused() => _options.Keys.Concat(_flags).Where(k => !_used.Contains(k) && !Common.Contains(k)).Select(k => "--" + k);
}

/// <summary>A command-line mistake: reported with a hint to the help, exit code 2.</summary>
public sealed class UsageException(string message) : Exception(message);
