// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using System.Globalization;
using System.Text;
using System.Text.Json;
using System.Text.Json.Nodes;
using System.Text.RegularExpressions;

namespace TomoStar.Cli;

/// <summary>
/// The interpreter of TomoSTAR scripts (<c>tomostar run FILE</c>): a whole processing sequence in one
/// text file, so that a study can be repeated exactly and its script published with the paper.
///
/// A script is a list of statements, one per line; a line ending in a backslash continues on the
/// next; lines starting with '#' are comments. Every tomostar command can be written as a line
/// (without the word tomostar), and these statements control the sequence:
///
/// <code>
/// set NAME = value          define a variable (used as ${NAME})
/// default NAME = value      define it only if not defined yet (so 'run --var NAME=...' wins)
/// config FILE               configuration file of the following commands ('config none' to clear)
/// options --opt value ...   options added to every following command ('options none' to clear)
/// echo text                 print a message
/// include FILE              run another script here, sharing the variables
/// foreach NAME in A B C     repeat the lines up to the matching 'end' for each value
/// if exists PATH            run the lines up to 'end' only if the file or folder exists
/// if missing PATH           ... only if it does not (to skip steps already done)
/// end                       closes foreach and if
/// exit                      stop the script here
/// </code>
///
/// Values are expanded before a line runs: <c>${NAME}</c> a variable, <c>${env:NAME}</c> an environment
/// variable, <c>${json:FILE:Path.To.Value}</c> a value read from a JSON file (for example the damping
/// recommended by an L-curve run, <c>${json:${OUT}/lcurve/summary.json:Recommended.Damping}</c>).
/// A command line starting with '-' may fail without stopping the script. Relative paths are relative
/// to the folder of the script. The predefined variables are SCRIPT_DIR, CALLER_DIR, DATE (yyyyMMdd)
/// and TIME (HHmmss).
/// </summary>
public sealed partial class ScriptRunner
{
    private readonly Dictionary<string, string> _vars = new(StringComparer.Ordinal);
    private readonly List<string> _options = [];
    private readonly CancellationToken _ct;
    private readonly bool _dryRun;
    private readonly bool _keepGoing;
    private string? _config;
    private bool _exit;
    private int _failures;

    public ScriptRunner(IEnumerable<string> assignments, bool dryRun, bool keepGoing, CancellationToken ct)
    {
        _dryRun = dryRun;
        _keepGoing = keepGoing;
        _ct = ct;
        _vars["CALLER_DIR"] = Environment.CurrentDirectory;
        _vars["DATE"] = DateTime.Now.ToString("yyyyMMdd", CultureInfo.InvariantCulture);
        _vars["TIME"] = DateTime.Now.ToString("HHmmss", CultureInfo.InvariantCulture);
        foreach (var a in assignments)
        {
            var eq = a.IndexOf('=');
            if (eq <= 0) throw new UsageException($"--var needs NAME=value, got '{a}'.");
            _vars[a[..eq].Trim()] = a[(eq + 1)..].Trim();
        }
    }

    private sealed record Line(string File, int Number, string Text);

    /// <summary>Runs a script file; returns 0 when every command succeeded.</summary>
    public int Run(string path)
    {
        var full = Path.GetFullPath(path);
        if (!File.Exists(full)) throw new UsageException($"Script not found: {path}");
        var previous = Environment.CurrentDirectory;
        var dir = Path.GetDirectoryName(full)!;
        _vars["SCRIPT_DIR"] = dir;
        Environment.CurrentDirectory = dir;
        try
        {
            Execute(Read(full), 0, int.MaxValue);
        }
        finally
        {
            Environment.CurrentDirectory = previous;
        }
        if (_failures > 0) Console.WriteLine($"Script finished with {_failures} failed command(s).");
        return _failures > 0 ? 1 : 0;
    }

    private static List<Line> Read(string file)
    {
        var lines = new List<Line>();
        var raw = File.ReadAllLines(file);
        for (var i = 0; i < raw.Length; i++)
        {
            var number = i + 1;
            var text = raw[i].TrimEnd();
            while (text.EndsWith('\\') && i + 1 < raw.Length) text = text[..^1] + " " + raw[++i].Trim();
            var t = text.Trim();
            if (t.Length == 0 || t[0] == '#') continue;
            lines.Add(new Line(file, number, t));
        }
        return lines;
    }

    /// <summary>Runs lines [from, to) of a block; handles foreach, if, end, exit.</summary>
    private void Execute(List<Line> lines, int from, int to)
    {
        to = Math.Min(to, lines.Count);
        for (var i = from; i < to && !_exit; i++)
        {
            _ct.ThrowIfCancellationRequested();
            var line = lines[i];
            var word = FirstWord(line.Text);
            switch (word)
            {
                case "foreach":
                {
                    var end = MatchingEnd(lines, i);
                    var m = ForeachPattern().Match(line.Text);
                    if (!m.Success) throw Error(line, "foreach NAME in VALUE VALUE ...");
                    var values = Tokenize(Expand(m.Groups[2].Value, line));
                    foreach (var v in values)
                    {
                        _vars[m.Groups[1].Value] = v;
                        Execute(lines, i + 1, end);
                        if (_exit) break;
                    }
                    i = end;
                    break;
                }
                case "if":
                {
                    var end = MatchingEnd(lines, i);
                    var parts = Tokenize(Expand(line.Text, line));
                    if (parts.Count != 3 || parts[1] is not ("exists" or "missing")) throw Error(line, "if exists PATH / if missing PATH");
                    var exists = File.Exists(parts[2]) || Directory.Exists(parts[2]);
                    if (exists == (parts[1] == "exists")) Execute(lines, i + 1, end);
                    i = end;
                    break;
                }
                case "end":
                    throw Error(line, "'end' without foreach or if.");
                default:
                    Statement(line);
                    break;
            }
        }
    }

    private static int MatchingEnd(List<Line> lines, int start)
    {
        var depth = 0;
        for (var i = start; i < lines.Count; i++)
        {
            var w = FirstWord(lines[i].Text);
            if (w is "foreach" or "if") depth++;
            else if (w == "end" && --depth == 0) return i;
        }
        throw Error(lines[start], $"'{FirstWord(lines[start].Text)}' has no matching 'end'.");
    }

    private void Statement(Line line)
    {
        var word = FirstWord(line.Text);
        var rest = line.Text[word.Length..].Trim();
        switch (word)
        {
            case "set":
            case "default":
            {
                var eq = rest.IndexOf('=');
                if (eq <= 0) throw Error(line, $"{word} NAME = value");
                var name = rest[..eq].Trim();
                if (word == "default" && _vars.ContainsKey(name)) return;
                var value = Expand(rest[(eq + 1)..].Trim(), line);
                if (value.Length >= 2 && value[0] == '"' && value[^1] == '"') value = value[1..^1];
                _vars[name] = value;
                return;
            }
            case "echo":
                Console.WriteLine(Expand(rest, line));
                return;
            case "config":
            {
                var v = Expand(rest, line).Trim('"');
                _config = v.Equals("none", StringComparison.OrdinalIgnoreCase) ? null : Path.GetFullPath(v);
                return;
            }
            case "options":
            {
                var v = Expand(rest, line);
                _options.Clear();
                if (!v.Trim().Equals("none", StringComparison.OrdinalIgnoreCase)) _options.AddRange(Tokenize(v));
                return;
            }
            case "include":
            {
                var file = Path.GetFullPath(Expand(rest, line).Trim('"'));
                if (!File.Exists(file)) throw Error(line, $"include: {file} not found.");
                Execute(Read(file), 0, int.MaxValue);
                return;
            }
            case "exit":
                _exit = true;
                return;
        }
        // A tomostar command.
        var tolerant = line.Text.StartsWith('-') && !line.Text.StartsWith("--", StringComparison.Ordinal);
        var tokens = Tokenize(Expand(tolerant ? line.Text[1..] : line.Text, line));
        if (tokens.Count > 0 && tokens[0].Equals("tomostar", StringComparison.OrdinalIgnoreCase)) tokens.RemoveAt(0);
        if (tokens.Count == 0) return;
        if (tokens[0] == "run") throw Error(line, "use 'include' to run another script from a script.");
        var args = new List<string> { tokens[0] };
        if (_config != null && !tokens.Contains("--config")) args.AddRange(["--config", _config]);
        args.AddRange(_options);
        args.AddRange(tokens.Skip(1));
        Console.WriteLine($"> tomostar {string.Join(' ', args.Select(Quote))}");
        if (_dryRun) return;
        var code = Program.Execute(args, _ct);
        if (code == 0) return;
        _failures++;
        if (!tolerant && !_keepGoing) throw new ScriptFailedException($"{Path.GetFileName(line.File)}, line {line.Number}: the command failed (exit code {code}); the script stops.");
    }

    // ---- Expansion and tokens --------------------------------------------------------------------

    /// <summary>Replaces ${...} references, innermost first (so a JSON path may itself contain ${OUT}).</summary>
    private string Expand(string text, Line line)
    {
        for (var guard = 0; guard < 50; guard++)
        {
            var m = ReferencePattern().Match(text);
            if (!m.Success) return text;
            text = text[..m.Index] + Resolve(m.Groups[1].Value, line) + text[(m.Index + m.Length)..];
        }
        throw Error(line, "too many nested ${...} references.");
    }

    private string Resolve(string reference, Line line)
    {
        if (reference.StartsWith("env:", StringComparison.Ordinal))
            return Environment.GetEnvironmentVariable(reference[4..]) ?? throw Error(line, $"environment variable {reference[4..]} is not set.");
        if (reference.StartsWith("json:", StringComparison.Ordinal))
        {
            var spec = reference[5..];
            var colon = spec.LastIndexOf(':');
            if (colon <= 0) throw Error(line, "${json:FILE:Path.To.Value}");
            var file = spec[..colon];
            if (_dryRun && !File.Exists(file)) return $"<{spec[(colon + 1)..]} from {Path.GetFileName(file)}>";
            if (!File.Exists(file)) throw Error(line, $"{file} not found (was the step that writes it run?).");
            JsonNode? node = JsonNode.Parse(File.ReadAllText(file), documentOptions: new JsonDocumentOptions { CommentHandling = JsonCommentHandling.Skip });
            foreach (var part in spec[(colon + 1)..].Split('.', StringSplitOptions.RemoveEmptyEntries))
            {
                node = node switch
                {
                    JsonObject o => o.FirstOrDefault(kv => string.Equals(kv.Key, part, StringComparison.OrdinalIgnoreCase)).Value,
                    JsonArray a when int.TryParse(part, out var k) && k >= 0 && k < a.Count => a[k],
                    _ => null
                };
                if (node == null) throw Error(line, $"{file} has no value at '{spec[(colon + 1)..]}'.");
            }
            return node is JsonValue v && v.TryGetValue<double>(out var d) ? d.ToString("R", CultureInfo.InvariantCulture)
                : node is JsonValue s && s.TryGetValue<string>(out var str) ? str : node!.ToJsonString();
        }
        return _vars.TryGetValue(reference, out var value) ? value : throw Error(line, $"variable {reference} is not defined.");
    }

    /// <summary>Splits a line into words; double or single quotes group words, a backslash escapes a quote.</summary>
    public static List<string> Tokenize(string text)
    {
        var tokens = new List<string>();
        var sb = new StringBuilder();
        char quote = '\0';
        var inToken = false;
        for (var i = 0; i < text.Length; i++)
        {
            var ch = text[i];
            if (quote != '\0')
            {
                if (ch == '\\' && i + 1 < text.Length && text[i + 1] == quote) { sb.Append(quote); i++; }
                else if (ch == quote) quote = '\0';
                else sb.Append(ch);
            }
            else if (ch is '"' or '\'') { quote = ch; inToken = true; }
            else if (char.IsWhiteSpace(ch))
            {
                if (inToken) { tokens.Add(sb.ToString()); sb.Clear(); inToken = false; }
            }
            else { sb.Append(ch); inToken = true; }
        }
        if (quote != '\0') throw new UsageException($"Unclosed quote in: {text}");
        if (inToken) tokens.Add(sb.ToString());
        return tokens;
    }

    private static string Quote(string s) => s.Length == 0 || s.Any(char.IsWhiteSpace) ? $"\"{s}\"" : s;

    private static string FirstWord(string text)
    {
        var i = 0;
        while (i < text.Length && !char.IsWhiteSpace(text[i])) i++;
        return text[..i];
    }

    private static UsageException Error(Line line, string message) => new($"{Path.GetFileName(line.File)}, line {line.Number}: {message}");

    [GeneratedRegex(@"\$\{([^${}]+)\}")]
    private static partial Regex ReferencePattern();

    [GeneratedRegex(@"^foreach\s+([A-Za-z_][A-Za-z0-9_]*)\s+in\s+(.+)$")]
    private static partial Regex ForeachPattern();
}

/// <summary>A command of a script failed and the script stopped.</summary>
public sealed class ScriptFailedException(string message) : Exception(message);
