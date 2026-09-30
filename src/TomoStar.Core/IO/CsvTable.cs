using System.Globalization;

namespace TomoStar.Core.IO;

/// <summary>
/// A small reader of delimited text tables, the input format of TomoSTAR. The first non-comment
/// line is the header; columns are found by name (without regard to case, with aliases), so extra
/// columns and any column order are accepted. Fields are separated by commas, semicolons or tabs
/// (whichever the header uses); lines starting with '#' and blank lines are skipped; a field may be
/// quoted with double quotes. Numbers use the invariant culture (a dot as decimal separator).
/// </summary>
public sealed class CsvTable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    private readonly Dictionary<string, int> _columns = new(StringComparer.OrdinalIgnoreCase);

    private CsvTable(string path) => Path = path;

    /// <summary>The file the table was read from (for error messages).</summary>
    public string Path { get; }

    /// <summary>Column names as written in the header.</summary>
    public IReadOnlyList<string> Header { get; private set; } = [];

    /// <summary>Data rows with the 1-based line number they came from.</summary>
    public List<(int Line, string[] Cells)> Rows { get; } = [];

    /// <summary>Reads a whole table.</summary>
    public static CsvTable Read(string path)
    {
        if (!File.Exists(path)) throw new FileNotFoundException($"Input file not found: {path}", path);
        var table = new CsvTable(path);
        char separator = ',';
        var lineNo = 0;
        foreach (var raw in File.ReadLines(path))
        {
            lineNo++;
            var line = raw.Trim().TrimStart('﻿');
            if (line.Length == 0 || line[0] == '#') continue;
            if (table.Header.Count == 0)
            {
                separator = line.Contains('\t') ? '\t' : line.Contains(';') && !line.Contains(',') ? ';' : ',';
                table.Header = Split(line, separator).Select(h => h.Trim()).ToArray();
                for (var i = 0; i < table.Header.Count; i++) table._columns.TryAdd(table.Header[i], i);
                continue;
            }
            table.Rows.Add((lineNo, Split(line, separator)));
        }
        if (table.Header.Count == 0) throw new InvalidDataException($"{path}: no header line.");
        return table;
    }

    /// <summary>Index of the first column with one of the given names, or -1.</summary>
    public int Find(params string[] names)
    {
        foreach (var n in names)
            if (_columns.TryGetValue(n, out var i)) return i;
        return -1;
    }

    /// <summary>Index of a required column; throws a message naming the accepted names otherwise.</summary>
    public int Require(string what, params string[] names)
    {
        var i = Find(names);
        if (i < 0) throw new InvalidDataException($"{Path}: no {what} column (expected one of: {string.Join(", ", names)}).");
        return i;
    }

    /// <summary>The text of a cell, empty when the row is shorter than the header or the column is absent.</summary>
    public static string Cell(string[] cells, int column) => column >= 0 && column < cells.Length ? cells[column].Trim() : "";

    /// <summary>A number in a cell; <paramref name="fallback"/> when the cell is empty or absent.</summary>
    public double Number(string[] cells, int column, int line, double fallback = double.NaN)
    {
        var s = Cell(cells, column);
        if (s.Length == 0) return fallback;
        if (double.TryParse(s, NumberStyles.Float, Inv, out var v)) return v;
        if (s.Equals("nan", StringComparison.OrdinalIgnoreCase)) return double.NaN;
        throw new InvalidDataException($"{Path}, line {line}: '{s}' is not a number.");
    }

    /// <summary>
    /// A UTC time in a cell: ISO 8601 (2016-10-30T06:40:17.32Z, with or without the Z, a space
    /// instead of the T accepted) or, when the text is a plain number, seconds since the Unix epoch.
    /// </summary>
    public DateTime Time(string[] cells, int column, int line)
    {
        var s = Cell(cells, column);
        if (double.TryParse(s, NumberStyles.Float, Inv, out var epoch))
            return DateTime.UnixEpoch.AddTicks((long)Math.Round(epoch * TimeSpan.TicksPerSecond));
        if (DateTime.TryParse(s, Inv, DateTimeStyles.AdjustToUniversal | DateTimeStyles.AssumeUniversal, out var t))
            return DateTime.SpecifyKind(t, DateTimeKind.Utc);
        throw new InvalidDataException($"{Path}, line {line}: '{s}' is not a time (use ISO 8601, e.g. 2016-10-30T06:40:17.32Z).");
    }

    /// <summary>A yes/no cell: 1, true, yes, y (any case) are true; empty is <paramref name="fallback"/>.</summary>
    public static bool Flag(string[] cells, int column, bool fallback = false)
    {
        var s = Cell(cells, column);
        if (s.Length == 0) return fallback;
        return s is "1" || s.Equals("true", StringComparison.OrdinalIgnoreCase) || s.Equals("yes", StringComparison.OrdinalIgnoreCase)
               || s.Equals("y", StringComparison.OrdinalIgnoreCase);
    }

    /// <summary>Splits a line on a separator, honouring double quotes ("" inside quotes is a quote).</summary>
    private static string[] Split(string line, char separator)
    {
        if (!line.Contains('"')) return line.Split(separator);
        var cells = new List<string>();
        var sb = new System.Text.StringBuilder();
        var quoted = false;
        for (var i = 0; i < line.Length; i++)
        {
            var c = line[i];
            if (c == '"')
            {
                if (quoted && i + 1 < line.Length && line[i + 1] == '"') { sb.Append('"'); i++; }
                else quoted = !quoted;
            }
            else if (c == separator && !quoted) { cells.Add(sb.ToString()); sb.Clear(); }
            else sb.Append(c);
        }
        cells.Add(sb.ToString());
        return cells.ToArray();
    }
}
