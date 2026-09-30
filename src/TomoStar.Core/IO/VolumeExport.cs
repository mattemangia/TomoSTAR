using System.Buffers.Binary;
using System.Globalization;
using System.Text;
using TomoStar.Core.Geo;

namespace TomoStar.Core.IO;

/// <summary>
/// Volumes as text and as VTK, in the exact layouts QUIVER exports and imports, so a model written
/// here opens in QUIVER (File, Import volume) as well as in ParaView, GMT or a spreadsheet.
///
/// CSV: comment lines "# QUIVER volume: name; quantity=...; units=..." and "# meta key=value", then
/// the header "lon,lat,depth_km,value" and one row per node, longitude fastest.
///
/// VTK: legacy STRUCTURED_GRID, big-endian binary, points in km in a local east-north-up frame at
/// the centre of the grid (so the true curved shape is shown); the second line carries the grid
/// definition, which is how QUIVER rebuilds the grid on import. Longitude, latitude and depth are
/// added as extra point scalars.
/// </summary>
public static class VolumeExport
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    private const string VtkTag = "QUIVER";

    /// <summary>Writes a node field as a QUIVER-readable CSV.</summary>
    public static void WriteCsv(string path, VolumeHeader header, ReadOnlySpan<double> values)
    {
        var g = new SphericalGrid(header.Grid);
        if (values.Length != g.Count) throw new ArgumentException($"Volume has {values.Length} values, the grid {g.Count}.");
        using var w = new StreamWriter(path, false, new UTF8Encoding(false), 1 << 20);
        w.WriteLine($"# QUIVER volume: {Escape(header.Name)}; quantity={Escape(header.Quantity)}; units={Escape(header.Units)}");
        foreach (var (key, value) in header.Metadata) w.WriteLine($"# meta {Escape(key)}={Escape(value)}");
        w.WriteLine("lon,lat,depth_km,value");
        for (var k = 0; k < g.Nz; k++)
        for (var j = 0; j < g.Ny; j++)
        for (var i = 0; i < g.Nx; i++)
        {
            var v = values[g.Index(i, j, k)];
            w.Write(g.LonDeg[i].ToString("0.######", Inv)); w.Write(',');
            w.Write(g.LatDeg[j].ToString("0.######", Inv)); w.Write(',');
            w.Write(g.DepthKm[k].ToString("0.####", Inv)); w.Write(',');
            w.WriteLine(double.IsFinite(v) ? ((float)v).ToString("R", Inv) : "NaN");
        }
    }

    /// <summary>Writes a node field as a QUIVER-readable legacy VTK structured grid.</summary>
    public static void WriteVtk(string path, VolumeHeader header, ReadOnlySpan<double> values)
    {
        var g = new SphericalGrid(header.Grid);
        if (values.Length != g.Count) throw new ArgumentException($"Volume has {values.Length} values, the grid {g.Count}.");
        var d = g.Definition;
        var frame = g.DisplayFrame();
        using var fs = new FileStream(path, FileMode.Create, FileAccess.Write, FileShare.None, 1 << 20);
        var title = string.Create(Inv,
            $"{VtkTag} volume;lon={d.MinLon},{d.MaxLon},{d.Nx};lat={d.MinLat},{d.MaxLat},{d.Ny};dep={d.MinDepthKm},{d.MaxDepthKm},{d.Nz};q={Escape(header.Quantity, true)};u={Escape(header.Units, true)};name={Escape(header.Name, true)}");
        Ascii(fs, "# vtk DataFile Version 3.0\n");
        Ascii(fs, (title.Length > 255 ? title[..255] : title) + "\n");
        Ascii(fs, "BINARY\nDATASET STRUCTURED_GRID\n");
        Ascii(fs, $"DIMENSIONS {g.Nx} {g.Ny} {g.Nz}\n");
        Ascii(fs, $"POINTS {g.Count} float\n");
        var buf = new byte[12 * g.Nx];
        for (var k = 0; k < g.Nz; k++)
        for (var j = 0; j < g.Ny; j++)
        {
            for (var i = 0; i < g.Nx; i++)
            {
                var p = frame.ToLocal(g.LonDeg[i], g.LatDeg[j], g.DepthKm[k]);
                BinaryPrimitives.WriteSingleBigEndian(buf.AsSpan(12 * i), (float)p.X);
                BinaryPrimitives.WriteSingleBigEndian(buf.AsSpan(12 * i + 4), (float)p.Y);
                BinaryPrimitives.WriteSingleBigEndian(buf.AsSpan(12 * i + 8), (float)p.Z);
            }
            fs.Write(buf);
        }
        Ascii(fs, $"\nPOINT_DATA {g.Count}\n");
        var copy = values.ToArray();
        Scalars(fs, VtkName(header.Quantity), g, (i, j, k) => (float)copy[g.Index(i, j, k)]);
        Scalars(fs, "depth_km", g, (i, j, k) => (float)g.DepthKm[k]);
        Scalars(fs, "lon", g, (i, j, k) => (float)g.LonDeg[i]);
        Scalars(fs, "lat", g, (i, j, k) => (float)g.LatDeg[j]);
    }

    private static void Scalars(Stream fs, string name, SphericalGrid g, Func<int, int, int, float> f)
    {
        Ascii(fs, $"SCALARS {name} float 1\nLOOKUP_TABLE default\n");
        var row = new byte[4 * g.Nx];
        for (var k = 0; k < g.Nz; k++)
        for (var j = 0; j < g.Ny; j++)
        {
            for (var i = 0; i < g.Nx; i++) BinaryPrimitives.WriteSingleBigEndian(row.AsSpan(4 * i), f(i, j, k));
            fs.Write(row);
        }
        Ascii(fs, "\n");
    }

    private static void Ascii(Stream s, string text) => s.Write(Encoding.ASCII.GetBytes(text));

    /// <summary>VTK names are ASCII without blanks: any other character becomes '_'.</summary>
    private static string VtkName(string quantity)
    {
        var s = new string(quantity.Select(c => c < 128 && char.IsLetterOrDigit(c) ? c : '_').ToArray());
        return s.Length == 0 ? "value" : s;
    }

    /// <summary>
    /// Percent-encodes the characters that separate fields in the header lines (';', '=', '%', line
    /// breaks) and, for VTK, anything outside ASCII; QUIVER decodes them on import.
    /// </summary>
    private static string Escape(string s, bool asciiOnly = false)
    {
        var sb = new StringBuilder(s.Length);
        foreach (var rune in s.EnumerateRunes())
        {
            if (rune.Value is ';' or '=' or '%' or '\n' or '\r' || asciiOnly && rune.Value > 127)
            {
                Span<byte> utf8 = stackalloc byte[4];
                var n = rune.EncodeToUtf8(utf8);
                for (var i = 0; i < n; i++) sb.Append('%').Append(utf8[i].ToString("X2", Inv));
            }
            else sb.Append(rune.ToString());
        }
        return sb.ToString();
    }
}
