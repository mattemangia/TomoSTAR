// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using System.Buffers.Binary;
using System.IO.MemoryMappedFiles;
using System.Text;
using System.Text.Json;
using TomoStar.Core.Geo;
using TomoStar.Core.Model;

namespace TomoStar.Core.IO;

/// <summary>Header of a <c>.qvol</c> file.</summary>
public sealed class VolumeHeader
{
    public string Name { get; set; } = "";
    public string Quantity { get; set; } = "";
    public string Units { get; set; } = "";
    public GridDefinition Grid { get; set; } = new();
    public float NoData { get; set; } = float.NaN;
    public double Min { get; set; } = double.NaN;
    public double Max { get; set; } = double.NaN;
    public string Source { get; set; } = "";
    public DateTime CreatedUtc { get; set; } = DateTime.UtcNow;

    /// <summary>Free-form provenance (tool, parameters) for a figure caption.</summary>
    public Dictionary<string, string> Metadata { get; set; } = [];
}

/// <summary>
/// The QUIVER volume format: an 8-byte magic, a little-endian int32 header length, a UTF-8 JSON
/// <see cref="VolumeHeader"/>, padding to a 64-byte boundary, then Nx·Ny·Nz little-endian float32
/// values with x fastest and depth slowest.
///
/// The layout is what makes streaming cheap: a depth slice is one contiguous read, a vertical
/// slice is Nz·Ny (or Nz·Nx) short reads, and nothing else is ever touched. Readers map the file
/// instead of loading it, so a multi-gigabyte model opens on a small workstation.
/// </summary>
public static class VolumeFile
{
    public const string Extension = ".qvol";
    private static readonly byte[] Magic = "QVOL\u0001\0\0\0"u8.ToArray();

    public static void Write(string path, VolumeHeader header, ReadOnlySpan<double> values)
    {
        var count = header.Grid.Count;
        if (values.Length != count) throw new ArgumentException($"Volume has {values.Length} values, the grid {count}.");
        var min = double.PositiveInfinity;
        var max = double.NegativeInfinity;
        foreach (var v in values)
        {
            if (!double.IsFinite(v)) continue;
            if (v < min) min = v;
            if (v > max) max = v;
        }
        header.Min = double.IsFinite(min) ? min : double.NaN;
        header.Max = double.IsFinite(max) ? max : double.NaN;

        var dir = Path.GetDirectoryName(Path.GetFullPath(path));
        if (dir != null) Directory.CreateDirectory(dir);
        var tmp = path + ".tmp";
        using (var fs = new FileStream(tmp, FileMode.Create, FileAccess.Write, FileShare.None, 1 << 20))
        {
            var json = JsonSerializer.SerializeToUtf8Bytes(header, TomoJson.Options);
            fs.Write(Magic);
            Span<byte> len = stackalloc byte[4];
            BinaryPrimitives.WriteInt32LittleEndian(len, json.Length);
            fs.Write(len);
            fs.Write(json);
            var pad = (int)((64 - fs.Position % 64) % 64);
            fs.Write(new byte[pad]);

            var buffer = new byte[4 * 65536];
            var n = 0;
            foreach (var v in values)
            {
                BinaryPrimitives.WriteSingleLittleEndian(buffer.AsSpan(n), double.IsFinite(v) ? (float)v : header.NoData);
                n += 4;
                if (n == buffer.Length) { fs.Write(buffer, 0, n); n = 0; }
            }
            if (n > 0) fs.Write(buffer, 0, n);
        }
        File.Move(tmp, path, overwrite: true);
    }

    public static VolumeHeader ReadHeader(string path) => ReadHeader(path, out _);

    internal static VolumeHeader ReadHeader(string path, out long dataOffset)
    {
        using var fs = File.OpenRead(path);
        Span<byte> head = stackalloc byte[12];
        if (fs.Read(head) != 12 || !head[..8].SequenceEqual(Magic)) throw new InvalidDataException($"{path} is not a QUIVER volume.");
        var len = BinaryPrimitives.ReadInt32LittleEndian(head[8..]);
        if (len <= 0 || len > 16 << 20) throw new InvalidDataException("Corrupt volume header length.");
        var json = new byte[len];
        fs.ReadExactly(json);
        var header = JsonSerializer.Deserialize<VolumeHeader>(json, TomoJson.Options) ?? throw new InvalidDataException("Empty volume header.");
        var end = 12 + len;
        dataOffset = end + (64 - end % 64) % 64;
        if (fs.Length < dataOffset + header.Grid.Count * 4) throw new InvalidDataException("Volume file is truncated.");
        return header;
    }
}

/// <summary>
/// Streaming reader of a <c>.qvol</c>: memory-mapped, read-only, thread-safe. Slices come out as
/// fresh arrays; the whole volume is never materialised unless <see cref="ReadAll"/> is called.
/// A small LRU keeps the last few slices so dragging a slider does not re-read the disk.
/// </summary>
public sealed class VolumeReader : IDisposable
{
    private readonly MemoryMappedFile _file;
    private readonly MemoryMappedViewAccessor _view;
    private readonly long _offset;
    private readonly object _cacheGate = new();
    private readonly LinkedList<(string Key, float[] Data)> _cache = new();
    private const int CacheSlots = 12;

    public VolumeReader(string path)
    {
        Path = System.IO.Path.GetFullPath(path);
        Header = VolumeFile.ReadHeader(path, out _offset);
        Grid = new SphericalGrid(Header.Grid);
        _file = MemoryMappedFile.CreateFromFile(Path, FileMode.Open, null, 0, MemoryMappedFileAccess.Read);
        _view = _file.CreateViewAccessor(0, 0, MemoryMappedFileAccess.Read);
    }

    public string Path { get; }
    public VolumeHeader Header { get; }
    public SphericalGrid Grid { get; }

    public float Value(int i, int j, int k)
    {
        var v = _view.ReadSingle(_offset + 4L * Grid.Index(i, j, k));
        return IsNoData(v) ? float.NaN : v;
    }

    private bool IsNoData(float v) => float.IsNaN(v) || (!float.IsNaN(Header.NoData) && v == Header.NoData);

    /// <summary>Depth slice k as [j·Nx + i]: one contiguous read.</summary>
    public float[] DepthSlice(int k) => Cached($"k{k}", () =>
    {
        var n = Grid.Nx * Grid.Ny;
        var data = new float[n];
        _view.ReadArray(_offset + 4L * k * n, data, 0, n);
        Clean(data);
        return data;
    });

    /// <summary>Longitude slice i as [k·Ny + j].</summary>
    public float[] LonSlice(int i) => Cached($"i{i}", () =>
    {
        var data = new float[Grid.Ny * Grid.Nz];
        for (var k = 0; k < Grid.Nz; k++)
        for (var j = 0; j < Grid.Ny; j++)
            data[k * Grid.Ny + j] = _view.ReadSingle(_offset + 4L * Grid.Index(i, j, k));
        Clean(data);
        return data;
    });

    /// <summary>Latitude slice j as [k·Nx + i].</summary>
    public float[] LatSlice(int j) => Cached($"j{j}", () =>
    {
        var data = new float[Grid.Nx * Grid.Nz];
        for (var k = 0; k < Grid.Nz; k++)
            _view.ReadArray(_offset + 4L * Grid.Index(0, j, k), data, k * Grid.Nx, Grid.Nx);
        Clean(data);
        return data;
    });

    /// <summary>
    /// Row (j, k) of nodes along longitude into <paramref name="buffer"/> (length ≥ Nx): one contiguous
    /// read, not cached, for callers that walk a large volume once (e.g. decimated surface extraction).
    /// </summary>
    public void ReadRow(int j, int k, float[] buffer)
    {
        _view.ReadArray(_offset + 4L * Grid.Index(0, j, k), buffer, 0, Grid.Nx);
        if (float.IsNaN(Header.NoData)) return;
        for (var n = 0; n < Grid.Nx; n++)
            if (buffer[n] == Header.NoData) buffer[n] = float.NaN;
    }

    /// <summary>Trilinear sample at a geographic point; NaN outside the grid or next to no-data.</summary>
    public double Sample(double lon, double lat, double depthKm)
    {
        if (!Grid.Contains(lon, lat, depthKm, 1e-6)) return double.NaN;
        Span<int> nodes = stackalloc int[8];
        Span<double> w = stackalloc double[8];
        Grid.TrilinearWeights(lon, lat, depthKm, nodes, w);
        double sum = 0, wsum = 0;
        for (var n = 0; n < 8; n++)
        {
            if (w[n] == 0) continue;
            var v = _view.ReadSingle(_offset + 4L * nodes[n]);
            if (IsNoData(v)) continue;
            sum += v * w[n];
            wsum += w[n];
        }
        return wsum > 0.5 ? sum / wsum : double.NaN;
    }

    /// <summary>Whole volume as doubles. Only for tools that genuinely need all of it.</summary>
    public double[] ReadAll()
    {
        var n = Grid.Count;
        var f = new float[n];
        _view.ReadArray(_offset, f, 0, n);
        var d = new double[n];
        for (var i = 0; i < n; i++) d[i] = IsNoData(f[i]) ? double.NaN : f[i];
        return d;
    }

    private void Clean(float[] data)
    {
        if (float.IsNaN(Header.NoData)) return;
        for (var n = 0; n < data.Length; n++)
            if (data[n] == Header.NoData) data[n] = float.NaN;
    }

    private float[] Cached(string key, Func<float[]> load)
    {
        lock (_cacheGate)
        {
            for (var node = _cache.First; node != null; node = node.Next)
            {
                if (node.Value.Key != key) continue;
                _cache.Remove(node);
                _cache.AddFirst(node);
                return node.Value.Data;
            }
        }
        var data = load();
        lock (_cacheGate)
        {
            _cache.AddFirst((key, data));
            while (_cache.Count > CacheSlots) _cache.RemoveLast();
        }
        return data;
    }

    public void Dispose()
    {
        _view.Dispose();
        _file.Dispose();
    }
}
