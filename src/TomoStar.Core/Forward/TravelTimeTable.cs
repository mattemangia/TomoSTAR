// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using System.Buffers.Binary;
using System.IO.MemoryMappedFiles;
using System.Security.Cryptography;
using System.Text.Json;
using TomoStar.Core.Compute;
using TomoStar.Core.Geo;
using TomoStar.Core.Model;

namespace TomoStar.Core.Forward;

public sealed class TravelTimeTableHeader
{
    public string StationId { get; set; } = "";
    public Phase Phase { get; set; }
    public double Lon { get; set; }
    public double Lat { get; set; }
    public double DepthKm { get; set; }
    public GridDefinition Grid { get; set; } = new();

    /// <summary>Hash of the slowness field, so a table is reused only for the model it was computed in.</summary>
    public string ModelKey { get; set; } = "";

    public string Solver { get; set; } = "";
}

/// <summary>
/// The first-arrival travel-time field from one station for one phase, stored on disk as float32
/// and read through a memory map. By reciprocity it gives the time from ANY source in the grid to
/// the station, so relocation, pick prediction and ray tracing all read the same table and none of
/// them holds hundreds of fields in memory (the approach of NonLinLoc time grids; Lomax et al.
/// 2000, in Advances in Seismic Event Location, Kluwer, 101-134).
/// </summary>
public sealed class TravelTimeTable : IDisposable
{
    private static readonly byte[] Magic = "QTTB\u0001\0\0\0"u8.ToArray();
    private readonly MemoryMappedFile _file;
    private readonly MemoryMappedViewAccessor _view;
    private readonly long _offset;

    private TravelTimeTable(string path, TravelTimeTableHeader header, long offset)
    {
        Path = path;
        Header = header;
        Grid = new SphericalGrid(header.Grid);
        _offset = offset;
        _file = MemoryMappedFile.CreateFromFile(path, FileMode.Open, null, 0, MemoryMappedFileAccess.Read);
        _view = _file.CreateViewAccessor(0, 0, MemoryMappedFileAccess.Read);
    }

    public string Path { get; }
    public TravelTimeTableHeader Header { get; }
    public SphericalGrid Grid { get; }

    public static void Write(string path, TravelTimeTableHeader header, ReadOnlySpan<double> times)
    {
        var tmp = path + ".tmp";
        using (var fs = new FileStream(tmp, FileMode.Create, FileAccess.Write, FileShare.None, 1 << 20))
        {
            var json = JsonSerializer.SerializeToUtf8Bytes(header, TomoJson.Options);
            fs.Write(Magic);
            Span<byte> len = stackalloc byte[4];
            BinaryPrimitives.WriteInt32LittleEndian(len, json.Length);
            fs.Write(len);
            fs.Write(json);
            fs.Write(new byte[(64 - (12 + json.Length) % 64) % 64]);
            var buf = new byte[1 << 18];
            var o = 0;
            foreach (var t in times)
            {
                BinaryPrimitives.WriteSingleLittleEndian(buf.AsSpan(o), double.IsFinite(t) ? (float)t : float.PositiveInfinity);
                o += 4;
                if (o == buf.Length) { fs.Write(buf); o = 0; }
            }
            if (o > 0) fs.Write(buf, 0, o);
        }
        File.Move(tmp, path, overwrite: true);
    }

    public static TravelTimeTableHeader? ReadHeader(string path)
    {
        try
        {
            using var fs = File.OpenRead(path);
            Span<byte> head = stackalloc byte[12];
            if (fs.Read(head) != 12 || !head[..8].SequenceEqual(Magic)) return null;
            var json = new byte[BinaryPrimitives.ReadInt32LittleEndian(head[8..])];
            fs.ReadExactly(json);
            return JsonSerializer.Deserialize<TravelTimeTableHeader>(json, TomoJson.Options);
        }
        catch (Exception ex) when (ex is IOException or JsonException)
        {
            return null;
        }
    }

    public static TravelTimeTable Open(string path)
    {
        using var fs = File.OpenRead(path);
        Span<byte> head = stackalloc byte[12];
        fs.ReadExactly(head);
        if (!head[..8].SequenceEqual(Magic)) throw new InvalidDataException("Not a travel-time table.");
        var len = BinaryPrimitives.ReadInt32LittleEndian(head[8..]);
        var json = new byte[len];
        fs.ReadExactly(json);
        var header = JsonSerializer.Deserialize<TravelTimeTableHeader>(json, TomoJson.Options)!;
        var offset = 12 + len + (64 - (12 + len) % 64) % 64;
        return new TravelTimeTable(path, header, offset);
    }

    public float Node(int index) => _view.ReadSingle(_offset + 4L * index);

    /// <summary>Travel time from a point to the station, s (trilinear; clamped to the grid).</summary>
    public double Time(double lon, double lat, double depthKm)
    {
        var (fx, fy, fz) = Grid.Fractional(lon, lat, depthKm);
        return TimeFractional(fx, fy, fz);
    }

    public double TimeFractional(double fx, double fy, double fz)
    {
        Span<int> nodes = stackalloc int[8];
        Span<double> w = stackalloc double[8];
        Grid.TrilinearWeightsFractional(fx, fy, fz, nodes, w);
        double t = 0;
        for (var n = 0; n < 8; n++)
        {
            if (w[n] == 0) continue;
            var v = Node(nodes[n]);
            if (!float.IsFinite(v)) return double.PositiveInfinity;
            t += v * w[n];
        }
        return t;
    }

    /// <summary>
    /// Gradient of the travel time at a point in the local (east, north, down) frame, s/km. Its
    /// magnitude is the slowness and it points away from the station along the ray; the hypocentre
    /// partial derivatives of a travel time are exactly this vector (Geiger 1912; Thurber 1983).
    /// </summary>
    public (double DtDe, double DtDn, double DtDz) Gradient(double lon, double lat, double depthKm)
    {
        var g = Grid;
        var (fx, fy, fz) = g.Fractional(lon, lat, depthKm);
        // Central differences of the trilinear interpolant, half a node each side.
        const double d = 0.5;
        var tx = (TimeFractional(fx + d, fy, fz) - TimeFractional(fx - d, fy, fz)) / (2 * d);
        var ty = (TimeFractional(fx, fy + d, fz) - TimeFractional(fx, fy - d, fz)) / (2 * d);
        var tz = (TimeFractional(fx, fy, fz + d) - TimeFractional(fx, fy, fz - d)) / (2 * d);
        var r = GeoMath.EarthRadiusKm - depthKm;
        var hx = r * Math.Cos(lat * GeoMath.Deg2Rad) * g.DLonRad;
        var hy = r * g.DLatRad;
        var hz = g.DDepthKm;
        return (tx / hx, ty / hy, tz / hz);
    }

    public void Dispose()
    {
        _view.Dispose();
        _file.Dispose();
    }
}

/// <summary>A velocity model on the forward grid, as slowness per phase.</summary>
public sealed class ForwardModel
{
    public required SphericalGrid Grid { get; init; }
    public required GridMetric Metric { get; init; }
    public required double[] SlownessP { get; init; }
    public required double[] SlownessS { get; init; }

    public double[] Slowness(Phase phase) => phase == Phase.P ? SlownessP : SlownessS;

    public string Key(Phase phase)
    {
        var s = Slowness(phase);
        var bytes = new byte[s.Length * 8];
        Buffer.BlockCopy(s, 0, bytes, 0, bytes.Length);
        var def = JsonSerializer.SerializeToUtf8Bytes(Grid.Definition, TomoJson.Options);
        return Convert.ToHexString(SHA256.HashData([.. bytes, .. def]))[..32];
    }

    /// <summary>Forward model from a 1-D model.</summary>
    public static ForwardModel From1D(GridDefinition forwardGrid, VelocityModel1D model)
    {
        var g = new SphericalGrid(forwardGrid);
        var sp = new double[g.Count];
        var ss = new double[g.Count];
        for (var k = 0; k < g.Nz; k++)
        {
            var vp = model.Vp(g.DepthKm[k]);
            var vs = model.Vs(g.DepthKm[k]);
            var n = g.Nx * g.Ny;
            for (var q = 0; q < n; q++)
            {
                sp[k * n + q] = 1 / vp;
                ss[k * n + q] = 1 / vs;
            }
        }
        return new ForwardModel { Grid = g, Metric = new GridMetric(g), SlownessP = sp, SlownessS = ss };
    }

    /// <summary>Forward model from velocities on the (coarser) inversion grid.</summary>
    public static ForwardModel FromInversionGrid(SphericalGrid inversion, double[] vp, double[] vs, GridDefinition forwardGrid) =>
        FromInversionGrid(inversion, vp, vs, new ForwardDomain(forwardGrid));

    /// <summary>
    /// Forward model on a domain that may extend beyond the inversion grid: interpolated from the
    /// inversion grid inside it, the domain's 1-D background outside.
    /// </summary>
    public static ForwardModel FromInversionGrid(SphericalGrid inversion, double[] vp, double[] vs, ForwardDomain domain)
    {
        var g = new SphericalGrid(domain.Forward);
        var fp = g.ResampleFrom(inversion, vp);
        var fs = g.ResampleFrom(inversion, vs);
        if (domain.Extended && domain.Background is { } bg)
        {
            var tol = 1e-6;
            Parallel.For(0, g.Nz, k =>
            {
                var bp = bg.Vp(g.DepthKm[k]);
                var bs = bg.Vs(g.DepthKm[k]);
                for (var j = 0; j < g.Ny; j++)
                for (var i = 0; i < g.Nx; i++)
                {
                    if (inversion.Contains(g.LonDeg[i], g.LatDeg[j], g.DepthKm[k], tol)) continue;
                    var n = g.Index(i, j, k);
                    fp[n] = bp;
                    fs[n] = bs;
                }
            });
        }
        for (var i = 0; i < fp.Length; i++)
        {
            fp[i] = 1 / fp[i];
            fs[i] = 1 / fs[i];
        }
        return new ForwardModel { Grid = g, Metric = new GridMetric(g), SlownessP = fp, SlownessS = fs };
    }
}

/// <summary>A station to compute a table for.</summary>
public sealed record TableSource(string StationId, double Lon, double Lat, double DepthKm);

/// <summary>
/// Builds (or reuses) the travel-time tables of a set of stations in a forward model, on the
/// OpenCL device when one passed its self-test, otherwise with fast marching on all CPU cores.
/// </summary>
public sealed class TravelTimeTableSet : IDisposable
{
    private readonly Dictionary<(string, Phase), TravelTimeTable> _tables = new();

    public string Folder { get; }

    private TravelTimeTableSet(string folder) => Folder = folder;

    public TravelTimeTable? Get(string stationId, Phase phase) =>
        _tables.TryGetValue((stationId.ToUpperInvariant(), phase), out var t) ? t : null;

    public IEnumerable<TravelTimeTable> All => _tables.Values;

    public static TravelTimeTableSet Build(
        ForwardModel model, IReadOnlyList<TableSource> sources, IReadOnlyList<Phase> phases, string folder,
        EikonalOpenCl? gpu, IProgress<(double, string)>? progress = null, Action<string>? log = null, CancellationToken ct = default)
    {
        Directory.CreateDirectory(folder);
        var set = new TravelTimeTableSet(folder);
        var jobs = new List<(TableSource Src, Phase Phase, string Path, string Key)>();
        foreach (var phase in phases)
        {
            var key = model.Key(phase);
            foreach (var s in sources)
            {
                var path = System.IO.Path.Combine(folder, $"{Safe(s.StationId)}.{phase}.qttb");
                var h = TravelTimeTable.ReadHeader(path);
                if (h != null && h.ModelKey == key && Math.Abs(h.Lon - s.Lon) < 1e-9 && Math.Abs(h.Lat - s.Lat) < 1e-9
                    && Math.Abs(h.DepthKm - s.DepthKm) < 1e-9) continue; // reuse: same model, same station
                jobs.Add((s, phase, path, key));
            }
        }
        var total = Math.Max(1, jobs.Count);
        log?.Invoke($"Travel-time tables: {jobs.Count} to compute, {sources.Count * phases.Count - jobs.Count} reused.");
        if (jobs.Count > 0)
        {
            if (MpiSession.Current is { CanDispatch: true } mpi)
            {
                ct.ThrowIfCancellationRequested();
                log?.Invoke($"MPI travel-time tables: {mpi.Size} ranks.");
                mpi.Dispatch(new MpiSession.Request
                {
                    Kind = "tables", Grid = model.Grid.Definition, P = model.SlownessP, S = model.SlownessS,
                    Sources = sources.ToList(), Phases = phases.ToList(), Folder = Path.GetFullPath(folder), UseGpu = gpu != null
                });
                ct.ThrowIfCancellationRequested();
            }
            else Compute(model, jobs, gpu, total, progress, log, ct);
        }

        foreach (var phase in phases)
        foreach (var s in sources)
        {
            var path = System.IO.Path.Combine(folder, $"{Safe(s.StationId)}.{phase}.qttb");
            if (File.Exists(path)) set._tables[(s.StationId.ToUpperInvariant(), phase)] = TravelTimeTable.Open(path);
        }
        return set;
    }

    /// <summary>
    /// Device and CPU cores work the same queue at once: the device takes batches of one phase
    /// (a batch shares one slowness field) from the front, each core one table from the back, so
    /// neither waits for the other and the faster of the two simply does more of the work.
    /// Whatever the device declines goes to the CPU.
    /// </summary>
    private static void Compute(ForwardModel model, List<(TableSource Src, Phase Phase, string Path, string Key)> jobs, EikonalOpenCl? gpu,
        int total, IProgress<(double, string)>? progress, Action<string>? log, CancellationToken ct)
    {
        jobs.Sort((x, y) => x.Phase.CompareTo(y.Phase));
        var gate = new object();
        int lo = 0, hi = jobs.Count, done = 0, onDevice = 0;
        var clock = System.Diagnostics.Stopwatch.StartNew();
        var batch = gpu != null ? Math.Max(1, gpu.SuggestBatch(model.Grid.Count)) : 0;
        var label = gpu != null ? "OpenCL + CPU" : "CPU";

        void Done(int count)
        {
            var d = Interlocked.Add(ref done, count);
            progress?.Report(((double)d / total, $"Travel-time tables {d}/{total} ({label})"));
        }

        void Cpu((TableSource Src, Phase Phase, string Path, string Key) job)
        {
            var s = model.Slowness(job.Phase);
            var seeds = SphericalEikonal.PointSource(model.Grid, s, job.Src.Lon, job.Src.Lat, job.Src.DepthKm);
            var times = SphericalEikonal.Solve(model.Metric, s, seeds, ct);
            WriteTable(job, model, times, "fast marching (CPU)");
            Done(1);
        }

        void CpuWorker()
        {
            while (true)
            {
                ct.ThrowIfCancellationRequested();
                (TableSource Src, Phase Phase, string Path, string Key) job;
                lock (gate)
                {
                    if (lo >= hi) return;
                    job = jobs[--hi];
                }
                Cpu(job);
            }
        }

        void DeviceWorker()
        {
            var first = true;
            while (true)
            {
                ct.ThrowIfCancellationRequested();
                List<(TableSource Src, Phase Phase, string Path, string Key)> items;
                lock (gate)
                {
                    if (lo >= hi) return;
                    var phase = jobs[lo].Phase;
                    var end = lo;
                    while (end < hi && end - lo < batch && jobs[end].Phase == phase) end++;
                    items = jobs.GetRange(lo, end - lo);
                    lo = end;
                }
                var s = model.Slowness(items[0].Phase);
                var seeds = items.Select(j => (IReadOnlyList<SeedNode>)SphericalEikonal.PointSource(model.Grid, s, j.Src.Lon, j.Src.Lat, j.Src.DepthKm)).ToList();
                var output = new double[(long)items.Count * model.Grid.Count];
                var t0 = clock.Elapsed;
                if (!gpu!.TrySolve(model.Metric, s, seeds, output, out var reason, ct))
                {
                    // From here on this thread is one more CPU worker.
                    log?.Invoke($"OpenCL declined the batch ({reason}); continuing on CPU.");
                    foreach (var job in items) Cpu(job);
                    CpuWorker();
                    return;
                }
                if (first)
                {
                    log?.Invoke($"OpenCL: {items.Count} tables in {(clock.Elapsed - t0).TotalSeconds:F1} s ({gpu.LastStats}); CPU cores take tables alongside.");
                    first = false;
                }
                for (var q = 0; q < items.Count; q++)
                    WriteTable(items[q], model, output.AsSpan(q * model.Grid.Count, model.Grid.Count), "OpenCL fast sweeping");
                Interlocked.Add(ref onDevice, items.Count);
                Done(items.Count);
            }
        }

        // One thread drives the device; it spends its time waiting on the queue, not computing.
        var cores = gpu != null ? Math.Max(1, ComputeSettings.Threads - 1) : ComputeSettings.Threads;
        var workers = new List<Task>();
        if (gpu != null) workers.Add(Task.Factory.StartNew(DeviceWorker, TaskCreationOptions.LongRunning));
        for (var w = 0; w < cores; w++) workers.Add(Task.Factory.StartNew(CpuWorker, TaskCreationOptions.LongRunning));
        try
        {
            Task.WaitAll(workers.ToArray());
        }
        catch (AggregateException ex)
        {
            if (ct.IsCancellationRequested) throw new OperationCanceledException(ct);
            System.Runtime.ExceptionServices.ExceptionDispatchInfo.Capture(ex.Flatten().InnerExceptions[0]).Throw();
        }
        if (gpu != null)
            log?.Invoke($"Travel-time tables: {onDevice} on OpenCL, {jobs.Count - onDevice} on CPU, in {clock.Elapsed.TotalSeconds:F0} s.");
    }

    private static void WriteTable((TableSource Src, Phase Phase, string Path, string Key) job, ForwardModel model, ReadOnlySpan<double> times, string solver)
    {
        TravelTimeTable.Write(job.Path, new TravelTimeTableHeader
        {
            StationId = job.Src.StationId, Phase = job.Phase, Lon = job.Src.Lon, Lat = job.Src.Lat, DepthKm = job.Src.DepthKm,
            Grid = model.Grid.Definition, ModelKey = job.Key, Solver = solver
        }, times);
    }

    private static string Safe(string s) => new(s.Select(c => char.IsLetterOrDigit(c) || c is '.' or '_' or '-' ? c : '_').ToArray());

    public void Dispose()
    {
        foreach (var t in _tables.Values) t.Dispose();
        _tables.Clear();
    }
}
