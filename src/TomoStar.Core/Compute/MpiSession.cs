// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0
using System.Runtime.InteropServices;
using System.Text.Json;
using System.Text.Json.Serialization;
using TomoStar.Core.Forward;
using TomoStar.Core.Geo;
using TomoStar.Core.Model;
using TomoStar.Core.Tomography;

namespace TomoStar.Core.Compute;

/// <summary>Optional MPI forward workers. Only rank zero runs the CLI and the inversion.</summary>
public sealed class MpiSession : IDisposable
{
    private static readonly JsonSerializerOptions Json = new()
    {
        PreferredObjectCreationHandling = JsonObjectCreationHandling.Populate,
        NumberHandling = JsonNumberHandling.AllowNamedFloatingPointLiterals
    };
    public static MpiSession? Current { get; private set; }
    private bool _local;
    private bool _disposed;
    public int Rank { get; }
    public int Size { get; }
    public bool CanDispatch => Rank == 0 && Size > 1 && !_local;

    public MpiSession()
    {
        if (Current != null) throw new InvalidOperationException("MPI is already initialized.");
        try { Check(Native.Init()); }
        catch (DllNotFoundException ex) { throw new InvalidOperationException("MPI requires libtomostar_mpi.so; build native/mpi/build.sh against the cluster MPI and set LD_LIBRARY_PATH.", ex); }
        Check(Native.Rank(out var rank));
        Check(Native.Size(out var size));
        Rank = rank; Size = size;
        Current = this;
    }

    internal sealed class Request
    {
        public string Kind { get; set; } = "";
        public GridDefinition Grid { get; set; } = new();
        public ForwardDomain? Domain { get; set; }
        public ObservationSet? Data { get; set; }
        public double[] P { get; set; } = [];
        public double[] S { get; set; } = [];
        public List<TableSource> Sources { get; set; } = [];
        public List<Phase> Phases { get; set; } = [];
        public string Folder { get; set; } = "";
        public RayMethod Method { get; set; }
        public int KeepPoints { get; set; }
        public int Threads { get; set; }
        public bool UseGpu { get; set; }
    }
    private sealed class Response
    {
        public string? Error { get; set; }
        public List<RayRow> Rows { get; set; } = [];
    }

    public void WorkerLoop()
    {
        if (Rank == 0) throw new InvalidOperationException("Rank zero cannot enter the worker loop.");
        try
        {
            while (true)
            {
                var payload = Broadcast(null, 0);
                if (payload.Length == 0) return;
                var request = JsonSerializer.Deserialize<Request>(payload, Json)!;
                var result = RunLocal(request);
                Collect(result);
            }
        }
        catch { Native.Abort(1); throw; }
    }

    internal List<RayRow> Dispatch(Request request)
    {
        if (!CanDispatch) throw new InvalidOperationException("Only rank zero may dispatch MPI work.");
        request.Threads = ComputeSettings.MaxThreads;
        Broadcast(JsonSerializer.SerializeToUtf8Bytes(request, Json), 0);
        var responses = Collect(RunLocal(request));
        var failures = responses.Select((r, i) => r.Error == null ? null : $"rank {i}: {r.Error}").Where(e => e != null).ToList();
        if (failures.Count != 0) throw new InvalidOperationException("MPI forward failed: " + string.Join("; ", failures));
        return responses.SelectMany(r => r.Rows).OrderBy(r => r.Arrival).ToList();
    }

    private Response RunLocal(Request request)
    {
        _local = true;
        try
        {
            ComputeSettings.MaxThreads = request.Threads;
            using var gpu = request.UseGpu ? EikonalOpenCl.TryCreate() : null;
            if (request.Kind == "tables")
            {
                var grid = new SphericalGrid(request.Grid);
                var model = new ForwardModel { Grid = grid, Metric = new GridMetric(grid), SlownessP = request.P, SlownessS = request.S };
                using var tables = TravelTimeTableSet.Build(model,
                    request.Sources.Where((_, i) => i % Size == Rank).ToList(), request.Phases, request.Folder, gpu);
                return new Response();
            }
            if (request.Kind != "rays") throw new InvalidOperationException("Unknown MPI request.");
            var full = request.Data!;
            var indices = Enumerable.Range(0, full.Arrivals.Count).Where(i => full.Arrivals[i].Station % Size == Rank).ToArray();
            if (indices.Length == 0) return new Response();
            var data = new ObservationSet();
            data.Events.AddRange(full.Events);
            data.Stations.AddRange(full.Stations);
            data.Arrivals.AddRange(indices.Select(i => full.Arrivals[i]));
            var rows = new RayKernels(new SphericalGrid(request.Grid), request.Domain!, request.Method,
                Path.Combine(request.Folder, $"mpi-rank-{Rank}")).Compute(data, request.P, request.S, gpu, request.KeepPoints);
            return new Response { Rows = rows.Select(r => r with
            {
                Arrival = indices[r.Arrival],
                Ray = r.Ray == null ? null : r.Ray with { Arrival = indices[r.Arrival] }
            }).ToList() };
        }
        catch (Exception ex) { return new Response { Error = ex.Message }; }
        finally { _local = false; }
    }

    private List<Response> Collect(Response local)
    {
        var all = new List<Response>();
        for (var rank = 0; rank < Size; rank++)
        {
            var payload = Broadcast(Rank == rank ? JsonSerializer.SerializeToUtf8Bytes(local, Json) : null, rank);
            if (Rank == 0) all.Add(JsonSerializer.Deserialize<Response>(payload, Json)!);
        }
        return all;
    }

    private byte[] Broadcast(byte[]? payload, int root)
    {
        long length = Rank == root ? payload!.LongLength : 0;
        Check(Native.Length(ref length, root));
        if (length < 0 || length > Array.MaxLength) { Native.Abort(1); throw new InvalidOperationException("MPI message exceeds the supported array size."); }
        var bytes = Rank == root ? payload! : new byte[(int)length];
        Check(Native.Bytes(bytes, bytes.Length, root));
        return bytes;
    }
    private static void Check(int code)
    {
        if (code == 0) return;
        Native.Abort(code);
        throw new InvalidOperationException($"MPI failed with code {code}.");
    }
    public void Dispose()
    {
        if (_disposed) return;
        _disposed = true;
        if (Rank == 0 && Size > 1) Broadcast([], 0);
        Check(Native.FinalizeMpi());
        Current = null;
    }
    private static class Native
    {
        private const string Library = "tomostar_mpi";
        [DllImport(Library, EntryPoint = "ts_mpi_init")] internal static extern int Init();
        [DllImport(Library, EntryPoint = "ts_mpi_rank")] internal static extern int Rank(out int rank);
        [DllImport(Library, EntryPoint = "ts_mpi_size")] internal static extern int Size(out int size);
        [DllImport(Library, EntryPoint = "ts_mpi_bcast_length")] internal static extern int Length(ref long length, int root);
        [DllImport(Library, EntryPoint = "ts_mpi_bcast_bytes")] internal static extern int Bytes([In, Out] byte[] bytes, int count, int root);
        [DllImport(Library, EntryPoint = "ts_mpi_abort")] internal static extern int Abort(int code);
        [DllImport(Library, EntryPoint = "ts_mpi_finalize")] internal static extern int FinalizeMpi();
    }
}
