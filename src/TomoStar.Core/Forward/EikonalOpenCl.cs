using TomoStar.Core.Compute;
using TomoStar.Core.Geo;
using Silk.NET.OpenCL;

namespace TomoStar.Core.Forward;

/// <summary>
/// The spherical eikonal solved on an OpenCL device for a batch of sources at once.
///
/// A GPU cannot run a priority queue, so the device sweeps instead of marching: the upwind update
/// of <see cref="SphericalEikonal"/> is applied Gauss-Seidel fashion in the eight diagonal sweep
/// orderings of the fast sweeping method (Zhao 2005, Math. Comp. 74(250), 603-627), each ordering
/// taken plane by plane along i+j+k = const so every node of a plane is updated at once: none of
/// them is another's neighbour (Detrixhe, Gibou &amp; Min 2013, J. Comput. Phys. 237, 46-55).
/// Information crosses the whole grid in one sweep instead of one node per pass, so a field
/// settles in a handful of rounds of eight sweeps where a Jacobi relaxation needs as many passes
/// as the grid is wide, each over every node. First order to its fixed point, then with the
/// second-order promotion, exactly as the host <see cref="SphericalEikonal.Relax"/>, which walks
/// the same planes in the same order; the kernel is a term-by-term transcription of that update,
/// with floating-point contraction disabled so device and host round the same way.
///
/// A device is not trusted on arrival: <see cref="TryCreate"/> solves a heterogeneous field on the
/// device and on the host and declines the device if they disagree. A batch that has not settled
/// within the round cap is refused rather than returned half-relaxed.
///
/// A device without double precision (Apple GPUs, most ARM GPUs) holds each time as two floats,
/// hi + lo (about 48 bits, in the 8 bytes a double would take), and computes the update in float
/// as the step above the smallest upwind neighbour: every term of the quadratic is then the size of
/// a cell's crossing time, where in absolute times its discriminant is the difference of two numbers
/// the size of t² and float would keep none of the digits that survive the subtraction. Storing the
/// time itself as one float is not enough: the rounding of each node is inherited by the next, and on
/// a 200 × 180 × 110 grid 1000 km deep the fields drifted 0.1 s from the double solution. Paired,
/// the median difference is 1.5e-8 of the time and 99 % of the nodes are within the 6e-8 a float32
/// table rounds to anyway, so such a device meets the same <see cref="SelfTestTolerance"/>.
///
/// A node is updated only when a node it reads (one or two steps away on each axis) has moved
/// since its last update; otherwise the update would return the value it already holds, so the
/// field is the same to the bit and the settled bulk of the grid costs a flag read per sweep.
/// </summary>
public sealed unsafe class EikonalOpenCl : IDisposable
{
    public const double SelfTestTolerance = 1e-6;

    /// <summary>Host memory the output of one batch may take; the device is busy enough well below it.</summary>
    private const long BatchHostBytes = 768L << 20;

    /// <summary>
    /// Field nodes per batch past which a bigger batch no longer buys throughput. A sweep is tens of
    /// thousands of launches whatever the batch, so a small batch is paid in launch overhead: on an
    /// Apple M1 Max a 4-million-node grid costs 450 ms a table in batches of 4 and 156 ms in 32.
    /// </summary>
    private const long BatchNodes = 128L << 20;

    private const int LocalSize = 64;
    private readonly OpenClContext _ctx;
    private readonly nint _program;
    private readonly nint _sweep;
    private readonly nint _seed;
    private readonly nint _clear;
    private readonly nint _dirty;
    private readonly object _gate = new();

    /// <summary>Double precision on the device; otherwise times are held as two floats (hi + lo).</summary>
    private readonly bool _fp64;

    /// <summary>Bytes of one field value on the device: a double, or the two floats (hi + lo) of one.</summary>
    private int Real => sizeof(double);

    private EikonalOpenCl(OpenClContext ctx, nint program, nint sweep, nint seed, nint clear, nint dirty)
    {
        _ctx = ctx;
        _fp64 = ctx.Device.Fp64;
        _clear = clear;
        _dirty = dirty;
        _program = program;
        _sweep = sweep;
        _seed = seed;
    }

    public string Device => _ctx.Device.ToString();

    /// <summary>Rounds (of eight sweeps) the last batch took in each phase, for the log.</summary>
    public string LastStats { get; private set; } = "";

    public static EikonalOpenCl? TryCreate(Action<string>? log = null, bool allowCpuDevice = false)
    {
        var ctx = OpenClContext.TryCreate(log, requireFp64: false, allowCpuDevice: allowCpuDevice);
        if (ctx == null) return null;
        var program = ctx.BuildProgram(Kernel.Replace("{{TOL}}", SphericalEikonal.SettleTolerance.ToString("R", System.Globalization.CultureInfo.InvariantCulture)),
            log, ctx.Device.Fp64 ? "-D QV_FP64" : "");
        var sweep = program != 0 ? ctx.CreateKernel(program, "eik_sweep") : 0;
        var seed = program != 0 ? ctx.CreateKernel(program, "eik_seed") : 0;
        var clear = program != 0 ? ctx.CreateKernel(program, "eik_clear") : 0;
        var dirty = program != 0 ? ctx.CreateKernel(program, "eik_dirty") : 0;
        if (sweep == 0 || seed == 0 || clear == 0 || dirty == 0)
        {
            if (sweep != 0) ctx.Cl.ReleaseKernel(sweep);
            if (seed != 0) ctx.Cl.ReleaseKernel(seed);
            if (clear != 0) ctx.Cl.ReleaseKernel(clear);
            if (dirty != 0) ctx.Cl.ReleaseKernel(dirty);
            if (program != 0) ctx.Cl.ReleaseProgram(program);
            ctx.Dispose();
            return null;
        }
        var solver = new EikonalOpenCl(ctx, program, sweep, seed, clear, dirty);
        var err = solver.SelfTest(out var message);
        if (!(err <= SelfTestTolerance))
        {
            log?.Invoke($"OpenCL eikonal self-test failed ({message}); using CPU.");
            solver.Dispose();
            return null;
        }
        log?.Invoke($"OpenCL eikonal ready on {ctx.Device}{(solver._fp64 ? "" : " in paired single precision")}: self-test {message}.");
        return solver;
    }

    /// <summary>How many fields fit in one batch on this device (and in host memory) for a grid of this size.</summary>
    public int SuggestBatch(int nodes)
    {
        var dev = _ctx.Device;
        var perField = (long)nodes * Real + 2L * nodes;
        var budget = (long)(dev.GlobalMemory * 0.5) - 2L * nodes * Real;
        var byAlloc = (long)dev.MaxAlloc / Math.Max(1, (long)nodes * Real);
        var byHost = BatchHostBytes / Math.Max(1, (long)nodes * sizeof(double));
        // Enough fields to amortise the launches; the CPU cores work the same queue from the other
        // end, so a larger batch only moves where the two meet.
        var enough = Math.Max(4, BatchNodes / Math.Max(1, nodes));
        return (int)Math.Clamp(Math.Min(Math.Min(budget / Math.Max(1, perField), byAlloc), Math.Min(byHost, enough)), 0, 32);
    }

    /// <summary>
    /// Solves one field per seed set into <paramref name="output"/> (field-major, grid order).
    /// Returns false (and the caller falls back to the CPU) when the device declines.
    /// </summary>
    public bool TrySolve(GridMetric metric, double[] slowness, IReadOnlyList<IReadOnlyList<SeedNode>> seedSets, double[] output,
        out string reason, CancellationToken ct = default)
    {
        var g = metric.Grid;
        var n = g.Count;
        var batch = seedSets.Count;
        reason = "";
        if ((long)output.Length < (long)batch * n) { reason = "output too small"; return false; }
        // Seeds only: the fields start at +inf on the device, so the host never holds a second copy.
        var at = new List<long>();
        var time = new List<double>();
        var merged = new Dictionary<int, double>();
        for (var e = 0; e < batch; e++)
        {
            merged.Clear();
            foreach (var s in seedSets[e])
            {
                if ((uint)s.Index >= (uint)n) continue;
                if (!merged.TryGetValue(s.Index, out var t) || s.Time < t) merged[s.Index] = s.Time;
            }
            foreach (var (index, t) in merged)
            {
                at.Add((long)index * batch + e); // node-major: the batch's fields of a node sit side by side
                time.Add(t);
            }
        }
        if (at.Count == 0) { reason = "no seed nodes"; return false; }
        lock (_gate) return Sweep(metric, slowness, at.ToArray(), time.ToArray(), batch, output, ct, out reason);
    }

    private bool Sweep(GridMetric m, double[] slowness, long[] seedAt, double[] seedTime, int batch, double[] output,
        CancellationToken ct, out string reason)
    {
        var cl = _ctx.Cl;
        var g = m.Grid;
        var n = g.Count;
        var total = (long)batch * n;
        reason = "";
        nint bS = 0, bHl = 0, bHt = 0, bF = 0, bM = 0, bC = 0, bAt = 0, bT = 0, bD = 0;
        var changed = new int[1];
        int err;
        var real = Real;
        try
        {
            fixed (long* pa = seedAt)
            {
                bS = UploadReal(slowness, out err); if (err != 0) { reason = $"buffer {err}"; return false; }
                bHl = UploadReal(m.HLon, out err); if (err != 0) { reason = $"buffer {err}"; return false; }
                bHt = UploadReal(m.HLat, out err); if (err != 0) { reason = $"buffer {err}"; return false; }
                bAt = cl.CreateBuffer(_ctx.Context, MemFlags.ReadOnly | MemFlags.CopyHostPtr, (nuint)(seedAt.Length * 8L), pa, &err); if (err != 0) { reason = $"buffer {err}"; return false; }
                bT = UploadTime(seedTime, out err); if (err != 0) { reason = $"buffer {err}"; return false; }
            }
            bF = cl.CreateBuffer(_ctx.Context, MemFlags.ReadWrite, (nuint)(total * real), null, &err); if (err != 0) { reason = $"field buffer {err} (too large for the device?)"; return false; }
            bM = cl.CreateBuffer(_ctx.Context, MemFlags.ReadWrite, (nuint)total, null, &err); if (err != 0) { reason = $"mask buffer {err}"; return false; }
            bD = cl.CreateBuffer(_ctx.Context, MemFlags.ReadWrite, (nuint)total, null, &err); if (err != 0) { reason = $"flag buffer {err}"; return false; }
            bC = cl.CreateBuffer(_ctx.Context, MemFlags.ReadWrite, sizeof(int), null, &err); if (err != 0) { reason = $"flag buffer {err}"; return false; }

            // A kernel, not clEnqueueFillBuffer: Apple's OpenCL fills a buffer with garbage.
            SetArg(_clear, 0, bF); SetArg(_clear, 1, bM); SetArg(_clear, 2, bD);
            cl.SetKernelArg(_clear, 3, sizeof(long), &total);
            var clearGlobal = (nuint)((total + LocalSize - 1) / LocalSize * LocalSize);
            err = cl.EnqueueNdrangeKernel(_ctx.Queue, _clear, 1, null, &clearGlobal, null, 0, null, null);
            if (err != 0) { reason = $"clear launch {err}"; return false; }
            var count = seedAt.Length;
            SetArg(_seed, 0, bF); SetArg(_seed, 1, bM); SetArg(_seed, 2, bAt); SetArg(_seed, 3, bT);
            cl.SetKernelArg(_seed, 4, sizeof(int), &count);
            var seedGlobal = (nuint)count;
            err = cl.EnqueueNdrangeKernel(_ctx.Queue, _seed, 1, null, &seedGlobal, null, 0, null, null);
            if (err != 0) { reason = $"seed launch {err}"; return false; }

            var hr = m.HR;
            var hr32 = (float)hr;
            int nx = g.Nx, ny = g.Ny, nz = g.Nz;
            SetArg(_sweep, 0, bF); SetArg(_sweep, 1, bS); SetArg(_sweep, 2, bHl); SetArg(_sweep, 3, bHt);
            if (_fp64) cl.SetKernelArg(_sweep, 4, sizeof(double), &hr);
            else cl.SetKernelArg(_sweep, 4, sizeof(float), &hr32);
            SetArg(_sweep, 5, bM);
            cl.SetKernelArg(_sweep, 6, sizeof(int), &nx);
            cl.SetKernelArg(_sweep, 7, sizeof(int), &ny);
            cl.SetKernelArg(_sweep, 8, sizeof(int), &nz);
            SetArg(_sweep, 10, bC);
            var fields = batch;
            cl.SetKernelArg(_sweep, 13, sizeof(int), &fields);
            SetArg(_sweep, 14, bD);
            SetArg(_dirty, 0, bD);
            cl.SetKernelArg(_dirty, 1, sizeof(long), &total);

            // Once, not per launch: a stackalloc inside the plane loop is only freed when the method
            // returns, and a large grid launches tens of thousands of planes.
            var global = stackalloc nuint[2];
            var local = stackalloc nuint[2];
            var cap = SphericalEikonal.SweepRoundCap(g);
            var planes = nx + ny + nz - 2;
            var rounds = new int[2];
            for (var phase = 0; phase < 2; phase++)
            {
                cl.SetKernelArg(_sweep, 9, sizeof(int), &phase);
                if (phase > 0)
                {
                    err = cl.EnqueueNdrangeKernel(_ctx.Queue, _dirty, 1, null, &clearGlobal, null, 0, null, null);
                    if (err != 0) { reason = $"flag launch {err}"; return false; }
                }
                var settled = false;
                for (var round = 0; round < cap && !settled; round++)
                {
                    ct.ThrowIfCancellationRequested();
                    changed[0] = 0;
                    fixed (int* pc = changed) err = cl.EnqueueWriteBuffer(_ctx.Queue, bC, false, 0, sizeof(int), pc, 0, null, null);
                    if (err != 0) { reason = $"flag write {err}"; return false; }
                    for (var dir = 0; dir < 8; dir++)
                    {
                        cl.SetKernelArg(_sweep, 11, sizeof(int), &dir);
                        for (var p = 0; p < planes; p++)
                        {
                            SphericalEikonal.PlaneBounds(nx, ny, nz, p, out var i0, out var i1, out var k0, out var k1);
                            int iw = i1 - i0 + 1, kw = k1 - k0 + 1;
                            cl.SetKernelArg(_sweep, 12, sizeof(int), &p);
                            var items = (long)iw * batch;
                            global[0] = (nuint)((items + LocalSize - 1) / LocalSize * LocalSize); global[1] = (nuint)kw;
                            local[0] = LocalSize; local[1] = 1;
                            err = cl.EnqueueNdrangeKernel(_ctx.Queue, _sweep, 2, null, global, local, 0, null, null);
                            if (err != 0) { reason = $"kernel launch {err}"; return false; }
                        }
                        // Let a stop land within a sweep, not a round, on a large grid.
                        if (ct.IsCancellationRequested) { cl.Finish(_ctx.Queue); ct.ThrowIfCancellationRequested(); }
                    }
                    fixed (int* pc = changed) err = cl.EnqueueReadBuffer(_ctx.Queue, bC, true, 0, sizeof(int), pc, 0, null, null);
                    if (err != 0) { reason = $"flag read {err}"; return false; }
                    rounds[phase] = round + 1;
                    settled = changed[0] == 0;
                }
                if (!settled) { reason = $"sweeping not settled after {cap} rounds"; return false; }
            }
            LastStats = $"{rounds[0]}+{rounds[1]} rounds of 8 sweeps";
            // Back to field-major in slabs of nodes, so the host never holds the batch twice.
            var slab = (int)Math.Min(n, Math.Max(1, (32L << 20) / batch));
            var tmp = _fp64 ? new double[(long)slab * batch] : [];
            var pairs = _fp64 ? [] : new float[2L * slab * batch];
            fixed (double* pt = tmp)
            fixed (float* pp = pairs)
            for (var n0 = 0; n0 < n; n0 += slab)
            {
                var len = Math.Min(slab, n - n0);
                err = cl.EnqueueReadBuffer(_ctx.Queue, bF, true, (nuint)((long)n0 * batch * real), (nuint)((long)len * batch * real),
                    _fp64 ? pt : pp, 0, null, null);
                if (err != 0) { reason = $"read back {err}"; return false; }
                Parallel.For(0, batch, e =>
                {
                    var dst = (long)e * n + n0;
                    if (_fp64) for (var q = 0; q < len; q++) output[dst + q] = tmp[(long)q * batch + e];
                    else
                        for (var q = 0; q < len; q++)
                        {
                            var at = 2 * ((long)q * batch + e);
                            output[dst + q] = (double)pairs[at] + pairs[at + 1];
                        }
                });
            }
            return true;
        }
        finally
        {
            cl.Finish(_ctx.Queue);
            foreach (var b in new[] { bS, bHl, bHt, bF, bM, bC, bAt, bT, bD })
                if (b != 0) cl.ReleaseMemObject(b);
        }

        void SetArg(nint kernel, uint index, nint buffer) => cl.SetKernelArg(kernel, index, (nuint)sizeof(nint), &buffer);

        // Times as the device holds them: doubles, or each as the two floats hi + lo.
        nint UploadTime(double[] data, out int e)
        {
            if (_fp64) return UploadReal(data, out e);
            var pairs = new float[2 * data.Length];
            for (var q = 0; q < data.Length; q++)
            {
                var hi = (float)data[q];
                pairs[2 * q] = hi;
                pairs[2 * q + 1] = float.IsFinite(hi) ? (float)(data[q] - hi) : 0f;
            }
            int code;
            nint b;
            fixed (float* p = pairs) b = cl.CreateBuffer(_ctx.Context, MemFlags.ReadOnly | MemFlags.CopyHostPtr, (nuint)(pairs.Length * 4L), p, &code);
            e = code;
            return b;
        }

        // A read-only copy of host doubles in the device's precision.
        nint UploadReal(double[] data, out int e)
        {
            int code;
            nint b;
            if (_fp64)
            {
                fixed (double* p = data) b = cl.CreateBuffer(_ctx.Context, MemFlags.ReadOnly | MemFlags.CopyHostPtr, (nuint)(data.Length * 8L), p, &code);
            }
            else
            {
                var narrow = new float[data.Length];
                for (var q = 0; q < data.Length; q++) narrow[q] = (float)data[q];
                fixed (float* p = narrow) b = cl.CreateBuffer(_ctx.Context, MemFlags.ReadOnly | MemFlags.CopyHostPtr, (nuint)(narrow.Length * 4L), p, &code);
            }
            e = code;
            return b;
        }
    }

    /// <summary>Largest device-host difference relative to the field's range on a heterogeneous test.</summary>
    public double SelfTest(out string message)
    {
        var grid = new SphericalGrid(new GridDefinition
        {
            MinLon = 10, MaxLon = 12, MinLat = 42, MaxLat = 44, MinDepthKm = -2, MaxDepthKm = 120, Nx = 22, Ny = 20, Nz = 16
        });
        var metric = new GridMetric(grid);
        var s = new double[grid.Count];
        for (var k = 0; k < grid.Nz; k++)
        for (var j = 0; j < grid.Ny; j++)
        for (var i = 0; i < grid.Nx; i++)
            s[grid.Index(i, j, k)] = 1 / (5.8 + 0.02 * grid.DepthKm[k] + 0.2 * Math.Sin(0.7 * i) + 0.15 * Math.Cos(0.5 * j));
        var seeds = SphericalEikonal.PointSource(grid, s, 10.37, 42.61, 8.3);
        var host = SphericalEikonal.Relax(metric, s, seeds);
        var dev = new double[grid.Count];
        if (!TrySolve(metric, s, [seeds], dev, out var reason)) { message = reason; return double.PositiveInfinity; }
        double lo = double.PositiveInfinity, hi = double.NegativeInfinity, worst = 0;
        for (var i = 0; i < host.Length; i++)
        {
            if (!double.IsFinite(host[i])) continue;
            if (!double.IsFinite(dev[i])) { message = "device left a node unreached"; return double.PositiveInfinity; }
            lo = Math.Min(lo, host[i]);
            hi = Math.Max(hi, host[i]);
            worst = Math.Max(worst, Math.Abs(dev[i] - host[i]));
        }
        var rel = worst / Math.Max(1e-12, hi - lo);
        message = $"max |Δt| {worst:E2} s over a {hi - lo:F2} s field ({rel:E2} relative)";
        return rel;
    }

    public void Dispose()
    {
        _ctx.Cl.ReleaseKernel(_sweep);
        _ctx.Cl.ReleaseKernel(_seed);
        _ctx.Cl.ReleaseKernel(_clear);
        _ctx.Cl.ReleaseKernel(_dirty);
        _ctx.Cl.ReleaseProgram(_program);
        _ctx.Dispose();
    }

    // ASCII only: some OpenCL compilers truncate or reject non-ASCII source.
    private const string Kernel = """
#ifdef QV_FP64
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
typedef double real;
typedef double tval;
#define TVAL_INF INFINITY
#define TOL {{TOL}}
#else
typedef float real;
typedef float2 tval;   // a time as hi + lo, two floats: about 48 bits, in the 8 bytes of a double
#define TVAL_INF ((float2)(INFINITY, 0.0f))
#define TOL {{TOL}}f
#endif
#pragma OPENCL FP_CONTRACT OFF

// Constants are float literals: every one is exact in float, so the double kernel rounds exactly as
// the host does and the float kernel never meets a double.

__kernel void eik_clear(__global tval* field, __global uchar* seeded, __global uchar* dirty, long total)
{
    long q = get_global_id(0);
    if (q >= total) return;
    field[q] = TVAL_INF;
    seeded[q] = 0;
    dirty[q] = 1;
}

// Every node due for an update: the start of a phase, whose update is a different function.
__kernel void eik_dirty(__global uchar* dirty, long total)
{
    long q = get_global_id(0);
    if (q < total) dirty[q] = 1;
}

__kernel void eik_seed(__global tval* field, __global uchar* seeded, __global const long* at, __global const tval* t, int count)
{
    int q = get_global_id(0);
    if (q >= count) return;
    field[at[q]] = t[q];
    seeded[at[q]] = 1;
}

#ifdef QV_FP64

// f points at this field's value of node 0; node x of the field is f[x * batch].
static void eik_axis(__global const real* f, long batch, int pos, int n, int idx, int stride, real h, int second_order,
                     real* te, real* we, real* tn, int* used)
{
    real minus = pos > 0 ? f[(idx - stride) * batch] : INFINITY;
    real plus = pos + 1 < n ? f[(idx + stride) * batch] : INFINITY;
    real first;
    int second;
    if (minus <= plus) { first = minus; second = pos > 1 ? idx - 2 * stride : -1; }
    else { first = plus; second = pos + 2 < n ? idx + 2 * stride : -1; }
    if (!isfinite(first)) return;
    real t2 = (second_order && second >= 0) ? f[second * batch] : INFINITY;
    int u = *used;
    real he;
    if (isfinite(t2) && t2 <= first) { te[u] = (4.0f * first - t2) / 3.0f; he = 2.0f * h / 3.0f; }
    else { te[u] = first; he = h; }
    we[u] = 1.0f / (he * he);
    tn[u] = first;
    *used = u + 1;
}

#else

// hi + lo of two times as one float: exact in the hi part when the two are within a factor of two
// (Sterbenz), which neighbouring times are; +inf when a is.
static float df_sub(float2 a, float2 b) { return (a.x - b.x) + (a.y - b.y); }

// r + x as hi + lo (Knuth's two-sum of r.x and r.y + x).
static float2 df_add(float2 r, float x)
{
    float b = r.y + x;
    float s = r.x + b;
    float bb = s - r.x;
    float e = (r.x - (s - bb)) + (b - bb);
    return (float2)(s, e);
}

static int df_less(float2 a, float2 b) { return a.x < b.x || (a.x == b.x && a.y < b.y); }

// The upwind neighbour of one axis, as a two-float time; second is its second-order partner.
static int eik_pick(__global const float2* f, long batch, int pos, int n, int idx, int stride, int second_order,
                    float2* first, float2* t2)
{
    float2 minus = pos > 0 ? f[(idx - stride) * batch] : TVAL_INF;
    float2 plus = pos + 1 < n ? f[(idx + stride) * batch] : TVAL_INF;
    int second;
    if (!df_less(plus, minus)) { *first = minus; second = pos > 1 ? idx - 2 * stride : -1; }
    else { *first = plus; second = pos + 2 < n ? idx + 2 * stride : -1; }
    if (!isfinite((*first).x)) return 0;
    *t2 = (second_order && second >= 0) ? f[second * batch] : TVAL_INF;
    return 1;
}

// One axis's term as offsets (float) from the reference time r, which is the smallest upwind
// neighbour: every quantity the update needs is then the size of a cell's crossing time.
static void eik_term(float2 first, float2 t2, float2 r, float h, float* te, float* we, float* tn, int* used)
{
    int u = *used;
    float of = df_sub(first, r);
    float he;
    if (isfinite(t2.x) && !df_less(first, t2)) { te[u] = of + df_sub(first, t2) / 3.0f; he = 2.0f * h / 3.0f; }
    else { te[u] = of; he = h; }
    we[u] = 1.0f / (he * he);
    tn[u] = of;
    *used = u + 1;
}

#endif

// One plane i'+j'+k' = p of one sweep ordering; (i', j', k') are (i, j, k) mirrored on the axes
// whose bit is set in dir. Work item (a * batch + e, c) is i' = i0 + a, k' = k0 + c of field e,
// and j' follows; items off the plane leave at once. No node of a plane is a neighbour of
// another, so the in-place update is race free and its result does not depend on the order of
// the items. The fields are stored node-major (node x of field e at x * batch + e), so the items
// of one node read and write adjacent words: a diagonal plane has no two nodes adjacent in memory.
__kernel void eik_sweep(__global tval* fields, __global const real* slowness, __global const real* hlon, __global const real* hlat,
                        real hr, __global const uchar* seeded, int nx, int ny, int nz, int second_order, __global int* changed,
                        int dir, int p, int batch, __global uchar* dirty)
{
    // The plane's i' and k' ranges, as SphericalEikonal.PlaneBounds gives them to the host.
    int i0 = max(0, p - (ny - 1) - (nz - 1));
    int iw = min(nx - 1, p) - i0 + 1;
    int k0 = max(0, p - (nx - 1) - (ny - 1));
    int t = get_global_id(0);
    int a = t / batch;
    if (a >= iw) return;
    int e = t - a * batch;
    int ip = i0 + a;
    int kp = k0 + (int)get_global_id(1);
    int jp = p - ip - kp;
    if (jp < 0 || jp >= ny) return;
    int i = (dir & 1) ? nx - 1 - ip : ip;
    int j = (dir & 2) ? ny - 1 - jp : jp;
    int k = (dir & 4) ? nz - 1 - kp : kp;
    int nxy = nx * ny;
    int node = (k * ny + j) * nx + i;
    long bl = batch;
    if (seeded[node * bl + e]) return;
    // A node none of whose inputs (the neighbours one and two nodes away on each axis) moved since
    // it was last updated would get the value it already has: skip it. The field is the same.
    if (!dirty[node * bl + e]) return;
    dirty[node * bl + e] = 0;
    real s = slowness[node];
    if (!(s > 0.0f)) return;
    __global tval* f = fields + e;
    tval cur = f[node * bl];

    real te[3], we[3], tn[3];
    int used = 0;
#ifdef QV_FP64
    eik_axis(f, bl, i, nx, node, 1, hlon[k * ny + j], second_order, te, we, tn, &used);
    eik_axis(f, bl, j, ny, node, nx, hlat[k], second_order, te, we, tn, &used);
    eik_axis(f, bl, k, nz, node, nxy, hr, second_order, te, we, tn, &used);
#else
    float2 f0, f1, f2, s0, s1, s2;
    int u0 = eik_pick(f, bl, i, nx, node, 1, second_order, &f0, &s0);
    int u1 = eik_pick(f, bl, j, ny, node, nx, second_order, &f1, &s1);
    int u2 = eik_pick(f, bl, k, nz, node, nxy, second_order, &f2, &s2);
    if (!(u0 | u1 | u2)) return;
    float2 r = TVAL_INF;
    if (u0 && df_less(f0, r)) r = f0;
    if (u1 && df_less(f1, r)) r = f1;
    if (u2 && df_less(f2, r)) r = f2;
    if (u0) eik_term(f0, s0, r, hlon[k * ny + j], te, we, tn, &used);
    if (u1) eik_term(f1, s1, r, hlat[k], te, we, tn, &used);
    if (u2) eik_term(f2, s2, r, hr, te, we, tn, &used);
#endif
    if (used == 0) return;

    for (int q = 1; q < used; q++) {
        real tt = te[q], ww = we[q], nn = tn[q];
        int b = q - 1;
        while (b >= 0 && te[b] > tt) { te[b + 1] = te[b]; we[b + 1] = we[b]; tn[b + 1] = tn[b]; b--; }
        te[b + 1] = tt; we[b + 1] = ww; tn[b + 1] = nn;
    }

    real cand = INFINITY;
#ifndef QV_FP64
    // sum w (t - te)^2 = s^2 solved for x = t - te[0]: the terms are cell crossing times, not t.
    real base = te[0];
#endif
    for (int terms = used; terms >= 1; terms--) {
#ifdef QV_FP64
        real sw = 0.0f, swt = 0.0f, swt2 = 0.0f;
        for (int q = 0; q < terms; q++) {
            real w = we[q];
            sw += w; swt += w * te[q]; swt2 += w * te[q] * te[q];
        }
        real disc = swt * swt - sw * (swt2 - s * s);
        if (disc < 0.0f) continue;
        real root = (swt + sqrt(disc)) / sw;
#else
        real sw = 0.0f, swd = 0.0f, swd2 = 0.0f;
        for (int q = 0; q < terms; q++) {
            real w = we[q], d = te[q] - base;
            sw += w; swd += w * d; swd2 += w * d * d;
        }
        real disc = swd * swd - sw * (swd2 - s * s);
        if (disc < 0.0f) continue;
        real root = base + (swd + sqrt(disc)) / sw;
#endif
        int ok = 1;
        for (int q = 0; q < terms; q++) if (root < tn[q]) { ok = 0; break; }
        if (ok) { cand = root; break; }
    }
#ifdef QV_FP64
    if (!(cand < cur)) return;
    f[node * bl] = cand;
    if (!isfinite(cur) || cur - cand > TOL) atomic_or(changed, 1);
#else
    if (!isfinite(cand)) return;
    float2 next = df_add(r, cand);
    float drop = df_sub(cur, next);   // +inf when cur is
    if (!(drop > 0.0f)) return;
    f[node * bl] = next;
    if (drop > TOL) atomic_or(changed, 1);
#endif
    // Wake the nodes that read this one. They lie on other planes of the sweep, so no item of
    // this launch reads or clears these flags; two items may set the same one, to the same value.
    __global uchar* d = dirty + e;
    if (i > 0) d[(node - 1) * bl] = 1;
    if (i > 1) d[(node - 2) * bl] = 1;
    if (i + 1 < nx) d[(node + 1) * bl] = 1;
    if (i + 2 < nx) d[(node + 2) * bl] = 1;
    if (j > 0) d[(node - nx) * bl] = 1;
    if (j > 1) d[(node - 2 * nx) * bl] = 1;
    if (j + 1 < ny) d[(node + nx) * bl] = 1;
    if (j + 2 < ny) d[(node + 2 * nx) * bl] = 1;
    if (k > 0) d[(node - nxy) * bl] = 1;
    if (k > 1) d[(node - 2 * nxy) * bl] = 1;
    if (k + 1 < nz) d[(node + nxy) * bl] = 1;
    if (k + 2 < nz) d[(node + 2 * nxy) * bl] = 1;
}
""";
}
