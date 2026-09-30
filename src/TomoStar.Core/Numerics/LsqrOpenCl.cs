// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using TomoStar.Core.Compute;
using Silk.NET.OpenCL;

namespace TomoStar.Core.Numerics;

/// <summary>
/// <see cref="Lsqr"/> run entirely on an OpenCL device: the matrix, its transpose and every LSQR vector
/// live in device memory, and only the three scalars of each iteration (‖u‖, ‖v‖, ‖w‖²) come back to
/// the host, which carries the plane rotations and the stopping rules exactly as <see cref="Lsqr.Solve"/>
/// does. The loop is a line-by-line transcription of the host solver.
///
/// Kernels: a row-per-work-item CSR product (applied to A and to the stored Aᵀ, so the transpose
/// product is a gather with no atomics), the vector updates of <see cref="SimdVector"/>, and a dot
/// product reduced in a fixed number of work-groups whose partial sums the host adds in order, so a
/// solve is repeatable to the bit on a given device. Results agree with the host solver to rounding,
/// not to the bit (the host reduces in SIMD lanes).
///
/// A device without double precision (Apple GPUs, most ARM GPUs) runs the same kernels on pairs of
/// floats, hi + lo, whose sum carries about 48 bits (<see cref="PairedSource"/>): the matrix, the
/// vectors and the dot products are paired, the host keeps the recurrences in double, and the
/// solution agrees with the host solver as a double device's does.
///
/// As every device path in TomoSTAR, a device is used only after <see cref="TryCreate"/> has solved a
/// test system on it and on the host and found them in agreement; any failure during a solve returns
/// null and the caller continues on the CPU (<see cref="SolveOrCpu"/>).
/// </summary>
public sealed unsafe class LsqrOpenCl : IDisposable
{
    /// <summary>Largest relative difference of the self-test solution from the host solver's.</summary>
    public const double SelfTestTolerance = 1e-8;

    /// <summary>Systems with fewer non-zeros than this stay on the CPU, where they are faster than the transfers.</summary>
    public static long MinNonZerosForDevice { get; set; } = 2_000_000;

    private const int Groups = 64;
    private const int GroupSize = 256;

    private const string Source = """
        #pragma OPENCL EXTENSION cl_khr_fp64 : enable
        #pragma OPENCL FP_CONTRACT OFF

        __kernel void spmv(const int rows, __global const int* rs, __global const int* ci, __global const double* va,
                           __global const double* x, __global double* y)
        {
            int i = get_global_id(0);
            if (i >= rows) return;
            double s = 0.0;
            for (int k = rs[i]; k < rs[i + 1]; k++) s += va[k] * x[ci[k]];
            y[i] = s;
        }

        __kernel void xmby(const int n, __global const double* x, const double b, __global double* y)
        {
            int i = get_global_id(0);
            if (i < n) y[i] = x[i] - b * y[i];
        }

        __kernel void xpby(const int n, __global const double* x, const double b, __global double* y)
        {
            int i = get_global_id(0);
            if (i < n) y[i] = x[i] + b * y[i];
        }

        __kernel void axpy(const int n, const double a, __global const double* x, __global double* y)
        {
            int i = get_global_id(0);
            if (i < n) y[i] += a * x[i];
        }

        __kernel void scal(const int n, __global double* x, const double s)
        {
            int i = get_global_id(0);
            if (i < n) x[i] *= s;
        }

        __kernel void dot_partial(const int n, __global const double* a, __global const double* b,
                                  __global double* partial, __local double* scratch)
        {
            int lid = get_local_id(0);
            int size = get_local_size(0);
            int stride = get_num_groups(0) * size;
            double s = 0.0;
            for (int i = get_group_id(0) * size + lid; i < n; i += stride) s += a[i] * b[i];
            scratch[lid] = s;
            barrier(CLK_LOCAL_MEM_FENCE);
            for (int h = size / 2; h > 0; h >>= 1)
            {
                if (lid < h) scratch[lid] += scratch[lid + h];
                barrier(CLK_LOCAL_MEM_FENCE);
            }
            if (lid == 0) partial[get_group_id(0)] = scratch[0];
        }
        """;

    /// <summary>
    /// <see cref="Source"/> on pairs of floats (hi + lo): Knuth's two-sum and an fma product give
    /// sums and products to about 48 bits, where one float would keep 24. Every buffer and scalar
    /// takes the 8 bytes a double would.
    /// </summary>
    private const string PairedSource = """
        #pragma OPENCL FP_CONTRACT OFF

        static float2 df_add(float2 a, float2 b)
        {
            float s = a.x + b.x;
            float bb = s - a.x;
            float e = (a.x - (s - bb)) + (b.x - bb);
            e += a.y + b.y;
            float hi = s + e;
            return (float2)(hi, e - (hi - s));
        }

        static float2 df_mul(float2 a, float2 b)
        {
            float p = a.x * b.x;
            float e = fma(a.x, b.x, -p);
            e += a.x * b.y + a.y * b.x;
            float hi = p + e;
            return (float2)(hi, e - (hi - p));
        }

        __kernel void spmv(const int rows, __global const int* rs, __global const int* ci, __global const float2* va,
                           __global const float2* x, __global float2* y)
        {
            int i = get_global_id(0);
            if (i >= rows) return;
            float2 s = (float2)(0.0f, 0.0f);
            for (int k = rs[i]; k < rs[i + 1]; k++) s = df_add(s, df_mul(va[k], x[ci[k]]));
            y[i] = s;
        }

        __kernel void xmby(const int n, __global const float2* x, const float2 b, __global float2* y)
        {
            int i = get_global_id(0);
            if (i < n) y[i] = df_add(x[i], df_mul(-b, y[i]));
        }

        __kernel void xpby(const int n, __global const float2* x, const float2 b, __global float2* y)
        {
            int i = get_global_id(0);
            if (i < n) y[i] = df_add(x[i], df_mul(b, y[i]));
        }

        __kernel void axpy(const int n, const float2 a, __global const float2* x, __global float2* y)
        {
            int i = get_global_id(0);
            if (i < n) y[i] = df_add(y[i], df_mul(a, x[i]));
        }

        __kernel void scal(const int n, __global float2* x, const float2 s)
        {
            int i = get_global_id(0);
            if (i < n) x[i] = df_mul(x[i], s);
        }

        __kernel void dot_partial(const int n, __global const float2* a, __global const float2* b,
                                  __global float2* partial, __local float2* scratch)
        {
            int lid = get_local_id(0);
            int size = get_local_size(0);
            int stride = get_num_groups(0) * size;
            float2 s = (float2)(0.0f, 0.0f);
            for (int i = get_group_id(0) * size + lid; i < n; i += stride) s = df_add(s, df_mul(a[i], b[i]));
            scratch[lid] = s;
            barrier(CLK_LOCAL_MEM_FENCE);
            for (int h = size / 2; h > 0; h >>= 1)
            {
                if (lid < h) scratch[lid] = df_add(scratch[lid], scratch[lid + h]);
                barrier(CLK_LOCAL_MEM_FENCE);
            }
            if (lid == 0) partial[get_group_id(0)] = scratch[0];
        }
        """;

    private readonly OpenClContext _ctx;
    private readonly nint _program;

    /// <summary>Double precision on the device; otherwise every value is a pair of floats (hi + lo).</summary>
    private readonly bool _fp64;
    private readonly nint _spmv, _xmby, _xpby, _axpy, _scal, _dot;
    private readonly object _gate = new();

    private LsqrOpenCl(OpenClContext ctx, nint program)
    {
        _ctx = ctx;
        _fp64 = ctx.Device.Fp64;
        _program = program;
        _spmv = ctx.CreateKernel(program, "spmv");
        _xmby = ctx.CreateKernel(program, "xmby");
        _xpby = ctx.CreateKernel(program, "xpby");
        _axpy = ctx.CreateKernel(program, "axpy");
        _scal = ctx.CreateKernel(program, "scal");
        _dot = ctx.CreateKernel(program, "dot_partial");
    }

    private bool Complete => _spmv != 0 && _xmby != 0 && _xpby != 0 && _axpy != 0 && _scal != 0 && _dot != 0;

    public string Device => _ctx.Device.ToString() + (_fp64 ? "" : " in paired single precision");

    /// <summary>
    /// The solver on the best double-precision device, or null (no device, build failure, or a self-test
    /// that disagrees with the host). CPU-type devices (POCL) are left out unless allowed: on them the
    /// native multithreaded solver is faster.
    /// </summary>
    public static LsqrOpenCl? TryCreate(Action<string>? log = null, bool allowCpuDevice = false)
    {
        var ctx = OpenClContext.TryCreate(log, requireFp64: false, allowCpuDevice: allowCpuDevice);
        if (ctx == null) return null;
        var program = ctx.BuildProgram(ctx.Device.Fp64 ? Source : PairedSource, log);
        if (program == 0) { ctx.Dispose(); return null; }
        var solver = new LsqrOpenCl(ctx, program);
        if (!solver.Complete)
        {
            log?.Invoke("OpenCL LSQR: kernels could not be created; using CPU.");
            solver.Dispose();
            return null;
        }
        var err = solver.SelfTest(out var message);
        if (!(err <= SelfTestTolerance))
        {
            log?.Invoke($"OpenCL LSQR self-test failed ({message}); using CPU.");
            solver.Dispose();
            return null;
        }
        log?.Invoke($"OpenCL LSQR ready on {solver.Device}: self-test {message}.");
        return solver;
    }

    /// <summary>
    /// Solves on the device when there is one and the system is large enough, otherwise (or when the
    /// device fails) on the CPU. <paramref name="where"/> says which ran.
    /// </summary>
    public static LsqrResult SolveOrCpu(LsqrOpenCl? device, CsrMatrix a, double[] b, double damp, int maxIterations, double atol, double btol,
        double conlim, CancellationToken ct, Action<int, double>? onIteration, out string where, Action<string>? log = null)
    {
        if (device != null && a.NonZeroCount >= MinNonZerosForDevice)
        {
            var r = device.TrySolve(a, b, damp, maxIterations, atol, btol, conlim, ct, onIteration, out var reason);
            if (r != null) { where = $"OpenCL ({device.Device})"; return r; }
            log?.Invoke($"OpenCL LSQR declined ({reason}); solving on the CPU.");
        }
        where = "CPU";
        return Lsqr.Solve(a, b, damp, maxIterations, atol, btol, conlim, ct, onIteration);
    }

    /// <summary>The LSQR solve on the device; null (with the reason) when the device cannot take it.</summary>
    public LsqrResult? TrySolve(CsrMatrix a, double[] b, double damp = 0, int maxIterations = 500, double atol = 1e-8, double btol = 1e-8,
        double conlim = 1e8, CancellationToken ct = default, Action<int, double>? onIteration = null) =>
        TrySolve(a, b, damp, maxIterations, atol, btol, conlim, ct, onIteration, out _);

    public LsqrResult? TrySolve(CsrMatrix a, double[] b, double damp, int maxIterations, double atol, double btol,
        double conlim, CancellationToken ct, Action<int, double>? onIteration, out string reason)
    {
        lock (_gate)
        {
            try
            {
                return Run(a, b, damp, maxIterations, atol, btol, conlim, ct, onIteration, out reason);
            }
            catch (OperationCanceledException) { throw; }
            catch (Exception ex)
            {
                reason = ex.Message;
                return null;
            }
        }
    }

    private LsqrResult? Run(CsrMatrix a, double[] b, double damp, int maxIterations, double atol, double btol,
        double conlim, CancellationToken ct, Action<int, double>? onIteration, out string reason)
    {
        reason = "";
        var m = a.RowCount;
        var n = a.ColumnCount;
        if (b.Length != m) throw new ArgumentException("Right-hand side length differs from the row count.");
        var x = new double[n];
        if (m == 0 || n == 0) return new LsqrResult(x, 0, 0, 0, 0, "empty system");
        var at = a.Transpose();
        var bytes = 8L * (2L * a.NonZeroCount + 3L * m + 4L * n) + 4L * (2L * a.NonZeroCount + m + n + 2);
        if (bytes > (long)(_ctx.Device.GlobalMemory * 0.8)) { reason = $"system needs {bytes >> 20} MB of the device's {_ctx.Device.GlobalMemory >> 20} MB"; return null; }
        if (8L * a.NonZeroCount > (long)_ctx.Device.MaxAlloc) { reason = "matrix larger than the device's largest allocation"; return null; }

        var cl = _ctx.Cl;
        var buffers = new List<nint>();
        nint Buffer(long size, void* host, bool readOnly)
        {
            int e;
            var flags = (readOnly ? MemFlags.ReadOnly : MemFlags.ReadWrite) | (host != null ? MemFlags.CopyHostPtr : MemFlags.None);
            var buf = cl.CreateBuffer(_ctx.Context, flags, (nuint)Math.Max(8, size), host, &e);
            if (e != 0) throw new InvalidOperationException($"buffer of {size} bytes ({e})");
            buffers.Add(buf);
            return buf;
        }
        try
        {
            nint aRs, aCi, aVa, tRs, tCi, tVa, bu, bv, bw, bx, bTm, bTn;
            // A paired device takes each double as the floats hi, lo: the same 8 bytes a value.
            var aValues = _fp64 ? a.Values : Pairs(a.Values);
            var tValues = _fp64 ? at.Values : Pairs(at.Values);
            var bValues = _fp64 ? b : Pairs(b);
            fixed (int* p1 = a.RowStart) fixed (int* p2 = a.ColumnIndex) fixed (double* p3 = aValues)
            fixed (int* q1 = at.RowStart) fixed (int* q2 = at.ColumnIndex) fixed (double* q3 = tValues)
            fixed (double* pb = bValues)
            {
                aRs = Buffer(4L * (m + 1), p1, true); aCi = Buffer(4L * a.NonZeroCount, p2, true); aVa = Buffer(8L * a.NonZeroCount, p3, true);
                tRs = Buffer(4L * (n + 1), q1, true); tCi = Buffer(4L * at.NonZeroCount, q2, true); tVa = Buffer(8L * at.NonZeroCount, q3, true);
                bu = Buffer(8L * m, pb, false);
            }
            // x starts at zero and w is written before it is read, but an uninitialised buffer may hold
            // NaN, and NaN·0 is NaN: every vector starts from zeros.
            fixed (double* pz = x)
            {
                bv = Buffer(8L * n, pz, false);
                bw = Buffer(8L * n, pz, false);
                bx = Buffer(8L * n, pz, false);
                bTn = Buffer(8L * n, pz, false);
            }
            bTm = Buffer(8L * m, null, false);

            // The host solver in device calls; comments name the host line each call replaces.
            var beta = Norm(bu, m);
            if (beta == 0) return new LsqrResult(x, 0, 0, 0, 0, "zero right-hand side");
            var bnorm = beta;
            Scal(bu, m, 1.0 / beta);
            Spmv(n, tRs, tCi, tVa, bu, bv);                         // v = Aᵀu
            var alpha = Norm(bv, n);
            if (alpha == 0) return new LsqrResult(x, 0, bnorm, bnorm, 0, "Aᵀb = 0: the data do not see the model");
            Scal(bv, n, 1.0 / alpha);
            Copy(bv, bw, n);                                        // w = v

            double rhobar = alpha, phibar = beta;
            double ddnorm = 0, xxnorm = 0, z = 0, cs2 = -1, sn2 = 0, anorm2 = 0, dampedRes2 = 0;
            double rnorm = beta, acond = 0;
            var ctol = conlim > 0 ? 1.0 / conlim : 0.0;
            var stop = "iteration budget exhausted";
            var iterations = 0;

            for (var it = 0; it < maxIterations; it++)
            {
                ct.ThrowIfCancellationRequested();
                iterations = it + 1;

                Spmv(m, aRs, aCi, aVa, bv, bTm);                    // tmpM = A v
                Xmby(bTm, alpha, bu, m);                            // u = tmpM − α u
                beta = Norm(bu, m);
                if (beta > 0)
                {
                    Scal(bu, m, 1.0 / beta);
                    anorm2 += alpha * alpha + beta * beta + damp * damp;
                    Spmv(n, tRs, tCi, tVa, bu, bTn);                // tmpN = Aᵀ u
                    Xmby(bTn, beta, bv, n);                         // v = tmpN − β v
                    alpha = Norm(bv, n);
                    if (alpha > 0) Scal(bv, n, 1.0 / alpha);
                }

                var rhobar1 = Math.Sqrt(rhobar * rhobar + damp * damp);
                var cs1 = rhobar / rhobar1;
                var sn1 = damp / rhobar1;
                var psi = sn1 * phibar;
                phibar = cs1 * phibar;

                var rho = Math.Sqrt(rhobar1 * rhobar1 + beta * beta);
                var cs = rhobar1 / rho;
                var sn = beta / rho;
                var theta = sn * alpha;
                rhobar = -cs * alpha;
                var phi = cs * phibar;
                phibar = sn * phibar;
                var tau = sn * phi;

                ddnorm += Dot(bw, bw, n) / (rho * rho);
                Axpy(phi / rho, bw, bx, n);                         // x += (φ/ρ) w
                Xpby(bv, -theta / rho, bw, n);                      // w = v − (θ/ρ) w

                var delta = sn2 * rho;
                var gambar = -cs2 * rho;
                var rhs = phi - delta * z;
                var gamma = Math.Sqrt(gambar * gambar + theta * theta);
                cs2 = gambar / gamma;
                sn2 = theta / gamma;
                z = rhs / gamma;
                xxnorm += z * z;

                dampedRes2 += psi * psi;
                var anorm = Math.Sqrt(anorm2);
                acond = anorm * Math.Sqrt(ddnorm);
                var arnorm = alpha * Math.Abs(tau);
                rnorm = Math.Sqrt(phibar * phibar + dampedRes2);
                var xnorm = Math.Sqrt(xxnorm);

                var test1 = rnorm / bnorm;
                var test2 = arnorm / Math.Max(1e-300, anorm * rnorm);
                var test3 = 1.0 / Math.Max(1e-300, acond);
                var t1 = test1 / (1 + anorm * xnorm / bnorm);
                var rtol = btol + atol * anorm * xnorm / bnorm;

                onIteration?.Invoke(iterations, rnorm);

                if (test1 <= rtol) { stop = "residual reduced to tolerance (btol)"; break; }
                if (test2 <= atol) { stop = "normal equations satisfied (atol)"; break; }
                if (ctol > 0 && test3 <= ctol) { stop = "condition-number limit reached (conlim)"; break; }
                if (1 + t1 <= 1) { stop = "residual at machine precision"; break; }
                if (1 + test2 <= 1) { stop = "normal equations at machine precision"; break; }
                if (1 + test3 <= 1) { stop = "condition number at machine precision"; break; }
                if (beta == 0 || alpha == 0) { stop = "bidiagonalisation terminated"; break; }
            }

            fixed (double* px = x) Check(cl.EnqueueReadBuffer(_ctx.Queue, bx, true, 0, (nuint)(8L * n), px, 0, null, null), "read x");
            if (!_fp64) Unpair(x);
            return new LsqrResult(x, iterations, bnorm, rnorm, acond, stop);
        }
        finally
        {
            foreach (var buf in buffers) cl.ReleaseMemObject(buf);
        }
    }

    // ---- Kernel launches ----

    private static void Check(int err, string what)
    {
        if (err != 0) throw new InvalidOperationException($"{what} failed ({err})");
    }

    private void Launch(nint kernel, int n)
    {
        var global = (nuint)(((n + GroupSize - 1) / GroupSize) * GroupSize);
        var local = (nuint)GroupSize;
        Check(_ctx.Cl.EnqueueNdrangeKernel(_ctx.Queue, kernel, 1, null, &global, &local, 0, null, null), "kernel launch");
    }

    private void Arg(nint kernel, uint i, nint buffer) => Check(_ctx.Cl.SetKernelArg(kernel, i, (nuint)sizeof(nint), &buffer), "argument");
    private void Arg(nint kernel, uint i, int v) => Check(_ctx.Cl.SetKernelArg(kernel, i, sizeof(int), &v), "argument");
    private void Arg(nint kernel, uint i, double v)
    {
        if (!_fp64) v = PairBits(v);
        Check(_ctx.Cl.SetKernelArg(kernel, i, sizeof(double), &v), "argument");
    }

    /// <summary>The floats hi, lo whose sum is <paramref name="v"/> to about 48 bits, in the 8 bytes of a double.</summary>
    private static double PairBits(double v)
    {
        var hi = (float)v;
        var lo = float.IsFinite(hi) ? (float)(v - hi) : 0f;
        return BitConverter.Int64BitsToDouble((long)(uint)BitConverter.SingleToInt32Bits(hi) | ((long)BitConverter.SingleToInt32Bits(lo) << 32));
    }

    /// <summary>A copy of <paramref name="data"/> with each value as its float pair.</summary>
    private static double[] Pairs(double[] data)
    {
        var r = new double[data.Length];
        for (var i = 0; i < data.Length; i++) r[i] = PairBits(data[i]);
        return r;
    }

    /// <summary>Float pairs read back into a double array, summed in place.</summary>
    private static void Unpair(double[] data)
    {
        for (var i = 0; i < data.Length; i++)
        {
            var bits = BitConverter.DoubleToInt64Bits(data[i]);
            data[i] = (double)BitConverter.Int32BitsToSingle((int)bits) + BitConverter.Int32BitsToSingle((int)(bits >> 32));
        }
    }

    private void Spmv(int rows, nint rs, nint ci, nint va, nint x, nint y)
    {
        Arg(_spmv, 0, rows); Arg(_spmv, 1, rs); Arg(_spmv, 2, ci); Arg(_spmv, 3, va); Arg(_spmv, 4, x); Arg(_spmv, 5, y);
        Launch(_spmv, rows);
    }

    private void Xmby(nint x, double b, nint y, int n)
    {
        Arg(_xmby, 0, n); Arg(_xmby, 1, x); Arg(_xmby, 2, b); Arg(_xmby, 3, y);
        Launch(_xmby, n);
    }

    private void Xpby(nint x, double b, nint y, int n)
    {
        Arg(_xpby, 0, n); Arg(_xpby, 1, x); Arg(_xpby, 2, b); Arg(_xpby, 3, y);
        Launch(_xpby, n);
    }

    private void Axpy(double a, nint x, nint y, int n)
    {
        Arg(_axpy, 0, n); Arg(_axpy, 1, a); Arg(_axpy, 2, x); Arg(_axpy, 3, y);
        Launch(_axpy, n);
    }

    private void Scal(nint x, int n, double s)
    {
        Arg(_scal, 0, n); Arg(_scal, 1, x); Arg(_scal, 2, s);
        Launch(_scal, n);
    }

    /// <summary>y ← x.</summary>
    private void Copy(nint x, nint y, int n) =>
        Check(_ctx.Cl.EnqueueCopyBuffer(_ctx.Queue, x, y, 0, 0, (nuint)(8L * n), 0, null, null), "copy");

    private readonly double[] _partial = new double[Groups];

    private double Dot(nint a, nint b, int n)
    {
        Arg(_dot, 0, n); Arg(_dot, 1, a); Arg(_dot, 2, b);
        var part = PartialBuffer;
        Arg(_dot, 3, part);
        Check(_ctx.Cl.SetKernelArg(_dot, 4, (nuint)(GroupSize * sizeof(double)), null), "local argument");
        var global = (nuint)(Groups * GroupSize);
        var local = (nuint)GroupSize;
        Check(_ctx.Cl.EnqueueNdrangeKernel(_ctx.Queue, _dot, 1, null, &global, &local, 0, null, null), "dot launch");
        fixed (double* p = _partial) Check(_ctx.Cl.EnqueueReadBuffer(_ctx.Queue, part, true, 0, (nuint)(8 * Groups), p, 0, null, null), "read partial sums");
        if (!_fp64) Unpair(_partial);
        double s = 0;
        for (var g = 0; g < Groups; g++) s += _partial[g];
        return s;
    }

    private nint _partialBuf;

    private nint PartialBuffer => _partialBuf != 0 ? _partialBuf : _partialBuf = CreatePartial();

    private nint CreatePartial()
    {
        int e;
        var buf = _ctx.Cl.CreateBuffer(_ctx.Context, MemFlags.ReadWrite, (nuint)(8 * Groups), null, &e);
        Check(e, "partial-sum buffer");
        return buf;
    }

    private double Norm(nint x, int n) => Math.Sqrt(Math.Max(0, Dot(x, x, n)));

    // ---- Self-test ----

    /// <summary>
    /// Largest difference between the device and host solutions of a random sparse damped system,
    /// relative to the largest component of the host solution.
    /// </summary>
    public double SelfTest(out string message)
    {
        var (a, b) = TestSystem(3000, 1200, 12, 17);
        var host = Lsqr.Solve(a, b, 0.05, 300, 1e-12, 1e-12, 1e10);
        var dev = TrySolve(a, b, 0.05, 300, 1e-12, 1e-12, 1e10, default, null, out var reason);
        if (dev == null) { message = reason; return double.PositiveInfinity; }
        double worst = 0, scale = 0;
        for (var i = 0; i < host.Solution.Length; i++)
        {
            worst = Math.Max(worst, Math.Abs(dev.Solution[i] - host.Solution[i]));
            scale = Math.Max(scale, Math.Abs(host.Solution[i]));
        }
        var rel = worst / Math.Max(1e-300, scale);
        message = $"max |Δx| {rel:E2} of max |x| after {dev.Iterations} iterations (host {host.Iterations})";
        return rel;
    }

    /// <summary>A random sparse system with a few non-zeros per row, as a ray matrix has.</summary>
    public static (CsrMatrix A, double[] B) TestSystem(int rows, int cols, int perRow, int seed)
    {
        var rnd = new Random(seed);
        var b = new CsrBuilder(cols);
        var used = new HashSet<int>();
        for (var i = 0; i < rows; i++)
        {
            used.Clear();
            while (used.Count < perRow) used.Add(rnd.Next(cols));
            var c = used.OrderBy(v => v).ToArray();
            var v = c.Select(_ => rnd.NextDouble() + 0.1).ToArray();
            b.AddRow(c, v, rnd.NextDouble() - 0.5);
        }
        return b.Build();
    }

    public void Dispose()
    {
        var cl = _ctx.Cl;
        foreach (var k in new[] { _spmv, _xmby, _xpby, _axpy, _scal, _dot })
            if (k != 0) cl.ReleaseKernel(k);
        if (_partialBuf != 0) cl.ReleaseMemObject(_partialBuf);
        if (_program != 0) cl.ReleaseProgram(_program);
        _ctx.Dispose();
    }
}
