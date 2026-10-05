// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using System.Text;
using Silk.NET.OpenCL;

namespace TomoStar.Core.Compute;

/// <summary>An OpenCL device as enumerated on this machine.</summary>
public sealed record OpenClDeviceInfo(nint Platform, nint Device, string Name, string Vendor, string PlatformName,
    ulong GlobalMemory, ulong MaxAlloc, bool Fp64, bool IsGpu)
{
    public override string ToString() => $"{Name} ({Vendor}, {GlobalMemory / (1024 * 1024)} MB{(Fp64 ? ", fp64" : "")})";
}

/// <summary>
/// One OpenCL context + queue on the best available device. Creating it never throws: a machine
/// with no OpenCL runtime (common on macOS since 10.14 deprecated it, and on headless Linux) simply
/// gets <see cref="TryCreate"/> = null and TomoSTAR stays on the CPU.
/// </summary>
public sealed unsafe class OpenClContext : IDisposable
{
    private bool _disposed;

    private OpenClContext(CL cl, OpenClDeviceInfo device, nint context, nint queue)
    {
        Cl = cl;
        Device = device;
        Context = context;
        Queue = queue;
    }

    public CL Cl { get; }
    public OpenClDeviceInfo Device { get; }
    public nint Context { get; }
    public nint Queue { get; }

    public static IReadOnlyList<OpenClDeviceInfo> Enumerate(out string? error)
    {
        error = null;
        var list = new List<OpenClDeviceInfo>();
        CL cl;
        try { cl = CL.GetApi(); }
        catch (Exception ex) { error = $"No OpenCL runtime: {ex.Message}"; return list; }
        try
        {
            uint np = 0;
            if (cl.GetPlatformIDs(0, null, &np) != 0 || np == 0) { error = "No OpenCL platform."; return list; }
            var platforms = new nint[np];
            fixed (nint* pp = platforms) cl.GetPlatformIDs(np, pp, null);
            foreach (var p in platforms)
            {
                uint nd = 0;
                if (cl.GetDeviceIDs(p, DeviceType.All, 0, null, &nd) != 0 || nd == 0) continue;
                var devices = new nint[nd];
                fixed (nint* pd = devices) cl.GetDeviceIDs(p, DeviceType.All, nd, pd, null);
                foreach (var d in devices)
                {
                    ulong mem = 0, alloc = 0, type = 0;
                    cl.GetDeviceInfo(d, DeviceInfo.GlobalMemSize, sizeof(ulong), &mem, null);
                    cl.GetDeviceInfo(d, DeviceInfo.MaxMemAllocSize, sizeof(ulong), &alloc, null);
                    cl.GetDeviceInfo(d, DeviceInfo.Type, sizeof(ulong), &type, null);
                    var ext = DeviceString(cl, d, DeviceInfo.Extensions);
                    list.Add(new OpenClDeviceInfo(p, d,
                        DeviceString(cl, d, DeviceInfo.Name).Trim(),
                        DeviceString(cl, d, DeviceInfo.Vendor).Trim(),
                        PlatformString(cl, p, PlatformInfo.Name).Trim(),
                        mem, alloc, ext.Contains("cl_khr_fp64", StringComparison.Ordinal),
                        (type & (ulong)DeviceType.Gpu) != 0));
                }
            }
        }
        catch (Exception ex)
        {
            error = $"OpenCL enumeration failed: {ex.Message}";
        }
        finally
        {
            cl.Dispose();
        }
        return list;
    }

    /// <summary>
    /// The best device: GPUs first, then (when single precision will do) those with fp64, then by
    /// global memory. Null when there is none.
    /// </summary>
    /// <param name="requireFp64">
    /// Only devices with cl_khr_fp64. Solvers that carry a single-precision kernel pass false and
    /// read <see cref="OpenClDeviceInfo.Fp64"/> to choose it: Apple GPUs and most ARM GPUs have no
    /// double precision at all.
    /// </param>
    /// <param name="allowCpuDevice">
    /// Accept an OpenCL device that is itself a CPU (e.g. POCL). Such a device runs the relaxation
    /// kernel far slower than the native multithreaded fast marching, so computations leave it out;
    /// self-tests and diagnostics use it to check the kernel.
    /// </param>
    public static OpenClContext? TryCreate(Action<string>? log = null, bool requireFp64 = true, bool allowCpuDevice = false)
    {
        var devices = Enumerate(out var error);
        if (error != null) log?.Invoke(error);
        var pick = devices.Where(d => (d.Fp64 || !requireFp64) && (d.IsGpu || allowCpuDevice))
            .OrderByDescending(d => d.IsGpu).ThenByDescending(d => d.Fp64).ThenByDescending(d => d.GlobalMemory).FirstOrDefault();
        if (pick == null)
        {
            log?.Invoke(devices.Count == 0 ? "OpenCL: no device; using CPU."
                : devices.All(d => !d.IsGpu) && !allowCpuDevice ? "OpenCL: only CPU devices, slower than the native solver; using CPU."
                : "OpenCL: no device with double precision (cl_khr_fp64); using CPU.");
            return null;
        }
        CL? cl = null;
        try
        {
            cl = CL.GetApi();
            int err;
            var dev = pick.Device;
            var ctx = cl.CreateContext(null, 1, &dev, null, null, &err);
            if (err != 0) throw new InvalidOperationException($"CreateContext failed ({err}).");
            var queue = cl.CreateCommandQueue(ctx, dev, CommandQueueProperties.None, &err);
            if (err != 0)
            {
                cl.ReleaseContext(ctx);
                throw new InvalidOperationException($"CreateCommandQueue failed ({err}).");
            }
            log?.Invoke($"OpenCL device: {pick}");
            return new OpenClContext(cl, pick, ctx, queue);
        }
        catch (Exception ex)
        {
            log?.Invoke($"OpenCL initialisation failed ({ex.Message}); using CPU.");
            cl?.Dispose();
            return null;
        }
    }

    /// <summary>Builds a program; returns 0 and logs the build log on failure.</summary>
    public nint BuildProgram(string source, Action<string>? log, string options = "")
    {
        int err;
        var bytes = Encoding.ASCII.GetBytes(source);
        nint program;
        fixed (byte* p = bytes)
        {
            var pp = p;
            var len = (nuint)bytes.Length;
            program = Cl.CreateProgramWithSource(Context, 1, &pp, &len, &err);
        }
        if (err != 0) { log?.Invoke($"OpenCL: CreateProgramWithSource failed ({err})."); return 0; }
        var dev = Device.Device;
        var opt = Encoding.ASCII.GetBytes(options + "\0");
        fixed (byte* po = opt) err = Cl.BuildProgram(program, 1, &dev, po, null, null);
        if (err != 0)
        {
            nuint size = 0;
            Cl.GetProgramBuildInfo(program, dev, ProgramBuildInfo.BuildLog, 0, null, &size);
            var buf = new byte[Math.Max(1, (int)size)];
            fixed (byte* pb = buf) Cl.GetProgramBuildInfo(program, dev, ProgramBuildInfo.BuildLog, size, pb, null);
            log?.Invoke($"OpenCL build failed ({err}): {Encoding.ASCII.GetString(buf).Trim('\0', ' ', '\n')}");
            Cl.ReleaseProgram(program);
            return 0;
        }
        return program;
    }

    public nint CreateKernel(nint program, string name)
    {
        int err;
        var bytes = Encoding.ASCII.GetBytes(name + "\0");
        nint k;
        fixed (byte* p = bytes) k = Cl.CreateKernel(program, p, &err);
        return err == 0 ? k : 0;
    }

    private static string DeviceString(CL cl, nint d, DeviceInfo info)
    {
        nuint size = 0;
        cl.GetDeviceInfo(d, info, 0, null, &size);
        if (size == 0) return "";
        var buf = new byte[size];
        fixed (byte* p = buf) cl.GetDeviceInfo(d, info, size, p, null);
        return Encoding.ASCII.GetString(buf).TrimEnd('\0');
    }

    private static string PlatformString(CL cl, nint p, PlatformInfo info)
    {
        nuint size = 0;
        cl.GetPlatformInfo(p, info, 0, null, &size);
        if (size == 0) return "";
        var buf = new byte[size];
        fixed (byte* b = buf) cl.GetPlatformInfo(p, info, size, b, null);
        return Encoding.ASCII.GetString(buf).TrimEnd('\0');
    }

    public void Dispose()
    {
        if (_disposed) return;
        _disposed = true;
        Cl.ReleaseCommandQueue(Queue);
        Cl.ReleaseContext(Context);
        Cl.Dispose();
    }
}
