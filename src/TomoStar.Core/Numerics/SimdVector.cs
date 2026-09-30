using System.Numerics;
using System.Runtime.InteropServices;

namespace TomoStar.Core.Numerics;

/// <summary>
/// Dense vector kernels written against <see cref="Vector{T}"/>, so they use the widest SIMD unit
/// the machine has (SSE/AVX2/AVX-512 on x64, NEON on ARM64 including Apple silicon) without a
/// platform-specific code path. The scalar tail handles lengths that are not a multiple of the
/// vector width.
///
/// Reductions accumulate in several SIMD lanes and are therefore summed in a different order from
/// a scalar loop: results agree to rounding, not bit for bit.
/// </summary>
public static class SimdVector
{
    public static double Dot(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
    {
        if (a.Length != b.Length) throw new ArgumentException("Length mismatch.");
        var va = MemoryMarshal.Cast<double, Vector<double>>(a);
        var vb = MemoryMarshal.Cast<double, Vector<double>>(b);
        var acc = Vector<double>.Zero;
        for (var i = 0; i < va.Length; i++) acc += va[i] * vb[i];
        var sum = Vector.Dot(acc, Vector<double>.One);
        for (var i = va.Length * Vector<double>.Count; i < a.Length; i++) sum += a[i] * b[i];
        return sum;
    }

    /// <summary>Euclidean norm, with scaling so it neither overflows nor underflows.</summary>
    public static double Norm(ReadOnlySpan<double> a)
    {
        var max = MaxAbs(a);
        if (max == 0 || !double.IsFinite(max)) return max;
        if (max > 1e150 || max < 1e-150)
        {
            double s = 0;
            var inv = 1.0 / max;
            foreach (var x in a) { var y = x * inv; s += y * y; }
            return max * Math.Sqrt(s);
        }
        return Math.Sqrt(Dot(a, a));
    }

    public static double MaxAbs(ReadOnlySpan<double> a)
    {
        var va = MemoryMarshal.Cast<double, Vector<double>>(a);
        var acc = Vector<double>.Zero;
        for (var i = 0; i < va.Length; i++) acc = Vector.Max(acc, Vector.Abs(va[i]));
        double m = 0;
        for (var k = 0; k < Vector<double>.Count; k++) m = Math.Max(m, acc[k]);
        for (var i = va.Length * Vector<double>.Count; i < a.Length; i++) m = Math.Max(m, Math.Abs(a[i]));
        return m;
    }

    /// <summary>a ← a·s</summary>
    public static void Scale(Span<double> a, double s)
    {
        var va = MemoryMarshal.Cast<double, Vector<double>>(a);
        var vs = new Vector<double>(s);
        for (var i = 0; i < va.Length; i++) va[i] *= vs;
        for (var i = va.Length * Vector<double>.Count; i < a.Length; i++) a[i] *= s;
    }

    /// <summary>y ← y + α·x</summary>
    public static void Axpy(double alpha, ReadOnlySpan<double> x, Span<double> y)
    {
        if (x.Length != y.Length) throw new ArgumentException("Length mismatch.");
        var vx = MemoryMarshal.Cast<double, Vector<double>>(x);
        var vy = MemoryMarshal.Cast<double, Vector<double>>(y);
        var va = new Vector<double>(alpha);
        for (var i = 0; i < vx.Length; i++) vy[i] += va * vx[i];
        for (var i = vx.Length * Vector<double>.Count; i < x.Length; i++) y[i] += alpha * x[i];
    }

    /// <summary>y ← x + β·y</summary>
    public static void Xpby(ReadOnlySpan<double> x, double beta, Span<double> y)
    {
        if (x.Length != y.Length) throw new ArgumentException("Length mismatch.");
        var vx = MemoryMarshal.Cast<double, Vector<double>>(x);
        var vy = MemoryMarshal.Cast<double, Vector<double>>(y);
        var vb = new Vector<double>(beta);
        for (var i = 0; i < vx.Length; i++) vy[i] = vx[i] + vb * vy[i];
        for (var i = vx.Length * Vector<double>.Count; i < x.Length; i++) y[i] = x[i] + beta * y[i];
    }

    /// <summary>y ← x − β·y</summary>
    public static void Xmby(ReadOnlySpan<double> x, double beta, Span<double> y)
    {
        if (x.Length != y.Length) throw new ArgumentException("Length mismatch.");
        var vx = MemoryMarshal.Cast<double, Vector<double>>(x);
        var vy = MemoryMarshal.Cast<double, Vector<double>>(y);
        var vb = new Vector<double>(beta);
        for (var i = 0; i < vx.Length; i++) vy[i] = vx[i] - vb * vy[i];
        for (var i = vx.Length * Vector<double>.Count; i < x.Length; i++) y[i] = x[i] - beta * y[i];
    }

    public static bool IsHardwareAccelerated => Vector.IsHardwareAccelerated;
    public static int Width => Vector<double>.Count;
}
