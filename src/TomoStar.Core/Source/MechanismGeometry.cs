// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

namespace TomoStar.Core.Source;

/// <summary>A fault plane and slip: strike, dip and rake in degrees (Aki &amp; Richards convention).</summary>
public readonly record struct NodalPlane(double Strike, double Dip, double Rake);

/// <summary>A principal axis: trend (azimuth, degrees from north) and plunge (degrees below the horizontal).</summary>
public readonly record struct Axis(double Trend, double Plunge);

/// <summary>
/// A moment tensor in the north–east–down frame (x = north, y = east, z = down), six independent
/// components, N·m (or normalised).
/// </summary>
public readonly record struct MomentTensor(double Nn, double Ee, double Dd, double Ne, double Nd, double Ed)
{
    public double[,] Matrix => new[,] { { Nn, Ne, Nd }, { Ne, Ee, Ed }, { Nd, Ed, Dd } };

    public static MomentTensor operator *(MomentTensor m, double k) => new(m.Nn * k, m.Ee * k, m.Dd * k, m.Ne * k, m.Nd * k, m.Ed * k);

    /// <summary>γᵀ M γ for a unit direction γ (north, east, down): the far-field P radiation.</summary>
    public double Radiation(double n, double e, double d) =>
        Nn * n * n + Ee * e * e + Dd * d * d + 2 * (Ne * n * e + Nd * n * d + Ed * e * d);

    /// <summary>
    /// The components in the up–south–east (r, θ, φ) convention of the Global CMT catalogue:
    /// Mrr = Mdd, Mθθ = Mnn, Mφφ = Mee, Mrθ = Mnd, Mrφ = −Med, Mθφ = −Mne.
    /// </summary>
    public (double Mrr, double Mtt, double Mpp, double Mrt, double Mrp, double Mtp) ToUse() => (Dd, Nn, Ee, Nd, -Ed, -Ne);
}

/// <summary>
/// Geometry of earthquake sources: double-couple planes, their tensor, principal axes, tensor
/// decomposition and the angle between two double couples. Formulae of Aki &amp; Richards
/// (Quantitative Seismology, 2nd ed., 2002, Box 4.4 and §4.2) in the north–east–down frame;
/// decomposition after Jost &amp; Herrmann (1989, Seismol. Res. Lett. 60(2), 37–57,
/// doi:10.1785/gssrl.60.2.37); minimum rotation angle after Kagan (1991, Geophys. J. Int. 106(3),
/// 709–716, doi:10.1111/j.1365-246X.1991.tb06343.x).
/// </summary>
public static class MechanismGeometry
{
    private const double D2R = Math.PI / 180, R2D = 180 / Math.PI;

    /// <summary>Unit-moment double-couple tensor of a fault plane and slip (Aki &amp; Richards Box 4.4).</summary>
    public static MomentTensor Tensor(NodalPlane p)
    {
        double f = p.Strike * D2R, d = p.Dip * D2R, l = p.Rake * D2R;
        double sd = Math.Sin(d), cd = Math.Cos(d), s2d = Math.Sin(2 * d), c2d = Math.Cos(2 * d);
        double sl = Math.Sin(l), cl = Math.Cos(l), sf = Math.Sin(f), cf = Math.Cos(f), s2f = Math.Sin(2 * f), c2f = Math.Cos(2 * f);
        return new MomentTensor(
            Nn: -(sd * cl * s2f + s2d * sl * sf * sf),
            Ee: sd * cl * s2f - s2d * sl * cf * cf,
            Dd: s2d * sl,
            Ne: sd * cl * c2f + 0.5 * s2d * sl * s2f,
            Nd: -(cd * cl * cf + c2d * sl * sf),
            Ed: -(cd * cl * sf - c2d * sl * cf));
    }

    /// <summary>Fault normal (pointing into the hanging wall, upwards) and slip vector of a plane, NED.</summary>
    public static (double[] Normal, double[] Slip) Vectors(NodalPlane p)
    {
        double f = p.Strike * D2R, d = p.Dip * D2R, l = p.Rake * D2R;
        double[] n = [-Math.Sin(d) * Math.Sin(f), Math.Sin(d) * Math.Cos(f), -Math.Cos(d)];
        double[] u = [Math.Cos(l) * Math.Cos(f) + Math.Sin(l) * Math.Cos(d) * Math.Sin(f),
            Math.Cos(l) * Math.Sin(f) - Math.Sin(l) * Math.Cos(d) * Math.Cos(f), -Math.Sin(l) * Math.Sin(d)];
        return (n, u);
    }

    /// <summary>The plane whose normal is <paramref name="n"/> and slip <paramref name="u"/> (NED unit vectors).</summary>
    public static NodalPlane FromVectors(double[] n, double[] u)
    {
        n = (double[])n.Clone();
        u = (double[])u.Clone();
        // The normal must point up (hanging wall); flipping both keeps the same double couple.
        if (n[2] > 0) { for (var i = 0; i < 3; i++) { n[i] = -n[i]; u[i] = -u[i]; } }
        var dip = Math.Acos(Math.Clamp(-n[2], -1, 1));
        double strike, rake;
        if (Math.Sin(dip) < 1e-9)
        {
            // Horizontal plane: strike is arbitrary; take the one that makes the rake zero.
            strike = Math.Atan2(u[1], u[0]);
            rake = 0;
        }
        else
        {
            strike = Math.Atan2(-n[0], n[1]);
            rake = Math.Atan2(-u[2] / Math.Sin(dip), u[0] * Math.Cos(strike) + u[1] * Math.Sin(strike));
        }
        return new NodalPlane(Wrap360(strike * R2D), dip * R2D, WrapRake(rake * R2D));
    }

    /// <summary>The other nodal plane of a double couple (normal and slip exchanged).</summary>
    public static NodalPlane Auxiliary(NodalPlane p)
    {
        var (n, u) = Vectors(p);
        return FromVectors(u, n);
    }

    /// <summary>P, T and B axes of a double-couple plane.</summary>
    public static (Axis P, Axis T, Axis B) Axes(NodalPlane p)
    {
        var (n, u) = Vectors(p);
        double[] t = [(n[0] + u[0]) / Math.Sqrt(2), (n[1] + u[1]) / Math.Sqrt(2), (n[2] + u[2]) / Math.Sqrt(2)];
        double[] pp = [(n[0] - u[0]) / Math.Sqrt(2), (n[1] - u[1]) / Math.Sqrt(2), (n[2] - u[2]) / Math.Sqrt(2)];
        return (ToAxis(pp), ToAxis(t), ToAxis(Cross(n, u)));
    }

    /// <summary>Trend and plunge of a direction (NED); the downward-pointing sense is used.</summary>
    public static Axis ToAxis(double[] v)
    {
        var x = v;
        if (x[2] < 0) x = [-v[0], -v[1], -v[2]];
        var len = Math.Sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
        return new Axis(Wrap360(Math.Atan2(x[1], x[0]) * R2D), Math.Asin(Math.Clamp(x[2] / len, -1, 1)) * R2D);
    }

    /// <summary>Unit vector (NED) of an axis.</summary>
    public static double[] FromAxis(Axis a)
    {
        double t = a.Trend * D2R, p = a.Plunge * D2R;
        return [Math.Cos(p) * Math.Cos(t), Math.Cos(p) * Math.Sin(t), Math.Sin(p)];
    }

    /// <summary>
    /// Eigen-decomposition of a moment tensor: eigenvalues ascending (P, B, T) with their vectors.
    /// </summary>
    public static (double[] Values, double[][] Vectors) Eigen(MomentTensor m)
    {
        var a = m.Matrix;
        var v = new double[,] { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
        // Cyclic Jacobi rotations: exact enough for a 3×3 symmetric matrix after a few sweeps.
        for (var sweep = 0; sweep < 50; sweep++)
        {
            var off = a[0, 1] * a[0, 1] + a[0, 2] * a[0, 2] + a[1, 2] * a[1, 2];
            if (off < 1e-30) break;
            for (var p = 0; p < 2; p++)
            for (var q = p + 1; q < 3; q++)
            {
                if (Math.Abs(a[p, q]) < 1e-300) continue;
                var theta = 0.5 * Math.Atan2(2 * a[p, q], a[q, q] - a[p, p]);
                double c = Math.Cos(theta), s = Math.Sin(theta);
                for (var k = 0; k < 3; k++)
                {
                    double akp = a[k, p], akq = a[k, q];
                    a[k, p] = c * akp - s * akq;
                    a[k, q] = s * akp + c * akq;
                }
                for (var k = 0; k < 3; k++)
                {
                    double apk = a[p, k], aqk = a[q, k];
                    a[p, k] = c * apk - s * aqk;
                    a[q, k] = s * apk + c * aqk;
                }
                for (var k = 0; k < 3; k++)
                {
                    double vkp = v[k, p], vkq = v[k, q];
                    v[k, p] = c * vkp - s * vkq;
                    v[k, q] = s * vkp + c * vkq;
                }
            }
        }
        var order = new[] { 0, 1, 2 }.OrderBy(i => a[i, i]).ToArray();
        return (order.Select(i => a[i, i]).ToArray(), order.Select(i => new[] { v[0, i], v[1, i], v[2, i] }).ToArray());
    }

    /// <summary>The best double couple of a tensor: its two nodal planes, from the P and T eigenvectors.</summary>
    public static (NodalPlane Plane1, NodalPlane Plane2) BestDoubleCouple(MomentTensor m)
    {
        var (_, vec) = Eigen(m);
        double[] p = vec[0], t = vec[2];
        double[] n = [(t[0] + p[0]) / Math.Sqrt(2), (t[1] + p[1]) / Math.Sqrt(2), (t[2] + p[2]) / Math.Sqrt(2)];
        double[] u = [(t[0] - p[0]) / Math.Sqrt(2), (t[1] - p[1]) / Math.Sqrt(2), (t[2] - p[2]) / Math.Sqrt(2)];
        return (FromVectors(n, u), FromVectors(u, n));
    }

    /// <summary>
    /// Scalar moment, moment magnitude and the isotropic / double-couple / CLVD percentages of a
    /// tensor. M₀ = √(Σ Mᵢⱼ²/2) (Silver &amp; Jordan 1982); Mw = (2/3)(log₁₀ M₀ − 9.1) for M₀ in N·m
    /// (the IASPEI standard form of Hanks &amp; Kanamori 1979). The deviatoric part is split with
    /// ε = −λ_min|/|λ_max| of its eigenvalues ordered by absolute value (Jost &amp; Herrmann 1989):
    /// CLVD % = 200|ε|, DC % = 100 − CLVD %, both of the deviatoric part; ISO % is the share of the
    /// isotropic part in |ISO| + |deviatoric|.
    /// </summary>
    public static (double M0, double Mw, double IsoPercent, double DcPercent, double ClvdPercent) Decompose(MomentTensor m)
    {
        var m0 = Math.Sqrt((m.Nn * m.Nn + m.Ee * m.Ee + m.Dd * m.Dd + 2 * (m.Ne * m.Ne + m.Nd * m.Nd + m.Ed * m.Ed)) / 2);
        var iso = (m.Nn + m.Ee + m.Dd) / 3;
        var dev = new MomentTensor(m.Nn - iso, m.Ee - iso, m.Dd - iso, m.Ne, m.Nd, m.Ed);
        var (ev, _) = Eigen(dev);
        var byAbs = ev.OrderBy(Math.Abs).ToArray();
        var eps = Math.Abs(byAbs[2]) > 0 ? -byAbs[0] / Math.Abs(byAbs[2]) : 0;
        var clvd = Math.Min(100, 200 * Math.Abs(eps));
        var devSize = Math.Abs(byAbs[2]);
        var isoPercent = Math.Abs(iso) + devSize > 0 ? 100 * Math.Abs(iso) / (Math.Abs(iso) + devSize) : 0;
        var mw = m0 > 0 ? 2.0 / 3.0 * (Math.Log10(m0) - 9.1) : double.NaN;
        return (m0, mw, isoPercent, 100 - clvd, clvd);
    }

    /// <summary>
    /// Kagan's minimum rotation angle between two double couples, degrees (0 to 120): the smallest
    /// rotation that takes one set of principal axes onto the other, over the four symmetries of a
    /// double couple.
    /// </summary>
    public static double KaganAngle(NodalPlane a, NodalPlane b)
    {
        var ra = Triad(a);
        var rb = Triad(b);
        var best = 180.0;
        // Symmetries: identity and 180° rotations about T, B and P (sign changes of two axes).
        int[][] signs = [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]];
        foreach (var s in signs)
        {
            // R = Rb · diag(s) · Raᵀ; angle = acos((tr R − 1)/2).
            double tr = 0;
            for (var k = 0; k < 3; k++)
            for (var i = 0; i < 3; i++)
                tr += rb[i][k] * s[k] * ra[i][k];
            var angle = Math.Acos(Math.Clamp((tr - 1) / 2, -1, 1)) * R2D;
            best = Math.Min(best, angle);
        }
        return best;
    }

    /// <summary>Columns T, B, P of the principal-axis frame of a double couple (unit vectors, NED).</summary>
    private static double[][] Triad(NodalPlane p)
    {
        var (n, u) = Vectors(p);
        double[] t = [(n[0] + u[0]) / Math.Sqrt(2), (n[1] + u[1]) / Math.Sqrt(2), (n[2] + u[2]) / Math.Sqrt(2)];
        double[] pp = [(n[0] - u[0]) / Math.Sqrt(2), (n[1] - u[1]) / Math.Sqrt(2), (n[2] - u[2]) / Math.Sqrt(2)];
        var b = Cross(t, pp);
        // Rows are components, columns the axes: m[i][k] = component i of axis k.
        return [[t[0], b[0], pp[0]], [t[1], b[1], pp[1]], [t[2], b[2], pp[2]]];
    }

    /// <summary>
    /// Faulting style from the plunges of P, T and B (Frohlich 1992, Phys. Earth Planet. Inter.
    /// 75(1–3), 193–198, doi:10.1016/0031-9201(92)90130-N): thrust when T is steep, normal when P is
    /// steep, strike-slip when B is steep, otherwise oblique.
    /// </summary>
    public static string Style(NodalPlane p)
    {
        var (pa, ta, ba) = Axes(p);
        double sp = Math.Sin(pa.Plunge * D2R), st = Math.Sin(ta.Plunge * D2R), sb = Math.Sin(ba.Plunge * D2R);
        if (st * st > 0.59) return "Thrust";
        if (sp * sp > 0.75) return "Normal";
        if (sb * sb > 0.75) return "Strike-slip";
        return "Oblique";
    }

    public static double[] Cross(double[] a, double[] b) =>
        [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];

    public static double Wrap360(double x) => ((x % 360) + 360) % 360;

    public static double WrapRake(double x)
    {
        x = ((x + 180) % 360 + 360) % 360 - 180;
        return x <= -180 ? x + 360 : x;
    }
}
