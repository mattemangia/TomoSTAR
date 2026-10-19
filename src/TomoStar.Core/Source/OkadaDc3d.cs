// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using System;

namespace TomoStar.Core.Source;

// Clean C# translation of the finite-rectangle part of Y. Okada's DC3D
// reference implementation (Sep. 1991, revisions through May 2002).
// Equations: Okada (1992), BSSA 82(2), 1018-1040,
// https://doi.org/10.1785/BSSA0820021018.
// Comparator/source distribution: cutde-org/okada_wrapper, MIT license,
// https://github.com/cutde-org/okada_wrapper/blob/master/okada_wrapper/DC3D.f
// The translation uses doubles throughout and keeps DC3D's x=strike,
// y=horizontal up-dip, z=up (half-space z <= 0) convention.

internal readonly record struct OkadaDc3dResult(
    int Status, double X, double Y, double Z,
    double Xx, double Yx, double Zx,
    double Xy, double Yy, double Zy,
    double Xz, double Yz, double Zz)
{
    public bool IsValid => Status == 0 &&
        double.IsFinite(X) && double.IsFinite(Y) && double.IsFinite(Z) &&
        double.IsFinite(Xx) && double.IsFinite(Yx) && double.IsFinite(Zx) &&
        double.IsFinite(Xy) && double.IsFinite(Yy) && double.IsFinite(Zy) &&
        double.IsFinite(Xz) && double.IsFinite(Yz) && double.IsFinite(Zz);

    public double[,] Gradient() => new[,]
    {
        { Xx, Xy, Xz },
        { Yx, Yy, Yz },
        { Zx, Zy, Zz }
    };
}

/// <summary>
/// Okada (1992) DC3D finite rectangular dislocation in a homogeneous elastic
/// half-space. Returns displacement and all nine analytic first derivatives.
/// </summary>
internal static class OkadaDc3d
{
    private const double Eps = 1e-6;
    private const double TwoPi = 2.0 * Math.PI;

    private readonly record struct Medium(
        double A1, double A2, double A3, double A4, double A5,
        double Sd, double Cd, double Sd2, double Cd2, double SdCd);

    private readonly record struct Geometry(
        double Xi2, double Et2, double Q2, double R, double R2, double R3,
        double R5, double Y, double D, double Tt, double Alx, double Ale,
        double X11, double Y11, double X32, double Y32,
        double Ey, double Ez, double Fy, double Fz,
        double Gy, double Gz, double Hy, double Hz);

    public static OkadaDc3dResult Evaluate(
        double alpha, double x, double y, double z, double depth, double dipDeg,
        double alongMin, double alongMax, double downMin, double downMax,
        double strikeSlip, double dipSlip, double opening)
    {
        if (z > 0) return Empty(2);
        var medium = Constants(alpha, dipDeg);
        var u = new double[12];
        double[] xi = { Snap(x - alongMin), Snap(x - alongMax) };

        if (!Accumulate(realSource: true, depth + z, x, y, z, xi,
                downMin, downMax, strikeSlip, dipSlip, opening, medium, u))
            return Empty(1);
        if (!Accumulate(realSource: false, depth - z, x, y, z, xi,
                downMin, downMax, strikeSlip, dipSlip, opening, medium, u))
            return Empty(1);

        return new OkadaDc3dResult(0,
            u[0], u[1], u[2], u[3], u[4], u[5],
            u[6], u[7], u[8], u[9], u[10], u[11]);
    }

    private static bool Accumulate(bool realSource, double d, double x, double y,
        double z, double[] xi, double aw1, double aw2,
        double disl1, double disl2, double disl3, Medium c, double[] u)
    {
        double p = y * c.Cd + d * c.Sd;
        double q = Snap(y * c.Sd - d * c.Cd);
        double[] et = { Snap(p - aw1), Snap(p - aw2) };
        if (q == 0 && ((xi[0] * xi[1] <= 0 && et[0] * et[1] == 0) ||
                       (et[0] * et[1] <= 0 && xi[0] * xi[1] == 0)))
            return false;

        int[] kxi = { 0, 0 }, ket = { 0, 0 };
        double r12 = Math.Sqrt(xi[0] * xi[0] + et[1] * et[1] + q * q);
        double r21 = Math.Sqrt(xi[1] * xi[1] + et[0] * et[0] + q * q);
        double r22 = Math.Sqrt(xi[1] * xi[1] + et[1] * et[1] + q * q);
        if (xi[0] < 0 && r21 + xi[1] < Eps) kxi[0] = 1;
        if (xi[0] < 0 && r22 + xi[1] < Eps) kxi[1] = 1;
        if (et[0] < 0 && r12 + et[1] < Eps) ket[0] = 1;
        if (et[0] < 0 && r22 + et[1] < Eps) ket[1] = 1;

        for (int k = 0; k < 2; k++)
        for (int j = 0; j < 2; j++)
        {
            Geometry g = GeometryConstants(xi[j], et[k], q, c.Sd, c.Cd,
                kxi[k], ket[j]);
            if (!(g.R > 0)) return false;
            double[] a = PartA(xi[j], et[k], q, disl1, disl2, disl3, c, g);
            var du = new double[12];
            if (realSource)
            {
                for (int i = 0; i < 12; i += 3)
                {
                    du[i] = -a[i];
                    du[i + 1] = -a[i + 1] * c.Cd + a[i + 2] * c.Sd;
                    du[i + 2] = -a[i + 1] * c.Sd - a[i + 2] * c.Cd;
                    if (i == 9)
                    {
                        du[i] = -du[i];
                        du[i + 1] = -du[i + 1];
                        du[i + 2] = -du[i + 2];
                    }
                }
            }
            else
            {
                double[] b = PartB(xi[j], et[k], q, disl1, disl2, disl3, c, g);
                double[] cc = PartC(xi[j], et[k], q, z, disl1, disl2, disl3, c, g);
                for (int i = 0; i < 12; i += 3)
                {
                    du[i] = a[i] + b[i] + z * cc[i];
                    du[i + 1] = (a[i + 1] + b[i + 1] + z * cc[i + 1]) * c.Cd
                              - (a[i + 2] + b[i + 2] + z * cc[i + 2]) * c.Sd;
                    du[i + 2] = (a[i + 1] + b[i + 1] - z * cc[i + 1]) * c.Sd
                              + (a[i + 2] + b[i + 2] - z * cc[i + 2]) * c.Cd;
                    if (i == 9)
                    {
                        du[9] += cc[0];
                        du[10] += cc[1] * c.Cd - cc[2] * c.Sd;
                        du[11] -= cc[1] * c.Sd + cc[2] * c.Cd;
                    }
                }
            }
            double sign = j + k == 1 ? -1.0 : 1.0;
            for (int i = 0; i < 12; i++) u[i] += sign * du[i];
        }
        return true;
    }

    private static double[] PartA(double xi, double et, double q,
        double d1, double d2, double d3, Medium c, Geometry g)
    {
        var u = new double[12];
        var v = new double[12];
        double xy = xi * g.Y11, qx = q * g.X11, qy = q * g.Y11;
        if (d1 != 0)
        {
            v[0] = g.Tt / 2 + c.A2 * xi * qy;
            v[1] = c.A2 * q / g.R;
            v[2] = c.A1 * g.Ale - c.A2 * q * qy;
            v[3] = -c.A1 * qy - c.A2 * g.Xi2 * q * g.Y32;
            v[4] = -c.A2 * xi * q / g.R3;
            v[5] = c.A1 * xy + c.A2 * xi * g.Q2 * g.Y32;
            v[6] = c.A1 * xy * c.Sd + c.A2 * xi * g.Fy + g.D / 2 * g.X11;
            v[7] = c.A2 * g.Ey;
            v[8] = c.A1 * (c.Cd / g.R + qy * c.Sd) - c.A2 * q * g.Fy;
            v[9] = c.A1 * xy * c.Cd + c.A2 * xi * g.Fz + g.Y / 2 * g.X11;
            v[10] = c.A2 * g.Ez;
            v[11] = -c.A1 * (c.Sd / g.R - qy * c.Cd) - c.A2 * q * g.Fz;
            Add(u, v, d1 / TwoPi);
        }
        if (d2 != 0)
        {
            v[0] = c.A2 * q / g.R;
            v[1] = g.Tt / 2 + c.A2 * et * qx;
            v[2] = c.A1 * g.Alx - c.A2 * q * qx;
            v[3] = -c.A2 * xi * q / g.R3;
            v[4] = -qy / 2 - c.A2 * et * q / g.R3;
            v[5] = c.A1 / g.R + c.A2 * g.Q2 / g.R3;
            v[6] = c.A2 * g.Ey;
            v[7] = c.A1 * g.D * g.X11 + xy / 2 * c.Sd + c.A2 * et * g.Gy;
            v[8] = c.A1 * g.Y * g.X11 - c.A2 * q * g.Gy;
            v[9] = c.A2 * g.Ez;
            v[10] = c.A1 * g.Y * g.X11 + xy / 2 * c.Cd + c.A2 * et * g.Gz;
            v[11] = -c.A1 * g.D * g.X11 - c.A2 * q * g.Gz;
            Add(u, v, d2 / TwoPi);
        }
        if (d3 != 0)
        {
            v[0] = -c.A1 * g.Ale - c.A2 * q * qy;
            v[1] = -c.A1 * g.Alx - c.A2 * q * qx;
            v[2] = g.Tt / 2 - c.A2 * (et * qx + xi * qy);
            v[3] = -c.A1 * xy + c.A2 * xi * g.Q2 * g.Y32;
            v[4] = -c.A1 / g.R + c.A2 * g.Q2 / g.R3;
            v[5] = -c.A1 * qy - c.A2 * q * g.Q2 * g.Y32;
            v[6] = -c.A1 * (c.Cd / g.R + qy * c.Sd) - c.A2 * q * g.Fy;
            v[7] = -c.A1 * g.Y * g.X11 - c.A2 * q * g.Gy;
            v[8] = c.A1 * (g.D * g.X11 + xy * c.Sd) + c.A2 * q * g.Hy;
            v[9] = c.A1 * (c.Sd / g.R - qy * c.Cd) - c.A2 * q * g.Fz;
            v[10] = c.A1 * g.D * g.X11 - c.A2 * q * g.Gz;
            v[11] = c.A1 * (g.Y * g.X11 + xy * c.Cd) + c.A2 * q * g.Hz;
            Add(u, v, d3 / TwoPi);
        }
        return u;
    }

    private static double[] PartB(double xi, double et, double q,
        double d1, double d2, double d3, Medium c, Geometry g)
    {
        var u = new double[12];
        var v = new double[12];
        double rd = g.R + g.D, d11 = 1 / (g.R * rd);
        double aj2 = xi * g.Y / rd * d11;
        double aj5 = -(g.D + g.Y * g.Y / rd) * d11;
        double ai3, ai4, ak1, ak3, aj3, aj6;
        if (c.Cd != 0)
        {
            if (xi == 0) ai4 = 0;
            else
            {
                double xx = Math.Sqrt(g.Xi2 + g.Q2);
                ai4 = (xi / rd * c.SdCd + 2 * Math.Atan(
                    (et * (xx + q * c.Cd) + xx * (g.R + xx) * c.Sd) /
                    (xi * (g.R + xx) * c.Cd))) / c.Cd2;
            }
            ai3 = (g.Y * c.Cd / rd - g.Ale + c.Sd * Math.Log(rd)) / c.Cd2;
            ak1 = xi * (d11 - g.Y11 * c.Sd) / c.Cd;
            ak3 = (q * g.Y11 - g.Y * d11) / c.Cd;
            aj3 = (ak1 - aj2 * c.Sd) / c.Cd;
            aj6 = (ak3 - aj5 * c.Sd) / c.Cd;
        }
        else
        {
            double rd2 = rd * rd;
            ai3 = (et / rd + g.Y * q / rd2 - g.Ale) / 2;
            ai4 = xi * g.Y / rd2 / 2;
            ak1 = xi * q / rd * d11;
            ak3 = c.Sd / rd * (g.Xi2 * d11 - 1);
            aj3 = -xi / rd2 * (g.Q2 * d11 - 0.5);
            aj6 = -g.Y / rd2 * (g.Xi2 * d11 - 0.5);
        }
        double xy = xi * g.Y11;
        double ai1 = -xi / rd * c.Cd - ai4 * c.Sd;
        double ai2 = Math.Log(rd) + ai3 * c.Sd;
        double ak2 = 1 / g.R + ak3 * c.Sd;
        double ak4 = xy * c.Cd - ak1 * c.Sd;
        double aj1 = aj5 * c.Cd - aj6 * c.Sd;
        double aj4 = -xy - aj2 * c.Cd + aj3 * c.Sd;
        double qx = q * g.X11, qy = q * g.Y11;
        if (d1 != 0)
        {
            v[0] = -xi * qy - g.Tt - c.A3 * ai1 * c.Sd;
            v[1] = -q / g.R + c.A3 * g.Y / rd * c.Sd;
            v[2] = q * qy - c.A3 * ai2 * c.Sd;
            v[3] = g.Xi2 * q * g.Y32 - c.A3 * aj1 * c.Sd;
            v[4] = xi * q / g.R3 - c.A3 * aj2 * c.Sd;
            v[5] = -xi * g.Q2 * g.Y32 - c.A3 * aj3 * c.Sd;
            v[6] = -xi * g.Fy - g.D * g.X11 + c.A3 * (xy + aj4) * c.Sd;
            v[7] = -g.Ey + c.A3 * (1 / g.R + aj5) * c.Sd;
            v[8] = q * g.Fy - c.A3 * (qy - aj6) * c.Sd;
            v[9] = -xi * g.Fz - g.Y * g.X11 + c.A3 * ak1 * c.Sd;
            v[10] = -g.Ez + c.A3 * g.Y * d11 * c.Sd;
            v[11] = q * g.Fz + c.A3 * ak2 * c.Sd;
            Add(u, v, d1 / TwoPi);
        }
        if (d2 != 0)
        {
            v[0] = -q / g.R + c.A3 * ai3 * c.SdCd;
            v[1] = -et * qx - g.Tt - c.A3 * xi / rd * c.SdCd;
            v[2] = q * qx + c.A3 * ai4 * c.SdCd;
            v[3] = xi * q / g.R3 + c.A3 * aj4 * c.SdCd;
            v[4] = et * q / g.R3 + qy + c.A3 * aj5 * c.SdCd;
            v[5] = -g.Q2 / g.R3 + c.A3 * aj6 * c.SdCd;
            v[6] = -g.Ey + c.A3 * aj1 * c.SdCd;
            v[7] = -et * g.Gy - xy * c.Sd + c.A3 * aj2 * c.SdCd;
            v[8] = q * g.Gy + c.A3 * aj3 * c.SdCd;
            v[9] = -g.Ez - c.A3 * ak3 * c.SdCd;
            v[10] = -et * g.Gz - xy * c.Cd - c.A3 * xi * d11 * c.SdCd;
            v[11] = q * g.Gz - c.A3 * ak4 * c.SdCd;
            Add(u, v, d2 / TwoPi);
        }
        if (d3 != 0)
        {
            v[0] = q * qy - c.A3 * ai3 * c.Sd2;
            v[1] = q * qx + c.A3 * xi / rd * c.Sd2;
            v[2] = et * qx + xi * qy - g.Tt - c.A3 * ai4 * c.Sd2;
            v[3] = -xi * g.Q2 * g.Y32 - c.A3 * aj4 * c.Sd2;
            v[4] = -g.Q2 / g.R3 - c.A3 * aj5 * c.Sd2;
            v[5] = q * g.Q2 * g.Y32 - c.A3 * aj6 * c.Sd2;
            v[6] = q * g.Fy - c.A3 * aj1 * c.Sd2;
            v[7] = q * g.Gy - c.A3 * aj2 * c.Sd2;
            v[8] = -q * g.Hy - c.A3 * aj3 * c.Sd2;
            v[9] = q * g.Fz + c.A3 * ak3 * c.Sd2;
            v[10] = q * g.Gz + c.A3 * xi * d11 * c.Sd2;
            v[11] = -q * g.Hz + c.A3 * ak4 * c.Sd2;
            Add(u, v, d3 / TwoPi);
        }
        return u;
    }

    private static double[] PartC(double xi, double et, double q, double z,
        double d1, double d2, double d3, Medium c, Geometry g)
    {
        var u = new double[12];
        var v = new double[12];
        double cc = g.D + z;
        double x53 = (8 * g.R2 + 9 * g.R * xi + 3 * g.Xi2) *
                     g.X11 * g.X11 * g.X11 / g.R2;
        double y53 = (8 * g.R2 + 9 * g.R * et + 3 * g.Et2) *
                     g.Y11 * g.Y11 * g.Y11 / g.R2;
        double h = q * c.Cd - z;
        double z32 = c.Sd / g.R3 - h * g.Y32;
        double z53 = 3 * c.Sd / g.R5 - h * y53;
        double y0 = g.Y11 - g.Xi2 * g.Y32;
        double z0 = z32 - g.Xi2 * z53;
        double ppy = c.Cd / g.R3 + q * g.Y32 * c.Sd;
        double ppz = c.Sd / g.R3 - q * g.Y32 * c.Cd;
        double qq = z * g.Y32 + z32 + z0;
        double qqy = 3 * cc * g.D / g.R5 - qq * c.Sd;
        double qqz = 3 * cc * g.Y / g.R5 - qq * c.Cd + q * g.Y32;
        double xy = xi * g.Y11, qy = q * g.Y11, qr = 3 * q / g.R5;
        double cdr = (cc + g.D) / g.R3;
        double yy0 = g.Y / g.R3 - y0 * c.Cd;
        if (d1 != 0)
        {
            v[0] = c.A4 * xy * c.Cd - c.A5 * xi * q * z32;
            v[1] = c.A4 * (c.Cd / g.R + 2 * qy * c.Sd) - c.A5 * cc * q / g.R3;
            v[2] = c.A4 * qy * c.Cd - c.A5 * (cc * et / g.R3 - z * g.Y11 + g.Xi2 * z32);
            v[3] = c.A4 * y0 * c.Cd - c.A5 * q * z0;
            v[4] = -c.A4 * xi * (c.Cd / g.R3 + 2 * q * g.Y32 * c.Sd) + c.A5 * cc * xi * qr;
            v[5] = -c.A4 * xi * q * g.Y32 * c.Cd + c.A5 * xi * (3 * cc * et / g.R5 - qq);
            v[6] = -c.A4 * xi * ppy * c.Cd - c.A5 * xi * qqy;
            v[7] = c.A4 * 2 * (g.D / g.R3 - y0 * c.Sd) * c.Sd - g.Y / g.R3 * c.Cd
                   - c.A5 * (cdr * c.Sd - et / g.R3 - cc * g.Y * qr);
            v[8] = -c.A4 * q / g.R3 + yy0 * c.Sd
                   + c.A5 * (cdr * c.Cd + cc * g.D * qr - (y0 * c.Cd + q * z0) * c.Sd);
            v[9] = c.A4 * xi * ppz * c.Cd - c.A5 * xi * qqz;
            v[10] = c.A4 * 2 * (g.Y / g.R3 - y0 * c.Cd) * c.Sd + g.D / g.R3 * c.Cd
                    - c.A5 * (cdr * c.Cd + cc * g.D * qr);
            v[11] = yy0 * c.Cd - c.A5 * (cdr * c.Sd - cc * g.Y * qr - y0 * c.Sd2 + q * z0 * c.Cd);
            Add(u, v, d1 / TwoPi);
        }
        if (d2 != 0)
        {
            v[0] = c.A4 * c.Cd / g.R - qy * c.Sd - c.A5 * cc * q / g.R3;
            v[1] = c.A4 * g.Y * g.X11 - c.A5 * cc * et * q * g.X32;
            v[2] = -g.D * g.X11 - xy * c.Sd - c.A5 * cc * (g.X11 - g.Q2 * g.X32);
            v[3] = -c.A4 * xi / g.R3 * c.Cd + c.A5 * cc * xi * qr + xi * q * g.Y32 * c.Sd;
            v[4] = -c.A4 * g.Y / g.R3 + c.A5 * cc * et * qr;
            v[5] = g.D / g.R3 - y0 * c.Sd + c.A5 * cc / g.R3 * (1 - 3 * g.Q2 / g.R2);
            v[6] = -c.A4 * et / g.R3 + y0 * c.Sd2 - c.A5 * (cdr * c.Sd - cc * g.Y * qr);
            v[7] = c.A4 * (g.X11 - g.Y * g.Y * g.X32)
                   - c.A5 * cc * ((g.D + 2 * q * c.Cd) * g.X32 - g.Y * et * q * x53);
            v[8] = xi * ppy * c.Sd + g.Y * g.D * g.X32
                   + c.A5 * cc * ((g.Y + 2 * q * c.Sd) * g.X32 - g.Y * g.Q2 * x53);
            v[9] = -q / g.R3 + y0 * c.SdCd - c.A5 * (cdr * c.Cd + cc * g.D * qr);
            v[10] = c.A4 * g.Y * g.D * g.X32
                    - c.A5 * cc * ((g.Y - 2 * q * c.Sd) * g.X32 + g.D * et * q * x53);
            v[11] = -xi * ppz * c.Sd + g.X11 - g.D * g.D * g.X32
                    - c.A5 * cc * ((g.D - 2 * q * c.Cd) * g.X32 - g.D * g.Q2 * x53);
            Add(u, v, d2 / TwoPi);
        }
        if (d3 != 0)
        {
            v[0] = -c.A4 * (c.Sd / g.R + qy * c.Cd) - c.A5 * (z * g.Y11 - g.Q2 * z32);
            v[1] = c.A4 * 2 * xy * c.Sd + g.D * g.X11 - c.A5 * cc * (g.X11 - g.Q2 * g.X32);
            v[2] = c.A4 * (g.Y * g.X11 + xy * c.Cd) + c.A5 * q * (cc * et * g.X32 + xi * z32);
            v[3] = c.A4 * xi / g.R3 * c.Sd + xi * q * g.Y32 * c.Cd
                   + c.A5 * xi * (3 * cc * et / g.R5 - 2 * z32 - z0);
            v[4] = c.A4 * 2 * y0 * c.Sd - g.D / g.R3 + c.A5 * cc / g.R3 * (1 - 3 * g.Q2 / g.R2);
            v[5] = -c.A4 * yy0 - c.A5 * (cc * et * qr - q * z0);
            v[6] = c.A4 * (q / g.R3 + y0 * c.SdCd)
                   + c.A5 * (z / g.R3 * c.Cd + cc * g.D * qr - q * z0 * c.Sd);
            v[7] = -c.A4 * 2 * xi * ppy * c.Sd - g.Y * g.D * g.X32
                   + c.A5 * cc * ((g.Y + 2 * q * c.Sd) * g.X32 - g.Y * g.Q2 * x53);
            v[8] = -c.A4 * (xi * ppy * c.Cd - g.X11 + g.Y * g.Y * g.X32)
                   + c.A5 * (cc * ((g.D + 2 * q * c.Cd) * g.X32 - g.Y * et * q * x53) + xi * qqy);
            v[9] = -et / g.R3 + y0 * c.Cd2
                   - c.A5 * (z / g.R3 * c.Sd - cc * g.Y * qr - y0 * c.Sd2 + q * z0 * c.Cd);
            v[10] = c.A4 * 2 * xi * ppz * c.Sd - g.X11 + g.D * g.D * g.X32
                    - c.A5 * cc * ((g.D - 2 * q * c.Cd) * g.X32 - g.D * g.Q2 * x53);
            v[11] = c.A4 * (xi * ppz * c.Cd + g.Y * g.D * g.X32)
                    + c.A5 * (cc * ((g.Y - 2 * q * c.Sd) * g.X32 + g.D * et * q * x53) + xi * qqz);
            Add(u, v, d3 / TwoPi);
        }
        return u;
    }

    private static Medium Constants(double alpha, double dipDeg)
    {
        double sd = Math.Sin(dipDeg * Math.PI / 180.0);
        double cd = Math.Cos(dipDeg * Math.PI / 180.0);
        if (Math.Abs(cd) < Eps)
        {
            cd = 0;
            sd = Math.Sign(sd);
        }
        return new Medium((1 - alpha) / 2, alpha / 2,
            (1 - alpha) / alpha, 1 - alpha, alpha,
            sd, cd, sd * sd, cd * cd, sd * cd);
    }

    private static Geometry GeometryConstants(double xi, double et, double q,
        double sd, double cd, int kxi, int ket)
    {
        xi = Snap(xi); et = Snap(et); q = Snap(q);
        double xi2 = xi * xi, et2 = et * et, q2 = q * q;
        double r2 = xi2 + et2 + q2, r = Math.Sqrt(r2);
        if (r == 0) return default;
        double r3 = r * r2, r5 = r3 * r2;
        double y = et * cd + q * sd, d = et * sd - q * cd;
        double tt = q == 0 ? 0 : Math.Atan(xi * et / (q * r));
        double alx, x11, x32;
        if (kxi == 1) { alx = -Math.Log(r - xi); x11 = x32 = 0; }
        else
        {
            double rxi = r + xi;
            alx = Math.Log(rxi); x11 = 1 / (r * rxi);
            x32 = (r + rxi) * x11 * x11 / r;
        }
        double ale, y11, y32;
        if (ket == 1) { ale = -Math.Log(r - et); y11 = y32 = 0; }
        else
        {
            double ret = r + et;
            ale = Math.Log(ret); y11 = 1 / (r * ret);
            y32 = (r + ret) * y11 * y11 / r;
        }
        double ey = sd / r - y * q / r3;
        double ez = cd / r + d * q / r3;
        double fy = d / r3 + xi2 * y32 * sd;
        double fz = y / r3 + xi2 * y32 * cd;
        double gy = 2 * x11 * sd - y * q * x32;
        double gz = 2 * x11 * cd + d * q * x32;
        double hy = d * q * x32 + xi * q * y32 * sd;
        double hz = y * q * x32 + xi * q * y32 * cd;
        return new Geometry(xi2, et2, q2, r, r2, r3, r5, y, d, tt,
            alx, ale, x11, y11, x32, y32, ey, ez, fy, fz, gy, gz, hy, hz);
    }

    private static void Add(double[] target, double[] values, double scale)
    {
        for (int i = 0; i < 12; i++) target[i] += values[i] * scale;
    }

    private static double Snap(double value) => Math.Abs(value) < Eps ? 0 : value;
    private static OkadaDc3dResult Empty(int status) => new(status,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
}
