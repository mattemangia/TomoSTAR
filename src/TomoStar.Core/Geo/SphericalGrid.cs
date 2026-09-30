// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using System.Text.Json.Serialization;

namespace TomoStar.Core.Geo;

/// <summary>
/// The serialisable definition of a regular grid in longitude, latitude and depth. Nodes sit ON
/// the bounds: node i of Nx is at <c>MinLon + i·(MaxLon − MinLon)/(Nx − 1)</c>, so the model is
/// defined everywhere inside the box and the volume a user asked for is the volume they get.
/// </summary>
public sealed record GridDefinition
{
    public double MinLon { get; init; }
    public double MaxLon { get; init; }
    public double MinLat { get; init; }
    public double MaxLat { get; init; }

    /// <summary>Top of the grid, km below sea level; negative to include topography.</summary>
    public double MinDepthKm { get; init; }

    public double MaxDepthKm { get; init; }
    public int Nx { get; init; } = 2;
    public int Ny { get; init; } = 2;
    public int Nz { get; init; } = 2;

    [JsonIgnore] public double DLon => (MaxLon - MinLon) / Math.Max(1, Nx - 1);
    [JsonIgnore] public double DLat => (MaxLat - MinLat) / Math.Max(1, Ny - 1);
    [JsonIgnore] public double DDepth => (MaxDepthKm - MinDepthKm) / Math.Max(1, Nz - 1);
    [JsonIgnore] public long Count => (long)Nx * Ny * Nz;

    public void Validate()
    {
        if (!double.IsFinite(MinLon) || !double.IsFinite(MaxLon) || !double.IsFinite(MinLat) || !double.IsFinite(MaxLat)
            || !double.IsFinite(MinDepthKm) || !double.IsFinite(MaxDepthKm))
            throw new ArgumentException("Grid bounds must be finite.");
        if (Nx < 2 || Ny < 2 || Nz < 2) throw new ArgumentException("A grid needs at least two nodes on every axis.");
        if (!(MaxLon > MinLon) || !(MaxLat > MinLat) || !(MaxDepthKm > MinDepthKm))
            throw new ArgumentException("Grid bounds are empty or inverted.");
        if (MaxLon - MinLon >= 360) throw new ArgumentException("A grid must span less than 360° of longitude.");
        if (MinLat < -89.9 || MaxLat > 89.9) throw new ArgumentException("Grids touching the poles are not supported.");
        if (MaxDepthKm >= GeoMath.EarthRadiusKm) throw new ArgumentException("The grid reaches the centre of the Earth.");
    }

    /// <summary>
    /// A grid over the same box with node spacing divided by <paramref name="factor"/>: the forward
    /// grid the eikonal is solved on, finer than the inversion grid so the travel times are not
    /// limited by the parameterisation.
    /// </summary>
    public GridDefinition Refined(int factor) => factor <= 1 ? this : this with
    {
        Nx = (Nx - 1) * factor + 1,
        Ny = (Ny - 1) * factor + 1,
        Nz = (Nz - 1) * factor + 1
    };

    /// <summary>A grid covering a box at approximately the requested node spacing in km.</summary>
    public static GridDefinition FromSpacing(
        double minLon, double maxLon, double minLat, double maxLat,
        double minDepthKm, double maxDepthKm, double horizontalKm, double verticalKm)
    {
        var midLat = 0.5 * (minLat + maxLat) * GeoMath.Deg2Rad;
        var widthKm = (maxLon - minLon) * GeoMath.Deg2Rad * GeoMath.EarthRadiusKm * Math.Cos(midLat);
        var heightKm = (maxLat - minLat) * GeoMath.Deg2Rad * GeoMath.EarthRadiusKm;
        return new GridDefinition
        {
            MinLon = minLon, MaxLon = maxLon, MinLat = minLat, MaxLat = maxLat,
            MinDepthKm = minDepthKm, MaxDepthKm = maxDepthKm,
            Nx = Math.Max(2, (int)Math.Round(widthKm / horizontalKm) + 1),
            Ny = Math.Max(2, (int)Math.Round(heightKm / horizontalKm) + 1),
            Nz = Math.Max(2, (int)Math.Round((maxDepthKm - minDepthKm) / verticalKm) + 1)
        };
    }
}

/// <summary>
/// A <see cref="GridDefinition"/> with its node coordinates precomputed: the lattice the eikonal is
/// solved on, the rays are traced through and the model is parameterised on.
///
/// Index order is x (longitude) fastest, then y (latitude), then z (depth), so a depth slice is a
/// contiguous block of the model vector and of a volume file.
/// </summary>
public sealed class SphericalGrid
{
    public SphericalGrid(GridDefinition definition)
    {
        definition.Validate();
        Definition = definition;
        Nx = definition.Nx;
        Ny = definition.Ny;
        Nz = definition.Nz;
        if ((long)Nx * Ny * Nz > int.MaxValue) throw new ArgumentException("Grid has more than 2^31 nodes.");
        Count = Nx * Ny * Nz;

        LonDeg = new double[Nx];
        for (var i = 0; i < Nx; i++) LonDeg[i] = definition.MinLon + i * definition.DLon;
        LatDeg = new double[Ny];
        CosLat = new double[Ny];
        for (var j = 0; j < Ny; j++)
        {
            LatDeg[j] = definition.MinLat + j * definition.DLat;
            CosLat[j] = Math.Cos(LatDeg[j] * GeoMath.Deg2Rad);
        }
        DepthKm = new double[Nz];
        RadiusKm = new double[Nz];
        for (var k = 0; k < Nz; k++)
        {
            DepthKm[k] = definition.MinDepthKm + k * definition.DDepth;
            RadiusKm[k] = GeoMath.EarthRadiusKm - DepthKm[k];
        }
        CentreLon = 0.5 * (definition.MinLon + definition.MaxLon);
        DLonRad = definition.DLon * GeoMath.Deg2Rad;
        DLatRad = definition.DLat * GeoMath.Deg2Rad;
        DDepthKm = definition.DDepth;
    }

    public GridDefinition Definition { get; }
    public int Nx { get; }
    public int Ny { get; }
    public int Nz { get; }
    public int Count { get; }
    public double[] LonDeg { get; }

    /// <summary>Longitude of the centre: longitudes are brought within 180° of it before indexing.</summary>
    public double CentreLon { get; }
    public double[] LatDeg { get; }
    public double[] CosLat { get; }
    public double[] DepthKm { get; }
    public double[] RadiusKm { get; }
    public double DLonRad { get; }
    public double DLatRad { get; }
    public double DDepthKm { get; }

    public int Index(int i, int j, int k) => (k * Ny + j) * Nx + i;

    public (int I, int J, int K) Decompose(int index)
    {
        var k = index / (Nx * Ny);
        var rest = index - k * Nx * Ny;
        var j = rest / Nx;
        return (rest - j * Nx, j, k);
    }

    public bool InRange(int i, int j, int k) => (uint)i < (uint)Nx && (uint)j < (uint)Ny && (uint)k < (uint)Nz;

    /// <summary>
    /// Physical node spacings at a node, km: the metric of the spherical eikonal. One longitude
    /// step is r·cosφ·Δλ long, so it shrinks with depth and towards the poles; treating the
    /// spacings as constants is how a "spherical" solver quietly becomes a Cartesian one.
    /// </summary>
    public (double HLon, double HLat, double HR) Spacing(int i, int j, int k)
    {
        var r = RadiusKm[k];
        return (Math.Max(1e-9, r * CosLat[j] * DLonRad), r * DLatRad, DDepthKm);
    }

    /// <summary>Smallest physical spacing anywhere on the grid, km.</summary>
    public double MinSpacingKm()
    {
        var rMin = RadiusKm[Nz - 1];
        var cosMin = CosLat.Min();
        return Math.Min(DDepthKm, Math.Min(rMin * DLatRad, rMin * cosMin * DLonRad));
    }

    public Vec3 Cartesian(int i, int j, int k) => GeoMath.ToCartesian(LonDeg[i], LatDeg[j], DepthKm[k]);

    public Vec3 Cartesian(int index)
    {
        var (i, j, k) = Decompose(index);
        return Cartesian(i, j, k);
    }

    /// <summary>
    /// Fractional node coordinates of a geographic point (unclamped). The longitude may be given in
    /// any turn (−178° and 182° are the same place), so a grid across the antimeridian works with
    /// positions that come back from Cartesian coordinates in [−180°, 180°].
    /// </summary>
    public (double Fx, double Fy, double Fz) Fractional(double lon, double lat, double depthKm)
    {
        var d = Definition;
        if (lon < CentreLon - 180 || lon > CentreLon + 180) lon = GeoMath.UnwrapLon(lon, CentreLon);
        return ((lon - d.MinLon) / d.DLon, (lat - d.MinLat) / d.DLat, (depthKm - d.MinDepthKm) / d.DDepth);
    }

    public bool Contains(double lon, double lat, double depthKm, double tolerance = 1e-9)
    {
        var (fx, fy, fz) = Fractional(lon, lat, depthKm);
        return fx >= -tolerance && fx <= Nx - 1 + tolerance
            && fy >= -tolerance && fy <= Ny - 1 + tolerance
            && fz >= -tolerance && fz <= Nz - 1 + tolerance;
    }

    /// <summary>
    /// The eight trilinear weights of a point, clamped to the grid. Returns the number of entries
    /// written (always 8; weights may be zero). This is the model's interpolation AND the kernel
    /// of the inversion: a ray segment of length l at this point contributes l·wₙ to node n.
    /// </summary>
    public int TrilinearWeights(double lon, double lat, double depthKm, Span<int> nodes, Span<double> weights)
    {
        var (fx, fy, fz) = Fractional(lon, lat, depthKm);
        return TrilinearWeightsFractional(fx, fy, fz, nodes, weights);
    }

    public int TrilinearWeightsFractional(double fx, double fy, double fz, Span<int> nodes, Span<double> weights)
    {
        fx = Math.Clamp(fx, 0, Nx - 1);
        fy = Math.Clamp(fy, 0, Ny - 1);
        fz = Math.Clamp(fz, 0, Nz - 1);
        var i0 = Math.Min((int)fx, Nx - 2);
        var j0 = Math.Min((int)fy, Ny - 2);
        var k0 = Math.Min((int)fz, Nz - 2);
        double tx = fx - i0, ty = fy - j0, tz = fz - k0;
        var n = 0;
        for (var dk = 0; dk < 2; dk++)
        for (var dj = 0; dj < 2; dj++)
        for (var di = 0; di < 2; di++)
        {
            nodes[n] = Index(i0 + di, j0 + dj, k0 + dk);
            weights[n] = (di == 0 ? 1 - tx : tx) * (dj == 0 ? 1 - ty : ty) * (dk == 0 ? 1 - tz : tz);
            n++;
        }
        return n;
    }

    /// <summary>Trilinear interpolation of a node field at a geographic point (clamped).</summary>
    public double Interpolate(ReadOnlySpan<double> field, double lon, double lat, double depthKm)
    {
        Span<int> nodes = stackalloc int[8];
        Span<double> w = stackalloc double[8];
        TrilinearWeights(lon, lat, depthKm, nodes, w);
        double sum = 0;
        for (var n = 0; n < 8; n++) sum += field[nodes[n]] * w[n];
        return sum;
    }

    /// <summary>
    /// Resamples a node field defined on <paramref name="source"/> onto this grid by trilinear
    /// interpolation: inversion grid → forward grid.
    /// </summary>
    public double[] ResampleFrom(SphericalGrid source, ReadOnlySpan<double> field)
    {
        var result = new double[Count];
        var src = field.ToArray();
        Parallel.For(0, Nz, k =>
        {
            for (var j = 0; j < Ny; j++)
            for (var i = 0; i < Nx; i++)
                result[Index(i, j, k)] = source.Interpolate(src, LonDeg[i], LatDeg[j], DepthKm[k]);
        });
        return result;
    }

    /// <summary>
    /// Local display frame centred on the grid: east, north, up at the centre of the top face. A
    /// deep grid drawn in this frame shows its true curvature (a truncated cone), a shallow one
    /// looks like a box because it nearly is one.
    /// </summary>
    public LocalFrame DisplayFrame() => new(
        0.5 * (Definition.MinLon + Definition.MaxLon),
        0.5 * (Definition.MinLat + Definition.MaxLat));
}

/// <summary>
/// A tangent-plane frame anchored at a surface point; converts geographic positions to local
/// (east, north, up) km in true spherical geometry, with an optional vertical exaggeration applied
/// to depth below/above the sphere (so the curvature itself is never exaggerated).
/// </summary>
public sealed class LocalFrame
{
    private readonly Vec3 _origin;
    private readonly Vec3 _east, _north, _up;

    public LocalFrame(double lon, double lat)
    {
        Lon = lon;
        Lat = lat;
        _origin = GeoMath.ToCartesian(lon, lat, 0);
        (_east, _north, _up) = GeoMath.LocalFrame(lon, lat);
    }

    public double Lon { get; }
    public double Lat { get; }

    /// <summary>Local (east, north, up) km of a geographic point, depth scaled by <paramref name="verticalExaggeration"/>.</summary>
    public Vec3 ToLocal(double lon, double lat, double depthKm, double verticalExaggeration = 1.0)
    {
        var p = GeoMath.ToCartesian(lon, lat, depthKm * verticalExaggeration) - _origin;
        return new Vec3(p.Dot(_east), p.Dot(_north), p.Dot(_up));
    }

    public Vec3 ToLocal(Vec3 cartesian)
    {
        var p = cartesian - _origin;
        return new Vec3(p.Dot(_east), p.Dot(_north), p.Dot(_up));
    }

    /// <summary>Inverse of <see cref="ToLocal(double,double,double,double)"/>: local (east, north, up) km to geography.</summary>
    public GeoPoint ToGeo(Vec3 local, double verticalExaggeration = 1.0)
    {
        var c = _origin + _east * local.X + _north * local.Y + _up * local.Z;
        var g = GeoMath.FromCartesian(c);
        return g with { DepthKm = g.DepthKm / verticalExaggeration };
    }

    /// <summary>Rotates a geocentric direction vector into the local frame.</summary>
    public Vec3 DirectionToLocal(Vec3 d) => new(d.Dot(_east), d.Dot(_north), d.Dot(_up));
}
