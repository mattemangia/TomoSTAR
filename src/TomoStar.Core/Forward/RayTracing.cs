// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using TomoStar.Core.Geo;

namespace TomoStar.Core.Forward;

/// <summary>A ray path from source to receiver, as geographic points.</summary>
public sealed class RayPath
{
    public required double[] Lon { get; init; }
    public required double[] Lat { get; init; }
    public required double[] Depth { get; init; }
    public int Count => Lon.Length;

    public double LengthKm()
    {
        double l = 0;
        for (var i = 1; i < Count; i++)
            l += (GeoMath.ToCartesian(Lon[i], Lat[i], Depth[i]) - GeoMath.ToCartesian(Lon[i - 1], Lat[i - 1], Depth[i - 1])).Length;
        return l;
    }

    /// <summary>A copy with at most <paramref name="max"/> points, for storage and display.</summary>
    public RayPath Decimated(int max)
    {
        if (Count <= max) return this;
        var idx = Enumerable.Range(0, max).Select(i => (int)Math.Round(i * (Count - 1.0) / (max - 1))).ToArray();
        return new RayPath { Lon = idx.Select(i => Lon[i]).ToArray(), Lat = idx.Select(i => Lat[i]).ToArray(), Depth = idx.Select(i => Depth[i]).ToArray() };
    }
}

public enum RayMethod
{
    /// <summary>Straight chords in geocentric space: fast, exact in a homogeneous Earth only.</summary>
    Straight,

    /// <summary>Curved rays traced back along the gradient of the fast-marching field.</summary>
    FastMarching
}

/// <summary>
/// Ray geometry and the travel-time kernel of a ray.
///
/// Curved rays are obtained by steepest descent on the travel-time field of the receiver: starting
/// at the source, step against ∇T until the receiver is reached. Because T is the first-arrival
/// field, the path followed is the first-arrival ray (the characteristic of the eikonal; e.g.
/// Rawlinson &amp; Sambridge 2004, §2.3; Podvin &amp; Lecomte 1991, Geophys. J. Int. 105, 271-284).
/// </summary>
public static class RayTracing
{
    public static RayPath Straight(double lon1, double lat1, double dep1, double lon2, double lat2, double dep2, double stepKm)
    {
        var a = GeoMath.ToCartesian(lon1, lat1, dep1);
        var b = GeoMath.ToCartesian(lon2, lat2, dep2);
        var len = (b - a).Length;
        var n = Math.Max(2, (int)Math.Ceiling(len / Math.Max(1e-6, stepKm)) + 1);
        var lon = new double[n];
        var lat = new double[n];
        var dep = new double[n];
        for (var i = 0; i < n; i++)
        {
            var p = GeoMath.FromCartesian(a + (b - a) * ((double)i / (n - 1)));
            lon[i] = p.Lon; lat[i] = p.Lat; dep[i] = p.DepthKm;
        }
        return new RayPath { Lon = lon, Lat = lat, Depth = dep };
    }

    /// <summary>
    /// Traces the ray from a source back to the station whose table is given. Returns null if the
    /// descent stalls (a flat or non-finite region of the field).
    /// </summary>
    public static RayPath? Backtrack(TravelTimeTable table, double lon, double lat, double depthKm, double stepKm)
    {
        var g = table.Grid;
        var h = table.Header;
        var station = GeoMath.ToCartesian(h.Lon, h.Lat, Math.Clamp(h.DepthKm, g.DepthKm[0], g.DepthKm[^1]));
        var p = GeoMath.ToCartesian(lon, lat, depthKm);
        var lons = new List<double> { lon };
        var lats = new List<double> { lat };
        var deps = new List<double> { depthKm };
        var t0 = table.Time(lon, lat, depthKm);
        if (!double.IsFinite(t0)) return null;
        var maxSteps = 20 * (int)Math.Ceiling((p - station).Length / stepKm + 10);
        var lastT = t0;
        var stalls = 0;
        for (var step = 0; step < maxSteps; step++)
        {
            if ((p - station).Length <= 1.5 * stepKm) break;
            var gp = GeoMath.FromCartesian(p);
            var (ge, gn, gz) = table.Gradient(gp.Lon, gp.Lat, gp.DepthKm);
            if (!double.IsFinite(ge) || !double.IsFinite(gn) || !double.IsFinite(gz)) return null;
            var (e, nn, up) = GeoMath.LocalFrame(gp.Lon, gp.Lat);
            var grad = e * ge + nn * gn + up * (-gz);
            var norm = grad.Length;
            if (norm < 1e-12) return null;
            p -= grad * (stepKm / norm);
            var q = GeoMath.FromCartesian(p);
            // Keep the path inside the grid (a ray grazing the side of the model).
            q = new GeoPoint(Math.Clamp(GeoMath.UnwrapLon(q.Lon, g.CentreLon), g.LonDeg[0], g.LonDeg[^1]), Math.Clamp(q.Lat, g.LatDeg[0], g.LatDeg[^1]),
                Math.Clamp(q.DepthKm, g.DepthKm[0], g.DepthKm[^1]));
            p = GeoMath.ToCartesian(q);
            lons.Add(q.Lon); lats.Add(q.Lat); deps.Add(q.DepthKm);
            var t = table.Time(q.Lon, q.Lat, q.DepthKm);
            if (!(t < lastT)) { if (++stalls > 20) return null; }
            else stalls = 0;
            lastT = Math.Min(lastT, t);
        }
        lons.Add(h.Lon); lats.Add(h.Lat); deps.Add(Math.Clamp(h.DepthKm, g.DepthKm[0], g.DepthKm[^1]));
        return new RayPath { Lon = lons.ToArray(), Lat = lats.ToArray(), Depth = deps.ToArray() };
    }

    /// <summary>
    /// Integrates a ray against the trilinear basis of a grid: kernel[n] = ∫ wₙ(x) ds, so that the
    /// travel time is Σ kernel[n]·s[n] for node slownesses s. Segments are subdivided so no piece is
    /// longer than half the node spacing. Returns the accumulated (node, length) pairs.
    /// </summary>
    public static Dictionary<int, double> Kernel(SphericalGrid grid, RayPath ray, double[]? weightByNode = null) =>
        Kernel(grid, ray, weightByNode, null, out _);

    /// <summary>
    /// Kernel of the part of a ray inside the grid. Pieces outside are not attributed to the edge
    /// nodes (as clamping would): their length times <paramref name="outsideSlowness"/>(depth) is
    /// summed in <paramref name="outsideTime"/>, the fixed time spent outside the model.
    /// </summary>
    public static Dictionary<int, double> Kernel(SphericalGrid grid, RayPath ray, double[]? weightByNode, Func<double, double>? outsideSlowness, out double outsideTime)
    {
        outsideTime = 0;
        var kernel = new Dictionary<int, double>();
        var h = 0.5 * grid.MinSpacingKm();
        Span<int> nodes = stackalloc int[8];
        Span<double> w = stackalloc double[8];
        for (var i = 1; i < ray.Count; i++)
        {
            var a = GeoMath.ToCartesian(ray.Lon[i - 1], ray.Lat[i - 1], ray.Depth[i - 1]);
            var b = GeoMath.ToCartesian(ray.Lon[i], ray.Lat[i], ray.Depth[i]);
            var len = (b - a).Length;
            if (len <= 0) continue;
            var pieces = Math.Max(1, (int)Math.Ceiling(len / h));
            var dl = len / pieces;
            for (var q = 0; q < pieces; q++)
            {
                var mid = GeoMath.FromCartesian(a + (b - a) * ((q + 0.5) / pieces));
                if (outsideSlowness != null && !grid.Contains(mid.Lon, mid.Lat, mid.DepthKm, 1e-6))
                {
                    outsideTime += dl * outsideSlowness(mid.DepthKm);
                    continue;
                }
                grid.TrilinearWeights(mid.Lon, mid.Lat, mid.DepthKm, nodes, w);
                for (var n = 0; n < 8; n++)
                {
                    if (w[n] == 0) continue;
                    var v = dl * w[n] * (weightByNode?[nodes[n]] ?? 1.0);
                    kernel[nodes[n]] = kernel.TryGetValue(nodes[n], out var old) ? old + v : v;
                }
            }
        }
        return kernel;
    }

    /// <summary>Travel time along a ray through a node slowness field: Σ kernel·s.</summary>
    public static double Time(Dictionary<int, double> kernel, double[] slowness)
    {
        double t = 0;
        foreach (var (n, l) in kernel) t += l * slowness[n];
        return t;
    }

    /// <summary>
    /// Take-off gradient of a straight ray at the source in the local (east, north, down) frame:
    /// −s·û where û is the unit vector from source to receiver.
    /// </summary>
    public static (double E, double N, double Z) StraightGradient(double lon, double lat, double depth, double rLon, double rLat, double rDepth, double slowness)
    {
        var src = GeoMath.ToCartesian(lon, lat, depth);
        var rec = GeoMath.ToCartesian(rLon, rLat, rDepth);
        var d = (rec - src).Normalized();
        var (e, n, up) = GeoMath.LocalFrame(lon, lat);
        // dT/dx_src = −s · (direction towards the receiver); "down" is −up.
        return (-slowness * d.Dot(e), -slowness * d.Dot(n), slowness * d.Dot(up));
    }
}
