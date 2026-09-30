using TomoStar.Core.Geo;
using TomoStar.Core.Model;

namespace TomoStar.Core.Forward;

/// <summary>
/// Where the forward problem is solved, which can be larger than the inversion grid.
///
/// Stations and events outside the study area are useful when their rays cross it (a regional
/// event recorded by the local network, a distant station of another network): the part of the ray
/// inside the area constrains the model there. The eikonal is then solved on a forward grid that
/// holds the inversion grid and every such station and event; inside the inversion grid the model is
/// the one being inverted, outside it is a fixed 1-D background. The inversion kernel is taken on the
/// inside part only, and the time spent outside enters the prediction but has no model parameter:
/// the classic treatment of paths leaving a regional model (e.g. Evans &amp; Achauer 1993, in Seismic
/// Tomography: Theory and Practice, 319-360, for teleseismic paths; Eberhart-Phillips 1990,
/// J. Geophys. Res. 95(B10), 15343-15363, for local events outside the model).
/// Errors of the background along the outside part are absorbed by the origin times of outside
/// events and by the station corrections, not by the model inside.
/// </summary>
public sealed class ForwardDomain
{
    public ForwardDomain(GridDefinition forward, VelocityModel1D? background = null, bool extended = false)
    {
        Forward = forward;
        Background = background;
        Extended = extended;
    }

    /// <summary>The grid the eikonal is solved on.</summary>
    public GridDefinition Forward { get; }

    /// <summary>Velocities outside the inversion grid; null when the forward grid equals the inversion box.</summary>
    public VelocityModel1D? Background { get; }

    /// <summary>The forward grid reaches beyond the inversion grid.</summary>
    public bool Extended { get; }

    /// <summary>The classic domain: the inversion box, refined.</summary>
    public static ForwardDomain Refined(GridDefinition inversion, int refinement) => new(inversion.Refined(refinement));

    /// <summary>
    /// The domain holding the inversion grid and every point given (stations and events), with a
    /// margin of two forward cells. Node spacing is the refined inversion spacing unless that would
    /// exceed <paramref name="maxNodes"/>, in which case the whole forward grid is coarsened
    /// uniformly (the log should say so: the accuracy of the times drops with it).
    /// </summary>
    public static ForwardDomain Covering(GridDefinition inversion, int refinement, IEnumerable<(double Lon, double Lat, double DepthKm)> points,
        VelocityModel1D background, long maxNodes = 12_000_000, Action<string>? log = null)
    {
        var fine = inversion.Refined(Math.Max(1, refinement));
        var centre = 0.5 * (inversion.MinLon + inversion.MaxLon);
        double minLon = inversion.MinLon, maxLon = inversion.MaxLon, minLat = inversion.MinLat, maxLat = inversion.MaxLat;
        double top = inversion.MinDepthKm, bottom = inversion.MaxDepthKm;
        var outside = 0;
        foreach (var (lon0, lat, dep) in points)
        {
            var lon = GeoMath.UnwrapLon(lon0, centre);
            var inside = lon >= inversion.MinLon && lon <= inversion.MaxLon && lat >= inversion.MinLat && lat <= inversion.MaxLat
                         && dep >= inversion.MinDepthKm && dep <= inversion.MaxDepthKm;
            if (inside) continue;
            outside++;
            minLon = Math.Min(minLon, lon); maxLon = Math.Max(maxLon, lon);
            minLat = Math.Min(minLat, lat); maxLat = Math.Max(maxLat, lat);
            // Stations above the grid top are placed on it (as inside), so only depth below extends.
            bottom = Math.Max(bottom, dep);
        }
        if (outside == 0) return new ForwardDomain(fine);

        // Margin of two forward cells so points are not on the boundary.
        minLon -= 2 * fine.DLon; maxLon += 2 * fine.DLon;
        minLat = Math.Max(-89, minLat - 2 * fine.DLat); maxLat = Math.Min(89, maxLat + 2 * fine.DLat);
        if (bottom > inversion.MaxDepthKm) bottom += 2 * fine.DDepth;
        if (maxLon - minLon >= 359) throw new InvalidOperationException("Stations and events outside the area span the whole globe: restrict the data.");

        double dLon = fine.DLon, dLat = fine.DLat, dDep = fine.DDepth;
        long Count(double f) =>
            (long)(Math.Ceiling((maxLon - minLon) / (dLon * f)) + 1) * (long)(Math.Ceiling((maxLat - minLat) / (dLat * f)) + 1)
            * (long)(Math.Ceiling((bottom - top) / (dDep * f)) + 1);
        var factor = 1.0;
        while (Count(factor) > maxNodes) factor *= 1.25;
        if (factor > 1)
            log?.Invoke($"Forward grid enlarged to hold {outside} stations/events outside the area; spacing coarsened ×{factor:0.##} to stay within {maxNodes / 1e6:0.#} M nodes.");
        else
            log?.Invoke($"Forward grid enlarged to hold {outside} stations/events outside the area (fixed 1-D model '{background.Name}' outside the inversion grid).");
        dLon *= factor; dLat *= factor; dDep *= factor;

        // Extend outward from the inversion box so its bounds stay on forward nodes when possible.
        int west = (int)Math.Ceiling((inversion.MinLon - minLon) / dLon), east = (int)Math.Ceiling((maxLon - inversion.MaxLon) / dLon);
        int south = (int)Math.Ceiling((inversion.MinLat - minLat) / dLat), north = (int)Math.Ceiling((maxLat - inversion.MaxLat) / dLat);
        int down = (int)Math.Ceiling((bottom - inversion.MaxDepthKm) / dDep);
        var nx = (int)Math.Round((inversion.MaxLon - inversion.MinLon) / dLon);
        var ny = (int)Math.Round((inversion.MaxLat - inversion.MinLat) / dLat);
        var nz = (int)Math.Round((inversion.MaxDepthKm - inversion.MinDepthKm) / dDep);
        var def = new GridDefinition
        {
            MinLon = inversion.MinLon - west * dLon, MaxLon = inversion.MinLon - west * dLon + (west + nx + east) * dLon,
            MinLat = inversion.MinLat - south * dLat, MaxLat = inversion.MinLat - south * dLat + (south + ny + north) * dLat,
            MinDepthKm = inversion.MinDepthKm, MaxDepthKm = inversion.MinDepthKm + (nz + down) * dDep,
            Nx = west + nx + east + 1, Ny = south + ny + north + 1, Nz = nz + down + 1
        };
        def.Validate();
        return new ForwardDomain(def, background, extended: true);
    }

    /// <summary>
    /// The domain an inversion of <paramref name="points"/> needs: the refined inversion grid when
    /// every point is inside it, otherwise a covering grid with <paramref name="background"/> (or
    /// the layer mean of vp/vs) outside.
    /// </summary>
    public static ForwardDomain For(SphericalGrid grid, int refinement, bool hasOutside, IEnumerable<(double Lon, double Lat, double DepthKm)> points,
        VelocityModel1D? background, double[] vp, double[] vs, Action<string>? log = null) =>
        hasOutside
            ? Covering(grid.Definition, refinement, points, background ?? LayerMean(grid, vp, vs), log: log)
            : Refined(grid.Definition, refinement);

    /// <summary>
    /// A 1-D background from a model on the inversion grid: the mean of each depth level, constant
    /// below the bottom. Used when no 1-D model is given (synthetic tests, a 3-D starting model).
    /// </summary>
    public static VelocityModel1D LayerMean(SphericalGrid grid, double[] vp, double[] vs)
    {
        var nodes = new List<VelocityNode>();
        var layer = grid.Nx * grid.Ny;
        for (var k = 0; k < grid.Nz; k++)
        {
            double p = 0, s = 0;
            for (var q = 0; q < layer; q++) { p += vp[k * layer + q]; s += vs[k * layer + q]; }
            nodes.Add(new VelocityNode(grid.DepthKm[k], p / layer, s / layer));
        }
        nodes.Add(new VelocityNode(GeoMath.EarthRadiusKm - 1, nodes[^1].Vp, nodes[^1].Vs));
        return new VelocityModel1D { Name = "layer mean of the starting model", Nodes = nodes };
    }
}
