// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using TomoStar.Core.Attenuation;
using TomoStar.Core.Forward;
using TomoStar.Core.Geo;
using TomoStar.Core.Model;
using TomoStar.Core.Tomography;

namespace TomoStar.Tests;

/// <summary>
/// The Q tomography's regularisation scaled by the sensitivity of the data: a high background Q and
/// realistic t* errors must not flatten the model onto the reference Q, and the station terms must
/// neither take up shallow structure that spans several stations nor miss genuine site effects.
/// </summary>
public class QShallowLayerTests
{
    private static readonly GridDefinition Grid = new()
    {
        MinLon = 13.0, MaxLon = 13.6, MinLat = 42.5, MaxLat = 43.0, MinDepthKm = -1, MaxDepthKm = 23, Nx = 13, Ny = 12, Nz = 9
    };

    /// <summary>
    /// Q 700 (Vp 6 km/s) under a slow 3 km layer whose western half has Q 100, 36 stations, 60 events
    /// at 4–21 km, 8 ms t* errors, and optionally a random t* offset at each station (site effect).
    /// </summary>
    private static (SphericalGrid Grid, double[] Vp, double[] Vs, double[] QTrue, double[] Site, ObservationSet Data) Layered(bool shallowLayer, double siteSpread)
    {
        var grid = new SphericalGrid(Grid);
        var vp = new double[grid.Count];
        var qTrue = new double[grid.Count];
        for (var i = 0; i < grid.Count; i++)
        {
            var (x, _, z) = grid.Decompose(i);
            var shallow = grid.DepthKm[z] < 3;
            vp[i] = shallow ? 3.5 : 6.0;
            qTrue[i] = shallowLayer && shallow && x < grid.Nx / 2 ? 1 / 100.0 : 1 / 700.0;
        }
        var vs = vp.Select(v => v / 1.73).ToArray();
        var data = new ObservationSet();
        var n = 0;
        for (var i = 0; i < 6; i++)
        for (var j = 0; j < 6; j++)
            data.Stations.Add(new StationState { Id = $"XX.S{n++:00}", Lon = 13.03 + 0.105 * i, Lat = 42.53 + 0.088 * j, DepthKm = -0.3 });
        var rnd = new Random(3);
        double G() => Math.Sqrt(-2 * Math.Log(1 - rnd.NextDouble())) * Math.Cos(2 * Math.PI * rnd.NextDouble());
        var site = data.Stations.Select(_ => siteSpread * G()).ToArray();
        var field = qTrue.Select((q, i) => q / vp[i]).ToArray();
        for (var e = 0; e < 60; e++)
        {
            var ev = new EventState { Id = $"e{e}", Lon = 13.03 + 0.54 * rnd.NextDouble(), Lat = 42.53 + 0.44 * rnd.NextDouble(), DepthKm = 4 + 17 * rnd.NextDouble(), Fixed = true };
            data.Events.Add(ev);
            for (var s = 0; s < data.Stations.Count; s++)
            {
                var st = data.Stations[s];
                var t = RayTracing.Time(RayTracing.Kernel(grid, RayTracing.Straight(ev.Lon, ev.Lat, ev.DepthKm, st.Lon, st.Lat, st.DepthKm, 0.5)), field);
                data.Arrivals.Add(new Observation { Event = e, Station = s, Phase = Phase.P, Time = t + site[s] + 0.008 * G(), Sigma = 0.008 });
            }
        }
        return (grid, vp, vs, qTrue, site, data);
    }

    private static QTomographyResult Run(SphericalGrid grid, ObservationSet data, double[] vp, double[] vs)
    {
        // Default regularisation: what is tested is the scaling, not tuned values.
        var settings = new QTomographySettings { RayMethod = RayMethod.Straight, UseOpenCl = false, StationTerms = true };
        var dir = Path.Combine(Path.GetTempPath(), $"quiver_{Guid.NewGuid():N}");
        try { return new QTomography(grid, settings, dir).Run(data, vp, vs); }
        finally { if (Directory.Exists(dir)) Directory.Delete(dir, true); }
    }

    /// <summary>
    /// The low-Q half of the shallow layer is imaged in the volume (it was flattened onto the
    /// reference Q, the contrast taken up by the station terms, when the regularisation was absolute).
    /// </summary>
    [Fact]
    public void ShallowLowQLayerIsImagedNotAbsorbedByStationTerms()
    {
        var (grid, vp, vs, qTrue, _, data) = Layered(shallowLayer: true, siteSpread: 0);
        var r = Run(grid, data, vp, vs);
        var surface = Enumerable.Range(0, grid.Count).Where(i => r.Dws[i] > 5 && grid.Decompose(i).Item3 <= 1).ToList();
        var low = Median(surface.Where(i => qTrue[i] > 1 / 200.0).Select(i => r.Q[i]));
        var high = Median(surface.Where(i => qTrue[i] < 1 / 200.0).Select(i => r.Q[i]));
        Assert.InRange(low, 70, 150);
        Assert.True(high > 3 * low, $"shallow Q {low:0} under the low-Q half, {high:0} under the other (reference {r.ReferenceQ:0})");
        // The ~9 ms the layer adds to the western stations' t* is in the model, not in their terms.
        Assert.True(r.StationTerms.Average(Math.Abs) < 0.003, $"mean |station term| {1000 * r.StationTerms.Average(Math.Abs):0.0} ms");
    }

    /// <summary>Uniform Q with a random t* offset at each station: the station terms recover the offsets.</summary>
    [Fact]
    public void StationTermsStillRecoverSiteEffects()
    {
        var (grid, vp, vs, _, site, data) = Layered(shallowLayer: false, siteSpread: 0.006);
        var r = Run(grid, data, vp, vs);
        var centred = site.Select(x => x - site.Average()).ToArray();
        var corr = Correlation(centred, r.StationTerms);
        Assert.True(corr > 0.9, $"station terms vs site offsets: correlation {corr:0.00}");
        Assert.True(r.StationTerms.Average(Math.Abs) > 0.5 * centred.Average(Math.Abs),
            $"mean |station term| {1000 * r.StationTerms.Average(Math.Abs):0.0} ms for {1000 * centred.Average(Math.Abs):0.0} ms of site offsets");
    }

    private static double Correlation(double[] x, double[] y)
    {
        double mx = x.Average(), my = y.Average();
        var d = Math.Sqrt(x.Sum(v => (v - mx) * (v - mx)) * y.Sum(v => (v - my) * (v - my)));
        return d == 0 ? 0 : x.Zip(y).Sum(p => (p.First - mx) * (p.Second - my)) / d;
    }

    private static double Median(IEnumerable<double> v)
    {
        var a = v.OrderBy(x => x).ToArray();
        return a.Length == 0 ? double.NaN : a[a.Length / 2];
    }
}
