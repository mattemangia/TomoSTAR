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

    /// <summary>
    /// Q 200 with a low-Q sphere (Q 80, radius 9 km) under the 36 stations, 60 events, t* from straight rays with 2 ms
    /// errors and 8 ms more noise added (σ 8 ms), and weak regularisation: the linear unknowns (fractional change of 1/Q)
    /// let the solve take 1/Q through zero at many nodes, which the bounds turn into Q = QMax next to low values. Solved
    /// for ln(1/Q) by Gauss–Newton, 1/Q stays positive, far fewer nodes sit on a bound and the model is closer to the truth.
    /// </summary>
    [Fact]
    public void LogParameterisationKeepsQOffTheBounds()
    {
        var grid = new SphericalGrid(Grid);
        var vp = Enumerable.Repeat(6.0, grid.Count).ToArray();
        var vs = Enumerable.Repeat(3.46, grid.Count).ToArray();
        var centre = GeoMath.ToCartesian(13.3, 42.75, 10);
        var qTrue = Enumerable.Range(0, grid.Count).Select(i =>
        {
            var (x, y, z) = grid.Decompose(i);
            return (grid.Cartesian(x, y, z) - centre).Length < 9 ? 1 / 80.0 : 1 / 200.0;
        }).ToArray();
        // The same noisy data for both solves (each run marks rejections on its own copy).
        ObservationSet Noisy()
        {
            var data = new ObservationSet();
            var n = 0;
            for (var i = 0; i < 6; i++)
            for (var j = 0; j < 6; j++)
                data.Stations.Add(new StationState { Id = $"XX.S{n++:00}", Lon = 13.03 + 0.105 * i, Lat = 42.53 + 0.088 * j, DepthKm = -0.3 });
            var rnd = new Random(8);
            double G() => Math.Sqrt(-2 * Math.Log(1 - rnd.NextDouble())) * Math.Cos(2 * Math.PI * rnd.NextDouble());
            var field = qTrue.Select((q, i) => q / vp[i]).ToArray();
            for (var e = 0; e < 60; e++)
            {
                var ev = new EventState { Id = $"e{e}", Lon = 13.03 + 0.54 * rnd.NextDouble(), Lat = 42.53 + 0.44 * rnd.NextDouble(), DepthKm = 2 + 19 * rnd.NextDouble(), Fixed = true };
                data.Events.Add(ev);
                for (var s = 0; s < data.Stations.Count; s++)
                {
                    var st = data.Stations[s];
                    var t = RayTracing.Time(RayTracing.Kernel(grid, RayTracing.Straight(ev.Lon, ev.Lat, ev.DepthKm, st.Lon, st.Lat, st.DepthKm, 0.5)), field);
                    data.Arrivals.Add(new Observation { Event = e, Station = s, Phase = Phase.P, Time = t + 0.002 * G(), Sigma = 0.002 });
                }
            }
            var noise = new Random(21);
            double N() => Math.Sqrt(-2 * Math.Log(1 - noise.NextDouble())) * Math.Cos(2 * Math.PI * noise.NextDouble());
            foreach (var a in data.Arrivals) { a.Time += 0.008 * N(); a.Sigma = 0.008; }
            return data;
        }
        var dir = Path.Combine(Path.GetTempPath(), $"quiver_{Guid.NewGuid():N}");
        try
        {
            QTomographyResult Solve(bool log) => new QTomography(grid, new QTomographySettings
            {
                RayMethod = RayMethod.Straight, UseOpenCl = false, Damping = 0.1, Smoothing = 0.3, StationTerms = true, LogParameterisation = log
            }, dir).Run(Noisy(), vp, vs);
            var linear = Solve(false);
            var logarithmic = Solve(true);
            var bounds = new QTomographySettings();
            var sampled = Enumerable.Range(0, grid.Count).Where(i => linear.Dws[i] > 0).ToList();
            int AtBounds(QTomographyResult r) => sampled.Count(i => r.Q[i] >= 0.999 * bounds.QMax || r.Q[i] <= 1.001 * bounds.QMin);
            double Err(QTomographyResult r) => sampled.Average(i => Math.Abs(Math.Log(qTrue[i] * r.Q[i])));
            var summary = $"nodes at a Q bound: linear {AtBounds(linear)}, logarithmic {AtBounds(logarithmic)} of {sampled.Count}; " +
                          $"mean |ln Q/Q_true|: linear {Err(linear):0.000}, logarithmic {Err(logarithmic):0.000}";
            Assert.True(AtBounds(linear) > 0, "the case must provoke the linear failure; " + summary);
            Assert.True(5 * AtBounds(logarithmic) < AtBounds(linear), summary);
            Assert.True(Err(logarithmic) < Err(linear), summary);
        }
        finally
        {
            if (Directory.Exists(dir)) Directory.Delete(dir, true);
        }
    }

    private static double Correlation(double[] x, double[] y)
    {
        double mx = x.Average(), my = y.Average();
        var d = Math.Sqrt(x.Sum(v => (v - mx) * (v - mx)) * y.Sum(v => (v - my) * (v - my)));
        return d == 0 ? 0 : x.Zip(y).Sum(p => (p.First - mx) * (p.Second - my)) / d;
    }

    /// <summary>
    /// The reference Q from the data: without station terms the median of the path averages T / t*,
    /// with station terms the uniform Q fitting all t* by weighted least squares, both computed here
    /// independently along the straight rays.
    /// </summary>
    [Fact]
    public void ReferenceQIsTheMedianPathAverageWithoutStationTerms()
    {
        var (grid, vp, vs, _, _, data) = Layered(shallowLayer: true, siteSpread: 0);
        var slow = vp.Select(v => 1 / v).ToArray();
        var ratios = new List<double>();
        double num = 0, den = 0;
        foreach (var a in data.Arrivals)
        {
            var e = data.Events[a.Event];
            var st = data.Stations[a.Station];
            var tt = RayTracing.Time(RayTracing.Kernel(grid, RayTracing.Straight(e.Lon, e.Lat, e.DepthKm, st.Lon, st.Lat, st.DepthKm, 0.5)), slow);
            ratios.Add(a.Time / tt);
            num += a.Time * tt / (a.Sigma * a.Sigma);
            den += tt * tt / (a.Sigma * a.Sigma);
        }
        ratios.Sort();
        var m = ratios.Count / 2;
        var median = 1 / (ratios.Count % 2 == 1 ? ratios[m] : 0.5 * (ratios[m - 1] + ratios[m]));
        var fit = den / num;
        Assert.True(Math.Abs(median / fit - 1) > 0.05, $"the test needs estimators that differ: median {median:0.0}, fit {fit:0.0}");
        foreach (var terms in new[] { false, true })
        {
            var settings = new QTomographySettings { RayMethod = RayMethod.Straight, UseOpenCl = false, StationTerms = terms };
            var dir = Path.Combine(Path.GetTempPath(), $"quiver_{Guid.NewGuid():N}");
            try
            {
                var r = new QTomography(grid, settings, dir).Run(data, vp, vs);
                var expected = terms ? fit : median;
                Assert.True(Math.Abs(r.ReferenceQ / expected - 1) < 0.01,
                    $"station terms {terms}: reference Q {r.ReferenceQ:0.0}, expected {expected:0.0} (median {median:0.0}, weighted fit {fit:0.0})");
            }
            finally
            {
                if (Directory.Exists(dir)) Directory.Delete(dir, true);
            }
        }
    }

    private static double Median(IEnumerable<double> v)
    {
        var a = v.OrderBy(x => x).ToArray();
        return a.Length == 0 ? double.NaN : a[a.Length / 2];
    }
}
