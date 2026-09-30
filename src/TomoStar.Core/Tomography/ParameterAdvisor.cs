using System.Globalization;
using TomoStar.Core.Forward;
using TomoStar.Core.Geo;
using TomoStar.Core.Model;

namespace TomoStar.Core.Tomography;

/// <summary>One recommendation with the reason behind it.</summary>
public sealed record Advice(string Parameter, string Value, string Reason);

/// <summary>The grid and settings proposed for a data set, with the reasons.</summary>
public sealed class AdvisorReport
{
    public required GridDefinition Grid { get; init; }
    public required TomographySettings Settings { get; init; }
    public List<Advice> Advice { get; } = [];
    public List<string> Warnings { get; } = [];

    /// <summary>The recommendations as text, one per line.</summary>
    public string Summary() => string.Join(Environment.NewLine,
        Advice.Select(a => $"{a.Parameter} = {a.Value}: {a.Reason}").Concat(Warnings.Select(w => "WARNING: " + w)));
}

/// <summary>What the user fixes in advance; the rest is proposed from the data.</summary>
public sealed class GridAdviceOptions
{
    /// <summary>Horizontal node spacing, km (null: from the station spacing).</summary>
    public double? HorizontalSpacingKm { get; set; }

    /// <summary>Vertical node spacing, km (null: about 12 layers, at most half the horizontal spacing).</summary>
    public double? VerticalSpacingKm { get; set; }

    /// <summary>Margin around the stations and events, km (null: the larger of the spacing and 5 % of the extent).</summary>
    public double? MarginKm { get; set; }

    /// <summary>Top of the grid, km below sea level (null: just above the highest station).</summary>
    public double? TopKm { get; set; }

    /// <summary>Bottom of the grid, km (null: below the deeper events and the turning depth of the long rays).</summary>
    public double? BottomKm { get; set; }
}

/// <summary>
/// Proposes an inversion grid and settings from the data. The user decides; this only says what the
/// data can support and why, following the usual rules of thumb of local earthquake tomography (node
/// spacing no finer than the station spacing near the surface; a data to parameter ratio well above
/// one; a grid deep enough for the rays of the longer paths; Kissling et al. 1994, JGR 99(B10),
/// 19635-19646, for the 1-D reference; Evans et al. 1994, USGS OFR 94-431, for SIMULPS practice).
/// </summary>
public static class ParameterAdvisor
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    public static AdvisorReport Advise(Catalogue catalogue, VelocityModel1D model, GridAdviceOptions? options = null)
    {
        options ??= new GridAdviceOptions();
        var stations = catalogue.StationIndex();
        var events = catalogue.Events.Where(e => e.Picks.Any(p => p.Usable && stations.ContainsKey(p.StationId))).ToList();
        var picks = events.SelectMany(e => e.Picks.Where(p => p.Usable && stations.ContainsKey(p.StationId)).Select(p => (Event: e, Pick: p))).ToList();
        var usedStations = picks.Select(p => stations[p.Pick.StationId]).Distinct().ToList();
        var report = new List<Advice>();
        var warnings = new List<string>();
        if (events.Count == 0 || usedStations.Count == 0)
        {
            // Without picks the grid can only cover what is there: stations and events.
            var lons = catalogue.Stations.Select(s => s.Lon).Concat(catalogue.Events.Select(e => e.Lon)).DefaultIfEmpty(0).ToArray();
            var lats = catalogue.Stations.Select(s => s.Lat).Concat(catalogue.Events.Select(e => e.Lat)).DefaultIfEmpty(0).ToArray();
            var g0 = GridDefinition.FromSpacing(lons.Min() - 0.1, lons.Max() + 0.1, lats.Min() - 0.1, lats.Max() + 0.1, -2, 30,
                options.HorizontalSpacingKm ?? 5, options.VerticalSpacingKm ?? 2);
            var r0 = new AdvisorReport { Grid = g0, Settings = new TomographySettings() };
            r0.Warnings.Add("No usable picks: the grid only covers the stations and events; pick arrivals first.");
            return r0;
        }

        // Lateral extent: envelope of the stations and events actually used, with a margin.
        var allLon = usedStations.Select(s => s.Lon).Concat(events.Select(e => e.Lon)).ToArray();
        var allLat = usedStations.Select(s => s.Lat).Concat(events.Select(e => e.Lat)).ToArray();
        double minLon = allLon.Min(), maxLon = allLon.Max(), minLat = allLat.Min(), maxLat = allLat.Max();
        var box = new GeoBox(minLon, maxLon, minLat, maxLat);
        var extentKm = Math.Max(1, Math.Max(box.WidthKm, box.HeightKm));

        // Horizontal spacing from the median nearest-neighbour station distance.
        var nn = usedStations.Select(a => usedStations.Where(b => b != a)
                .Select(b => GeoMath.SurfaceDistanceKm(a.Lon, a.Lat, b.Lon, b.Lat)).DefaultIfEmpty(extentKm).Min())
            .OrderBy(x => x).ToArray();
        var medianNn = nn[nn.Length / 2];
        var h = options.HorizontalSpacingKm ?? Nice(Math.Clamp(medianNn, Math.Max(0.5, extentKm / 60), Math.Max(2, extentKm / 6)));
        report.Add(new Advice("Horizontal spacing", string.Create(Inv, $"{h:0.##} km"), options.HorizontalSpacingKm != null ? "given"
            : string.Create(Inv, $"median distance between neighbouring stations is {medianNn:0.#} km; structure finer than the station spacing is not resolved near the surface")));

        var margin = options.MarginKm ?? Math.Max(h, 0.05 * extentKm);
        var dLat = margin / 111.195;
        var dLon = margin / (111.195 * Math.Cos(0.5 * (minLat + maxLat) * GeoMath.Deg2Rad));
        minLon -= dLon; maxLon += dLon; minLat -= dLat; maxLat += dLat;
        report.Add(new Advice("Lateral extent", $"{GeoMath.FormatLonRange(minLon, maxLon, "0.###")}, {GeoMath.FormatLatRange(minLat, maxLat, "0.###")}",
            string.Create(Inv, $"envelope of {usedStations.Count} stations and {events.Count} events plus {margin:0.#} km, so no ray leaves the grid")));

        // Depth range: from the highest station to below the deepest (95th percentile) event and below
        // the depth where the rays of the longer source-receiver paths turn.
        var topKm = options.TopKm ?? Math.Floor(-usedStations.Max(s => s.ElevationM) / 1000.0 * 2) / 2 - 0.5;
        var depths = events.Select(e => e.DepthKm).OrderBy(x => x).ToArray();
        var deep = depths[(int)(0.95 * (depths.Length - 1))];
        var offsets = picks.Select(p => GeoMath.SurfaceDistanceKm(p.Event.Lon, p.Event.Lat, stations[p.Pick.StationId].Lon, stations[p.Pick.StationId].Lat))
            .OrderBy(x => x).ToArray();
        var far = offsets[(int)(0.95 * (offsets.Length - 1))];
        var turning = TurningDepthKm(model, far);
        var bottom = options.BottomKm ?? Math.Max(10, Math.Max(deep * 1.2 + h, turning * 1.1));
        var dz = options.VerticalSpacingKm ?? Nice(Math.Clamp(Math.Min(h / 2, (bottom - topKm) / 12), 0.5, 50));
        bottom = topKm + Math.Ceiling((bottom - topKm) / dz - 1e-9) * dz;
        report.Add(new Advice("Depth range", string.Create(Inv, $"{topKm:0.#} to {bottom:0.#} km"),
            string.Create(Inv, $"top above the highest station ({usedStations.Max(s => s.ElevationM):0} m); bottom below 95 % of the events ({deep:0.#} km) and below the turning depth ({turning:0.#} km in '{model.Name}') of the rays to 95 % of the picked stations (up to {far:0} km away)")));
        report.Add(new Advice("Vertical spacing", string.Create(Inv, $"{dz:0.##} km"), options.VerticalSpacingKm != null ? "given"
            : "about 12 layers over the depth range, at most half the horizontal spacing: velocity varies faster with depth than laterally in the crust"));

        var grid = GridDefinition.FromSpacing(minLon, maxLon, minLat, maxLat, topKm, bottom, h, dz);
        var nodes = grid.Count;
        var nP = picks.Count(p => p.Pick.Phase == Phase.P);
        var nS = picks.Count(p => p.Pick.Phase == Phase.S);

        var settings = new TomographySettings();
        // The data to parameter ratio decides how hard to regularise.
        var ratio = (double)(nP + nS) / Math.Max(1, nodes * (nS > 0 ? 2 : 1));
        settings.DampingVelocity = ratio < 1 ? 50 : ratio < 5 ? 20 : 10;
        settings.Smoothing = ratio < 1 ? 50 : ratio < 5 ? 20 : 10;
        report.Add(new Advice("Damping / smoothing", string.Create(Inv, $"{settings.DampingVelocity} / {settings.Smoothing}"),
            string.Create(Inv, $"{nP + nS} arrivals for {nodes} nodes per phase (ratio {ratio:0.0}); refine with 'tomostar lcurve' before the final run")));
        if (ratio < 1) warnings.Add(string.Create(Inv, $"Fewer arrivals than model parameters (ratio {ratio:0.00}): coarsen the grid, use the adaptive grid, or add data."));

        if (nS == 0)
        {
            settings.InvertS = false;
            report.Add(new Advice("Phases", "P only", "no S picks"));
        }
        else if (nS < 0.3 * nP)
        {
            settings.Parameterization = VelocityParameterization.VpVpVs;
            report.Add(new Advice("Phases", "P + S as Vp and Vp/Vs", string.Create(Inv, $"S picks are {100.0 * nS / nP:0}% of P: invert the ratio, which S times constrain more stably than Vs itself (Thurber 1993)")));
        }
        else report.Add(new Advice("Phases", "P + S jointly (Vp, Vs)", $"{nS} S picks for {nP} P picks"));

        // Starting Vp/Vs against the Wadati estimate of the data.
        if (nS > 0)
        {
            var obs = ObservationSet.FromCatalogue(catalogue, new SphericalGrid(grid), true, true);
            var w = Wadati.Estimate(obs);
            if (double.IsFinite(w.VpVs))
            {
                var ratios = Enumerable.Range(0, 41).Select(i => topKm + (deep - topKm) * i / 40.0).Select(z => model.Vp(z) / model.Vs(z)).ToArray();
                double lo = ratios.Min(), hi = ratios.Max();
                var off = Math.Max(Math.Abs(lo / w.VpVs - 1), Math.Abs(hi / w.VpVs - 1));
                var what = string.Create(Inv, $"the data's Vp/Vs is {w.VpVs:0.000} +/- {w.StandardError:0.000} (Wadati diagram of {w.Pairs} S-P pairs, {w.Events} events); '{model.Name}' has {lo:0.00} to {hi:0.00} between {topKm:0.#} and {deep:0.#} km");
                if (off > 0.02)
                {
                    settings.StartVpVs = StartingVpVs.FromData;
                    report.Add(new Advice("Starting Vp/Vs", string.Create(Inv, $"from the data ({w.VpVs:0.000})"),
                        string.Create(Inv, $"{what}, up to {100 * off:0} % away: where S rays are few the start would remain and be shown as structure")));
                }
                else report.Add(new Advice("Starting Vp/Vs", "from the starting model", string.Create(Inv, $"{what}, within {100 * off:0.#} %")));
            }
        }

        settings.RayMethod = RayMethod.FastMarching;
        report.Add(new Advice("Rays", "fast marching", extentKm < 20
            ? "the area is small enough that straight rays would be a fair first approximation, but curved rays cost little here"
            : "rays bend in a velocity gradient; straight rays bias the model at these distances"));
        var hLatKm = (grid.MaxLat - grid.MinLat) / (grid.Ny - 1) * 111.195;
        var aniso = Math.Max(h, hLatKm) / dz;
        settings.ForwardRefinement = aniso > 2.5 ? 3 : 2;
        report.Add(new Advice("Forward refinement", $"{settings.ForwardRefinement}x",
            string.Create(Inv, $"forward cells {h / settings.ForwardRefinement:0.##} x {dz / settings.ForwardRefinement:0.##} km keep the eikonal error well below the pick errors")));
        settings.JointHypocentres = true;
        report.Add(new Advice("Hypocentres", "inverted jointly", "catalogue locations were computed in a different model; holding them fixed maps their errors into velocity"));
        settings.StationCorrections = usedStations.Count >= 8;
        if (settings.StationCorrections)
            report.Add(new Advice("Station corrections", "on", "absorb near-surface structure under each station that the grid cannot represent"));
        report.Add(new Advice("Iterations", "5", "the nonlinear update usually stabilises in 3 to 6 iterations; iterations.csv shows it"));

        var outside = catalogue.Events.Count(e => !(e.Lon >= grid.MinLon && e.Lon <= grid.MaxLon && e.Lat >= grid.MinLat && e.Lat <= grid.MaxLat
                                                    && e.DepthKm >= grid.MinDepthKm && e.DepthKm <= grid.MaxDepthKm));
        if (outside > 0) warnings.Add($"{outside} events fall outside the proposed grid and would be excluded (or used only where their rays cross it, with OutsideData = WhenRaysCross).");
        var fewPhase = events.Count(e => e.Picks.Count(p => p.Usable) < 6);
        if (fewPhase > 0.3 * events.Count) warnings.Add($"{fewPhase} events have fewer than 6 picks; consider requiring more phases per event.");
        if (depths[^1] > 100 || extentKm > 500)
            warnings.Add("Deep or regional data set: check that the starting model extends below the deepest ray.");

        var result = new AdvisorReport { Grid = grid, Settings = settings };
        result.Advice.AddRange(report);
        result.Warnings.AddRange(warnings);
        return result;
    }

    /// <summary>
    /// Depth at which a ray from the surface turns to reach the surface again
    /// <paramref name="distanceKm"/> away in a 1-D model (P velocities): the shallowest turning depth
    /// z with X(p) at least the distance, X(p) = 2 integral of p v / sqrt(1 - p^2 v^2) dz and
    /// p = 1 / max v above z. Rays of that length sample the model down to at least this depth.
    /// </summary>
    public static double TurningDepthKm(VelocityModel1D model, double distanceKm, double maxDepthKm = 200)
    {
        if (distanceKm <= 0) return 0;
        const double dz = 0.1;
        var n = (int)(maxDepthKm / dz);
        var v = new double[n + 1];
        for (var i = 0; i <= n; i++) v[i] = model.Vp(i * dz);
        var vmax = v[0];
        for (var k = 1; k <= n; k++)
        {
            vmax = Math.Max(vmax, v[k]);
            var p = 1 / vmax;
            // Linear velocity within each 0.1 km slab: horizontal distance (c1 - c2) / (p b), with b
            // the gradient and c = sqrt(1 - p^2 v^2), or dz p v / c where it is constant.
            double x = 0;
            for (var i = 0; i < k; i++)
            {
                double v1 = v[i], v2 = Math.Min(v[i + 1], vmax), b = (v2 - v1) / dz;
                double c1 = Math.Sqrt(Math.Max(0, 1 - p * p * v1 * v1)), c2 = Math.Sqrt(Math.Max(0, 1 - p * p * v2 * v2));
                x += Math.Abs(b) > 1e-9 ? (c1 - c2) / (p * b) : c1 > 1e-6 ? dz * p * v1 / c1 : 0;
            }
            if (2 * x >= distanceKm) return k * dz;
        }
        return maxDepthKm;
    }

    /// <summary>Rounds to 1, 2, 2.5 or 5 times a power of ten.</summary>
    public static double Nice(double x)
    {
        if (!(x > 0)) return 1;
        var e = Math.Pow(10, Math.Floor(Math.Log10(x)));
        var m = x / e;
        var n = m < 1.5 ? 1 : m < 2.25 ? 2 : m < 3.75 ? 2.5 : m < 7.5 ? 5 : 10;
        return n * e;
    }
}
