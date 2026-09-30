// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using TomoStar.Core.Attenuation;
using TomoStar.Core.Forward;
using TomoStar.Core.Geo;
using TomoStar.Core.Model;

namespace TomoStar.Core.Tomography;

/// <summary>A spike: a Gaussian anomaly.</summary>
public sealed record Spike(double Lon, double Lat, double DepthKm, double RadiusKm, double AmplitudePercent);

/// <summary>
/// A tabular body: a plate of the given thickness centred on (lon, lat, depth), long
/// <paramref name="LengthKm"/> along the strike and wide <paramref name="WidthKm"/> down the dip.
/// Dip 0° is a horizontal layer, 90° a vertical dike, anything between a dipping slab (dipping to
/// the right of the strike, Aki &amp; Richards). Its edges are smooth (the amplitude falls to half at
/// the nominal half-size), so the forward solver sees no infinite gradients. It tests whether the
/// data resolve structure elongated in a given direction, which a checkerboard does not show.
/// </summary>
public sealed record TabularBody(double Lon, double Lat, double DepthKm, double StrikeDeg, double DipDeg,
    double LengthKm, double WidthKm, double ThicknessKm, double AmplitudePercent)
{
    /// <summary>Fractional perturbation at a point.</summary>
    public double At(double lon, double lat, double depthKm)
    {
        var east = (lon - Lon) * GeoMath.Deg2Rad * GeoMath.EarthRadiusKm * Math.Cos(Lat * GeoMath.Deg2Rad);
        var north = (lat - Lat) * GeoMath.Deg2Rad * GeoMath.EarthRadiusKm;
        var down = depthKm - DepthKm;
        double st = StrikeDeg * GeoMath.Deg2Rad, dp = DipDeg * GeoMath.Deg2Rad;
        // Unit vectors (east, north, down): along strike, down the dip, normal to the plate.
        var alongStrike = east * Math.Sin(st) + north * Math.Cos(st);
        var downDip = (east * Math.Cos(st) - north * Math.Sin(st)) * Math.Cos(dp) + down * Math.Sin(dp);
        var normal = -(east * Math.Cos(st) - north * Math.Sin(st)) * Math.Sin(dp) + down * Math.Cos(dp);
        return AmplitudePercent / 100 * Edge(alongStrike, LengthKm) * Edge(downDip, WidthKm) * Edge(normal, ThicknessKm);
    }

    // Super-Gaussian of order 6: 1 inside, 0.5 at ±size/2, 0 well outside.
    private static double Edge(double u, double size) => Math.Exp(-Math.Log(2) * Math.Pow(Math.Abs(u) / (0.5 * size), 6));

    public string Describe() => string.Create(System.Globalization.CultureInfo.InvariantCulture,
        $"{(DipDeg <= 0.5 ? "layer" : DipDeg >= 89.5 ? "vertical body" : $"slab dipping {DipDeg:0}°")} {ThicknessKm:0.#} km thick, {AmplitudePercent:+0.#;-0.#}%");
}

/// <summary>A synthetic anomaly pattern for a resolution test.</summary>
public sealed class ResolutionPattern
{
    public bool IsCheckerboard { get; init; } = true;

    /// <summary>Checkerboard cell size, km (horizontal), one full sign change per cell.</summary>
    public double CellHorizontalKm { get; init; } = 20;

    public double CellVerticalKm { get; init; } = 10;
    public double AmplitudePercent { get; init; } = 5;
    public List<Spike> Spikes { get; init; } = [];

    /// <summary>Tabular bodies (layers, dikes, dipping slabs), added to the spikes when not a checkerboard.</summary>
    public List<TabularBody> Bodies { get; init; } = [];

    /// <summary>Fractional perturbation at every node of the grid.</summary>
    public double[] Evaluate(SphericalGrid g)
    {
        var p = new double[g.Count];
        var def = g.Definition;
        for (var k = 0; k < g.Nz; k++)
        for (var j = 0; j < g.Ny; j++)
        for (var i = 0; i < g.Nx; i++)
        {
            var idx = g.Index(i, j, k);
            if (IsCheckerboard)
            {
                // Distances from the south-west top corner along the surface, km.
                var x = GeoMath.SurfaceDistanceKm(def.MinLon, g.LatDeg[j], g.LonDeg[i], g.LatDeg[j]);
                var y = GeoMath.SurfaceDistanceKm(g.LonDeg[i], def.MinLat, g.LonDeg[i], g.LatDeg[j]);
                var z = g.DepthKm[k] - def.MinDepthKm;
                // Smooth sine-product cells (no infinite gradients for the ray tracer); sign
                // alternates every cell in all three directions.
                p[idx] = AmplitudePercent / 100 * Math.Sin(Math.PI * x / CellHorizontalKm) * Math.Sin(Math.PI * y / CellHorizontalKm)
                         * Math.Sin(Math.PI * (z + 0.5 * CellVerticalKm) / CellVerticalKm);
            }
            else
            {
                var c = g.Cartesian(i, j, k);
                foreach (var sp in Spikes)
                {
                    var d = (c - GeoMath.ToCartesian(sp.Lon, sp.Lat, sp.DepthKm)).Length;
                    p[idx] += sp.AmplitudePercent / 100 * Math.Exp(-d * d / (2 * sp.RadiusKm * sp.RadiusKm));
                }
                foreach (var b in Bodies) p[idx] += b.At(g.LonDeg[i], g.LatDeg[j], g.DepthKm[k]);
            }
        }
        return p;
    }

    public string Describe() => IsCheckerboard
        ? $"checkerboard {CellHorizontalKm:0.#}×{CellHorizontalKm:0.#}×{CellVerticalKm:0.#} km, ±{AmplitudePercent:0.#}%"
        : Bodies.Count == 0 ? $"spike test, {Spikes.Count} spike(s)"
        : string.Join("; ", Bodies.Select(b => b.Describe())) + (Spikes.Count > 0 ? $"; {Spikes.Count} spike(s)" : "");
}

/// <summary>The quantity a velocity resolution test perturbs and checks.</summary>
public enum ResolutionTarget
{
    Vp,
    Vs,

    /// <summary>Vp/Vs alone: Vs carries the pattern inverted, Vp is unchanged.</summary>
    VpVs,

    /// <summary>The same pattern in Vp and Vs, Vp/Vs unchanged: any recovered Vp/Vs is an artefact of the P and S coverage.</summary>
    VpVsLeakage
}

public sealed class ResolutionTestResult
{
    public required double[] TruePerturbation { get; init; }
    public required double[] RecoveredPerturbation { get; init; }
    public required double[] Difference { get; init; }
    public required string Quantity { get; init; }
    public required string Description { get; init; }

    /// <summary>Correlation between true and recovered perturbation over sampled nodes.</summary>
    public double Correlation { get; init; }

    /// <summary>
    /// The same over the well-sampled nodes (derivative weight sum above the median of the sampled
    /// ones): nodes touched by a few rays drag the first figure down in sparse geometries.
    /// </summary>
    public double CorrelationWellSampled { get; init; }

    /// <summary>RMS of recovered minus true over the sampled nodes (fraction).</summary>
    public double ErrorRms { get; init; } = double.NaN;

    /// <summary>Derivative weight sum of each node in the synthetic inversion (target phase; for Vp/Vs the smaller of P and S).</summary>
    public double[]? Dws { get; init; }

    public TomographyResult? Velocity { get; init; }
    public QTomographyResult? Attenuation { get; init; }
}

/// <summary>
/// Checkerboard and spike tests (Spakman &amp; Nolet 1988; Humphreys &amp; Clayton 1988): a known pattern
/// is added to the background model, synthetic data are computed for the SAME source-receiver
/// geometry with the same forward solver, Gaussian noise is added, and the data are inverted from
/// the background with the same settings as the real inversion. What comes back shows where, and
/// how faithfully, the data can image structure of that size.
/// </summary>
public static class ResolutionTests
{
    public static ResolutionTestResult Velocity(
        SphericalGrid grid, TomographySettings settings, ObservationSet geometry,
        double[] bgVp, double[] bgVs, ResolutionPattern pattern, Phase target, double noiseSeconds, int seed,
        string workFolder, IProgress<(double, string)>? progress = null, Action<string>? log = null, CancellationToken ct = default,
        VelocityModel1D? background = null) =>
        Velocity(grid, settings, geometry, bgVp, bgVs, pattern, target == Phase.P ? ResolutionTarget.Vp : ResolutionTarget.Vs,
            noiseSeconds, noiseSeconds, seed, workFolder, progress, log, ct, background);

    /// <summary>
    /// The test for any of the quantities of <see cref="ResolutionTarget"/>. For
    /// <see cref="ResolutionTarget.VpVs"/> only Vs changes, so that Vp/Vs carries the pattern and Vp
    /// does not; for <see cref="ResolutionTarget.VpVsLeakage"/> Vp and Vs change alike, the true
    /// Vp/Vs is the background one, and whatever Vp/Vs comes back is made by the different P and S
    /// coverage (the correlation is then undefined and <see cref="ResolutionTestResult.ErrorRms"/> is
    /// the measure). P and S times take their own noise.
    /// </summary>
    public static ResolutionTestResult Velocity(
        SphericalGrid grid, TomographySettings settings, ObservationSet geometry,
        double[] bgVp, double[] bgVs, ResolutionPattern pattern, ResolutionTarget target, double noiseP, double noiseS, int seed,
        string workFolder, IProgress<(double, string)>? progress = null, Action<string>? log = null, CancellationToken ct = default,
        VelocityModel1D? background = null)
    {
        var pert = pattern.Evaluate(grid);
        var ratio = target is ResolutionTarget.VpVs or ResolutionTarget.VpVsLeakage;
        var vpTrue = bgVp.Select((v, i) => target is ResolutionTarget.Vp or ResolutionTarget.VpVsLeakage ? v * (1 + pert[i]) : v).ToArray();
        var vsTrue = bgVs.Select((v, i) => target switch
        {
            ResolutionTarget.Vs or ResolutionTarget.VpVsLeakage => v * (1 + pert[i]),
            ResolutionTarget.VpVs => v / (1 + pert[i]),
            _ => v
        }).ToArray();
        var truth = target == ResolutionTarget.VpVsLeakage ? new double[pert.Length] : pert;

        // Synthetic data in the true model, same geometry and forward method.
        var synth = geometry.CloneGeometry();
        foreach (var e in synth.Events) { e.Fixed = true; e.T0 = 0; }
        progress?.Report((0, "Synthetic data"));
        // Same domain as the inversion below: outside the grid the fixed background, never the pattern.
        var domain = ForwardDomain.For(grid, settings.ForwardRefinement, synth.HasOutsidePoints, synth.Points(), background, bgVp, bgVs);
        var fwd = new RayKernels(grid, domain, settings.RayMethod, Path.Combine(workFolder, "synthetic"), log);
        var rows = fwd.Compute(synth, vpTrue.Select(v => 1 / v).ToArray(), vsTrue.Select(v => 1 / v).ToArray(), null, 0, null, 0, ct);
        var rnd = new Random(seed);
        var keep = new HashSet<int>();
        foreach (var r in rows)
        {
            var a = synth.Arrivals[r.Arrival];
            var noise = a.Phase == Phase.S ? noiseS : noiseP;
            a.Time = r.TCalc + noise * Gaussian(rnd);
            a.Sigma = Math.Max(0.01, noise);
            keep.Add(r.Arrival);
        }
        var data = Subset(synth, keep);

        var s = settings.Clone();
        s.JointHypocentres = false; // hypocentres are exact in a synthetic test
        s.StationCorrections = false;
        s.InvertP = target != ResolutionTarget.Vs || settings.InvertP;
        s.InvertS = target != ResolutionTarget.Vp || (settings.InvertS && data.Arrivals.Any(a => a.Phase == Phase.S));
        var result = new TravelTimeTomography(grid, s, workFolder, log) { Background = background }.Run(data, bgVp, bgVs, progress, ct);
        var rec = target switch
        {
            ResolutionTarget.Vp => result.Vp.Select((v, i) => v / bgVp[i] - 1).ToArray(),
            ResolutionTarget.Vs => result.Vs.Select((v, i) => v / bgVs[i] - 1).ToArray(),
            _ => result.Vp.Select((v, i) => v / result.Vs[i] / (bgVp[i] / bgVs[i]) - 1).ToArray()
        };
        var dws = target == ResolutionTarget.Vp ? result.DwsP : result.DwsS;
        if (ratio) dws = dws.Select((d, i) => Math.Min(d, result.DwsP[i])).ToArray(); // a ratio needs both phases
        var diff = rec.Select((x, i) => x - truth[i]).ToArray();
        var sampled = Enumerable.Range(0, dws.Length).Where(i => dws[i] > 0).ToArray();
        return new ResolutionTestResult
        {
            TruePerturbation = truth, RecoveredPerturbation = rec, Difference = diff,
            Quantity = target switch { ResolutionTarget.Vp => "dVp", ResolutionTarget.Vs => "dVs", _ => "d(Vp/Vs)" },
            Description = pattern.Describe() + (target == ResolutionTarget.VpVsLeakage ? " in Vp and Vs alike (Vp/Vs unchanged)" : ""),
            Correlation = target == ResolutionTarget.VpVsLeakage ? double.NaN : Correlation(pert, rec, dws),
            CorrelationWellSampled = target == ResolutionTarget.VpVsLeakage ? double.NaN : Correlation(pert, rec, dws, wellSampled: true),
            ErrorRms = sampled.Length > 0 ? Math.Sqrt(sampled.Average(i => diff[i] * diff[i])) : double.NaN,
            Dws = dws, Velocity = result
        };
    }

    /// <summary>
    /// Attenuation test: the pattern in 1/Q, synthetic t* along the rays of the velocity model,
    /// inverted with the same settings (parameterisation included). With
    /// <paramref name="velocityLeakage"/> the pattern goes instead into the velocity of the phase
    /// (1/Q uniform): the t* are computed along the rays and slowness of that true model and inverted
    /// in the given one, as when the velocity model of a real Q inversion is wrong. Any 1/Q that comes
    /// back is then leaked from velocity (true 1/Q perturbation zero, measured by
    /// <see cref="ResolutionTestResult.ErrorRms"/>).
    /// </summary>
    public static ResolutionTestResult Attenuation(
        SphericalGrid grid, QTomographySettings settings, ObservationSet geometry, double[] vp, double[] vs,
        ResolutionPattern pattern, double noiseSeconds, int seed, string workFolder,
        IProgress<(double, string)>? progress = null, Action<string>? log = null, CancellationToken ct = default,
        VelocityModel1D? background = null, bool velocityLeakage = false)
    {
        var pert = pattern.Evaluate(grid);
        var q0 = 1 / settings.Q0;
        var qTrue = pert.Select(p => velocityLeakage ? q0 : q0 * (1 + p)).ToArray();
        var vpTrue = velocityLeakage && settings.Phase == Phase.P ? vp.Select((v, i) => v * (1 + pert[i])).ToArray() : vp;
        var vsTrue = velocityLeakage && settings.Phase == Phase.S ? vs.Select((v, i) => v * (1 + pert[i])).ToArray() : vs;
        var truth = velocityLeakage ? new double[pert.Length] : pert;
        var synth = geometry.CloneGeometry();
        var sP = (settings.Phase == Phase.S ? vsTrue : vpTrue).Select(v => 1 / v).ToArray(); // slowness of the phase inverted
        var domain = ForwardDomain.For(grid, settings.ForwardRefinement, synth.HasOutsidePoints, synth.Points(), background, vp, vs);
        var rows = new RayKernels(grid, domain, settings.RayMethod, Path.Combine(workFolder, "synthetic"), log)
            .Compute(synth, sP, vsTrue.Select(v => 1 / v).ToArray(), null, 0, null, 0, ct);
        var rnd = new Random(seed);
        var keep = new HashSet<int>();
        foreach (var r in rows)
        {
            // Outside the grid the attenuation is the reference, as the inversion assumes.
            var t = q0 * r.OutsideTime;
            for (var k = 0; k < r.Nodes.Length; k++) t += r.Length[k] * sP[r.Nodes[k]] * qTrue[r.Nodes[k]];
            var a = synth.Arrivals[r.Arrival];
            a.Time = t + noiseSeconds * Gaussian(rnd);
            a.Sigma = Math.Max(0.002, noiseSeconds);
            keep.Add(r.Arrival);
        }
        var data = Subset(synth, keep);
        var s = settings.Clone();
        s.StationTerms = false;
        s.EstimateQ0 = false;
        s.OutlierMads = 100;
        var result = new QTomography(grid, s, workFolder, log) { Background = background }.Run(data, vp, vs, progress, ct);
        var rec = result.Q.Select(q => (1 / q) / q0 - 1).ToArray();
        var diff = rec.Select((x, i) => x - truth[i]).ToArray();
        var sampled = Enumerable.Range(0, diff.Length).Where(i => result.Dws[i] > 0).ToArray();
        return new ResolutionTestResult
        {
            TruePerturbation = truth, RecoveredPerturbation = rec, Difference = diff,
            Quantity = "d(1/Q)", Description = pattern.Describe() + (velocityLeakage ? $" in V{(settings.Phase == Phase.P ? "p" : "s")}, 1/Q uniform" : ""),
            Correlation = velocityLeakage ? double.NaN : Correlation(pert, rec, result.Dws),
            CorrelationWellSampled = velocityLeakage ? double.NaN : Correlation(pert, rec, result.Dws, wellSampled: true),
            ErrorRms = sampled.Length > 0 ? Math.Sqrt(sampled.Average(i => diff[i] * diff[i])) : double.NaN,
            Dws = result.Dws, Attenuation = result
        };
    }

    private static ObservationSet Subset(ObservationSet s, HashSet<int> keep)
    {
        var o = new ObservationSet();
        o.Events.AddRange(s.Events);
        o.Stations.AddRange(s.Stations);
        o.Arrivals.AddRange(s.Arrivals.Where((_, i) => keep.Contains(i)));
        return o;
    }

    private static double Gaussian(Random r) =>
        Math.Sqrt(-2 * Math.Log(1 - r.NextDouble())) * Math.Cos(2 * Math.PI * r.NextDouble());

    /// <summary>Pearson correlation over nodes that at least one ray samples.</summary>
    private static double Correlation(double[] a, double[] b, double[] dws, bool wellSampled = false)
    {
        var idx = Enumerable.Range(0, a.Length).Where(i => dws[i] > 0).ToArray();
        if (wellSampled && idx.Length > 0)
        {
            var sorted = idx.Select(i => dws[i]).OrderBy(x => x).ToArray();
            var median = sorted[sorted.Length / 2];
            idx = idx.Where(i => dws[i] >= median).ToArray();
        }
        if (idx.Length < 3) return double.NaN;
        var ma = idx.Average(i => a[i]);
        var mb = idx.Average(i => b[i]);
        double sab = 0, saa = 0, sbb = 0;
        foreach (var i in idx)
        {
            sab += (a[i] - ma) * (b[i] - mb);
            saa += (a[i] - ma) * (a[i] - ma);
            sbb += (b[i] - mb) * (b[i] - mb);
        }
        return saa > 0 && sbb > 0 ? sab / Math.Sqrt(saa * sbb) : double.NaN;
    }
}
