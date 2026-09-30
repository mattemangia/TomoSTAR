using System.Globalization;

namespace TomoStar.Core.Model;

/// <summary>One node of a 1-D model. Two nodes at the same depth make a discontinuity.</summary>
public sealed record VelocityNode(double DepthKm, double Vp, double Vs);

/// <summary>
/// A 1-D velocity model, piecewise linear in depth. At a repeated depth the model is
/// discontinuous: a point exactly on the interface takes the upper value, anything below it the
/// lower one, so Moho, 410 and 660 stay steps instead of being smeared into ramps.
/// Above the first node (topography) the first node's value is used.
/// </summary>
public sealed class VelocityModel1D
{
    public string Name { get; set; } = "custom";
    public string Reference { get; set; } = "";

    /// <summary>DOI of the paper the model comes from (empty for a custom model).</summary>
    public string Doi { get; set; } = "";

    public List<VelocityNode> Nodes { get; set; } = [];

    public double Vp(double depthKm) => Evaluate(depthKm, n => n.Vp);
    public double Vs(double depthKm) => Evaluate(depthKm, n => n.Vs);
    public double Velocity(Phase phase, double depthKm) => phase == Phase.P ? Vp(depthKm) : Vs(depthKm);

    private double Evaluate(double depth, Func<VelocityNode, double> value)
    {
        if (Nodes.Count == 0) throw new InvalidOperationException("The 1-D model has no nodes.");
        if (depth <= Nodes[0].DepthKm) return value(Nodes[0]);
        for (var i = 1; i < Nodes.Count; i++)
        {
            var a = Nodes[i - 1];
            var b = Nodes[i];
            if (depth > b.DepthKm) continue;
            if (b.DepthKm - a.DepthKm <= 1e-12) continue; // discontinuity: move below it
            var t = (depth - a.DepthKm) / (b.DepthKm - a.DepthKm);
            return value(a) + t * (value(b) - value(a));
        }
        return value(Nodes[^1]);
    }

    public void Validate()
    {
        if (Nodes.Count == 0) throw new InvalidOperationException("The 1-D model has no nodes.");
        for (var i = 0; i < Nodes.Count; i++)
        {
            if (!(Nodes[i].Vp > 0) || !(Nodes[i].Vs > 0)) throw new InvalidOperationException($"Node {i}: velocities must be positive.");
            if (Nodes[i].Vs >= Nodes[i].Vp) throw new InvalidOperationException($"Node {i}: Vs must be smaller than Vp.");
            if (i > 0 && Nodes[i].DepthKm < Nodes[i - 1].DepthKm) throw new InvalidOperationException($"Node {i}: depths must not decrease.");
        }
    }

    /// <summary>Layer-cake model from a table of layer tops: each layer is constant down to the next top.</summary>
    public static VelocityModel1D FromLayers(string name, IReadOnlyList<(double TopKm, double Vp, double Vs)> layers, double bottomKm)
    {
        var m = new VelocityModel1D { Name = name };
        for (var i = 0; i < layers.Count; i++)
        {
            var bottom = i + 1 < layers.Count ? layers[i + 1].TopKm : bottomKm;
            m.Nodes.Add(new VelocityNode(layers[i].TopKm, layers[i].Vp, layers[i].Vs));
            m.Nodes.Add(new VelocityNode(bottom, layers[i].Vp, layers[i].Vs));
        }
        return m;
    }

    /// <summary>Parses "depth vp vs" lines (whitespace or comma separated, '#' comments).</summary>
    public static VelocityModel1D Parse(string name, string text)
    {
        var m = new VelocityModel1D { Name = name };
        foreach (var raw in text.Split('\n'))
        {
            var line = raw.Split('#')[0].Trim();
            if (line.Length == 0) continue;
            var parts = line.Split([' ', '\t', ',', ';'], StringSplitOptions.RemoveEmptyEntries);
            if (parts.Length < 3) continue;
            if (!double.TryParse(parts[0], NumberStyles.Float, CultureInfo.InvariantCulture, out var d)) continue;
            var vp = double.Parse(parts[1], CultureInfo.InvariantCulture);
            var vs = double.Parse(parts[2], CultureInfo.InvariantCulture);
            m.Nodes.Add(new VelocityNode(d, vp, vs));
        }
        m.Validate();
        return m;
    }

    public string Format() =>
        string.Join('\n', Nodes.Select(n => string.Create(CultureInfo.InvariantCulture, $"{n.DepthKm:0.###} {n.Vp:0.####} {n.Vs:0.####}")));

    /// <summary>
    /// ak135, continental travel-time edition, surface to the core-mantle boundary at 2891.5 km
    /// (Kennett, Engdahl &amp; Buland 1995, Geophys. J. Int. 122, 108-124; table at
    /// https://rses.anu.edu.au/seismology/ak135/, as distributed with ObsPy/TauP in ak135.tvel).
    /// The liquid outer core (Vs = 0) is left out: no travel-time tomography grid reaches it, and a
    /// zero shear velocity has no slowness.
    /// </summary>
    public static VelocityModel1D Ak135()
    {
        double[] t =
        [
            0.000, 5.8000, 3.4600,
            20.000, 5.8000, 3.4600,
            20.000, 6.5000, 3.8500,
            35.000, 6.5000, 3.8500,
            35.000, 8.0400, 4.4800,
            77.500, 8.0450, 4.4900,
            120.000, 8.0500, 4.5000,
            165.000, 8.1750, 4.5090,
            210.000, 8.3000, 4.5180,
            210.000, 8.3000, 4.5230,
            260.000, 8.4825, 4.6090,
            310.000, 8.6650, 4.6960,
            360.000, 8.8475, 4.7830,
            410.000, 9.0300, 4.8700,
            410.000, 9.3600, 5.0800,
            460.000, 9.5280, 5.1860,
            510.000, 9.6960, 5.2920,
            560.000, 9.8640, 5.3980,
            610.000, 10.0320, 5.5040,
            660.000, 10.2000, 5.6100,
            660.000, 10.7900, 5.9600,
            710.000, 10.9229, 6.0897,
            760.000, 11.0558, 6.2095,
            809.500, 11.1353, 6.2426,
            859.000, 11.2221, 6.2798,
            908.500, 11.3068, 6.3160,
            958.000, 11.3896, 6.3512,
            1007.500, 11.4705, 6.3854,
            1057.000, 11.5495, 6.4187,
            1106.500, 11.6269, 6.4510,
            1156.000, 11.7026, 6.4828,
            1205.500, 11.7766, 6.5138,
            1255.000, 11.8491, 6.5439,
            1304.500, 11.9200, 6.5727,
            1354.000, 11.9895, 6.6008,
            1403.500, 12.0577, 6.6285,
            1453.000, 12.1245, 6.6555,
            1502.500, 12.1912, 6.6815,
            1552.000, 12.2550, 6.7073,
            1601.500, 12.3185, 6.7326,
            1651.000, 12.3819, 6.7573,
            1700.500, 12.4426, 6.7815,
            1750.000, 12.5031, 6.8052,
            1799.500, 12.5631, 6.8286,
            1849.000, 12.6221, 6.8515,
            1898.500, 12.6804, 6.8742,
            1948.000, 12.7382, 6.8972,
            1997.500, 12.7956, 6.9194,
            2047.000, 12.8526, 6.9418,
            2096.500, 12.9096, 6.9627,
            2146.000, 12.9668, 6.9855,
            2195.500, 13.0222, 7.0063,
            2245.000, 13.0783, 7.0281,
            2294.500, 13.1336, 7.0500,
            2344.000, 13.1894, 7.0720,
            2393.500, 13.2465, 7.0931,
            2443.000, 13.3018, 7.1144,
            2492.500, 13.3585, 7.1369,
            2542.000, 13.4156, 7.1586,
            2591.500, 13.4741, 7.1807,
            2640.000, 13.5312, 7.2031,
            2690.000, 13.5900, 7.2258,
            2740.000, 13.6494, 7.2490,
            2789.670, 13.6530, 7.2597,
            2839.330, 13.6566, 7.2704,
            2891.500, 13.6602, 7.2811
        ];
        var m = new VelocityModel1D
        {
            Name = "ak135",
            Reference = "Kennett, B. L. N., Engdahl, E. R., & Buland, R. (1995). Constraints on seismic velocities in the Earth from traveltimes. Geophysical Journal International, 122(1), 108-124."
        };
        for (var i = 0; i < t.Length; i += 3) m.Nodes.Add(new VelocityNode(t[i], t[i + 1], t[i + 2]));
        return m;
    }

    /// <summary>
    /// PREM, isotropic version, evaluated from the published polynomial coefficients in normalised
    /// radius x = r/6371 (Dziewonski &amp; Anderson 1981, Phys. Earth Planet. Inter. 25, 297-356,
    /// Table I) from the surface to the core-mantle boundary at 2891 km, sampled every ~10 km with
    /// both sides of each discontinuity kept. The liquid outer core (Vs = 0) is left out, as for
    /// <see cref="Ak135"/>.
    ///
    /// The 3 km ocean layer of PREM is replaced by the upper crust (the "continental" variant used
    /// for regional work): a water layer with Vs = 0 has no place in a tomography over land.
    /// </summary>
    public static VelocityModel1D Prem()
    {
        const double a = 6371.0;
        // (top depth, bottom depth, vp coefficients, vs coefficients), polynomials in x.
        (double Top, double Bottom, double[] Vp, double[] Vs)[] layers =
        [
            (0, 15, [5.8], [3.2]),
            (15, 24.4, [6.8], [3.9]),
            (24.4, 80, [4.1875, 3.9382], [2.1519, 2.3481]),
            (80, 220, [4.1875, 3.9382], [2.1519, 2.3481]),
            (220, 400, [20.3926, -12.2569], [8.9496, -4.4597]),
            (400, 600, [39.7027, -32.6166], [22.3512, -18.5856]),
            (600, 670, [19.0957, -9.8672], [9.9839, -4.9324]),
            (670, 771, [29.2766, -23.6027, 5.5242, -2.5514], [22.3459, -17.2473, -2.0834, 0.9783]),
            (771, 2741, [24.9520, -40.4673, 51.4832, -26.6419], [11.1671, -13.7818, 17.4575, -9.2777]),
            (2741, 2891, [15.3891, -5.3181, 5.5242, -2.5514], [6.9254, 1.4672, -2.0834, 0.9783])
        ];
        static double Poly(double[] c, double x)
        {
            double v = 0, p = 1;
            foreach (var k in c) { v += k * p; p *= x; }
            return v;
        }
        var m = new VelocityModel1D
        {
            Name = "PREM",
            Reference = "Dziewonski, A. M., & Anderson, D. L. (1981). Preliminary reference Earth model. Physics of the Earth and Planetary Interiors, 25(4), 297-356. Isotropic, ocean layer replaced by upper crust."
        };
        foreach (var (top, bottom, vp, vs) in layers)
        {
            var steps = Math.Max(1, (int)Math.Ceiling((bottom - top) / 10.0));
            for (var s = 0; s <= steps; s++)
            {
                var d = top + (bottom - top) * s / steps;
                var x = (a - d) / a;
                m.Nodes.Add(new VelocityNode(d, Math.Round(Poly(vp, x), 4), Math.Round(Poly(vs, x), 4)));
            }
        }
        // Consecutive layers without a discontinuity (80 km, 771 km, 2741 km) produce duplicate equal nodes;
        // drop exact duplicates so they are not read as zero-thickness steps.
        m.Nodes = m.Nodes.Where((n, i) => i == 0 || n != m.Nodes[i - 1]).ToList();
        return m;
    }

    /// <summary>The published models of <see cref="ReferenceModelLibrary"/>, global and regional.</summary>
    public static IReadOnlyList<VelocityModel1D> BuiltIn() => ReferenceModelLibrary.All.Select(m => m.Model()).ToList();

    /// <summary>
    /// <paramref name="upper"/> down to <paramref name="depthKm"/>, then <paramref name="below"/>. The
    /// join is placed at the first depth from <paramref name="depthKm"/> down where <paramref name="below"/>
    /// is at least as fast as the upper model's bottom, which is held down to it: a regional model whose
    /// half-space is faster than the global crust at that depth is not followed by a slow layer.
    /// </summary>
    public static VelocityModel1D ContinueBelow(VelocityModel1D upper, double depthKm, VelocityModel1D below)
    {
        upper.Validate();
        below.Validate();
        var m = new VelocityModel1D { Name = upper.Name, Reference = upper.Reference, Doi = upper.Doi };
        m.Nodes.AddRange(upper.Nodes.Where(n => n.DepthKm < depthKm));
        // Bottom of the upper model (the upper side of a discontinuity at depthKm).
        var vp = upper.Vp(depthKm);
        var vs = upper.Vs(depthKm);
        m.Nodes.Add(new VelocityNode(depthKm, vp, vs));

        const double eps = 1e-6;
        var candidates = new[] { depthKm }.Concat(below.Nodes.Select(n => n.DepthKm).Where(d => d > depthKm)).Distinct().Order();
        var join = candidates.FirstOrDefault(d => below.Vp(d + eps) >= vp - 1e-9, depthKm);
        if (join > depthKm) m.Nodes.Add(new VelocityNode(join, vp, vs));
        m.Nodes.Add(new VelocityNode(join, below.Vp(join + eps), below.Vs(join + eps)));
        m.Nodes.AddRange(below.Nodes.Where(n => n.DepthKm > join));
        m.Nodes = m.Nodes.Where((n, i) => i == 0 || n != m.Nodes[i - 1]).ToList();
        m.Validate();
        return m;
    }
}
