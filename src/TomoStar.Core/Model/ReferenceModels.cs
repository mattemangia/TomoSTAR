using System.Globalization;

namespace TomoStar.Core.Model;

/// <summary>How far a published model reaches: the whole Earth, a region (hundreds of km) or a single volcano or fault zone.</summary>
public enum ModelScope
{
    Global,
    Regional,
    Local
}

/// <summary>
/// A published 1-D velocity model with its source. Every entry of <see cref="ReferenceModelLibrary"/>
/// comes from a paper whose DOI is registered with Crossref (<see cref="Doi"/>); when the numbers were
/// taken from a data set deposited with the paper, the data set's DOI is <see cref="DataDoi"/>.
///
/// <see cref="Published"/> holds the model exactly as published, down to <see cref="ValidToKm"/> (the
/// last layer of a layered model is the published half-space and is held down to that depth). Below
/// it <see cref="Model"/> continues with another published model (<see cref="ContinueWithId"/>, ak135
/// by default), so that a grid deeper than the local study is not filled with a crustal half-space.
/// </summary>
public sealed class PublishedVelocityModel
{
    public required string Id { get; init; }

    /// <summary>Name shown in the application and stored in the project (e.g. "Carannante et al. 2013: Umbria-Marche").</summary>
    public required string Name { get; init; }

    public required string Region { get; init; }
    public required ModelScope Scope { get; init; }

    /// <summary>Where the model applies. Global models cover the whole Earth.</summary>
    public required GeoBox Bounds { get; init; }

    /// <summary>Among the models that cover an area, the higher priority is chosen first (then the smaller region).</summary>
    public int Priority { get; init; } = 1;

    /// <summary>Full reference (APA).</summary>
    public required string Citation { get; init; }

    /// <summary>DOI of the paper, registered with Crossref.</summary>
    public required string Doi { get; init; }

    /// <summary>Title of the paper as registered with Crossref, used to check the DOI.</summary>
    public required string Title { get; init; }

    /// <summary>DOI of the data set (Zenodo, KIT RADAR, publisher supplement...) the numbers were read from, when there is one.</summary>
    public string DataDoi { get; init; } = "";

    /// <summary>Where the numbers come from and how they were turned into nodes (table, file, figure).</summary>
    public required string Source { get; init; }

    /// <summary>How Vs was obtained when the paper gives only Vp (a Vp/Vs ratio and who published it).</summary>
    public string VsNote { get; init; } = "";

    public string Notes { get; init; } = "";

    /// <summary>Depth down to which the published model is used; below it <see cref="ContinueWithId"/> takes over.</summary>
    public required double ValidToKm { get; init; }

    /// <summary>Model used below <see cref="ValidToKm"/>; empty for the global models.</summary>
    public string ContinueWithId { get; init; } = "";

    /// <summary>Builds the published nodes.</summary>
    public required Func<VelocityModel1D> Build { get; init; }

    /// <summary>The model as published, to <see cref="ValidToKm"/>.</summary>
    public VelocityModel1D Published()
    {
        var m = Build();
        m.Name = Name;
        m.Reference = Citation;
        m.Doi = Doi;
        return m;
    }

    /// <summary>The published model continued below <see cref="ValidToKm"/> by <see cref="ContinueWithId"/>.</summary>
    public VelocityModel1D Model()
    {
        var m = Published();
        if (ContinueWithId.Length == 0) return m;
        var below = ReferenceModelLibrary.Require(ContinueWithId);
        var joined = VelocityModel1D.ContinueBelow(m, ValidToKm, below.Model());
        joined.Name = Name;
        joined.Reference = Citation + $" Below {ValidToKm.ToString("0.#", CultureInfo.InvariantCulture)} km: {below.Name}.";
        joined.Doi = Doi;
        return joined;
    }

    /// <summary>The model's velocities as strictly increasing depth samples (discontinuities moved by 1 m), from the surface down.</summary>
    public (double[] DepthKm, double[] Vp, double[] Vs) Samples(double maxDepthKm = double.PositiveInfinity)
    {
        var nodes = Model().Nodes;
        var d = new List<double>();
        var p = new List<double>();
        var s = new List<double>();
        var first = true;
        foreach (var n in nodes)
        {
            if (n.DepthKm > maxDepthKm) break;
            var z = Math.Max(0, n.DepthKm);
            if (!first && z <= d[^1])
            {
                // Above the surface only the node nearest to it counts; a discontinuity keeps both sides 1 m apart.
                if (n.DepthKm <= 0) { p[^1] = n.Vp; s[^1] = n.Vs; continue; }
                z = d[^1] + 0.001;
            }
            d.Add(z); p.Add(n.Vp); s.Add(n.Vs);
            first = false;
        }
        return (d.ToArray(), p.ToArray(), s.ToArray());
    }

    public override string ToString() => Name;
}

/// <summary>
/// The published 1-D velocity models known to TomoSTAR (the same library as QUIVER), used as starting
/// models of the location, of the minimum 1-D model and of the travel-time and Q tomography.
/// Only models published in peer-reviewed papers are listed; each DOI was checked against Crossref
/// when the library was compiled.
/// </summary>
public static class ReferenceModelLibrary
{
    /// <summary>Default for areas no regional model covers.</summary>
    public const string GlobalDefaultId = "ak135";

    private static readonly Lazy<IReadOnlyList<PublishedVelocityModel>> Entries = new(Create);

    public static IReadOnlyList<PublishedVelocityModel> All => Entries.Value;

    public static bool TryGet(string id, out PublishedVelocityModel model)
    {
        model = All.FirstOrDefault(m => string.Equals(m.Id, id, StringComparison.OrdinalIgnoreCase)
                                        || string.Equals(m.Name, id, StringComparison.OrdinalIgnoreCase))!;
        return model is not null;
    }

    public static PublishedVelocityModel Require(string id) =>
        TryGet(id, out var m) ? m : throw new KeyNotFoundException($"No published 1-D model '{id}'.");

    /// <summary>
    /// The models that apply to an area, best first. A model applies when its region contains the centre
    /// of the area and the area is not much larger than the region (a grid over the whole of southern
    /// Italy centred on Naples is not a Campi Flegrei study). The global models always apply and come last.
    /// </summary>
    public static IReadOnlyList<PublishedVelocityModel> ForArea(GeoBox area)
    {
        var lat = area.CenterLat;
        var lon = area.CenterLon;
        var areaSize = Size(area);
        return All
            .Where(m => m.Scope == ModelScope.Global
                        || (Contains(m.Bounds, lat, lon) && areaSize <= 4.0 * Size(m.Bounds)))
            .OrderBy(m => m.Scope == ModelScope.Global ? 1 : 0)
            .ThenByDescending(m => m.Scope == ModelScope.Global ? (m.Id == GlobalDefaultId ? 1 : 0) : m.Priority)
            .ThenBy(m => m.Scope == ModelScope.Global ? 0 : Size(m.Bounds))
            .ToList();
    }

    /// <summary>The best model for an area, with the reason for the choice.</summary>
    public static (PublishedVelocityModel Model, string Reason) Select(GeoBox area)
    {
        var list = ForArea(area);
        var best = list[0];
        var where = string.Create(CultureInfo.InvariantCulture, $"({area.CenterLat:0.##}°, {area.CenterLon:0.##}°)");
        if (best.Scope == ModelScope.Global)
            return (best, $"No regional published model covers the area centred at {where}; global reference {best.Name}.");
        var others = list.Where(m => m.Scope != ModelScope.Global && m != best).Select(m => m.Name).ToList();
        return (best, $"{best.Name} covers the area centred at {where} ({best.Region})."
                      + (others.Count > 0 ? " Also applicable: " + string.Join("; ", others) + "." : ""));
    }

    private static bool Contains(GeoBox b, double lat, double lon) =>
        lat >= b.MinLat && lat <= b.MaxLat && lon >= b.MinLon && lon <= b.MaxLon;

    private static double Size(GeoBox b) =>
        (b.MaxLat - b.MinLat) * (b.MaxLon - b.MinLon) * Math.Cos(0.5 * (b.MinLat + b.MaxLat) * Math.PI / 180);

    private static readonly GeoBox World = new(-180, 180, -90, 90);

    /// <summary>Nodes from (depth, vp, vs) triplets; a repeated depth is a discontinuity.</summary>
    private static VelocityModel1D Nodes(params double[] t)
    {
        var m = new VelocityModel1D();
        for (var i = 0; i < t.Length; i += 3) m.Nodes.Add(new VelocityNode(t[i], t[i + 1], t[i + 2]));
        m.Validate();
        return m;
    }

    /// <summary>Layered model from (top, vp, vs) triplets, each layer constant down to the next top, the last one down to <paramref name="bottomKm"/>.</summary>
    private static VelocityModel1D Layers(double bottomKm, params double[] t)
    {
        var layers = new List<(double, double, double)>();
        for (var i = 0; i < t.Length; i += 3) layers.Add((t[i], t[i + 1], t[i + 2]));
        var m = VelocityModel1D.FromLayers("", layers, bottomKm);
        m.Validate();
        return m;
    }

    /// <summary>Layered P model from (top, vp) pairs with Vs = Vp / <paramref name="vpVs"/>.</summary>
    private static VelocityModel1D LayersVp(double bottomKm, double vpVs, params double[] t)
    {
        var triplets = new List<double>();
        for (var i = 0; i < t.Length; i += 2) triplets.AddRange([t[i], t[i + 1], Math.Round(t[i + 1] / vpVs, 4)]);
        return Layers(bottomKm, [.. triplets]);
    }

    /// <summary>Layered model from (top, vp, vp/vs) triplets.</summary>
    private static VelocityModel1D LayersVpVs(double bottomKm, params double[] t)
    {
        var triplets = new List<double>();
        for (var i = 0; i < t.Length; i += 3) triplets.AddRange([t[i], t[i + 1], Math.Round(t[i + 1] / t[i + 2], 4)]);
        return Layers(bottomKm, [.. triplets]);
    }

    private static IReadOnlyList<PublishedVelocityModel> Create() =>
    [
        // ---- Global references ----
        new()
        {
            Id = "ak135", Name = "ak135", Region = "Global", Scope = ModelScope.Global, Bounds = World, Priority = 0,
            Citation = "Kennett, B. L. N., Engdahl, E. R., & Buland, R. (1995). Constraints on seismic velocities in the Earth from traveltimes. Geophysical Journal International, 122(1), 108-124.",
            Doi = "10.1111/j.1365-246x.1995.tb03540.x",
            Title = "Constraints on seismic velocities in the Earth from traveltimes",
            Source = "Continental ak135 table (RSES/ANU), surface to the core-mantle boundary (2891.5 km).",
            ValidToKm = 2891.5,
            Build = VelocityModel1D.Ak135
        },
        new()
        {
            Id = "iasp91", Name = "iasp91", Region = "Global", Scope = ModelScope.Global, Bounds = World, Priority = 0,
            Citation = "Kennett, B. L. N., & Engdahl, E. R. (1991). Traveltimes for global earthquake location and phase identification. Geophysical Journal International, 105(2), 429-465.",
            Doi = "10.1111/j.1365-246x.1991.tb06724.x",
            Title = "Traveltimes for global earthquake location and phase identification",
            Source = "IASP91 table as distributed with ObsPy/TauP (iasp91.tvel), surface to the core-mantle boundary (2889 km).",
            ValidToKm = 2889,
            Build = () => Nodes(
                0, 5.8000, 3.3600, 20, 5.8000, 3.3600, 20, 6.5000, 3.7500, 35, 6.5000, 3.7500,
                35, 8.0400, 4.4700, 77.5, 8.0450, 4.4850, 120, 8.0500, 4.5000, 165, 8.1750, 4.5090,
                210, 8.3000, 4.5180, 210, 8.3000, 4.5220, 260, 8.4825, 4.6090, 310, 8.6650, 4.6960,
                360, 8.8475, 4.7830, 410, 9.0300, 4.8700, 410, 9.3600, 5.0700, 460, 9.5280, 5.1760,
                510, 9.6960, 5.2820, 560, 9.8640, 5.3880, 610, 10.0320, 5.4940, 660, 10.2000, 5.6000,
                660, 10.7900, 5.9500, 710, 10.9229, 6.0797, 760, 11.0558, 6.2095, 809.5, 11.1440, 6.2474,
                859, 11.2300, 6.2841, 908.5, 11.3140, 6.3199, 958, 11.3960, 6.3546, 1007.5, 11.4761, 6.3883,
                1057, 11.5543, 6.4211, 1106.5, 11.6308, 6.4530, 1156, 11.7056, 6.4841, 1205.5, 11.7787, 6.5143,
                1255, 11.8504, 6.5438, 1304.5, 11.9205, 6.5725, 1354, 11.9893, 6.6006, 1403.5, 12.0568, 6.6280,
                1453, 12.1231, 6.6547, 1502.5, 12.1881, 6.6809, 1552, 12.2521, 6.7066, 1601.5, 12.3151, 6.7317,
                1651, 12.3772, 6.7564, 1700.5, 12.4383, 6.7807, 1750, 12.4987, 6.8046, 1799.5, 12.5584, 6.8282,
                1849, 12.6174, 6.8514, 1898.5, 12.6759, 6.8745, 1948, 12.7339, 6.8972, 1997.5, 12.7915, 6.9199,
                2047, 12.8487, 6.9423, 2096.5, 12.9057, 6.9647, 2146, 12.9625, 6.9870, 2195.5, 13.0192, 7.0093,
                2245, 13.0758, 7.0316, 2294.5, 13.1325, 7.0540, 2344, 13.1892, 7.0765, 2393.5, 13.2462, 7.0991,
                2443, 13.3034, 7.1218, 2492.5, 13.3610, 7.1449, 2542, 13.4190, 7.1681, 2591.5, 13.4774, 7.1917,
                2641, 13.5364, 7.2156, 2690.5, 13.5961, 7.2398, 2740, 13.6564, 7.2645, 2789.67, 13.6679, 7.2768,
                2839.33, 13.6793, 7.2892, 2889, 13.6908, 7.3015)
        },
        new()
        {
            Id = "PREM", Name = "PREM", Region = "Global", Scope = ModelScope.Global, Bounds = World, Priority = 0,
            Citation = "Dziewonski, A. M., & Anderson, D. L. (1981). Preliminary reference Earth model. Physics of the Earth and Planetary Interiors, 25(4), 297-356.",
            Doi = "10.1016/0031-9201(81)90046-7",
            Title = "Preliminary reference Earth model",
            Source = "Isotropic polynomials of Table I evaluated every ~10 km to the core-mantle boundary (2891 km); the ocean layer is replaced by the upper crust.",
            ValidToKm = 2891,
            Build = VelocityModel1D.Prem
        },

        // ---- Italy ----
        new()
        {
            Id = "Carannante2013", Name = "Carannante et al. 2013: Umbria-Marche / northern-central Apennines",
            Region = "Northern-central Apennines (Umbria-Marche)", Scope = ModelScope.Regional,
            Bounds = new GeoBox(11.5, 15.8, 41.2, 44.7), Priority = 3,
            Citation = "Carannante, S., Monachesi, G., Cattaneo, M., Amato, A., & Chiarabba, C. (2013). Deep structure and tectonics of the northern-central Apennines as seen by regional-scale tomography and 3-D located earthquakes. Journal of Geophysical Research: Solid Earth, 118(10), 5391-5403.",
            Doi = "10.1002/jgrb.50371",
            Title = "Deep structure and tectonics of the northern-central Apennines as seen by regional-scale tomography and 3-D located earthquakes",
            DataDoi = "10.5281/zenodo.16535187",
            Source = "Files ModP9_1D_Marche and ModS9_1D_Marche of the Zenodo data set (Carannante, Cattaneo & Monachesi 2025): SIMULPS nodes, linear between nodes; the -100 km guard node is left out.",
            ValidToKm = 400,
            Build = () => Nodes(
                0, 5.63, 2.80, 4, 6.22, 3.34, 8, 6.23, 3.386, 12, 6.24, 3.39,
                20, 6.26, 3.37, 30, 6.62, 3.64, 80, 7.92, 4.33, 400, 8.50, 4.59),
            ContinueWithId = "ak135"
        },
        new()
        {
            Id = "ChiarabbaFrepoli1997North", Name = "Chiarabba & Frepoli 1997, northern Apennines",
            Region = "Northern Apennines (north of 42.7°N)", Scope = ModelScope.Regional,
            Bounds = new GeoBox(9.5, 14.0, 42.7, 44.6), Priority = 1,
            Citation = "Chiarabba, C., & Frepoli, A. (1997). Minimum 1D velocity models in Central and Southern Italy: a contribution to better constrain hypocentral determinations. Annals of Geophysics, 40(4), 937-954.",
            Doi = "10.4401/ag-3888", Title = "Minimum 1D velocity models in Central and Southern Italy: a contribution to better constrain hypocentral determinations",
            Source = "Minimum 1-D P model of Fig. 3 (VELEST): layer tops 0, 12, 18, 26, 35 km.",
            VsNote = "P model only; Vs = Vp / 1.73 (not from the paper).",
            ValidToKm = 40,
            Build = () => LayersVp(40, 1.73, 0, 5.72, 12, 5.91, 18, 6.03, 26, 6.28, 35, 7.52),
            ContinueWithId = "ak135"
        },
        new()
        {
            Id = "ChiarabbaFrepoli1997South", Name = "Chiarabba & Frepoli 1997, southern Apennines",
            Region = "Southern Apennines (38.8-42.7°N)", Scope = ModelScope.Regional,
            Bounds = new GeoBox(11.5, 18.5, 38.8, 42.7), Priority = 1,
            Citation = "Chiarabba, C., & Frepoli, A. (1997). Minimum 1D velocity models in Central and Southern Italy: a contribution to better constrain hypocentral determinations. Annals of Geophysics, 40(4), 937-954.",
            Doi = "10.4401/ag-3888", Title = "Minimum 1D velocity models in Central and Southern Italy: a contribution to better constrain hypocentral determinations",
            Source = "Minimum 1-D P model of Fig. 6 (VELEST): layer tops 0, 12, 18, 26, 35 km.",
            VsNote = "P model only; Vs = Vp / 1.73 (not from the paper).",
            ValidToKm = 40,
            Build = () => LayersVp(40, 1.73, 0, 5.97, 12, 6.08, 18, 6.34, 26, 6.70, 35, 7.56),
            ContinueWithId = "ak135"
        },
        new()
        {
            Id = "ChiarabbaFrepoli1997Sicily", Name = "Chiarabba & Frepoli 1997: Sicily and Calabria",
            Region = "Sicily and southern Calabria (south of 38.8°N)", Scope = ModelScope.Regional,
            Bounds = new GeoBox(12.0, 16.8, 36.0, 38.8), Priority = 1,
            Citation = "Chiarabba, C., & Frepoli, A. (1997). Minimum 1D velocity models in Central and Southern Italy: a contribution to better constrain hypocentral determinations. Annals of Geophysics, 40(4), 937-954.",
            Doi = "10.4401/ag-3888", Title = "Minimum 1D velocity models in Central and Southern Italy: a contribution to better constrain hypocentral determinations",
            Source = "Minimum 1-D P model of Fig. 9 (VELEST): layer tops 0, 12, 18, 26, 35 km.",
            VsNote = "P model only; Vs = Vp / 1.73 (not from the paper).",
            ValidToKm = 40,
            Build = () => LayersVp(40, 1.73, 0, 5.16, 12, 5.66, 18, 6.25, 26, 6.79, 35, 7.30),
            ContinueWithId = "ak135"
        },
        new()
        {
            Id = "Carannante2015", Name = "Carannante et al. 2015: Po Plain (Emilia)",
            Region = "Po Plain, Emilia and the Ferrara arc", Scope = ModelScope.Regional,
            Bounds = new GeoBox(8.5, 12.6, 44.2, 45.8), Priority = 3,
            Citation = "Carannante, S., Argnani, A., Massa, M., D'Alema, E., Lovati, S., Moretti, M., Cattaneo, M., & Augliera, P. (2015). The May 20 (MW 6.1) and 29 (MW 6.0), 2012, Emilia (Po Plain, northern Italy) earthquakes: New seismotectonic implications from subsurface geology and high-quality hypocenter location. Tectonophysics, 655, 107-123.",
            Doi = "10.1016/j.tecto.2015.05.015",
            Title = "The May 20 (MW 6.1) and 29 (MW 6.0), 2012, Emilia (Po Plain, northern Italy) earthquakes: New seismotectonic implications from subsurface geology and high-quality hypocenter location",
            Source = "Final minimum 1-D Vp and Vs of Fig. 3b (SIMULPS14 nodes every 4 km to 20 km, then 30 and 40 km), digitised from the figure (±0.05 km/s); the Vp/Vs of the digitised nodes matches the published Vp/Vs panel.",
            Notes = "Calibrated on 650 events of the 2012 Emilia sequence; the slow top (Vp 3.5 km/s) is the Plio-Quaternary fill of the Po Plain.",
            ValidToKm = 40,
            Build = () => Nodes(
                0, 3.49, 1.83, 4, 4.06, 2.14, 8, 5.74, 3.16, 12, 5.90, 3.28,
                16, 6.37, 3.74, 20, 6.40, 3.62, 30, 7.68, 4.23, 40, 7.69, 4.23),
            ContinueWithId = "ak135"
        },
        new()
        {
            Id = "Matrullo2013", Name = "Matrullo et al. 2013: Campania-Lucania (Irpinia)",
            Region = "Campania-Lucania Apennines (Irpinia)", Scope = ModelScope.Regional,
            Bounds = new GeoBox(14.3, 16.3, 39.9, 41.6), Priority = 3,
            Citation = "Matrullo, E., De Matteis, R., Satriano, C., Amoroso, O., & Zollo, A. (2013). An improved 1-D seismic velocity model for seismological studies in the Campania-Lucania region (Southern Italy). Geophysical Journal International, 195(1), 460-473.",
            Doi = "10.1093/gji/ggt224",
            Title = "An improved 1-D seismic velocity model for seismological studies in the Campania-Lucania region (Southern Italy)",
            Source = "Minimum 1-D P model of Table 1 (VELEST), layer tops from -1.5 km.",
            VsNote = "Vs = Vp / 1.85, the ratio the paper uses for locations.",
            ValidToKm = 40,
            Build = () => LayersVp(40, 1.85, -1.5, 3.25, 2, 4.72, 6, 5.51, 8, 6.20, 12, 6.55, 14, 6.70, 22, 6.80, 30, 6.85, 35, 7.03),
            ContinueWithId = "ak135"
        },
        new()
        {
            Id = "CampiFlegrei2022", Name = "Giacomuzzi et al. 2025: Campi Flegrei (2022 model)",
            Region = "Campi Flegrei caldera", Scope = ModelScope.Local,
            Bounds = new GeoBox(13.95, 14.30, 40.70, 40.95), Priority = 4,
            Citation = "Giacomuzzi, G., Fonzetti, R., Govoni, A., De Gori, P., & Chiarabba, C. (2025). Causal processes of shallow and deep seismicity at Campi Flegrei caldera. Communications Earth & Environment, 6, 70.",
            Doi = "10.1038/s43247-025-02045-2",
            Title = "Causal processes of shallow and deep seismicity at Campi Flegrei caldera",
            Source = "Supplementary Table S1, 'Vp 2022' and 'Vp/Vs 2022': the 2022 epoch of the 4-D tomography (Giacomuzzi et al. 2024, doi:10.1016/j.epsl.2024.118744) averaged over 8 × 8 km centred on Pozzuoli; layer tops 0-7 km.",
            ValidToKm = 10,
            Build = () => LayersVpVs(10, 0, 1.45, 1.96, 0.5, 2.09, 1.91, 1, 2.73, 1.86, 2, 3.81, 1.62, 3, 3.67, 1.40,
                4, 3.64, 1.34, 5, 4.50, 1.53, 6, 4.80, 1.72, 7, 4.97, 1.82),
            ContinueWithId = "Matrullo2013"
        },
        new()
        {
            Id = "CampiFlegreiOV", Name = "INGV-OV routine model: Campi Flegrei",
            Region = "Campi Flegrei caldera", Scope = ModelScope.Local,
            Bounds = new GeoBox(13.95, 14.30, 40.70, 40.95), Priority = 2,
            Citation = "Giacomuzzi, G., Fonzetti, R., Govoni, A., De Gori, P., & Chiarabba, C. (2025). Causal processes of shallow and deep seismicity at Campi Flegrei caldera. Communications Earth & Environment, 6, 70 (Supplementary Table S1, model of the INGV Osservatorio Vesuviano).",
            Doi = "10.1038/s43247-025-02045-2",
            Title = "Causal processes of shallow and deep seismicity at Campi Flegrei caldera",
            Source = "Supplementary Table S1, 'Vp OV' and 'Vp/Vs OV': the model used routinely by INGV-Osservatorio Vesuviano; layer tops 0, 0.5, 1, 3 km.",
            Notes = "Tuned to the shallow north-eastern seismicity; Giacomuzzi et al. find it places events deeper than 3 km up to 1.5 km too shallow.",
            ValidToKm = 10,
            Build = () => LayersVpVs(10, 0, 1.70, 1.78, 0.5, 2.00, 1.78, 1, 3.00, 1.78, 3, 5.50, 1.78),
            ContinueWithId = "Matrullo2013"
        },

        // ---- Europe ----
        new()
        {
            Id = "Braszus2024", Name = "Braszus et al. 2024, greater Alpine region (GAR1D_PS)",
            Region = "Greater Alpine region", Scope = ModelScope.Regional,
            Bounds = new GeoBox(4.5, 17.0, 43.5, 49.0), Priority = 1,
            Citation = "Braszus, B., Rietbrock, A., Haberland, C., & Ryberg, T. (2024). AI based 1-D P- and S-wave velocity models for the greater alpine region from local earthquake data. Geophysical Journal International, 237(2), 916-930.",
            Doi = "10.1093/gji/ggae077",
            Title = "AI based 1-D P- and S-wave velocity models for the greater alpine region from local earthquake data",
            DataDoi = "10.35097/1965",
            Source = "File GAR1DPS_VEL.mod of the KIT data set (GAR1D_PS_VELEST, lime in Fig. 6): VELEST layer tops -5, 7.5, 16, 32, 46, 52.5, 70 km.",
            ValidToKm = 80,
            Build = () => Layers(80, -5, 5.29, 2.89, 7.5, 5.92, 3.46, 16, 6.20, 3.62, 32, 7.43, 4.24,
                46, 7.65, 4.45, 52.5, 7.87, 4.57, 70, 8.17, 4.69),
            ContinueWithId = "ak135"
        },
        new()
        {
            Id = "vanLaaten2023", Name = "van Laaten et al. 2023: Leipzig-Regensburg fault zone (north)",
            Region = "Northern Leipzig-Regensburg fault zone, Germany", Scope = ModelScope.Regional,
            Bounds = new GeoBox(11.5, 13.3, 50.3, 51.9), Priority = 2,
            Citation = "van Laaten, M., Wegler, U., & Eulenfeld, T. (2023). On the trail of fluids in the northernmost intracontinental earthquake swarm areas of the Leipzig-Regensburg fault zone, Germany. Journal of Seismology, 27, 573-597.",
            Doi = "10.1007/s10950-023-10146-8",
            Title = "On the trail of fluids in the northernmost intracontinental earthquake swarm areas of the Leipzig-Regensburg fault zone, Germany",
            DataDoi = "10.5281/zenodo.6511543",
            Source = "VELEST output model velest_vel.dat of the Zenodo data set: 1 km layers from 0 to 30 km (the -1 km guard layer, with Vs > Vp, is left out).",
            ValidToKm = 35,
            Build = () => Layers(35,
                0, 5.50, 3.18, 1, 5.75, 3.28, 2, 5.82, 3.37, 3, 5.83, 3.37, 4, 5.83, 3.43, 5, 5.83, 3.43,
                6, 5.83, 3.43, 7, 5.84, 3.43, 8, 5.84, 3.46, 9, 5.84, 3.48, 10, 5.86, 3.49, 11, 5.93, 3.52,
                12, 6.04, 3.53, 13, 6.09, 3.55, 14, 6.10, 3.55, 15, 6.12, 3.58, 16, 6.12, 3.58, 17, 6.12, 3.58,
                18, 6.12, 3.65, 19, 6.33, 3.67, 20, 6.36, 3.67, 21, 6.44, 3.71, 22, 6.48, 3.73, 23, 6.48, 3.76,
                24, 6.56, 3.77, 25, 6.63, 3.81, 26, 6.64, 3.81, 27, 7.06, 3.97, 28, 7.20, 4.02, 29, 7.26, 4.07,
                30, 7.36, 4.14),
            ContinueWithId = "ak135"
        },
        new()
        {
            Id = "GlastonburySouthern2025", Name = "Glastonbury-Southern et al. 2025: Fagradalsfjall, Reykjanes Peninsula",
            Region = "Reykjanes Peninsula (Fagradalsfjall), Iceland", Scope = ModelScope.Local,
            Bounds = new GeoBox(-22.6, -21.9, 63.75, 64.1), Priority = 3,
            Citation = "Glastonbury-Southern, E., Winder, T., Rawlinson, N., White, R. S., Greenfield, T., Bacon, C. A., et al. (2025). Pre-existing structures control the orientation of strike-slip faulting during the 2021 dike intrusion at Fagradalsfjall, Iceland. Journal of Geophysical Research: Solid Earth, 130(6), e2024JB030162.",
            Doi = "10.1029/2024jb030162",
            Title = "Pre‐Existing Structures Control the Orientation of Strike‐Slip Faulting During the 2021 Dike Intrusion at Fagradalsfjall, Iceland",
            DataDoi = "10.5281/zenodo.13225186",
            Source = "File velmod.vm of the Zenodo data set (depth, Vp, Vs nodes, linear between nodes).",
            ValidToKm = 25,
            Build = () => Nodes(
                -1, 2.159, 1.220, 0, 2.336, 1.320, 0.75, 2.529, 1.429, 1, 2.546, 1.438,
                1.25, 2.745, 1.551, 1.5, 3.176, 1.794, 1.75, 3.778, 2.134, 2, 4.419, 2.497,
                2.25, 4.842, 2.736, 2.5, 5.047, 2.851, 2.75, 5.283, 2.985, 3, 5.471, 3.091,
                3.25, 5.688, 3.214, 3.5, 5.898, 3.332, 3.75, 6.086, 3.438, 4, 6.280, 3.548,
                4.25, 6.378, 3.603, 4.5, 6.470, 3.655, 4.75, 6.514, 3.680, 5, 6.561, 3.707,
                5.25, 6.566, 3.710, 5.5, 6.624, 3.742, 5.75, 6.681, 3.775, 6, 6.723, 3.798,
                9, 6.764, 3.821, 12, 7.008, 3.959, 15, 7.376, 4.167, 15.1, 7.665, 4.331,
                20, 7.835, 4.427),
            ContinueWithId = "ak135"
        },
        new()
        {
            Id = "Hicks2026", Name = "Hicks et al. 2026: São Jorge, Azores",
            Region = "São Jorge island, Azores", Scope = ModelScope.Local,
            Bounds = new GeoBox(-28.6, -27.6, 38.35, 38.95), Priority = 3,
            Citation = "Hicks, S. P., González, P. J., Lomax, A., Ferreira, A. M. G., Ramalho, R. S., Mitchell, N. C., et al. (2026). Fault-mediated magma propagation and triggered seismicity revealed by the 2022 São Jorge Azores unrest. Nature Communications, 17, 3531.",
            Doi = "10.1038/s41467-026-71668-6",
            Title = "Fault-mediated magma propagation and triggered seismicity revealed by the 2022 São Jorge Azores unrest",
            DataDoi = "10.5281/zenodo.18223639",
            Source = "File SaoJorge_1D_velocitymodel_NLLOC.dat of the Zenodo data set (NonLinLoc layers, tops from -2 km).",
            ValidToKm = 40,
            Build = () => Layers(40, -2, 3.98, 2.44, 1.5, 4.44, 2.48, 7.71, 6.46, 3.72, 9.14, 6.61, 3.74,
                11.80, 7.46, 4.28, 12.77, 7.74, 4.43, 16.77, 7.75, 4.55, 17.39, 8.13, 4.64, 18.90, 8.14, 4.72,
                22.64, 8.14, 4.77, 23.78, 8.14, 4.79),
            ContinueWithId = "ak135"
        },

        // ---- Africa and the Indian Ocean ----
        new()
        {
            Id = "Lavayssiere2022", Name = "Lavayssière et al. 2022: Mayotte (ALav)",
            Region = "Mayotte and the submarine volcano east of it", Scope = ModelScope.Local,
            Bounds = new GeoBox(44.8, 45.9, -13.3, -12.3), Priority = 3,
            Citation = "Lavayssière, A., Crawford, W. C., Saurel, J.-M., Satriano, C., Feuillet, N., Jacques, E., & Komorowski, J.-C. (2022). A new 1D velocity model and absolute locations image the Mayotte seismo-volcanic region. Journal of Volcanology and Geothermal Research, 421, 107440.",
            Doi = "10.1016/j.jvolgeores.2021.107440",
            Title = "A new 1D velocity model and absolute locations image the Mayotte seismo-volcanic region",
            DataDoi = "10.18715/IPGP.2021.kwdfof3p",
            Source = "File ALav_Mayotte1DModel.txt of the IPGP data set (VELEST, NonLinLoc layers).",
            ValidToKm = 70,
            Build = () => Layers(70,
                0, 3.63, 2.61, 4.08, 4.86, 2.71, 6.03, 4.91, 2.71, 8.17, 5.18, 2.90, 10.11, 5.95, 3.93,
                11.08, 6.13, 3.95, 13.03, 6.21, 3.96, 15.17, 6.45, 3.97, 16.14, 6.59, 4.08, 17.11, 6.82, 4.27,
                18.08, 7.00, 4.50, 19.06, 7.15, 4.50, 20.03, 7.32, 4.53, 22.17, 7.79, 4.53, 24.11, 7.88, 4.90,
                26.05, 8.15, 4.95, 28.19, 8.15, 5.00, 30.14, 8.15, 5.06, 45.11, 8.15, 5.07, 50.17, 8.16, 5.07),
            ContinueWithId = "ak135"
        },

        // ---- Americas ----
        new()
        {
            Id = "HadleyKanamori1977", Name = "Hadley & Kanamori 1977, southern California (SCSN)",
            Region = "Southern California", Scope = ModelScope.Regional,
            Bounds = new GeoBox(-121.0, -114.0, 32.0, 36.5), Priority = 2,
            Citation = "Hadley, D., & Kanamori, H. (1977). Seismic structure of the Transverse Ranges, California. Geological Society of America Bulletin, 88(10), 1469-1478.",
            Doi = "10.1130/0016-7606(1977)88<1469:ssottr>2.0.co;2",
            Title = "Seismic structure of the Transverse Ranges, California",
            Source = "Layer tops 0, 5.5, 16, 32 km with Vp 5.5, 6.3, 6.7, 7.8 km/s, as used by the Southern California Seismic Network.",
            VsNote = "Vs = Vp / 1.73, the ratio of the SCSN locations (Hutton, Woessner & Hauksson 2010, doi:10.1785/0120090130).",
            ValidToKm = 40,
            Build = () => LayersVp(40, 1.73, 0, 5.5, 5.5, 6.3, 16, 6.7, 32, 7.8),
            ContinueWithId = "ak135"
        },
        new()
        {
            Id = "Munchmeyer2025", Name = "Münchmeyer et al. 2025: Atacama segment, Chile",
            Region = "Atacama segment of the Chilean subduction zone (24-31°S)", Scope = ModelScope.Regional,
            Bounds = new GeoBox(-72.5, -68.0, -31.0, -24.0), Priority = 2,
            Citation = "Münchmeyer, J., Molina-Ormazabal, D., Marsan, D., Langlais, M., Baez, J. C., Heit, B., et al. (2025). Characterizing the Atacama segment of the Chile subduction margin (24°S-31°S) with >165,000 earthquakes. Journal of Geophysical Research: Solid Earth, 130(7), e2025JB031256.",
            Doi = "10.1029/2025jb031256",
            Title = "Characterizing the Atacama Segment of the Chile Subduction Margin (24°S-31°S) With >165,000 Earthquakes",
            DataDoi = "10.5281/zenodo.15083298",
            Source = "File chile_temp_velest.csv of the Zenodo data set (VELEST layers of the temporary-network model, mean of the ensemble), tops from -10 km, rounded to 1 m/s.",
            ValidToKm = 400,
            Build = () => Layers(400,
                -10, 6.133, 3.577, 0.05, 6.167, 3.590, 3.33, 6.179, 3.597, 6.67, 6.182, 3.598,
                10, 6.186, 3.599, 13.33, 6.263, 3.688, 16.67, 6.358, 3.733, 20, 6.442, 3.792,
                23.33, 6.472, 3.841, 26.67, 6.548, 3.863, 30, 6.612, 3.884, 33.33, 6.749, 3.933,
                36.67, 6.889, 4.043, 40, 7.168, 4.187, 43.33, 7.225, 4.236, 46.67, 7.265, 4.254,
                50, 7.311, 4.269, 53.33, 7.338, 4.287, 56.67, 7.396, 4.321, 60, 7.448, 4.354,
                63.33, 7.600, 4.452, 66.67, 7.818, 4.574, 70, 7.970, 4.701, 73.33, 8.017, 4.712,
                76.67, 8.079, 4.728, 80, 8.131, 4.739, 83.33, 8.296, 4.891, 86.67, 8.458, 4.985,
                90, 8.515, 5.106, 96.67, 8.536, 5.112, 103.33, 8.555, 5.117, 110, 8.580, 5.122,
                116.67, 8.604, 5.132, 123.33, 8.635, 5.163, 130, 8.661, 5.190, 220, 8.666, 5.192,
                310, 8.675, 5.198),
            ContinueWithId = "ak135"
        }
    ];
}
