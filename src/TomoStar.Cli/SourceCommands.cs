// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using System.Globalization;
using System.Text;
using TomoStar.Core.IO;
using TomoStar.Core.Location;
using TomoStar.Core.Model;
using TomoStar.Core.Source;

namespace TomoStar.Cli;

/// <summary>
/// Source studies: first-motion focal mechanisms, stress inversion of mechanisms, waveform
/// cross-correlation times for the double-difference relocation, and Coulomb stress changes.
/// </summary>
public static class SourceCommands
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    // ---- mechanism ---------------------------------------------------------------------------------

    /// <summary>
    /// First-motion focal mechanisms of every event with enough polarities: grid search over strike,
    /// dip and rake with the acceptable set and quality classes of Hardebeck &amp; Shearer (2002,
    /// doi:10.1785/0120010200). The polarity table gives event, station, polarity (+1 up, −1 down),
    /// azimuth and takeoff (degrees from the downward vertical) of the ray at the source.
    /// </summary>
    public static int Mechanism(CommandContext c)
    {
        var polPath = c.Args.Require("polarities", "table of P first-motion polarities");
        var events = c.Catalogue.Events.ToDictionary(e => e.Id, StringComparer.Ordinal);
        var settings = new FirstMotionSettings
        {
            GridStepDeg = c.Args.Double("grid-step") ?? 5,
            ExpectedErrorFraction = c.Args.Double("error-fraction") ?? 0.1,
            MinPolarities = c.Args.Int("min-polarities") ?? 10
        };
        var pols = ReadPolarities(polPath, c.Args.Double("automatic-weight") ?? 0.5);
        c.Log($"Mechanisms: {pols.Count} events with polarities in {polPath}; grid step {settings.GridStepDeg}°, expected wrong polarities {settings.ExpectedErrorFraction:P0}, at least {settings.MinPolarities}.");
        var sb = new StringBuilder("id,time,lon,lat,depth_km,magnitude,npol,strike1,dip1,rake1,strike2,dip2,rake2,p_trend,p_plunge,t_trend,t_plunge,misfit,uncert_deg,prob,stdr,az_gap,to_gap,quality,style,n_acceptable\n");
        var acc = new StringBuilder("id,strike,dip,rake\n");
        var count = new Dictionary<string, int>();
        foreach (var (id, obs) in pols.OrderBy(k => k.Key, StringComparer.Ordinal))
        {
            if (!events.TryGetValue(id, out var e) || obs.Count < settings.MinPolarities) continue;
            var r = FirstMotionInversion.Invert(obs, settings);
            if (r == null) continue;
            var (p, t, _) = MechanismGeometry.Axes(r.Plane1);
            sb.AppendLine(string.Create(Inv,
                $"{id},{e.OriginTime:yyyy-MM-ddTHH:mm:ss.fffZ},{e.Lon:0.######},{e.Lat:0.######},{e.DepthKm:0.###},{e.Magnitude:0.00},{r.Polarities}," +
                $"{r.Plane1.Strike:0},{r.Plane1.Dip:0},{r.Plane1.Rake:0},{r.Plane2.Strike:0},{r.Plane2.Dip:0},{r.Plane2.Rake:0}," +
                $"{p.Trend:0},{p.Plunge:0},{t.Trend:0},{t.Plunge:0},{r.MisfitFraction:0.###},{r.UncertaintyDeg:0.#},{r.Probability:0.##},{r.Stdr:0.##}," +
                $"{r.AzimuthalGapDeg:0},{r.TakeoffGapDeg:0},{r.Quality},{MechanismGeometry.Style(r.Plane1)},{r.Acceptable.Count}"));
            foreach (var a in r.Acceptable.Take(400)) acc.AppendLine(string.Create(Inv, $"{id},{a.Strike:0},{a.Dip:0},{a.Rake:0}"));
            count[r.Quality] = count.GetValueOrDefault(r.Quality) + 1;
        }
        Directory.CreateDirectory(c.RunFolder);
        File.WriteAllText(Path.Combine(c.RunFolder, "mechanisms.csv"), sb.ToString());
        File.WriteAllText(Path.Combine(c.RunFolder, "acceptable_planes.csv"), acc.ToString());
        c.Log($"Mechanisms: {count.Values.Sum()} ({string.Join(", ", count.OrderBy(k => k.Key).Select(k => $"{k.Value} {k.Key}"))}); mechanisms.csv and acceptable_planes.csv in {c.RunFolder}.");
        RunWriter.WriteSummary(c.RunFolder, new { Tool = "First-motion mechanisms", Software = "TomoSTAR", Method = "Hardebeck & Shearer (2002) grid search", Settings = settings, Qualities = count });
        RunWriter.WriteLog(c.RunFolder, c.LogLines);
        return 0;
    }

    private static Dictionary<string, List<PolarityObservation>> ReadPolarities(string path, double automaticWeight)
    {
        var t = CsvTable.Read(path);
        int ce = t.Require("event", "event"), cs = t.Require("station", "station"), cp = t.Require("polarity", "polarity");
        int ca = t.Require("azimuth", "azimuth", "azimuth_deg"), ct = t.Require("takeoff", "takeoff", "takeoff_deg");
        int cm = t.Find("manual");
        var res = new Dictionary<string, List<PolarityObservation>>(StringComparer.Ordinal);
        foreach (var (line, cells) in t.Rows)
        {
            double az = t.Number(cells, ca, line), to = t.Number(cells, ct, line), pol = t.Number(cells, cp, line);
            if (!double.IsFinite(az) || !double.IsFinite(to) || !double.IsFinite(pol) || pol == 0) continue;
            var ev = CsvTable.Cell(cells, ce);
            var st = CsvTable.Cell(cells, cs);
            if (!res.TryGetValue(ev, out var list)) res[ev] = list = [];
            if (list.Any(o => o.StationId.Equals(st, StringComparison.OrdinalIgnoreCase))) continue; // one reading per station
            var w = cm >= 0 && CsvTable.Cell(cells, cm) == "0" ? automaticWeight : 1.0;
            list.Add(new PolarityObservation(st, az, to, Math.Sign(pol), w));
        }
        return res;
    }

    // ---- kagan ------------------------------------------------------------------------------------

    /// <summary>Kagan (1991, doi:10.1111/j.1365-246X.1991.tb06343.x) rotation angle between two double couples.</summary>
    public static int Kagan(CommandContext c)
    {
        if (c.Args.Positional.Count < 2) throw new UsageException("Give two mechanisms as strike/dip/rake.");
        static NodalPlane Parse(string s)
        {
            var f = s.Split('/').Select(x => double.Parse(x, Inv)).ToArray();
            if (f.Length != 3) throw new UsageException($"{s}: expected strike/dip/rake.");
            return new NodalPlane(f[0], f[1], f[2]);
        }
        var (a, b) = (Parse(c.Args.Positional[0]), Parse(c.Args.Positional[1]));
        Console.WriteLine(string.Create(Inv, $"{c.Args.Positional[0]} vs {c.Args.Positional[1]}: Kagan angle {MechanismGeometry.KaganAngle(a, b):0.0}°"));
        return 0;
    }

    // ---- stress -----------------------------------------------------------------------------------

    /// <summary>
    /// Reduced stress tensor from focal mechanisms: linear inversion of Michael (1984,
    /// doi:10.1029/JB089iB13p11517) with the iterative fault-plane selection by fault instability of
    /// Vavryčuk (2014, doi:10.1093/gji/ggu224), friction grid search and bootstrap confidence;
    /// optionally the spatially damped inversion of Hardebeck &amp; Michael (2006,
    /// doi:10.1029/2005JB004144) on a map grid.
    /// </summary>
    public static int Stress(CommandContext c)
    {
        var mechPath = c.Args.Positional.FirstOrDefault() ?? c.Args.Require("mechanisms", "mechanisms.csv of 'tomostar mechanism'");
        var qualities = c.Args.Get("quality") ?? "ABC";
        var centre = c.Args.Numbers("centre");
        var radius = c.Args.Double("radius") ?? double.PositiveInfinity;
        var from = c.Args.Get("from") ?? "0000";
        var to = c.Args.Get("to") ?? "9999";
        var obs = new List<FocalMechanismObservation>();
        var t = CsvTable.Read(mechPath);
        int cid = t.Require("id", "id"), ctime = t.Require("time", "time"), clon = t.Require("lon", "lon"), clat = t.Require("lat", "lat"), cz = t.Require("depth", "depth_km");
        int cq = t.Require("quality", "quality"), cu = t.Find("uncert_deg");
        int[] pl = new[] { "strike1", "dip1", "rake1", "strike2", "dip2", "rake2" }.Select(n => t.Require(n, n)).ToArray();
        foreach (var (line, cells) in t.Rows)
        {
            if (!qualities.Contains(CsvTable.Cell(cells, cq), StringComparison.OrdinalIgnoreCase)) continue;
            var time = CsvTable.Cell(cells, ctime);
            if (string.CompareOrdinal(time, from) < 0 || string.CompareOrdinal(time, to) >= 0) continue;
            double lon = t.Number(cells, clon, line), lat = t.Number(cells, clat, line);
            if (centre is { Length: 2 } && Core.Geo.GeoMath.SurfaceDistanceKm(centre[0], centre[1], lon, lat) > radius) continue;
            var v = pl.Select(k => t.Number(cells, k, line)).ToArray();
            obs.Add(new FocalMechanismObservation
            {
                EventId = CsvTable.Cell(cells, cid), LongitudeDeg = lon, LatitudeDeg = lat, DepthKm = t.Number(cells, cz, line),
                NodalPlane1StrikeDeg = v[0], NodalPlane1DipDeg = v[1], NodalPlane1RakeDeg = v[2],
                NodalPlane2StrikeDeg = v[3], NodalPlane2DipDeg = v[4], NodalPlane2RakeDeg = v[5],
                AngularUncertaintyDeg = cu >= 0 ? t.Number(cells, cu, line, 0) : 0,
                Provider = "first motion", SourceId = CsvTable.Cell(cells, cid), Provenance = Path.GetFileName(mechPath)
            });
        }
        var label = $"quality {qualities}" + (centre is { Length: 2 } ? string.Create(Inv, $", within {radius:0.#} km of {centre[0]:0.###},{centre[1]:0.###}") : "") +
                    (from != "0000" || to != "9999" ? $", {from} to {to}" : "");
        c.Log($"Stress inversion: {obs.Count} mechanisms ({label}).");
        var options = new StressInversionOptions { BootstrapRealizations = c.Args.Int("bootstrap") ?? 2000, MinimumMechanisms = 5 };
        var service = new FocalMechanismStressInversionService();
        var res = service.Invert(obs, options);
        Directory.CreateDirectory(c.RunFolder);
        var sb = new StringBuilder($"# Stress inversion ({res.Method}); {label}\n");
        sb.AppendLine("n,sigma1_az,sigma1_pl,sigma1_cone,sigma2_az,sigma2_pl,sigma2_cone,sigma3_az,sigma3_pl,sigma3_cone,R,R_lo,R_hi,shmax,regime,friction,misfit_deg,instability,quality");
        if (res.Tensor is { } s)
        {
            sb.AppendLine(string.Create(Inv,
                $"{res.MechanismCount},{s.Sigma1.AzimuthDeg:0},{s.Sigma1.PlungeDeg:0},{s.Sigma1.ConfidenceCone95Deg:0},{s.Sigma2.AzimuthDeg:0},{s.Sigma2.PlungeDeg:0},{s.Sigma2.ConfidenceCone95Deg:0}," +
                $"{s.Sigma3.AzimuthDeg:0},{s.Sigma3.PlungeDeg:0},{s.Sigma3.ConfidenceCone95Deg:0},{s.ShapeRatio:0.##},{s.ShapeRatioLower95:0.##},{s.ShapeRatioUpper95:0.##}," +
                $"{StressField.SHmaxAzimuth(s):0},{StressField.Regime(s)},{res.FrictionCoefficient:0.##},{res.MeanSlipMisfitDeg:0.#},{res.MeanFaultInstability:0.##},{res.QualityStatus}"));
            c.Log(string.Create(Inv,
                $"σ1 {s.Sigma1.AzimuthDeg:0}/{s.Sigma1.PlungeDeg:0} (±{s.Sigma1.ConfidenceCone95Deg:0}°), σ2 {s.Sigma2.AzimuthDeg:0}/{s.Sigma2.PlungeDeg:0}, σ3 {s.Sigma3.AzimuthDeg:0}/{s.Sigma3.PlungeDeg:0} (±{s.Sigma3.ConfidenceCone95Deg:0}°), " +
                $"R {s.ShapeRatio:0.00} [{s.ShapeRatioLower95:0.00}–{s.ShapeRatioUpper95:0.00}], S_Hmax {StressField.SHmaxAzimuth(s):0}°, {StressField.RegimeName(StressField.Regime(s))}, friction {res.FrictionCoefficient:0.00}, misfit {res.MeanSlipMisfitDeg:0.0}°."));
        }
        foreach (var m in res.QualityMessages) c.Log("  " + m);
        File.WriteAllText(Path.Combine(c.RunFolder, "stress.csv"), sb.ToString());
        var planes = new StringBuilder("id,plane,instability,slip_misfit_deg\n");
        foreach (var p in res.SelectedFaultPlanes) planes.AppendLine(string.Create(Inv, $"{p.EventId},{p.NodalPlane},{p.Instability:0.###},{p.SlipMisfitDeg:0.#}"));
        File.WriteAllText(Path.Combine(c.RunFolder, "selected_planes.csv"), planes.ToString());

        // Spatially damped inversion on a map grid: --map LON1,LON2,DLON,LAT1,LAT2,DLAT
        if (c.Args.Numbers("map") is { Length: 6 } g)
        {
            var depth = c.Args.Double("map-depth") ?? obs.Select(o => o.DepthKm).DefaultIfEmpty(5).Average();
            var grid = service.InvertSpatial(new StressGridRequest
            {
                Observations = obs,
                MinimumLongitudeDeg = g[0], MaximumLongitudeDeg = g[1], LongitudeStepDeg = g[2],
                MinimumLatitudeDeg = g[3], MaximumLatitudeDeg = g[4], LatitudeStepDeg = g[5],
                MinimumDepthKm = depth, MaximumDepthKm = depth, DepthStepKm = 1,
                HorizontalSearchRadiusKm = c.Args.Double("map-radius") ?? 15, VerticalSearchRadiusKm = c.Args.Double("map-vertical-radius") ?? 15,
                NeighborDamping = c.Args.Double("map-damping") ?? 1,
                InversionOptions = new StressInversionOptions { BootstrapRealizations = c.Args.Int("map-bootstrap") ?? 300, MinimumMechanisms = 5 }
            });
            var gs = new StringBuilder("lon,lat,depth_km,n,n_eff,sigma1_az,sigma1_pl,sigma3_az,sigma3_pl,R,shmax,regime,misfit_deg,quality\n");
            foreach (var nd in grid.Nodes)
            {
                if (nd.Tensor is not { } ns) continue;
                gs.AppendLine(string.Create(Inv,
                    $"{nd.LongitudeDeg:0.###},{nd.LatitudeDeg:0.###},{nd.DepthKm:0.#},{nd.MechanismCount},{nd.EffectiveMechanismCount:0.#},{ns.Sigma1.AzimuthDeg:0},{ns.Sigma1.PlungeDeg:0},{ns.Sigma3.AzimuthDeg:0},{ns.Sigma3.PlungeDeg:0},{ns.ShapeRatio:0.##},{StressField.SHmaxAzimuth(ns):0},{StressField.Regime(ns)},{nd.MeanSlipMisfitDeg:0.#},{nd.QualityStatus}"));
            }
            File.WriteAllText(Path.Combine(c.RunFolder, "stress_map.csv"), gs.ToString());
            c.Log($"{grid.Method}: {grid.Nodes.Count(n => n.Tensor != null)} nodes with a tensor, converged {grid.Converged}, {grid.QualityStatus}; stress_map.csv.");
        }
        RunWriter.WriteSummary(c.RunFolder, new { Tool = "Stress inversion", Software = "TomoSTAR", res.Method, Filters = label, Mechanisms = obs.Count, res.QualityStatus });
        RunWriter.WriteLog(c.RunFolder, c.LogLines);
        return res.Tensor == null ? 2 : 0;
    }

    // ---- crosscorr --------------------------------------------------------------------------------

    /// <summary>
    /// Differential arrival times of neighbouring events from the cross-correlation of their waveforms
    /// (Poupinet et al. 1984, doi:10.1029/JB089iB07p05719; Schaff et al. 2004, doi:10.1785/0120020238),
    /// for 'relocate --method dd --cc dt_cc.csv' (Waldhauser &amp; Ellsworth 2000, doi:10.1785/0120000006).
    /// </summary>
    public static int CrossCorr(CommandContext c)
    {
        var cat = c.Catalogue;
        var dir = c.Args.Require("waveforms", "folder of miniSEED files");
        var archive = WaveformArchive.FromFolder(dir);
        archive.BeforeSeconds = c.Args.Double("before") ?? 10;
        archive.AfterSeconds = c.Args.Double("after") ?? 40;
        var band = c.Args.Numbers("band") ?? [2, 15];
        var s = new CrossCorrelationSettings
        {
            MaxSeparationKm = c.Args.Double("max-separation") ?? 3,
            MaxNeighbours = c.Args.Int("neighbours") ?? 40,
            MaxStationDistanceKm = c.Args.Double("max-distance") ?? 60,
            LowHz = band[0], HighHz = band[1],
            MaxLagSeconds = c.Args.Double("max-lag") ?? 0.3,
            MinCoefficient = c.Args.Double("min-cc") ?? 0.7,
            IncludeS = !c.Args.Flag("p-only"),
            AutomaticPicks = c.Config.AutomaticPicks
        };
        c.Log(s.Describe());
        // Waveforms of every event once, regrouped by station.
        var byStation = new Dictionary<string, List<(string, Trace)>>(StringComparer.OrdinalIgnoreCase);
        foreach (var e in cat.Events)
            foreach (var (sid, traces) in archive.ForEvent(e, c.Log))
            {
                if (!byStation.TryGetValue(sid, out var l)) byStation[sid] = l = [];
                l.AddRange(traces.Select(tr => (e.Id, tr)));
            }
        var stations = cat.StationIndex();
        var list = WaveformCrossCorrelation.Measure(cat.Events, stations,
            sid => byStation.TryGetValue(sid, out var l) ? l : Enumerable.Empty<(string, Trace)>(), s, c.Progress(), c.Log, c.Cancel);
        Directory.CreateDirectory(c.RunFolder);
        WaveformCrossCorrelation.Write(Path.Combine(c.RunFolder, WaveformCrossCorrelation.FileName), list);
        c.Log($"{list.Count} cross-correlation times in {Path.Combine(c.RunFolder, WaveformCrossCorrelation.FileName)}.");
        RunWriter.WriteSummary(c.RunFolder, new { Tool = "Waveform cross-correlation", Software = "TomoSTAR", Settings = s, Times = list.Count });
        RunWriter.WriteLog(c.RunFolder, c.LogLines);
        return 0;
    }

    // ---- coulomb ----------------------------------------------------------------------------------

    /// <summary>
    /// Static Coulomb stress change ΔCFS = Δτ + μ′Δσn (King et al. 1994, doi:10.1785/BSSA0840030935) of a
    /// rectangular source with tapered slip, from the finite-rectangle solution of Okada (1992,
    /// doi:10.1785/BSSA0820021018) in an elastic half-space, on receivers given in a table (--receivers)
    /// and/or on map grids at given depths (--map with --receiver STRIKE/DIP/RAKE).
    /// </summary>
    public static int Coulomb(CommandContext c)
    {
        var src = c.Args.Numbers("source") ?? throw new UsageException("--source LAT,LON,DEPTH,STRIKE,DIP,RAKE,MW,LENGTH_KM,WIDTH_KM is required.");
        if (src.Length != 9) throw new UsageException("--source needs 9 numbers: LAT,LON,DEPTH,STRIKE,DIP,RAKE,MW,LENGTH_KM,WIDTH_KM.");
        var mu = (c.Args.Double("shear-gpa") ?? 30) * 1e9;
        var nu = c.Args.Double("poisson") ?? 0.25;
        var friction = c.Args.Double("friction") ?? 0.4;
        var n = c.Args.Int("patches") ?? 7;
        var mask = (c.Args.Double("mask-km") ?? 0.1) * 1000;
        var patches = CoulombStress.TaperedSource(src[0], src[1], src[2], src[3], src[4], src[5], src[6], src[7], src[8], n, n, mu);
        Directory.CreateDirectory(c.RunFolder);
        c.Log(string.Create(Inv,
            $"Source: Mw {src[6]:0.0} (M0 {CoulombStress.Moment(patches, mu):0.###e+0} N m), {src[3]:0}/{src[4]:0}/{src[5]:0}, {src[7]:0.##} × {src[8]:0.##} km, {n}×{n} patches, cosine-tapered slip up to {patches.Max(p => p.SlipM) * 100:0.0} cm; μ {mu / 1e9:0} GPa, ν {nu:0.##}, μ′ {friction:0.##}."));
        using (var w = new StreamWriter(Path.Combine(c.RunFolder, "source_patches.csv")))
        {
            w.WriteLine("lat,lon,depth_km,strike,dip,rake,length_km,width_km,slip_m");
            foreach (var p in patches)
                w.WriteLine(string.Create(Inv, $"{p.Lat:0.######},{p.Lon:0.######},{p.DepthKm:0.####},{p.StrikeDeg},{p.DipDeg},{p.RakeDeg},{p.LengthKm:0.####},{p.WidthKm:0.####},{p.SlipM:0.######}"));
        }
        if (c.Args.Get("receivers") is { } recPath)
        {
            var t = CsvTable.Read(recPath);
            int ci = t.Require("id", "id"), cla = t.Require("lat", "lat"), clo = t.Require("lon", "lon"), cz = t.Require("depth", "depth", "depth_km");
            int cs = t.Require("strike", "strike"), cd = t.Require("dip", "dip"), cr = t.Require("rake", "rake");
            var receivers = t.Rows.Select(r => new CoulombStress.Receiver(CsvTable.Cell(r.Cells, ci), t.Number(r.Cells, cla, r.Line), t.Number(r.Cells, clo, r.Line),
                t.Number(r.Cells, cz, r.Line), t.Number(r.Cells, cs, r.Line), t.Number(r.Cells, cd, r.Line), t.Number(r.Cells, cr, r.Line))).ToList();
            var res = CoulombStress.Compute(patches, receivers, friction, mu, nu, additionalMaskDistanceM: mask);
            using var w = new StreamWriter(Path.Combine(c.RunFolder, "coulomb_receivers.csv"));
            w.WriteLine("id,lat,lon,depth_km,strike,dip,rake,dcfs_bar,shear_bar,unclamping_bar,mask");
            for (var i = 0; i < res.Count; i++)
            {
                var r = res[i];
                var rc = receivers[i];
                w.WriteLine(string.Create(Inv, $"{r.Id},{r.Lat:0.######},{r.Lon:0.######},{r.DepthKm:0.####},{rc.StrikeDeg},{rc.DipDeg},{rc.RakeDeg},{r.DeltaCfsPa / 1e5:0.#####},{r.ShearPa / 1e5:0.#####},{r.UnclampingPa / 1e5:0.#####},{r.Mask}"));
            }
            c.Log($"ΔCFS on {res.Count} receivers: coulomb_receivers.csv.");
        }
        // --map HALF_WIDTH_KM,STEP_KM --depths Z1,Z2 --receiver STRIKE/DIP/RAKE (grid centred on the source)
        if (c.Args.Numbers("map") is { Length: 2 } m)
        {
            var rec = (c.Args.Get("receiver") ?? $"{src[3]}/{src[4]}/{src[5]}").Split('/').Select(x => double.Parse(x, Inv)).ToArray();
            var depths = c.Args.Numbers("depths") ?? [src[2]];
            int k = (int)Math.Round(2 * m[0] / m[1]) + 1;
            var tag = string.Create(Inv, $"{rec[0]:0}_{rec[1]:0}_{rec[2]:0}");
            using var w = new StreamWriter(Path.Combine(c.RunFolder, $"coulomb_map_{tag}.csv"));
            w.WriteLine("lon,lat,depth_km,dcfs_bar,shear_bar,unclamping_bar");
            foreach (var z in depths)
            {
                var res = CoulombStress.SampleVolume(patches, new CoulombStress.VolumeSamplingRequest("map", src[0], src[1],
                    -m[0], m[0], k, -m[0], m[0], k, z, z, 1, rec[0], rec[1], rec[2]), friction, mu, nu, additionalMaskDistanceM: mask);
                foreach (var r in res)
                    w.WriteLine(string.Create(Inv, $"{r.Lon:0.######},{r.Lat:0.######},{r.DepthKm:0.###},{r.DeltaCfsPa / 1e5:G6},{r.ShearPa / 1e5:G6},{r.UnclampingPa / 1e5:G6}"));
            }
            c.Log(string.Create(Inv, $"ΔCFS maps on receivers {rec[0]:0}/{rec[1]:0}/{rec[2]:0} at {string.Join(", ", depths)} km ({k}×{k} nodes): coulomb_map_{tag}.csv."));
        }
        RunWriter.WriteSummary(c.RunFolder, new
        {
            Tool = "Coulomb stress change", Software = "TomoSTAR", Model = CoulombStress.Provenance, Source = src, ShearModulusPa = mu, Poisson = nu, Friction = friction, Patches = n * n
        });
        RunWriter.WriteLog(c.RunFolder, c.LogLines);
        return 0;
    }
}
