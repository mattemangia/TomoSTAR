// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using System.Globalization;
using System.Text;
using System.Text.Json;
using System.Text.Json.Nodes;
using TomoStar.Core.Attenuation;
using TomoStar.Core.Geo;
using TomoStar.Core.Model;
using TomoStar.Core.Numerics;
using TomoStar.Core.Tomography;

namespace TomoStar.Core.IO;

/// <summary>Which volume formats a run writes besides the native <c>.qvol</c>.</summary>
[Flags]
public enum VolumeFormats
{
    /// <summary>Only the QUIVER volume (<c>.qvol</c>), always written.</summary>
    Qvol = 0,

    /// <summary>Also a CSV with the QUIVER comment header.</summary>
    Csv = 1,

    /// <summary>Also a legacy VTK structured grid.</summary>
    Vtk = 2
}

/// <summary>A volume written by a run: its file and what it holds.</summary>
public sealed record WrittenVolume(string File, string Name, string Quantity, string Units, string Source);

/// <summary>
/// Everything a run leaves behind, in a folder of its own and in the layout QUIVER uses for its own
/// runs, so that the folder can be dropped into a project's <c>runs</c> folder and every diagnostic
/// of QUIVER redrawn from it: <c>summary.json</c>, <c>residuals.csv</c>, <c>hypocentres.csv</c>,
/// <c>stations.csv</c>, <c>rays.qray</c>, <c>G.qcsr</c> (+ <c>G.layout.json</c>), and the volumes under
/// <c>volumes/</c>. Plain tables for scripts and papers are added (<c>iterations.csv</c>,
/// <c>events_relocated.csv</c>, which is again a valid event table for the next step).
/// </summary>
public static class RunWriter
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>Formats written for every volume.</summary>
    public static VolumeFormats Formats { get; set; } = VolumeFormats.Csv;

    // ---- Volumes ---------------------------------------------------------------------------------

    /// <summary>Writes one node field as <c>volumes/&lt;run&gt;_&lt;quantity&gt;.qvol</c> (and the other formats asked for).</summary>
    public static WrittenVolume SaveVolume(string runFolder, string runName, string quantity, string units, SphericalGrid g,
        double[] values, Dictionary<string, string>? meta = null)
    {
        var dir = Path.Combine(runFolder, "volumes");
        Directory.CreateDirectory(dir);
        var file = Path.Combine(dir, $"{runName}_{quantity}{VolumeFile.Extension}");
        var header = new VolumeHeader
        {
            Name = $"{quantity} ({runName})", Quantity = quantity, Units = units, Grid = g.Definition, Source = runName,
            Metadata = meta != null ? new Dictionary<string, string>(meta) : []
        };
        header.Metadata["software"] = "TomoSTAR";
        VolumeFile.Write(file, header, values);
        if (Formats.HasFlag(VolumeFormats.Csv)) VolumeExport.WriteCsv(Path.ChangeExtension(file, ".csv"), header, values);
        if (Formats.HasFlag(VolumeFormats.Vtk)) VolumeExport.WriteVtk(Path.ChangeExtension(file, ".vtk"), header, values);
        return new WrittenVolume(file, header.Name, quantity, units, runName);
    }

    // ---- Velocity tomography -----------------------------------------------------------------------

    /// <summary>
    /// A velocity run: Vp, dVp, DWS and hits for P, and Vs, dVs, Vp/Vs, d(Vp/Vs) and DWS for S when S
    /// was inverted; the cell size of an adaptive run; tables, rays and the last system.
    /// </summary>
    public static List<WrittenVolume> SaveVelocity(string runFolder, string runName, TomographyResult r, TomographySettings settings,
        Catalogue? catalogue, VelocityModel1D? startingModel, IEnumerable<string> log, object? extra = null)
    {
        Directory.CreateDirectory(runFolder);
        var g = r.Grid;
        var volumes = new List<WrittenVolume>();
        var meta = new Dictionary<string, string>
        {
            ["run"] = runName,
            ["forward"] = r.ForwardSolver,
            ["iterations"] = settings.Iterations.ToString(Inv),
            ["damping"] = settings.DampingVelocity.ToString(Inv),
            ["smoothing"] = settings.Smoothing.ToString(Inv),
            ["final RMS (s)"] = r.Iterations[^1].Rms.ToString("0.0000", Inv)
        };
        if (startingModel != null) meta["starting model"] = startingModel.Name;
        if (r.Mesh != null) meta["adaptive cells"] = r.Mesh.CellCount.ToString(Inv);
        var hasS = r.Data.Arrivals.Any(a => a.Phase == Phase.S) && settings.InvertS;
        volumes.Add(SaveVolume(runFolder, runName, "Vp", "km/s", g, r.Vp, meta));
        // Perturbations and ratios cover the whole grid, like Vp and Vs: the smoothing moves the model
        // where no ray passes too. To see only the sampled part, mask with DWS_P or DWS_S.
        volumes.Add(SaveVolume(runFolder, runName, "dVp", "%", g, r.Vp.Select((v, i) => 100 * (v / r.StartVp[i] - 1)).ToArray(), meta));
        volumes.Add(SaveVolume(runFolder, runName, "DWS_P", "km", g, r.DwsP, meta));
        volumes.Add(SaveVolume(runFolder, runName, "Hits_P", "rays", g, r.HitsP.Select(x => (double)x).ToArray(), meta));
        if (hasS)
        {
            volumes.Add(SaveVolume(runFolder, runName, "Vs", "km/s", g, r.Vs, meta));
            volumes.Add(SaveVolume(runFolder, runName, "dVs", "%", g, r.Vs.Select((v, i) => 100 * (v / r.StartVs[i] - 1)).ToArray(), meta));
            volumes.Add(SaveVolume(runFolder, runName, "VpVs", "", g, r.Vp.Select((p, i) => p / r.Vs[i]).ToArray(), meta));
            volumes.Add(SaveVolume(runFolder, runName, "dVpVs", "%", g,
                r.Vp.Select((v, i) => 100 * ((v / r.Vs[i]) / (r.StartVp[i] / r.StartVs[i]) - 1)).ToArray(), meta));
            volumes.Add(SaveVolume(runFolder, runName, "DWS_S", "km", g, r.DwsS, meta));
            volumes.Add(SaveVolume(runFolder, runName, "Hits_S", "rays", g, r.HitsS.Select(x => (double)x).ToArray(), meta));
        }
        if (r.Mesh != null)
        {
            r.Mesh.Save(Path.Combine(runFolder, "mesh.json"));
            volumes.Add(SaveVolume(runFolder, runName, "CellSize", "km", g, r.Mesh.CellSizeKm(), meta));
        }
        // Formal resolution: the diagonal, and the standard deviation in % of the unknown (slowness, or Vp/Vs).
        var sName = settings.Parameterization == VelocityParameterization.VpVpVs ? "VpVs" : "S";
        foreach (var (tag, f) in new[] { ("P", r.FormalP), (sName, r.FormalS) })
        {
            if (f == null) continue;
            volumes.Add(SaveVolume(runFolder, runName, $"Resolution_{tag}", "", g, f.Resolution, meta));
            volumes.Add(SaveVolume(runFolder, runName, $"Sigma_{tag}", "%", g, f.Sigma.Select(x => 100 * x).ToArray(), meta));
            volumes.Add(SaveVolume(runFolder, runName, $"ResolutionLength_{tag}", "km", g, f.LengthKm, meta));
        }
        WriteSummary(runFolder, new
        {
            Tool = "Travel-time tomography",
            Software = "TomoSTAR",
            Run = runName,
            r.ForwardSolver,
            Grid = g.Definition,
            StartingModel = startingModel?.Name,
            Settings = settings,
            r.Iterations,
            Warnings = r.Data.Warnings,
            Volumes = volumes.Select(v => Relative(runFolder, v.File)),
            Extra = extra
        });
        WriteIterations(Path.Combine(runFolder, "iterations.csv"), r.Iterations);
        WriteResiduals(Path.Combine(runFolder, "residuals.csv"), r.Data);
        WriteHypocentres(Path.Combine(runFolder, "hypocentres.csv"), catalogue, r.Data);
        WriteRelocatedEvents(Path.Combine(runFolder, "events_relocated.csv"), catalogue, r.Data);
        WriteStations(Path.Combine(runFolder, "stations.csv"), r.Data);
        WriteRays(Path.Combine(runFolder, "rays.qray"), r.Rays);
        if (r.LastMatrix != null && r.LastLayout != null)
        {
            WriteMatrix(Path.Combine(runFolder, "G.qcsr"), r.LastMatrix);
            File.WriteAllText(Path.Combine(runFolder, "G.layout.json"), JsonSerializer.Serialize(r.LastLayout, TomoJson.Options));
        }
        WriteLog(runFolder, log);
        return volumes;
    }

    // ---- Q tomography ----------------------------------------------------------------------------

    /// <summary>A Q run: Qp (or Qs), 1000/Q, DWS and, for an adaptive run, the cell size.</summary>
    public static List<WrittenVolume> SaveQ(string runFolder, string runName, QTomographyResult r, QTomographySettings s, IEnumerable<string> log, object? extra = null)
    {
        Directory.CreateDirectory(runFolder);
        var meta = new Dictionary<string, string>
        {
            ["run"] = runName, ["reference Q"] = r.ReferenceQ.ToString("0", Inv), ["RMS after (s)"] = r.RmsAfter.ToString("0.00000", Inv),
            ["damping"] = s.Damping.ToString(Inv), ["smoothing"] = s.Smoothing.ToString(Inv)
        };
        var qName = s.Phase == Phase.S ? "Qs" : "Qp";
        var list = new List<WrittenVolume>
        {
            SaveVolume(runFolder, runName, qName, "", r.Grid, r.Q, meta),
            SaveVolume(runFolder, runName, $"1000_over_{qName}", "", r.Grid, r.Q.Select(q => 1000 / q).ToArray(), meta),
            SaveVolume(runFolder, runName, "DWS_Q", "km", r.Grid, r.Dws, meta)
        };
        if (r.Mesh != null)
        {
            r.Mesh.Save(Path.Combine(runFolder, "mesh.json"));
            list.Add(SaveVolume(runFolder, runName, "CellSize", "km", r.Grid, r.Mesh.CellSizeKm(), meta));
        }
        if (r.Formal != null)
        {
            // The standard deviation in % of 1/Q (the unknowns are its fractional changes from the reference).
            list.Add(SaveVolume(runFolder, runName, "Resolution_Q", "", r.Grid, r.Formal.Resolution, meta));
            list.Add(SaveVolume(runFolder, runName, "Sigma_Q", "%", r.Grid, r.Formal.Sigma.Select(x => 100 * x).ToArray(), meta));
            list.Add(SaveVolume(runFolder, runName, "ResolutionLength_Q", "km", r.Grid, r.Formal.LengthKm, meta));
        }
        WriteSummary(runFolder, new
        {
            Tool = "Q tomography", Software = "TomoSTAR", Run = runName, Grid = r.Grid.Definition, Settings = s, r.ReferenceQ, r.RmsBefore, r.RmsAfter,
            Lsqr = new { r.Lsqr.Iterations, r.Lsqr.StopReason, r.Lsqr.ResidualNorm, r.Lsqr.EstimatedConditionNumber },
            LsqrHistory = r.LsqrHistory, StationTerms = r.Data.Stations.Select((st, i) => new { st.Id, Term = i < r.StationTerms.Length ? r.StationTerms[i] : 0 }),
            Warnings = r.Data.Warnings,
            Volumes = list.Select(v => Relative(runFolder, v.File)),
            Extra = extra
        });
        WriteResiduals(Path.Combine(runFolder, "residuals.csv"), r.Data);
        using (var w = new StreamWriter(Path.Combine(runFolder, "station_terms.csv"), false, new UTF8Encoding(false)))
        {
            w.WriteLine("station,lon,lat,tstar_term_s");
            for (var i = 0; i < r.Data.Stations.Count; i++)
            {
                var st = r.Data.Stations[i];
                w.WriteLine(string.Create(Inv, $"{st.Id},{st.Lon:0.#####},{st.Lat:0.#####},{(i < r.StationTerms.Length ? r.StationTerms[i] : 0):0.######}"));
            }
        }
        WriteRays(Path.Combine(runFolder, "rays.qray"), r.Rays);
        if (r.Matrix != null) WriteMatrix(Path.Combine(runFolder, "G.qcsr"), r.Matrix);
        WriteLog(runFolder, log);
        return list;
    }

    // ---- Resolution tests ----------------------------------------------------------------------------

    /// <summary>A checkerboard, spike or body test: true, recovered and difference (in %), and the DWS.</summary>
    public static List<WrittenVolume> SaveResolutionTest(string runFolder, string runName, ResolutionTestResult t, SphericalGrid g, IEnumerable<string> log)
    {
        Directory.CreateDirectory(runFolder);
        var meta = new Dictionary<string, string>
        {
            ["run"] = runName, ["pattern"] = t.Description,
            ["correlation"] = t.Correlation.ToString("0.000", Inv), ["correlation_well_sampled"] = t.CorrelationWellSampled.ToString("0.000", Inv)
        };
        var q = t.Quantity.Replace("(", "").Replace(")", "").Replace("/", "");
        var list = new List<WrittenVolume>
        {
            SaveVolume(runFolder, runName, $"{q}_true", "%", g, t.TruePerturbation.Select(x => 100 * x).ToArray(), meta),
            SaveVolume(runFolder, runName, $"{q}_recovered", "%", g, t.RecoveredPerturbation.Select(x => 100 * x).ToArray(), meta),
            SaveVolume(runFolder, runName, $"{q}_difference", "%", g, t.Difference.Select(x => 100 * x).ToArray(), meta)
        };
        // The node sampling, read by the diagnostics of QUIVER as dws.qvol in the run folder.
        if (t.Dws != null)
        {
            VolumeFile.Write(Path.Combine(runFolder, "dws" + VolumeFile.Extension), new VolumeHeader
            {
                Name = $"DWS ({runName})", Quantity = "DWS", Units = "km", Grid = g.Definition, Source = runName, Metadata = meta
            }, t.Dws);
            list.Add(SaveVolume(runFolder, runName, "DWS", "km", g, t.Dws, meta));
        }
        WriteSummary(runFolder, new
        {
            Tool = "Resolution test", Software = "TomoSTAR", Run = runName, t.Description, t.Quantity,
            Correlation = Finite(t.Correlation), CorrelationWellSampled = Finite(t.CorrelationWellSampled), ErrorRms = Finite(t.ErrorRms),
            Grid = g.Definition, Iterations = t.Velocity?.Iterations, Volumes = list.Select(v => Relative(runFolder, v.File))
        });
        if (t.Velocity != null) WriteRays(Path.Combine(runFolder, "rays.qray"), t.Velocity.Rays);
        if (t.Attenuation != null) WriteRays(Path.Combine(runFolder, "rays.qray"), t.Attenuation.Rays);
        WriteLog(runFolder, log);
        return list;
    }

    private static double? Finite(double x) => double.IsFinite(x) ? x : null;

    // ---- Trade-off (L-curve) -------------------------------------------------------------------------

    /// <summary>
    /// The points of a trade-off run as a table (one row per damping and smoothing pair), the corner of
    /// each L-curve and the recommended pair in <c>summary.json</c>.
    /// </summary>
    public static void SaveTradeOff(string runFolder, string tool, IReadOnlyList<TravelTimeTomography.TradeOffSummary> points,
        IReadOnlyList<(double DataVariance, double ModelVariance)> variances, (double Damping, double Smoothing) chosen,
        IReadOnlyList<(double Smoothing, double Damping)> corners, object settings, IEnumerable<string> log)
    {
        Directory.CreateDirectory(runFolder);
        using (var w = new StreamWriter(Path.Combine(runFolder, "lcurve.csv"), false, new UTF8Encoding(false)))
        {
            w.WriteLine("damping,smoothing,residual_norm,model_norm,roughness,data_variance,model_variance,corner,chosen");
            for (var i = 0; i < points.Count; i++)
            {
                var p = points[i];
                var corner = corners.Any(c => c.Smoothing == p.Smoothing && c.Damping == p.Damping) ? 1 : 0;
                var isChosen = p.Damping == chosen.Damping && p.Smoothing == chosen.Smoothing ? 1 : 0;
                w.WriteLine(string.Create(Inv,
                    $"{p.Damping:G6},{p.Smoothing:G6},{p.ResidualNorm:G8},{p.ModelNorm:G8},{p.Roughness:G8},{variances[i].DataVariance:G8},{variances[i].ModelVariance:G8},{corner},{isChosen}"));
            }
        }
        WriteSummary(runFolder, new
        {
            Tool = tool, Software = "TomoSTAR", Recommended = new { chosen.Damping, chosen.Smoothing },
            Corners = corners.Select(c => new { c.Smoothing, c.Damping }), Settings = settings, Points = points
        });
        WriteLog(runFolder, log);
    }

    // ---- Tables ----------------------------------------------------------------------------------

    /// <summary>Settings and results of a run as indented JSON (the file QUIVER's diagnostics read first).</summary>
    public static void WriteSummary(string runFolder, object summary) =>
        File.WriteAllText(Path.Combine(runFolder, "summary.json"), JsonSerializer.Serialize(summary, TomoJson.Options));

    /// <summary>The log of a run, one message per line.</summary>
    public static void WriteLog(string runFolder, IEnumerable<string> log) =>
        File.WriteAllLines(Path.Combine(runFolder, "log.txt"), log);

    /// <summary>The statistics of each nonlinear iteration.</summary>
    public static void WriteIterations(string path, IReadOnlyList<IterationStats> it)
    {
        using var w = new StreamWriter(path, false, new UTF8Encoding(false));
        w.WriteLine("iteration,used,rejected,rms_s,rms_p_s,rms_s_s,weighted_rms_s,variance_reduction_percent,lsqr_iterations,lsqr_stop,update_norm,step_scale,roughness,mean_hypocentre_shift_km,rows,columns,nonzeros");
        foreach (var s in it)
            w.WriteLine(string.Create(Inv,
                $"{s.Iteration},{s.Used},{s.Rejected},{s.Rms:0.#####},{s.RmsP:0.#####},{s.RmsS:0.#####},{s.WeightedRms:0.#####},{s.VarianceReductionPercent:0.##},{s.LsqrIterations},{Csv(s.LsqrStop)},{s.UpdateNorm:G6},{s.StepScale:0.###},{s.Roughness:G6},{s.MeanHypocentreShiftKm:0.###},{s.MatrixRows},{s.MatrixColumns},{s.MatrixNonZeros}"));
    }

    /// <summary>Every arrival with its initial and final residual and, when rejected, the reason.</summary>
    public static void WriteResiduals(string path, ObservationSet d)
    {
        using var w = new StreamWriter(path, false, new UTF8Encoding(false));
        w.WriteLine("event,station,phase,distance_km,azimuth_deg,event_depth_km,sigma_s,observed_s,initial_residual_s,final_residual_s,rejected,reason");
        foreach (var a in d.Arrivals)
        {
            var e = d.Events[a.Event];
            var s = d.Stations[a.Station];
            var dist = GeoMath.SurfaceDistanceKm(e.Lon, e.Lat, s.Lon, s.Lat);
            var az = GeoMath.Azimuth(e.Lon, e.Lat, s.Lon, s.Lat);
            w.WriteLine(string.Create(Inv,
                $"{e.Id},{s.Id},{a.Phase},{dist:0.###},{az:0.#},{e.DepthKm:0.###},{a.Sigma:0.####},{a.Time:0.####},{a.InitialResidual:0.####},{a.Residual:0.####},{(a.Rejected ? 1 : 0)},{Csv(a.RejectReason)}"));
        }
    }

    /// <summary>A free-text CSV field: commas would split it, so they become semicolons.</summary>
    private static string Csv(string s) => s.Replace(',', ';').ReplaceLineEndings(" ");

    /// <summary>Starting and final hypocentres with the shifts (the table QUIVER's diagnostics read).</summary>
    public static void WriteHypocentres(string path, Catalogue? catalogue, ObservationSet d)
    {
        var events = catalogue?.EventIndex();
        using var w = new StreamWriter(path, false, new UTF8Encoding(false));
        w.WriteLine("event,catalog_lon,catalog_lat,catalog_depth_km,lon,lat,depth_km,origin_shift_s,horizontal_shift_km,depth_shift_km,fixed");
        foreach (var e in d.Events)
        {
            double cLon = e.Lon, cLat = e.Lat, cDep = e.DepthKm;
            if (events != null && events.TryGetValue(e.Id, out var ev)) (cLon, cLat, cDep) = (ev.Lon, ev.Lat, ev.DepthKm);
            w.WriteLine(string.Create(Inv,
                $"{e.Id},{cLon:0.#####},{cLat:0.#####},{cDep:0.###},{e.Lon:0.#####},{e.Lat:0.#####},{e.DepthKm:0.###},{e.T0:0.####},{GeoMath.SurfaceDistanceKm(cLon, cLat, e.Lon, e.Lat):0.###},{e.DepthKm - cDep:0.###},{(e.Fixed ? 1 : 0)}"));
        }
    }

    /// <summary>
    /// The final hypocentres as an event table (id, time, lon, lat, depth_km, magnitude, rms_s,
    /// phases): the input of the next step (t* measurement, Q tomography, a second run).
    /// </summary>
    public static void WriteRelocatedEvents(string path, Catalogue? catalogue, ObservationSet d)
    {
        var events = catalogue?.EventIndex();
        var count = new int[d.Events.Count];
        var sum = new double[d.Events.Count];
        foreach (var a in d.Arrivals)
        {
            if (a.Rejected || double.IsNaN(a.Residual)) continue;
            count[a.Event]++;
            sum[a.Event] += a.Residual * a.Residual;
        }
        using var w = new StreamWriter(path, false, new UTF8Encoding(false));
        w.WriteLine("id,time,lon,lat,depth_km,magnitude,rms_s,phases,fixed");
        for (var i = 0; i < d.Events.Count; i++)
        {
            var e = d.Events[i];
            var mag = events != null && events.TryGetValue(e.Id, out var ev) ? ev.Magnitude : double.NaN;
            var t = e.Reference.AddSeconds(e.T0);
            w.WriteLine(string.Create(Inv,
                $"{e.Id},{t:yyyy-MM-ddTHH:mm:ss.fffZ},{e.Lon:0.######},{e.Lat:0.######},{e.DepthKm:0.####},{(double.IsFinite(mag) ? mag.ToString("0.0#", Inv) : "")},{(count[i] > 0 ? Math.Sqrt(sum[i] / count[i]) : double.NaN):0.####},{count[i]},{(e.Fixed ? 1 : 0)}"));
        }
    }

    /// <summary>Stations with their corrections and mean residuals (also a valid station table).</summary>
    public static void WriteStations(string path, ObservationSet d)
    {
        var p = new List<double>[d.Stations.Count];
        var q = new List<double>[d.Stations.Count];
        for (var i = 0; i < p.Length; i++) { p[i] = []; q[i] = []; }
        foreach (var a in d.Arrivals)
            if (!a.Rejected && !double.IsNaN(a.Residual)) (a.Phase == Phase.P ? p : q)[a.Station].Add(a.Residual);
        using var w = new StreamWriter(path, false, new UTF8Encoding(false));
        w.WriteLine("station,lon,lat,depth_km,correction_p_s,correction_s_s,n_p,n_s,mean_residual_p_s,mean_residual_s_s");
        for (var i = 0; i < d.Stations.Count; i++)
        {
            var s = d.Stations[i];
            w.WriteLine(string.Create(Inv,
                $"{s.Id},{s.Lon:0.#####},{s.Lat:0.#####},{s.DepthKm:0.###},{s.CorrectionP:0.####},{s.CorrectionS:0.####},{p[i].Count},{q[i].Count},{(p[i].Count > 0 ? p[i].Average() : double.NaN):0.####},{(q[i].Count > 0 ? q[i].Average() : double.NaN):0.####}"));
        }
    }

    // ---- Rays and matrix (binary formats of QUIVER) ------------------------------------------------------

    /// <summary>
    /// Rays as "QRAY", count, then per ray: arrival, event, station (int32), phase (byte), point
    /// count (int32), and the longitudes, latitudes and depths (float32 each).
    /// </summary>
    public static void WriteRays(string path, IReadOnlyList<StoredRay> rays)
    {
        using var fs = new BufferedStream(File.Create(path), 1 << 20);
        using var w = new BinaryWriter(fs);
        w.Write("QRAY"u8);
        w.Write(rays.Count);
        foreach (var r in rays)
        {
            w.Write(r.Arrival); w.Write(r.Event); w.Write(r.Station); w.Write((byte)r.Phase); w.Write(r.Lon.Length);
            foreach (var x in r.Lon) w.Write(x);
            foreach (var x in r.Lat) w.Write(x);
            foreach (var x in r.Depth) w.Write(x);
        }
    }

    /// <summary>Reads a ray file written by <see cref="WriteRays"/>.</summary>
    public static List<StoredRay> ReadRays(string path)
    {
        var list = new List<StoredRay>();
        using var r = new BinaryReader(new BufferedStream(File.OpenRead(path), 1 << 20));
        if (Encoding.ASCII.GetString(r.ReadBytes(4)) != "QRAY") throw new InvalidDataException("Not a ray file.");
        var n = r.ReadInt32();
        for (var i = 0; i < n; i++)
        {
            int arr = r.ReadInt32(), ev = r.ReadInt32(), st = r.ReadInt32();
            var ph = (Phase)r.ReadByte();
            var m = r.ReadInt32();
            var lon = new float[m]; var lat = new float[m]; var dep = new float[m];
            for (var k = 0; k < m; k++) lon[k] = r.ReadSingle();
            for (var k = 0; k < m; k++) lat[k] = r.ReadSingle();
            for (var k = 0; k < m; k++) dep[k] = r.ReadSingle();
            list.Add(new StoredRay(arr, ev, st, ph, lon, lat, dep));
        }
        return list;
    }

    /// <summary>A CSR matrix as "QCSR", rows, columns, nonzeros (int32), row starts, column indices (int32) and values (float32).</summary>
    public static void WriteMatrix(string path, CsrMatrix m)
    {
        using var w = new BinaryWriter(new BufferedStream(File.Create(path), 1 << 20));
        w.Write("QCSR"u8);
        w.Write(m.RowCount); w.Write(m.ColumnCount); w.Write(m.NonZeroCount);
        foreach (var x in m.RowStart) w.Write(x);
        foreach (var x in m.ColumnIndex) w.Write(x);
        foreach (var x in m.Values) w.Write((float)x);
    }

    private static string Relative(string folder, string file) => Path.GetRelativePath(folder, file).Replace('\\', '/');

    // ---- QUIVER project registration -------------------------------------------------------------------

    /// <summary>
    /// Adds a finished run to a QUIVER project: the run folder is copied into the project's
    /// <c>runs</c> folder (where the Diagnostics tab lists it), its volumes into <c>volumes</c>, and each
    /// volume is appended to the "Volumes" list of <c>project.json</c>, so it appears in the project
    /// tree the next time the project is opened. Every other field of project.json is kept as it is.
    /// QUIVER must not have the project open meanwhile: it would write its own list back on saving.
    /// </summary>
    public static int RegisterInQuiverProject(string projectPath, string runFolder, string runName, IReadOnlyList<WrittenVolume> volumes, string group = "")
    {
        var project = CatalogueReader.ProjectFolder(projectPath);
        var file = Path.Combine(project, "project.json");
        if (!File.Exists(file)) throw new FileNotFoundException($"{project} is not a QUIVER project (no project.json).");
        var runsDir = Path.Combine(project, "runs", runName);
        CopyFolder(runFolder, runsDir, skip: "volumes");
        var volDir = Path.Combine(project, "volumes");
        Directory.CreateDirectory(volDir);
        var root = JsonNode.Parse(File.ReadAllText(file))!.AsObject();
        var key = root.Select(kv => kv.Key).FirstOrDefault(k => k.Equals("Volumes", StringComparison.OrdinalIgnoreCase)) ?? "Volumes";
        if (root[key] is not JsonArray list) root[key] = list = [];
        var added = 0;
        foreach (var v in volumes)
        {
            var target = Path.Combine(volDir, Path.GetFileName(v.File));
            File.Copy(v.File, target, overwrite: true);
            var rel = Path.GetRelativePath(project, target).Replace('\\', '/');
            // Replace an entry of the same file (a run registered twice), otherwise append.
            for (var i = list.Count - 1; i >= 0; i--)
                if (list[i] is JsonObject o && o.TryGetPropertyValue("File", out var f) && f?.GetValue<string>() == rel) list.RemoveAt(i);
            list.Add(new JsonObject
            {
                ["Id"] = Guid.NewGuid().ToString("N"), ["Name"] = v.Name, ["File"] = rel, ["Quantity"] = v.Quantity, ["Units"] = v.Units,
                ["CreatedUtc"] = DateTime.UtcNow, ["Source"] = $"TomoSTAR {v.Source}", ["Group"] = group
            });
            added++;
        }
        var tmp = file + ".tomostar.tmp";
        File.WriteAllText(tmp, root.ToJsonString(new JsonSerializerOptions { WriteIndented = true }));
        File.Copy(file, file + ".bak", overwrite: true);
        File.Move(tmp, file, overwrite: true);
        return added;
    }

    private static void CopyFolder(string from, string to, string? skip = null)
    {
        Directory.CreateDirectory(to);
        foreach (var f in Directory.GetFiles(from)) File.Copy(f, Path.Combine(to, Path.GetFileName(f)), overwrite: true);
        foreach (var d in Directory.GetDirectories(from))
            if (!string.Equals(Path.GetFileName(d), skip, StringComparison.OrdinalIgnoreCase))
                CopyFolder(d, Path.Combine(to, Path.GetFileName(d)));
    }
}
