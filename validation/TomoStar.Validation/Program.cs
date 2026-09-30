// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

// The TomoSTAR side of the validation cases. Every command computes one thing with the TomoSTAR
// library and writes it to plain files (JSON for small results, raw little-endian float64 or
// float32 arrays for large ones) that the Python scripts of this folder compare with the reference
// codes. Nothing here is tuned to the references: the library is called as the program calls it.
//
//     tomostar-validation eikonal-taup OUT.json
//     tomostar-validation eikonal-gradient OUT.json
//     tomostar-validation eikonal-3d OUT_FOLDER
//     tomostar-validation lsqr G.qcsr b.f64 DAMP ITERATIONS OUT.f64
//     tomostar-validation filter IN.f32 RATE LOW HIGH ORDER ZEROPHASE(0|1) OUT.f32
//     tomostar-validation dpss N NW K OUT.f64
//     tomostar-validation stalta IN.f32 RATE STA LTA OUT.f32
//     tomostar-validation aic IN.f32 OUT.f64
//     tomostar-validation mseed FILE OUT_FOLDER
//     tomostar-validation quakeml FILE_OR_FOLDER OUT.json
//     tomostar-validation stationxml FILE FREQUENCIES(comma) OUT.json
//     tomostar-validation remove-response MSEED STATIONXML F1 F2 F3 F4 WATERLEVEL OUT_FOLDER

using System.Globalization;
using System.Text;
using System.Text.Json;
using TomoStar.Core.Forward;
using TomoStar.Core.Geo;
using TomoStar.Core.IO;
using TomoStar.Core.Model;
using TomoStar.Core.Numerics;
using TomoStar.Core.Signal;
using TomoStar.Core.Tomography;

CultureInfo.DefaultThreadCurrentCulture = CultureInfo.InvariantCulture;
CultureInfo.CurrentCulture = CultureInfo.InvariantCulture;
var json = new JsonSerializerOptions { WriteIndented = true, NumberHandling = System.Text.Json.Serialization.JsonNumberHandling.AllowNamedFloatingPointLiterals };
double D(string s) => double.Parse(s, CultureInfo.InvariantCulture);

switch (args.Length > 0 ? args[0] : "")
{
    case "eikonal-taup": EikonalTauP(args[1]); break;
    case "eikonal-gradient": EikonalGradient(args[1]); break;
    case "eikonal-3d": Eikonal3D(args[1]); break;
    case "eikonal-homogeneous": EikonalHomogeneous(args[1]); break;
    case "lsqr": LsqrCase(args[1], args[2], D(args[3]), int.Parse(args[4]), args[5]); break;
    case "filter":
    {
        var x = ReadF32(args[1]);
        Butterworth.Filter(x, D(args[2]), D(args[3]), D(args[4]), int.Parse(args[5]), args[6] == "1");
        WriteF32(args[7], x);
        break;
    }
    case "dpss":
    {
        var (tapers, _) = Multitaper.Tapers(int.Parse(args[1]), D(args[2]), int.Parse(args[3]));
        WriteF64(args[4], tapers.SelectMany(t => t).ToArray());
        break;
    }
    case "stalta": WriteF32(args[5], AutoPicker.StaLta(ReadF32(args[1]), D(args[2]), D(args[3]), D(args[4]))); break;
    case "aic":
    {
        var x = ReadF32(args[1]);
        WriteF64(args[2], AutoPicker.Aic(x, 0, x.Length - 1));
        break;
    }
    case "mseed": MSeed(args[1], args[2]); break;
    case "quakeml": QuakeMlCase(args[1], args[2]); break;
    case "stationxml": StationXmlCase(args[1], args[2], args[3]); break;
    case "remove-response": RemoveResponse(args[1], args[2], D(args[3]), D(args[4]), D(args[5]), D(args[6]), D(args[7]), args[8]); break;
    default:
        Console.Error.WriteLine("Usage: see the header of Program.cs.");
        return 2;
}
return 0;

// ---- Eikonal ----------------------------------------------------------------------------------------

// First-arrival P and S times in ak135 from sources at 5, 15 and 30 km to surface receivers 10 to
// 500 km away, on a spherical grid of 0.02 deg x 0.02 deg x 1 km (the reference: TauP).
void EikonalTauP(string output)
{
    var def = new GridDefinition { MinLon = 12.0, MaxLon = 18.4, MinLat = 41.4, MaxLat = 42.3, MinDepthKm = 0, MaxDepthKm = 100, Nx = 321, Ny = 46, Nz = 101 };
    var model = VelocityModel1D.Ak135();
    var fm = ForwardModel.From1D(def, model);
    var results = new List<object>();
    var clock = System.Diagnostics.Stopwatch.StartNew();
    foreach (var depth in new[] { 5.0, 15.0, 30.0 })
    foreach (var phase in new[] { Phase.P, Phase.S })
    {
        var s = fm.Slowness(phase);
        var t = SphericalEikonal.Solve(fm.Metric, s, SphericalEikonal.PointSource(fm.Grid, s, 12.2, 42.0, depth));
        for (var d = 10.0; d <= 500; d += 10)
        {
            var (lon, lat) = GeoMath.Destination(12.2, 42.0, 90, d);
            results.Add(new { Phase = phase.ToString(), SourceDepthKm = depth, DistanceKm = d, DistanceDeg = d / GeoMath.EarthRadiusKm * GeoMath.Rad2Deg, Time = fm.Grid.Interpolate(t, lon, lat, 0) });
        }
    }
    File.WriteAllText(output, JsonSerializer.Serialize(new { Grid = def, Model = model.Name, Seconds = clock.Elapsed.TotalSeconds, Times = results }, json));
}

// A linear increase of velocity with depth, v = 4 + 0.05 z km/s, for which the travel time is known
// in closed form; the times at every node of four depth levels are written with the analytic ones.
void EikonalGradient(string output)
{
    var def = new GridDefinition { MinLon = 13.0, MaxLon = 13.8, MinLat = 42.5, MaxLat = 43.1, MinDepthKm = 0, MaxDepthKm = 40, Nx = 161, Ny = 121, Nz = 81 };
    var model = new VelocityModel1D { Name = "gradient", Nodes = [new(0, 4.0, 4.0 / 1.75), new(100, 9.0, 9.0 / 1.75)] };
    var fm = ForwardModel.From1D(def, model);
    var g = fm.Grid;
    var (lon0, lat0, z0) = (13.4, 42.8, 15.0);
    var t = SphericalEikonal.Solve(fm.Metric, fm.SlownessP, SphericalEikonal.PointSource(g, fm.SlownessP, lon0, lat0, z0));
    var src = GeoMath.ToCartesian(lon0, lat0, z0);
    var rows = new List<double[]>();
    foreach (var z in new[] { 0.0, 5.0, 20.0, 30.0 })
    {
        var k = (int)Math.Round((z - def.MinDepthKm) / def.DDepth);
        for (var j = 0; j < g.Ny; j += 4)
        for (var i = 0; i < g.Nx; i += 4)
        {
            // Analytic: T = arccosh(1 + b^2 r^2 / (2 v1 v2)) / b, b the gradient, r the distance.
            var r = (g.Cartesian(i, j, k) - src).Length;
            if (r < 3) continue; // the seeded neighbourhood of the source is imposed, not computed
            double v1 = model.Vp(z0), v2 = model.Vp(g.DepthKm[k]), b = 0.05;
            var exact = Math.Acosh(1 + b * b * r * r / (2 * v1 * v2)) / b;
            rows.Add([g.LonDeg[i], g.LatDeg[j], g.DepthKm[k], r, t[g.Index(i, j, k)], exact]);
        }
    }
    File.WriteAllText(output, JsonSerializer.Serialize(new { Grid = def, Columns = new[] { "lon", "lat", "depth_km", "distance_km", "fmm_s", "exact_s" }, Rows = rows }, json));
}

// The grid of the 3-D case with a uniform 6 km/s: the times of a surface source at every node, whose
// exact values are the distances divided by 6 km/s (the near-source accuracy of both solvers).
void EikonalHomogeneous(string output)
{
    var def = new GridDefinition { MinLon = 12.8, MaxLon = 13.7, MinLat = 42.3, MaxLat = 43.2, MinDepthKm = -2, MaxDepthKm = 30, Nx = 91, Ny = 91, Nz = 65 };
    var g = new SphericalGrid(def);
    var s = Enumerable.Repeat(1 / 6.0, g.Count).ToArray();
    var t = SphericalEikonal.Solve(new GridMetric(g), s, SphericalEikonal.PointSource(g, s, 13.05, 42.55, 0));
    WriteF64(output, t);
}

// A 3-D model (the regional 1-D model of Carannante et al. 2013 with an 8 % checkerboard of 16 km
// cells): the velocity, the travel-time fields of three surface sources, and the rays from ten
// hypocentres to each source, for the comparison with PyKonal.
void Eikonal3D(string folder)
{
    Directory.CreateDirectory(folder);
    var def = new GridDefinition { MinLon = 12.8, MaxLon = 13.7, MinLat = 42.3, MaxLat = 43.2, MinDepthKm = -2, MaxDepthKm = 30, Nx = 91, Ny = 91, Nz = 65 };
    var g = new SphericalGrid(def);
    var model = ReferenceModelLibrary.Require("Carannante2013").Model();
    var (vp1, _) = TravelTimeTomography.StartingModel(g, model);
    var pattern = new ResolutionPattern { IsCheckerboard = true, CellHorizontalKm = 16, CellVerticalKm = 10, AmplitudePercent = 8 };
    var p = pattern.Evaluate(g);
    var vp = vp1.Select((v, i) => v * (1 + p[i])).ToArray();
    WriteF64(Path.Combine(folder, "vp.f64"), vp);
    var s = vp.Select(v => 1 / v).ToArray();
    var metric = new GridMetric(g);
    var sources = new[] { (Lon: 13.05, Lat: 42.55, Depth: 0.0), (Lon: 13.25, Lat: 42.95, Depth: 0.0), (Lon: 13.5, Lat: 42.7, Depth: 0.0) };
    var rnd = new Random(5);
    var hypos = Enumerable.Range(0, 10).Select(_ => (Lon: 12.95 + 0.6 * rnd.NextDouble(), Lat: 42.45 + 0.6 * rnd.NextDouble(), Depth: 4 + 12 * rnd.NextDouble())).ToArray();
    var rays = new List<object>();
    for (var q = 0; q < sources.Length; q++)
    {
        var src = sources[q];
        var t = SphericalEikonal.Solve(metric, s, SphericalEikonal.PointSource(g, s, src.Lon, src.Lat, src.Depth));
        WriteF64(Path.Combine(folder, $"t{q}.f64"), t);
        var tablePath = Path.Combine(folder, $"t{q}.qttb");
        TravelTimeTable.Write(tablePath, new TravelTimeTableHeader { StationId = $"S{q}", Phase = Phase.P, Lon = src.Lon, Lat = src.Lat, DepthKm = src.Depth, Grid = def }, t);
        using var table = TravelTimeTable.Open(tablePath);
        foreach (var h in hypos)
        {
            var ray = RayTracing.Backtrack(table, h.Lon, h.Lat, h.Depth, 0.25 * g.MinSpacingKm());
            if (ray == null) continue;
            // The ray's time as the tomography computes it: the path integral of the trilinear basis.
            var kernel = RayTracing.Kernel(g, ray);
            rays.Add(new
            {
                Source = q, h.Lon, h.Lat, h.Depth, FieldTime = table.Time(h.Lon, h.Lat, h.Depth), KernelTime = RayTracing.Time(kernel, s),
                LengthKm = ray.LengthKm(), RayLon = ray.Lon, RayLat = ray.Lat, RayDepth = ray.Depth
            });
        }
    }
    File.WriteAllText(Path.Combine(folder, "case.json"), JsonSerializer.Serialize(new { Grid = def, Sources = sources.Select(x => new { x.Lon, x.Lat, x.Depth }), Rays = rays }, json));
}

// ---- LSQR ----------------------------------------------------------------------------------------

void LsqrCase(string matrix, string rhs, double damp, int iterations, string output)
{
    var a = ReadQcsr(matrix);
    var b = ReadF64(rhs);
    var clock = System.Diagnostics.Stopwatch.StartNew();
    var r = Lsqr.Solve(a, b, damp, iterations, 1e-12, 1e-12, 1e12);
    WriteF64(output, r.Solution);
    File.WriteAllText(output + ".json", JsonSerializer.Serialize(new { r.Iterations, r.StopReason, r.ResidualNorm, Seconds = clock.Elapsed.TotalSeconds, SimdWidth = SimdVector.Width }, json));
}

CsrMatrix ReadQcsr(string path)
{
    using var r = new BinaryReader(File.OpenRead(path));
    if (Encoding.ASCII.GetString(r.ReadBytes(4)) != "QCSR") throw new InvalidDataException("not a QCSR file");
    int rows = r.ReadInt32(), cols = r.ReadInt32(), nnz = r.ReadInt32();
    var start = new int[rows + 1];
    for (var i = 0; i <= rows; i++) start[i] = r.ReadInt32();
    var ci = new int[nnz];
    for (var i = 0; i < nnz; i++) ci[i] = r.ReadInt32();
    var v = new double[nnz];
    for (var i = 0; i < nnz; i++) v[i] = r.ReadSingle();
    return new CsrMatrix(rows, cols, start, ci, v);
}

// ---- Data formats and responses ------------------------------------------------------------------

void MSeed(string file, string folder)
{
    Directory.CreateDirectory(folder);
    var list = new List<object>();
    var n = 0;
    foreach (var t in MiniSeed.ReadFile(file))
    {
        var name = $"trace{n++}.f32";
        WriteF32(Path.Combine(folder, name), t.Data);
        list.Add(new { t.Network, t.Station, t.Location, t.Channel, Start = t.StartTime.ToString("o"), t.SampleRate, Samples = t.Data.Length, File = name });
    }
    File.WriteAllText(Path.Combine(folder, "traces.json"), JsonSerializer.Serialize(list, json));
}

void QuakeMlCase(string input, string output)
{
    var files = Directory.Exists(input) ? Directory.GetFiles(input, "*.xml").Order().ToArray() : [input];
    var events = files.SelectMany(QuakeMl.ReadFile).Select(e => new
    {
        e.Id, Time = e.OriginTime.ToString("o"), e.Lon, e.Lat, e.DepthKm, e.Magnitude,
        Picks = e.Picks.Select(p => new { p.StationId, Phase = p.Phase.ToString(), Time = p.Time.ToString("o"), p.Sigma, Origin = p.Origin.ToString(), p.Quality })
    });
    File.WriteAllText(output, JsonSerializer.Serialize(events, json));
}

void StationXmlCase(string file, string frequencies, string output)
{
    var r = StationXml.ReadFile(file);
    var freqs = frequencies.Split(',').Select(D).ToArray();
    var responses = new List<object>();
    foreach (var (nslc, epochs) in r.Responses)
    foreach (var e in epochs)
    {
        if (e.Response == null) continue;
        var values = freqs.Select(f => e.Response.Evaluate(f)).ToArray();
        responses.Add(new
        {
            Nslc = nslc, Start = e.Start?.ToString("o"), End = e.End?.ToString("o"), e.Response.InputUnits,
            Amplitude = values.Select(v => v.Magnitude), Phase = values.Select(v => v.Phase)
        });
    }
    File.WriteAllText(output, JsonSerializer.Serialize(new
    {
        Frequencies = freqs, Stations = r.Stations.Select(s => new { s.Id, s.Lon, s.Lat, s.ElevationM }), Responses = responses
    }, json));
}

void RemoveResponse(string mseed, string stationXml, double f1, double f2, double f3, double f4, double waterLevel, string folder)
{
    Directory.CreateDirectory(folder);
    var inventory = StationXml.ReadFile(stationXml);
    var list = new List<object>();
    var n = 0;
    foreach (var t in MiniSeed.ReadFile(mseed))
    {
        var response = StationXml.ResponseAt(inventory.Responses, t.Nslc, t.StartTime);
        if (response == null || t.Data.Length < 1000) continue;
        var y = ResponseRemoval.Remove(t.Data, t.SampleRate, response, GroundMotion.Velocity, f1, f2, f3, f4, waterLevel);
        var name = $"vel{n++}.f32";
        WriteF32(Path.Combine(folder, name), y);
        list.Add(new { Nslc = t.Nslc, Start = t.StartTime.ToString("o"), t.SampleRate, File = name });
    }
    File.WriteAllText(Path.Combine(folder, "traces.json"), JsonSerializer.Serialize(list, json));
}

// ---- Raw arrays --------------------------------------------------------------------------------

static float[] ReadF32(string path)
{
    var b = File.ReadAllBytes(path);
    var x = new float[b.Length / 4];
    Buffer.BlockCopy(b, 0, x, 0, x.Length * 4);
    return x;
}

static double[] ReadF64(string path)
{
    var b = File.ReadAllBytes(path);
    var x = new double[b.Length / 8];
    Buffer.BlockCopy(b, 0, x, 0, x.Length * 8);
    return x;
}

static void WriteF32(string path, float[] x)
{
    var b = new byte[x.Length * 4];
    Buffer.BlockCopy(x, 0, b, 0, b.Length);
    File.WriteAllBytes(path, b);
}

static void WriteF64(string path, double[] x)
{
    var b = new byte[x.Length * 8];
    Buffer.BlockCopy(x, 0, b, 0, b.Length);
    File.WriteAllBytes(path, b);
}
