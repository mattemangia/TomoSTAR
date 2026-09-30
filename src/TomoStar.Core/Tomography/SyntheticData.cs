using System.Numerics;
using TomoStar.Core.Forward;
using TomoStar.Core.Geo;
using TomoStar.Core.IO;
using TomoStar.Core.Model;
using TomoStar.Core.Signal;

namespace TomoStar.Core.Tomography;

/// <summary>What a synthetic data set contains and how it is made.</summary>
public sealed class SyntheticOptions
{
    public double MinLon { get; set; } = 13.0;
    public double MaxLon { get; set; } = 13.6;
    public double MinLat { get; set; } = 42.6;
    public double MaxLat { get; set; } = 43.1;

    /// <summary>Node spacing of the grid the true model lives on, km.</summary>
    public double HorizontalSpacingKm { get; set; } = 4;

    public double VerticalSpacingKm { get; set; } = 2;
    public double TopKm { get; set; } = -2;
    public double BottomKm { get; set; } = 24;

    public int Stations { get; set; } = 25;
    public int Events { get; set; } = 150;
    public double MinEventDepthKm { get; set; } = 2;
    public double MaxEventDepthKm { get; set; } = 16;
    public double MinMagnitude { get; set; } = 1.5;
    public double MaxMagnitude { get; set; } = 3.5;
    public int Seed { get; set; } = 12345;

    /// <summary>1-D background (a name from the model library or a file).</summary>
    public string Model { get; set; } = "ak135";

    /// <summary>Vp anomaly: a checkerboard of this cell size (km) and amplitude (%); 0 amplitude for none.</summary>
    public double CheckerCellKm { get; set; } = 16;

    public double CheckerDepthCellKm { get; set; } = 10;
    public double VpAnomalyPercent { get; set; } = 5;

    /// <summary>Vs anomaly, the same pattern (so Vp/Vs varies when the two amplitudes differ).</summary>
    public double VsAnomalyPercent { get; set; } = 5;

    /// <summary>Uniform background Qp, and the anomaly of 1/Qp in % (same pattern, opposite sign to Vp).</summary>
    public double Q0 { get; set; } = 300;

    public double QAnomalyPercent { get; set; } = 30;

    /// <summary>Qs / Qp, for the S waveforms.</summary>
    public double QsOverQp { get; set; } = 0.6;

    /// <summary>Gaussian noise of the picks and of the t*, s.</summary>
    public double NoiseP { get; set; } = 0.03;

    public double NoiseS { get; set; } = 0.06;
    public double NoiseTStar { get; set; } = 0.003;

    /// <summary>Error of the catalogue hypocentres given as the starting locations (km, and s for the origin time).</summary>
    public double LocationErrorKm { get; set; } = 1.5;

    public double OriginTimeErrorS { get; set; } = 0.2;

    /// <summary>Also write three-component waveforms (miniSEED) in which the picks and t* can be measured.</summary>
    public bool Waveforms { get; set; }

    public double SampleRateHz { get; set; } = 100;

    /// <summary>Peak P signal over the RMS of the noise on the vertical.</summary>
    public double WaveformSnr { get; set; } = 200;

    /// <summary>Brune stress drop of the synthetic sources, MPa.</summary>
    public double StressDropMPa { get; set; } = 3;
}

/// <summary>The files and the truth of a synthetic data set.</summary>
public sealed record SyntheticResult(Catalogue Catalogue, GridDefinition Grid, double[] TrueVp, double[] TrueVs, double[] TrueQp,
    List<EventRecord> TrueEvents, int Arrivals, int TStar);

/// <summary>
/// Synthetic local-earthquake data sets: stations and earthquakes in a box, a true 3-D model (a 1-D
/// background with a smooth checkerboard in Vp, Vs and 1/Q), travel times and t* along the curved
/// fast-marching rays of that model with Gaussian noise, starting hypocentres perturbed from the true
/// ones, and optionally three-component waveforms (Brune omega-squared sources attenuated with the
/// t* of each path, P on the vertical, S on the horizontals, white noise). They exercise every step
/// of the pipeline, picking included, against a known answer, and are the example data of TomoSTAR.
/// </summary>
public static class SyntheticData
{
    public static SyntheticResult Generate(SyntheticOptions o, string outputFolder, Action<string>? log = null, CancellationToken ct = default)
    {
        Directory.CreateDirectory(outputFolder);
        var rnd = new Random(o.Seed);
        var model = CatalogueReader.ReadModel1D(o.Model);
        var def = GridDefinition.FromSpacing(o.MinLon, o.MaxLon, o.MinLat, o.MaxLat, o.TopKm, o.BottomKm, o.HorizontalSpacingKm, o.VerticalSpacingKm);
        var grid = new SphericalGrid(def);
        var (vp0, vs0) = TravelTimeTomography.StartingModel(grid, model);
        var pattern = new ResolutionPattern { IsCheckerboard = true, CellHorizontalKm = o.CheckerCellKm, CellVerticalKm = o.CheckerDepthCellKm, AmplitudePercent = 100 };
        var shape = pattern.Evaluate(grid);
        var vp = vp0.Select((v, i) => v * (1 + o.VpAnomalyPercent / 100 * shape[i])).ToArray();
        var vs = vs0.Select((v, i) => v * (1 + o.VsAnomalyPercent / 100 * shape[i])).ToArray();
        // Slow rock attenuates more: the 1/Q anomaly has the opposite sign of the velocity anomaly.
        var q = shape.Select(x => 1 / o.Q0 * (1 - o.QAnomalyPercent / 100 * x)).ToArray();

        // Stations on a jittered regular pattern over the box, events uniformly inside it.
        var c = new Catalogue();
        var nSide = (int)Math.Ceiling(Math.Sqrt(o.Stations));
        for (var k = 0; k < o.Stations; k++)
        {
            var u = (k % nSide + 0.5 + 0.6 * (rnd.NextDouble() - 0.5)) / nSide;
            var v = (k / nSide + 0.5 + 0.6 * (rnd.NextDouble() - 0.5)) / nSide;
            c.Stations.Add(new StationRecord
            {
                Id = $"SY.S{k + 1:D3}", Lon = Lerp(o.MinLon, o.MaxLon, 0.05 + 0.9 * u), Lat = Lerp(o.MinLat, o.MaxLat, 0.05 + 0.9 * v),
                ElevationM = Math.Round(rnd.NextDouble() * 1200)
            });
        }
        var t0 = new DateTime(2026, 1, 1, 0, 0, 0, DateTimeKind.Utc);
        var truth = new List<EventRecord>();
        for (var k = 0; k < o.Events; k++)
        {
            truth.Add(new EventRecord
            {
                Id = $"ev{k + 1:D4}", OriginTime = t0.AddMinutes(17 * k).AddSeconds(Math.Round(rnd.NextDouble() * 60, 2)),
                Lon = Lerp(o.MinLon, o.MaxLon, 0.1 + 0.8 * rnd.NextDouble()), Lat = Lerp(o.MinLat, o.MaxLat, 0.1 + 0.8 * rnd.NextDouble()),
                DepthKm = Lerp(o.MinEventDepthKm, o.MaxEventDepthKm, rnd.NextDouble()),
                Magnitude = Math.Round(Lerp(o.MinMagnitude, o.MaxMagnitude, Math.Pow(rnd.NextDouble(), 2)), 1)
            });
        }

        // Travel times and ray kernels in the true model, every event at every station, P and S.
        var set = new ObservationSet();
        set.Stations.AddRange(c.Stations.Select(s => new StationState { Id = s.Id, Lon = s.Lon, Lat = s.Lat, DepthKm = Math.Max(def.MinDepthKm, s.DepthKm) }));
        set.Events.AddRange(truth.Select(e => new EventState { Id = e.Id, Lon = e.Lon, Lat = e.Lat, DepthKm = e.DepthKm, Reference = e.OriginTime, Fixed = true }));
        for (var e = 0; e < truth.Count; e++)
        for (var s = 0; s < c.Stations.Count; s++)
        {
            set.Arrivals.Add(new Observation { Event = e, Station = s, Phase = Phase.P });
            set.Arrivals.Add(new Observation { Event = e, Station = s, Phase = Phase.S });
        }
        log?.Invoke($"Synthetic: {c.Stations.Count} stations, {truth.Count} events, grid {def.Nx}x{def.Ny}x{def.Nz}; tracing {set.Arrivals.Count} rays in the true model.");
        var sP = vp.Select(v => 1 / v).ToArray();
        var sS = vs.Select(v => 1 / v).ToArray();
        var rows = new RayKernels(grid, def.Refined(2), RayMethod.FastMarching, Path.Combine(outputFolder, ".work"), log)
            .Compute(set, sP, sS, null, 0, null, 0, ct);
        var tt = new double[set.Arrivals.Count];
        var ts = new double[set.Arrivals.Count];
        Array.Fill(tt, double.NaN);
        foreach (var r in rows)
        {
            tt[r.Arrival] = r.TCalc;
            var slow = set.Arrivals[r.Arrival].Phase == Phase.P ? sP : sS;
            var qf = set.Arrivals[r.Arrival].Phase == Phase.P ? 1.0 : 1 / o.QsOverQp;
            double t = 0;
            for (var k = 0; k < r.Nodes.Length; k++) t += r.Length[k] * slow[r.Nodes[k]] * q[r.Nodes[k]] * qf;
            ts[r.Arrival] = t;
        }
        try { Directory.Delete(Path.Combine(outputFolder, ".work"), true); } catch (IOException) { }

        // Catalogue: perturbed starting locations, noisy absolute picks, noisy P t*.
        var arrivals = 0;
        for (var e = 0; e < truth.Count; e++)
        {
            var te = truth[e];
            var (lon, lat) = GeoMath.Destination(te.Lon, te.Lat, 360 * rnd.NextDouble(), o.LocationErrorKm * Math.Abs(Gauss(rnd)));
            var ev = new EventRecord
            {
                Id = te.Id, Magnitude = te.Magnitude,
                OriginTime = te.OriginTime.AddSeconds(o.OriginTimeErrorS * Gauss(rnd)),
                Lon = lon, Lat = lat, DepthKm = Math.Clamp(te.DepthKm + o.LocationErrorKm * Gauss(rnd), 0, def.MaxDepthKm - 1)
            };
            for (var s = 0; s < c.Stations.Count; s++)
            {
                foreach (var ph in new[] { Phase.P, Phase.S })
                {
                    var ai = 2 * (e * c.Stations.Count + s) + (ph == Phase.P ? 0 : 1);
                    if (!double.IsFinite(tt[ai])) continue;
                    var noise = ph == Phase.P ? o.NoiseP : o.NoiseS;
                    ev.Picks.Add(new PickRecord
                    {
                        StationId = c.Stations[s].Id, Phase = ph, Sigma = Math.Max(0.01, noise),
                        Time = te.OriginTime.AddTicks((long)Math.Round((tt[ai] + noise * Gauss(rnd)) * TimeSpan.TicksPerSecond))
                    });
                    arrivals++;
                    if (ph == Phase.P)
                        c.TStar.Add(new TStarRecord
                        {
                            EventId = te.Id, StationId = c.Stations[s].Id, Phase = Phase.P,
                            TStar = ts[ai] + o.NoiseTStar * Gauss(rnd), Sigma = Math.Max(0.001, o.NoiseTStar)
                        });
                }
            }
            c.Events.Add(ev);
        }

        WriteFiles(o, outputFolder, c, truth, def, grid, model, vp, vs, q);
        if (o.Waveforms)
        {
            log?.Invoke("Synthetic: writing three-component waveforms.");
            WriteWaveforms(o, Path.Combine(outputFolder, "waveforms"), c, truth, tt, ts, ct);
        }
        return new SyntheticResult(c, def, vp, vs, q.Select(x => 1 / x).ToArray(), truth, arrivals, c.TStar.Count);
    }

    private static void WriteFiles(SyntheticOptions o, string folder, Catalogue c, List<EventRecord> truth, GridDefinition def,
        SphericalGrid grid, VelocityModel1D model, double[] vp, double[] vs, double[] q)
    {
        var inv = System.Globalization.CultureInfo.InvariantCulture;
        File.WriteAllLines(Path.Combine(folder, "stations.csv"),
            ["id,lon,lat,elevation_m", .. c.Stations.Select(s => string.Create(inv, $"{s.Id},{s.Lon:0.#####},{s.Lat:0.#####},{s.ElevationM:0}"))]);
        static string EventLine(EventRecord e) => string.Create(System.Globalization.CultureInfo.InvariantCulture,
            $"{e.Id},{e.OriginTime:yyyy-MM-ddTHH:mm:ss.fffZ},{e.Lon:0.#####},{e.Lat:0.#####},{e.DepthKm:0.###},{e.Magnitude:0.0}");
        File.WriteAllLines(Path.Combine(folder, "events.csv"), ["id,time,lon,lat,depth_km,magnitude", .. c.Events.Select(EventLine)]);
        File.WriteAllLines(Path.Combine(folder, "events_true.csv"), ["id,time,lon,lat,depth_km,magnitude", .. truth.Select(EventLine)]);
        File.WriteAllLines(Path.Combine(folder, "picks.csv"),
            ["event,station,phase,time,sigma_s", .. c.Events.SelectMany(e => e.Picks.Select(p =>
                string.Create(inv, $"{e.Id},{p.StationId},{p.Phase},{p.Time:yyyy-MM-ddTHH:mm:ss.ffffZ},{p.Sigma:0.###}")))]);
        File.WriteAllLines(Path.Combine(folder, "tstar.csv"),
            ["event,station,phase,tstar_s,sigma_s", .. c.TStar.Select(t => string.Create(inv, $"{t.EventId},{t.StationId},{t.Phase},{t.TStar:0.######},{t.Sigma:0.####}"))]);
        File.WriteAllText(Path.Combine(folder, "model1d.txt"), $"# {model.Name}\n# depth_km vp_km_s vs_km_s\n" + model.Format() + "\n");
        File.WriteAllText(Path.Combine(folder, "grid.json"), System.Text.Json.JsonSerializer.Serialize(def, TomoJson.Options));
        File.WriteAllText(Path.Combine(folder, "synthetic.json"), System.Text.Json.JsonSerializer.Serialize(o, TomoJson.Options));
        var truthDir = Path.Combine(folder, "truth");
        Directory.CreateDirectory(truthDir);
        void Vol(string quantity, string units, double[] values) => VolumeFile.Write(Path.Combine(truthDir, $"true_{quantity}.qvol"),
            new VolumeHeader { Name = $"true {quantity}", Quantity = quantity, Units = units, Grid = def, Source = "TomoSTAR synthetic" }, values);
        Vol("Vp", "km/s", vp);
        Vol("Vs", "km/s", vs);
        Vol("VpVs", "", vp.Zip(vs, (a, b) => a / b).ToArray());
        Vol("Qp", "", q.Select(x => 1 / x).ToArray());
        var (vp0, vs0) = TravelTimeTomography.StartingModel(grid, model);
        Vol("dVp", "%", vp.Select((v, i) => 100 * (v / vp0[i] - 1)).ToArray());
        Vol("dVs", "%", vs.Select((v, i) => 100 * (v / vs0[i] - 1)).ToArray());
    }

    /// <summary>
    /// Three components per station and event. Displacement of each phase: the causal Brune pulse
    /// Omega0 a^2 t exp(-a t) H(t), a = 2 pi fc, whose amplitude spectrum is Omega0 / (1 + (f/fc)^2),
    /// with fc from the magnitude and the stress drop, attenuated by exp(-pi f t*) and delayed to the
    /// arrival; recorded as ground velocity (times 2 pi i f), plus white noise. The attenuation operator
    /// is zero phase (no dispersion), so the onset stays at the ray travel time within about t*/2.
    /// </summary>
    private static void WriteWaveforms(SyntheticOptions o, string folder, Catalogue c, List<EventRecord> truth, double[] tt, double[] ts, CancellationToken ct)
    {
        var fs = o.SampleRateHz;
        Parallel.For(0, truth.Count, Compute.ComputeSettings.Options(ct), e =>
        {
            var rnd = new Random(o.Seed * 7919 + e);
            var ev = truth[e];
            var m0 = Math.Pow(10, 1.5 * ev.Magnitude + 9.1);
            var fc = 0.37 * 3500 * Math.Cbrt(16 * o.StressDropMPa * 1e6 / (7 * m0));
            var traces = new List<Trace>();
            var before = 20.0;
            for (var s = 0; s < c.Stations.Count; s++)
            {
                var ai = 2 * (e * c.Stations.Count + s);
                if (!double.IsFinite(tt[ai]) || !double.IsFinite(tt[ai + 1])) continue;
                var st = c.Stations[s];
                var start = ev.OriginTime.AddSeconds(-before);
                var length = before + tt[ai + 1] + 25;
                var n = Fft.NextPow2((int)Math.Ceiling(length * fs));
                var dist = Math.Max(1, (GeoMath.ToCartesian(ev.Lon, ev.Lat, ev.DepthKm) - GeoMath.ToCartesian(st.Lon, st.Lat, st.DepthKm)).Length);
                var omega = m0 * 1e-8 / dist;
                var p = Pulse(n, fs, before + tt[ai], ts[ai], fc, omega);
                var sw = Pulse(n, fs, before + tt[ai + 1], ts[ai + 1], fc / 1.3, 3 * omega);
                var az = GeoMath.Azimuth(ev.Lon, ev.Lat, st.Lon, st.Lat) * GeoMath.Deg2Rad;
                var peak = p.Max(Math.Abs);
                var noise = peak / Math.Max(1, o.WaveformSnr);
                float[] Comp(double pw, double sw1) => Enumerable.Range(0, (int)Math.Ceiling(length * fs))
                    .Select(i => (float)(pw * p[i] + sw1 * sw[i] + noise * Gauss(rnd))).ToArray();
                var parts = st.Id.Split('.');
                Trace T(string ch, float[] d) => new()
                {
                    Network = parts[0], Station = parts[^1], Location = "", Channel = ch, StartTime = start, SampleRate = fs, Data = d, Units = "nm/s"
                };
                // P radial motion projects on N and E; S shared between the horizontals.
                traces.Add(T("HHZ", Comp(1, 0.2)));
                traces.Add(T("HHN", Comp(0.3 * Math.Cos(az), 0.8 * Math.Sin(az + 0.7))));
                traces.Add(T("HHE", Comp(0.3 * Math.Sin(az), 0.8 * Math.Cos(az + 0.7))));
            }
            MiniSeedWriter.Write(Path.Combine(folder, WaveformArchive.Safe(ev.Id), $"{ev.Id}.mseed"), traces);
        });
    }

    /// <summary>Ground velocity of one attenuated Brune pulse arriving at <paramref name="arrival"/> seconds.</summary>
    private static double[] Pulse(int n, double fs, double arrival, double tstar, double fc, double omega)
    {
        var a = new Complex[n];
        for (var k = 0; k <= n / 2; k++)
        {
            var f = k * fs / n;
            var w = 2 * Math.PI * f;
            var a0 = 2 * Math.PI * fc;
            // Causal source: the Fourier transform of a^2 t exp(-a t) H(t) is a^2 / (a + i w)^2.
            var source = omega * a0 * a0 / Complex.Pow(new Complex(a0, w), 2);
            var spec = new Complex(0, w) * source * Math.Exp(-Math.PI * f * tstar) * Complex.Exp(new Complex(0, -w * arrival));
            a[k] = spec;
            if (k > 0 && k < n / 2) a[n - k] = Complex.Conjugate(spec);
        }
        a[n / 2] = new Complex(a[n / 2].Real, 0);
        Fft.Transform(a, inverse: true);
        // Continuous spectrum to samples: the inverse DFT times the sampling rate.
        return a.Select(x => x.Real * fs).ToArray();
    }

    private static double Lerp(double a, double b, double t) => a + (b - a) * t;

    private static double Gauss(Random r) => Math.Sqrt(-2 * Math.Log(1 - r.NextDouble())) * Math.Cos(2 * Math.PI * r.NextDouble());
}
