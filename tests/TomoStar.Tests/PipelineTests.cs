// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using System.Text.Json.Nodes;
using TomoStar.Cli;
using TomoStar.Core.Attenuation;
using TomoStar.Core.Geo;
using TomoStar.Core.IO;
using TomoStar.Core.Location;
using TomoStar.Core.Model;
using TomoStar.Core.Signal;
using TomoStar.Core.Tomography;

namespace TomoStar.Tests;

/// <summary>
/// The processing steps on a small synthetic data set with a known answer: every step must move the
/// result towards the truth. One data set is shared by the tests of the class (made once).
/// </summary>
public class PipelineTests : IClassFixture<PipelineTests.SmallDataSet>
{
    /// <summary>12 stations, 30 events, waveforms, on a 5 km grid.</summary>
    public sealed class SmallDataSet : IDisposable
    {
        public string Folder { get; } = Directory.CreateTempSubdirectory("tomostar-pipeline-").FullName;
        public SyntheticResult Result { get; }
        public SyntheticOptions Options { get; } = new()
        {
            Stations = 12, Events = 30, HorizontalSpacingKm = 5, VerticalSpacingKm = 2.5, BottomKm = 20, Waveforms = true,
            MinMagnitude = 2.5, MaxMagnitude = 3.5, NoiseTStar = 0.0005, Seed = 99
        };

        public SmallDataSet() => Result = SyntheticData.Generate(Options, Path.Combine(Folder, "data"));

        public void Dispose()
        {
            try { Directory.Delete(Folder, true); } catch (IOException) { }
        }
    }

    private readonly SmallDataSet _set;

    public PipelineTests(SmallDataSet set) => _set = set;

    private string Data => Path.Combine(_set.Folder, "data");

    private static double MedianHorizontalError(IEnumerable<EventRecord> events, IReadOnlyDictionary<string, EventRecord> truth) =>
        Median(events.Select(e => GeoMath.SurfaceDistanceKm(e.Lon, e.Lat, truth[e.Id].Lon, truth[e.Id].Lat)));

    private static double Median(IEnumerable<double> v)
    {
        var a = v.OrderBy(x => x).ToArray();
        return a[a.Length / 2];
    }

    [Fact]
    public void SyntheticDataSetIsComplete()
    {
        foreach (var f in new[] { "stations.csv", "events.csv", "events_true.csv", "picks.csv", "tstar.csv", "grid.json", "model1d.txt", "truth/true_Vp.qvol" })
            Assert.True(File.Exists(Path.Combine(Data, f)), f);
        Assert.Equal(30, Directory.GetDirectories(Path.Combine(Data, "waveforms")).Length);
        var c = CatalogueReader.ReadCsv(Path.Combine(Data, "stations.csv"), Path.Combine(Data, "events.csv"), Path.Combine(Data, "picks.csv"), Path.Combine(Data, "tstar.csv"));
        Assert.Equal(12 * 30 * 2, c.Events.Sum(e => e.Picks.Count));
        Assert.Equal(12 * 30, c.TStar.Count);
    }

    [Fact]
    public void AutomaticPicksFindTheOnsets()
    {
        // Every P pick of the first events against the travel time of the true model (the synthetic
        // picks carry noise, the waveforms are what the picker sees).
        var archive = WaveformArchive.FromFolder(Path.Combine(Data, "waveforms"));
        var c = CatalogueReader.ReadCsv(Path.Combine(Data, "stations.csv"), Path.Combine(Data, "events_true.csv"), Path.Combine(Data, "picks.csv"));
        var errors = new List<double>();
        foreach (var ev in c.Events.Take(5))
        {
            var reference = ev.Picks.Where(p => p.Phase == Phase.P).ToDictionary(p => p.StationId, p => p.Time);
            foreach (var (station, traces) in archive.ForEvent(ev))
            {
                var picks = AutoPicker.PickStation(traces, reference[station], null, new AutoPickerSettings { PickS = false });
                var p = Assert.Single(picks);
                errors.Add(Math.Abs((p.Time - reference[station]).TotalSeconds));
            }
        }
        Assert.Equal(5 * 12, errors.Count);
        // The reference picks have 30 ms of noise themselves.
        Assert.True(Median(errors) < 0.06, $"median pick error {Median(errors):0.000} s");
    }

    [Fact]
    public void RelocationMovesEventsTowardsTheTruth()
    {
        var c = CatalogueReader.ReadCsv(Path.Combine(Data, "stations.csv"), Path.Combine(Data, "events.csv"), Path.Combine(Data, "picks.csv"));
        var truth = new Catalogue();
        CatalogueReader.ReadEvents(Path.Combine(Data, "events_true.csv"), truth);
        var byId = truth.EventIndex();
        var grid = new SphericalGrid(_set.Result.Grid);
        var model = CatalogueReader.ReadModel1D(Path.Combine(Data, "model1d.txt"));
        var (vp, vs) = TravelTimeTomography.StartingModel(grid, model);
        var located = Relocation.Absolute(c, grid, vp, vs, 2, new LocatorSettings(), Path.Combine(_set.Folder, "loc"), false);
        Assert.All(located, r => Assert.True(r.Ok, r.Message));
        var before = MedianHorizontalError(c.Events, byId);
        var after = MedianHorizontalError(Relocation.Apply(c, located).Events, byId);
        Assert.True(after < 0.5 * before, $"median horizontal error {before:0.00} km before, {after:0.00} km after");
    }

    [Fact]
    public void TomographyRecoversTheCheckerboardAndReducesTheMisfit()
    {
        var c = CatalogueReader.ReadCsv(Path.Combine(Data, "stations.csv"), Path.Combine(Data, "events.csv"), Path.Combine(Data, "picks.csv"));
        var grid = new SphericalGrid(_set.Result.Grid);
        var model = CatalogueReader.ReadModel1D(Path.Combine(Data, "model1d.txt"));
        var (vp, vs) = TravelTimeTomography.StartingModel(grid, model);
        var s = new TomographySettings { UseOpenCl = false, Iterations = 4, StationCorrections = false };
        var data = ObservationSet.FromCatalogue(c, grid, true, true);
        var r = new TravelTimeTomography(grid, s, Path.Combine(_set.Folder, "tomo")) { Background = model }.Run(data, vp, vs);
        Assert.True(r.Iterations[^1].Rms < 0.5 * r.Iterations[0].Rms, $"RMS {r.Iterations[0].Rms:0.000} to {r.Iterations[^1].Rms:0.000} s");
        // Correlation of the recovered and true dVp over the well-sampled nodes.
        var truth = _set.Result.TrueVp.Select((v, i) => v / vp[i] - 1).ToArray();
        var rec = r.Vp.Select((v, i) => v / vp[i] - 1).ToArray();
        var med = Median(r.DwsP.Where(x => x > 0));
        var idx = Enumerable.Range(0, rec.Length).Where(i => r.DwsP[i] >= med).ToArray();
        Assert.True(Correlation(truth, rec, idx) > 0.5, $"correlation {Correlation(truth, rec, idx):0.00}");
    }

    [Fact]
    public void CheckerboardTestRecoversItsPattern()
    {
        var c = CatalogueReader.ReadCsv(Path.Combine(Data, "stations.csv"), Path.Combine(Data, "events_true.csv"), Path.Combine(Data, "picks.csv"));
        var grid = new SphericalGrid(_set.Result.Grid);
        var model = CatalogueReader.ReadModel1D(Path.Combine(Data, "model1d.txt"));
        var (vp, vs) = TravelTimeTomography.StartingModel(grid, model);
        var geometry = ObservationSet.FromCatalogue(c, grid, true, true);
        var pattern = new ResolutionPattern { IsCheckerboard = true, CellHorizontalKm = 16, CellVerticalKm = 10, AmplitudePercent = 5 };
        var t = ResolutionTests.Velocity(grid, new TomographySettings { UseOpenCl = false, Iterations = 3 }, geometry, vp, vs, pattern, Phase.P, 0.02, 1,
            Path.Combine(_set.Folder, "chk"), background: model);
        Assert.True(t.CorrelationWellSampled > 0.6, $"correlation {t.CorrelationWellSampled:0.00}");
    }

    [Fact]
    public void TStarFromWaveformsMatchesThePathIntegral()
    {
        var c = CatalogueReader.ReadCsv(Path.Combine(Data, "stations.csv"), Path.Combine(Data, "events_true.csv"), Path.Combine(Data, "picks.csv"), Path.Combine(Data, "tstar.csv"));
        var exact = c.TStar.ToDictionary(t => (t.EventId, t.StationId), t => t.TStar);
        var archive = WaveformArchive.FromFolder(Path.Combine(Data, "waveforms"));
        var stations = c.StationIndex();
        var diffs = new List<double>();
        foreach (var ev in c.Events)
        {
            var (results, _, _) = TStarEstimator.MeasureEvent(ev, archive.ForEvent(ev), new TStarSettings(), null, id => stations.GetValueOrDefault(id));
            diffs.AddRange(results.Where(m => !m.Disabled).Select(m => m.TStar - exact[(m.EventId, m.StationId)]));
        }
        Assert.True(diffs.Count > 50, $"{diffs.Count} resolved t*");
        Assert.True(Math.Abs(Median(diffs)) < 0.003, $"median t* bias {Median(diffs) * 1000:0.0} ms");
    }

    [Fact]
    public void QTomographyRecoversTheReferenceQ()
    {
        var c = CatalogueReader.ReadCsv(Path.Combine(Data, "stations.csv"), Path.Combine(Data, "events_true.csv"), null, Path.Combine(Data, "tstar.csv"));
        var grid = new SphericalGrid(_set.Result.Grid);
        var data = QTomography.FromCatalogue(c, grid);
        var q = new QTomography(grid, new Core.Attenuation.QTomographySettings { UseOpenCl = false, StationTerms = false }, Path.Combine(_set.Folder, "q"))
            .Run(data, _set.Result.TrueVp, _set.Result.TrueVs);
        Assert.InRange(q.ReferenceQ, 0.85 * _set.Options.Q0, 1.15 * _set.Options.Q0);
        Assert.True(q.RmsAfter < q.RmsBefore);
    }

    private static double Correlation(double[] a, double[] b, int[] idx)
    {
        var ma = idx.Average(i => a[i]);
        var mb = idx.Average(i => b[i]);
        double sab = 0, saa = 0, sbb = 0;
        foreach (var i in idx)
        {
            sab += (a[i] - ma) * (b[i] - mb);
            saa += (a[i] - ma) * (a[i] - ma);
            sbb += (b[i] - mb) * (b[i] - mb);
        }
        return sab / Math.Sqrt(saa * sbb);
    }

    [Fact]
    public void CommandLineRunsTheLCurveAndWritesTheChoiceIntoTheConfiguration()
    {
        var outDir = Path.Combine(_set.Folder, "lcurve");
        var cfg = Path.Combine(_set.Folder, "cfg.json");
        var code = Program.Execute(["lcurve", Data, "--grid", Path.Combine(Data, "grid.json"), "--model", Path.Combine(Data, "model1d.txt"),
            "--hypocentres", Path.Combine(Data, "events_true.csv"), "--dampings", "5,20,80", "--smoothings", "20",
            "--out", outDir, "--update-config", cfg, "--no-opencl", "--quiet"], CancellationToken.None);
        Assert.Equal(0, code);
        Assert.True(File.Exists(Path.Combine(outDir, "lcurve.csv")));
        Assert.True(File.Exists(Path.Combine(outDir, "lcurve.svg")));
        var written = JsonNode.Parse(File.ReadAllText(cfg))!;
        Assert.Equal(20, written["Tomography"]!["Smoothing"]!.GetValue<double>());
    }

    /// <summary>The library example of the README, as written there.</summary>
    [Fact]
    public void ReadmeLibraryExampleRuns()
    {
        var previous = Environment.CurrentDirectory;
        Environment.CurrentDirectory = Data;
        try
        {
            var data  = CatalogueReader.ReadCsv("stations.csv", "events.csv", "picks.csv");
            var grid  = new SphericalGrid(CatalogueReader.ReadGrid("grid.json"));
            var model = CatalogueReader.ReadModel1D("ak135");
            var (vp, vs) = TravelTimeTomography.StartingModel(grid, model);
            var obs   = ObservationSet.FromCatalogue(data, grid, includeP: true, includeS: true);
            var run   = new TravelTimeTomography(grid, new TomographySettings { UseOpenCl = false, Iterations = 1 }, "work", Console.WriteLine) { Background = model }
                        .Run(obs, vp, vs);
            RunWriter.SaveVelocity("runs/vel", "vel", run, new TomographySettings(), data, model, []);
            Assert.True(File.Exists("runs/vel/volumes/vel_Vp.qvol"));
        }
        finally
        {
            Environment.CurrentDirectory = previous;
        }
    }
}
