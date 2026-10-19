// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using TomoStar.Core.IO;
using TomoStar.Core.Location;
using TomoStar.Core.Model;
using TomoStar.Core.Source;

namespace TomoStar.Tests;

/// <summary>Source studies against known answers: mechanisms, stress, Coulomb stress change, cross-correlation.</summary>
public class SourceTests
{
    private static readonly DateTime T0 = new(2020, 1, 1, 0, 0, 0, DateTimeKind.Utc);

    [Fact]
    public void FirstMotionRecoversAKnownMechanismDespiteWrongPolarities()
    {
        var truth = new NodalPlane(160, 80, -10);
        var m = MechanismGeometry.Tensor(truth);
        var rnd = new Random(5);
        var obs = new List<PolarityObservation>();
        for (var k = 0; k < 120; k++)
        {
            double az = 360 * rnd.NextDouble(), to = 20 + 140 * rnd.NextDouble();
            var o = new PolarityObservation($"S{k}", az, to, 1);
            var (n, e, d) = o.Direction;
            var pol = Math.Sign(m.Radiation(n, e, d));
            if (rnd.NextDouble() < 0.05) pol = -pol;
            obs.Add(o with { Polarity = pol });
        }
        var r = FirstMotionInversion.Invert(obs, new FirstMotionSettings { GridStepDeg = 5 });
        Assert.NotNull(r);
        var angle = MechanismGeometry.KaganAngle(r!.Plane1, truth);
        Assert.True(angle < 15, $"Kagan angle {angle:0.0}°");
        Assert.True(r.MisfitFraction < 0.1, $"misfit {r.MisfitFraction:0.00}");
    }

    [Fact]
    public void StressInversionRecoversAStrikeSlipTensor()
    {
        // σ1 horizontal towards N130°E, σ3 horizontal towards N40°E, σ2 vertical, R = 0.5; slip along the
        // resolved shear traction on random planes (Wallace-Bott), one nodal plane given as the fault.
        double a1 = 130 * Math.PI / 180, a3 = 40 * Math.PI / 180;
        double[] e1 = [Math.Sin(a1), Math.Cos(a1), 0], e3 = [Math.Sin(a3), Math.Cos(a3), 0], e2 = [0, 0, 1]; // east, north, down
        double s1 = 1, s2 = 0.5, s3 = 0;
        var S = new double[3, 3];
        for (var i = 0; i < 3; i++)
        for (var j = 0; j < 3; j++)
            S[i, j] = s1 * e1[i] * e1[j] + s2 * e2[i] * e2[j] + s3 * e3[i] * e3[j];
        var rnd = new Random(11);
        var obs = new List<FocalMechanismObservation>();
        while (obs.Count < 60)
        {
            double strike = 360 * rnd.NextDouble(), dip = 45 + 45 * rnd.NextDouble();
            double sr = strike * Math.PI / 180, dr = dip * Math.PI / 180;
            // normal (east, north, down) of the plane, pointing up into the hanging wall
            double[] n = [Math.Cos(sr) * Math.Sin(dr), -Math.Sin(sr) * Math.Sin(dr), -Math.Cos(dr)];
            var t = new double[3];
            for (var i = 0; i < 3; i++) t[i] = S[i, 0] * n[0] + S[i, 1] * n[1] + S[i, 2] * n[2];
            var tn = t[0] * n[0] + t[1] * n[1] + t[2] * n[2];
            double[] shear = [t[0] - tn * n[0], t[1] - tn * n[1], t[2] - tn * n[2]];
            var len = Math.Sqrt(shear.Sum(x => x * x));
            if (len < 0.2) continue;
            // Compression positive: the hanging wall slips opposite to the shear traction it feels.
            double[] slip = [-shear[0] / len, -shear[1] / len, -shear[2] / len];
            double[] along = [Math.Sin(sr), Math.Cos(sr), 0];
            double[] down = [Math.Cos(sr) * Math.Cos(dr), -Math.Sin(sr) * Math.Cos(dr), Math.Sin(dr)];
            var rake = Math.Atan2(-(slip[0] * down[0] + slip[1] * down[1] + slip[2] * down[2]), slip[0] * along[0] + slip[1] * along[1]) * 180 / Math.PI;
            var p1 = new NodalPlane(strike, dip, rake);
            var p2 = MechanismGeometry.Auxiliary(p1);
            obs.Add(new FocalMechanismObservation
            {
                EventId = $"e{obs.Count}", LatitudeDeg = 46.9, LongitudeDeg = 8.9, DepthKm = 5,
                NodalPlane1StrikeDeg = p1.Strike, NodalPlane1DipDeg = p1.Dip, NodalPlane1RakeDeg = p1.Rake,
                NodalPlane2StrikeDeg = p2.Strike, NodalPlane2DipDeg = p2.Dip, NodalPlane2RakeDeg = p2.Rake,
                Provider = "synthetic", SourceId = "synthetic", Provenance = "test"
            });
        }
        var res = new FocalMechanismStressInversionService().Invert(obs, new StressInversionOptions { BootstrapRealizations = 100 });
        Assert.NotNull(res.Tensor);
        static double AxisDiff(double az, double target) { var d = Math.Abs(((az - target) % 180 + 180) % 180); return Math.Min(d, 180 - d); }
        Assert.True(AxisDiff(res.Tensor!.Sigma1.AzimuthDeg, 130) < 10, $"σ1 azimuth {res.Tensor.Sigma1.AzimuthDeg:0}");
        Assert.True(AxisDiff(res.Tensor.Sigma3.AzimuthDeg, 40) < 10, $"σ3 azimuth {res.Tensor.Sigma3.AzimuthDeg:0}");
        Assert.True(Math.Abs(StressField.SHmaxAzimuth(res.Tensor) - 130) < 10, $"S_Hmax {StressField.SHmaxAzimuth(res.Tensor):0}");
        Assert.Equal("SS", StressField.Regime(res.Tensor));
    }

    [Fact]
    public void CoulombStressLoadsTheTipsAndUnloadsTheSidesOfAStrikeSlipFault()
    {
        // A N-S sinistral fault, 2 km long, at 8 km; receivers parallel to it.
        var patches = CoulombStress.TaperedSource(46.0, 8.0, 8.0, 0, 90, 0, 5.0, 2.0, 2.0, 9, 9);
        var kmLat = 1 / 111.195;
        var rec = new List<CoulombStress.Receiver>
        {
            new("tip", 46.0 + 1.6 * kmLat, 8.0, 8.0, 0, 90, 0),
            new("side", 46.0, 8.0 + 0.5 / (111.195 * Math.Cos(46 * Math.PI / 180)), 8.0, 0, 90, 0)
        };
        var r = CoulombStress.Compute(patches, rec, 0.4);
        Assert.True(r[0].DeltaCfsPa > 1e4, $"tip ΔCFS {r[0].DeltaCfsPa / 1e5:0.00} bar");
        Assert.True(r[1].DeltaCfsPa < -1e4, $"side ΔCFS {r[1].DeltaCfsPa / 1e5:0.00} bar");
        Assert.Equal(Math.Pow(10, 1.5 * 5.0 + 9.1), CoulombStress.Moment(patches), 1e10);
    }

    [Fact]
    public void CoulombStressVanishesOnTheFreeSurfaceTractions()
    {
        var patches = CoulombStress.TaperedSource(46.0, 8.0, 3.0, 30, 60, 90, 5.0, 2.0, 2.0, 5, 5);
        var t = CoulombStress.Tensors(patches, [(46.01, 8.02, 0.0)])[0].Stress;
        var scale = Math.Max(Math.Abs(t.EePa), Math.Abs(t.NnPa));
        Assert.True(scale > 1e3);
        Assert.True(Math.Abs(t.UuPa) < 1e-6 * scale && Math.Abs(t.EuPa) < 1e-6 * scale && Math.Abs(t.NuPa) < 1e-6 * scale,
            $"surface tractions {t.UuPa:0.###}, {t.EuPa:0.###}, {t.NuPa:0.###} Pa against {scale:0} Pa");
    }

    [Fact]
    public void CrossCorrelationRecoversASubsampleDelay()
    {
        const double fs = 100;
        double Wavelet(double t) => t < 0 ? 0 : (1 - Math.Exp(-t / 0.03)) * Math.Exp(-t / 0.6) *
                                                 (Math.Sin(2 * Math.PI * 6 * t) + 0.6 * Math.Sin(2 * Math.PI * 9.3 * t + 1) + 0.4 * Math.Sin(2 * Math.PI * 3.1 * t + 2));
        Trace Record(double arrival, int seed)
        {
            var rnd = new Random(seed);
            var data = new float[3000];
            for (var i = 0; i < data.Length; i++) data[i] = (float)(Wavelet(i / fs - arrival) + 0.01 * (rnd.NextDouble() - 0.5));
            return new Trace { Network = "XX", Station = "A", Channel = "HHZ", StartTime = T0, SampleRate = fs, Data = data };
        }
        double[] arrivals = [10.0, 10.2374];
        double[] pickErrors = [0.03, -0.04];
        var events = new List<EventRecord>();
        var records = new List<(string, Trace)>();
        for (var e = 0; e < 2; e++)
        {
            var ev = new EventRecord { Id = $"e{e}", OriginTime = T0, Lon = 8.0, Lat = 46.0 + 0.001 * e, DepthKm = 5 };
            ev.Picks.Add(new PickRecord { StationId = "XX.A", Phase = Phase.P, Time = T0.AddSeconds(arrivals[e] + pickErrors[e]), Sigma = 0.05, Origin = PickOrigin.Manual });
            events.Add(ev);
            records.Add((ev.Id, Record(arrivals[e], e + 1)));
        }
        var stations = new Dictionary<string, StationRecord>(StringComparer.OrdinalIgnoreCase) { ["XX.A"] = new() { Id = "XX.A", Lon = 8.1, Lat = 46.1 } };
        var list = WaveformCrossCorrelation.Measure(events, stations, _ => records, new CrossCorrelationSettings { IncludeS = false, MaxLagSeconds = 0.2 });
        Assert.Single(list);
        Assert.True(Math.Abs(list[0].ArrivalDifference - (arrivals[0] - arrivals[1])) < 0.002, $"Δt {list[0].ArrivalDifference:0.0000} s");
        Assert.True(list[0].Coefficient > 0.95);
    }
}
