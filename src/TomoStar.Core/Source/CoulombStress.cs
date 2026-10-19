// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

namespace TomoStar.Core.Source;

/// <summary>
/// Multi-plane static stress transfer from finite rectangular slip patches.
/// Each rectangle is evaluated with the analytic DC3D displacement-gradient
/// solution in a homogeneous, isotropic elastic half-space, including the
/// free-surface image terms of Okada (1992).
///
/// Stress is tensile-positive ENU in Pa. Receiver resolution follows
/// King, Stein &amp; Lin (1994): unclamping = n·Δσ·n, resolved shear =
/// s·Δσ·n, and ΔCFS = Δτ + μ'Δσn.
/// References: Okada (1992), DOI 10.1785/BSSA0820021018; King et al.
/// (1994), DOI 10.1785/BSSA0840030935.
/// </summary>
public static class CoulombStress
{
    public const string ModelName = "okada-1992-dc3d-half-space-v1";
    public const string Units = "Pa";
    public const string Confidence = "VERIFIED";
    public const string Provenance =
        "Okada (1992) DC3D finite rectangular dislocations in a homogeneous isotropic elastic half-space; " +
        "Coulomb resolution after King, Stein & Lin (1994), doi:10.1785/BSSA0840030935. " +
        "Okada (1992), doi:10.1785/BSSA0820021018.";

    /// <summary>
    /// One uniformly slipping rectangle: centre (lat, lon, depth below sea level), orientation
    /// (Aki &amp; Richards strike/dip/rake), size and slip.
    /// </summary>
    public readonly record struct SlipPatch(
        double Lat, double Lon, double DepthKm, double StrikeDeg, double DipDeg, double RakeDeg,
        double LengthKm, double WidthKm, double SlipM);

    /// <summary>
    /// A rectangular source centred on a hypocentre, divided into <paramref name="alongStrike"/> ×
    /// <paramref name="downDip"/> patches with a cosine taper of slip that vanishes at the edges,
    /// scaled so that the seismic moment equals that of <paramref name="mw"/>
    /// (M0 = 10^(1.5 Mw + 9.1) N m, Hanks &amp; Kanamori 1979, doi:10.1029/JB084iB05p02348).
    /// </summary>
    public static List<SlipPatch> TaperedSource(double lat, double lon, double depthKm, double strikeDeg, double dipDeg, double rakeDeg,
        double mw, double lengthKm, double widthKm, int alongStrike = 7, int downDip = 7, double shearModulusPa = 30e9)
    {
        double m0 = Math.Pow(10, 1.5 * mw + 9.1);
        double s = Rad(strikeDeg), d = Rad(dipDeg);
        double[] along = [Math.Sin(s), Math.Cos(s), 0];
        double[] down = [Math.Cos(s) * Math.Cos(d), -Math.Sin(s) * Math.Cos(d), -Math.Sin(d)];
        double dl = lengthKm / alongStrike, dw = widthKm / downDip;
        var taper = new double[alongStrike * downDip];
        for (int j = 0; j < downDip; j++)
        for (int i = 0; i < alongStrike; i++)
            taper[j * alongStrike + i] = Math.Cos(Math.PI * ((i + 0.5) / alongStrike - 0.5)) * Math.Cos(Math.PI * ((j + 0.5) / downDip - 0.5));
        double scale = m0 / (shearModulusPa * dl * dw * 1e6 * taper.Sum());
        var patches = new List<SlipPatch>();
        for (int j = 0; j < downDip; j++)
        for (int i = 0; i < alongStrike; i++)
        {
            double a = (-lengthKm / 2 + (i + 0.5) * dl) * 1000, b = (-widthKm / 2 + (j + 0.5) * dw) * 1000;
            var ll = ToLatLon(a * along[0] + b * down[0], a * along[1] + b * down[1], lat, lon);
            patches.Add(new SlipPatch(ll.Lat, ll.Lon, depthKm - b * down[2] / 1000, strikeDeg, dipDeg, rakeDeg, dl, dw, scale * taper[j * alongStrike + i]));
        }
        return patches;
    }

    /// <summary>Seismic moment of a set of patches, N m.</summary>
    public static double Moment(IEnumerable<SlipPatch> patches, double shearModulusPa = 30e9) =>
        patches.Sum(p => shearModulusPa * p.LengthKm * p.WidthKm * 1e6 * Math.Abs(p.SlipM));

    private const double EarthRadiusM = 6371000.0;

    /// <summary>Latitude and longitude (degrees) to east and north metres in the tangent plane at (lat0, lon0).</summary>
    public static (double East, double North) ToEastNorth(double lat, double lon, double lat0, double lon0) =>
        ((lon - lon0) * Math.PI / 180 * EarthRadiusM * Math.Cos(lat0 * Math.PI / 180), (lat - lat0) * Math.PI / 180 * EarthRadiusM);

    /// <summary>East and north metres in the tangent plane at (lat0, lon0) to latitude and longitude, degrees.</summary>
    public static (double Lat, double Lon) ToLatLon(double east, double north, double lat0, double lon0) =>
        (lat0 + north / EarthRadiusM * 180 / Math.PI, lon0 + east / (EarthRadiusM * Math.Cos(lat0 * Math.PI / 180)) * 180 / Math.PI);

    public enum SampleMask
    {
        None,
        NoActiveSources,
        SourceNearField,
        NonFinite
    }

    public readonly record struct Tensor(
        double EePa, double NnPa, double UuPa,
        double EnPa, double EuPa, double NuPa)
    {
        public bool IsFinite => double.IsFinite(EePa) && double.IsFinite(NnPa) &&
            double.IsFinite(UuPa) && double.IsFinite(EnPa) &&
            double.IsFinite(EuPa) && double.IsFinite(NuPa);

        public static Tensor operator *(Tensor a, double scale) => new(
            a.EePa * scale, a.NnPa * scale, a.UuPa * scale,
            a.EnPa * scale, a.EuPa * scale, a.NuPa * scale);

        public double[,] Matrix() => new[,]
        {
            { EePa, EnPa, EuPa }, { EnPa, NnPa, NuPa }, { EuPa, NuPa, UuPa }
        };
    }

    public readonly record struct Receiver(
        string Id, double Lat, double Lon, double DepthKm,
        double StrikeDeg, double DipDeg, double RakeDeg);

    public readonly record struct Result(
        string Id, double Lat, double Lon, double DepthKm,
        double DeltaCfsPa, double ShearPa, double UnclampingPa,
        double EffectiveFriction, Tensor Stress, SampleMask Mask,
        int SourcePatches, string Model, string StressUnits,
        string ConfidenceLabel, string DataProvenance)
    {
        public bool IsValid => Mask == SampleMask.None && double.IsFinite(DeltaCfsPa);
    }

    public readonly record struct PlaneSamplingRequest(
        string Id, double OriginLat, double OriginLon, double TopDepthKm,
        double StrikeDeg, double DipDeg, double RakeDeg,
        double LengthKm, double WidthKm, int AlongStrikeSamples,
        int DownDipSamples);

    public readonly record struct VolumeSamplingRequest(
        string IdPrefix, double OriginLat, double OriginLon,
        double MinEastKm, double MaxEastKm, int EastSamples,
        double MinNorthKm, double MaxNorthKm, int NorthSamples,
        double MinDepthKm, double MaxDepthKm, int DepthSamples,
        double StrikeDeg, double DipDeg, double RakeDeg);

    private readonly record struct Source(
        double EastM, double NorthM, double UpM,
        double StrikeDeg, double DipDeg, double RakeDeg,
        double LengthM, double WidthM, double SlipM);

    public static List<Result> Compute(
        IReadOnlyList<SlipPatch> patches, IReadOnlyList<Receiver> receivers,
        double effectiveFriction,
        double shearModulusPa = 30e9, double poisson = 0.25,
        double targetMomentNm = 0, int quadratureOrder = 4,
        double additionalMaskDistanceM = 0)
    {
        Validate(patches, receivers, effectiveFriction, shearModulusPa, poisson,
            quadratureOrder, additionalMaskDistanceM);
        var reference = Reference(patches);
        var sources = BuildSources(patches, reference);
        double scale = MomentScale(sources, shearModulusPa, targetMomentNm);
        var output = new Result[receivers.Count];
        // Receivers are independent: evaluated in parallel (a map or a volume has 10^4-10^6 of them).
        Parallel.For(0, receivers.Count, idx =>
        {
            var receiver = receivers[idx];
            var (tensor, mask) = TensorAt(receiver.Lat, receiver.Lon, receiver.DepthKm, reference, sources,
                shearModulusPa, poisson, scale, quadratureOrder, additionalMaskDistanceM);
            if (mask != SampleMask.None)
            {
                output[idx] = Masked(receiver, effectiveFriction, sources.Count, mask);
                return;
            }
            var (shear, unclamping) = Resolve(tensor, receiver.StrikeDeg, receiver.DipDeg, receiver.RakeDeg);
            output[idx] = new Result(receiver.Id, receiver.Lat, receiver.Lon,
                receiver.DepthKm, shear + effectiveFriction * unclamping,
                shear, unclamping, effectiveFriction, tensor, SampleMask.None,
                sources.Count, ModelName, Units, Confidence, Provenance);
        });
        return output.ToList();
    }

    /// <summary>
    /// The stress-change tensor (ENU, tensile positive, Pa) at points, without resolving it on a
    /// receiver: for receivers whose orientation depends on the total stress (optimally oriented
    /// planes) or that are resolved on several planes (both nodal planes of a mechanism).
    /// </summary>
    public static (Tensor Stress, SampleMask Mask)[] Tensors(
        IReadOnlyList<SlipPatch> patches, IReadOnlyList<(double Lat, double Lon, double DepthKm)> points,
        double shearModulusPa = 30e9, double poisson = 0.25,
        double targetMomentNm = 0, int quadratureOrder = 4, double additionalMaskDistanceM = 0)
    {
        Validate(patches, Array.Empty<Receiver>(), 0, shearModulusPa, poisson, quadratureOrder, additionalMaskDistanceM);
        ArgumentNullException.ThrowIfNull(points);
        var reference = Reference(patches);
        var sources = BuildSources(patches, reference);
        double scale = MomentScale(sources, shearModulusPa, targetMomentNm);
        var output = new (Tensor, SampleMask)[points.Count];
        Parallel.For(0, points.Count, i =>
        {
            var (lat, lon, depth) = points[i];
            output[i] = TensorAt(lat, lon, depth, reference, sources, shearModulusPa, poisson, scale, quadratureOrder, additionalMaskDistanceM);
        });
        return output;
    }

    /// <summary>Resolved shear (along the receiver's slip) and normal (unclamping positive) stress on a plane.</summary>
    public static (double Shear, double Unclamping) Resolve(Tensor tensor, double strikeDeg, double dipDeg, double rakeDeg)
    {
        var basis = Basis(strikeDeg, dipDeg, rakeDeg);
        double[] traction = MatVec(tensor.Matrix(), basis.Normal);
        return (Dot(traction, basis.Slip), Dot(traction, basis.Normal));
    }

    private static (Tensor, SampleMask) TensorAt(double lat, double lon, double depthKm, (double Lat, double Lon) reference,
        List<Source> sources, double mu, double nu, double scale, int order, double extraMask)
    {
        var nan = new Tensor(double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN);
        if (sources.Count == 0) return (nan, SampleMask.NoActiveSources);
        var en = ToEastNorth(lat, lon, reference.Lat, reference.Lon);
        var point = new[] { en.East, en.North, -depthKm * 1000.0 };
        var mask = MaskAt(point, sources, order, extraMask);
        if (mask != SampleMask.None) return (nan, mask);
        var tensor = StressAt(point, sources, mu, nu) * scale;
        return tensor.IsFinite ? (tensor, SampleMask.None) : (nan, SampleMask.NonFinite);
    }

    public static List<Result> SamplePlane(
        IReadOnlyList<SlipPatch> patches, PlaneSamplingRequest plane, double effectiveFriction,
        double shearModulusPa = 30e9, double poisson = 0.25,
        int quadratureOrder = 4, double additionalMaskDistanceM = 0)
    {
        if (plane.AlongStrikeSamples < 1 || plane.DownDipSamples < 1 ||
            !(plane.LengthKm > 0) || !(plane.WidthKm > 0))
            throw new ArgumentOutOfRangeException(nameof(plane));
        var receivers = new List<Receiver>(plane.AlongStrikeSamples * plane.DownDipSamples);
        double strike = Rad(plane.StrikeDeg), dip = Rad(plane.DipDeg);
        for (int row = 0; row < plane.DownDipSamples; row++)
        for (int col = 0; col < plane.AlongStrikeSamples; col++)
        {
            double w = (row + 0.5) / plane.DownDipSamples * plane.WidthKm * 1000.0;
            double s = ((col + 0.5) / plane.AlongStrikeSamples - 0.5) *
                       plane.LengthKm * 1000.0;
            double east = s * Math.Sin(strike) + w * Math.Cos(dip) * Math.Cos(strike);
            double north = s * Math.Cos(strike) - w * Math.Cos(dip) * Math.Sin(strike);
            var ll = ToLatLon(east, north, plane.OriginLat, plane.OriginLon);
            receivers.Add(new Receiver($"{plane.Id}:{row}:{col}", ll.Lat, ll.Lon,
                plane.TopDepthKm + w * Math.Sin(dip) / 1000.0,
                plane.StrikeDeg, plane.DipDeg, plane.RakeDeg));
        }
        return Compute(patches, receivers, effectiveFriction,
            shearModulusPa, poisson, quadratureOrder: quadratureOrder,
            additionalMaskDistanceM: additionalMaskDistanceM);
    }

    public static List<Result> SampleVolume(
        IReadOnlyList<SlipPatch> patches, VolumeSamplingRequest volume, double effectiveFriction,
        double shearModulusPa = 30e9, double poisson = 0.25,
        int quadratureOrder = 4, double additionalMaskDistanceM = 0)
    {
        if (volume.EastSamples < 1 || volume.NorthSamples < 1 ||
            volume.DepthSamples < 1 || volume.MaxEastKm < volume.MinEastKm ||
            volume.MaxNorthKm < volume.MinNorthKm ||
            volume.MaxDepthKm < volume.MinDepthKm)
            throw new ArgumentOutOfRangeException(nameof(volume));
        static double At(double min, double max, int count, int index) =>
            count == 1 ? 0.5 * (min + max) : min + (max - min) * index / (count - 1.0);
        var receivers = new List<Receiver>(volume.EastSamples * volume.NorthSamples *
                                           volume.DepthSamples);
        for (int k = 0; k < volume.DepthSamples; k++)
        for (int j = 0; j < volume.NorthSamples; j++)
        for (int i = 0; i < volume.EastSamples; i++)
        {
            double east = At(volume.MinEastKm, volume.MaxEastKm,
                volume.EastSamples, i) * 1000.0;
            double north = At(volume.MinNorthKm, volume.MaxNorthKm,
                volume.NorthSamples, j) * 1000.0;
            var ll = ToLatLon(east, north, volume.OriginLat, volume.OriginLon);
            receivers.Add(new Receiver($"{volume.IdPrefix}:{k}:{j}:{i}", ll.Lat, ll.Lon,
                At(volume.MinDepthKm, volume.MaxDepthKm, volume.DepthSamples, k),
                volume.StrikeDeg, volume.DipDeg, volume.RakeDeg));
        }
        return Compute(patches, receivers, effectiveFriction,
            shearModulusPa, poisson, quadratureOrder: quadratureOrder,
            additionalMaskDistanceM: additionalMaskDistanceM);
    }

    private static (double Lat, double Lon) Reference(IReadOnlyList<SlipPatch> patches) =>
        (patches.Average(p => p.Lat), patches.Average(p => p.Lon));

    private static List<Source> BuildSources(IReadOnlyList<SlipPatch> patches, (double Lat, double Lon) reference)
    {
        var sources = new List<Source>();
        foreach (var p in patches.Where(p => p.SlipM != 0 && double.IsFinite(p.SlipM) && p.LengthKm > 0 && p.WidthKm > 0))
        {
            var en = ToEastNorth(p.Lat, p.Lon, reference.Lat, reference.Lon);
            sources.Add(new Source(en.East, en.North, -p.DepthKm * 1000.0, p.StrikeDeg, p.DipDeg, p.RakeDeg,
                p.LengthKm * 1000.0, p.WidthKm * 1000.0, p.SlipM));
        }
        return sources;
    }

    private static double MomentScale(IReadOnlyList<Source> sources, double mu,
        double targetMomentNm)
    {
        if (!(targetMomentNm > 0)) return 1.0;
        double current = sources.Sum(s => mu * s.LengthM * s.WidthM * Math.Abs(s.SlipM));
        return current > 0 ? targetMomentNm / current : 1.0;
    }

    private static SampleMask MaskAt(double[] point, IReadOnlyList<Source> sources,
        int order, double extra)
    {
        foreach (var source in sources)
        {
            double singularTolerance = Math.Max(source.LengthM, source.WidthM) * 1e-6;
            if (DistanceToRectangle(point, source) <= Math.Max(extra, singularTolerance))
                return SampleMask.SourceNearField;
        }
        return SampleMask.None;
    }

    private static double DistanceToRectangle(double[] point, Source source)
    {
        var b = Basis(source.StrikeDeg, source.DipDeg, source.RakeDeg);
        var d = new[] { point[0] - source.EastM, point[1] - source.NorthM,
            point[2] - source.UpM };
        double da = Math.Max(0, Math.Abs(Dot(d, b.Along)) - source.LengthM / 2.0);
        double dw = Math.Max(0, Math.Abs(Dot(d, b.Down)) - source.WidthM / 2.0);
        double dn = Dot(d, b.Normal);
        return Math.Sqrt(da * da + dw * dw + dn * dn);
    }

    private static Tensor StressAt(double[] point, IReadOnlyList<Source> sources,
        double mu, double nu)
    {
        var gradient = new double[3, 3];
        double alpha = 1.0 / (2.0 * (1.0 - nu));
        foreach (var source in sources)
        {
            var b = Basis(source.StrikeDeg, source.DipDeg, source.RakeDeg);
            double de = point[0] - source.EastM;
            double dn = point[1] - source.NorthM;
            double localX = de * b.Along[0] + dn * b.Along[1];
            double localY = de * b.UpDip[0] + dn * b.UpDip[1];
            double rake = Rad(source.RakeDeg);
            var field = OkadaDc3d.Evaluate(alpha, localX, localY, point[2],
                -source.UpM, source.DipDeg,
                -source.LengthM / 2, source.LengthM / 2,
                -source.WidthM / 2, source.WidthM / 2,
                source.SlipM * Math.Cos(rake), source.SlipM * Math.Sin(rake), 0);
            if (!field.IsValid)
                return new Tensor(double.NaN, double.NaN, double.NaN,
                    double.NaN, double.NaN, double.NaN);

            double[,] local = field.Gradient();
            // Columns of R are local x=strike, y=up-dip-horizontal, z=up.
            double[,] rotation =
            {
                { b.Along[0], b.UpDip[0], 0 },
                { b.Along[1], b.UpDip[1], 0 },
                { 0, 0, 1 }
            };
            for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            for (int a = 0; a < 3; a++)
            for (int d = 0; d < 3; d++)
                gradient[i, j] += rotation[i, a] * local[a, d] * rotation[j, d];
        }
        double lambda = 2.0 * mu * nu / (1.0 - 2.0 * nu);
        double trace = gradient[0, 0] + gradient[1, 1] + gradient[2, 2];
        var sigma = new double[3, 3];
        for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            sigma[i, j] = (i == j ? lambda * trace : 0.0) +
                          mu * (gradient[i, j] + gradient[j, i]);
        return new Tensor(sigma[0, 0], sigma[1, 1], sigma[2, 2], sigma[0, 1],
            sigma[0, 2], sigma[1, 2]);
    }

    private static (double[] Along, double[] Down, double[] UpDip,
        double[] Normal, double[] Slip)
        Basis(double strikeDeg, double dipDeg, double rakeDeg)
    {
        double strike = Rad(strikeDeg), dip = Rad(dipDeg), rake = Rad(rakeDeg);
        var along = new[] { Math.Sin(strike), Math.Cos(strike), 0.0 };
        var down = new[] { Math.Cos(strike) * Math.Cos(dip),
            -Math.Sin(strike) * Math.Cos(dip), -Math.Sin(dip) };
        var upDip = new[] { -Math.Cos(strike), Math.Sin(strike), 0.0 };
        // Normal from the footwall into the hanging wall (upward for a dipping fault): the traction
        // σ·n is then what the hanging wall's side exerts, and s·σ·n is the shear that drives the
        // hanging wall along s. With the opposite normal the resolved shear,
        // and so ΔCFS, came out with the wrong sign: a receiver on strike beyond a strike-slip
        // rupture's tip, the classic loaded lobe (King et al. 1994), read as a stress shadow.
        var normal = new[] { Math.Cos(strike) * Math.Sin(dip),
            -Math.Sin(strike) * Math.Sin(dip), Math.Cos(dip) };
        var slip = new[] { Math.Cos(rake) * along[0] - Math.Sin(rake) * down[0],
            Math.Cos(rake) * along[1] - Math.Sin(rake) * down[1],
            Math.Cos(rake) * along[2] - Math.Sin(rake) * down[2] };
        return (along, down, upDip, normal, slip);
    }

    private static Result Masked(Receiver r, double friction, int count, SampleMask mask) =>
        new(r.Id, r.Lat, r.Lon, r.DepthKm, double.NaN, double.NaN, double.NaN,
            friction, new Tensor(double.NaN, double.NaN, double.NaN, double.NaN,
                double.NaN, double.NaN), mask, count, ModelName, Units, Confidence, Provenance);

    private static void Validate(IReadOnlyList<SlipPatch> patches,
        IReadOnlyList<Receiver> receivers, double friction, double mu, double nu,
        int order, double maskDistance)
    {
        if (patches == null || patches.Count == 0) throw new ArgumentException("No slip patches.", nameof(patches));
        if (receivers == null) throw new ArgumentNullException(nameof(receivers));
        if (!(mu > 0) || !double.IsFinite(mu)) throw new ArgumentOutOfRangeException(nameof(mu));
        if (!(nu > -1 && nu < 0.5) || !double.IsFinite(nu))
            throw new ArgumentOutOfRangeException(nameof(nu));
        if (!(friction >= 0) || !double.IsFinite(friction))
            throw new ArgumentOutOfRangeException(nameof(friction));
        if (order < 1 || order > 4) throw new ArgumentOutOfRangeException(nameof(order));
        if (!(maskDistance >= 0) || !double.IsFinite(maskDistance))
            throw new ArgumentOutOfRangeException(nameof(maskDistance));
    }

    private static double Rad(double degrees) => degrees * Math.PI / 180.0;
    private static double Dot(double[] a, double[] b) =>
        a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    private static double[] MatVec(double[,] a, double[] x) => new[]
    {
        a[0, 0] * x[0] + a[0, 1] * x[1] + a[0, 2] * x[2],
        a[1, 0] * x[0] + a[1, 1] * x[1] + a[1, 2] * x[2],
        a[2, 0] * x[0] + a[2, 1] * x[1] + a[2, 2] * x[2]
    };
}
