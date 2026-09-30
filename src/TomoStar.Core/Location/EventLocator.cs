using TomoStar.Core.Forward;
using TomoStar.Core.Geo;
using TomoStar.Core.Model;

namespace TomoStar.Core.Location;

public sealed class LocatorSettings
{
    /// <summary>Search the whole grid first (true) or only around the current location (false).</summary>
    public bool GlobalSearch { get; set; } = true;

    /// <summary>Maximum number of candidate nodes in the coarse search.</summary>
    public int CoarseCandidates { get; set; } = 8000;

    public int MaxGeigerIterations { get; set; } = 20;

    /// <summary>Levenberg-Marquardt damping of the Geiger step, relative to the diagonal.</summary>
    public double Damping { get; set; } = 0.01;

    public int MinPhases { get; set; } = 4;

    /// <summary>Residuals beyond this many scaled MADs are rejected (and the event relocated once more).</summary>
    public double OutlierMads { get; set; } = 4;

    /// <summary>Never reject a residual smaller than this, s.</summary>
    public double OutlierFloor { get; set; } = 0.5;

    public bool FixDepth { get; set; }
    public bool UseStationCorrections { get; set; } = true;

    /// <summary>When automatic picks are used (see <see cref="PickSelection"/>).</summary>
    public AutomaticPicks AutomaticPicks { get; set; } = AutomaticPicks.IfFewAnalystPicks;
}

public sealed record LocationResult(Hypocentre Hypocentre, int Used, int Rejected, int Iterations, string Message);

/// <summary>
/// Hypocentre location in a 3-D (or 1-D) model through precomputed travel-time tables.
///
/// 1. Grid search over the forward-grid nodes: at each node the origin time that best fits the
///    picks is the weighted median of (observed − predicted), and the misfit is the weighted L1
///    norm of what remains. L1 with a median origin time is insensitive to the odd wrong pick
///    (the approach of NonLinLoc's EDT-like misfits; Lomax et al. 2000).
/// 2. Geiger's method (Geiger 1912; linearised least squares) from the best node: the travel-time
///    derivatives with respect to the source are the gradient of the station's table at the source,
///    solved with Levenberg-Marquardt damping (Lee &amp; Stewart 1981, Principles and Applications of
///    Microearthquake Networks).
/// 3. Errors from the scaled covariance (GᵀWG)⁻¹·max(1, χ²/ν): the 1σ horizontal error ellipse
///    (semi-axes and azimuth of the major axis, from the eigen-decomposition of the east-north
///    block), depth and origin-time errors; azimuthal gap from the stations used.
/// </summary>
public sealed class EventLocator(TravelTimeTableSet tables, IReadOnlyDictionary<string, StationRecord> stations, LocatorSettings settings)
{
    private sealed record Obs(PickRecord Pick, TravelTimeTable Table, double ObsSeconds, double Weight, double Correction, StationRecord Station);

    public LocationResult Locate(EventRecord ev, CancellationToken ct = default)
    {
        var start = new Hypocentre { OriginTime = ev.OriginTime, Lon = ev.Lon, Lat = ev.Lat, DepthKm = ev.DepthKm };
        var reference = start.OriginTime;
        var all = new List<Obs>();
        foreach (var p in PickSelection.Best(ev.Picks, settings.AutomaticPicks))
        {
            if (!stations.TryGetValue(p.StationId, out var st)) continue;
            var table = tables.Get(st.Id, p.Phase);
            if (table == null) continue;
            var corr = settings.UseStationCorrections ? (p.Phase == Phase.P ? st.CorrectionP : st.CorrectionS) : 0;
            all.Add(new Obs(p, table, (p.Time - reference).TotalSeconds, p.Weight, corr, st));
        }
        if (all.Count < settings.MinPhases)
            return new LocationResult(start, all.Count, 0, 0, $"only {all.Count} usable phases (minimum {settings.MinPhases})");

        var used = all;
        var rejected = 0;
        Hypocentre? best = null;
        var iterations = 0;
        for (var round = 0; round < 3; round++)
        {
            ct.ThrowIfCancellationRequested();
            var (lon, lat, dep, t0) = GridSearch(used, start);
            (best, iterations) = Geiger(used, lon, lat, dep, t0, reference);
            // Reject outliers once the solution is stable, then relocate without them.
            var res = used.Select(o => o.ObsSeconds - (best.OriginTime - reference).TotalSeconds - o.Table.Time(best.Lon, best.Lat, best.DepthKm) - o.Correction).ToArray();
            var med = Median(res);
            var mad = 1.4826 * Median(res.Select(r => Math.Abs(r - med)).ToArray());
            var limit = Math.Max(settings.OutlierFloor, settings.OutlierMads * mad);
            var keep = used.Where((o, i) => Math.Abs(res[i] - med) <= limit).ToList();
            if (keep.Count == used.Count || keep.Count < settings.MinPhases) break;
            rejected += used.Count - keep.Count;
            used = keep;
        }

        // Store the residual of every pick (rejected ones included, so the diagnostics can show them).
        foreach (var o in all)
        {
            var tt = o.Table.Time(best!.Lon, best.Lat, best.DepthKm);
            o.Pick.Residual = o.ObsSeconds - (best.OriginTime - reference).TotalSeconds - tt - o.Correction;
        }
        best!.PhaseCount = used.Count;
        best.AzimuthalGapDeg = Gap(used.Select(o => GeoMath.Azimuth(best.Lon, best.Lat, o.Station.Lon, o.Station.Lat)));

        return new LocationResult(best, used.Count, rejected, iterations, "ok");
    }

    private (double Lon, double Lat, double Dep, double T0) GridSearch(List<Obs> obs, Hypocentre start)
    {
        var g = obs[0].Table.Grid;
        var weights = obs.Select(o => o.Weight).ToArray();
        var r = new double[obs.Count];

        double Misfit(int idx, out double t0)
        {
            for (var q = 0; q < obs.Count; q++)
            {
                var t = obs[q].Table.Node(idx);
                r[q] = float.IsFinite(t) ? obs[q].ObsSeconds - t - obs[q].Correction : double.NaN;
            }
            t0 = WeightedMedian(r, weights);
            double m = 0;
            for (var q = 0; q < obs.Count; q++) m += double.IsNaN(r[q]) ? 10 * weights[q] : weights[q] * Math.Abs(r[q] - t0);
            return m;
        }

        int bestIdx;
        if (settings.GlobalSearch)
        {
            var stride = Math.Max(1, (int)Math.Ceiling(Math.Cbrt((double)g.Count / settings.CoarseCandidates)));
            bestIdx = Search(0, g.Nx - 1, 0, g.Ny - 1, 0, g.Nz - 1, stride);
            var (bi, bj, bk) = g.Decompose(bestIdx);
            bestIdx = Search(bi - 2 * stride, bi + 2 * stride, bj - 2 * stride, bj + 2 * stride, bk - 2 * stride, bk + 2 * stride, 1);
        }
        else
        {
            var (fx, fy, fz) = g.Fractional(start.Lon, start.Lat, start.DepthKm);
            int ci = (int)Math.Round(fx), cj = (int)Math.Round(fy), ck = (int)Math.Round(fz);
            bestIdx = Search(ci - 6, ci + 6, cj - 6, cj + 6, ck - 6, ck + 6, 1);
        }
        var (i, j, k) = g.Decompose(bestIdx);
        Misfit(bestIdx, out var bestT0);
        return (g.LonDeg[i], g.LatDeg[j], settings.FixDepth ? start.DepthKm : g.DepthKm[k], bestT0);

        int Search(int i0, int i1, int j0, int j1, int k0, int k1, int stride)
        {
            i0 = Math.Max(0, i0); j0 = Math.Max(0, j0); k0 = Math.Max(0, k0);
            i1 = Math.Min(g.Nx - 1, i1); j1 = Math.Min(g.Ny - 1, j1); k1 = Math.Min(g.Nz - 1, k1);
            var best = -1;
            var bestM = double.PositiveInfinity;
            if (settings.FixDepth)
            {
                var (_, _, fz) = g.Fractional(start.Lon, start.Lat, start.DepthKm);
                k0 = k1 = Math.Clamp((int)Math.Round(fz), 0, g.Nz - 1);
            }
            for (var kk = k0; kk <= k1; kk += stride)
            for (var jj = j0; jj <= j1; jj += stride)
            for (var ii = i0; ii <= i1; ii += stride)
            {
                var idx = g.Index(ii, jj, kk);
                var m = Misfit(idx, out _);
                if (m < bestM) { bestM = m; best = idx; }
            }
            return best < 0 ? g.Index(i0, j0, k0) : best;
        }
    }

    private (Hypocentre Hyp, int Iterations) Geiger(List<Obs> obs, double lon, double lat, double dep, double t0, DateTime reference)
    {
        var g = obs[0].Table.Grid;
        var nPar = settings.FixDepth ? 3 : 4;
        var ata = new double[4, 4];
        var atb = new double[4];

        // Weighted misfit at a trial source and, when asked, the normal equations there.
        double Misfit(double x, double y, double z, double t, bool normal)
        {
            if (normal) { Array.Clear(ata); Array.Clear(atb); }
            double sum = 0;
            var rows = 0;
            foreach (var o in obs)
            {
                var tt = o.Table.Time(x, y, z);
                if (!double.IsFinite(tt)) continue;
                var res = o.ObsSeconds - t - tt - o.Correction;
                var w = o.Weight * o.Weight;
                sum += w * res * res;
                rows++;
                if (!normal) continue;
                var (ge, gn, gz) = o.Table.Gradient(x, y, z);
                Span<double> row = [1, ge, gn, gz];
                for (var a = 0; a < nPar; a++)
                {
                    atb[a] += w * row[a] * res;
                    for (var b = 0; b < nPar; b++) ata[a, b] += w * row[a] * row[b];
                }
            }
            return rows >= nPar ? sum : double.PositiveInfinity;
        }

        // Levenberg-Marquardt: a step is taken only when it lowers the misfit; otherwise the damping
        // grows (shorter steps towards the gradient direction) and the step is tried again. A fixed
        // damping could let a shallow event, where the travel-time field bends strongly, jump to the
        // bottom of the grid and stay there.
        var lambda = Math.Max(1e-6, settings.Damping);
        var chi2 = Misfit(lon, lat, dep, t0, true);
        var it = 0;
        for (; it < settings.MaxGeigerIterations && double.IsFinite(chi2); it++)
        {
            var accepted = false;
            double[]? step = null;
            while (lambda < 1e8)
            {
                var lhs = (double[,])ata.Clone();
                for (var a = 0; a < nPar; a++) lhs[a, a] *= 1 + lambda;
                step = Solve(lhs, atb, nPar);
                if (step == null) break;
                var (nlon, nlat) = GeoMath.Destination(lon, lat, Math.Atan2(step[1], step[2]) * GeoMath.Rad2Deg, Math.Sqrt(step[1] * step[1] + step[2] * step[2]));
                var tLon = Math.Clamp(GeoMath.UnwrapLon(nlon, g.CentreLon), g.LonDeg[0], g.LonDeg[^1]);
                var tLat = Math.Clamp(nlat, g.LatDeg[0], g.LatDeg[^1]);
                var tDep = settings.FixDepth ? dep : Math.Clamp(dep + step[3], g.DepthKm[0], g.DepthKm[^1]);
                var tT0 = t0 + step[0];
                var trial = Misfit(tLon, tLat, tDep, tT0, false);
                if (trial <= chi2)
                {
                    (lon, lat, dep, t0) = (tLon, tLat, tDep, tT0);
                    chi2 = Misfit(lon, lat, dep, t0, true);
                    lambda = Math.Max(1e-6, lambda / 10);
                    accepted = true;
                    break;
                }
                lambda *= 10;
            }
            if (!accepted || step == null) break;
            var moved = Math.Sqrt(step[1] * step[1] + step[2] * step[2] + (nPar > 3 ? step[3] * step[3] : 0));
            if (moved < 0.005 && Math.Abs(step[0]) < 0.001) { it++; break; }
        }
        // Covariance below from the undamped normal matrix at the solution.
        Misfit(lon, lat, dep, t0, true);

        var hyp = new Hypocentre
        {
            OriginTime = reference.AddSeconds(t0), Lon = lon, Lat = lat, DepthKm = dep,
            Rms = Math.Sqrt(obs.Sum(o =>
            {
                var r = o.ObsSeconds - t0 - o.Table.Time(lon, lat, dep) - o.Correction;
                return r * r;
            }) / obs.Count)
        };
        // Covariance from the undamped normal matrix, scaled by the reduced χ² when the fit is worse
        // than the pick uncertainties say.
        var cov = Invert(ata, nPar);
        if (cov != null)
        {
            var dof = Math.Max(1, obs.Count - nPar);
            var scale = Math.Max(1, chi2 / dof);
            double cee = cov[1, 1] * scale, cnn = cov[2, 2] * scale, cen = cov[1, 2] * scale;
            var tr = 0.5 * (cee + cnn);
            var det = cee * cnn - cen * cen;
            var root = Math.Sqrt(Math.Max(0, tr * tr - det));
            hyp.ErrorHorizontalKm = Math.Sqrt(Math.Max(0, tr + root));
            hyp.ErrorMinorKm = Math.Sqrt(Math.Max(0, tr - root));
            // Major axis at θ = ½·atan2(2σen, σee − σnn) from east; as an azimuth from north in [0, 180).
            var theta = 0.5 * Math.Atan2(2 * cen, cee - cnn) * GeoMath.Rad2Deg;
            hyp.ErrorAzimuthDeg = ((90 - theta) % 180 + 180) % 180;
            hyp.ErrorDepthKm = nPar > 3 ? Math.Sqrt(Math.Max(0, cov[3, 3] * scale)) : 0;
            hyp.ErrorTimeS = Math.Sqrt(Math.Max(0, cov[0, 0] * scale));
        }
        return (hyp, it);
    }

    private static double[]? Solve(double[,] a, double[] b, int n)
    {
        var m = new double[n, n + 1];
        for (var i = 0; i < n; i++)
        {
            for (var j = 0; j < n; j++) m[i, j] = a[i, j];
            m[i, n] = b[i];
        }
        for (var c = 0; c < n; c++)
        {
            var p = c;
            for (var r = c + 1; r < n; r++) if (Math.Abs(m[r, c]) > Math.Abs(m[p, c])) p = r;
            if (Math.Abs(m[p, c]) < 1e-30) return null;
            if (p != c) for (var j = 0; j <= n; j++) (m[c, j], m[p, j]) = (m[p, j], m[c, j]);
            for (var r = 0; r < n; r++)
            {
                if (r == c) continue;
                var f = m[r, c] / m[c, c];
                for (var j = c; j <= n; j++) m[r, j] -= f * m[c, j];
            }
        }
        var x = new double[n];
        for (var i = 0; i < n; i++) x[i] = m[i, n] / m[i, i];
        return x;
    }

    private static double[,]? Invert(double[,] a, int n)
    {
        var inv = new double[n, n];
        for (var c = 0; c < n; c++)
        {
            var e = new double[n];
            e[c] = 1;
            var col = Solve(a, e, n);
            if (col == null) return null;
            for (var r = 0; r < n; r++) inv[r, c] = col[r];
        }
        return inv;
    }

    internal static double WeightedMedian(double[] values, double[] weights)
    {
        var pairs = values.Select((v, i) => (v, w: weights[i])).Where(p => !double.IsNaN(p.v)).OrderBy(p => p.v).ToArray();
        if (pairs.Length == 0) return 0;
        var half = pairs.Sum(p => p.w) / 2;
        double acc = 0;
        foreach (var (v, w) in pairs)
        {
            acc += w;
            if (acc >= half) return v;
        }
        return pairs[^1].v;
    }

    internal static double Median(double[] v)
    {
        if (v.Length == 0) return 0;
        var s = v.Where(x => !double.IsNaN(x)).OrderBy(x => x).ToArray();
        if (s.Length == 0) return 0;
        return s.Length % 2 == 1 ? s[s.Length / 2] : 0.5 * (s[s.Length / 2 - 1] + s[s.Length / 2]);
    }

    /// <summary>Largest azimuthal gap between consecutive station azimuths, degrees.</summary>
    public static double Gap(IEnumerable<double> azimuths)
    {
        var a = azimuths.Distinct().OrderBy(x => x).ToArray();
        if (a.Length < 2) return 360;
        var gap = 360 - a[^1] + a[0];
        for (var i = 1; i < a.Length; i++) gap = Math.Max(gap, a[i] - a[i - 1]);
        return gap;
    }
}
