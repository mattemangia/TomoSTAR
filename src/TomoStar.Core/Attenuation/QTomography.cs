// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using TomoStar.Core.Forward;
using TomoStar.Core.Geo;
using TomoStar.Core.Model;
using TomoStar.Core.Numerics;
using TomoStar.Core.Tomography;

namespace TomoStar.Core.Attenuation;

public sealed class QTomographySettings
{
    public RayMethod RayMethod { get; set; } = RayMethod.FastMarching;
    public bool UseOpenCl { get; set; } = true;
    public int ForwardRefinement { get; set; } = 2;

    /// <summary>P t* gives Qp (rays and slowness from Vp), S t* gives Qs (from Vs).</summary>
    public Phase Phase { get; set; } = Phase.P;

    /// <summary>Starting (and reference) quality factor, uniform; replaced by the data's when <see cref="EstimateQ0"/>.</summary>
    public double Q0 { get; set; } = 300;

    /// <summary>
    /// Take the reference Q from the data: the uniform Q that best fits all t* (weighted least
    /// squares). With station terms a uniform change of Q is nearly indistinguishable from the
    /// terms, so the damping, not the data, would otherwise set the background.
    /// </summary>
    public bool EstimateQ0 { get; set; } = true;

    /// <summary>
    /// Damping and smoothing of the fractional change of 1/Q, relative to the typical weighted
    /// sensitivity of a t* to a node (1: a unit of regularisation costs as much as one datum), so the
    /// same values mean the same whatever the background Q, the t* uncertainties and the grid.
    /// </summary>
    public double Damping { get; set; } = 2;
    public double Smoothing { get; set; } = 10;
    public double VerticalSmoothingWeight { get; set; } = 0.1;

    /// <summary>Second differences in km (default) or between neighbouring nodes.</summary>
    public SmoothingScale SmoothingScale { get; set; } = SmoothingScale.Kilometres;

    public bool StationTerms { get; set; } = true;
    /// <summary>Damping of the station terms, relative to their sensitivity (1 / σ) as <see cref="Damping"/> is to the nodes'.</summary>
    public double DampingStation { get; set; } = 3;
    public int LsqrIterations { get; set; } = 400;
    public double QMin { get; set; } = 10;
    public double QMax { get; set; } = 5000;
    public double OutlierMads { get; set; } = 4;
    public int StoredRayPoints { get; set; } = 48;

    /// <summary>Adaptive parameterisation (octree cells refined where the t* paths are dense), see <see cref="AdaptiveMesh"/>.</summary>
    public AdaptiveGridSettings Adaptive { get; set; } = new();

    /// <summary>How roughness is penalised (see <see cref="Tomography.SmoothingMethod"/>); the adaptive grid and lattices use the Laplacian.</summary>
    public SmoothingMethod SmoothingMethod { get; set; } = SmoothingMethod.Laplacian;

    /// <summary>For the reweighted methods: the roughness beyond which a contrast is kept sharp, in robust spreads.</summary>
    public double EdgeScale { get; set; } = 1;

    /// <summary>For the reweighted methods: solves after the first, each with the weights of the previous model.</summary>
    public int ReweightingPasses { get; set; } = 4;

    /// <summary>Unknowns on a rotated, shifted lattice instead of the nodes (see <see cref="ParameterLattice"/>); ignored with the adaptive grid.</summary>
    public LatticeSettings Lattice { get; set; } = new();

    public QTomographySettings Clone()
    {
        var c = (QTomographySettings)MemberwiseClone();
        c.Adaptive = Adaptive.Clone();
        c.Lattice = Lattice.Clone();
        return c;
    }
}

public sealed class QTomographyResult
{
    public required SphericalGrid Grid { get; init; }
    public required double[] Q { get; init; }
    public required double[] Dws { get; init; }
    public required ObservationSet Data { get; init; }
    public required double RmsBefore { get; init; }
    public required double RmsAfter { get; init; }
    public required LsqrResult Lsqr { get; init; }
    public required List<double> LsqrHistory { get; init; }
    public required List<StoredRay> Rays { get; init; }
    public required double[] StationTerms { get; init; }

    /// <summary>The reference Q used (estimated from the data or as set).</summary>
    public double ReferenceQ { get; init; }
    public CsrMatrix? Matrix { get; init; }

    /// <summary>The adaptive cells (null on the regular grid).</summary>
    public AdaptiveMesh? Mesh { get; init; }
}

/// <summary>
/// Attenuation tomography from t*. Along a ray, t* = ∫ ds / (v Q) = Σₙ Lₙ sₙ qₙ with q = 1/Q on the
/// nodes, Lₙ the ray's trilinear kernel and sₙ the slowness of the velocity model the rays were
/// traced in, a linear problem in q once the velocity model is fixed (e.g. Rietbrock 2001, J. Geophys.
/// Res. 106(B3), 4141-4154; Eberhart-Phillips &amp; Chadwick 2002). Unknowns are fractional changes of q
/// relative to 1/Q₀; damping and smoothing are scaled by the typical sensitivity of the data to them
/// (see <see cref="QTomographySettings.Damping"/>), and optional station terms, damped the same way,
/// absorb near-surface attenuation under each site.
/// </summary>
public sealed class QTomography(SphericalGrid grid, QTomographySettings settings, string workFolder, Action<string>? log = null)
{
    /// <summary>
    /// The t* measurements of one phase (P for Qp, S for Qs) as observations. Events and stations outside the grid are
    /// used when <paramref name="outside"/> allows it and the path crosses the grid (the attenuation
    /// outside is then the reference Q, fixed). Hypocentres are those of the catalogue, which should
    /// be the relocated ones of the velocity inversion the rays are traced in.
    /// </summary>
    public static ObservationSet FromCatalogue(Catalogue catalogue, SphericalGrid grid, OutsideData outside = OutsideData.Exclude, Phase phase = Phase.P)
    {
        var set = new ObservationSet();
        var stations = catalogue.StationIndex();
        var evIndex = new Dictionary<string, int>();
        var stIndex = new Dictionary<string, int>(StringComparer.OrdinalIgnoreCase);
        var events = catalogue.EventIndex();
        var top = grid.DepthKm[0];
        int left = 0, unknown = 0;
        foreach (var m in catalogue.TStar.Where(m => !m.Disabled && m.Phase == phase))
        {
            if (!events.TryGetValue(m.EventId, out var ev) || !stations.TryGetValue(m.StationId, out var st)) { unknown++; continue; }
            var depth = Math.Max(top, st.DepthKm);
            var evInside = grid.Contains(ev.Lon, ev.Lat, ev.DepthKm);
            var stInside = grid.Contains(st.Lon, st.Lat, depth);
            if (!evInside || !stInside)
            {
                if (outside == OutsideData.Exclude || !ObservationSet.PathCrosses(grid, ev.Lon, ev.Lat, ev.DepthKm, st.Lon, st.Lat, depth))
                {
                    left++;
                    continue;
                }
            }
            if (!evIndex.TryGetValue(ev.Id, out var ei))
            {
                ei = set.Events.Count;
                evIndex[ev.Id] = ei;
                set.Events.Add(new EventState { Id = ev.Id, Lon = ev.Lon, Lat = ev.Lat, DepthKm = ev.DepthKm, Reference = ev.OriginTime, Fixed = true, Outside = !evInside });
            }
            if (!stIndex.TryGetValue(st.Id, out var si))
            {
                si = set.Stations.Count;
                stIndex[st.Id] = si;
                set.Stations.Add(new StationState { Id = st.Id, Lon = st.Lon, Lat = st.Lat, DepthKm = depth, Outside = !stInside });
            }
            set.Arrivals.Add(new Observation { Event = ei, Station = si, Phase = phase, Time = m.TStar, Sigma = Math.Max(0.002, m.Sigma) });
        }
        if (unknown > 0) set.Warnings.Add($"{unknown} t* measurements refer to an event or a station that is not in the catalogue and were left out.");
        if (left > 0) set.Warnings.Add($"{left} t* measurements with the event or the station outside the grid were left out" +
                                       (outside == OutsideData.Exclude ? "." : " (path not crossing the grid)."));
        return set;
    }

    /// <summary>Sensitivity scales of the last solve (1/Q block, station terms) the regularisation was multiplied by.</summary>
    private (double Model, double Station) _scales = (1, 1);

    /// <summary>1-D velocities outside the grid when outside data are used; null: layer mean of the model.</summary>
    public VelocityModel1D? Background { get; init; }

    /// <summary>What every solve shares: the rays, the slowness they were traced in, the reference q.</summary>
    private sealed class Prepared
    {
        public required List<RayRow> Rows { get; init; }
        /// <summary>Slowness of the phase inverted (P or S).</summary>
        public required double[] Slowness { get; init; }
        public required double Q0 { get; init; }
        public required Dictionary<int, double> StartResidual { get; init; }
        public required double RmsBefore { get; init; }
        public AdaptiveMesh? Mesh { get; init; }
        public ParameterLattice? Lattice { get; init; }
    }

    /// <summary>Traces the rays, sets the reference Q (from the data when asked) and the starting residuals.</summary>
    private Prepared Prepare(ObservationSet data, double[] vp, double[] vs, IProgress<(double, string)>? progress, CancellationToken ct)
    {
        var s = settings;
        var n = grid.Count;
        var sP = vp.Select(v => 1 / v).ToArray();
        var sS = vs.Select(v => 1 / v).ToArray();
        var q0 = 1 / Math.Max(1, s.Q0);
        EikonalOpenCl? gpu = s.RayMethod == RayMethod.FastMarching && s.UseOpenCl ? EikonalOpenCl.TryCreate(log) : null;
        List<RayRow> rows;
        try
        {
            progress?.Report((0.1, "Q tomography: tracing rays in the velocity model"));
            var domain = ForwardDomain.For(grid, s.ForwardRefinement, data.HasOutsidePoints, data.Points(), Background, vp, vs, log);
            rows = new RayKernels(grid, domain, s.RayMethod, workFolder, log)
                .Compute(data, sP, sS, gpu, s.StoredRayPoints, progress, 0.1, ct);
        }
        finally
        {
            gpu?.Dispose();
        }
        foreach (var a in data.Arrivals)
            if (a.Phase != s.Phase) throw new InvalidOperationException($"Q tomography for {s.Phase} was given {a.Phase} t*.");
        var slow = s.Phase == Phase.S ? sS : sP;
        if (s.EstimateQ0 && rows.Count > 0)
        {
            // t* = q · T* with T* = Σ L s the travel time: weighted least squares for a uniform q.
            double num = 0, den = 0;
            foreach (var r in rows)
            {
                var a = data.Arrivals[r.Arrival];
                var tt = r.OutsideTime;
                for (var k = 0; k < r.Nodes.Length; k++) tt += r.Length[k] * slow[r.Nodes[k]];
                var w = 1 / (a.Sigma * a.Sigma);
                num += w * a.Time * tt;
                den += w * tt * tt;
            }
            if (den > 0 && num > 0)
            {
                q0 = Math.Clamp(num / den, 1 / s.QMax, 1 / s.QMin);
                log?.Invoke($"Q tomography: reference Q from the data {1 / q0:0} (uniform model fitting all t*).");
            }
        }
        var qStart = Enumerable.Repeat(q0, n).ToArray();
        var startResidual = new Dictionary<int, double>();
        foreach (var a in data.Arrivals) { a.Rejected = true; a.Residual = double.NaN; a.RejectReason = "no ray: the forward pass found no path"; }
        foreach (var r in rows)
        {
            var a = data.Arrivals[r.Arrival];
            a.Residual = a.InitialResidual = a.Time - Predict(r, slow, qStart, q0);
            startResidual[r.Arrival] = a.Residual;
            a.Rejected = false;
            a.RejectReason = "";
        }
        AdaptiveMesh? mesh = null;
        if (s.Adaptive.Enabled)
        {
            var coverage = new double[n];
            var byLength = s.Adaptive.Coverage == AdaptiveCoverage.Dws;
            foreach (var r in rows)
                for (var k = 0; k < r.Nodes.Length; k++) coverage[r.Nodes[k]] += byLength ? r.Length[k] : 1;
            mesh = AdaptiveMesh.Build(grid, coverage, s.Adaptive);
            log?.Invoke($"Q tomography: adaptive grid {mesh.Summary()}.");
        }
        ParameterLattice? lattice = null;
        if (mesh == null && s.Lattice.Enabled)
        {
            lattice = new ParameterLattice(grid, s.Lattice);
            log?.Invoke($"Q tomography: {lattice.Summary()}.");
        }
        return new Prepared { Rows = rows, Slowness = slow, Q0 = q0, StartResidual = startResidual, RmsBefore = Rms(data), Mesh = mesh, Lattice = lattice };
    }

    /// <summary>Predicted t*: inside the grid Σ L s q, outside the fixed reference q₀ times the time spent there.</summary>
    private static double Predict(RayRow r, double[] slowness, double[] q, double q0)
    {
        var t = q0 * r.OutsideTime;
        for (var k = 0; k < r.Nodes.Length; k++) t += r.Length[k] * slowness[r.Nodes[k]] * q[r.Nodes[k]];
        return t;
    }

    /// <summary>
    /// One damped, smoothed least-squares solve on the arrivals not rejected; leaves the residuals of
    /// the fitted model on the arrivals. Returns the system, the solution, q on the nodes and the
    /// station terms.
    /// </summary>
    private (CsrMatrix G, LsqrResult Sol, List<double> History, double[] Q, double[] Terms) Solve(
        Prepared pr, ObservationSet data, double damping, double smoothing, CancellationToken ct, double[]? weightModel = null)
    {
        var s = settings;
        var mesh = pr.Mesh;
        var lattice = pr.Lattice;
        var n = mesh?.CellCount ?? lattice?.Count ?? grid.Count;
        var acc = lattice != null ? new Dictionary<int, double>() : null;
        var offSt = s.StationTerms ? n : -1;
        var b = new CsrBuilder(n + (s.StationTerms ? data.Stations.Count : 0));
        var cols = new List<int>();
        var vals = new List<double>();
        double sumModel = 0, sumStation = 0;
        long countModel = 0, countStation = 0;
        foreach (var r in pr.Rows)
        {
            var a = data.Arrivals[r.Arrival];
            if (a.Rejected) continue;
            var nodeVals = r.Nodes.Select((nd, k) => (double)r.Length[k] * pr.Slowness[nd] * pr.Q0).ToArray();
            if (mesh != null) mesh.MergeRow(r.Nodes, nodeVals, 0, cols, vals);
            else if (lattice != null)
            {
                acc!.Clear();
                for (var k = 0; k < r.Nodes.Length; k++)
                    lattice.Spread(r.Nodes[k], nodeVals[k], (c, v) => acc[c] = acc.GetValueOrDefault(c) + v);
                cols.Clear(); vals.Clear();
                foreach (var (c, v) in acc.OrderBy(x => x.Key)) { cols.Add(c); vals.Add(v); }
            }
            else { cols.Clear(); vals.Clear(); cols.AddRange(r.Nodes); vals.AddRange(nodeVals); }
            if (offSt >= 0) { cols.Add(offSt + a.Station); vals.Add(1); }
            for (var k = 0; k < cols.Count; k++)
            {
                var w = vals[k] / a.Sigma;
                if (cols[k] < n) { sumModel += w * w; countModel++; }
                else { sumStation += w * w; countStation++; }
            }
            b.AddRow(cols.ToArray(), vals.ToArray(), pr.StartResidual[r.Arrival], 1 / a.Sigma);
        }
        // Damping and smoothing are relative to the typical weighted sensitivity of a data row to a
        // parameter of their block: a node's entry is (time in the node) × q₀ / σ, a station term's 1 / σ,
        // so absolute weights would regularise a high-Q₀ or noisy data set far more than a low-Q₀ one,
        // and would leave the station terms almost free next to the nodes (they would take up all the
        // shallow attenuation). Scaled, a unit of regularisation costs as much as one row of data.
        var scaleModel = countModel > 0 ? Math.Sqrt(sumModel / countModel) : 1;
        var scaleStation = countStation > 0 ? Math.Sqrt(sumStation / countStation) : 1;
        _scales = (scaleModel, scaleStation);
        if (lattice != null)
        {
            lattice.AddDamping(b, 0, damping * scaleModel);
            lattice.AddLaplacian(b, 0, smoothing * scaleModel, s.VerticalSmoothingWeight, s.SmoothingScale, null);
        }
        else if (mesh != null)
        {
            mesh.AddDamping(b, 0, damping * scaleModel);
            mesh.AddLaplacian(b, 0, smoothing * scaleModel, s.VerticalSmoothingWeight, s.SmoothingScale, null);
        }
        else
        {
            Regularization.AddDamping(b, 0, n, damping * scaleModel);
            Regularization.AddSmoothing(b, grid, 0, smoothing * scaleModel, s.VerticalSmoothingWeight, null, s.SmoothingMethod, s.EdgeScale, weightModel, s.SmoothingScale);
        }
        if (offSt >= 0)
        {
            Regularization.AddDamping(b, offSt, data.Stations.Count, s.DampingStation * scaleStation);
            var zc = Enumerable.Range(offSt, data.Stations.Count).ToArray();
            b.AddRow(zc, Enumerable.Repeat(1.0, zc.Length).ToArray(), 0, 10 * scaleStation);
        }
        var (g, rhs) = b.Build();
        var history = new List<double>();
        LsqrResult sol;
        if (s.UseOpenCl && g.NonZeroCount >= LsqrOpenCl.MinNonZerosForDevice)
        {
            using var device = LsqrOpenCl.TryCreate(log);
            sol = LsqrOpenCl.SolveOrCpu(device, g, rhs, 0, s.LsqrIterations, 1e-6, 1e-6, 1e8, ct, (_, rn) => history.Add(rn), out _, log);
        }
        else sol = Lsqr.Solve(g, rhs, 0, s.LsqrIterations, 1e-6, 1e-6, 1e8, ct, (_, rn) => history.Add(rn));
        var q = new double[grid.Count];
        var onNodes = lattice?.Prolong(sol.Solution);
        for (var i = 0; i < q.Length; i++) q[i] = Math.Clamp(pr.Q0 * (1 + (onNodes?[i] ?? sol.Solution[mesh?.LeafOfNode[i] ?? i])), 1 / s.QMax, 1 / s.QMin);
        var terms = offSt >= 0 ? sol.Solution.Skip(offSt).Take(data.Stations.Count).ToArray() : [];
        foreach (var r in pr.Rows)
        {
            var a = data.Arrivals[r.Arrival];
            a.Residual = a.Time - Predict(r, pr.Slowness, q, pr.Q0) - (terms.Length > 0 ? terms[a.Station] : 0);
        }
        return (g, sol, history, q, terms);
    }

    /// <summary>
    /// Rejects the arrivals whose residual on the fitted model lies beyond OutlierMads scaled MADs
    /// (never tighter than twice the typical t* uncertainty). Judged on the fitted model, not on the
    /// uniform start: against the start the rays that cross the strongest anomalies have the largest
    /// residuals and would be thrown away although they are the most informative data.
    /// </summary>
    private int RejectOutliers(ObservationSet data)
    {
        var s = settings;
        var fitted = data.Arrivals.Where(a => !a.Rejected).Select(a => a.Residual).OrderBy(x => x).ToArray();
        var med = fitted.Length > 0 ? fitted[fitted.Length / 2] : 0;
        var mad = 1.4826 * (fitted.Length > 0 ? fitted.Select(x => Math.Abs(x - med)).OrderBy(x => x).ElementAt(fitted.Length / 2) : 0);
        var limit = Math.Max(s.OutlierMads * mad, 2 * data.Arrivals.Where(a => !a.Rejected).Select(a => a.Sigma).DefaultIfEmpty(0.002).Average());
        var rejected = 0;
        foreach (var a in data.Arrivals.Where(a => !a.Rejected))
            if (Math.Abs(a.Residual - med) > limit)
            {
                a.Rejected = true;
                rejected++;
                a.RejectReason = string.Create(System.Globalization.CultureInfo.InvariantCulture,
                    $"t* residual {a.Residual - med:+0.00000;-0.00000} s from the median of the fitted model, beyond {limit:0.00000} s ({s.OutlierMads:0.#} × MAD)");
            }
        return rejected;
    }

    public QTomographyResult Run(ObservationSet data, double[] vp, double[] vs, IProgress<(double, string)>? progress = null, CancellationToken ct = default)
    {
        var s = settings;
        var n = grid.Count;
        var pr = Prepare(data, vp, vs, progress, ct);
        progress?.Report((0.6, "Q tomography: LSQR"));
        var pass = Solve(pr, data, s.Damping, s.Smoothing, ct);
        var rejected = RejectOutliers(data);
        if (rejected > 0)
        {
            progress?.Report((0.8, "Q tomography: LSQR without outliers"));
            pass = Solve(pr, data, s.Damping, s.Smoothing, ct);
        }
        if (s.SmoothingMethod is SmoothingMethod.TotalVariation or SmoothingMethod.EdgePreserving && pr.Mesh == null && pr.Lattice == null)
            for (var k = 0; k < s.ReweightingPasses; k++)
            {
                // Reweighted least squares: the roughness weights of the fractional change of q of the last solve.
                var previous = pass.Q.Select(q => q / pr.Q0 - 1).ToArray();
                progress?.Report((0.85, $"Q tomography: reweighted solve {k + 1}/{s.ReweightingPasses}"));
                pass = Solve(pr, data, s.Damping, s.Smoothing, ct, previous);
            }
        if (s.SmoothingMethod != SmoothingMethod.Laplacian)
            log?.Invoke(pr.Mesh != null || pr.Lattice != null
                ? $"Q tomography: {s.SmoothingMethod} acts on the node grid; the {(pr.Mesh != null ? "adaptive cells" : "parameter lattice")} used the Laplacian."
                : $"Q tomography: regularisation {s.SmoothingMethod}.");
        var (g, sol, history, q, terms) = pass;
        var used = data.Arrivals.Count(a => !a.Rejected);
        var dws = new double[n];
        foreach (var r in pr.Rows.Where(r => !data.Arrivals[r.Arrival].Rejected))
            for (var k = 0; k < r.Nodes.Length; k++) dws[r.Nodes[k]] += r.Length[k];
        log?.Invoke($"Q tomography: {used} t* used, {rejected} rejected on the fitted model, RMS {pr.RmsBefore * 1000:0.0} → {Rms(data) * 1000:0.0} ms, LSQR {sol.Iterations} it. ({sol.StopReason}); " +
                     $"regularisation × {_scales.Model:0.###} (1/Q){(s.StationTerms ? $", × {_scales.Station:0.#} (station terms)" : "")}, the typical weighted sensitivity of a t* to each parameter.");
        return new QTomographyResult
        {
            Grid = grid, Q = q.Select(x => 1 / x).ToArray(), Dws = dws, Data = data, RmsBefore = pr.RmsBefore, RmsAfter = Rms(data),
            Lsqr = sol, LsqrHistory = history, Rays = pr.Rows.Where(r => r.Ray != null).Select(r => r.Ray!).ToList(),
            StationTerms = terms, Matrix = g, ReferenceQ = 1 / pr.Q0, Mesh = pr.Mesh
        };
    }

    /// <summary>One point of the Q trade-off: a complete solve for one damping and smoothing.</summary>
    public sealed record QTradeOffPoint(double Damping, double Smoothing, double ResidualNorm, double ModelNorm, double Roughness,
        double DataVariance, double ModelVariance, double[] Q);

    /// <summary>
    /// The Q problem solved for every (damping, smoothing) pair with the rays traced once: the
    /// material of L-curves (Hansen 1992) and of data-variance/model-variance trade-offs
    /// (Eberhart-Phillips 1986). The problem is linear in q, so each point is the complete inversion
    /// for its pair, not only a first step. Outliers are rejected once, on a fit with the current
    /// settings, and the same data enter every point, so the points differ only in regularisation.
    /// </summary>
    public List<QTradeOffPoint> TradeOff(ObservationSet data, double[] vp, double[] vs, IReadOnlyList<(double Damping, double Smoothing)> pairs,
        IProgress<(double, string)>? progress = null, CancellationToken ct = default)
    {
        var n = grid.Count;
        var pr = Prepare(data, vp, vs, progress, ct);
        Solve(pr, data, settings.Damping, settings.Smoothing, ct);
        var rejected = RejectOutliers(data);
        log?.Invoke($"Q trade-off: {data.Arrivals.Count(a => !a.Rejected)} t* used, {rejected} rejected on a fit with damping {settings.Damping} and smoothing {settings.Smoothing}; reference Q {1 / pr.Q0:0}.");
        var points = new List<QTradeOffPoint>();
        for (var i = 0; i < pairs.Count; i++)
        {
            ct.ThrowIfCancellationRequested();
            var (d, sm) = pairs[i];
            progress?.Report((0.2 + 0.8 * i / pairs.Count, $"Q trade-off: damping {d:0.###}, smoothing {sm:0.###}"));
            var (_, sol, _, q, _) = Solve(pr, data, d, sm, ct);
            var used = data.Arrivals.Where(a => !a.Rejected && !double.IsNaN(a.Residual)).ToList();
            var residualNorm = Math.Sqrt(used.Sum(a => Math.Pow(a.Residual / a.Sigma, 2)));
            var model = pr.Mesh != null ? pr.Mesh.Prolong(sol.Solution.AsSpan(0, pr.Mesh.CellCount)) : sol.Solution.AsSpan(0, n).ToArray();
            var mean = model.Average();
            points.Add(new QTradeOffPoint(d, sm, residualNorm, SimdVector.Norm(model),
                Regularization.Roughness(grid, model, settings.VerticalSmoothingWeight, settings.SmoothingScale),
                used.Count > 0 ? used.Average(a => a.Residual * a.Residual) : 0, model.Average(x => (x - mean) * (x - mean)),
                q.Select(x => 1 / x).ToArray()));
        }
        return points;
    }

    private static double Rms(ObservationSet d)
    {
        var used = d.Arrivals.Where(a => !a.Rejected && !double.IsNaN(a.Residual)).ToArray();
        return used.Length == 0 ? 0 : Math.Sqrt(used.Average(a => a.Residual * a.Residual));
    }
}
