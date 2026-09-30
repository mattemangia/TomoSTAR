// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using TomoStar.Core.Forward;
using TomoStar.Core.Geo;
using TomoStar.Core.Model;
using TomoStar.Core.Numerics;

namespace TomoStar.Core.Tomography;

public enum VelocityParameterization
{
    /// <summary>P slowness and S slowness as independent node fields.</summary>
    VpVs,

    /// <summary>P slowness and the Vp/Vs ratio (Thurber 1993): S times constrain the ratio directly.</summary>
    VpVpVs
}

/// <summary>Where the starting Vs comes from.</summary>
public enum StartingVpVs
{
    /// <summary>The Vs of the starting model.</summary>
    FromModel,

    /// <summary>Vp divided by the Vp/Vs of the Wadati diagram of the picks, the same at every depth.</summary>
    FromData,

    /// <summary>Vp divided by <see cref="TomographySettings.StartVpVsValue"/>.</summary>
    Constant
}

public sealed class TomographySettings
{
    public RayMethod RayMethod { get; set; } = RayMethod.FastMarching;
    public bool UseOpenCl { get; set; } = true;
    public int ForwardRefinement { get; set; } = 2;
    public int Iterations { get; set; } = 5;
    public bool InvertP { get; set; } = true;
    public bool InvertS { get; set; } = true;
    public VelocityParameterization Parameterization { get; set; } = VelocityParameterization.VpVs;
    public bool JointHypocentres { get; set; } = true;
    public bool StationCorrections { get; set; }

    /// <summary>Damping of the velocity update (dimensionless: acts on fractional slowness change).</summary>
    public double DampingVelocity { get; set; } = 20;

    /// <summary>Laplacian smoothing weight of the total deviation from the starting model.</summary>
    public double Smoothing { get; set; } = 20;

    /// <summary>
    /// Factor on the damping and smoothing of the second field (Vs, or Vp/Vs with
    /// <see cref="VelocityParameterization.VpVpVs"/>). S picks are fewer and noisier than P picks, and
    /// the noise of two independently regularised fields adds up in their ratio: a Vp/Vs model usually
    /// wants more regularisation than Vp. Choose it with the Vp/Vs leakage and resolution tests.
    /// </summary>
    public double SRegularisationFactor { get; set; } = 1;

    /// <summary>How roughness is penalised: the Laplacian (default) or one of the other <see cref="Tomography.SmoothingMethod"/>s. The adaptive grid and parameter lattices use the Laplacian.</summary>
    public SmoothingMethod SmoothingMethod { get; set; } = SmoothingMethod.Laplacian;

    /// <summary>For the reweighted methods: the roughness beyond which a contrast is kept sharp, in robust spreads of the model's roughness.</summary>
    public double EdgeScale { get; set; } = 1;

    /// <summary>
    /// Vertical weight of the Laplacian relative to horizontal (1 = isotropic). With the smoothing in
    /// km, 0.1 lets velocity change with depth over about a third of the distance it may change
    /// laterally; on the central Italy geometry it recovers checkerboards as well as the node-based
    /// smoothing did, where 0.5 smoothed 8 km layers away (docs/validation.md).
    /// </summary>
    public double VerticalSmoothingWeight { get; set; } = 0.1;

    /// <summary>Second differences in km (default) or between neighbouring nodes.</summary>
    public SmoothingScale SmoothingScale { get; set; } = SmoothingScale.Kilometres;

    /// <summary>
    /// Invert for a layered (1-D) model: one Vp (and Vs) per grid depth, with the hypocentres and
    /// the station corrections: the "minimum 1-D model" that fits the data best (Kissling et al.
    /// 1994, JGR 99(B10), 19635-19646), the starting model a 3-D inversion should use.
    /// </summary>
    public bool Layered { get; set; }

    /// <summary>
    /// Largest change of slowness (as a fraction) one iteration may make at any node: a longer
    /// Gauss-Newton step is shortened (the whole update scaled), as the linearisation holds only near
    /// the current model. Without it a 1-D inversion from a poor starting model oscillates.
    /// </summary>
    public double MaxSlownessStep { get; set; } = 0.1;

    /// <summary>Bounds of Vp/Vs kept at every node whatever the parameterisation.</summary>
    public double VpVsMin { get; set; } = 1.4;

    /// <summary>
    /// Line search: when an update makes the weighted misfit worse, it is undone and half of it is
    /// applied instead, up to this many times (0 = accept every step).
    /// </summary>
    public int MaxStepHalvings { get; set; } = 3;

    public double VpVsMax { get; set; } = 2.5;

    /// <summary>Damping of hypocentre updates (per km and per s).</summary>
    public double DampingHypocentre { get; set; } = 1;

    public double DampingStation { get; set; } = 5;
    public int LsqrIterations { get; set; } = 400;
    public double LsqrTolerance { get; set; } = 1e-6;

    /// <summary>Residuals beyond this many scaled MADs (after removing each event's median) are rejected.</summary>
    public double OutlierMads { get; set; } = 4;

    public double OutlierFloorSeconds { get; set; } = 0.3;
    public double MaxAbsResidualSeconds { get; set; } = 5;
    public double VpMin { get; set; } = 1.5;
    public double VpMax { get; set; } = 14.5;
    public double VsMin { get; set; } = 0.7;
    public double VsMax { get; set; } = 8.5;
    public int StoredRayPoints { get; set; } = 48;
    public int MaxStoredRays { get; set; } = 50000;

    /// <summary>
    /// Adaptive parameterisation: the unknowns are octree cells, fine where rays are dense and coarse
    /// where they are few, rebuilt from the coverage of every iteration (<see cref="AdaptiveMesh"/>).
    /// Ignored by a layered inversion.
    /// </summary>
    public AdaptiveGridSettings Adaptive { get; set; } = new();

    /// <summary>
    /// Double-difference tomography (Zhang &amp; Thurber 2003): differential times of event pairs inverted with the
    /// absolute times; when enabled its iteration sets replace <see cref="Iterations"/>. Ignored by a
    /// layered inversion.
    /// </summary>
    public DoubleDifferenceSettings DoubleDifference { get; set; } = new();

    /// <summary>
    /// Velocity unknowns on a rotated, shifted lattice instead of the grid nodes (as in Zaharia et al.
    /// 2025; the forward and the output stay on the grid). Ignored by a layered or an
    /// adaptive inversion.
    /// </summary>
    public LatticeSettings Lattice { get; set; } = new();

    /// <summary>Starting Vs: the model's, or Vp over the data's (Wadati) or a given Vp/Vs. Applied by the tools that build the starting model.</summary>
    public StartingVpVs StartVpVs { get; set; } = StartingVpVs.FromModel;

    public double StartVpVsValue { get; set; } = 1.80;

    public TomographySettings Clone()
    {
        var c = (TomographySettings)MemberwiseClone();
        c.Adaptive = Adaptive.Clone();
        c.DoubleDifference = DoubleDifference.Clone();
        c.Lattice = Lattice.Clone();
        return c;
    }
}

public sealed class IterationStats
{
    public int Iteration { get; set; }
    public int Used { get; set; }
    public int Rejected { get; set; }
    public double RmsP { get; set; } = double.NaN;
    public double RmsS { get; set; } = double.NaN;
    public double Rms { get; set; }
    public double WeightedRms { get; set; }
    public double VarianceReductionPercent { get; set; }
    public int LsqrIterations { get; set; }
    public string LsqrStop { get; set; } = "";
    public double LsqrConditionEstimate { get; set; }
    public List<double> LsqrResidualHistory { get; set; } = [];
    public double UpdateNorm { get; set; }

    /// <summary>Factor the update was shortened by (1 = full Gauss-Newton step).</summary>
    public double StepScale { get; set; } = 1;
    public double Roughness { get; set; }
    public int MatrixRows { get; set; }
    public int MatrixColumns { get; set; }
    public long MatrixNonZeros { get; set; }
    public double MeanHypocentreShiftKm { get; set; }

    /// <summary>Double difference: RMS of the differential residuals in use, s (NaN without differential data).</summary>
    public double DifferentialRms { get; set; } = double.NaN;

    /// <summary>Double difference: differential times in use (weight above the minimum).</summary>
    public int DifferentialUsed { get; set; }
}

/// <summary>A traced ray kept for display (decimated).</summary>
public sealed record StoredRay(int Arrival, int Event, int Station, Phase Phase, float[] Lon, float[] Lat, float[] Depth);

/// <summary>Row metadata of the stored G matrix: which arrival each data row belongs to.</summary>
public sealed record MatrixLayout(int VelocityPOffset, int VelocitySOffset, int HypoOffset, int StationPOffset, int StationSOffset,
    int NodeCount, int DataRows, int[] RowArrival)
{
    /// <summary>Double difference: first row of the differential times (after the regularisation) and their number.</summary>
    public int DifferentialFirstRow { get; init; } = -1;

    public int DifferentialRows { get; init; }

    /// <summary>The velocity column of every grid node when the system was built (adaptive cells change between iterations).</summary>
    [System.Text.Json.Serialization.JsonIgnore]
    public int[]? ColumnOfNode { get; init; }
}

public sealed class TomographyResult
{
    public required SphericalGrid Grid { get; init; }
    public required double[] Vp { get; init; }
    public required double[] Vs { get; init; }
    public required double[] StartVp { get; init; }
    public required double[] StartVs { get; init; }
    public required double[] DwsP { get; init; }
    public required double[] DwsS { get; init; }
    public required int[] HitsP { get; init; }
    public required int[] HitsS { get; init; }
    public required ObservationSet Data { get; init; }
    public required List<IterationStats> Iterations { get; init; }
    public required List<StoredRay> Rays { get; init; }
    public CsrMatrix? LastMatrix { get; init; }
    public MatrixLayout? LastLayout { get; init; }

    /// <summary>The adaptive cells of the last update (null on the regular grid).</summary>
    public AdaptiveMesh? Mesh { get; init; }
    public List<string> Log { get; } = [];
    public string ForwardSolver { get; init; } = "";
}

/// <summary>
/// Local-earthquake travel-time tomography, P and S jointly with hypocentres and station
/// corrections (the coupled problem of Aki &amp; Lee 1976, JGR 81(23), 4381-4399; Thurber 1983, JGR
/// 88(B10), 8226-8236; parameterised on nodes with trilinear interpolation as in SIMULPS).
///
/// Each nonlinear iteration: forward (straight chords, or curved rays back-traced through
/// fast-marching tables of every station); residuals with robust outlier rejection; one damped,
/// smoothed least-squares update solved by LSQR (Paige &amp; Saunders 1982). Velocity unknowns are
/// fractional slowness changes relative to the starting model, so damping and smoothing are
/// dimensionless and comparable between P and S. Smoothing acts on the total deviation from the
/// starting model, not on each update, so it does not weaken with the number of iterations.
///
/// The same engine inverts synthetic data for checkerboard and spike tests.
/// </summary>
public sealed class TravelTimeTomography(SphericalGrid grid, TomographySettings settings, string workFolder, Action<string>? log = null)
{
    private ForwardDomain? _domain;

    /// <summary>
    /// Velocities outside the inversion grid, used when stations or events outside it are part of
    /// the data (see <see cref="ForwardDomain"/>). Null: the layer mean of the starting model.
    /// </summary>
    public VelocityModel1D? Background { get; init; }

    /// <summary>The forward domain of the last run.</summary>
    public ForwardDomain? Domain => _domain;

    private ForwardDomain MakeDomain(ObservationSet data, double[] startVp, double[] startVs) =>
        ForwardDomain.For(Grid, Settings.ForwardRefinement, data.HasOutsidePoints, data.Points(), Background, startVp, startVs, log);

    /// <summary>
    /// Sets the forward domain from <paramref name="data"/> for engines driven from outside
    /// <see cref="Run"/> (the time-lapse inversion traces each epoch with the domain of all epochs).
    /// </summary>
    internal void PrepareDomain(ObservationSet data, double[] startVp, double[] startVs) =>
        _domain = MakeDomain(data, startVp, startVs);

    public SphericalGrid Grid { get; } = grid;
    public TomographySettings Settings { get; } = settings;
    public string WorkFolder { get; } = workFolder;

    /// <summary>The model, hypocentres and station terms before an update, to undo it.</summary>
    private sealed class State
    {
        private double[] _sP = [], _sS = [];
        private (double Lon, double Lat, double Depth, double T0)[] _events = [];
        private (double P, double S)[] _stations = [];

        public static State Take(ObservationSet data, double[] sP, double[] sS) => new()
        {
            _sP = (double[])sP.Clone(), _sS = (double[])sS.Clone(),
            _events = data.Events.Select(e => (e.Lon, e.Lat, e.DepthKm, e.T0)).ToArray(),
            _stations = data.Stations.Select(x => (x.CorrectionP, x.CorrectionS)).ToArray()
        };

        public void Restore(ObservationSet data, double[] sP, double[] sS)
        {
            _sP.CopyTo(sP, 0);
            _sS.CopyTo(sS, 0);
            for (var i = 0; i < _events.Length; i++)
                (data.Events[i].Lon, data.Events[i].Lat, data.Events[i].DepthKm, data.Events[i].T0) = _events[i];
            for (var i = 0; i < _stations.Length; i++)
                (data.Stations[i].CorrectionP, data.Stations[i].CorrectionS) = _stations[i];
        }
    }

    /// <summary>Velocity unknowns per phase: one per node, per adaptive cell, or per layer for a 1-D inversion.</summary>
    private int VelocityColumns => Settings.Layered ? Grid.Nz : _lattice?.Count ?? _mesh?.CellCount ?? Grid.Count;

    private ParameterLattice? _lattice;

    /// <summary>The parameter lattice in use (null: the unknowns are the nodes, cells or layers).</summary>
    public ParameterLattice? Lattice => _lattice;

    private void PrepareLattice()
    {
        if (_lattice != null || !Settings.Lattice.Enabled || Settings.Layered || Adaptive) return;
        _lattice = new ParameterLattice(Grid, Settings.Lattice);
        log?.Invoke($"Parameterisation: {_lattice.Summary()}.");
    }

    /// <summary>The node field of the velocity block of a solution (P·x on a lattice, x[column of node] otherwise).</summary>
    private double[] NodeUpdate(double[] x, int offset, MatrixLayout l)
    {
        if (_lattice != null) return _lattice.Prolong(x, offset);
        var colOf = l.ColumnOfNode ?? ColumnOfNode;
        var f = new double[Grid.Count];
        for (var i = 0; i < f.Length; i++) f[i] = x[offset + colOf[i]];
        return f;
    }

    private int[]? _columnOfNode;
    private AdaptiveMesh? _mesh;

    /// <summary>The adaptive cells in use (null on the regular grid).</summary>
    public AdaptiveMesh? Mesh => _mesh;

    private bool Adaptive => Settings.Adaptive.Enabled && !Settings.Layered;

    /// <summary>The velocity column of each node (itself, its adaptive cell, or its layer).</summary>
    private int[] ColumnOfNode => _mesh != null && !Settings.Layered ? _mesh.LeafOfNode : _columnOfNode ??=
        Enumerable.Range(0, Grid.Count).Select(i => Settings.Layered ? Grid.Decompose(i).K : i).ToArray();

    /// <summary>
    /// Builds (first call) or refines the adaptive cells from the coverage of the current rays: rays
    /// per node or ray length per node, P and S together, over the arrivals in use.
    /// </summary>
    private void UpdateMesh(ObservationSet data, List<RayRow> rows)
    {
        if (!Adaptive || (_mesh != null && !Settings.Adaptive.RefineEachIteration)) return;
        var coverage = new double[Grid.Count];
        var byLength = Settings.Adaptive.Coverage == AdaptiveCoverage.Dws;
        foreach (var r in rows)
        {
            if (data.Arrivals[r.Arrival].Rejected) continue;
            for (var q = 0; q < r.Nodes.Length; q++) coverage[r.Nodes[q]] += byLength ? r.Length[q] : 1;
        }
        var previous = _mesh;
        _mesh = AdaptiveMesh.Build(Grid, coverage, Settings.Adaptive, previous);
        log?.Invoke($"  Adaptive grid{(previous == null ? "" : " refined")}: {_mesh.Summary()}.");
    }

    /// <summary>Mean of a node field over each layer (the current state of a layered model).</summary>
    private double[] LayerMeans(double[] field)
    {
        var sum = new double[Grid.Nz];
        var per = Grid.Nx * Grid.Ny;
        for (var i = 0; i < field.Length; i++) sum[Grid.Decompose(i).K] += field[i];
        return sum.Select(x => x / per).ToArray();
    }

    public TomographyResult Run(
        ObservationSet data, double[] startVp, double[] startVs,
        IProgress<(double, string)>? progress = null, CancellationToken ct = default)
    {
        var s = Settings;
        var n = Grid.Count;
        if (s.SmoothingMethod != SmoothingMethod.Laplacian)
            log?.Invoke(s.Layered || Adaptive || s.Lattice.Enabled
                ? $"Regularisation: {s.SmoothingMethod} acts on the node grid; the {(s.Layered ? "layered model" : Adaptive ? "adaptive cells" : "parameter lattice")} uses the Laplacian."
                : string.Create(System.Globalization.CultureInfo.InvariantCulture, $"Regularisation: {s.SmoothingMethod}{(s.SmoothingMethod is SmoothingMethod.TotalVariation or SmoothingMethod.EdgePreserving ? $", edge scale {s.EdgeScale:0.##} robust spreads (weights from the model of the previous iteration)" : "")}."));
        var sRefP = startVp.Select(v => 1 / v).ToArray();
        var sRefS = startVs.Select(v => 1 / v).ToArray();
        var kRef = startVp.Zip(startVs, (p, q) => p / q).ToArray();
        var sP = (double[])sRefP.Clone();
        var sS = (double[])sRefS.Clone();
        var stats = new List<IterationStats>();
        CsrMatrix? lastG = null;
        MatrixLayout? lastLayout = null;
        State? snapshot = null;
        double previousWrms = double.NaN;
        double[]? lastStep = null;
        MatrixLayout? lastLayoutApplied = null;
        var halvings = 0;
        List<RayRow> rows = [];
        double initialVariance = double.NaN;

        EikonalOpenCl? gpu = null;
        var solverName = s.RayMethod == RayMethod.Straight ? "straight rays" : "fast marching (CPU)";
        if (s.RayMethod == RayMethod.FastMarching && s.UseOpenCl)
        {
            gpu = EikonalOpenCl.TryCreate(log);
            if (gpu != null) solverName = $"OpenCL fast sweeping on {gpu.Device} + CPU";
        }
        if (!s.StationCorrections && data.Stations.Any(x => x.CorrectionP != 0 || x.CorrectionS != 0))
        {
            // Corrections of an earlier run belong to its model: without station terms none is applied.
            foreach (var x in data.Stations) { x.CorrectionP = 0; x.CorrectionS = 0; }
            log?.Invoke("Station corrections of earlier runs are not used (station corrections are off).");
        }
        _domain = MakeDomain(data, startVp, startVs);
        PrepareLattice();
        using var lsqrDevice = s.UseOpenCl ? LsqrOpenCl.TryCreate(log) : null;
        log?.Invoke($"Tomography: {data.Events.Count} events, {data.Stations.Count} stations, {data.Arrivals.Count} arrivals; grid {Grid.Nx}×{Grid.Ny}×{Grid.Nz}; forward {solverName}.");
        foreach (var w in data.Warnings) log?.Invoke("Warning: " + w);
        // Double difference: the event pairs are chosen once, from the starting hypocentres.
        var dd = s.DoubleDifference is { Enabled: true } && !s.Layered ? s.DoubleDifference : null;
        _differential = dd != null ? DoubleDifferencePairs.Build(data, dd, log) : null;
        var iterations = dd?.TotalIterations ?? s.Iterations;

        try
        {
            for (var it = 0; it <= iterations; it++)
            {
                ct.ThrowIfCancellationRequested();
                var final = it == iterations;
                var frac0 = (double)it / (iterations + 1);
                progress?.Report((frac0, final ? "Final forward pass" : $"Iteration {it + 1}/{iterations}: forward"));

                rows = Forward(data, sP, sS, gpu, keepRays: final, progress, frac0, ct, (final ? 1.0 : 0.5) / (iterations + 1));
                var st = Residuals(data, rows, it, ref initialVariance);
                if (snapshot != null && st.WeightedRms > previousWrms * 1.001 && halvings < s.MaxStepHalvings)
                {
                    // The last step made things worse: go back and take half of it.
                    halvings++;
                    snapshot.Restore(data, sP, sS);
                    var shorter = lastStep!.Select(v => v * Math.Pow(0.5, halvings)).ToArray();
                    log?.Invoke($"  Weighted RMS {st.WeightedRms:0.000} s is worse than {previousWrms:0.000} s: the step is halved ({halvings}/{s.MaxStepHalvings}).");
                    var baseScale = stats[^1].StepScale / (halvings > 1 ? Math.Pow(0.5, halvings - 1) : 1);
                    ApplyUpdate(data, shorter, lastLayoutApplied!, sP, sS, sRefP, sRefS, kRef, stats[^1]);
                    stats[^1].StepScale = baseScale * Math.Pow(0.5, halvings);
                    it--;
                    continue;
                }
                halvings = 0;
                stats.Add(st);
                var set = dd?.SetOf(Math.Min(it, iterations - 1));
                var ddRows = set != null ? DifferentialResiduals(data, rows, dd!, set, st) : null;
                var ddText = ddRows != null ? $" Differential: RMS {st.DifferentialRms:0.000} s, {st.DifferentialUsed} in use." : "";
                log?.Invoke(final
                    ? $"Final: RMS {st.Rms:0.000} s (P {st.RmsP:0.000}, S {st.RmsS:0.000}), {st.Used} used, {st.Rejected} rejected, variance reduction {st.VarianceReductionPercent:0.0}%.{ddText}"
                    : $"Iteration {it + 1}: RMS {st.Rms:0.000} s (P {st.RmsP:0.000}, S {st.RmsS:0.000}), {st.Used} used, {st.Rejected} rejected.{ddText}");
                if (final) break;

                progress?.Report((frac0 + 0.5 / (iterations + 1), $"Iteration {it + 1}/{iterations}: LSQR"));
                UpdateMesh(data, rows);
                var (g, rhs, layout) = BuildSystem(data, rows, sP, sS, sRefP, sRefS, kRef, s.DampingVelocity, s.Smoothing, set?.AbsoluteWeight ?? 1, ddRows);
                var history = new List<double>();
                var sol = LsqrOpenCl.SolveOrCpu(lsqrDevice, g, rhs, 0, s.LsqrIterations, s.LsqrTolerance, s.LsqrTolerance, 1e8, ct, (_, r) => history.Add(r), out _, log);
                st.LsqrIterations = sol.Iterations;
                st.LsqrStop = sol.StopReason;
                st.LsqrConditionEstimate = sol.EstimatedConditionNumber;
                st.LsqrResidualHistory = history;
                st.MatrixRows = g.RowCount;
                st.MatrixColumns = g.ColumnCount;
                st.MatrixNonZeros = g.NonZeroCount;
                st.UpdateNorm = SimdVector.Norm(sol.Solution.AsSpan(0, layout.HypoOffset >= 0 ? layout.HypoOffset : g.ColumnCount));
                snapshot = State.Take(data, sP, sS);
                previousWrms = st.WeightedRms;
                (lastStep, st.StepScale) = Capped(sol.Solution, layout);
                lastLayoutApplied = layout;
                ApplyUpdate(data, lastStep, layout, sP, sS, sRefP, sRefS, kRef, st);
                st.Roughness = Regularization.Roughness(Grid, sP.Select((x, i) => x / sRefP[i] - 1).ToArray(), s.VerticalSmoothingWeight, s.SmoothingScale);
                log?.Invoke($"  LSQR {sol.Iterations} it. ({sol.StopReason}), |Δm| {st.UpdateNorm:0.000}{(st.StepScale < 1 ? $" (step shortened to {st.StepScale:P0})" : "")}, mean hypocentre shift {st.MeanHypocentreShiftKm:0.00} km.");
                lastG = g;
                lastLayout = layout;
            }
        }
        finally
        {
            gpu?.Dispose();
        }

        var (dwsP, hitsP) = Coverage(rows, data, Phase.P);
        var (dwsS, hitsS) = Coverage(rows, data, Phase.S);
        return new TomographyResult
        {
            Grid = Grid,
            Vp = sP.Select(x => 1 / x).ToArray(),
            Vs = sS.Select(x => 1 / x).ToArray(),
            StartVp = startVp, StartVs = startVs,
            DwsP = dwsP, DwsS = dwsS, HitsP = hitsP, HitsS = hitsS,
            Data = data, Iterations = stats,
            Rays = rows.Where(r => r.Ray != null).Select(r => r.Ray!).Take(s.MaxStoredRays).ToList(),
            LastMatrix = lastG, LastLayout = lastLayout, ForwardSolver = solverName, Mesh = _mesh
        };
    }

    internal List<RayRow> Forward(ObservationSet data, double[] sP, double[] sS, EikonalOpenCl? gpu, bool keepRays,
        IProgress<(double, string)>? progress, double frac0, CancellationToken ct, double span = 0) =>
        new RayKernels(Grid, _domain ?? ForwardDomain.Refined(Grid.Definition, Settings.ForwardRefinement), Settings.RayMethod, WorkFolder, log)
            .Compute(data, sP, sS, gpu, keepRays ? Settings.StoredRayPoints : 0, progress, frac0, ct, span);

    /// <summary>One point of a regularisation trade-off curve.</summary>
    public sealed record TradeOffPoint(double Damping, double Smoothing, double ResidualNorm, double ModelNorm, double Roughness,
        double DataVariance, double ModelVariance, double[] Vp, double[] Vs);

    /// <summary>
    /// The first linearised step of the inversion solved for every (damping, smoothing) pair, with
    /// the forward computed once: the material of L-curves (Hansen 1992, SIAM Rev. 34(4), 561-580)
    /// and of data-variance/model-variance trade-off curves (Eberhart-Phillips 1986, BSSA 76(4),
    /// 1025-1052). Each point carries the model it produced, so any of them can be opened.
    /// </summary>
    public List<TradeOffPoint> TradeOff(ObservationSet data, double[] startVp, double[] startVs,
        IReadOnlyList<(double Damping, double Smoothing)> pairs, IProgress<(double, string)>? progress = null, CancellationToken ct = default)
    {
        var s = Settings;
        var sRefP = startVp.Select(v => 1 / v).ToArray();
        var sRefS = startVs.Select(v => 1 / v).ToArray();
        var kRef = startVp.Zip(startVs, (p, q) => p / q).ToArray();
        _domain = MakeDomain(data, startVp, startVs);
        PrepareLattice();
        EikonalOpenCl? gpu = s.RayMethod == RayMethod.FastMarching && s.UseOpenCl ? EikonalOpenCl.TryCreate(log) : null;
        List<RayRow> rows;
        try
        {
            progress?.Report((0, "Trade-off: forward"));
            rows = Forward(data, sRefP, sRefS, gpu, false, progress, 0, ct, 1.0 / (pairs.Count + 1));
        }
        finally
        {
            gpu?.Dispose();
        }
        var nan = double.NaN;
        Residuals(data, rows, 0, ref nan);
        UpdateMesh(data, rows);
        var points = new List<TradeOffPoint>();
        var n = Grid.Count;
        for (var q = 0; q < pairs.Count; q++)
        {
            ct.ThrowIfCancellationRequested();
            var (d, sm) = pairs[q];
            progress?.Report(((q + 1.0) / (pairs.Count + 1), $"Trade-off: damping {d:0.###}, smoothing {sm:0.###}"));
            var (g, rhs, layout) = BuildSystem(data, rows, sRefP, sRefS, sRefP, sRefS, kRef, d, sm);
            var sol = Lsqr.Solve(g, rhs, 0, s.LsqrIterations, s.LsqrTolerance, s.LsqrTolerance, 1e8, ct);
            // Data misfit after the linear update: ‖W(Gx − b)‖ over the data rows only.
            var pred = new double[g.RowCount];
            g.MultiplyInto(sol.Solution, pred);
            double rr = 0, dv = 0;
            for (var i = 0; i < layout.DataRows; i++)
            {
                var r = pred[i] - rhs[i];
                rr += r * r;
                var unweighted = r * data.Arrivals[layout.RowArrival[i]].Sigma;
                dv += unweighted * unweighted;
            }
            var mP = layout.VelocityPOffset >= 0 ? NodeUpdate(sol.Solution, layout.VelocityPOffset, layout) : new double[n];
            var mS = layout.VelocitySOffset >= 0 ? NodeUpdate(sol.Solution, layout.VelocitySOffset, layout) : new double[n];
            var model = mP.Concat(mS).ToArray();
            var mean = model.Average();
            var vp = startVp.Select((_, i) => 1 / (sRefP[i] * (1 + mP[i]))).ToArray();
            var vs = s.Parameterization == VelocityParameterization.VpVpVs
                ? vp.Select((v, i) => v / (kRef[i] * (1 + mS[i]))).ToArray()
                : startVs.Select((_, i) => 1 / (sRefS[i] * (1 + mS[i]))).ToArray();
            points.Add(new TradeOffPoint(d, sm, Math.Sqrt(rr), SimdVector.Norm(model),
                Regularization.Roughness(Grid, mP, s.VerticalSmoothingWeight, s.SmoothingScale),
                layout.DataRows > 0 ? dv / layout.DataRows : 0, model.Average(x => (x - mean) * (x - mean)), vp, vs));
        }
        return points;
    }

    /// <summary>
    /// The regularisation a trade-off run recommends. For each smoothing value, the damping at the
    /// corner of its L-curve (misfit against model norm). With several smoothing values, those corner
    /// points form a second curve, misfit against roughness, ordered by smoothing, and its corner gives
    /// the smoothing; with two values, the one with the smaller product of misfit and roughness.
    /// </summary>
    public static (double Damping, double Smoothing) ChooseRegularisation(IReadOnlyList<TradeOffPoint> points) =>
        ChooseRegularisation(points.Select(p => new TradeOffSummary(p.Damping, p.Smoothing, p.ResidualNorm, p.ModelNorm, p.Roughness)).ToList());

    /// <summary>The norms of one trade-off point, whatever the inversion (velocity or attenuation).</summary>
    public sealed record TradeOffSummary(double Damping, double Smoothing, double ResidualNorm, double ModelNorm, double Roughness);

    /// <inheritdoc cref="ChooseRegularisation(IReadOnlyList{TradeOffPoint})"/>
    public static (double Damping, double Smoothing) ChooseRegularisation(IReadOnlyList<TradeOffSummary> points)
    {
        if (points.Count == 0) throw new ArgumentException("No trade-off points.");
        var corners = points.GroupBy(p => p.Smoothing).OrderBy(g => g.Key).Select(g =>
        {
            var list = g.OrderBy(p => p.Damping).ToList();
            return list[Math.Clamp(Corner(list.Select(p => (p.ResidualNorm, p.ModelNorm)).ToList()), 0, list.Count - 1)];
        }).ToList();
        if (corners.Count == 1) return (corners[0].Damping, corners[0].Smoothing);
        if (corners.Count == 2)
        {
            var best = corners.OrderBy(p => p.ResidualNorm * Math.Max(1e-30, p.Roughness)).First();
            return (best.Damping, best.Smoothing);
        }
        var k = Math.Clamp(Corner(corners.Select(p => (p.ResidualNorm, Math.Max(1e-30, p.Roughness))).ToList()), 0, corners.Count - 1);
        return (corners[k].Damping, corners[k].Smoothing);
    }

    /// <summary>
    /// Index of the L-curve corner: the point of maximum curvature of (log ‖r‖, log ‖x‖), the
    /// curvature of the parametric curve estimated by finite differences along the point order.
    /// </summary>
    public static int Corner(IReadOnlyList<(double Residual, double Model)> curve)
    {
        if (curve.Count < 3) return curve.Count - 1;
        var x = curve.Select(c => Math.Log(Math.Max(1e-300, c.Residual))).ToArray();
        var y = curve.Select(c => Math.Log(Math.Max(1e-300, c.Model))).ToArray();
        var best = 1;
        var bestK = double.NegativeInfinity;
        for (var i = 1; i < curve.Count - 1; i++)
        {
            double dx = (x[i + 1] - x[i - 1]) / 2, dy = (y[i + 1] - y[i - 1]) / 2;
            double ddx = x[i + 1] - 2 * x[i] + x[i - 1], ddy = y[i + 1] - 2 * y[i] + y[i - 1];
            var k = Math.Abs(dx * ddy - dy * ddx) / Math.Pow(dx * dx + dy * dy, 1.5);
            if (double.IsFinite(k) && k > bestK) { bestK = k; best = i; }
        }
        return best;
    }

    // ---- Residuals and outliers ----

    internal IterationStats Residuals(ObservationSet data, List<RayRow> rows, int iteration, ref double initialVariance)
    {
        var s = Settings;
        foreach (var a in data.Arrivals) { a.Residual = double.NaN; a.Rejected = true; a.RejectReason = "no ray: the forward pass found no path (station or event outside the traced field)"; }
        foreach (var r in rows)
        {
            var a = data.Arrivals[r.Arrival];
            var ev = data.Events[a.Event];
            a.Residual = a.Time - ev.T0 - r.TCalc - data.Stations[a.Station].Correction(a.Phase);
            a.Rejected = false;
            a.RejectReason = "";
            if (iteration == 0) a.InitialResidual = a.Residual;
        }
        // Robust rejection on residuals with each event's median removed (that part is origin time).
        var byEvent = data.Arrivals.Where(a => !double.IsNaN(a.Residual)).GroupBy(a => a.Event)
            .ToDictionary(g => g.Key, g => Median(g.Select(a => a.Residual)));
        var centred = data.Arrivals.Where(a => !double.IsNaN(a.Residual)).Select(a => a.Residual - byEvent[a.Event]).ToArray();
        var mad = 1.4826 * Median(centred.Select(Math.Abs));
        var limit = Math.Max(s.OutlierFloorSeconds, s.OutlierMads * mad);
        var used = 0;
        var rejected = 0;
        double sum2 = 0, wsum2 = 0, wsum = 0, sumP = 0, sumS = 0;
        int nP = 0, nS = 0;
        foreach (var a in data.Arrivals)
        {
            if (double.IsNaN(a.Residual)) { rejected++; continue; }
            var deviation = a.Residual - byEvent[a.Event];
            if (Math.Abs(a.Residual) > s.MaxAbsResidualSeconds || Math.Abs(deviation) > limit)
            {
                a.Rejected = true;
                a.RejectReason = Math.Abs(a.Residual) > s.MaxAbsResidualSeconds
                    ? string.Create(System.Globalization.CultureInfo.InvariantCulture, $"|residual| {Math.Abs(a.Residual):0.000} s above the {s.MaxAbsResidualSeconds:0.###} s limit")
                    : string.Create(System.Globalization.CultureInfo.InvariantCulture, $"{deviation:+0.000;-0.000} s from the event median, beyond {limit:0.000} s ({s.OutlierMads:0.#} × MAD {mad:0.000} s, floor {s.OutlierFloorSeconds:0.###} s)");
                rejected++;
                continue;
            }
            used++;
            sum2 += a.Residual * a.Residual;
            var w = 1 / (a.Sigma * a.Sigma);
            wsum2 += w * a.Residual * a.Residual;
            wsum += w;
            if (a.Phase == Phase.P) { sumP += a.Residual * a.Residual; nP++; }
            else { sumS += a.Residual * a.Residual; nS++; }
        }
        var variance = used > 0 ? sum2 / used : 0;
        if (double.IsNaN(initialVariance)) initialVariance = variance;
        return new IterationStats
        {
            Iteration = iteration,
            Used = used,
            Rejected = rejected,
            Rms = Math.Sqrt(variance),
            RmsP = nP > 0 ? Math.Sqrt(sumP / nP) : double.NaN,
            RmsS = nS > 0 ? Math.Sqrt(sumS / nS) : double.NaN,
            WeightedRms = wsum > 0 ? Math.Sqrt(wsum2 / wsum) : 0,
            VarianceReductionPercent = initialVariance > 0 ? 100 * (1 - variance / initialVariance) : 0
        };
    }

    private static double Median(IEnumerable<double> values)
    {
        var a = values.Where(v => !double.IsNaN(v)).OrderBy(v => v).ToArray();
        if (a.Length == 0) return 0;
        return a.Length % 2 == 1 ? a[a.Length / 2] : 0.5 * (a[a.Length / 2 - 1] + a[a.Length / 2]);
    }

    // ---- Double difference ----

    private List<DifferentialTime>? _differential;

    /// <summary>A differential time in use: its rays, residual and weight.</summary>
    private sealed record DifferentialRow(DifferentialTime Time, RayRow RowA, RayRow RowB, double Residual, double Weight);

    /// <summary>
    /// Residuals and weights of the differential times for one iteration: (tA − tB)obs − (tA − tB)calc,
    /// i.e. the difference of the two absolute residuals (origin times included, station terms cancel);
    /// a priori weight of the phase and the pick errors, then reweighting by event separation
    /// and by residual (tricube tapers). Times whose arrivals were rejected are left out.
    /// </summary>
    private List<DifferentialRow>? DifferentialResiduals(ObservationSet data, List<RayRow> rows, DoubleDifferenceSettings dd,
        DoubleDifferenceIterationSet set, IterationStats st)
    {
        if (_differential == null) return null;
        var rowOf = new RayRow?[data.Arrivals.Count];
        foreach (var r in rows) rowOf[r.Arrival] = r;
        var pos = data.Events.Select(e => GeoMath.ToCartesian(e.Lon, e.Lat, e.DepthKm)).ToArray();
        var all = new (DifferentialTime T, RayRow A, RayRow B, double Res, double W, double Sep)?[_differential.Count];
        Parallel.For(0, _differential.Count, k =>
        {
            var t = _differential[k];
            var a = data.Arrivals[t.ArrivalA];
            var b = data.Arrivals[t.ArrivalB];
            if (a.Rejected || b.Rejected || rowOf[t.ArrivalA] is not { } ra || rowOf[t.ArrivalB] is not { } rb) return;
            if (double.IsNaN(a.Residual) || double.IsNaN(b.Residual)) return;
            var prior = a.Phase == Phase.P ? set.DifferentialWeightP : set.DifferentialWeightS;
            if (prior <= 0) return;
            all[k] = (t, ra, rb, a.Residual - b.Residual, prior, (pos[t.EventA] - pos[t.EventB]).Length);
        });
        var candidates = all.Where(x => x.HasValue).Select(x => x!.Value).ToList();
        var cutoff = set.ResidualCutoff;
        if (cutoff >= 1 && candidates.Count > 0)
            cutoff *= DoubleDifferencePairs.RobustSpread(candidates.Select(c => c.Res).ToList()).Sigma;
        var list = new List<DifferentialRow>();
        double sum2 = 0;
        foreach (var c in candidates)
        {
            // Relative weight (a priori × separation × residual), compared with the minimum weight,
            // then divided by the combined pick error of the two readings.
            var w = c.W;
            if (set.MaxSeparationKm > 0) w *= DoubleDifferencePairs.Tricube(c.Sep, set.MaxSeparationKm);
            if (cutoff > 0) w *= DoubleDifferencePairs.Tricube(c.Res, cutoff);
            if (w < dd.MinWeight) continue;
            var ea = data.Arrivals[c.T.ArrivalA].Sigma;
            var eb = data.Arrivals[c.T.ArrivalB].Sigma;
            list.Add(new DifferentialRow(c.T, c.A, c.B, c.Res, w / Math.Sqrt(ea * ea + eb * eb)));
            sum2 += c.Res * c.Res;
        }
        st.DifferentialUsed = list.Count;
        st.DifferentialRms = list.Count > 0 ? Math.Sqrt(sum2 / list.Count) : double.NaN;
        return list;
    }

    // ---- Linear system ----

    private (CsrMatrix G, double[] Rhs, MatrixLayout Layout) BuildSystem(
        ObservationSet data, List<RayRow> rows, double[] sP, double[] sS, double[] sRefP, double[] sRefS, double[] kRef,
        double damping, double smoothing, double absoluteWeight = 1, List<DifferentialRow>? differential = null)
    {
        var s = Settings;
        var n = VelocityColumns;
        var colOf = ColumnOfNode;
        var hasP = s.InvertP && data.Arrivals.Any(a => a.Phase == Phase.P);
        var hasS = s.InvertS && data.Arrivals.Any(a => a.Phase == Phase.S);
        var col = 0;
        var offP = -1;
        var offS = -1;
        var vpvs = s.Parameterization == VelocityParameterization.VpVpVs;
        // Vp/Vs parameterisation needs the P block even for S rows (t_S depends on s_P).
        if (hasP || (hasS && vpvs)) { offP = col; col += n; }
        if (hasS) { offS = col; col += n; }
        var offH = -1;
        if (s.JointHypocentres) { offH = col; col += 4 * data.Events.Count; }
        int offCP = -1, offCS = -1;
        if (s.StationCorrections)
        {
            if (data.Arrivals.Any(a => a.Phase == Phase.P)) { offCP = col; col += data.Stations.Count; }
            if (data.Arrivals.Any(a => a.Phase == Phase.S)) { offCS = col; col += data.Stations.Count; }
        }

        var b = new CsrBuilder(col);
        var rowArrival = new List<int>();
        var cols = new List<int>();
        var vals = new List<double>();
        var kCur = sS.Select((x, i) => x / sP[i]).ToArray(); // κ = Vp/Vs = sS/sP
        // Velocity and hypocentre terms of one arrival's row, added with a sign: + for an absolute
        // time or the first event of a pair, − for the second event of a pair.
        var lattice = _lattice;
        void Terms(RayRow r, Observation a, double sign, Dictionary<int, double> acc)
        {
            void Add(int c, double v) => acc[c] = acc.GetValueOrDefault(c) + sign * v;
            // A node's entry: its own column, or spread over the lattice parameters it interpolates.
            void AddNode(int off, int node, double v)
            {
                if (lattice == null) { Add(off + colOf[node], v); return; }
                var (cs, ws) = lattice.Of(node);
                for (var q = 0; q < cs.Length; q++) Add(off + cs.Span[q], v * ws.Span[q]);
            }
            var ev = data.Events[a.Event];
            if (a.Phase == Phase.P)
            {
                if (offP >= 0 && s.InvertP)
                    for (var q = 0; q < r.Nodes.Length; q++) AddNode(offP, r.Nodes[q], r.Length[q] * sRefP[r.Nodes[q]]);
            }
            else if (offS >= 0)
            {
                if (vpvs)
                {
                    for (var q = 0; q < r.Nodes.Length; q++)
                    {
                        var nd = r.Nodes[q];
                        if (offP >= 0 && s.InvertP) AddNode(offP, nd, r.Length[q] * kCur[nd] * sRefP[nd]);
                        AddNode(offS, nd, r.Length[q] * sP[nd] * kRef[nd]);
                    }
                }
                else
                    for (var q = 0; q < r.Nodes.Length; q++) AddNode(offS, r.Nodes[q], r.Length[q] * sRefS[r.Nodes[q]]);
            }
            if (offH >= 0 && !ev.Fixed)
            {
                var h = offH + 4 * a.Event;
                Add(h, 1);
                if (!ev.Outside) { Add(h + 1, r.Ge); Add(h + 2, r.Gn); Add(h + 3, r.Gz); }
            }
        }
        // Entries below 0.1 ms per unit fractional slowness are left out: in a differential row the
        // two rays share most of their path and cancel there, leaving values far below any pick error.
        static (int[] Cols, double[] Vals) Packed(Dictionary<int, double> acc)
        {
            var keys = acc.Where(x => Math.Abs(x.Value) >= 1e-4).Select(x => x.Key).ToArray();
            Array.Sort(keys);
            return (keys, keys.Select(k => acc[k]).ToArray());
        }
        // Rows are formed in parallel and added in order, so the system is the same on any number of
        // threads; in blocks, so the formed rows waiting to be added never hold more than a block.
        const int block = 32768;
        void AddRows(int count, Action<int, Dictionary<int, double>> form, Action<int, int[], double[]> add)
        {
            var packed = new (int[] Cols, double[] Vals)[Math.Min(block, count)];
            for (var start = 0; start < count; start += block)
            {
                var len = Math.Min(block, count - start);
                Parallel.For(0, len, () => new Dictionary<int, double>(), (i, _, acc) =>
                {
                    acc.Clear();
                    form(start + i, acc);
                    packed[i] = Packed(acc);
                    return acc;
                }, _ => { });
                for (var i = 0; i < len; i++)
                    if (packed[i].Cols.Length > 0) add(start + i, packed[i].Cols, packed[i].Vals);
                Array.Clear(packed);
            }
        }
        var used = rows.Where(r => !data.Arrivals[r.Arrival].Rejected).ToArray();
        if (absoluteWeight > 0)
            AddRows(used.Length, (i, acc) =>
            {
                var r = used[i];
                var a = data.Arrivals[r.Arrival];
                Terms(r, a, 1, acc);
                var cOff = a.Phase == Phase.P ? offCP : offCS;
                if (cOff >= 0) acc[cOff + a.Station] = acc.GetValueOrDefault(cOff + a.Station) + 1;
            }, (i, cols, vals) =>
            {
                var a = data.Arrivals[used[i].Arrival];
                b.AddRow(cols, vals, a.Residual, absoluteWeight / a.Sigma);
                rowArrival.Add(used[i].Arrival);
            });
        var dataRows = b.RowCount;

        // Regularisation.
        // A layer column stands for every node of the layer: the same penalty as the 3-D model
        // constrained to layers is the node penalty times √(nodes per layer). Without this a layered
        // inversion is in effect undamped.
        var regScale = s.Layered ? Math.Sqrt(Grid.Nx * Grid.Ny) : 1;
        damping *= regScale;
        smoothing *= regScale;
        var mesh = s.Layered ? null : _mesh;
        void Smooth(int off, double[] deviation, double f)
        {
            var w = smoothing * f;
            if (s.Layered) Regularization.AddLaplacian1D(b, n, off, w, LayerMeans(deviation));
            else if (lattice != null) lattice.AddLaplacian(b, off, w, s.VerticalSmoothingWeight, s.SmoothingScale, lattice.Restrict(deviation));
            else if (mesh != null) mesh.AddLaplacian(b, off, w, s.VerticalSmoothingWeight, s.SmoothingScale, mesh.Restrict(deviation));
            else Regularization.AddSmoothing(b, Grid, off, w, s.VerticalSmoothingWeight, deviation, s.SmoothingMethod, s.EdgeScale, deviation, s.SmoothingScale);
        }
        void Damp(int off, double f)
        {
            if (lattice != null) lattice.AddDamping(b, off, damping * f);
            else if (mesh != null) mesh.AddDamping(b, off, damping * f);
            else Regularization.AddDamping(b, off, n, damping * f);
        }
        if (offP >= 0)
        {
            Damp(offP, 1);
            Smooth(offP, sP.Select((x, i) => x / sRefP[i] - 1).ToArray(), 1);
        }
        if (offS >= 0)
        {
            var f = s.SRegularisationFactor > 0 ? s.SRegularisationFactor : 1;
            Damp(offS, f);
            Smooth(offS, vpvs ? kCur.Select((x, i) => x / kRef[i] - 1).ToArray() : sS.Select((x, i) => x / sRefS[i] - 1).ToArray(), f);
        }
        if (offH >= 0) Regularization.AddDamping(b, offH, 4 * data.Events.Count, s.DampingHypocentre);
        foreach (var off in new[] { offCP, offCS })
        {
            if (off < 0) continue;
            Regularization.AddDamping(b, off, data.Stations.Count, s.DampingStation);
            // Σ corrections = 0: otherwise a constant shift trades off with every origin time.
            var zc = Enumerable.Range(off, data.Stations.Count).ToArray();
            b.AddRow(zc, Enumerable.Repeat(1.0, zc.Length).ToArray(), 0, 10);
        }
        // Differential times (after the regularisation, so the absolute rows keep their numbering):
        // the difference of the two arrivals' rows; station terms cancel.
        var ddFirst = b.RowCount;
        if (differential != null)
            AddRows(differential.Count, (i, acc) =>
            {
                var d = differential[i];
                Terms(d.RowA, data.Arrivals[d.Time.ArrivalA], 1, acc);
                Terms(d.RowB, data.Arrivals[d.Time.ArrivalB], -1, acc);
            }, (i, cols, vals) => b.AddRow(cols, vals, differential[i].Residual, differential[i].Weight));
        var differentialRows = b.RowCount - ddFirst;
        var (g, rhs) = b.Build();
        return (g, rhs, new MatrixLayout(offP, offS, offH, offCP, offCS, n, dataRows, rowArrival.ToArray())
        {
            ColumnOfNode = lattice == null ? colOf : null, DifferentialFirstRow = differential != null ? ddFirst : -1, DifferentialRows = differentialRows
        });
    }

    /// <summary>
    /// Step-length control: the whole update (velocities, hypocentres and station terms alike, so
    /// its direction is kept) is scaled down when a velocity column would move more than
    /// <see cref="TomographySettings.MaxSlownessStep"/>. Returns the step to apply and its scale.
    /// </summary>
    private (double[] Step, double Scale) Capped(double[] x, MatrixLayout l)
    {
        double largest = 0;
        foreach (var off in new[] { l.VelocityPOffset, l.VelocitySOffset })
        {
            if (off < 0) continue;
            // On a lattice the limit is on what the nodes receive: a parameter reaching the grid
            // with a small weight may take a large value without changing any node much.
            if (_lattice != null)
                foreach (var v in NodeUpdate(x, off, l)) largest = Math.Max(largest, Math.Abs(v));
            else
                for (var c = 0; c < l.NodeCount; c++) largest = Math.Max(largest, Math.Abs(x[off + c]));
        }
        if (Settings.MaxSlownessStep <= 0 || largest <= Settings.MaxSlownessStep) return (x, 1);
        var f = Settings.MaxSlownessStep / largest;
        return (x.Select(v => v * f).ToArray(), f);
    }

    private void ApplyUpdate(ObservationSet data, double[] x, MatrixLayout l, double[] sP, double[] sS,
        double[] sRefP, double[] sRefS, double[] kRef, IterationStats st)
    {
        var s = Settings;
        var n = Grid.Count;
        var vpvs = s.Parameterization == VelocityParameterization.VpVpVs;
        var kCur = sS.Select((v, i) => v / sP[i]).ToArray();
        if (l.VelocityPOffset >= 0 && s.InvertP)
        {
            var dp = NodeUpdate(x, l.VelocityPOffset, l);
            for (var i = 0; i < n; i++)
                sP[i] = Math.Clamp(sP[i] + dp[i] * sRefP[i], 1 / s.VpMax, 1 / s.VpMin);
        }
        if (l.VelocitySOffset >= 0)
        {
            var ds = NodeUpdate(x, l.VelocitySOffset, l);
            for (var i = 0; i < n; i++)
            {
                if (vpvs)
                {
                    var k = Math.Clamp(kCur[i] + ds[i] * kRef[i], 1.3, 3.0);
                    sS[i] = sP[i] * k;
                }
                else sS[i] = sS[i] + ds[i] * sRefS[i];
                sS[i] = Math.Clamp(sS[i], 1 / s.VsMax, 1 / s.VsMin);
            }
        }
        else if (vpvs)
        {
            for (var i = 0; i < n; i++) sS[i] = Math.Clamp(sP[i] * kCur[i], 1 / s.VsMax, 1 / s.VsMin);
        }
        // Vp/Vs within physical bounds (sS/sP = Vp/Vs).
        for (var i = 0; i < n; i++) sS[i] = Math.Clamp(sS[i], sP[i] * s.VpVsMin, sP[i] * s.VpVsMax);
        if (l.HypoOffset >= 0)
        {
            double shift = 0;
            var moved = 0;
            for (var e = 0; e < data.Events.Count; e++)
            {
                var ev = data.Events[e];
                if (ev.Fixed) continue;
                var h = l.HypoOffset + 4 * e;
                ev.T0 += x[h];
                if (ev.Outside) continue;
                var de = x[h + 1];
                var dn = x[h + 2];
                var dist = Math.Sqrt(de * de + dn * dn);
                if (dist > 0)
                {
                    var (lon, lat) = GeoMath.Destination(ev.Lon, ev.Lat, Math.Atan2(de, dn) * GeoMath.Rad2Deg, dist);
                    ev.Lon = Math.Clamp(GeoMath.UnwrapLon(lon, Grid.CentreLon), Grid.LonDeg[0], Grid.LonDeg[^1]);
                    ev.Lat = Math.Clamp(lat, Grid.LatDeg[0], Grid.LatDeg[^1]);
                }
                ev.DepthKm = Math.Clamp(ev.DepthKm + x[h + 3], Grid.DepthKm[0], Grid.DepthKm[^1]);
                shift += Math.Sqrt(dist * dist + x[h + 3] * x[h + 3]);
                moved++;
            }
            st.MeanHypocentreShiftKm = moved > 0 ? shift / moved : 0;
        }
        if (l.StationPOffset >= 0) for (var i = 0; i < data.Stations.Count; i++) data.Stations[i].CorrectionP += x[l.StationPOffset + i];
        if (l.StationSOffset >= 0) for (var i = 0; i < data.Stations.Count; i++) data.Stations[i].CorrectionS += x[l.StationSOffset + i];
    }

    /// <summary>Derivative weight sum (ray length per node, Toomey &amp; Foulger 1989) and ray hit count.</summary>
    private (double[] Dws, int[] Hits) Coverage(List<RayRow> rows, ObservationSet data, Phase phase)
    {
        var dws = new double[Grid.Count];
        var hits = new int[Grid.Count];
        foreach (var r in rows)
        {
            if (data.Arrivals[r.Arrival].Phase != phase || data.Arrivals[r.Arrival].Rejected) continue;
            for (var q = 0; q < r.Nodes.Length; q++)
            {
                dws[r.Nodes[q]] += r.Length[q];
                hits[r.Nodes[q]]++;
            }
        }
        return (dws, hits);
    }

    /// <summary>Samples a 1-D model on the grid nodes (starting model).</summary>
    /// <summary>
    /// The 1-D model of a layered inversion: one node per grid depth with the layer's velocities, and
    /// below the grid the model it started from. Layers no ray sampled keep their starting values.
    /// </summary>
    public static VelocityModel1D LayeredModel(TomographyResult r, VelocityModel1D below, string name)
    {
        var g = r.Grid;
        var per = g.Nx * g.Ny;
        var m = new VelocityModel1D { Name = name, Reference = "minimum 1-D model (Kissling et al. 1994) from the picks" };
        for (var k = 0; k < g.Nz; k++)
        {
            double vp = 0, vs = 0;
            for (var n = 0; n < per; n++) { vp += r.Vp[k * per + n]; vs += r.Vs[k * per + n]; }
            m.Nodes.Add(new VelocityNode(g.DepthKm[k], vp / per, vs / per));
        }
        var bottom = g.DepthKm[^1];
        m.Nodes.Add(new VelocityNode(bottom, below.Vp(bottom + 1e-6), below.Vs(bottom + 1e-6)));
        m.Nodes.AddRange(below.Nodes.Where(x => x.DepthKm > bottom));
        return m;
    }

    public static (double[] Vp, double[] Vs) StartingModel(SphericalGrid g, VelocityModel1D model)
    {
        var vp = new double[g.Count];
        var vs = new double[g.Count];
        for (var k = 0; k < g.Nz; k++)
        {
            var p = model.Vp(g.DepthKm[k]);
            var q = model.Vs(g.DepthKm[k]);
            for (var j = 0; j < g.Ny; j++)
            for (var i = 0; i < g.Nx; i++)
            {
                vp[g.Index(i, j, k)] = p;
                vs[g.Index(i, j, k)] = q;
            }
        }
        return (vp, vs);
    }
}
