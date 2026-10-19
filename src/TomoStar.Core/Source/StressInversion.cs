// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

namespace TomoStar.Core.Source;

/// <summary>
/// Clean-room reduced-stress inversion from double-couple focal mechanisms.
/// The linear system follows Michael (1984, doi:10.1029/JB089iB13p11517),
/// nodal-plane selection follows Vavrycuk (2014, doi:10.1093/gji/ggu224),
/// and the optional spatial coupling follows Hardebeck &amp; Michael (2006,
/// doi:10.1029/2005JB004144). No absolute stress magnitude is inferred.
/// </summary>
public sealed class FocalMechanismStressInversionService
{
    private const double DegToRad = Math.PI / 180.0;

    public StressInversionResult Invert(
        IReadOnlyList<FocalMechanismObservation> observations,
        StressInversionOptions? options = null)
        => Invert(observations, options, progress: null, CancellationToken.None);

    public StressInversionResult Invert(
        IReadOnlyList<FocalMechanismObservation> observations,
        StressInversionOptions? options,
        IProgress<StressInversionProgress>? progress,
        CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(observations);
        cancellationToken.ThrowIfCancellationRequested();
        options ??= new StressInversionOptions();
        var valid = observations.Where(IsValid).ToArray();
        if (valid.Length < Math.Max(5, options.MinimumMechanisms))
            return Failure(valid.Length, $"At least {Math.Max(5, options.MinimumMechanisms)} valid mechanisms are required for five independent tensor components.");
        if (valid.Select(o => o.EventId).Distinct(StringComparer.Ordinal).Count() != valid.Length)
            return Failure(valid.Length, "EventId values must be unique within one stress inversion.");

        var planes = valid.Select(ToPlanes).ToArray();
        CoreSolution? best = null;
        double[] frictions = (options.OptimizeFriction
            ? InclusiveRange(options.MinimumFriction, options.MaximumFriction, options.FrictionStep)
            : new[] { options.FrictionCoefficient }).ToArray();

        for (int frictionIndex = 0; frictionIndex < frictions.Length; frictionIndex++)
        {
            cancellationToken.ThrowIfCancellationRequested();
            double friction = frictions[frictionIndex];
            progress?.Report(new StressInversionProgress(StressInversionProgressStage.FrictionSearch,
                frictionIndex + 1, frictions.Length, $"Evaluating friction coefficient {friction:0.###}.", friction));
            var candidate = SolveIterative(valid, planes, friction, options, options.RandomSeed,
                progress, cancellationToken);
            if (candidate is null)
                continue;
            if (best is null || candidate.MeanInstability > best.MeanInstability + 1e-12 ||
                (Math.Abs(candidate.MeanInstability - best.MeanInstability) <= 1e-12 &&
                 candidate.MeanMisfitDeg < best.MeanMisfitDeg))
                best = candidate;
        }

        if (best is null)
            return Failure(valid.Length, "The weighted Michael system is rank deficient or exceeds the configured condition-number limit.");

        var bootstrap = Bootstrap(valid, options, best, progress, cancellationToken);
        var tensor = ToResult(best.Tensor, bootstrap);
        var messages = new List<string>();
        string status;
        if (!best.Converged)
        {
            status = "FAILED-QC";
            messages.Add("Nodal-plane selection did not converge within the iteration limit.");
        }
        else if (valid.Length < options.CertificationMinimumMechanisms)
        {
            status = "RESEARCH-GRADE";
            messages.Add($"Only {valid.Length} mechanisms are available; Vavrycuk (2014) demonstrates robust noisy-data recovery with sets of 20 or more, but this is not a universal resolution threshold.");
        }
        else if (options.BootstrapRealizations <= 0)
        {
            status = "VERIFIED-WITH-CAVEAT";
            messages.Add("Uncertainty resampling was disabled; enable bootstrap realizations for a VERIFIED global tensor.");
        }
        else if (bootstrap.Count < options.BootstrapRealizations)
        {
            status = "VERIFIED-WITH-CAVEAT";
            messages.Add($"Only {bootstrap.Count} of {options.BootstrapRealizations} uncertainty realizations passed the rank/conditioning gates.");
        }
        else
        {
            status = "VERIFIED";
            messages.Add("The normalized reduced-stress tensor passed convergence, conditioning, sample-size, provenance and bootstrap-completion gates.");
        }
        if (best.CycleResolved)
            messages.Add("Nodal-plane iteration entered a finite instability-selection cycle; the last internally consistent tensor/plane state was retained deterministically.");
        if (valid.Any(o => string.IsNullOrWhiteSpace(o.Provenance)))
            messages.Add("One or more mechanisms have no provenance; the numerical solution must not be used as a certified project result until provenance is supplied.");

        var result = new StressInversionResult
        {
            Tensor = tensor,
            MechanismCount = valid.Length,
            Iterations = best.Iterations,
            BootstrapCompleted = bootstrap.Count,
            FrictionCoefficient = best.Friction,
            MeanFaultInstability = best.MeanInstability,
            MeanSlipMisfitDeg = best.MeanMisfitDeg,
            ConditionNumber = best.ConditionNumber,
            SelectedFaultPlanes = SelectedResults(valid, best),
            QualityStatus = status,
            QualityMessages = messages
        };
        progress?.Report(new StressInversionProgress(StressInversionProgressStage.Completed, 1, 1,
            $"Stress inversion completed with status {result.QualityStatus}.", result.FrictionCoefficient));
        return result;
    }

    public StressGridResult InvertSpatial(StressGridRequest request)
        => InvertSpatial(request, progress: null, CancellationToken.None);

    public StressGridResult InvertSpatial(
        StressGridRequest request,
        IProgress<StressInversionProgress>? progress,
        CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(request);
        cancellationToken.ThrowIfCancellationRequested();
        ValidateGrid(request);
        double[] xs = AxisValues(request.MinimumLongitudeDeg, request.MaximumLongitudeDeg, request.LongitudeStepDeg);
        double[] ys = AxisValues(request.MinimumLatitudeDeg, request.MaximumLatitudeDeg, request.LatitudeStepDeg);
        double[] zs = AxisValues(request.MinimumDepthKm, request.MaximumDepthKm, request.DepthStepKm);
        int total = checked(xs.Length * ys.Length * zs.Length);
        var work = new GridWork?[total];

        for (int iz = 0; iz < zs.Length; iz++)
        for (int iy = 0; iy < ys.Length; iy++)
        for (int ix = 0; ix < xs.Length; ix++)
        {
            cancellationToken.ThrowIfCancellationRequested();
            int flat = Flat(ix, iy, iz, xs.Length, ys.Length);
            progress?.Report(new StressInversionProgress(StressInversionProgressStage.SpatialGridNode,
                flat + 1, total, $"Preparing spatial stress node {flat + 1}/{total}.",
                LongitudeDeg: xs[ix], LatitudeDeg: ys[iy], DepthKm: zs[iz]));
            var local = request.Observations
                .Select(o => (Observation: o, KernelWeight: SpatialKernelWeight(o, xs[ix], ys[iy], zs[iz], request)))
                .Where(x => x.KernelWeight > 0)
                .Select(x => WithWeight(x.Observation, x.Observation.Weight * x.KernelWeight))
                .ToArray();
            int minimumMechanisms = Math.Max(5, request.InversionOptions.MinimumMechanisms);
            double effectiveCount = EffectiveSampleCount(local);
            if (local.Length < minimumMechanisms || effectiveCount < minimumMechanisms)
                continue;
            var localOptions = CloneWithoutBootstrap(request.InversionOptions);
            var inversion = Invert(local, localOptions, progress: null, cancellationToken);
            if (inversion.Tensor is null || inversion.QualityStatus == "FAILED-QC")
                continue;
            var selected = inversion.SelectedFaultPlanes.ToDictionary(p => p.EventId, StringComparer.Ordinal);
            var normal = BuildNormalEquations(local, selected);
            if (normal is null)
                continue;
            work[flat] = new GridWork(local, selected, normal.Value.H, normal.Value.B,
                Components(inversion.Tensor), inversion.ConditionNumber, inversion.QualityStatus, effectiveCount);
        }

        var active = Enumerable.Range(0, total).Where(i => work[i] is not null).ToArray();
        var activeIndex = active.Select((flat, i) => (flat, i)).ToDictionary(x => x.flat, x => x.i);
        var edges = new List<(int A, int B)>();
        foreach (int flat in active)
        {
            Unflat(flat, xs.Length, ys.Length, out int ix, out int iy, out int iz);
            AddEdge(ix + 1, iy, iz);
            AddEdge(ix, iy + 1, iz);
            AddEdge(ix, iy, iz + 1);

            void AddEdge(int nx, int ny, int nz)
            {
                if (nx >= xs.Length || ny >= ys.Length || nz >= zs.Length) return;
                int other = Flat(nx, ny, nz, xs.Length, ys.Length);
                if (activeIndex.TryGetValue(other, out int b)) edges.Add((activeIndex[flat], b));
            }
        }

        int dimension = active.Length * 5;
        var x = new double[dimension];
        var rhs = new double[dimension];
        for (int a = 0; a < active.Length; a++)
        {
            var w = work[active[a]]!;
            Array.Copy(w.Initial, 0, x, a * 5, 5);
            Array.Copy(w.B, 0, rhs, a * 5, 5);
        }

        double damping = Math.Max(0, request.NeighborDamping);
        double[] Apply(double[] v)
        {
            var y = new double[v.Length];
            for (int a = 0; a < active.Length; a++)
            {
                var h = work[active[a]]!.H;
                for (int r = 0; r < 5; r++)
                for (int c = 0; c < 5; c++)
                    y[a * 5 + r] += h[r, c] * v[a * 5 + c];
            }
            foreach (var (a, b) in edges)
            for (int c = 0; c < 5; c++)
            {
                double d = damping * (v[a * 5 + c] - v[b * 5 + c]);
                y[a * 5 + c] += d;
                y[b * 5 + c] -= d;
            }
            return y;
        }

        var (solution, iterations, converged) = dimension == 0
            ? (Array.Empty<double>(), 0, false)
            : ConjugateGradient(Apply, rhs, x, request.MaximumConjugateGradientIterations, 1e-9,
                progress, cancellationToken);

        var nodes = new StressGridNodeResult[total];
        for (int flat = 0; flat < total; flat++)
        {
            Unflat(flat, xs.Length, ys.Length, out int ix, out int iy, out int iz);
            if (!activeIndex.TryGetValue(flat, out int ai))
            {
                nodes[flat] = new StressGridNodeResult
                {
                    LongitudeDeg = xs[ix], LatitudeDeg = ys[iy], DepthKm = zs[iz],
                    QualityStatus = "FAILED-QC"
                };
                continue;
            }
            var c = solution.AsSpan(ai * 5, 5).ToArray();
            var tensor = TensorFromComponents(c, null);
            var local = work[flat]!;
            double misfit = MeanMisfit(local.Observations, local.Selected, tensor);
            bool nodeVerified = converged && request.SpatialWeighting == StressSpatialWeighting.Gaussian &&
                                local.EffectiveMechanismCount >= Math.Max(5, request.InversionOptions.MinimumMechanisms) &&
                                local.Observations.All(o => !string.IsNullOrWhiteSpace(o.Provenance)) &&
                                double.IsFinite(misfit) && misfit <= 45;
            nodes[flat] = new StressGridNodeResult
            {
                LongitudeDeg = xs[ix], LatitudeDeg = ys[iy], DepthKm = zs[iz],
                MechanismCount = local.Observations.Length, Tensor = tensor,
                EffectiveMechanismCount = local.EffectiveMechanismCount,
                MeanSlipMisfitDeg = misfit,
                QualityStatus = nodeVerified ? "VERIFIED" : converged ? local.QualityStatus : "FAILED-QC"
            };
        }

        var messages = new List<string>();
        if (active.Length == 0) messages.Add("No grid node has enough well-conditioned focal mechanisms.");
        if (!converged && active.Length > 0) messages.Add("The simultaneous damped system did not meet the conjugate-gradient tolerance.");
        bool verified = converged && active.Length > 0 && request.SpatialWeighting == StressSpatialWeighting.Gaussian &&
                        active.All(flat => nodes[flat].QualityStatus == "VERIFIED");
        if (verified)
            messages.Add("The spatial normalized reduced-stress solution passed convergence, conditioning, effective-sample, provenance, slip-misfit and Gaussian-localization resolution gates.");
        else if (converged)
            messages.Add("The spatial solve converged, but at least one resolved node did not pass every certification gate.");
        messages.Add(request.SpatialWeighting == StressSpatialWeighting.Gaussian
            ? $"Data were localized with a Gaussian kernel truncated at {request.SearchRadiusStandardDeviations:0.###} standard deviations; search radii are cutoffs, not uniform averaging windows."
            : "Legacy uniform spatial weighting was requested; adjacent nodes may be nearly identical when search radii greatly exceed grid spacing.");
        var result = new StressGridResult
        {
            Nodes = nodes, LongitudeCount = xs.Length, LatitudeCount = ys.Length, DepthCount = zs.Length,
            Damping = damping, ConjugateGradientIterations = iterations, Converged = converged,
            QualityStatus = verified ? "VERIFIED" : converged ? "VERIFIED-WITH-CAVEAT" : "FAILED-QC", QualityMessages = messages
        };
        progress?.Report(new StressInversionProgress(StressInversionProgressStage.Completed, 1, 1,
            $"Spatial stress inversion completed with status {result.QualityStatus}."));
        return result;
    }

    private static CoreSolution? SolveIterative(
        FocalMechanismObservation[] observations, PlanePair[] planes, double friction,
        StressInversionOptions options, int seed,
        IProgress<StressInversionProgress>? progress, CancellationToken cancellationToken)
    {
        if (!double.IsFinite(friction) || friction < 0) return null;
        var random = new Random(seed);
        var selected = Enumerable.Range(0, observations.Length).Select(_ => random.Next(2)).ToArray();
        TensorState? tensor = null;
        double condition = double.PositiveInfinity;
        int iterations = 0;
        bool converged = false;
        bool cycleResolved = false;
        var visited = new HashSet<string>(StringComparer.Ordinal) { SelectionKey(selected) };
        for (int iteration = 1; iteration <= Math.Max(1, options.MaximumIterations); iteration++)
        {
            cancellationToken.ThrowIfCancellationRequested();
            iterations = iteration;
            progress?.Report(new StressInversionProgress(StressInversionProgressStage.NodalPlaneIteration,
                iteration, Math.Max(1, options.MaximumIterations),
                $"Selecting unstable nodal planes, iteration {iteration}.", friction));
            var normal = BuildNormalEquations(observations, planes, selected);
            if (normal is null) return null;
            var solved = SolveSymmetric(normal.Value.H, normal.Value.B);
            if (solved.X is null || solved.ConditionNumber > options.MaximumConditionNumber) return null;
            tensor = TensorFromFive(solved.X);
            condition = solved.ConditionNumber;
            var next = new int[selected.Length];
            for (int i = 0; i < next.Length; i++)
            {
                double i1 = Instability(planes[i].First.Normal, tensor, friction);
                double i2 = Instability(planes[i].Second.Normal, tensor, friction);
                next[i] = i2 > i1 ? 1 : 0;
            }
            if (next.SequenceEqual(selected))
            {
                converged = true;
                selected = next;
                break;
            }
            if (!visited.Add(SelectionKey(next)))
            {
                // Real catalogues can enter a short discrete plane-selection cycle while the
                // continuous stress tensor is already stable. Keep the current selection,
                // whose tensor was solved on the same iteration, rather than pairing a tensor
                // with the following cycle state. This is deterministic and explicitly
                // reported in QC instead of exhausting an arbitrary iteration limit.
                converged = true;
                cycleResolved = true;
                break;
            }
            selected = next;
        }
        if (tensor is null) return null;
        double meanInstability = selected.Select((p, i) => Instability(planes[i][p].Normal, tensor, friction)).Average();
        double meanMisfit = selected.Select((p, i) => SlipMisfit(planes[i][p], tensor)).Average();
        return new CoreSolution(tensor, selected, friction, iterations, converged, cycleResolved,
            condition, meanInstability, meanMisfit, planes);
    }

    private static string SelectionKey(IReadOnlyList<int> selected) =>
        string.Create(selected.Count, selected, static (span, values) =>
        {
            for (int i = 0; i < values.Count; i++) span[i] = values[i] == 0 ? '0' : '1';
        });

    private static BootstrapSummary Bootstrap(
        FocalMechanismObservation[] observations, StressInversionOptions options, CoreSolution reference,
        IProgress<StressInversionProgress>? progress, CancellationToken cancellationToken)
    {
        if (options.BootstrapRealizations <= 0) return new BootstrapSummary();
        var random = new Random(options.RandomSeed ^ 0x51A7_2024);
        var tensors = new List<TensorState>(options.BootstrapRealizations);
        for (int realization = 0; realization < options.BootstrapRealizations; realization++)
        {
            cancellationToken.ThrowIfCancellationRequested();
            progress?.Report(new StressInversionProgress(StressInversionProgressStage.Bootstrap,
                realization + 1, options.BootstrapRealizations,
                $"Running uncertainty realization {realization + 1}/{options.BootstrapRealizations}.",
                reference.Friction));
            var sample = new FocalMechanismObservation[observations.Length];
            for (int i = 0; i < sample.Length; i++)
            {
                var source = observations[random.Next(observations.Length)];
                sample[i] = Perturb(source, random);
            }
            var planes = sample.Select(ToPlanes).ToArray();
            var core = SolveIterative(sample, planes, reference.Friction,
                CloneWithoutBootstrap(options, optimizeFriction: false), options.RandomSeed + realization + 1,
                progress: null, cancellationToken);
            if (core is not null && core.Converged) tensors.Add(core.Tensor);
        }
        return new BootstrapSummary(tensors, reference.Tensor);
    }

    private static FocalMechanismObservation Perturb(FocalMechanismObservation source, Random random)
    {
        double sigma = Math.Max(0, source.AngularUncertaintyDeg);
        Plane original = ToPlanes(source).First;
        Plane perturbed = sigma == 0 ? original : RotatePlane(original, sigma * DegToRad * Gaussian(random), random);
        var (strike, dip, rake) = ToStrikeDipRake(perturbed);
        return new FocalMechanismObservation
        {
            EventId = source.EventId, LatitudeDeg = source.LatitudeDeg, LongitudeDeg = source.LongitudeDeg,
            DepthKm = source.DepthKm, Weight = source.Weight, AngularUncertaintyDeg = source.AngularUncertaintyDeg,
            DoubleCoupleFraction = source.DoubleCoupleFraction, Provider = source.Provider,
            SourceId = source.SourceId, Provenance = source.Provenance,
            NodalPlane1StrikeDeg = strike, NodalPlane1DipDeg = dip, NodalPlane1RakeDeg = rake
        };
    }

    private static FocalMechanismObservation WithWeight(FocalMechanismObservation source, double weight) => new()
    {
        EventId = source.EventId,
        LatitudeDeg = source.LatitudeDeg,
        LongitudeDeg = source.LongitudeDeg,
        DepthKm = source.DepthKm,
        NodalPlane1StrikeDeg = source.NodalPlane1StrikeDeg,
        NodalPlane1DipDeg = source.NodalPlane1DipDeg,
        NodalPlane1RakeDeg = source.NodalPlane1RakeDeg,
        NodalPlane2StrikeDeg = source.NodalPlane2StrikeDeg,
        NodalPlane2DipDeg = source.NodalPlane2DipDeg,
        NodalPlane2RakeDeg = source.NodalPlane2RakeDeg,
        Weight = weight,
        AngularUncertaintyDeg = source.AngularUncertaintyDeg,
        DoubleCoupleFraction = source.DoubleCoupleFraction,
        Provider = source.Provider,
        SourceId = source.SourceId,
        Provenance = source.Provenance
    };

    private static Plane RotatePlane(Plane plane, double angleRad, Random random)
    {
        var axis = new V3(Gaussian(random), Gaussian(random), Gaussian(random)).Unit;
        V3 Rotate(V3 v)
        {
            double c = Math.Cos(angleRad), s = Math.Sin(angleRad);
            return v * c + Cross(axis, v) * s + axis * (Dot(axis, v) * (1 - c));
        }
        V3 normal = Rotate(plane.Normal).Unit;
        V3 slip = Rotate(plane.Slip).Unit;
        if (normal.Z > 0) { normal = -normal; slip = -slip; }
        return new Plane(normal, slip);
    }

    private static (double Strike, double Dip, double Rake) ToStrikeDipRake(Plane plane)
    {
        V3 normal = plane.Normal, slip = plane.Slip;
        if (normal.Z > 0) { normal = -normal; slip = -slip; }
        double dip = Math.Acos(Math.Clamp(-normal.Z, -1, 1));
        double strike = Math.Atan2(-normal.X, normal.Y);
        var alongStrike = new V3(Math.Cos(strike), Math.Sin(strike), 0);
        var rake90 = new V3(Math.Cos(dip) * Math.Sin(strike), -Math.Cos(dip) * Math.Cos(strike), -Math.Sin(dip));
        double rake = Math.Atan2(Dot(slip, rake90), Dot(slip, alongStrike));
        return (Normalize360(strike / DegToRad), dip / DegToRad, Normalize180(rake / DegToRad));
    }

    private static (double[,] H, double[] B)? BuildNormalEquations(
        IReadOnlyList<FocalMechanismObservation> observations, PlanePair[] planes, int[] selected)
    {
        var h = new double[5, 5];
        var b = new double[5];
        for (int i = 0; i < observations.Count; i++)
            Accumulate(h, b, planes[i][selected[i]], EffectiveWeight(observations[i]));
        return HasFinite(h, b) ? (h, b) : null;
    }

    private static (double[,] H, double[] B)? BuildNormalEquations(
        IReadOnlyList<FocalMechanismObservation> observations,
        IReadOnlyDictionary<string, SelectedFaultPlaneResult> selected)
    {
        var h = new double[5, 5];
        var b = new double[5];
        foreach (var observation in observations)
        {
            if (!selected.TryGetValue(observation.EventId, out var choice)) continue;
            var pair = ToPlanes(observation);
            Accumulate(h, b, pair[Math.Clamp(choice.NodalPlane - 1, 0, 1)], EffectiveWeight(observation));
        }
        return HasFinite(h, b) ? (h, b) : null;
    }

    private static void Accumulate(double[,] h, double[] b, Plane plane, double weight)
    {
        for (int component = 0; component < 3; component++)
        {
            var g = new double[5];
            for (int k = 0; k < 5; k++)
            {
                var basis = Basis(k);
                var traction = Multiply(basis, plane.Normal);
                var shear = traction - plane.Normal * Dot(traction, plane.Normal);
                g[k] = shear[component];
            }
            // Wallace–Bott: slip is parallel to the shear traction of the tension-positive stress on the
            // footwall-to-hanging-wall normal, i.e. antiparallel to that of the compression-positive
            // tensor solved for here. (A +slip convention and then called the largest eigenvalue σ1,
            // which swapped σ1 and σ3.)
            double datum = -plane.Slip[component];
            for (int r = 0; r < 5; r++)
            {
                b[r] += weight * g[r] * datum;
                for (int c = 0; c < 5; c++) h[r, c] += weight * g[r] * g[c];
            }
        }
    }

    private static IReadOnlyList<SelectedFaultPlaneResult> SelectedResults(
        FocalMechanismObservation[] observations, CoreSolution solution)
        => observations.Select((o, i) => new SelectedFaultPlaneResult
        {
            EventId = o.EventId, NodalPlane = solution.Selected[i] + 1,
            Instability = Instability(solution.Planes[i][solution.Selected[i]].Normal, solution.Tensor, solution.Friction),
            SlipMisfitDeg = SlipMisfit(solution.Planes[i][solution.Selected[i]], solution.Tensor)
        }).ToArray();

    private static StressTensorResult ToResult(TensorState tensor, BootstrapSummary bootstrap)
        => TensorFromState(tensor, bootstrap);

    private static StressTensorResult TensorFromComponents(double[] c, BootstrapSummary? bootstrap)
        => TensorFromState(TensorFromFive(c), bootstrap);

    private static StressTensorResult TensorFromState(TensorState tensor, BootstrapSummary? bootstrap)
    {
        double r = tensor.ShapeRatio;
        double[] eigenvalues = DeviatoricEigenvalues(r);
        double[,] canonical = Recompose(eigenvalues, tensor.Axes);
        return new StressTensorResult
        {
            // Internal frame is north–east–down (the frame of the strike/dip/rake formulas); the
            // result is reported east–north–down as documented.
            Sxx = canonical[1, 1], Syy = canonical[0, 0], Szz = canonical[2, 2],
            Sxy = canonical[0, 1], Sxz = canonical[1, 2], Syz = canonical[0, 2],
            Sigma1 = AxisResult(tensor.Axis(0), bootstrap?.Cone(0) ?? double.NaN),
            Sigma2 = AxisResult(tensor.Axis(1), bootstrap?.Cone(1) ?? double.NaN),
            Sigma3 = AxisResult(tensor.Axis(2), bootstrap?.Cone(2) ?? double.NaN),
            ShapeRatio = r, ShapeRatioLower95 = bootstrap?.RLower ?? double.NaN,
            ShapeRatioUpper95 = bootstrap?.RUpper ?? double.NaN
        };
    }

    private static StressAxisResult AxisResult(V3 v, double cone)
    {
        if (v.Z < 0) v = -v;
        // v is north–east–down.
        double az = Normalize360(Math.Atan2(v.Y, v.X) / DegToRad);
        double plunge = Math.Asin(Math.Clamp(v.Z, -1, 1)) / DegToRad;
        return new StressAxisResult { East = v.Y, North = v.X, Down = v.Z, AzimuthDeg = az, PlungeDeg = plunge, ConfidenceCone95Deg = cone };
    }

    private static TensorState TensorFromFive(IReadOnlyList<double> c)
    {
        var matrix = new double[,]
        {
            { c[0], c[2], c[3] },
            { c[2], c[1], c[4] },
            { c[3], c[4], -c[0] - c[1] }
        };
        var (values, vectors) = JacobiEigen(matrix);
        int[] order = Enumerable.Range(0, 3).OrderByDescending(i => values[i]).ToArray();
        double l1 = values[order[0]], l2 = values[order[1]], l3 = values[order[2]];
        double range = l1 - l3;
        double shape = range > 1e-14 ? Math.Clamp((l1 - l2) / range, 0, 1) : 0.5;
        var axes = new double[3, 3];
        for (int column = 0; column < 3; column++)
        for (int row = 0; row < 3; row++) axes[row, column] = vectors[row, order[column]];
        return new TensorState(matrix, axes, shape);
    }

    private static double Instability(V3 normal, TensorState tensor, double friction)
    {
        double[] eigenvalues = DeviatoricEigenvalues(tensor.ShapeRatio);
        var canonical = Recompose(eigenvalues, tensor.Axes);
        V3 traction = Multiply(canonical, normal);
        double sigma = Dot(normal, traction);
        double tau = (traction - normal * sigma).Norm;
        double denominator = friction + Math.Sqrt(1 + friction * friction);
        // Vavryčuk (2014) instability, compression positive, σ1 − σ3 = 2: I = (τ − μ(σn − σ1)) / (μ + √(1+μ²)).
        return Math.Clamp((tau - friction * (sigma - eigenvalues[0])) / denominator, 0, 1);
    }

    private static double SlipMisfit(Plane plane, TensorState tensor)
    {
        double[] eigenvalues = DeviatoricEigenvalues(tensor.ShapeRatio);
        var stress = Recompose(eigenvalues, tensor.Axes);
        V3 traction = Multiply(stress, plane.Normal);
        V3 shear = traction - plane.Normal * Dot(traction, plane.Normal);
        if (shear.Norm < 1e-14) return 90;
        // Compression-positive tensor: the predicted slip is along −shear (see Accumulate).
        return Math.Acos(Math.Clamp(-Dot(shear.Unit, plane.Slip), -1, 1)) / DegToRad;
    }

    private static double MeanMisfit(
        IReadOnlyList<FocalMechanismObservation> observations,
        IReadOnlyDictionary<string, SelectedFaultPlaneResult> selected, StressTensorResult tensor)
    {
        var state = TensorFromFive(Components(tensor));
        var values = new List<double>();
        foreach (var observation in observations)
            if (selected.TryGetValue(observation.EventId, out var choice))
                values.Add(SlipMisfit(ToPlanes(observation)[choice.NodalPlane - 1], state));
        return values.Count == 0 ? double.NaN : values.Average();
    }

    private static PlanePair ToPlanes(FocalMechanismObservation observation)
    {
        Plane first = FaultVectors(observation.NodalPlane1StrikeDeg, observation.NodalPlane1DipDeg, observation.NodalPlane1RakeDeg);
        Plane second = observation.NodalPlane2StrikeDeg is double strike && observation.NodalPlane2DipDeg is double dip && observation.NodalPlane2RakeDeg is double rake
            ? FaultVectors(strike, dip, rake)
            : new Plane(first.Slip, first.Normal);
        return new PlanePair(first, second);
    }

    private static Plane FaultVectors(double strikeDeg, double dipDeg, double rakeDeg)
    {
        double strike = strikeDeg * DegToRad, dip = dipDeg * DegToRad, rake = rakeDeg * DegToRad;
        var normal = new V3(-Math.Sin(dip) * Math.Sin(strike), Math.Sin(dip) * Math.Cos(strike), -Math.Cos(dip)).Unit;
        var slip = new V3(
            Math.Cos(rake) * Math.Cos(strike) + Math.Sin(rake) * Math.Cos(dip) * Math.Sin(strike),
            Math.Cos(rake) * Math.Sin(strike) - Math.Sin(rake) * Math.Cos(dip) * Math.Cos(strike),
            -Math.Sin(rake) * Math.Sin(dip)).Unit;
        return new Plane(normal, slip);
    }

    private static bool IsValid(FocalMechanismObservation o)
        => !string.IsNullOrWhiteSpace(o.EventId) && double.IsFinite(o.NodalPlane1StrikeDeg) &&
           o.NodalPlane1DipDeg > 0 && o.NodalPlane1DipDeg <= 90 && double.IsFinite(o.NodalPlane1RakeDeg) &&
           double.IsFinite(o.Weight) && o.Weight > 0 && double.IsFinite(o.DoubleCoupleFraction) &&
           o.DoubleCoupleFraction > 0 && o.DoubleCoupleFraction <= 1 &&
           double.IsFinite(o.AngularUncertaintyDeg) && o.AngularUncertaintyDeg >= 0 &&
           ((o.NodalPlane2StrikeDeg is null && o.NodalPlane2DipDeg is null && o.NodalPlane2RakeDeg is null) ||
            (o.NodalPlane2StrikeDeg is double s && double.IsFinite(s) &&
             o.NodalPlane2DipDeg is double d && d > 0 && d <= 90 &&
             o.NodalPlane2RakeDeg is double r && double.IsFinite(r)));

    private static double EffectiveWeight(FocalMechanismObservation o)
        => o.Weight * Math.Clamp(o.DoubleCoupleFraction, 0, 1);

    private static StressInversionResult Failure(int count, string message) => new()
    {
        MechanismCount = count, QualityStatus = "FAILED-QC", QualityMessages = new[] { message }
    };

    private static StressInversionOptions CloneWithoutBootstrap(StressInversionOptions o, bool? optimizeFriction = null) => new()
    {
        MaximumIterations = o.MaximumIterations, FrictionCoefficient = o.FrictionCoefficient,
        OptimizeFriction = optimizeFriction ?? o.OptimizeFriction, MinimumFriction = o.MinimumFriction,
        MaximumFriction = o.MaximumFriction, FrictionStep = o.FrictionStep, BootstrapRealizations = 0,
        RandomSeed = o.RandomSeed, MinimumMechanisms = o.MinimumMechanisms,
        CertificationMinimumMechanisms = o.CertificationMinimumMechanisms, MaximumConditionNumber = o.MaximumConditionNumber
    };

    private static IEnumerable<double> InclusiveRange(double minimum, double maximum, double step)
    {
        if (!double.IsFinite(minimum) || !double.IsFinite(maximum) || !double.IsFinite(step) || step <= 0 || maximum < minimum)
            yield break;
        for (int i = 0; ; i++)
        {
            double value = minimum + i * step;
            if (value > maximum + step * 1e-8) break;
            yield return value;
        }
    }

    private static (double[]? X, double ConditionNumber) SolveSymmetric(double[,] a, double[] b)
    {
        var (eigenvalues, _) = JacobiEigen(a);
        double max = eigenvalues.Select(Math.Abs).Max();
        double min = eigenvalues.Select(Math.Abs).Where(v => v > max * 1e-13).DefaultIfEmpty(0).Min();
        double condition = min > 0 ? max / min : double.PositiveInfinity;
        int n = b.Length;
        var m = new double[n, n + 1];
        for (int r = 0; r < n; r++)
        {
            for (int c = 0; c < n; c++) m[r, c] = a[r, c];
            m[r, n] = b[r];
        }
        for (int p = 0; p < n; p++)
        {
            int pivot = p;
            for (int r = p + 1; r < n; r++) if (Math.Abs(m[r, p]) > Math.Abs(m[pivot, p])) pivot = r;
            if (Math.Abs(m[pivot, p]) <= Math.Max(1e-14, max * 1e-13)) return (null, double.PositiveInfinity);
            if (pivot != p) for (int c = p; c <= n; c++) (m[p, c], m[pivot, c]) = (m[pivot, c], m[p, c]);
            double diagonal = m[p, p];
            for (int c = p; c <= n; c++) m[p, c] /= diagonal;
            for (int r = 0; r < n; r++)
            {
                if (r == p) continue;
                double factor = m[r, p];
                for (int c = p; c <= n; c++) m[r, c] -= factor * m[p, c];
            }
        }
        var x = new double[n];
        for (int i = 0; i < n; i++) x[i] = m[i, n];
        return (x, condition);
    }

    // Shared by the source-tensor decomposition so one audited symmetric
    // eigensolver instead of introducing a second numerical implementation.
    internal static (double[] Values, double[,] Vectors) JacobiEigen(double[,] source)
    {
        int n = source.GetLength(0);
        var a = (double[,])source.Clone();
        var v = new double[n, n];
        for (int i = 0; i < n; i++) v[i, i] = 1;
        for (int iteration = 0; iteration < 100 * n * n; iteration++)
        {
            int p = 0, q = 1;
            double largest = 0;
            for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
                if (Math.Abs(a[i, j]) > largest) { largest = Math.Abs(a[i, j]); p = i; q = j; }
            if (largest < 1e-13) break;
            double phi = 0.5 * Math.Atan2(2 * a[p, q], a[q, q] - a[p, p]);
            double c = Math.Cos(phi), s = Math.Sin(phi);
            for (int k = 0; k < n; k++)
            {
                double apk = a[p, k], aqk = a[q, k];
                a[p, k] = c * apk - s * aqk;
                a[q, k] = s * apk + c * aqk;
            }
            for (int k = 0; k < n; k++)
            {
                double akp = a[k, p], akq = a[k, q];
                a[k, p] = c * akp - s * akq;
                a[k, q] = s * akp + c * akq;
                double vkp = v[k, p], vkq = v[k, q];
                v[k, p] = c * vkp - s * vkq;
                v[k, q] = s * vkp + c * vkq;
            }
        }
        return (Enumerable.Range(0, n).Select(i => a[i, i]).ToArray(), v);
    }

    private static (double[] Solution, int Iterations, bool Converged) ConjugateGradient(
        Func<double[], double[]> apply, double[] b, double[] initial, int maximumIterations, double tolerance,
        IProgress<StressInversionProgress>? progress, CancellationToken cancellationToken)
    {
        var x = (double[])initial.Clone();
        var ax = apply(x);
        var r = b.Zip(ax, (bi, ai) => bi - ai).ToArray();
        var p = (double[])r.Clone();
        double rr = Dot(r, r), bNorm = Math.Sqrt(Math.Max(Dot(b, b), 1e-30));
        if (Math.Sqrt(rr) / bNorm <= tolerance) return (x, 0, true);
        for (int iteration = 1; iteration <= Math.Max(1, maximumIterations); iteration++)
        {
            cancellationToken.ThrowIfCancellationRequested();
            var ap = apply(p);
            double denominator = Dot(p, ap);
            if (!(denominator > 0) || !double.IsFinite(denominator)) return (x, iteration, false);
            double alpha = rr / denominator;
            for (int i = 0; i < x.Length; i++) { x[i] += alpha * p[i]; r[i] -= alpha * ap[i]; }
            double next = Dot(r, r);
            progress?.Report(new StressInversionProgress(StressInversionProgressStage.SpatialConjugateGradient,
                iteration, Math.Max(1, maximumIterations), $"Solving damped spatial system, CG iteration {iteration}.",
                ResidualNorm: Math.Sqrt(next) / bNorm));
            if (Math.Sqrt(next) / bNorm <= tolerance) return (x, iteration, true);
            double beta = next / rr;
            for (int i = 0; i < p.Length; i++) p[i] = r[i] + beta * p[i];
            rr = next;
        }
        return (x, maximumIterations, false);
    }

    private static double[,] Basis(int component) => component switch
    {
        0 => new double[,] { { 1, 0, 0 }, { 0, 0, 0 }, { 0, 0, -1 } },
        1 => new double[,] { { 0, 0, 0 }, { 0, 1, 0 }, { 0, 0, -1 } },
        2 => new double[,] { { 0, 1, 0 }, { 1, 0, 0 }, { 0, 0, 0 } },
        3 => new double[,] { { 0, 0, 1 }, { 0, 0, 0 }, { 1, 0, 0 } },
        _ => new double[,] { { 0, 0, 0 }, { 0, 0, 1 }, { 0, 1, 0 } }
    };

    private static double[,] Recompose(double[] eigenvalues, double[,] axes)
    {
        var result = new double[3, 3];
        for (int k = 0; k < 3; k++)
        for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) result[i, j] += eigenvalues[k] * axes[i, k] * axes[j, k];
        return result;
    }

    private static double[] DeviatoricEigenvalues(double shapeRatio)
    {
        // Preserve sigma1-sigma3=2 while removing the hydrostatic component
        // from the convenient (1, 1-2R, -1) representation.
        double sigma1 = 2 * (1 + shapeRatio) / 3;
        return new[] { sigma1, sigma1 - 2 * shapeRatio, sigma1 - 2 };
    }

    private static V3 Multiply(double[,] matrix, V3 v) => new(
        matrix[0, 0] * v.X + matrix[0, 1] * v.Y + matrix[0, 2] * v.Z,
        matrix[1, 0] * v.X + matrix[1, 1] * v.Y + matrix[1, 2] * v.Z,
        matrix[2, 0] * v.X + matrix[2, 1] * v.Y + matrix[2, 2] * v.Z);

    /// <summary>The five internal (north–east–down) components [nn, ee, ne, nd, ed] of a reported (east–north–down) tensor.</summary>
    private static double[] Components(StressTensorResult t) => new[] { t.Syy, t.Sxx, t.Sxy, t.Syz, t.Sxz };
    private static bool HasFinite(double[,] h, double[] b) => h.Cast<double>().All(double.IsFinite) && b.All(double.IsFinite);
    private static double Dot(V3 a, V3 b) => a.X * b.X + a.Y * b.Y + a.Z * b.Z;
    private static V3 Cross(V3 a, V3 b) => new(a.Y * b.Z - a.Z * b.Y, a.Z * b.X - a.X * b.Z, a.X * b.Y - a.Y * b.X);
    private static double Dot(double[] a, double[] b) => a.Zip(b, (x, y) => x * y).Sum();
    private static double Gaussian(Random random)
    {
        double u1 = Math.Max(double.Epsilon, random.NextDouble());
        return Math.Sqrt(-2 * Math.Log(u1)) * Math.Cos(2 * Math.PI * random.NextDouble());
    }

    private static double SpatialKernelWeight(FocalMechanismObservation o, double longitude, double latitude, double depth, StressGridRequest r)
    {
        if (!double.IsFinite(o.LongitudeDeg) || !double.IsFinite(o.LatitudeDeg) || !double.IsFinite(o.DepthKm)) return 0;
        double north = (o.LatitudeDeg - latitude) * 111.195;
        double east = (o.LongitudeDeg - longitude) * 111.195 * Math.Cos(latitude * DegToRad);
        double horizontal = Math.Sqrt(east * east + north * north);
        double vertical = Math.Abs(o.DepthKm - depth);
        if (horizontal > r.HorizontalSearchRadiusKm || vertical > r.VerticalSearchRadiusKm) return 0;
        if (r.SpatialWeighting == StressSpatialWeighting.Uniform) return 1;
        double cutoff = r.SearchRadiusStandardDeviations;
        double normalizedSquared = Math.Pow(cutoff * horizontal / r.HorizontalSearchRadiusKm, 2)
                                   + Math.Pow(cutoff * vertical / r.VerticalSearchRadiusKm, 2);
        return Math.Exp(-0.5 * normalizedSquared);
    }

    private static double EffectiveSampleCount(IReadOnlyList<FocalMechanismObservation> observations)
    {
        double sum = 0, sumSquared = 0;
        foreach (var observation in observations)
        {
            double weight = EffectiveWeight(observation);
            sum += weight;
            sumSquared += weight * weight;
        }
        return sumSquared > 0 ? sum * sum / sumSquared : 0;
    }

    private static void ValidateGrid(StressGridRequest r)
    {
        if (r.Observations is null) throw new ArgumentNullException(nameof(r.Observations));
        if (!(r.LongitudeStepDeg > 0) || !(r.LatitudeStepDeg > 0) || !(r.DepthStepKm > 0))
            throw new ArgumentOutOfRangeException(nameof(r), "All grid steps must be positive.");
        if (r.MaximumLongitudeDeg < r.MinimumLongitudeDeg || r.MaximumLatitudeDeg < r.MinimumLatitudeDeg || r.MaximumDepthKm < r.MinimumDepthKm)
            throw new ArgumentException("Grid maxima must not be smaller than minima.", nameof(r));
        if (!(r.HorizontalSearchRadiusKm > 0) || !(r.VerticalSearchRadiusKm > 0))
            throw new ArgumentOutOfRangeException(nameof(r), "Search radii must be positive.");
        if (!double.IsFinite(r.SearchRadiusStandardDeviations) || r.SearchRadiusStandardDeviations <= 0)
            throw new ArgumentOutOfRangeException(nameof(r), "SearchRadiusStandardDeviations must be positive and finite.");
    }

    private static double[] AxisValues(double minimum, double maximum, double step)
    {
        int count = checked((int)Math.Floor((maximum - minimum) / step + 1e-9) + 1);
        if (count > 1_000) throw new ArgumentOutOfRangeException(nameof(step), "A grid axis may contain at most 1000 nodes.");
        return Enumerable.Range(0, count).Select(i => minimum + i * step).ToArray();
    }

    private static int Flat(int x, int y, int z, int nx, int ny) => (z * ny + y) * nx + x;
    private static void Unflat(int flat, int nx, int ny, out int x, out int y, out int z)
    {
        x = flat % nx; int q = flat / nx; y = q % ny; z = q / ny;
    }
    private static double Normalize360(double value) { value %= 360; return value < 0 ? value + 360 : value; }
    private static double Normalize180(double value) { value = Normalize360(value); return value > 180 ? value - 360 : value; }

    private readonly record struct V3(double X, double Y, double Z)
    {
        public double Norm => Math.Sqrt(X * X + Y * Y + Z * Z);
        public V3 Unit => Norm > 0 ? this * (1 / Norm) : this;
        public double this[int i] => i == 0 ? X : i == 1 ? Y : Z;
        public static V3 operator +(V3 a, V3 b) => new(a.X + b.X, a.Y + b.Y, a.Z + b.Z);
        public static V3 operator -(V3 a, V3 b) => new(a.X - b.X, a.Y - b.Y, a.Z - b.Z);
        public static V3 operator -(V3 a) => new(-a.X, -a.Y, -a.Z);
        public static V3 operator *(V3 a, double k) => new(a.X * k, a.Y * k, a.Z * k);
    }

    private readonly record struct Plane(V3 Normal, V3 Slip);
    private readonly record struct PlanePair(Plane First, Plane Second)
    {
        public Plane this[int i] => i == 0 ? First : Second;
    }
    private sealed record TensorState(double[,] Matrix, double[,] Axes, double ShapeRatio)
    {
        public V3 Axis(int i) => new(Axes[0, i], Axes[1, i], Axes[2, i]);
    }
    private sealed record CoreSolution(TensorState Tensor, int[] Selected, double Friction, int Iterations,
        bool Converged, bool CycleResolved, double ConditionNumber, double MeanInstability,
        double MeanMisfitDeg, PlanePair[] Planes);
    private sealed record GridWork(FocalMechanismObservation[] Observations,
        IReadOnlyDictionary<string, SelectedFaultPlaneResult> Selected, double[,] H, double[] B,
        double[] Initial, double ConditionNumber, string QualityStatus, double EffectiveMechanismCount);

    private sealed class BootstrapSummary
    {
        private readonly List<TensorState> _tensors;
        private readonly TensorState? _reference;
        public BootstrapSummary() { _tensors = new(); }
        public BootstrapSummary(List<TensorState> tensors, TensorState reference) { _tensors = tensors; _reference = reference; }
        public int Count => _tensors.Count;
        public double RLower => Percentile(_tensors.Select(t => t.ShapeRatio), 0.025);
        public double RUpper => Percentile(_tensors.Select(t => t.ShapeRatio), 0.975);
        public double Cone(int axis)
        {
            if (_reference is null) return double.NaN;
            V3 reference = _reference.Axis(axis);
            return Percentile(_tensors.Select(t => Math.Acos(Math.Clamp(Math.Abs(Dot(reference, t.Axis(axis))), -1, 1)) / DegToRad), 0.95);
        }
        private static double Percentile(IEnumerable<double> source, double probability)
        {
            double[] values = source.Where(double.IsFinite).OrderBy(v => v).ToArray();
            if (values.Length == 0) return double.NaN;
            double position = probability * (values.Length - 1);
            int lower = (int)Math.Floor(position), upper = (int)Math.Ceiling(position);
            return values[lower] + (values[upper] - values[lower]) * (position - lower);
        }
    }
}
