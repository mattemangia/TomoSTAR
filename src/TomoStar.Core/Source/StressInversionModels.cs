// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

namespace TomoStar.Core.Source;

public enum StressInversionProgressStage
{
    FrictionSearch,
    NodalPlaneIteration,
    Bootstrap,
    SpatialGridNode,
    SpatialConjugateGradient,
    Completed
}

public sealed record StressInversionProgress(
    StressInversionProgressStage Stage,
    int Current,
    int Total,
    string Message,
    double? FrictionCoefficient = null,
    double? LongitudeDeg = null,
    double? LatitudeDeg = null,
    double? DepthKm = null,
    double? ResidualNorm = null)
{
    public double Fraction => Total > 0 ? Math.Clamp((double)Current / Total, 0, 1) : 0;
}

/// <summary>
/// One double-couple focal mechanism used by the stress inversion.
/// Angles follow the right-hand-rule convention; depth is positive down.
/// Both nodal planes are explicit so imported QuakeML solutions retain their
/// published geometry and provenance.
/// </summary>
public sealed class FocalMechanismObservation
{
    public string EventId { get; init; } = string.Empty;
    public double LatitudeDeg { get; init; } = double.NaN;
    public double LongitudeDeg { get; init; } = double.NaN;
    public double DepthKm { get; init; } = double.NaN;
    public double NodalPlane1StrikeDeg { get; init; }
    public double NodalPlane1DipDeg { get; init; }
    public double NodalPlane1RakeDeg { get; init; }
    public double? NodalPlane2StrikeDeg { get; init; }
    public double? NodalPlane2DipDeg { get; init; }
    public double? NodalPlane2RakeDeg { get; init; }
    public double Weight { get; init; } = 1.0;
    public double AngularUncertaintyDeg { get; init; }
    public double DoubleCoupleFraction { get; init; } = 1.0;
    public string Provider { get; init; } = string.Empty;
    public string SourceId { get; init; } = string.Empty;
    public string Provenance { get; init; } = string.Empty;
}

public sealed class StressInversionOptions
{
    public int MaximumIterations { get; init; } = 20;
    public double FrictionCoefficient { get; init; } = 0.6;
    public bool OptimizeFriction { get; init; } = true;
    public double MinimumFriction { get; init; } = 0.2;
    public double MaximumFriction { get; init; } = 1.2;
    public double FrictionStep { get; init; } = 0.05;
    public int BootstrapRealizations { get; init; } = 1_000;
    public int RandomSeed { get; init; } = 17_031;
    public int MinimumMechanisms { get; init; } = 5;
    public int CertificationMinimumMechanisms { get; init; } = 20;
    public double MaximumConditionNumber { get; init; } = 1e10;
}

public sealed class StressAxisResult
{
    public double East { get; init; }
    public double North { get; init; }
    public double Down { get; init; }
    public double AzimuthDeg { get; init; }
    public double PlungeDeg { get; init; }
    public double ConfidenceCone95Deg { get; init; } = double.NaN;
}

public sealed class SelectedFaultPlaneResult
{
    public string EventId { get; init; } = string.Empty;
    public int NodalPlane { get; init; }
    public double Instability { get; init; }
    public double SlipMisfitDeg { get; init; }
}

public sealed class StressTensorResult
{
    /// <summary>Reduced deviatoric tensor in END coordinates, compression positive.</summary>
    public double Sxx { get; init; }
    public double Syy { get; init; }
    public double Szz { get; init; }
    public double Sxy { get; init; }
    public double Sxz { get; init; }
    public double Syz { get; init; }
    public required StressAxisResult Sigma1 { get; init; }
    public required StressAxisResult Sigma2 { get; init; }
    public required StressAxisResult Sigma3 { get; init; }
    /// <summary>R=(sigma1-sigma2)/(sigma1-sigma3), in [0,1].</summary>
    public double ShapeRatio { get; init; }
    public double ShapeRatioLower95 { get; init; } = double.NaN;
    public double ShapeRatioUpper95 { get; init; } = double.NaN;
}

public sealed class StressInversionResult
{
    public StressTensorResult? Tensor { get; init; }
    public int MechanismCount { get; init; }
    public int Iterations { get; init; }
    public int BootstrapCompleted { get; init; }
    public double FrictionCoefficient { get; init; }
    public double MeanFaultInstability { get; init; } = double.NaN;
    public double MeanSlipMisfitDeg { get; init; } = double.NaN;
    public double ConditionNumber { get; init; } = double.PositiveInfinity;
    public IReadOnlyList<SelectedFaultPlaneResult> SelectedFaultPlanes { get; init; }
        = Array.Empty<SelectedFaultPlaneResult>();
    public string QualityStatus { get; init; } = "FAILED-QC";
    public IReadOnlyList<string> QualityMessages { get; init; } = Array.Empty<string>();
    public string Method { get; init; } = "Michael (1984) + Vavrycuk (2014)";
    public string CoordinateConvention { get; init; }
        = "END Cartesian; azimuth clockwise from north; plunge positive down; compression positive";
}

public sealed class StressGridRequest
{
    public required IReadOnlyList<FocalMechanismObservation> Observations { get; init; }
    public double MinimumLongitudeDeg { get; init; }
    public double MaximumLongitudeDeg { get; init; }
    public double LongitudeStepDeg { get; init; }
    public double MinimumLatitudeDeg { get; init; }
    public double MaximumLatitudeDeg { get; init; }
    public double LatitudeStepDeg { get; init; }
    public double MinimumDepthKm { get; init; }
    public double MaximumDepthKm { get; init; }
    public double DepthStepKm { get; init; }
    public double HorizontalSearchRadiusKm { get; init; } = 25;
    public double VerticalSearchRadiusKm { get; init; } = 9;
    /// <summary>
    /// Spatial data localization. Gaussian is the defensible default; Uniform
    /// reproduces the legacy top-hat window and can strongly correlate adjacent
    /// depth nodes when the radius is large relative to the grid step.
    /// </summary>
    public StressSpatialWeighting SpatialWeighting { get; init; } = StressSpatialWeighting.Gaussian;
    /// <summary>Number of Gaussian standard deviations at each search-radius cutoff.</summary>
    public double SearchRadiusStandardDeviations { get; init; } = 3.0;
    public double NeighborDamping { get; init; } = 1.0;
    public int MaximumConjugateGradientIterations { get; init; } = 500;
    public StressInversionOptions InversionOptions { get; init; } = new();
}

public enum StressSpatialWeighting
{
    Gaussian,
    Uniform
}

public sealed class StressGridNodeResult
{
    public double LongitudeDeg { get; init; }
    public double LatitudeDeg { get; init; }
    public double DepthKm { get; init; }
    public int MechanismCount { get; init; }
    public double EffectiveMechanismCount { get; init; }
    public StressTensorResult? Tensor { get; init; }
    public double MeanSlipMisfitDeg { get; init; } = double.NaN;
    public string QualityStatus { get; init; } = "FAILED-QC";
}

public sealed class StressGridResult
{
    public IReadOnlyList<StressGridNodeResult> Nodes { get; init; }
        = Array.Empty<StressGridNodeResult>();
    public int LongitudeCount { get; init; }
    public int LatitudeCount { get; init; }
    public int DepthCount { get; init; }
    public double Damping { get; init; }
    public int ConjugateGradientIterations { get; init; }
    public bool Converged { get; init; }
    public string QualityStatus { get; init; } = "FAILED-QC";
    public IReadOnlyList<string> QualityMessages { get; init; } = Array.Empty<string>();
    public string Method { get; init; }
        = "Spatially damped simultaneous reduced-stress inversion (Hardebeck and Michael, 2006)";
}
