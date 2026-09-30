// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

namespace TomoStar.Core.Compute;

/// <summary>
/// Process-wide compute limits set from the user preferences: how many CPU threads the heavy
/// parallel loops (eikonal tables, ray tracing, relocation) may use, so a workstation stays usable
/// while an inversion runs.
/// </summary>
public static class ComputeSettings
{
    /// <summary>0 = all logical processors.</summary>
    public static int MaxThreads { get; set; }

    public static int Threads => MaxThreads > 0 ? Math.Min(MaxThreads, Environment.ProcessorCount) : Environment.ProcessorCount;

    public static ParallelOptions Options(CancellationToken ct = default) => new() { CancellationToken = ct, MaxDegreeOfParallelism = Threads };
}
