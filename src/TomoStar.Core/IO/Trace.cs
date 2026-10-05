// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

namespace TomoStar.Core.IO;

/// <summary>A continuous, evenly sampled seismogram segment of one channel.</summary>
public sealed class Trace
{
    public string Network { get; set; } = "";
    public string Station { get; set; } = "";
    public string Location { get; set; } = "";
    public string Channel { get; set; } = "";
    public DateTime StartTime { get; set; }
    public double SampleRate { get; set; }
    public float[] Data { get; set; } = [];
    public string Units { get; set; } = "counts";

    public double Delta => 1.0 / SampleRate;
    public DateTime EndTime => StartTime.AddSeconds((Data.Length - 1) / SampleRate);
    public string Nslc => $"{Network}.{Station}.{Location}.{Channel}";
    public string StationId => $"{Network}.{Station}";

    /// <summary>Component letter: the last character of the channel code (Z, N, E, 1, 2, …).</summary>
    public char Component => Channel.Length > 0 ? char.ToUpperInvariant(Channel[^1]) : '?';

    /// <summary>
    /// Vertical by SEED orientation code. Channels coded 1/2/3 are not assumed vertical: their
    /// orientation lives in the station metadata (dip), which callers consult when it matters.
    /// </summary>
    public bool IsVertical => Component == 'Z';

    /// <summary>Second letter of the SEED channel code: H/L/G/N = instrument type.</summary>
    public char InstrumentCode => Channel.Length > 1 ? char.ToUpperInvariant(Channel[1]) : '?';

    /// <summary>True for accelerometers (SEED instrument code N or G).</summary>
    public bool IsAccelerometer => InstrumentCode is 'N' or 'G' || Units.Contains("s^2", StringComparison.Ordinal) || Units.Contains("S**2", StringComparison.OrdinalIgnoreCase);

    /// <summary>Sample index nearest to an absolute time (may be outside the trace).</summary>
    public int IndexOf(DateTime t) => (int)Math.Round((t - StartTime).TotalSeconds * SampleRate);

    public DateTime TimeOf(int index) => StartTime.AddSeconds(index / SampleRate);

    public Trace Clone(float[]? data = null) => new()
    {
        Network = Network, Station = Station, Location = Location, Channel = Channel,
        StartTime = StartTime, SampleRate = SampleRate, Units = Units,
        Data = data ?? (float[])Data.Clone()
    };

    /// <summary>Copy of the samples between two times (clamped to the trace).</summary>
    public Trace Slice(DateTime from, DateTime to)
    {
        var a = Math.Clamp(IndexOf(from), 0, Data.Length);
        var b = Math.Clamp(IndexOf(to), a, Data.Length);
        var t = Clone(Data[a..b]);
        t.StartTime = TimeOf(a);
        return t;
    }
}
