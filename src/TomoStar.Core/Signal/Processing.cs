using System.Numerics;
using System.Text.Json.Serialization;
using TomoStar.Core.IO;

namespace TomoStar.Core.Signal;

public enum FilterKind
{
    None,
    BandPass,
    HighPass,
    LowPass
}

/// <summary>
/// Butterworth filters of any kind, causal or zero-phase. A causal filter runs forward only: it has
/// the phase delay of a minimum-phase filter but puts no energy before an onset, which is what an
/// onset picker needs. A zero-phase filter runs forward and backward: no delay, but the impulse
/// response is symmetric, so energy leaks before a sharp onset and a pick made on it can be early.
/// </summary>
public static class Filters
{
    public static void Apply(Span<float> x, double sampleRate, FilterKind kind, double lowHz, double highHz, int order, bool zeroPhase)
    {
        switch (kind)
        {
            case FilterKind.BandPass: Butterworth.Filter(x, sampleRate, lowHz, highHz, order, zeroPhase); break;
            case FilterKind.HighPass: Butterworth.Filter(x, sampleRate, lowHz, 0, order, zeroPhase); break;
            case FilterKind.LowPass: Butterworth.Filter(x, sampleRate, 0, highHz, order, zeroPhase); break;
        }
    }

    /// <summary>
    /// Second-order notch at <paramref name="frequencyHz"/> with quality factor Q (bandwidth f₀/Q),
    /// from the bilinear-transform biquad of R. Bristow-Johnson's "Audio EQ cookbook"; applied
    /// forward and backward so it adds no phase shift. For mains hum (50 or 60 Hz).
    /// </summary>
    public static void Notch(Span<float> x, double sampleRate, double frequencyHz, double q)
    {
        if (frequencyHz <= 0 || frequencyHz >= 0.5 * sampleRate || q <= 0) return;
        var w0 = 2 * Math.PI * frequencyHz / sampleRate;
        var alpha = Math.Sin(w0) / (2 * q);
        var a0 = 1 + alpha;
        var cos = Math.Cos(w0);
        var f = new Butterworth.Section(1 / a0, -2 * cos / a0, 1 / a0, -2 * cos / a0, (1 - alpha) / a0);
        Butterworth.Run(x, f, reverse: false);
        Butterworth.Run(x, f, reverse: true);
    }
}

/// <summary>
/// Change of sampling rate. Going down, the signal is first low-passed (zero-phase Butterworth,
/// order 8, corner at 40% of the new rate, i.e. 80% of the new Nyquist) so nothing aliases; then it is
/// interpolated at the new instants with a Lanczos-windowed sinc kernel of 8 lobes per side.
/// </summary>
public static class Resampler
{
    private const int Lobes = 8;

    public static float[] Resample(ReadOnlySpan<float> x, double fromHz, double toHz)
    {
        if (x.Length == 0 || fromHz <= 0 || toHz <= 0 || Math.Abs(fromHz - toHz) < 1e-9 * fromHz) return x.ToArray();
        var src = x.ToArray();
        if (toHz < fromHz) Butterworth.Filter(src, fromHz, 0, 0.4 * toHz, 8, zeroPhase: true);
        var duration = (src.Length - 1) / fromHz;
        var n = (int)Math.Floor(duration * toHz) + 1;
        var y = new float[n];
        // When decimating the kernel is stretched to the output rate (it is a low-pass there too).
        var scale = Math.Min(1.0, toHz / fromHz);
        var support = Lobes / scale;
        Parallel.For(0, n, j =>
        {
            var t = j * fromHz / toHz; // position in input samples
            var i0 = (int)Math.Ceiling(t - support);
            var i1 = (int)Math.Floor(t + support);
            double sum = 0, wsum = 0;
            for (var i = Math.Max(0, i0); i <= Math.Min(src.Length - 1, i1); i++)
            {
                var d = (t - i) * scale;
                var w = Lanczos(d);
                sum += w * src[i];
                wsum += w;
            }
            y[j] = (float)(wsum != 0 ? sum / wsum : 0);
        });
        return y;
    }

    private static double Lanczos(double d)
    {
        if (Math.Abs(d) < 1e-12) return 1;
        if (Math.Abs(d) >= Lobes) return 0;
        var pd = Math.PI * d;
        return Lobes * Math.Sin(pd) * Math.Sin(pd / Lobes) / (pd * pd);
    }

    public static Trace Resample(Trace t, double toHz)
    {
        var r = t.Clone(Resample(t.Data, t.SampleRate, toHz));
        r.SampleRate = toHz;
        return r;
    }
}

/// <summary>
/// Rotation of horizontal components. Horizontals of any orientation (N/E, or 1/2 with azimuths from
/// the station metadata) are first brought to north and east, then to radial and transverse for a
/// given back azimuth (station to event, degrees from north). Conventions of ObsPy's rotate_ne_rt:
/// R = −N·cos(baz) − E·sin(baz) points away from the source, T = N·sin(baz) − E·cos(baz).
/// </summary>
public static class Rotation
{
    /// <summary>North and east from two horizontals of azimuths a1 and a2 (degrees from north).</summary>
    public static (float[] N, float[] E) ToNorthEast(ReadOnlySpan<float> h1, double azimuth1, ReadOnlySpan<float> h2, double azimuth2)
    {
        var n = Math.Min(h1.Length, h2.Length);
        var north = new float[n];
        var east = new float[n];
        double c1 = Math.Cos(azimuth1 * Math.PI / 180), s1 = Math.Sin(azimuth1 * Math.PI / 180);
        double c2 = Math.Cos(azimuth2 * Math.PI / 180), s2 = Math.Sin(azimuth2 * Math.PI / 180);
        // Solve [c1 c2; s1 s2]·[?] … i.e. h1 = N·c1 + E·s1, h2 = N·c2 + E·s2 for N and E (works for
        // non-orthogonal pairs too, as long as they are not parallel).
        var det = c1 * s2 - s1 * c2;
        if (Math.Abs(det) < 1e-6) throw new InvalidOperationException("The two horizontal components are parallel.");
        for (var i = 0; i < n; i++)
        {
            north[i] = (float)((h1[i] * s2 - h2[i] * s1) / det);
            east[i] = (float)((h2[i] * c1 - h1[i] * c2) / det);
        }
        return (north, east);
    }

    public static (float[] R, float[] T) NorthEastToRadialTransverse(ReadOnlySpan<float> north, ReadOnlySpan<float> east, double backAzimuthDeg)
    {
        var n = Math.Min(north.Length, east.Length);
        var r = new float[n];
        var t = new float[n];
        var b = backAzimuthDeg * Math.PI / 180;
        double cb = Math.Cos(b), sb = Math.Sin(b);
        for (var i = 0; i < n; i++)
        {
            r[i] = (float)(-north[i] * cb - east[i] * sb);
            t[i] = (float)(north[i] * sb - east[i] * cb);
        }
        return (r, t);
    }

    /// <summary>
    /// Radial and transverse traces of a station from its horizontals (at least two, same sample
    /// rate), cut to their common time span. Azimuths come from <paramref name="azimuthOf"/>
    /// (channel code → degrees from north), with N = 0° and E = 90° when the metadata say nothing.
    /// </summary>
    public static (Trace R, Trace T)? RadialTransverse(IReadOnlyList<Trace> horizontals, double backAzimuthDeg, Func<Trace, double?> azimuthOf)
    {
        if (horizontals.Count < 2) return null;
        var a = horizontals[0];
        var b = horizontals.Skip(1).FirstOrDefault(h => Math.Abs(h.SampleRate - a.SampleRate) < 1e-6 * a.SampleRate);
        if (b == null) return null;
        double Az(Trace t) => azimuthOf(t) ?? t.Component switch { 'N' => 0, 'E' => 90, '1' => 0, '2' => 90, _ => double.NaN };
        double azA = Az(a), azB = Az(b);
        if (!double.IsFinite(azA) || !double.IsFinite(azB)) return null;
        var start = a.StartTime > b.StartTime ? a.StartTime : b.StartTime;
        int oa = a.IndexOf(start), ob = b.IndexOf(start);
        var len = Math.Min(a.Data.Length - oa, b.Data.Length - ob);
        if (len < 2) return null;
        var (n, e) = ToNorthEast(a.Data.AsSpan(oa, len), azA, b.Data.AsSpan(ob, len), azB);
        var (r, t) = NorthEastToRadialTransverse(n, e, backAzimuthDeg);
        var rt = a.Clone(r);
        rt.StartTime = a.TimeOf(oa);
        rt.Channel = a.Channel.Length > 0 ? a.Channel[..^1] + "R" : "R";
        var tt = a.Clone(t);
        tt.StartTime = rt.StartTime;
        tt.Channel = a.Channel.Length > 0 ? a.Channel[..^1] + "T" : "T";
        return (rt, tt);
    }
}

/// <summary>Short-time Fourier power of a trace: time × frequency, in dB.</summary>
public sealed record SpectrogramResult(double[] Times, double[] Frequencies, float[,] PowerDb, double MinDb, double MaxDb);

public static class Spectrogram
{
    /// <summary>
    /// Hann-windowed STFT with the given window length and overlap; times are window centres in s
    /// from the first sample. Power is |X|² in dB relative to the maximum.
    /// </summary>
    public static SpectrogramResult Compute(ReadOnlySpan<float> x, double sampleRate, double windowSeconds, double overlap = 0.75)
    {
        var w = Math.Clamp((int)Math.Round(windowSeconds * sampleRate), 16, Math.Max(16, x.Length));
        var nfft = Fft.NextPow2(w);
        var hop = Math.Max(1, (int)Math.Round(w * (1 - Math.Clamp(overlap, 0, 0.95))));
        var frames = x.Length >= w ? 1 + (x.Length - w) / hop : 0;
        var half = nfft / 2;
        var power = new float[frames, half];
        var times = new double[frames];
        var freqs = Enumerable.Range(0, half).Select(k => k * sampleRate / nfft).ToArray();
        var window = Enumerable.Range(0, w).Select(i => 0.5 * (1 - Math.Cos(2 * Math.PI * i / (w - 1)))).ToArray();
        var data = x.ToArray();
        double max = double.NegativeInfinity;
        var maxLock = new object();
        Parallel.For(0, frames, f =>
        {
            var a = new Complex[nfft];
            var o = f * hop;
            double mean = 0;
            for (var i = 0; i < w; i++) mean += data[o + i];
            mean /= w;
            for (var i = 0; i < w; i++) a[i] = new Complex((data[o + i] - mean) * window[i], 0);
            Fft.Transform(a);
            var localMax = double.NegativeInfinity;
            for (var k = 0; k < half; k++)
            {
                var p = a[k].Real * a[k].Real + a[k].Imaginary * a[k].Imaginary;
                var db = 10 * Math.Log10(p + 1e-30);
                power[f, k] = (float)db;
                if (db > localMax) localMax = db;
            }
            times[f] = (o + 0.5 * w) / sampleRate;
            lock (maxLock) max = Math.Max(max, localMax);
        });
        if (!double.IsFinite(max)) max = 0;
        for (var f = 0; f < frames; f++)
        for (var k = 0; k < half; k++) power[f, k] -= (float)max;
        return new SpectrogramResult(times, freqs, power, -80, 0);
    }
}

public enum ProcessingStepKind
{
    Demean,
    Detrend,
    Taper,
    Filter,
    Notch,
    Resample,
    RemoveResponse
}

/// <summary>One step of a processing chain. Only the fields of its kind matter.</summary>
public sealed class ProcessingStep
{
    public ProcessingStepKind Kind { get; set; }
    public bool Enabled { get; set; } = true;
    public FilterKind Filter { get; set; } = FilterKind.BandPass;
    public double LowHz { get; set; } = 1;
    public double HighHz { get; set; } = 20;
    public int Order { get; set; } = 4;
    public bool ZeroPhase { get; set; } = true;
    public double NotchHz { get; set; } = 50;
    public double NotchQ { get; set; } = 30;
    public double TaperFraction { get; set; } = 0.05;
    public double SampleRateHz { get; set; } = 50;

    /// <summary>Response removal: output ground motion, pre-filter corners (0 = automatic) and water level.</summary>
    public GroundMotion Output { get; set; } = GroundMotion.Velocity;

    public double PreF1 { get; set; }
    public double PreF2 { get; set; }
    public double PreF3 { get; set; }
    public double PreF4 { get; set; }
    public double WaterLevelDb { get; set; } = 60;

    [JsonIgnore]
    public string Description => Kind switch
    {
        ProcessingStepKind.Demean => "Remove mean",
        ProcessingStepKind.Detrend => "Remove linear trend",
        ProcessingStepKind.Taper => $"Cosine taper {100 * TaperFraction:0.#}% each end",
        ProcessingStepKind.Filter => Filter switch
        {
            FilterKind.BandPass => $"Band-pass {LowHz:0.###}-{HighHz:0.###} Hz",
            FilterKind.HighPass => $"High-pass {LowHz:0.###} Hz",
            FilterKind.LowPass => $"Low-pass {HighHz:0.###} Hz",
            _ => "No filter"
        } + $", Butterworth order {Order}, {(ZeroPhase ? "zero-phase" : "causal")}",
        ProcessingStepKind.Notch => $"Notch {NotchHz:0.##} Hz (Q {NotchQ:0.#})",
        ProcessingStepKind.Resample => $"Resample to {SampleRateHz:0.###} Hz",
        ProcessingStepKind.RemoveResponse => $"Remove instrument response → {Output.ToString().ToLowerInvariant()} ({ResponseRemoval.Units(Output)}), water level {WaterLevelDb:0} dB" +
                                             (PreF1 > 0 ? $", pre-filter {PreF1:0.###}/{PreF2:0.###}-{PreF3:0.###}/{PreF4:0.###} Hz" : ", automatic pre-filter"),
        _ => Kind.ToString()
    };

    public ProcessingStep Clone() => (ProcessingStep)MemberwiseClone();
}

/// <summary>
/// An ordered list of processing steps, saved with the project and applied on the fly to copies of
/// the waveforms (the viewer, the autopicker, exports). The recorded data are never modified.
/// </summary>
public sealed class ProcessingChain
{
    public List<ProcessingStep> Steps { get; set; } = [];

    public static ProcessingChain Default() => new()
    {
        Steps =
        [
            new ProcessingStep { Kind = ProcessingStepKind.Demean },
            new ProcessingStep { Kind = ProcessingStepKind.Detrend },
            new ProcessingStep { Kind = ProcessingStepKind.Taper, TaperFraction = 0.02 }
        ]
    };

    [JsonIgnore] public bool IsEmpty => !Steps.Any(s => s.Enabled);

    public string Describe() => IsEmpty ? "(no processing)" : string.Join(" → ", Steps.Where(s => s.Enabled).Select(s => s.Description));

    /// <summary>A processed copy of the trace (the sample rate may change).</summary>
    public Trace Apply(Trace trace) => Apply(trace, null);

    /// <summary>
    /// A processed copy of the trace. <paramref name="responseOf"/> gives the instrument response of
    /// a trace for the response-removal step; without one that step is skipped.
    /// </summary>
    public Trace Apply(Trace trace, Func<Trace, InstrumentResponse?>? responseOf)
    {
        var data = (float[])trace.Data.Clone();
        var fs = trace.SampleRate;
        var units = trace.Units;
        foreach (var s in Steps.Where(s => s.Enabled))
        {
            switch (s.Kind)
            {
                case ProcessingStepKind.Demean: SignalProcessing.Demean(data); break;
                case ProcessingStepKind.Detrend: SignalProcessing.Detrend(data); break;
                case ProcessingStepKind.Taper: SignalProcessing.Taper(data, s.TaperFraction); break;
                case ProcessingStepKind.Filter:
                    Filters.Apply(data, fs, s.Filter, s.LowHz, Math.Min(s.HighHz, 0.49 * fs), s.Order, s.ZeroPhase);
                    break;
                case ProcessingStepKind.Notch: Filters.Notch(data, fs, s.NotchHz, s.NotchQ); break;
                case ProcessingStepKind.Resample:
                    if (s.SampleRateHz > 0 && Math.Abs(s.SampleRateHz - fs) > 1e-9 * fs)
                    {
                        data = Resampler.Resample(data, fs, s.SampleRateHz);
                        fs = s.SampleRateHz;
                    }
                    break;
                case ProcessingStepKind.RemoveResponse:
                    if (responseOf?.Invoke(trace) is not { } r) break;
                    var (f1, f2, f3, f4) = s.PreF1 > 0 && s.PreF4 > s.PreF1 ? (s.PreF1, s.PreF2, s.PreF3, s.PreF4) : ResponseRemoval.DefaultPreFilter(fs, data.Length / fs);
                    data = ResponseRemoval.Remove(data, fs, r, s.Output, f1, f2, Math.Min(f3, 0.49 * fs), Math.Min(f4, 0.5 * fs), s.WaterLevelDb);
                    units = ResponseRemoval.Units(s.Output);
                    break;
            }
        }
        var t = trace.Clone(data);
        t.SampleRate = fs;
        t.Units = units;
        return t;
    }
}
