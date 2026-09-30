// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using TomoStar.Core.IO;
using TomoStar.Core.Model;
using TomoStar.Core.Signal;

namespace TomoStar.Core.Attenuation;

/// <summary>How t* is measured.</summary>
public enum TStarMethod
{
    /// <summary>Tapered FFT spectra, ln Ω₀ and t* fitted per record at a common corner (Eberhart-Phillips &amp; Chadwick 2002).</summary>
    SingleTaper,

    /// <summary>
    /// Multitaper spectra, the records of an event inverted jointly with a shared source level,
    /// optional frequency-dependent attenuation, t* ≥ 0 and S t* on the horizontals
    /// (after Stachnik et al. 2004 and Wei &amp; Wiens 2018; see <see cref="MultitaperTStar"/>).
    /// </summary>
    MultitaperJoint
}

public sealed class TStarSettings
{
    public TStarMethod Method { get; set; } = TStarMethod.SingleTaper;

    /// <summary>Length of the P window after the pick, s (shortened so it ends before S).</summary>
    public double WindowSeconds { get; set; } = 2.56;

    /// <summary>
    /// The window starts this long before the P pick. It must hold the whole onset even when the pick
    /// is late: with dispersion the high frequencies of an attenuated pulse arrive before the pick
    /// made on the dominant frequency, and a window that cuts the onset puts a step into the signal
    /// whose flat spectrum biases t* low (by half for t* ≈ 0.09 s with a 0.1 s margin on the
    /// synthetic benchmark, exact with 0.3 s even for picks 0.1 s late).
    /// </summary>
    public double PreSeconds { get; set; } = 0.3;
    public double MinFrequencyHz { get; set; } = 1;
    public double MaxFrequencyHz { get; set; } = 30;
    public double MinCornerHz { get; set; } = 0.3;
    public double MaxCornerHz { get; set; } = 40;

    /// <summary>Minimum spectral signal-to-noise ratio for a frequency to enter the fit.</summary>
    public double MinSnr { get; set; } = 3;

    /// <summary>Minimum usable bandwidth as a ratio fmax/fmin.</summary>
    public double MinBandRatio { get; set; } = 3;

    public int MinStationsPerEvent { get; set; } = 3;

    /// <summary>
    /// Stress-drop range (MPa) that bounds the corner-frequency search when the event has a
    /// magnitude: f_c = 0.37 β (16 Δσ / 7 M₀)^(1/3) (Brune 1970) with β = 3.5 km/s and
    /// M₀ = 10^(1.5 M + 9.1) N·m (Hanks &amp; Kanamori 1979). Without it the search can give a large
    /// earthquake a corner above the fit band (3 Hz for the 2016 Amatrice Mw 6.0), which the t* of
    /// every record then absorbs. Set the low end to 0 to disable.
    /// </summary>
    public double MinStressDropMPa { get; set; } = 0.1;

    public double MaxStressDropMPa { get; set; } = 100;

    /// <summary>
    /// An event's corner frequency counts as resolved when it lies at least this factor below the
    /// highest fitted frequency of its records (median): corners within a factor of three of the
    /// maximum signal frequency are underestimated (Abercrombie 2015; Abercrombie et al. 2017, GJI 208,
    /// 306-320, for spectral ratios, where attenuation is already removed; fitting t* as well only
    /// widens the trade-off, Ko et al. 2012).
    /// </summary>
    public double CornerBandFactor { get; set; } = 3;

    /// <summary>Flag (disable) the t* of events whose corner frequency is not resolved.</summary>
    public bool ExcludeUnresolvedCorners { get; set; } = true;

    /// <summary>Events below this magnitude are not measured (NaN or ≤ 0: no limit; events without magnitude are kept).</summary>
    public double MinMagnitude { get; set; } = double.NaN;

    /// <summary>Remove the noise power (from the pre-event window) from the signal spectrum before fitting.</summary>
    public bool SubtractNoise { get; set; } = true;

    // ---- Multitaper joint method ----

    /// <summary>
    /// Time-bandwidth product NW of the Slepian tapers (2NW − 1 tapers); 0 = automatic, NW = f_min × window
    /// length within 1.5-4. The spectrum is averaged over ±W = NW / T: when W exceeds the lowest fitted
    /// frequency the low end of the spectrum is flattened and mimics a higher corner, so small events get
    /// their corner at the top of the search range and t* tens of milliseconds too high (NW = 4 on a
    /// 2.56 s window does this to 12 of the 70 benchmark events; NW = 1.5-3 to none).
    /// </summary>
    public double TimeBandwidth { get; set; }

    /// <summary>The time-bandwidth product used for a window of <paramref name="lengthSeconds"/>.</summary>
    public double TimeBandwidthFor(double lengthSeconds) =>
        TimeBandwidth > 0 ? TimeBandwidth : Math.Clamp(MinFrequencyHz * lengthSeconds, 1.5, 4);

    /// <summary>Frequency dependence of attenuation, Q ∝ f^α (t*(f) = t*₁ f^−α); 0 = constant Q.</summary>
    public double Alpha { get; set; }

    /// <summary>Search α (0 to <see cref="MaxAlpha"/>) jointly with the corner frequency.</summary>
    public bool SearchAlpha { get; set; }

    public double MaxAlpha { get; set; } = 0.6;

    /// <summary>One source level per event (joint inversion) rather than one per record.</summary>
    public bool SharedSourceLevel { get; set; } = true;

    /// <summary>
    /// Keep t* ≥ 0 (active-set least squares, as the published method). Off by default: with a shared
    /// source level, a record held at zero moves that level and with it every other t* of the event
    /// (+5 ms on the test event whose nearest path has t* = 0); clearly negative values are flagged
    /// instead, as in the single-taper method.
    /// </summary>
    public bool NonNegative { get; set; }

    /// <summary>
    /// Fit only the longest contiguous run of frequencies above the SNR threshold. Off (default): all
    /// of them, as the single-taper method; on the 2016 Amatrice records the contiguous band gave t*
    /// half as consistent (18.5 against 12.6 ms scatter), on the synthetic benchmark the same result.
    /// </summary>
    public bool ContiguousBand { get; set; }

    /// <summary>ln A must decrease with f over the band with at least this correlation.</summary>
    public double MinSpectralCorrelation { get; set; } = 0.5;

    /// <summary>Records used to find the corner frequency: SNR at least this many times the threshold.</summary>
    public double HighQualitySnrFactor { get; set; } = 2;

    /// <summary>Records whose misfit exceeds this multiple of the event's median are flagged.</summary>
    public double MaxRelativeMisfit { get; set; } = 3;

    /// <summary>Measure S t* on the horizontal components too.</summary>
    public bool MeasureS { get; set; }

    /// <summary>
    /// Search the S corner frequency on the S spectra of the event (default), between f_c(P)/2 and
    /// f_c(P). Off: f_c(S) = f_c(P) / <see cref="CornerRatioPToS"/>.
    /// The ratio is not a constant of nature (1.5 for Madariaga's 1976 model; about 1-1.3 observed
    /// for small earthquakes), and when it is wrong t*_S absorbs the difference: on a synthetic with a
    /// true ratio of 1.1, the fixed 1.5 gave t*_S/t*_P = 1.76 for a true 2.0 (std 9 ms), the search 1.98
    /// (std 2 ms).
    /// </summary>
    public bool SearchCornerS { get; set; } = true;

    /// <summary>f_c(P) / f_c(S) when the S corner is not searched (Madariaga 1976: about 1.5).</summary>
    public double CornerRatioPToS { get; set; } = 1.5;

    /// <summary>
    /// The S noise window lies between the P onset and the S pick (the P coda, which is what overlaps
    /// the S wave), not before the P wave: frequencies where the P coda is comparable to S are then left
    /// out instead of flattening the S spectrum with less attenuated P energy.
    /// </summary>
    public bool SNoiseFromPCoda { get; set; } = true;

    /// <summary>Length of the P-coda noise window before S (s).</summary>
    public double SNoiseSeconds { get; set; } = 1.0;
}

/// <summary>A measured spectrum with its fit, kept for the diagnostics plot.</summary>
public sealed record SpectrumFit(
    string EventId, string StationId, double[] Freq, double[] Signal, double[] Noise, bool[] Used,
    double TStar, double TStarError, double CornerHz, double LogOmega0, double Snr, double Alpha = 0);

/// <summary>
/// Whole-path P attenuation t* from displacement spectra.
///
/// The displacement amplitude spectrum of a P wave is modelled as a Brune (1970, JGR 75(26),
/// 4997-5009) ω² source times a frequency-independent-Q path term,
///
///     A(f) = Ω₀ · exp(−π f t*) / (1 + (f/f_c)²),
///
/// so ln A = ln Ω₀ − π f t* − ln(1 + (f/f_c)²) is linear in (ln Ω₀, t*) for a given corner
/// frequency f_c. The corner frequency is a property of the source, so it is searched once per event
/// over all its stations (grid search in log f_c), and t* and Ω₀ are then solved per record by least
/// squares; this is the procedure of Eberhart-Phillips &amp; Chadwick (2002, JGR 107(B2), 2033).
/// Only frequencies where the signal exceeds the pre-event noise by <see cref="TStarSettings.MinSnr"/>
/// enter the fit, and the noise power is subtracted from the signal power in each of them. Spectra
/// are smoothed in log-spaced bins (1/6 octave) before fitting.
///
/// When the channel's response is known (StationXML at response level), the whole record is first
/// deconvolved to ground displacement (spectral division with a water level and a pre-filter outside
/// the fit band) and the signal and noise windows are cut from it. Without it the response is
/// assumed flat over the fit band (broadband velocity or accelerometer channels), and velocity or
/// acceleration is converted to displacement by dividing by (2πf) or (2πf)².
/// </summary>
public static class TStarEstimator
{
    public static (List<TStarMeasurement> Results, List<SpectrumFit> Fits, string Message) MeasureEvent(
        EventRecord ev, IReadOnlyDictionary<string, List<Trace>> tracesByStation, TStarSettings s,
        Func<Trace, InstrumentResponse?>? responseOf = null, Func<string, StationRecord?>? stationOf = null)
    {
        if (s.MinMagnitude > 0 && double.IsFinite(ev.Magnitude) && ev.Magnitude < s.MinMagnitude)
            return ([], [], $"{ev.Id}: M {ev.Magnitude:0.0} below the minimum {s.MinMagnitude:0.0}");
        if (s.Method == TStarMethod.MultitaperJoint) return MultitaperTStar.MeasureEvent(ev, tracesByStation, s, responseOf, stationOf);
        var prepared = new List<(string Station, double[] F, double[] A, bool[] Used, double[] Noise, double Snr)>();
        var chosen = PickSelection.Best(ev.Picks).ToList();
        foreach (var p in chosen.Where(p => p.Phase == Phase.P))
        {
            if (!tracesByStation.TryGetValue(p.StationId, out var traces)) continue;
            var tr = traces.FirstOrDefault(t => t.IsVertical) ?? traces.FirstOrDefault();
            if (tr == null) continue;
            var sPick = chosen.FirstOrDefault(q => q.Phase == Phase.S && q.StationId == p.StationId);
            var spec = Spectrum(tr, p.Time, sPick?.Time, s, responseOf?.Invoke(tr));
            if (spec != null) prepared.Add((p.StationId, spec.Value.F, spec.Value.A, spec.Value.Used, spec.Value.Noise, spec.Value.Snr));
        }
        if (prepared.Count < s.MinStationsPerEvent) return ([], [], $"{ev.Id}: {prepared.Count} usable spectra (minimum {s.MinStationsPerEvent})");

        // Corner-frequency grid search on the summed misfit of all stations of the event, within
        // the range the magnitude allows.
        var (fcLo, fcHi) = CornerRange(ev.Magnitude, s);
        var best = (Fc: double.NaN, Misfit: double.PositiveInfinity);
        var grid = new double[61];
        var rss = new double[61];
        var bestIndex = -1;
        var nData = 0;
        for (var q = 0; q <= 60; q++)
        {
            var fc = fcLo * Math.Pow(fcHi / fcLo, q / 60.0);
            grid[q] = fc;
            double total = 0;
            var n = 0;
            foreach (var x in prepared)
            {
                var fit = Fit(x.F, x.A, x.Used, fc);
                if (fit == null) continue;
                total += fit.Value.Rss;
                n += fit.Value.N;
            }
            rss[q] = n == 0 ? double.PositiveInfinity : total;
            if (n == 0) continue;
            var m = total / n;
            if (m < best.Misfit) { best = (fc, m); bestIndex = q; nData = n; }
        }
        if (double.IsNaN(best.Fc)) return ([], [], $"{ev.Id}: no corner frequency fits");
        var (iLo, iHi) = CornerTradeOff.Interval(rss, bestIndex, nData, 2 * prepared.Count + 1);

        var results = new List<TStarMeasurement>();
        var fits = new List<SpectrumFit>();
        var spread = new List<double>();
        foreach (var x in prepared)
        {
            var fit = Fit(x.F, x.A, x.Used, best.Fc);
            if (fit == null) continue;
            var (lnO, tstar, err, _, _) = fit.Value;
            fits.Add(new SpectrumFit(ev.Id, x.Station, x.F, x.A, x.Noise, x.Used, tstar, err, best.Fc, lnO, x.Snr));
            // A clearly negative t* is unphysical (site amplification or a bad corner); keep it flagged.
            results.Add(new TStarMeasurement
            {
                EventId = ev.Id, StationId = x.Station, Phase = Phase.P, TStar = tstar, Uncertainty = Math.Max(0.002, err),
                CornerFrequencyHz = best.Fc, LogOmega0 = lnO, FitMinHz = x.F.Where((_, i) => x.Used[i]).Min(),
                FitMaxHz = x.F.Where((_, i) => x.Used[i]).Max(), Snr = x.Snr, Disabled = tstar < -2 * err
            });
            spread.Add(Math.Max(Math.Abs((Fit(x.F, x.A, x.Used, grid[iLo])?.TStar ?? tstar) - tstar),
                                Math.Abs((Fit(x.F, x.A, x.Used, grid[iHi])?.TStar ?? tstar) - tstar)));
        }
        CornerTradeOff.Apply(results, spread, grid[iLo], grid[iHi], best.Fc, s);
        return (results, fits, $"{ev.Id}: fc {best.Fc:0.00} Hz, {results.Count} t*");
    }

    /// <summary>The corner-frequency search range: the settings, narrowed by the magnitude when known.</summary>
    public static (double Lo, double Hi) CornerRange(double magnitude, TStarSettings s)
    {
        double lo = s.MinCornerHz, hi = s.MaxCornerHz;
        if (double.IsFinite(magnitude) && s.MinStressDropMPa > 0 && s.MaxStressDropMPa > s.MinStressDropMPa)
        {
            var m0 = Math.Pow(10, 1.5 * magnitude + 9.1);
            double Fc(double mpa) => 0.37 * 3500 * Math.Cbrt(16 * mpa * 1e6 / (7 * m0));
            var a = Math.Max(lo, Fc(s.MinStressDropMPa));
            var b = Math.Min(hi, Fc(s.MaxStressDropMPa));
            // A magnitude inconsistent with the settings' range leaves the settings alone.
            if (b > a * 1.05) (lo, hi) = (a, b);
            else if (Fc(s.MaxStressDropMPa) < lo) hi = lo * 1.5;
        }
        return (lo, hi);
    }

    /// <summary>
    /// The signal window (from <see cref="TStarSettings.PreSeconds"/> before the pick, at most
    /// <see cref="TStarSettings.WindowSeconds"/> long and ending before <paramref name="endBefore"/>)
    /// and a noise window of the same length before <paramref name="noiseBefore"/>, demeaned, in
    /// displacement when the response is known (power 0) or else in the recorded quantity, whose
    /// spectrum must be divided by (2πf)^power. With <paramref name="noiseNotBefore"/> the noise window
    /// may not start earlier (the S noise is the P coda between the P onset and the S pick); it is then
    /// shorter than the signal when the two picks are close, down to <see cref="MinNoiseSeconds"/>,
    /// and the caller scales its spectrum by √(signal length / noise length). Below that the noise is
    /// taken before <paramref name="noiseNotBefore"/> instead.
    /// </summary>
    /// <summary>Shortest noise window between two picks (s).</summary>
    internal const double MinNoiseSeconds = 0.64;

    internal static (float[] Signal, float[] Noise, int Power, double Length)? Windows(Trace tr, DateTime pick, DateTime? endBefore, TStarSettings s,
        InstrumentResponse? response, DateTime? noiseBefore = null, DateTime? noiseNotBefore = null)
    {
        var fs = tr.SampleRate;
        var len = s.WindowSeconds;
        if (endBefore is { } st) len = Math.Min(len, (st - pick).TotalSeconds - 0.05);
        if (len < 0.5) return null;
        var i0 = tr.IndexOf(pick.AddSeconds(-s.PreSeconds));
        var n = (int)(len * fs);
        var noiseEnd = noiseBefore is { } nb ? tr.IndexOf(nb.AddSeconds(-s.PreSeconds)) : i0;
        var j0 = noiseEnd - n - (int)(0.2 * fs);
        var nNoise = n;
        if (noiseNotBefore is { } nnb)
        {
            var first = tr.IndexOf(nnb);
            if (j0 < first)
            {
                nNoise = noiseEnd - (int)(0.2 * fs) - first;
                j0 = first;
                if (nNoise < MinNoiseSeconds * fs)
                {
                    // Picks too close for a noise window between them: noise before the earlier one.
                    nNoise = n;
                    j0 = tr.IndexOf(nnb.AddSeconds(-s.PreSeconds)) - n - (int)(0.2 * fs);
                }
            }
        }
        if (j0 < 0 || i0 < 0 || i0 + n > tr.Data.Length) return null;

        // With a known response, the whole record is converted to ground displacement first and the
        // windows are cut from it. Dividing the spectrum of a short window by the response instead
        // would amplify the window's spectral leakage wherever the instrument is insensitive (below
        // a geophone's natural frequency, for example), inflating the low-frequency amplitudes.
        float[] record;
        int power;
        if (response != null)
        {
            var nyq = 0.5 * fs;
            var f4 = Math.Min(0.95 * nyq, Math.Max(1.25 * s.MaxFrequencyHz, 0.8 * nyq));
            record = ResponseRemoval.Remove(tr.Data, fs, response, GroundMotion.Displacement,
                0.25 * s.MinFrequencyHz, 0.5 * s.MinFrequencyHz, Math.Min(0.9 * f4, s.MaxFrequencyHz * 1.1), f4, 80);
            power = 0;
        }
        else
        {
            record = tr.Data;
            power = tr.IsAccelerometer ? 2 : 1; // velocity or acceleration assumed flat over the band
        }
        var sig = record.AsSpan(i0, n).ToArray();
        var noi = record.AsSpan(j0, nNoise).ToArray();
        SignalProcessing.Demean(sig);
        SignalProcessing.Demean(noi);
        return (sig, noi, power, len);
    }

    private static (double[] F, double[] A, bool[] Used, double[] Noise, double Snr)? Spectrum(Trace tr, DateTime pTime, DateTime? sTime, TStarSettings s,
        InstrumentResponse? response)
    {
        if (Windows(tr, pTime, sTime, s, response) is not { } cut) return null;
        var (sig, noi, power, len) = cut;
        var fs = tr.SampleRate;
        // Cosine ramps of 0.05 s at the window ends: the window starts at the onset, so a taper
        // proportional to its length (let alone a Hann window) would suppress the onset itself.
        var taper = Math.Min(1, 2 * 0.05 / len);
        var (f, a) = Fft.AmplitudeSpectrum(sig, fs, 256, taper);
        var (_, b) = Fft.AmplitudeSpectrum(noi, fs, 256, taper);
        for (var k = 1; k < f.Length; k++)
        {
            var w = Math.Pow(2 * Math.PI * f[k], power);
            a[k] /= w;
            b[k] /= w;
        }
        // 1/6-octave log bins between fmin and min(fmax, 0.8 Nyquist).
        var fmax = Math.Min(s.MaxFrequencyHz, 0.4 * fs);
        var bins = new List<(double F, double A, double N)>();
        for (var fc = s.MinFrequencyHz; fc <= fmax; fc *= Math.Pow(2, 1.0 / 6))
        {
            var lo = fc / Math.Pow(2, 1.0 / 12);
            var hi = fc * Math.Pow(2, 1.0 / 12);
            double sa = 0, sn = 0;
            var c = 0;
            for (var k = 1; k < f.Length; k++)
            {
                if (f[k] < lo || f[k] >= hi) continue;
                sa += a[k] * a[k];
                sn += b[k] * b[k];
                c++;
            }
            if (c == 0) continue;
            bins.Add((fc, Math.Sqrt(sa / c), Math.Sqrt(sn / c)));
        }
        if (bins.Count < 6) return null;
        var used = bins.Select(x => x.N > 0 && x.A / x.N >= s.MinSnr).ToArray();
        if (s.SubtractNoise)
        {
            // Signal and noise add in power: |S + N|² ≈ |S|² + |N|² on average. Left in, the noise
            // flattens the high-frequency end of the spectrum, where the SNR approaches the
            // threshold, and biases t* (and the corner) low. The pre-event window gives the noise
            // power; it is removed from each bin before the fit (bins used have SNR ≥ MinSnr, so
            // at least 1 − 1/MinSnr² of the power remains).
            for (var i = 0; i < bins.Count; i++)
            {
                var (bf, ba, bn) = bins[i];
                bins[i] = (bf, Math.Sqrt(Math.Max(ba * ba - bn * bn, 0.01 * ba * ba)), bn);
            }
        }
        var usedF = bins.Where((_, i) => used[i]).Select(x => x.F).ToArray();
        if (usedF.Length < 5 || usedF.Max() / usedF.Min() < s.MinBandRatio) return null;
        var snr = bins.Where((_, i) => used[i]).Average(x => x.A / Math.Max(1e-30, x.N));
        return (bins.Select(x => x.F).ToArray(), bins.Select(x => x.A).ToArray(), used, bins.Select(x => x.N).ToArray(), snr);
    }

    /// <summary>Least squares for (ln Ω₀, t*) at a fixed corner frequency.</summary>
    private static (double LnOmega, double TStar, double Err, double Rss, int N)? Fit(double[] f, double[] a, bool[] used, double fc)
    {
        double s1 = 0, sx = 0, sxx = 0, sy = 0, sxy = 0;
        var n = 0;
        var ys = new List<(double X, double Y)>();
        for (var k = 0; k < f.Length; k++)
        {
            if (!used[k] || !(a[k] > 0)) continue;
            var x = -Math.PI * f[k];
            var y = Math.Log(a[k]) + Math.Log(1 + f[k] * f[k] / (fc * fc));
            s1++; sx += x; sxx += x * x; sy += y; sxy += x * y;
            ys.Add((x, y));
            n++;
        }
        if (n < 4) return null;
        var det = s1 * sxx - sx * sx;
        if (Math.Abs(det) < 1e-30) return null;
        var tstar = (s1 * sxy - sx * sy) / det;
        var lnO = (sy - tstar * sx) / s1;
        var rss = ys.Sum(p => Math.Pow(p.Y - lnO - tstar * p.X, 2));
        var sigma2 = rss / Math.Max(1, n - 2);
        var err = Math.Sqrt(sigma2 * s1 / det);
        return (lnO, tstar, err, rss, n);
    }
}
