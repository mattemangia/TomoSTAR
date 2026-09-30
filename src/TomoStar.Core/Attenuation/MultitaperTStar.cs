// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using TomoStar.Core.Geo;
using TomoStar.Core.IO;
using TomoStar.Core.Model;
using TomoStar.Core.Signal;

namespace TomoStar.Core.Attenuation;

/// <summary>
/// t* from multitaper spectra inverted jointly for all the records of an event: the approach of
/// Stachnik et al. (2004, JGR 109, B10304) as developed by Wei &amp; Wiens (2018, EPSL 502, 187-199),
/// written for QUIVER and TomoSTAR from the published method.
///
/// For every record of the event the displacement spectrum, corrected for geometrical spreading
/// (1/R with R the hypocentral distance), is modelled as
///
///     ln A_i(f) + ln R_i = ln Ω₀ − ln(1 + (f/f_c)²) − π f^(1−α) t*_i,
///
/// that is a Brune ω² source common to the event, its seismic moment term ln Ω₀ shared by every
/// station (the joint inversion that constrains the trade-off between source level and t*), and a
/// path term with t*(f) = t*_i f^(−α) (frequency-dependent Q ∝ f^α; α = 0 gives constant Q). The
/// source level can instead be left free per record, which absorbs site amplification and radiation
/// pattern at the cost of a weaker constraint. For a given f_c (and α) the problem is linear in
/// (ln Ω₀, t*_i); t* can be kept non-negative with an active-set least squares. f_c (and α when asked) is
/// found by a grid search on the summed misfit of the better records only (high SNR, wide band), then
/// every record is inverted with it.
///
/// What this implementation adds to the published procedure, from the synthetic t* benchmark it was validated on:
/// the window starts before the pick so a late pick does not cut the onset, the noise power measured
/// before the P wave is removed from each frequency, and the corner-frequency search is bounded by
/// the event's magnitude (Brune stress drop 0.1-100 MPa), without which the search runs to the edge
/// of its range for part of the events and moves their t* by tens of milliseconds.
///
/// Spectra are multitaper estimates (Thomson 1982) with 2NW − 1 Slepian tapers, averaged in 1/6-octave
/// bins. NW is kept small enough that the averaging band ±NW/T stays below the lowest fitted frequency:
/// a wider band flattens the low end of the spectrum, which the corner search reads as a higher corner
/// (the benchmark's small events, and the likely origin of the corners the published implementation
/// placed at the top of its range there). A record needs a continuous band where the signal exceeds
/// the noise by the SNR threshold, at least MinBandRatio wide, and every bin above the threshold is
/// fitted (only that band with ContiguousBand); records whose spectrum does not decay with frequency
/// (|correlation of ln A with f| below a threshold, or rising) are left out. S t* is measured on the
/// horizontal components (their powers summed), with its own corner frequency searched on the S
/// spectra (or f_c(P) / 1.5, Madariaga 1976, BSSA 66(3), 639-666) and its noise taken from the P coda
/// just before the S pick.
/// </summary>
public static class MultitaperTStar
{
    private sealed record Record(string Station, Phase Phase, double[] F, double[] A, double[] Noise, bool[] Used, double Snr, double LnR, double BandRatio);

    public static (List<TStarMeasurement> Results, List<SpectrumFit> Fits, string Message) MeasureEvent(
        EventRecord ev, IReadOnlyDictionary<string, List<Trace>> tracesByStation, TStarSettings s,
        Func<Trace, InstrumentResponse?>? responseOf, Func<string, StationRecord?>? stationOf)
    {
        var chosen = PickSelection.Best(ev.Picks).ToList();
        var h = ev;
        var records = new List<Record>();
        var skipped = 0;
        foreach (var p in chosen.Where(p => p.Phase == Phase.P || (s.MeasureS && p.Phase == Phase.S)))
        {
            if (!tracesByStation.TryGetValue(p.StationId, out var traces) || traces.Count == 0) continue;
            var station = stationOf?.Invoke(p.StationId);
            var r = station != null
                ? Math.Max(1, (GeoMath.ToCartesian(h.Lon, h.Lat, h.DepthKm) - GeoMath.ToCartesian(station.Lon, station.Lat, -station.ElevationM / 1000)).Length)
                : double.NaN;
            var pPick = chosen.FirstOrDefault(q => q.Phase == Phase.P && q.StationId == p.StationId);
            var sPick = chosen.FirstOrDefault(q => q.Phase == Phase.S && q.StationId == p.StationId);
            var rec = p.Phase == Phase.P
                ? SpectrumOf([traces.FirstOrDefault(t => t.IsVertical) ?? traces[0]], p.Time, sPick?.Time, null, s, responseOf)
                : SpectrumOf(traces.Where(t => !t.IsVertical).Take(2).ToList(), p.Time, null, s.SNoiseFromPCoda ? p.Time : pPick?.Time,
                    s, responseOf, s.SNoiseFromPCoda ? SNoiseStart(pPick?.Time, p.Time, s) : null);
            if (rec == null) { skipped++; continue; }
            var (f, a, noise, used, snr, ratio) = rec.Value;
            records.Add(new Record(p.StationId, p.Phase, f, a, noise, used, snr, double.IsFinite(r) ? Math.Log(r) : 0, ratio));
        }
        var results = new List<TStarMeasurement>();
        var fits = new List<SpectrumFit>();
        var pRecords = records.Where(x => x.Phase == Phase.P).ToList();
        if (pRecords.Count < s.MinStationsPerEvent)
            return ([], [], $"{ev.Id}: {pRecords.Count} usable P spectra (minimum {s.MinStationsPerEvent}){(skipped > 0 ? $", {skipped} rejected" : "")}");

        // Corner frequency (and α) from the better records: high SNR and a wide band.
        var good = pRecords.Where(x => x.Snr >= s.HighQualitySnrFactor * s.MinSnr && x.BandRatio >= 1.5 * s.MinBandRatio).ToList();
        if (good.Count < s.MinStationsPerEvent) good = pRecords;
        var (fcLo, fcHi) = TStarEstimator.CornerRange(ev.Magnitude, s);
        var alphas = s.SearchAlpha
            ? Enumerable.Range(0, (int)Math.Round(s.MaxAlpha / 0.05) + 1).Select(i => i * 0.05).ToArray()
            : [s.Alpha];
        var best = (Fc: double.NaN, Alpha: s.Alpha, Misfit: double.PositiveInfinity);
        var gridP = Enumerable.Range(0, 61).Select(q => fcLo * Math.Pow(fcHi / fcLo, q / 60.0)).ToArray();
        double[] rssP = [];
        int bestIndex = -1, nData = 0, nPar = 0;
        foreach (var alpha in alphas)
        {
            var rss = new double[61];
            var improved = false;
            for (var q = 0; q <= 60; q++)
            {
                var sol = Solve(good, gridP[q], alpha, s);
                rss[q] = sol == null ? double.PositiveInfinity : sol.Misfit * sol.N;
                if (sol != null && sol.Misfit < best.Misfit) { best = (gridP[q], alpha, sol.Misfit); bestIndex = q; nData = sol.N; nPar = sol.Parameters; improved = true; }
            }
            if (improved) rssP = rss;
        }
        if (double.IsNaN(best.Fc)) return ([], [], $"{ev.Id}: no corner frequency fits");
        var (iLoP, iHiP) = CornerTradeOff.Interval(rssP, bestIndex, nData, nPar + 1);

        var pResolved = true;
        foreach (var phase in s.MeasureS ? new[] { Phase.P, Phase.S } : [Phase.P])
        {
            var set = records.Where(x => x.Phase == phase).ToList();
            if (set.Count == 0 || (phase == Phase.S && set.Count < 2)) continue;
            var fc = phase == Phase.P ? best.Fc : best.Fc / s.CornerRatioPToS;
            // Corner interval of this phase (for S with a fixed ratio: the P interval scaled).
            double cLo = phase == Phase.P ? gridP[iLoP] : gridP[iLoP] / s.CornerRatioPToS;
            double cHi = phase == Phase.P ? gridP[iHiP] : gridP[iHiP] / s.CornerRatioPToS;
            if (phase == Phase.S && s.SearchCornerS)
            {
                // The S corner from the S spectra, on the better records as for P, between f_c(P)/2
                // and f_c(P) (P/S corner ratios of 1-2). Left free over the magnitude's range, it traded
                // with t* on records whose band the P coda cuts short and ran to the range's end for
                // whole events (t*_S +15 to +20 ms on the synthetic).
                var goodS = set.Where(x => x.Snr >= s.HighQualitySnrFactor * s.MinSnr && x.BandRatio >= 1.5 * s.MinBandRatio).ToList();
                if (goodS.Count < Math.Min(s.MinStationsPerEvent, set.Count)) goodS = set;
                var lo = best.Fc / 2;
                var bestS = double.PositiveInfinity;
                var gridS = Enumerable.Range(0, 31).Select(q => lo * Math.Pow(2, q / 30.0)).ToArray();
                var rssS = new double[31];
                int iS = -1, nS = 0, pS = 0;
                for (var q = 0; q <= 30; q++)
                {
                    var trial = Solve(goodS, gridS[q], best.Alpha, s);
                    rssS[q] = trial == null ? double.PositiveInfinity : trial.Misfit * trial.N;
                    if (trial != null && trial.Misfit < bestS) { bestS = trial.Misfit; fc = gridS[q]; iS = q; nS = trial.N; pS = trial.Parameters; }
                }
                if (iS >= 0)
                {
                    var (a, b) = CornerTradeOff.Interval(rssS, iS, nS, pS + 1);
                    (cLo, cHi) = (gridS[a], gridS[b]);
                }
            }
            var sol = Solve(set, fc, best.Alpha, s);
            if (sol == null) continue;
            var solLo = Solve(set, cLo, best.Alpha, s);
            var solHi = Solve(set, cHi, best.Alpha, s);
            var spread = new List<double>();
            var phaseResults = new List<TStarMeasurement>();
            // Records the common model fits much worse than the rest are flagged, not used.
            var misfits = sol.RecordMisfit.OrderBy(x => x).ToArray();
            var median = misfits[misfits.Length / 2];
            for (var i = 0; i < set.Count; i++)
            {
                var x = set[i];
                var lnO = sol.LnOmega[i] - x.LnR;
                fits.Add(new SpectrumFit(ev.Id, phase == Phase.P ? x.Station : x.Station + " (S)", x.F, x.A, x.Noise, x.Used,
                    sol.TStar[i], sol.Error[i], fc, lnO, x.Snr, best.Alpha));
                spread.Add(Math.Max(solLo == null ? 0 : Math.Abs(solLo.TStar[i] - sol.TStar[i]), solHi == null ? 0 : Math.Abs(solHi.TStar[i] - sol.TStar[i])));
                phaseResults.Add(new TStarMeasurement
                {
                    EventId = ev.Id, StationId = x.Station, Phase = phase, TStar = sol.TStar[i], Uncertainty = Math.Max(0.002, sol.Error[i]),
                    CornerFrequencyHz = fc, LogOmega0 = lnO, Alpha = best.Alpha,
                    FitMinHz = x.F.Where((_, k) => x.Used[k]).Min(), FitMaxHz = x.F.Where((_, k) => x.Used[k]).Max(), Snr = x.Snr,
                    // Flagged: a fit far worse than the event's other records, or a clearly negative t*.
                    Disabled = sol.RecordMisfit[i] > s.MaxRelativeMisfit * Math.Max(median, 1e-6) || sol.TStar[i] < -2 * sol.Error[i]
                });
            }
            // The S corner's range is set by the P corner: an unresolved P corner leaves S unresolved.
            CornerTradeOff.Apply(phaseResults, spread, cLo, cHi, fc, s, forceUnresolved: phase == Phase.S && !pResolved);
            if (phase == Phase.P) pResolved = phaseResults.Count == 0 || phaseResults[0].CornerResolved;
            results.AddRange(phaseResults);
        }
        return (results, fits, $"{ev.Id}: fc {best.Fc:0.00} Hz{(s.SearchAlpha || best.Alpha != 0 ? $", α {best.Alpha:0.00}" : "")}, " +
                               $"{results.Count(m => m.Phase == Phase.P)} P{(s.MeasureS ? $" and {results.Count(m => m.Phase == Phase.S)} S" : "")} t*");
    }

    /// <summary>
    /// Multitaper spectrum of one or more components (powers summed), converted to displacement, in
    /// 1/6-octave bins, with the noise removed and the bins above the SNR threshold marked.
    /// </summary>
    private static (double[] F, double[] A, double[] Noise, bool[] Used, double Snr, double BandRatio)? SpectrumOf(
        List<Trace> traces, DateTime pick, DateTime? endBefore, DateTime? noiseBefore, TStarSettings s, Func<Trace, InstrumentResponse?>? responseOf,
        DateTime? noiseNotBefore = null)
    {
        if (traces.Count == 0) return null;
        double[]? freq = null;
        double[]? sigPow = null, noiPow = null;
        var fs = traces[0].SampleRate;
        foreach (var tr in traces)
        {
            if (Math.Abs(tr.SampleRate - fs) > 1e-6) continue;
            if (TStarEstimator.Windows(tr, pick, endBefore, s, responseOf?.Invoke(tr), noiseBefore, noiseNotBefore) is not { } w) continue;
            var nw = s.TimeBandwidthFor(w.Length);
            var nfft = Math.Max(256, w.Signal.Length);
            var (f, a) = Multitaper.AmplitudeSpectrum(w.Signal, fs, nw, minLength: nfft);
            // A shorter noise window (between two picks) is brought to the signal's length: the
            // amplitude of stationary noise grows as the square root of the duration.
            var (_, b) = Multitaper.AmplitudeSpectrum(w.Noise, fs, s.TimeBandwidthFor(w.Noise.Length / fs), minLength: nfft);
            var noiseScale = Math.Sqrt((double)w.Signal.Length / w.Noise.Length);
            for (var k = 0; k < b.Length; k++) b[k] *= noiseScale;
            freq ??= f;
            sigPow ??= new double[f.Length];
            noiPow ??= new double[f.Length];
            for (var k = 1; k < f.Length; k++)
            {
                var d = Math.Pow(2 * Math.PI * f[k], w.Power);
                sigPow[k] += a[k] * a[k] / (d * d);
                noiPow[k] += b[k] * b[k] / (d * d);
            }
        }
        if (freq == null) return null;
        var fmax = Math.Min(s.MaxFrequencyHz, 0.4 * fs);
        var bins = new List<(double F, double A, double N)>();
        for (var fc = s.MinFrequencyHz; fc <= fmax; fc *= Math.Pow(2, 1.0 / 6))
        {
            var lo = fc / Math.Pow(2, 1.0 / 12);
            var hi = fc * Math.Pow(2, 1.0 / 12);
            double sa = 0, sn = 0;
            var c = 0;
            for (var k = 1; k < freq.Length; k++)
            {
                if (freq[k] < lo || freq[k] >= hi) continue;
                sa += sigPow![k];
                sn += noiPow![k];
                c++;
            }
            if (c == 0) continue;
            bins.Add((fc, Math.Sqrt(sa / c), Math.Sqrt(sn / c)));
        }
        if (bins.Count < 6) return null;
        // Longest continuous run of bins above the SNR threshold.
        var above = bins.Select(x => x.N > 0 && x.A / x.N >= s.MinSnr).ToArray();
        int bestStart = -1, bestLen = 0;
        for (var i = 0; i < above.Length;)
        {
            if (!above[i]) { i++; continue; }
            var j = i;
            while (j < above.Length && above[j]) j++;
            if (j - i > bestLen) { bestStart = i; bestLen = j - i; }
            i = j;
        }
        if (bestLen < 5) return null;
        var used = new bool[bins.Count];
        for (var i = bestStart; i < bestStart + bestLen; i++) used[i] = true;
        // By default every bin above the threshold is fitted: on real records a notch in the SNR
        // (a spectral hole of the source or a site resonance) otherwise cuts the band short.
        if (!s.ContiguousBand) used = above;
        var ratio = bins[bestStart + bestLen - 1].F / bins[bestStart].F;
        if (ratio < s.MinBandRatio) return null;
        if (s.SubtractNoise)
            for (var i = 0; i < bins.Count; i++)
            {
                var (bf, ba, bn) = bins[i];
                bins[i] = (bf, Math.Sqrt(Math.Max(ba * ba - bn * bn, 0.01 * ba * ba)), bn);
            }
        // The spectrum must decay with frequency over the band (Wei & Wiens 2018 check this with the
        // linear correlation of ln A and f): flat or rising spectra are noise, resonance or clipping.
        var xs = Enumerable.Range(bestStart, bestLen).Select(i => bins[i].F).ToArray();
        var ys = Enumerable.Range(bestStart, bestLen).Select(i => Math.Log(bins[i].A)).ToArray();
        var corr = Correlation(xs, ys);
        if (!(corr < -s.MinSpectralCorrelation)) return null;
        var snr = Enumerable.Range(bestStart, bestLen).Average(i => bins[i].A / Math.Max(1e-30, bins[i].N));
        return (bins.Select(x => x.F).ToArray(), bins.Select(x => x.A).ToArray(), bins.Select(x => x.N).ToArray(), used, snr, ratio);
    }

    /// <summary>
    /// Start of the S noise window: the last <see cref="TStarSettings.SNoiseSeconds"/> of P coda before
    /// the S window, and never within 0.5 s of the P onset. The coda decays, so a window reaching back
    /// to the P wave overstates the noise under S and, once subtracted, the high-frequency S amplitude
    /// (t*_S +4 to +11 ms on the synthetic).
    /// </summary>
    private static DateTime? SNoiseStart(DateTime? pTime, DateTime sTime, TStarSettings s)
    {
        if (pTime is not { } p) return null;
        var start = sTime.AddSeconds(-s.PreSeconds - 0.2 - s.SNoiseSeconds);
        var early = p.AddSeconds(0.5);
        return start > early ? start : early;
    }

    private sealed record Solution(double[] TStar, double[] Error, double[] LnOmega, double[] RecordMisfit, double Misfit, int N, int Parameters);

    /// <summary>
    /// Joint least squares for one event at fixed f_c and α: ln Ω₀ shared (or one per record) and one
    /// t* per record, t* ≥ 0 when asked. Returns null when there are too few data.
    /// </summary>
    private static Solution? Solve(List<Record> recs, double fc, double alpha, TStarSettings s)
    {
        var shared = s.SharedSourceLevel;
        var nRec = recs.Count;
        var nSrc = shared ? 1 : nRec;
        var nPar = nSrc + nRec;
        var rows = new List<(int Rec, double[] G, double Y)>();
        for (var i = 0; i < nRec; i++)
        {
            var x = recs[i];
            for (var k = 0; k < x.F.Length; k++)
            {
                if (!x.Used[k] || !(x.A[k] > 0)) continue;
                var g = new double[nPar];
                g[shared ? 0 : i] = 1;
                g[nSrc + i] = -Math.PI * Math.Pow(x.F[k], 1 - alpha);
                rows.Add((i, g, Math.Log(x.A[k]) + x.LnR + Math.Log(1 + x.F[k] * x.F[k] / (fc * fc))));
            }
        }
        if (rows.Count < nPar + 2) return null;
        var free = new bool[nPar];
        for (var j = 0; j < nPar; j++) free[j] = true;
        double[] m;
        // Active set: t* that come out negative are fixed at zero and the rest re-solved; a fixed one
        // is released when its gradient says the misfit would drop by raising it (Lawson & Hanson 1974).
        for (var pass = 0; ; pass++)
        {
            m = LeastSquares(rows, free, nPar);
            if (!s.NonNegative || pass > 2 * nRec) break;
            var worst = -1;
            for (var j = nSrc; j < nPar; j++)
                if (free[j] && m[j] < 0 && (worst < 0 || m[j] < m[worst])) worst = j;
            if (worst >= 0) { free[worst] = false; continue; }
            var release = -1;
            double bestGrad = 0;
            for (var j = nSrc; j < nPar; j++)
            {
                if (free[j]) continue;
                double grad = 0;
                foreach (var (_, g, y) in rows) grad += g[j] * (y - Dot(g, m));
                if (grad > bestGrad) { bestGrad = grad; release = j; }
            }
            if (release < 0) break;
            free[release] = true;
        }
        var recMisfit = new double[nRec];
        var recCount = new int[nRec];
        double rss = 0;
        foreach (var (rec, g, y) in rows)
        {
            var r = y - Dot(g, m);
            rss += r * r;
            recMisfit[rec] += r * r;
            recCount[rec]++;
        }
        for (var i = 0; i < nRec; i++) recMisfit[i] = Math.Sqrt(recMisfit[i] / Math.Max(1, recCount[i]));
        // Standard errors from the covariance σ²(GᵀG)⁻¹ of the free parameters.
        var sigma2 = rss / Math.Max(1, rows.Count - free.Count(x => x));
        var cov = InverseNormal(rows, free, nPar);
        var err = new double[nRec];
        for (var i = 0; i < nRec; i++) err[i] = free[nSrc + i] ? Math.Sqrt(Math.Max(0, sigma2 * cov[nSrc + i, nSrc + i])) : Math.Sqrt(sigma2 / Math.Max(1, recCount[i]));
        var lnO = Enumerable.Range(0, nRec).Select(i => m[shared ? 0 : i]).ToArray();
        return new Solution(m.Skip(nSrc).ToArray(), err, lnO, recMisfit, rss / rows.Count, rows.Count, free.Count(x => x));
    }

    private static double Dot(double[] a, double[] b)
    {
        double s = 0;
        for (var i = 0; i < a.Length; i++) s += a[i] * b[i];
        return s;
    }

    /// <summary>Normal-equation least squares on the free parameters (fixed ones are zero).</summary>
    private static double[] LeastSquares(List<(int Rec, double[] G, double Y)> rows, bool[] free, int n)
    {
        var idx = Enumerable.Range(0, n).Where(j => free[j]).ToArray();
        var k = idx.Length;
        var ata = new double[k, k];
        var atb = new double[k];
        foreach (var (_, g, y) in rows)
            for (var a = 0; a < k; a++)
            {
                var ga = g[idx[a]];
                if (ga == 0) continue;
                atb[a] += ga * y;
                for (var b = 0; b < k; b++) ata[a, b] += ga * g[idx[b]];
            }
        var x = Cholesky(ata, atb);
        var m = new double[n];
        for (var a = 0; a < k; a++) m[idx[a]] = x[a];
        return m;
    }

    private static double[,] InverseNormal(List<(int Rec, double[] G, double Y)> rows, bool[] free, int n)
    {
        var idx = Enumerable.Range(0, n).Where(j => free[j]).ToArray();
        var k = idx.Length;
        var ata = new double[k, k];
        foreach (var (_, g, _) in rows)
            for (var a = 0; a < k; a++)
            {
                var ga = g[idx[a]];
                if (ga == 0) continue;
                for (var b = 0; b < k; b++) ata[a, b] += ga * g[idx[b]];
            }
        var inv = new double[n, n];
        for (var c = 0; c < k; c++)
        {
            var e = new double[k];
            e[c] = 1;
            var col = Cholesky(ata, e);
            for (var r = 0; r < k; r++) inv[idx[r], idx[c]] = col[r];
        }
        return inv;
    }

    /// <summary>Solves a symmetric positive-definite system (a tiny ridge keeps it definite).</summary>
    private static double[] Cholesky(double[,] a, double[] b)
    {
        var n = b.Length;
        var l = new double[n, n];
        var trace = 0.0;
        for (var i = 0; i < n; i++) trace += a[i, i];
        var ridge = 1e-12 * Math.Max(1e-30, trace / Math.Max(1, n));
        for (var i = 0; i < n; i++)
        for (var j = 0; j <= i; j++)
        {
            var sum = a[i, j] + (i == j ? ridge : 0);
            for (var k = 0; k < j; k++) sum -= l[i, k] * l[j, k];
            l[i, j] = i == j ? Math.Sqrt(Math.Max(sum, 1e-300)) : sum / l[j, j];
        }
        var y = new double[n];
        for (var i = 0; i < n; i++)
        {
            var sum = b[i];
            for (var k = 0; k < i; k++) sum -= l[i, k] * y[k];
            y[i] = sum / l[i, i];
        }
        var x = new double[n];
        for (var i = n - 1; i >= 0; i--)
        {
            var sum = y[i];
            for (var k = i + 1; k < n; k++) sum -= l[k, i] * x[k];
            x[i] = sum / l[i, i];
        }
        return x;
    }

    private static double Correlation(double[] x, double[] y)
    {
        double mx = x.Average(), my = y.Average(), sxy = 0, sxx = 0, syy = 0;
        for (var i = 0; i < x.Length; i++)
        {
            sxy += (x[i] - mx) * (y[i] - my);
            sxx += (x[i] - mx) * (x[i] - mx);
            syy += (y[i] - my) * (y[i] - my);
        }
        return sxx > 0 && syy > 0 ? sxy / Math.Sqrt(sxx * syy) : 0;
    }
}
