// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using TomoStar.Core.Geo;
using TomoStar.Core.Model;

namespace TomoStar.Core.Tomography;

/// <summary>Vp/Vs read from the data, with its spread and the number of pairs behind it.</summary>
public sealed record WadatiEstimate(double VpVs, double StandardError, double ScatterSeconds, int Pairs, int Events, int Rejected);

/// <summary>
/// Vp/Vs from the S and P times alone (Wadati 1933), in the form free of origin times (Chatelain
/// 1978): within each event, tS − ⟨tS⟩ = (Vp/Vs)(tP − ⟨tP⟩) over its stations, whatever the origin
/// time and the velocity model, so long as P and S follow similar paths. One slope through all the
/// events, fitted by least squares with pairs off the line by more than 3 robust sigmas left out.
/// This is the Vp/Vs the data carry on average; it is the starting value of a local inversion
/// when the regional model's ratio is not constrained for the area (as near the surface).
/// </summary>
public static class Wadati
{
    public static WadatiEstimate Estimate(ObservationSet data, double maxDistanceKm = 100, int minStations = 3)
    {
        var xs = new List<double>();
        var ys = new List<double>();
        var events = 0;
        foreach (var byEvent in data.Arrivals.GroupBy(a => a.Event))
        {
            var ev = data.Events[byEvent.Key];
            var pairs = new List<(double Tp, double Ts)>();
            foreach (var st in byEvent.GroupBy(a => a.Station))
            {
                var p = st.Where(a => a.Phase == Phase.P).MinBy(a => a.Sigma);
                var s = st.Where(a => a.Phase == Phase.S).MinBy(a => a.Sigma);
                if (p == null || s == null) continue;
                var sta = data.Stations[st.Key];
                if (GeoMath.SurfaceDistanceKm(ev.Lon, ev.Lat, sta.Lon, sta.Lat) > maxDistanceKm) continue;
                pairs.Add((p.Time, s.Time));
            }
            if (pairs.Count < minStations) continue;
            events++;
            double mp = pairs.Average(q => q.Tp), ms = pairs.Average(q => q.Ts);
            foreach (var (tp, ts) in pairs) { xs.Add(tp - mp); ys.Add(ts - ms); }
        }
        if (xs.Count < 3) return new WadatiEstimate(double.NaN, double.NaN, double.NaN, xs.Count, events, 0);

        var use = Enumerable.Repeat(true, xs.Count).ToArray();
        double k = 0, scatter = 0;
        for (var pass = 0; pass < 5; pass++)
        {
            double sxy = 0, sxx = 0;
            for (var i = 0; i < xs.Count; i++) if (use[i]) { sxy += xs[i] * ys[i]; sxx += xs[i] * xs[i]; }
            k = sxy / sxx;
            var res = Enumerable.Range(0, xs.Count).Select(i => ys[i] - k * xs[i]).ToArray();
            var kept = Enumerable.Range(0, xs.Count).Where(i => use[i]).Select(i => Math.Abs(res[i])).OrderBy(v => v).ToArray();
            scatter = 1.4826 * kept[kept.Length / 2];
            var changed = false;
            for (var i = 0; i < xs.Count; i++)
            {
                var keep = Math.Abs(res[i]) <= 3 * scatter;
                changed |= keep != use[i];
                use[i] = keep;
            }
            if (!changed) break;
        }
        double sx = 0, rss = 0;
        var n = 0;
        for (var i = 0; i < xs.Count; i++)
            if (use[i]) { sx += xs[i] * xs[i]; rss += Math.Pow(ys[i] - k * xs[i], 2); n++; }
        var se = n > 1 ? Math.Sqrt(rss / (n - 1) / sx) : double.NaN;
        return new WadatiEstimate(k, se, scatter, n, events, xs.Count - n);
    }

    /// <summary>
    /// The starting Vs of <paramref name="settings"/> applied to a starting model on the grid and to
    /// the 1-D model outside it: unchanged for <see cref="StartingVpVs.FromModel"/>, else Vs = Vp / k
    /// with k from the data's Wadati diagram or the given value. <paramref name="vs"/> is changed in place.
    /// </summary>
    public static VelocityModel1D ApplyStart(TomographySettings settings, double[] vp, double[] vs, VelocityModel1D background, ObservationSet data, Action<string>? log = null)
    {
        if (settings.StartVpVs == StartingVpVs.FromModel) return background;
        var k = settings.StartVpVsValue;
        if (settings.StartVpVs == StartingVpVs.FromData)
        {
            var w = Estimate(data);
            if (!double.IsFinite(w.VpVs)) throw new InvalidOperationException("No event has P and S picks at three stations: the Wadati diagram cannot give Vp/Vs.");
            log?.Invoke(string.Create(System.Globalization.CultureInfo.InvariantCulture,
                $"Wadati diagram: Vp/Vs {w.VpVs:0.000} ± {w.StandardError:0.000} from {w.Pairs} S-P pairs of {w.Events} events ({w.Rejected} off the line), scatter {w.ScatterSeconds:0.000} s."));
            k = w.VpVs;
        }
        for (var i = 0; i < vs.Length; i++) vs[i] = vp[i] / k;
        return WithVpVs(background, k);
    }

    /// <summary>The model with its Vp kept and Vs set to Vp / <paramref name="vpVs"/> at every node.</summary>
    public static VelocityModel1D WithVpVs(VelocityModel1D model, double vpVs)
    {
        var m = new VelocityModel1D { Name = $"{model.Name}, Vp/Vs {vpVs:0.000}", Reference = model.Reference, Doi = model.Doi };
        m.Nodes.AddRange(model.Nodes.Select(n => n with { Vs = n.Vp / vpVs }));
        return m;
    }
}
