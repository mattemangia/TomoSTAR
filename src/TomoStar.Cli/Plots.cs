// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using System.Globalization;
using System.Text;
using TomoStar.Core.Tomography;

namespace TomoStar.Cli;

/// <summary>
/// Small vector figures written next to the results (SVG, readable in any browser and importable in
/// Inkscape or Illustrator for a publication figure). Only the L-curve is drawn here; the volumes
/// are meant for QUIVER, ParaView or GMT.
/// </summary>
public static class Plots
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>
    /// The L-curves of a trade-off run: log residual norm against log model norm, one curve per
    /// smoothing value with the damping growing along it, each point labelled with its damping, the
    /// corner of each curve filled and the recommended pair circled.
    /// </summary>
    public static void LCurveSvg(string path, string title, IReadOnlyList<TravelTimeTomography.TradeOffSummary> points, (double Damping, double Smoothing) chosen)
    {
        const double w = 720, h = 520, left = 80, right = 170, top = 50, bottom = 60;
        var pts = points.Where(p => p.ResidualNorm > 0 && p.ModelNorm > 0).ToList();
        if (pts.Count == 0) return;
        double X(double r) => Math.Log10(r);
        double Y(double m) => Math.Log10(m);
        double x0 = pts.Min(p => X(p.ResidualNorm)), x1 = pts.Max(p => X(p.ResidualNorm));
        double y0 = pts.Min(p => Y(p.ModelNorm)), y1 = pts.Max(p => Y(p.ModelNorm));
        if (x1 - x0 < 1e-9) { x0 -= 0.5; x1 += 0.5; }
        if (y1 - y0 < 1e-9) { y0 -= 0.5; y1 += 0.5; }
        var padX = 0.06 * (x1 - x0);
        var padY = 0.06 * (y1 - y0);
        x0 -= padX; x1 += padX; y0 -= padY; y1 += padY;
        double Px(double v) => left + (X(v) - x0) / (x1 - x0) * (w - left - right);
        double Py(double v) => h - bottom - (Y(v) - y0) / (y1 - y0) * (h - top - bottom);
        string[] colours = ["#1f5fa8", "#c2410c", "#15803d", "#7e22ce", "#b91c1c", "#0e7490", "#a16207"];

        var sb = new StringBuilder();
        sb.Append(string.Create(Inv, $"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{w}\" height=\"{h}\" viewBox=\"0 0 {w} {h}\" font-family=\"Helvetica, Arial, sans-serif\" font-size=\"12\">\n"));
        sb.Append("<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n");
        sb.Append(string.Create(Inv, $"<text x=\"{left}\" y=\"28\" font-size=\"15\">{Escape(title)}</text>\n"));
        sb.Append(string.Create(Inv, $"<rect x=\"{left}\" y=\"{top}\" width=\"{w - left - right}\" height=\"{h - top - bottom}\" fill=\"none\" stroke=\"#333\"/>\n"));
        // Axes: ticks at every power of ten and at 2 and 5 when the range is short.
        foreach (var (lo, hi, horizontal) in new[] { (x0, x1, true), (y0, y1, false) })
        {
            var decades = hi - lo;
            var mantissas = decades < 1.5 ? new[] { 1.0, 2, 5 } : new[] { 1.0 };
            for (var e = Math.Floor(lo); e <= Math.Ceiling(hi); e++)
            foreach (var m in mantissas)
            {
                var v = Math.Log10(m) + e;
                if (v < lo || v > hi) continue;
                var label = (m * Math.Pow(10, e)).ToString("G3", Inv);
                if (horizontal)
                {
                    var px = left + (v - x0) / (x1 - x0) * (w - left - right);
                    sb.Append(string.Create(Inv, $"<line x1=\"{px:0.#}\" y1=\"{h - bottom}\" x2=\"{px:0.#}\" y2=\"{h - bottom + 5}\" stroke=\"#333\"/><text x=\"{px:0.#}\" y=\"{h - bottom + 18}\" text-anchor=\"middle\">{label}</text>\n"));
                }
                else
                {
                    var py = h - bottom - (v - y0) / (y1 - y0) * (h - top - bottom);
                    sb.Append(string.Create(Inv, $"<line x1=\"{left - 5}\" y1=\"{py:0.#}\" x2=\"{left}\" y2=\"{py:0.#}\" stroke=\"#333\"/><text x=\"{left - 8}\" y=\"{py + 4:0.#}\" text-anchor=\"end\">{label}</text>\n"));
                }
            }
        }
        sb.Append(string.Create(Inv, $"<text x=\"{left + (w - left - right) / 2}\" y=\"{h - 18}\" text-anchor=\"middle\">weighted residual norm ||W(Gm - d)||</text>\n"));
        sb.Append(string.Create(Inv, $"<text transform=\"translate(22 {top + (h - top - bottom) / 2}) rotate(-90)\" text-anchor=\"middle\">model norm ||m||</text>\n"));

        var k = 0;
        foreach (var group in pts.GroupBy(p => p.Smoothing).OrderBy(g => g.Key))
        {
            var colour = colours[k % colours.Length];
            var list = group.OrderBy(p => p.Damping).ToList();
            var corner = list[Math.Clamp(TravelTimeTomography.Corner(list.Select(p => (p.ResidualNorm, p.ModelNorm)).ToList()), 0, list.Count - 1)];
            sb.Append($"<polyline fill=\"none\" stroke=\"{colour}\" stroke-width=\"1.5\" points=\"");
            foreach (var p in list) sb.Append(string.Create(Inv, $"{Px(p.ResidualNorm):0.#},{Py(p.ModelNorm):0.#} "));
            sb.Append("\"/>\n");
            foreach (var p in list)
            {
                var isCorner = p == corner;
                sb.Append(string.Create(Inv, $"<circle cx=\"{Px(p.ResidualNorm):0.#}\" cy=\"{Py(p.ModelNorm):0.#}\" r=\"{(isCorner ? 5 : 3)}\" fill=\"{(isCorner ? colour : "white")}\" stroke=\"{colour}\"/>"));
                sb.Append(string.Create(Inv, $"<text x=\"{Px(p.ResidualNorm) + 6:0.#}\" y=\"{Py(p.ModelNorm) - 5:0.#}\" font-size=\"9\" fill=\"{colour}\">{p.Damping:G4}</text>\n"));
                if (p.Damping == chosen.Damping && p.Smoothing == chosen.Smoothing)
                    sb.Append(string.Create(Inv, $"<circle cx=\"{Px(p.ResidualNorm):0.#}\" cy=\"{Py(p.ModelNorm):0.#}\" r=\"10\" fill=\"none\" stroke=\"black\" stroke-width=\"1.5\"/>\n"));
            }
            var ly = top + 20 + 20 * k;
            sb.Append(string.Create(Inv, $"<line x1=\"{w - right + 20}\" y1=\"{ly}\" x2=\"{w - right + 45}\" y2=\"{ly}\" stroke=\"{colour}\" stroke-width=\"2\"/><text x=\"{w - right + 52}\" y=\"{ly + 4}\">smoothing {group.Key:G4}</text>\n"));
            k++;
        }
        var cy = top + 20 + 20 * k + 10;
        sb.Append(string.Create(Inv, $"<circle cx=\"{w - right + 32}\" cy=\"{cy}\" r=\"8\" fill=\"none\" stroke=\"black\" stroke-width=\"1.5\"/><text x=\"{w - right + 52}\" y=\"{cy + 4}\">chosen</text>\n"));
        sb.Append(string.Create(Inv, $"<text x=\"{w - right + 20}\" y=\"{cy + 24}\" font-size=\"10\">filled: corner</text>\n"));
        sb.Append(string.Create(Inv, $"<text x=\"{w - right + 20}\" y=\"{cy + 38}\" font-size=\"10\">labels: damping</text>\n"));
        sb.Append("</svg>\n");
        File.WriteAllText(path, sb.ToString());
    }

    private static string Escape(string s) => s.Replace("&", "&amp;").Replace("<", "&lt;").Replace(">", "&gt;");
}
