// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

namespace TomoStar.Core.Source;

/// <summary>
/// What is read off a reduced stress tensor for maps: the azimuth of the maximum horizontal
/// compression and the faulting regime of the World Stress Map.
/// </summary>
public static class StressField
{
    /// <summary>
    /// Azimuth of S_Hmax, degrees in [0, 180): the horizontal direction of the largest normal
    /// compression, i.e. the maximum of n·S·n over horizontal unit vectors n (Lund &amp; Townend 2007,
    /// doi:10.1111/j.1365-246X.2007.03468.x). Tensor in east–north–down, compression positive.
    /// </summary>
    public static double SHmaxAzimuth(StressTensorResult t)
    {
        double see = t.Sxx, snn = t.Syy, sen = t.Sxy;
        // σ(α) = Snn cos²α + 2 Sne sinα cosα + See sin²α, α clockwise from north.
        var alpha = 0.5 * Math.Atan2(2 * sen, snn - see) * 180 / Math.PI;
        return ((alpha % 180) + 180) % 180;
    }

    /// <summary>
    /// Faulting regime from the plunges of the principal axes (Zoback 1992,
    /// doi:10.1029/92JB00132): NF, NS, SS, TS, TF, or U when the axes fit none of the classes.
    /// </summary>
    public static string Regime(double sigma1Plunge, double sigma2Plunge, double sigma3Plunge)
    {
        double p1 = sigma1Plunge, p2 = sigma2Plunge, p3 = sigma3Plunge;
        if (p1 >= 52 && p3 <= 35) return "NF";
        if (p1 >= 40 && p1 < 52 && p3 <= 20) return "NS";
        if (p1 < 40 && p3 <= 20 && p2 >= 45) return "SS";
        if (p1 <= 20 && p3 < 40 && p2 >= 45) return "SS";
        if (p1 <= 20 && p3 >= 40 && p3 < 52) return "TS";
        if (p1 <= 35 && p3 >= 52) return "TF";
        return "U";
    }

    public static string Regime(StressTensorResult t) => Regime(t.Sigma1.PlungeDeg, t.Sigma2.PlungeDeg, t.Sigma3.PlungeDeg);

    public static string RegimeName(string code) => code switch
    {
        "NF" => "Normal faulting",
        "NS" => "Normal / strike-slip",
        "SS" => "Strike-slip",
        "TS" => "Thrust / strike-slip",
        "TF" => "Thrust faulting",
        _ => "Unknown"
    };
}
