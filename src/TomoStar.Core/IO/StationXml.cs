// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using System.Globalization;
using System.Xml;
using System.Xml.Linq;
using TomoStar.Core.Model;
using TomoStar.Core.Signal;

namespace TomoStar.Core.IO;

/// <summary>
/// FDSN StationXML 1.x reader (https://www.fdsn.org/xml/station/). Reads networks, stations and,
/// at channel level, each channel's orientation, sample rate and overall instrument sensitivity; at
/// response level also the analogue poles-and-zeros stages that give the shape of the response.
/// Stations that appear several times (several epochs, or several data centres) are merged by
/// network and code, keeping the most recent coordinates. The instrument responses are returned by
/// channel (NET.STA.LOC.CHA) for the t* measurement, which needs ground displacement.
/// </summary>
public static class StationXml
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>The stations of a StationXML document and the responses of their channels.</summary>
    public sealed record Result(List<StationRecord> Stations, Dictionary<string, List<ChannelResponse>> Responses);

    /// <summary>One channel epoch and its response (null when only the sensitivity is known and it is 1).</summary>
    public sealed record ChannelResponse(string Nslc, DateTime? Start, DateTime? End, double Azimuth, double Dip, InstrumentResponse? Response);

    public static Result Parse(Stream stream)
    {
        var doc = XDocument.Load(XmlReader.Create(stream, new XmlReaderSettings { DtdProcessing = DtdProcessing.Prohibit }));
        var byId = new Dictionary<string, (double Lat, double Lon, double Elev, DateTime Start)>(StringComparer.OrdinalIgnoreCase);
        var responses = new Dictionary<string, List<ChannelResponse>>(StringComparer.OrdinalIgnoreCase);
        foreach (var net in doc.Descendants().Where(e => e.Name.LocalName == "Network"))
        {
            var netCode = net.Attribute("code")?.Value ?? "";
            foreach (var sta in net.Elements().Where(e => e.Name.LocalName == "Station"))
            {
                var code = sta.Attribute("code")?.Value ?? "";
                var id = $"{netCode}.{code}";
                var start = Date(sta.Attribute("startDate")?.Value) ?? DateTime.MinValue;
                var lat = Num(Text(sta, "Latitude")) ?? 0;
                var lon = Num(Text(sta, "Longitude")) ?? 0;
                var elev = Num(Text(sta, "Elevation")) ?? 0;
                // Coordinates of the latest station epoch.
                if (!byId.TryGetValue(id, out var old) || start >= old.Start) byId[id] = (lat, lon, elev, start);
                foreach (var ch in sta.Elements().Where(e => e.Name.LocalName == "Channel"))
                {
                    var cha = ch.Attribute("code")?.Value ?? "";
                    var loc = (ch.Attribute("locationCode")?.Value ?? "").Trim();
                    InstrumentResponse? response = null;
                    var sens = Child(Child(ch, "Response"), "InstrumentSensitivity");
                    if (sens != null)
                    {
                        var stageElements = Child(ch, "Response")!.Elements().Where(e => e.Name.LocalName == "Stage").ToList();
                        var stages = stageElements.Where(st => Child(st, "PolesZeros") != null)
                            .Select(st =>
                            {
                                var pz = PolesZeros(Child(st, "PolesZeros")!);
                                if (pz != null) pz.GainFrequency = Num(Text(Child(st, "StageGain"), "Frequency")) ?? double.NaN;
                                return pz;
                            }).Where(pz => pz != null).ToList();
                        var digital = stageElements.Select(DigitalStageOf).Where(d => d != null).Select(d => d!).ToList();
                        // The product of the stage gains, when every stage declares one (as evalresp uses them).
                        var gains = stageElements.Select(st => Num(Text(Child(st, "StageGain"), "Value"))).ToList();
                        response = new InstrumentResponse
                        {
                            InputUnits = Text(Child(sens, "InputUnits"), "Name") ?? "",
                            Sensitivity = Num(Text(sens, "Value")) ?? 1,
                            SensitivityFrequency = Num(Text(sens, "Frequency")) ?? 1,
                            StageGain = gains.Count > 0 && gains.All(g => g is > 0) ? gains.Aggregate(1.0, (a, g) => a * g!.Value) : double.NaN,
                            Stages = stages!,
                            DigitalStages = digital
                        };
                        if (response.Stages.Count == 0 && (response.Sensitivity == 1 || response.InputUnits.Length == 0)) response = null;
                    }
                    var nslc = $"{netCode}.{code}.{loc}.{cha}";
                    if (!responses.TryGetValue(nslc, out var list)) responses[nslc] = list = [];
                    list.Add(new ChannelResponse(nslc, Date(ch.Attribute("startDate")?.Value), Date(ch.Attribute("endDate")?.Value),
                        Num(Text(ch, "Azimuth")) ?? 0, Num(Text(ch, "Dip")) ?? 0, response));
                }
            }
        }
        var stations = byId.Select(kv => new StationRecord { Id = kv.Key, Lat = kv.Value.Lat, Lon = kv.Value.Lon, ElevationM = kv.Value.Elev }).ToList();
        return new Result(stations, responses);
    }

    /// <summary>Reads a StationXML file.</summary>
    public static Result ReadFile(string path)
    {
        using var fs = File.OpenRead(path);
        return Parse(fs);
    }

    /// <summary>The response of a channel at a time (the epoch that contains it, else the last one).</summary>
    public static InstrumentResponse? ResponseAt(Dictionary<string, List<ChannelResponse>> responses, string nslc, DateTime time)
    {
        if (!responses.TryGetValue(nslc, out var list) || list.Count == 0) return null;
        var hit = list.FirstOrDefault(c => (c.Start == null || c.Start <= time) && (c.End == null || c.End >= time));
        return (hit ?? list[^1]).Response;
    }

    /// <summary>An analogue poles-and-zeros stage; digital (Z-transform) stages are skipped.</summary>
    private static PolesZerosStage? PolesZeros(XElement pz)
    {
        var type = Text(pz, "PzTransferFunctionType") ?? "";
        if (type.Contains("DIGITAL", StringComparison.OrdinalIgnoreCase)) return null;
        double[] Pair(XElement e) => [Num(Text(e, "Real")) ?? 0, Num(Text(e, "Imaginary")) ?? 0];
        return new PolesZerosStage
        {
            TransferFunction = type.Length > 0 ? type : "LAPLACE (RADIANS/SECOND)",
            NormalizationFactor = Num(Text(pz, "NormalizationFactor")) ?? 1,
            NormalizationFrequency = Num(Text(pz, "NormalizationFrequency")) ?? 1,
            Zeros = pz.Elements().Where(e => e.Name.LocalName == "Zero").Select(Pair).ToList(),
            Poles = pz.Elements().Where(e => e.Name.LocalName == "Pole").Select(Pair).ToList()
        };
    }

    /// <summary>A FIR or coefficients stage with its coefficients (symmetric halves expanded), or null.</summary>
    private static DigitalStage? DigitalStageOf(XElement stage)
    {
        var dec = Child(stage, "Decimation");
        var rate = Num(Text(dec, "InputSampleRate")) ?? 0;
        var delay = Num(Text(dec, "Delay")) ?? 0;
        var correction = Num(Text(dec, "Correction")) ?? 0;
        if (Child(stage, "FIR") is { } fir)
        {
            var c = fir.Elements().Where(e => e.Name.LocalName == "NumeratorCoefficient").Select(e => Num(e.Value.Trim()) ?? 0).ToList();
            var symmetry = (Text(fir, "Symmetry") ?? "NONE").ToUpperInvariant();
            // EVEN: the coefficients are the first half of an even-length filter; ODD: of an odd-length one.
            if (symmetry == "EVEN") c.AddRange(Enumerable.Reverse(c).ToList());
            else if (symmetry == "ODD") c.AddRange(Enumerable.Reverse(c).Skip(1).ToList());
            return c.Count > 1 ? new DigitalStage { InputSampleRate = rate, Coefficients = c.ToArray(), Delay = delay, Correction = correction } : null;
        }
        if (Child(stage, "Coefficients") is { } cf && !cf.Elements().Any(e => e.Name.LocalName == "Denominator"))
        {
            var c = cf.Elements().Where(e => e.Name.LocalName == "Numerator").Select(e => Num(e.Value.Trim()) ?? 0).ToArray();
            return c.Length > 1 ? new DigitalStage { InputSampleRate = rate, Coefficients = c, Delay = delay, Correction = correction } : null;
        }
        return null;
    }

    private static XElement? Child(XElement? e, string name) => e?.Elements().FirstOrDefault(x => x.Name.LocalName == name);
    private static string? Text(XElement? e, string name) => Child(e, name)?.Value.Trim();
    private static double? Num(string? s) => double.TryParse(s, NumberStyles.Float, Inv, out var v) ? v : null;

    private static DateTime? Date(string? s) =>
        DateTime.TryParse(s, Inv, DateTimeStyles.AdjustToUniversal | DateTimeStyles.AssumeUniversal, out var d)
            ? DateTime.SpecifyKind(d, DateTimeKind.Utc) : null;
}
