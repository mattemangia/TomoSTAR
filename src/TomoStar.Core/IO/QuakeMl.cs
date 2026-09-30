using System.Globalization;
using System.Xml;
using System.Xml.Linq;
using TomoStar.Core.Model;

namespace TomoStar.Core.IO;

/// <summary>
/// QuakeML 1.2 (BED) reader. Reads the preferred origin and magnitude of each event and, when the
/// service returned them (<c>includearrivals=true</c>), its picks with the phase taken from the
/// arrival that references them. Elements are matched by local name, so the several namespaces in
/// use by FDSN data centres all read.
/// </summary>
public static class QuakeMl
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    public static List<EventRecord> Parse(Stream stream, string source)
    {
        var settings = new XmlReaderSettings { DtdProcessing = DtdProcessing.Prohibit, IgnoreWhitespace = true };
        using var reader = XmlReader.Create(stream, settings);
        var events = new List<EventRecord>();
        // Stream event by event: a catalogue response can be hundreds of MB.
        // XNode.ReadFrom leaves the reader ON the node after the element, so the loop must not call
        // Read() after it: doing so would skip the next sibling <event>.
        reader.MoveToContent();
        while (!reader.EOF)
        {
            if (reader.NodeType != XmlNodeType.Element || reader.LocalName != "event")
            {
                reader.Read();
                continue;
            }
            if (XNode.ReadFrom(reader) is not XElement ev) continue;
            try
            {
                var parsed = ParseEvent(ev, source);
                if (parsed != null) events.Add(parsed);
            }
            catch (FormatException)
            {
                // A malformed event is skipped; the rest of the catalogue still reads.
            }
        }
        return events;
    }

    public static List<EventRecord> Parse(string xml, string source)
    {
        using var ms = new MemoryStream(System.Text.Encoding.UTF8.GetBytes(xml));
        return Parse(ms, source);
    }

    private static EventRecord? ParseEvent(XElement ev, string source)
    {
        var id = Attr(ev, "publicID") ?? Guid.NewGuid().ToString("N");
        var preferredOrigin = Text(ev, "preferredOriginID");
        var origins = Children(ev, "origin").ToList();
        if (origins.Count == 0) return null;
        var origin = origins.FirstOrDefault(o => Attr(o, "publicID") == preferredOrigin) ?? origins[0];

        var lat = Double(Value(origin, "latitude"));
        var lon = Double(Value(origin, "longitude"));
        var time = Time(Value(origin, "time"));
        if (lat is null || lon is null || time is null) return null;
        var depthM = Double(Value(origin, "depth"));

        var hyp = new Hypocentre
        {
            OriginTime = time.Value,
            Lat = lat.Value,
            Lon = lon.Value,
            DepthKm = (depthM ?? 0) / 1000.0
        };
        var quality = Child(origin, "quality");
        if (quality != null)
        {
            hyp.Rms = Double(Text(quality, "standardError")) ?? double.NaN;
            hyp.AzimuthalGapDeg = Double(Text(quality, "azimuthalGap")) ?? double.NaN;
            hyp.PhaseCount = (int)(Double(Text(quality, "usedPhaseCount")) ?? 0);
        }
        var herr = Child(origin, "originUncertainty");
        if (herr != null)
        {
            // An ellipse (max, min, azimuth of max) when given; otherwise the circular uncertainty.
            var maxH = Double(Text(herr, "maxHorizontalUncertainty"));
            var minH = Double(Text(herr, "minHorizontalUncertainty"));
            var azH = Double(Text(herr, "azimuthMaxHorizontalUncertainty"));
            if (maxH != null && minH != null && azH != null)
            {
                hyp.ErrorHorizontalKm = maxH.Value / 1000.0;
                hyp.ErrorMinorKm = minH.Value / 1000.0;
                hyp.ErrorAzimuthDeg = (azH.Value % 180 + 180) % 180;
            }
            else hyp.ErrorHorizontalKm = (Double(Text(herr, "horizontalUncertainty")) ?? double.NaN) / 1000.0;
        }
        var depthUnc = Double(Child(origin, "depth") is { } de ? Text(de, "uncertainty") : null);
        if (depthUnc != null) hyp.ErrorDepthKm = depthUnc.Value / 1000.0;
        var timeUnc = Double(Child(origin, "time") is { } ot ? Text(ot, "uncertainty") : null);
        if (timeUnc != null) hyp.ErrorTimeS = timeUnc.Value;

        var e = new EventRecord { Id = CleanId(id), OriginTime = hyp.OriginTime, Lon = hyp.Lon, Lat = hyp.Lat, DepthKm = hyp.DepthKm };

        var preferredMag = Text(ev, "preferredMagnitudeID");
        var mags = Children(ev, "magnitude").ToList();
        var mag = mags.FirstOrDefault(m => Attr(m, "publicID") == preferredMag) ?? mags.FirstOrDefault();
        if (mag != null)
        {
            e.Magnitude = Double(Value(mag, "mag")) ?? double.NaN;
        }

        // Phase comes from the arrival (the associated, interpreted label); fall back to the
        // pick's own phaseHint. An event can carry the picks of several solutions (INGV returns
        // the automatic and the revised ones: up to three P picks per station, up to a second
        // apart), so when the preferred origin has arrivals only the picks it associates are
        // kept; the others belong to superseded solutions.
        var preferredArrivals = Children(origin, "arrival").ToList();
        var phaseByPick = new Dictionary<string, string>(StringComparer.Ordinal);
        var weightByPick = new Dictionary<string, double>(StringComparer.Ordinal);
        foreach (var arr in preferredArrivals.Count > 0 ? preferredArrivals : origins.SelectMany(o => Children(o, "arrival")))
        {
            var pickId = Text(arr, "pickID");
            var phase = Text(arr, "phase");
            if (!string.IsNullOrEmpty(pickId) && !string.IsNullOrEmpty(phase)) phaseByPick.TryAdd(pickId, phase);
            if (!string.IsNullOrEmpty(pickId) && Double(Text(arr, "timeWeight")) is { } w) weightByPick.TryAdd(pickId, w);
        }
        // Arrival weights run 0-1 in QuakeML; INGV writes them in per cent (0-100).
        var weightScale = weightByPick.Count > 0 && weightByPick.Values.Max() > 1 ? 100.0 : 1.0;
        foreach (var pk in Children(ev, "pick"))
        {
            var pickId = Attr(pk, "publicID") ?? "";
            if (preferredArrivals.Count > 0 && !phaseByPick.ContainsKey(pickId)) continue;
            var wf = Child(pk, "waveformID");
            var t = Time(Value(pk, "time"));
            if (wf == null || t == null) continue;
            phaseByPick.TryGetValue(pickId, out var label);
            label ??= Text(pk, "phaseHint");
            var phase = MapPhase(label);
            if (phase == null) continue;
            var unc = Double(Child(pk, "time") is { } te ? Text(te, "uncertainty") : null);
            var lower = Double(Child(pk, "time") is { } tl ? Text(tl, "lowerUncertainty") : null);
            var upper = Double(Child(pk, "time") is { } tu ? Text(tu, "upperUncertainty") : null);
            if (unc == null && lower != null && upper != null) unc = 0.5 * (lower + upper);
            // The agency's own arrival weight: a pick its locator left out (weight 0: a far station, a
            // late phase or a misread arrival, often seconds off) is kept but rejected (quality 4, no
            // weight in an inversion), the others get the HYPO71 class of their weight.
            var pickQuality = weightByPick.TryGetValue(pickId, out var aw)
                ? (aw / weightScale) switch { <= 0 => 4, >= 0.75 => 0, >= 0.5 => 1, >= 0.25 => 2, _ => 3 }
                : 0;
            e.Picks.Add(new PickRecord
            {
                StationId = $"{Attr(wf, "networkCode") ?? ""}.{Attr(wf, "stationCode") ?? ""}",
                Channel = Attr(wf, "channelCode") ?? "",
                Phase = phase.Value,
                Time = t.Value,
                Sigma = unc is > 0 ? unc.Value : (phase == Phase.P ? 0.1 : 0.2),
                Quality = pickQuality,
                Origin = Text(pk, "evaluationMode") == "automatic" ? PickOrigin.Automatic : PickOrigin.Catalog
            });
        }
        return e;
    }

    /// <summary>
    /// Maps a catalogue phase label to the first-arrival phase TomoSTAR inverts. Crustal (Pg/Sg),
    /// Moho-refracted (Pn/Sn) and generic (P/S, Pb/Sb, P*) first arrivals are kept; depth phases,
    /// reflections and core phases are not first arrivals on a regional grid and are dropped.
    /// </summary>
    public static Phase? MapPhase(string? label)
    {
        if (string.IsNullOrWhiteSpace(label)) return null;
        var l = label.Trim();
        return l switch
        {
            "P" or "Pg" or "Pn" or "Pb" or "P*" or "Pp" or "PG" or "PN" or "p" => Phase.P,
            "S" or "Sg" or "Sn" or "Sb" or "S*" or "SG" or "SN" or "s" => Phase.S,
            _ => null
        };
    }

    private static string CleanId(string id)
    {
        // smi:org.gfz-potsdam.de/geofon/gfz2016qhoa → gfz2016qhoa. A query-style id names the event in
        // one of its parameters: quakeml:earthquake.usgs.gov/fdsnws/event/1/query?eventid=us7000d18q&format=quakeml
        // → us7000d18q (the last parameter is the format, the same for every event).
        var s = id;
        var q = s.IndexOf('?');
        if (q >= 0)
        {
            string? fallback = null;
            foreach (var part in s[(q + 1)..].Split('&', ';'))
            {
                var eq = part.IndexOf('=');
                if (eq <= 0 || eq == part.Length - 1) continue;
                var key = part[..eq].Trim().ToLowerInvariant();
                var value = Uri.UnescapeDataString(part[(eq + 1)..].Trim());
                if (key is "eventid" or "evid" or "id" or "originid" or "pickid" or "amplitudeid") return value;
                if (key is not ("format" or "includeallorigins" or "includearrivals" or "includeallmagnitudes")) fallback ??= value;
            }
            if (fallback != null) return fallback;
            s = s[..q];
        }
        var slash = s.LastIndexOf('/');
        if (slash >= 0 && slash < s.Length - 1) s = s[(slash + 1)..];
        var eqIdx = s.IndexOf('=');
        if (eqIdx >= 0 && eqIdx < s.Length - 1)
        {
            var key = s[..eqIdx].Trim().ToLowerInvariant();
            if (key is "eventid" or "evid" or "id" or "originid" or "origid" or "pickid" or "amplitudeid" or "arrival" or "arrivalid")
                return s[(eqIdx + 1)..].Trim();
        }
        return s;
    }

    private static IEnumerable<XElement> Children(XElement e, string name) => e.Elements().Where(x => x.Name.LocalName == name);
    private static XElement? Child(XElement e, string name) => e.Elements().FirstOrDefault(x => x.Name.LocalName == name);
    private static string? Text(XElement e, string name) => Child(e, name)?.Value.Trim();
    private static string? Value(XElement e, string name) => Child(e, name) is { } c ? Text(c, "value") : null;

    private static string? Attr(XElement e, string name) =>
        e.Attribute(name)?.Value ?? e.Attributes().FirstOrDefault(a => a.Name.LocalName == name)?.Value;

    private static double? Double(string? s) =>
        double.TryParse(s, NumberStyles.Float, Inv, out var v) ? v : null;

    internal static DateTime? Time(string? s)
    {
        if (string.IsNullOrWhiteSpace(s)) return null;
        return DateTime.TryParse(s, Inv, DateTimeStyles.AdjustToUniversal | DateTimeStyles.AssumeUniversal, out var t)
            ? DateTime.SpecifyKind(t, DateTimeKind.Utc)
            : null;
    }

    /// <summary>Reads a QuakeML file.</summary>
    public static List<EventRecord> ReadFile(string path)
    {
        using var fs = File.OpenRead(path);
        return Parse(fs, Path.GetFileName(path));
    }
}
