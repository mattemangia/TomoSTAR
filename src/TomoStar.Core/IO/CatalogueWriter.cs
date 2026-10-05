// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using System.Globalization;
using System.Text;
using TomoStar.Core.Location;
using TomoStar.Core.Model;

namespace TomoStar.Core.IO;

/// <summary>
/// Writers of the input tables (stations, events, picks, t*), in the columns
/// <see cref="CatalogueReader"/> reads, so that the output of one step is the input of the next:
/// automatic picks feed the relocation, relocated events feed the tomography and the t* measurement,
/// t* feed the Q tomography.
/// </summary>
public static class CatalogueWriter
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    private static StreamWriter Open(string path)
    {
        var dir = Path.GetDirectoryName(Path.GetFullPath(path));
        if (dir != null) Directory.CreateDirectory(dir);
        return new StreamWriter(path, false, new UTF8Encoding(false));
    }

    /// <summary>Time as ISO 8601 UTC with 0.1 ms.</summary>
    public static string Time(DateTime t) => t.ToUniversalTime().ToString("yyyy-MM-ddTHH:mm:ss.ffffZ", Inv);

    public static void WriteStations(string path, IEnumerable<StationRecord> stations)
    {
        using var w = Open(path);
        w.WriteLine("id,lon,lat,elevation_m,correction_p_s,correction_s_s");
        foreach (var s in stations)
            w.WriteLine(string.Create(Inv, $"{s.Id},{s.Lon:0.######},{s.Lat:0.######},{s.ElevationM:0.#},{s.CorrectionP:0.####},{s.CorrectionS:0.####}"));
    }

    public static void WriteEvents(string path, IEnumerable<EventRecord> events)
    {
        using var w = Open(path);
        w.WriteLine("id,time,lon,lat,depth_km,magnitude,fixed");
        foreach (var e in events)
            w.WriteLine(string.Create(Inv,
                $"{e.Id},{Time(e.OriginTime)},{e.Lon:0.######},{e.Lat:0.######},{e.DepthKm:0.####},{(double.IsFinite(e.Magnitude) ? e.Magnitude.ToString("0.0#", Inv) : "")},{(e.Fixed ? 1 : 0)}"));
    }

    /// <summary>Every pick of every event with its absolute time, origin and, for automatic picks, channel and SNR.</summary>
    public static void WritePicks(string path, IEnumerable<EventRecord> events)
    {
        using var w = Open(path);
        w.WriteLine("event,station,phase,time,sigma_s,quality,origin,disabled,channel,snr");
        foreach (var e in events)
        foreach (var p in e.Picks)
            w.WriteLine(string.Create(Inv,
                $"{e.Id},{p.StationId},{p.Phase},{Time(p.Time)},{p.Sigma:0.####},{p.Quality},{p.Origin.ToString().ToLowerInvariant()},{(p.Disabled ? 1 : 0)},{p.Channel},{(double.IsFinite(p.Snr) ? p.Snr.ToString("0.#", Inv) : "")}"));
    }

    /// <summary>t* measurements with the details of their fit (the extra columns are ignored on reading).</summary>
    public static void WriteTStar(string path, IEnumerable<TStarMeasurement> measurements)
    {
        using var w = Open(path);
        w.WriteLine("event,station,phase,tstar_s,sigma_s,disabled,corner_hz,corner_low_hz,corner_high_hz,corner_resolved,fit_min_hz,fit_max_hz,snr,alpha");
        foreach (var m in measurements)
            w.WriteLine(string.Create(Inv,
                $"{m.EventId},{m.StationId},{m.Phase},{m.TStar:0.######},{Math.Max(0.001, m.Uncertainty):0.######},{(m.Disabled ? 1 : 0)},{m.CornerFrequencyHz:0.###},{m.CornerLowHz:0.###},{m.CornerHighHz:0.###},{(m.CornerResolved ? 1 : 0)},{m.FitMinHz:0.##},{m.FitMaxHz:0.##},{m.Snr:0.#},{m.Alpha:0.##}"));
    }

    /// <summary>
    /// Absolute locations with their errors, one row per event (failed events keep their start and
    /// say why), and the events as an event table for the next step.
    /// </summary>
    public static void WriteLocations(string path, IEnumerable<RelocatedEvent> located)
    {
        using var w = Open(path);
        w.WriteLine("id,time,lon,lat,depth_km,rms_s,phases,rejected,err_h_major_km,err_h_minor_km,err_h_azimuth_deg,err_depth_km,err_time_s,gap_deg,shift_km,status");
        foreach (var r in located)
        {
            var h = r.Result;
            w.WriteLine(string.Create(Inv,
                $"{r.Start.Id},{Time(h.OriginTime)},{h.Lon:0.######},{h.Lat:0.######},{h.DepthKm:0.####},{h.Rms:0.####},{r.Used},{r.Rejected},{h.ErrorHorizontalKm:0.###},{h.ErrorMinorKm:0.###},{h.ErrorAzimuthDeg:0.#},{h.ErrorDepthKm:0.###},{h.ErrorTimeS:0.####},{h.AzimuthalGapDeg:0.#},{Relocation.Shift(r):0.###},{r.Message.Replace(',', ';')}"));
        }
    }
}
