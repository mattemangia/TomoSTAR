// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using System.Text.Json;
using System.Text.Json.Serialization;
using System.Text.Json.Serialization.Metadata;
using TomoStar.Core.Geo;

namespace TomoStar.Core.Model;

/// <summary>Seismic phase of an arrival or of a t* measurement.</summary>
public enum Phase
{
    /// <summary>Compressional wave.</summary>
    P,

    /// <summary>Shear wave.</summary>
    S
}

/// <summary>What to do with events and stations outside the inversion grid.</summary>
public enum OutsideData
{
    /// <summary>Leave them out.</summary>
    Exclude,

    /// <summary>Use their data whose path crosses the grid; the forward grid is enlarged to hold them.</summary>
    WhenRaysCross
}

/// <summary>A geographic bounding box, degrees.</summary>
public sealed record GeoBox(double MinLon, double MaxLon, double MinLat, double MaxLat)
{
    [JsonIgnore] public double CenterLon => 0.5 * (MinLon + MaxLon);
    [JsonIgnore] public double CenterLat => 0.5 * (MinLat + MaxLat);

    /// <summary>East-west extent at the central latitude, km.</summary>
    [JsonIgnore]
    public double WidthKm => (MaxLon - MinLon) * GeoMath.Deg2Rad * GeoMath.EarthRadiusKm * Math.Cos(CenterLat * GeoMath.Deg2Rad);

    /// <summary>North-south extent, km.</summary>
    [JsonIgnore]
    public double HeightKm => (MaxLat - MinLat) * GeoMath.Deg2Rad * GeoMath.EarthRadiusKm;

    /// <summary>Whether a point lies in the box (its longitude taken in the box's continuous range).</summary>
    public bool Contains(double lon, double lat)
    {
        lon = GeoMath.UnwrapLon(lon, CenterLon);
        return lon >= MinLon && lon <= MaxLon && lat >= MinLat && lat <= MaxLat;
    }
}

/// <summary>
/// The JSON conventions of every file TomoSTAR reads or writes. They are the conventions of QUIVER
/// (enums as strings, NaN and infinities allowed, property names matched without regard to case), so
/// grid definitions, 1-D models, settings and volume headers pass between the two programs unchanged.
/// Comments and trailing commas are accepted in input files, which makes configuration files
/// easier to annotate by hand.
/// </summary>
public static class TomoJson
{
    /// <summary>Options for reading and writing.</summary>
    public static readonly JsonSerializerOptions Options = new()
    {
        WriteIndented = true,
        NumberHandling = JsonNumberHandling.AllowNamedFloatingPointLiterals,
        Converters = { new JsonStringEnumConverter() },
        PropertyNameCaseInsensitive = true,
        ReadCommentHandling = JsonCommentHandling.Skip,
        AllowTrailingCommas = true,
        // Explicit, so serialisation works even where reflection defaults are switched off
        // (trimmed or single-file builds).
        TypeInfoResolver = new DefaultJsonTypeInfoResolver()
    };
}
