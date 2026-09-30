namespace TomoStar.Core.Geo;

/// <summary>A geocentric Cartesian position or vector, km.</summary>
public readonly record struct Vec3(double X, double Y, double Z)
{
    public static Vec3 operator +(Vec3 a, Vec3 b) => new(a.X + b.X, a.Y + b.Y, a.Z + b.Z);
    public static Vec3 operator -(Vec3 a, Vec3 b) => new(a.X - b.X, a.Y - b.Y, a.Z - b.Z);
    public static Vec3 operator *(Vec3 a, double s) => new(a.X * s, a.Y * s, a.Z * s);
    public static Vec3 operator *(double s, Vec3 a) => new(a.X * s, a.Y * s, a.Z * s);
    public static Vec3 operator /(Vec3 a, double s) => new(a.X / s, a.Y / s, a.Z / s);
    public static Vec3 operator -(Vec3 a) => new(-a.X, -a.Y, -a.Z);

    public double Dot(Vec3 b) => X * b.X + Y * b.Y + Z * b.Z;
    public Vec3 Cross(Vec3 b) => new(Y * b.Z - Z * b.Y, Z * b.X - X * b.Z, X * b.Y - Y * b.X);
    public double Length => Math.Sqrt(X * X + Y * Y + Z * Z);

    public Vec3 Normalized()
    {
        var l = Length;
        return l > 0 ? this / l : this;
    }
}

/// <summary>A geographic position: degrees and km below sea level (negative above it).</summary>
public readonly record struct GeoPoint(double Lon, double Lat, double DepthKm = 0);

/// <summary>
/// Spherical-Earth geometry. Every position in TomoSTAR goes through here, so a local survey and a
/// 1000 km deep model use one set of formulas: a local model is simply one where the curvature
/// terms are small, never one where they are dropped.
///
/// The Earth is a sphere of radius <see cref="EarthRadiusKm"/>. Geographic latitude is used as
/// geocentric latitude; the ellipticity correction (up to ~0.19° in latitude, ~21 km in radius)
/// is below the resolution of travel-time tomography at the scales TomoSTAR targets and is not
/// applied, as in the spherical tomography codes it follows (e.g. Rawlinson &amp; Sambridge 2004).
/// </summary>
public static class GeoMath
{
    /// <summary>Mean Earth radius, km (the ak135 value).</summary>
    public const double EarthRadiusKm = 6371.0;

    public const double Deg2Rad = Math.PI / 180.0;
    public const double Rad2Deg = 180.0 / Math.PI;

    /// <summary>Geocentric Cartesian position of a geographic point, km.</summary>
    public static Vec3 ToCartesian(double lonDeg, double latDeg, double depthKm)
    {
        var r = EarthRadiusKm - depthKm;
        var lat = latDeg * Deg2Rad;
        var lon = lonDeg * Deg2Rad;
        var cl = Math.Cos(lat);
        return new Vec3(r * cl * Math.Cos(lon), r * cl * Math.Sin(lon), r * Math.Sin(lat));
    }

    public static Vec3 ToCartesian(GeoPoint p) => ToCartesian(p.Lon, p.Lat, p.DepthKm);

    /// <summary>Inverse of <see cref="ToCartesian(double,double,double)"/>.</summary>
    public static GeoPoint FromCartesian(Vec3 v)
    {
        var r = v.Length;
        if (r <= 0) return new GeoPoint(0, 0, EarthRadiusKm);
        var lat = Math.Asin(Math.Clamp(v.Z / r, -1, 1)) * Rad2Deg;
        var lon = Math.Atan2(v.Y, v.X) * Rad2Deg;
        return new GeoPoint(lon, lat, EarthRadiusKm - r);
    }

    /// <summary>Local unit vectors (east, north, up) at a geographic point.</summary>
    public static (Vec3 East, Vec3 North, Vec3 Up) LocalFrame(double lonDeg, double latDeg)
    {
        var lat = latDeg * Deg2Rad;
        var lon = lonDeg * Deg2Rad;
        double sl = Math.Sin(lat), cl = Math.Cos(lat), so = Math.Sin(lon), co = Math.Cos(lon);
        var east = new Vec3(-so, co, 0);
        var north = new Vec3(-sl * co, -sl * so, cl);
        var up = new Vec3(cl * co, cl * so, sl);
        return (east, north, up);
    }

    /// <summary>Great-circle angular distance, radians (haversine form: exact at small angles).</summary>
    public static double AngularDistance(double lon1, double lat1, double lon2, double lat2)
    {
        var p1 = lat1 * Deg2Rad;
        var p2 = lat2 * Deg2Rad;
        var dp = p2 - p1;
        var dl = (lon2 - lon1) * Deg2Rad;
        var a = Math.Sin(dp / 2) * Math.Sin(dp / 2) + Math.Cos(p1) * Math.Cos(p2) * Math.Sin(dl / 2) * Math.Sin(dl / 2);
        return 2 * Math.Asin(Math.Min(1, Math.Sqrt(a)));
    }

    /// <summary>Great-circle distance on the surface, km.</summary>
    public static double SurfaceDistanceKm(double lon1, double lat1, double lon2, double lat2)
        => AngularDistance(lon1, lat1, lon2, lat2) * EarthRadiusKm;

    /// <summary>Initial azimuth from point 1 to point 2, degrees clockwise from north in [0, 360).</summary>
    public static double Azimuth(double lon1, double lat1, double lon2, double lat2)
    {
        var p1 = lat1 * Deg2Rad;
        var p2 = lat2 * Deg2Rad;
        var dl = (lon2 - lon1) * Deg2Rad;
        var y = Math.Sin(dl) * Math.Cos(p2);
        var x = Math.Cos(p1) * Math.Sin(p2) - Math.Sin(p1) * Math.Cos(p2) * Math.Cos(dl);
        var az = Math.Atan2(y, x) * Rad2Deg;
        return (az + 360) % 360;
    }

    /// <summary>Point at fraction <paramref name="f"/> of the great circle from 1 to 2 (slerp).</summary>
    public static (double Lon, double Lat) Intermediate(double lon1, double lat1, double lon2, double lat2, double f)
    {
        var a = ToCartesian(lon1, lat1, 0).Normalized();
        var b = ToCartesian(lon2, lat2, 0).Normalized();
        var omega = Math.Acos(Math.Clamp(a.Dot(b), -1, 1));
        if (omega < 1e-12) return (lon1, lat1);
        var s = Math.Sin(omega);
        var v = a * (Math.Sin((1 - f) * omega) / s) + b * (Math.Sin(f * omega) / s);
        var g = FromCartesian(v * EarthRadiusKm);
        return (g.Lon, g.Lat);
    }

    /// <summary>Destination point from a start, an azimuth (deg) and a surface distance (km).</summary>
    public static (double Lon, double Lat) Destination(double lon, double lat, double azimuthDeg, double distanceKm)
    {
        var d = distanceKm / EarthRadiusKm;
        var p1 = lat * Deg2Rad;
        var l1 = lon * Deg2Rad;
        var az = azimuthDeg * Deg2Rad;
        var p2 = Math.Asin(Math.Sin(p1) * Math.Cos(d) + Math.Cos(p1) * Math.Sin(d) * Math.Cos(az));
        var l2 = l1 + Math.Atan2(Math.Sin(az) * Math.Sin(d) * Math.Cos(p1), Math.Cos(d) - Math.Sin(p1) * Math.Sin(p2));
        return (NormalizeLon(l2 * Rad2Deg), p2 * Rad2Deg);
    }

    public static double NormalizeLon(double lon)
    {
        lon = (lon + 180) % 360;
        if (lon < 0) lon += 360;
        return lon - 180;
    }

    /// <summary>
    /// The longitude equal to <paramref name="lon"/> within 180° of <paramref name="reference"/>. A
    /// project across the antimeridian keeps one continuous range (e.g. 175° to 185°), and every
    /// longitude that comes in from outside, or from a Cartesian position, is brought into it.
    /// </summary>
    public static double UnwrapLon(double lon, double reference)
    {
        if (lon < reference - 180) lon += 360 * Math.Ceiling((reference - 180 - lon) / 360);
        else if (lon > reference + 180) lon -= 360 * Math.Ceiling((lon - reference - 180) / 360);
        return lon;
    }

    /// <summary>
    /// A continuous longitude range as one or two ranges inside [−180°, 180°], for services that
    /// take a box and do not accept one across the antimeridian.
    /// </summary>
    public static IReadOnlyList<(double Min, double Max)> LonRanges(double min, double max)
    {
        if (!double.IsFinite(min) || !double.IsFinite(max)) throw new ArgumentException("Longitude bounds must be finite.");
        if (max < min) max += 360 * Math.Ceiling((min - max) / 360);
        if (max - min >= 360) return [(-180, 180)];
        var a = NormalizeLon(min);
        var b = a + (max - min);
        return b <= 180 ? [(a, b)] : [(a, 180), (-180, b - 360)];
    }

    /// <summary>A longitude with its hemisphere: "70.6°W", "12.3°E".</summary>
    public static string FormatLon(double lon, string format = "0.####")
    {
        lon = NormalizeLon(lon);
        if (lon == -180) lon = 180;
        return Math.Abs(lon).ToString(format, System.Globalization.CultureInfo.InvariantCulture) + "°" + (lon < 0 ? "W" : "E");
    }

    /// <summary>A latitude with its hemisphere: "33.4°S", "42.1°N".</summary>
    public static string FormatLat(double lat, string format = "0.####") =>
        Math.Abs(lat).ToString(format, System.Globalization.CultureInfo.InvariantCulture) + "°" + (lat < 0 ? "S" : "N");

    /// <summary>"33.4°S 70.6°W".</summary>
    public static string FormatPosition(double lon, double lat, string format = "0.####") => FormatLat(lat, format) + " " + FormatLon(lon, format);

    /// <summary>A longitude range with its hemispheres: "175°E to 175°W".</summary>
    public static string FormatLonRange(double min, double max, string format = "0.###") => FormatLon(min, format) + " to " + FormatLon(max, format);

    /// <summary>A latitude range with its hemispheres.</summary>
    public static string FormatLatRange(double min, double max, string format = "0.###") => FormatLat(min, format) + " to " + FormatLat(max, format);

    /// <summary>
    /// Compass label for an azimuth, 16-point rose ("N", "NNE", ...). Used to orient slices
    /// ("SW → NE") so a section can never be read mirrored.
    /// </summary>
    public static string CompassPoint(double azimuthDeg)
    {
        string[] points = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE", "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"];
        var i = (int)Math.Round(((azimuthDeg % 360) + 360) % 360 / 22.5) % 16;
        return points[i];
    }
}
