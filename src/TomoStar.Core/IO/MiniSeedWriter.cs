using System.Buffers.Binary;
using System.Text;

namespace TomoStar.Core.IO;

/// <summary>
/// A minimal miniSEED 2.4 writer: 512-byte data records, FLOAT32 samples (encoding 4), big-endian,
/// one blockette 1000 per record (FDSN SEED Reference Manual v2.4, 2012). It exists so that the
/// synthetic data sets of TomoSTAR are complete (waveforms included) and can be read back by any
/// miniSEED reader (ObsPy, libmseed, QUIVER); it is not meant as a general archive writer.
/// </summary>
public static class MiniSeedWriter
{
    private const int RecordLength = 512;
    private const int DataOffset = 64;
    private const int SamplesPerRecord = (RecordLength - DataOffset) / 4;

    /// <summary>Writes the traces to one file, each trace as a run of consecutive records.</summary>
    public static void Write(string path, IEnumerable<Trace> traces)
    {
        var dir = Path.GetDirectoryName(Path.GetFullPath(path));
        if (dir != null) Directory.CreateDirectory(dir);
        using var fs = new FileStream(path, FileMode.Create, FileAccess.Write, FileShare.None, 1 << 16);
        var sequence = 1;
        var rec = new byte[RecordLength];
        foreach (var t in traces)
        {
            var (factor, multiplier) = RateFactors(t.SampleRate);
            for (var start = 0; start < t.Data.Length; start += SamplesPerRecord)
            {
                var n = Math.Min(SamplesPerRecord, t.Data.Length - start);
                Array.Clear(rec);
                Ascii(rec, 0, (sequence++ % 1000000).ToString("D6"), 6);
                rec[6] = (byte)'D';
                rec[7] = (byte)' ';
                Ascii(rec, 8, t.Station, 5);
                Ascii(rec, 13, t.Location, 2);
                Ascii(rec, 15, t.Channel, 3);
                Ascii(rec, 18, t.Network, 2);
                // Record start time as BTIME, rounded to 0.1 ms.
                var time = t.StartTime.AddTicks((long)Math.Round(start / t.SampleRate * TimeSpan.TicksPerSecond));
                BinaryPrimitives.WriteUInt16BigEndian(rec.AsSpan(20), (ushort)time.Year);
                BinaryPrimitives.WriteUInt16BigEndian(rec.AsSpan(22), (ushort)time.DayOfYear);
                rec[24] = (byte)time.Hour;
                rec[25] = (byte)time.Minute;
                rec[26] = (byte)time.Second;
                BinaryPrimitives.WriteUInt16BigEndian(rec.AsSpan(28), (ushort)(time.Ticks % TimeSpan.TicksPerSecond / 1000));
                BinaryPrimitives.WriteUInt16BigEndian(rec.AsSpan(30), (ushort)n);
                BinaryPrimitives.WriteInt16BigEndian(rec.AsSpan(32), factor);
                BinaryPrimitives.WriteInt16BigEndian(rec.AsSpan(34), multiplier);
                rec[39] = 1; // one blockette
                BinaryPrimitives.WriteUInt16BigEndian(rec.AsSpan(44), DataOffset);
                BinaryPrimitives.WriteUInt16BigEndian(rec.AsSpan(46), 48);
                // Blockette 1000: encoding FLOAT32, big-endian, record length 2^9.
                BinaryPrimitives.WriteUInt16BigEndian(rec.AsSpan(48), 1000);
                BinaryPrimitives.WriteUInt16BigEndian(rec.AsSpan(50), 0);
                rec[52] = 4;
                rec[53] = 1;
                rec[54] = 9;
                for (var i = 0; i < n; i++) BinaryPrimitives.WriteSingleBigEndian(rec.AsSpan(DataOffset + 4 * i), t.Data[start + i]);
                fs.Write(rec);
            }
        }
    }

    /// <summary>Sample-rate factor and multiplier: integer rates directly, others as the inverse of a period.</summary>
    private static (short Factor, short Multiplier) RateFactors(double rate)
    {
        if (Math.Abs(rate - Math.Round(rate)) < 1e-9 && rate >= 1 && rate <= short.MaxValue) return ((short)Math.Round(rate), 1);
        var period = 1 / rate;
        if (Math.Abs(period - Math.Round(period)) < 1e-9 && period <= short.MaxValue) return ((short)-Math.Round(period), 1);
        throw new ArgumentException($"Sample rate {rate} Hz cannot be written as an integer rate or period.");
    }

    private static void Ascii(byte[] rec, int at, string text, int width)
    {
        var s = (text ?? "").PadRight(width)[..width];
        Encoding.ASCII.GetBytes(s, 0, width, rec, at);
    }
}
