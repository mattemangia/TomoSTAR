// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using System.Buffers.Binary;
using System.Text;

namespace TomoStar.Core.IO;

/// <summary>
/// miniSEED 2.4 reader (FDSN SEED Reference Manual v2.4, 2012): the fixed 48-byte data header,
/// blockette 1000 (encoding, byte order, record length) and blockette 100 (exact sample rate), and
/// the encodings found in FDSN data centres: INT16 (1), INT24 (2), INT32 (3), FLOAT32 (4),
/// FLOAT64 (5), STEIM-1 (10) and STEIM-2 (11).
///
/// Records are grouped by channel and joined where they are contiguous; a gap or an unmergeable
/// overlap starts a new segment, as libmseed and ObsPy do, so no data after a gap is lost.
/// Record lengths may vary within a file (each record's own blockette 1000 is read).
/// </summary>
public static class MiniSeed
{
    private const int FixedHeader = 48;

    public static List<Trace> ReadFile(string path)
    {
        using var fs = File.OpenRead(path);
        return Read(fs);
    }

    public static List<Trace> Read(Stream stream)
    {
        using var ms = new MemoryStream();
        stream.CopyTo(ms);
        return Read(ms.GetBuffer().AsSpan(0, (int)ms.Length));
    }

    public static List<Trace> Read(ReadOnlySpan<byte> data)
    {
        var records = new List<Record>();
        var offset = 0;
        while (offset + FixedHeader <= data.Length)
        {
            var length = RecordLength(data[offset..]);
            if (length < FixedHeader || offset + length > data.Length)
            {
                // Not a data record (or garbage): resynchronise on the next 256-byte boundary.
                offset += 256;
                continue;
            }
            try
            {
                var r = ParseRecord(data.Slice(offset, length));
                if (r != null && r.Samples.Length > 0) records.Add(r);
            }
            catch (Exception ex) when (ex is NotSupportedException or IndexOutOfRangeException or ArgumentOutOfRangeException)
            {
                // A malformed or unsupported record is skipped; the rest of the stream still reads.
            }
            offset += length;
        }
        return Stitch(records);
    }

    private static int RecordLength(ReadOnlySpan<byte> rec)
    {
        var type = (char)rec[6];
        if (type is not ('D' or 'R' or 'Q' or 'M')) return -1;
        int blk = BinaryPrimitives.ReadUInt16BigEndian(rec[46..]);
        for (var guard = 0; blk != 0 && blk + 8 <= rec.Length && guard < 16; guard++)
        {
            int kind = BinaryPrimitives.ReadUInt16BigEndian(rec[blk..]);
            int next = BinaryPrimitives.ReadUInt16BigEndian(rec[(blk + 2)..]);
            if (kind == 1000) return 1 << rec[blk + 6];
            if (next <= blk) break;
            blk = next;
        }
        return -1;
    }

    private sealed class Record
    {
        public string Net = "", Sta = "", Loc = "", Cha = "";
        public DateTime Start;
        public double Rate;
        public float[] Samples = [];
        public string Key => $"{Net}.{Sta}.{Loc}.{Cha}";
    }

    private static Record? ParseRecord(ReadOnlySpan<byte> rec)
    {
        var r = new Record
        {
            Sta = Ascii(rec, 8, 5), Loc = Ascii(rec, 13, 2), Cha = Ascii(rec, 15, 3), Net = Ascii(rec, 18, 2)
        };
        int year = BinaryPrimitives.ReadUInt16BigEndian(rec[20..]);
        int jday = BinaryPrimitives.ReadUInt16BigEndian(rec[22..]);
        int hour = rec[24], minute = rec[25], second = rec[26];
        int frac = BinaryPrimitives.ReadUInt16BigEndian(rec[28..]);
        int nSamples = BinaryPrimitives.ReadUInt16BigEndian(rec[30..]);
        short factor = BinaryPrimitives.ReadInt16BigEndian(rec[32..]);
        short mult = BinaryPrimitives.ReadInt16BigEndian(rec[34..]);
        byte activity = rec[36];
        int nBlockettes = rec[39];
        int timeCorrection = BinaryPrimitives.ReadInt32BigEndian(rec[40..]);
        int dataOffset = BinaryPrimitives.ReadUInt16BigEndian(rec[44..]);
        int blk = BinaryPrimitives.ReadUInt16BigEndian(rec[46..]);

        if (year < 1 || year > 9999 || jday < 1 || jday > 366) return null;
        r.Rate = SampleRate(factor, mult);
        var encoding = 11;
        var littleEndian = false;
        var microseconds = 0;
        for (var b = 0; b < nBlockettes && blk != 0 && blk + 4 <= rec.Length; b++)
        {
            int kind = BinaryPrimitives.ReadUInt16BigEndian(rec[blk..]);
            int next = BinaryPrimitives.ReadUInt16BigEndian(rec[(blk + 2)..]);
            switch (kind)
            {
                case 1000:
                    encoding = rec[blk + 4];
                    littleEndian = rec[blk + 5] == 0;
                    break;
                case 100:
                    r.Rate = BinaryPrimitives.ReadSingleBigEndian(rec[(blk + 4)..]);
                    break;
                case 1001:
                    microseconds = (sbyte)rec[blk + 5];
                    break;
            }
            if (next <= blk) break;
            blk = next;
        }
        if (!(r.Rate > 0) || nSamples == 0) return null;

        var start = new DateTime(year, 1, 1, 0, 0, 0, DateTimeKind.Utc)
            .AddDays(jday - 1).AddHours(hour).AddMinutes(minute).AddSeconds(second)
            .AddTicks(frac * 1000L + microseconds * 10L);
        // Bit 1 of the activity flags: the time correction has already been applied.
        if ((activity & 0x02) == 0 && timeCorrection != 0) start = start.AddTicks(timeCorrection * 1000L);
        r.Start = start;
        r.Samples = Decode(rec, dataOffset, nSamples, encoding, littleEndian);
        return r;
    }

    private static string Ascii(ReadOnlySpan<byte> rec, int at, int n) => Encoding.ASCII.GetString(rec.Slice(at, n)).Trim();

    private static double SampleRate(short f, short m)
    {
        if (f > 0 && m > 0) return (double)f * m;
        if (f > 0 && m < 0) return -(double)f / m;
        if (f < 0 && m > 0) return -(double)m / f;
        if (f < 0 && m < 0) return 1.0 / ((double)f * m);
        return 0;
    }

    private static float[] Decode(ReadOnlySpan<byte> rec, int off, int n, int encoding, bool le)
    {
        var x = new float[n];
        switch (encoding)
        {
            case 1:
                for (var i = 0; i < n; i++) x[i] = le ? BinaryPrimitives.ReadInt16LittleEndian(rec[(off + 2 * i)..]) : BinaryPrimitives.ReadInt16BigEndian(rec[(off + 2 * i)..]);
                return x;
            case 2:
                for (var i = 0; i < n; i++)
                {
                    var a = off + 3 * i;
                    var v = le ? rec[a] | (rec[a + 1] << 8) | (rec[a + 2] << 16) : (rec[a] << 16) | (rec[a + 1] << 8) | rec[a + 2];
                    if ((v & 0x800000) != 0) v |= unchecked((int)0xFF000000);
                    x[i] = v;
                }
                return x;
            case 3:
                for (var i = 0; i < n; i++) x[i] = le ? BinaryPrimitives.ReadInt32LittleEndian(rec[(off + 4 * i)..]) : BinaryPrimitives.ReadInt32BigEndian(rec[(off + 4 * i)..]);
                return x;
            case 4:
                for (var i = 0; i < n; i++) x[i] = le ? BinaryPrimitives.ReadSingleLittleEndian(rec[(off + 4 * i)..]) : BinaryPrimitives.ReadSingleBigEndian(rec[(off + 4 * i)..]);
                return x;
            case 5:
                for (var i = 0; i < n; i++) x[i] = (float)(le ? BinaryPrimitives.ReadDoubleLittleEndian(rec[(off + 8 * i)..]) : BinaryPrimitives.ReadDoubleBigEndian(rec[(off + 8 * i)..]));
                return x;
            case 10:
            case 11:
                return DecodeSteim(rec, off, n, encoding == 11);
            default:
                throw new NotSupportedException($"miniSEED encoding {encoding} is not supported.");
        }
    }

    /// <summary>
    /// STEIM-1/2 difference decompression. A frame is 16 big-endian 32-bit words; word 0 holds a
    /// 2-bit code per word; in frame 0, words 1 and 2 are the forward (X0) and reverse (Xn)
    /// integration constants. The first difference refers to the previous record's last sample and
    /// is discarded: the first sample is X0 by definition (as in libmseed).
    /// </summary>
    private static float[] DecodeSteim(ReadOnlySpan<byte> rec, int off, int n, bool steim2)
    {
        var output = new float[n];
        var count = 0;
        var previous = 0;
        var skipFirstDifference = false;
        var first = true;

        void Emit(int d)
        {
            if (skipFirstDifference) { skipFirstDifference = false; return; }
            if (count >= n) return;
            previous += d;
            output[count++] = previous;
        }

        for (var frame = off; frame + 64 <= rec.Length && count < n; frame += 64)
        {
            var nibbles = BinaryPrimitives.ReadUInt32BigEndian(rec[frame..]);
            var w0 = 1;
            if (first)
            {
                var x0 = BinaryPrimitives.ReadInt32BigEndian(rec[(frame + 4)..]);
                previous = x0;
                output[count++] = x0;
                skipFirstDifference = true;
                w0 = 3;
                first = false;
            }
            for (var w = w0; w < 16 && count < n; w++)
            {
                var code = (int)((nibbles >> (30 - 2 * w)) & 3);
                var at = frame + 4 * w;
                if (code == 0) continue;
                if (code == 1)
                {
                    for (var i = 0; i < 4; i++) Emit((sbyte)rec[at + i]);
                    continue;
                }
                if (!steim2)
                {
                    if (code == 2)
                        for (var i = 0; i < 2; i++) Emit(BinaryPrimitives.ReadInt16BigEndian(rec[(at + 2 * i)..]));
                    else
                        Emit(BinaryPrimitives.ReadInt32BigEndian(rec[at..]));
                    continue;
                }
                var word = BinaryPrimitives.ReadUInt32BigEndian(rec[at..]);
                var dnib = (int)(word >> 30);
                if (code == 2)
                {
                    switch (dnib)
                    {
                        case 1: Emit(Sign(word & 0x3FFFFFFF, 30)); break;
                        case 2: for (var i = 0; i < 2; i++) Emit(Sign((word >> (15 * (1 - i))) & 0x7FFF, 15)); break;
                        case 3: for (var i = 0; i < 3; i++) Emit(Sign((word >> (10 * (2 - i))) & 0x3FF, 10)); break;
                    }
                }
                else
                {
                    switch (dnib)
                    {
                        case 0: for (var i = 0; i < 5; i++) Emit(Sign((word >> (6 * (4 - i))) & 0x3F, 6)); break;
                        case 1: for (var i = 0; i < 6; i++) Emit(Sign((word >> (5 * (5 - i))) & 0x1F, 5)); break;
                        case 2: for (var i = 0; i < 7; i++) Emit(Sign((word >> (4 * (6 - i))) & 0xF, 4)); break;
                    }
                }
            }
        }
        return count == n ? output : output[..count];
    }

    private static int Sign(uint value, int bits)
    {
        var shift = 32 - bits;
        return (int)(value << shift) >> shift;
    }

    private static List<Trace> Stitch(List<Record> records)
    {
        var result = new List<Trace>();
        foreach (var group in records.GroupBy(r => r.Key))
        {
            var ordered = group.OrderBy(r => r.Start).ToList();
            var parts = new List<float[]>();
            Record? head = null;
            long total = 0;

            void Flush()
            {
                if (head == null) return;
                var data = new float[total];
                var o = 0;
                foreach (var p in parts) { p.CopyTo(data, o); o += p.Length; }
                result.Add(new Trace
                {
                    Network = head.Net, Station = head.Sta, Location = head.Loc, Channel = head.Cha,
                    StartTime = head.Start, SampleRate = head.Rate, Data = data
                });
            }

            foreach (var r in ordered)
            {
                if (head == null)
                {
                    head = r; parts = [r.Samples]; total = r.Samples.Length;
                    continue;
                }
                var dt = 1.0 / head.Rate;
                var expected = head.Start.AddSeconds(total * dt);
                var gap = (r.Start - expected).TotalSeconds;
                if (Math.Abs(gap) < 0.5 * dt && Math.Abs(r.Rate - head.Rate) < 1e-6 * head.Rate)
                {
                    parts.Add(r.Samples);
                    total += r.Samples.Length;
                    continue;
                }
                if (gap < 0 && r.Start.AddSeconds(r.Samples.Length * dt) <= expected) continue; // duplicate
                Flush();
                head = r; parts = [r.Samples]; total = r.Samples.Length;
            }
            Flush();
        }
        return result;
    }

    // ---- Writer (INT32 / FLOAT32, uncompressed) ----

    /// <summary>
    /// Writes traces as 4096-byte FLOAT32 records with blockette 1000. Used to store processed or
    /// imported waveforms inside a project; any miniSEED reader can open them.
    /// </summary>
    public static void Write(string path, IEnumerable<Trace> traces)
    {
        const int recLen = 4096;
        const int dataStart = 64;
        const int perRecord = (recLen - dataStart) / 4;
        using var fs = File.Create(path);
        var seq = 1;
        var rec = new byte[recLen];
        foreach (var t in traces)
        {
            for (var start = 0; start < t.Data.Length; start += perRecord)
            {
                Array.Clear(rec);
                var n = Math.Min(perRecord, t.Data.Length - start);
                var time = t.TimeOf(start);
                Encoding.ASCII.GetBytes(seq++.ToString("D6")[^6..]).CopyTo(rec, 0);
                rec[6] = (byte)'D';
                rec[7] = (byte)' ';
                Pad(t.Station, 5).CopyTo(rec, 8);
                Pad(t.Location, 2).CopyTo(rec, 13);
                Pad(t.Channel, 3).CopyTo(rec, 15);
                Pad(t.Network, 2).CopyTo(rec, 18);
                BinaryPrimitives.WriteUInt16BigEndian(rec.AsSpan(20), (ushort)time.Year);
                BinaryPrimitives.WriteUInt16BigEndian(rec.AsSpan(22), (ushort)time.DayOfYear);
                rec[24] = (byte)time.Hour;
                rec[25] = (byte)time.Minute;
                rec[26] = (byte)time.Second;
                BinaryPrimitives.WriteUInt16BigEndian(rec.AsSpan(28), (ushort)(time.Ticks % TimeSpan.TicksPerSecond / 1000));
                BinaryPrimitives.WriteUInt16BigEndian(rec.AsSpan(30), (ushort)n);
                var (factor, mult) = RateFactors(t.SampleRate);
                BinaryPrimitives.WriteInt16BigEndian(rec.AsSpan(32), factor);
                BinaryPrimitives.WriteInt16BigEndian(rec.AsSpan(34), mult);
                rec[39] = 2;
                BinaryPrimitives.WriteUInt16BigEndian(rec.AsSpan(44), dataStart);
                BinaryPrimitives.WriteUInt16BigEndian(rec.AsSpan(46), 48);
                // Blockette 1000 at 48, blockette 100 at 56.
                BinaryPrimitives.WriteUInt16BigEndian(rec.AsSpan(48), 1000);
                BinaryPrimitives.WriteUInt16BigEndian(rec.AsSpan(50), 56);
                rec[52] = 4;   // FLOAT32
                rec[53] = 1;   // big-endian
                rec[54] = 12;  // 2^12 = 4096
                BinaryPrimitives.WriteUInt16BigEndian(rec.AsSpan(56), 100);
                BinaryPrimitives.WriteUInt16BigEndian(rec.AsSpan(58), 0);
                BinaryPrimitives.WriteSingleBigEndian(rec.AsSpan(60), (float)t.SampleRate);
                for (var i = 0; i < n; i++)
                    BinaryPrimitives.WriteSingleBigEndian(rec.AsSpan(dataStart + 4 * i), t.Data[start + i]);
                fs.Write(rec);
            }
        }
    }

    private static byte[] Pad(string s, int n) => Encoding.ASCII.GetBytes((s ?? "").PadRight(n)[..n]);

    private static (short Factor, short Mult) RateFactors(double rate)
    {
        if (rate >= 1 && Math.Abs(rate - Math.Round(rate)) < 1e-9 && rate < 32767) return ((short)Math.Round(rate), 1);
        if (rate < 1 && rate > 0 && Math.Abs(1 / rate - Math.Round(1 / rate)) < 1e-9) return ((short)-Math.Round(1 / rate), 1);
        return ((short)Math.Round(rate * 100), -100); // blockette 100 carries the exact value
    }
}
