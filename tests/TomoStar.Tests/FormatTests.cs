// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using Microsoft.Data.Sqlite;
using TomoStar.Core.Geo;
using TomoStar.Core.IO;
using TomoStar.Core.Model;

namespace TomoStar.Tests;

/// <summary>Every file TomoSTAR reads or writes, by round trip or against a hand-made file.</summary>
public class FormatTests : IDisposable
{
    private readonly string _dir = Directory.CreateTempSubdirectory("tomostar-test-").FullName;

    public void Dispose()
    {
        try { Directory.Delete(_dir, true); } catch (IOException) { }
    }

    private string P(string name) => Path.Combine(_dir, name);

    [Fact]
    public void CsvTablesRoundTrip()
    {
        File.WriteAllText(P("stations.csv"), "network;station;lon;lat;elevation_m\nIV;CAMP;13.40;42.53;1250\nIV;ARRO;12.77;42.58;253\n");
        File.WriteAllText(P("events.csv"), "# a comment\nid,time,lon,lat,depth_km,magnitude\ne1,2016-10-30T06:40:17.32Z,13.11,42.83,9.2,6.5\n");
        File.WriteAllText(P("picks.csv"), "event,station,phase,travel_time_s,sigma_s,origin\ne1,IV.CAMP,Pg,5.25,0.05,manual\ne1,IV.ARRO,Sn,12.5,,automatic\ne1,IV.ARRO,PmP,9,,\n");
        var c = CatalogueReader.ReadCsv(P("stations.csv"), P("events.csv"), P("picks.csv"));
        Assert.Equal(2, c.Stations.Count);
        Assert.Equal("IV.CAMP", c.Stations[0].Id);
        Assert.Equal(-1.25, c.Stations[0].DepthKm, 9);
        var e = Assert.Single(c.Events);
        Assert.Equal(new DateTime(2016, 10, 30, 6, 40, 17, 320, DateTimeKind.Utc), e.OriginTime);
        // PmP is a later phase, not a first arrival: skipped with a warning.
        Assert.Equal(2, e.Picks.Count);
        Assert.Contains(c.Warnings, w => w.Contains("other than P and S"));
        Assert.Equal(PickOrigin.Manual, e.Picks[0].Origin);
        Assert.Equal(5.25, (e.Picks[0].Time - e.OriginTime).TotalSeconds, 6);
        Assert.Equal(0.2, e.Picks[1].Sigma, 9); // the S default

        CatalogueWriter.WriteStations(P("s2.csv"), c.Stations);
        CatalogueWriter.WriteEvents(P("e2.csv"), c.Events);
        CatalogueWriter.WritePicks(P("p2.csv"), c.Events);
        var d = CatalogueReader.ReadCsv(P("s2.csv"), P("e2.csv"), P("p2.csv"));
        Assert.Equal(c.Stations.Select(s => (s.Id, s.Lon, s.Lat, s.ElevationM)), d.Stations.Select(s => (s.Id, s.Lon, s.Lat, s.ElevationM)));
        Assert.Equal(c.Events[0].Picks.Select(p => (p.StationId, p.Phase, p.Time, p.Origin)), d.Events[0].Picks.Select(p => (p.StationId, p.Phase, p.Time, p.Origin)));
    }

    [Fact]
    public void BadInputIsReportedWithFileAndLine()
    {
        File.WriteAllText(P("stations.csv"), "id,lon,lat\nA,13.1,oops\n");
        var ex = Assert.Throws<InvalidDataException>(() => CatalogueReader.ReadCsv(P("stations.csv"), null, null));
        Assert.Contains("line 2", ex.Message);
    }

    [Fact]
    public void VolumesRoundTripInEveryFormat()
    {
        var def = new GridDefinition { MinLon = 12, MaxLon = 13, MinLat = 42, MaxLat = 43, MinDepthKm = -1, MaxDepthKm = 20, Nx = 4, Ny = 3, Nz = 5 };
        var values = Enumerable.Range(0, (int)def.Count).Select(i => i * 0.5).ToArray();
        values[7] = double.NaN;
        var header = new VolumeHeader { Name = "Vp (run; test)", Quantity = "Vp", Units = "km/s", Grid = def };
        VolumeFile.Write(P("v.qvol"), header, values);
        using (var r = new VolumeReader(P("v.qvol")))
        {
            var back = r.ReadAll();
            Assert.True(double.IsNaN(back[7]));
            Assert.Equal(values[8], back[8], 6);
            Assert.Equal("Vp (run; test)", r.Header.Name);
        }
        VolumeExport.WriteCsv(P("v.csv"), header, values);
        var lines = File.ReadAllLines(P("v.csv"));
        // The header QUIVER reads back: separators inside the name are percent-encoded.
        Assert.Equal("# QUIVER volume: Vp (run%3B test); quantity=Vp; units=km/s", lines[0]);
        Assert.Equal("lon,lat,depth_km,value", lines[1]);
        Assert.Equal(def.Count + 2, lines.Length);
        VolumeExport.WriteVtk(P("v.vtk"), header, values);
        using var vtk = new StreamReader(P("v.vtk"));
        Assert.Equal("# vtk DataFile Version 3.0", vtk.ReadLine());
        Assert.StartsWith("QUIVER volume;lon=12,13,4;lat=42,43,3;dep=-1,20,5;q=Vp", vtk.ReadLine());
    }

    [Fact]
    public void MiniSeedWriterIsReadBack()
    {
        var t = new Trace
        {
            Network = "SY", Station = "S001", Location = "", Channel = "HHZ", SampleRate = 100,
            StartTime = new DateTime(2026, 3, 4, 5, 6, 7, 890, DateTimeKind.Utc),
            Data = Enumerable.Range(0, 1000).Select(i => (float)Math.Sin(i * 0.05)).ToArray()
        };
        MiniSeedWriter.Write(P("a.mseed"), [t]);
        var back = Assert.Single(MiniSeed.ReadFile(P("a.mseed")));
        Assert.Equal("SY.S001..HHZ", back.Nslc);
        Assert.Equal(t.StartTime, back.StartTime);
        Assert.Equal(100, back.SampleRate);
        Assert.Equal(t.Data, back.Data);
    }

    [Fact]
    public void ModelsAreReadFromTheLibraryAndFromFiles()
    {
        Assert.Equal("ak135", CatalogueReader.ReadModel1D("ak135").Name);
        File.WriteAllText(P("m.txt"), "# depth vp vs\n0 5.0 2.9\n10 6.0 3.5\n10 6.5 3.7\n40 8.0 4.5\n");
        var m = CatalogueReader.ReadModel1D(P("m.txt"));
        Assert.Equal(5.5, m.Vp(5), 9);
        Assert.Equal(6.0, m.Vp(10), 9); // on a discontinuity: the upper value
        Assert.Throws<FileNotFoundException>(() => CatalogueReader.ReadModel1D("no-such-model"));
    }

    /// <summary>
    /// A project folder laid out as QUIVER writes it (project.json, stations.json, catalog.sqlite with
    /// its tables and pick flags), read without QUIVER.
    /// </summary>
    [Fact]
    public void QuiverProjectsAreReadInPlace()
    {
        var project = P("study.quiver");
        Directory.CreateDirectory(project);
        File.WriteAllText(Path.Combine(project, "project.json"), """
            {
              "SchemaVersion": 1, "Name": "study",
              "Grid": { "MinLon": 13, "MaxLon": 13.5, "MinLat": 42.5, "MaxLat": 43, "MinDepthKm": -2, "MaxDepthKm": 20, "Nx": 6, "Ny": 6, "Nz": 12 },
              "StartingModel": { "Name": "mine", "Nodes": [ { "DepthKm": 0, "Vp": 5.5, "Vs": 3.1 }, { "DepthKm": 30, "Vp": 6.8, "Vs": 3.9 } ] },
              "OutsideData": "WhenRaysCross",
              "Tomography": { "Smoothing": 12, "SomeFutureSetting": 1 },
              "Volumes": []
            }
            """);
        File.WriteAllText(Path.Combine(project, "stations.json"), """
            [ { "Network": "IV", "Code": "AAA", "Latitude": 42.7, "Longitude": 13.2, "ElevationM": 500, "Channels": [] } ]
            """);
        using (var c = new SqliteConnection($"Data Source={Path.Combine(project, "catalog.sqlite")};Pooling=False"))
        {
            c.Open();
            using var cmd = c.CreateCommand();
            cmd.CommandText = """
                CREATE TABLE events(rid INTEGER PRIMARY KEY, ord INTEGER NOT NULL, id TEXT NOT NULL, origin_time INTEGER NOT NULL, lat REAL, lon REAL, depth_km REAL, magnitude REAL, print INTEGER NOT NULL, data BLOB NOT NULL);
                CREATE TABLE channels(cid INTEGER PRIMARY KEY, network TEXT NOT NULL, station TEXT NOT NULL, location TEXT NOT NULL, channel TEXT NOT NULL);
                CREATE TABLE sources(sid INTEGER PRIMARY KEY, text TEXT NOT NULL UNIQUE);
                CREATE TABLE picks(event_rid INTEGER NOT NULL, seq INTEGER NOT NULL, cid INTEGER NOT NULL, time INTEGER NOT NULL, uncertainty REAL NOT NULL, residual REAL, flags INTEGER NOT NULL, sid INTEGER NOT NULL, pick_id TEXT, PRIMARY KEY(event_rid, seq));
                CREATE TABLE tstar(rid INTEGER PRIMARY KEY, event_id TEXT NOT NULL, station_id TEXT NOT NULL, data BLOB NOT NULL);
                INSERT INTO channels VALUES (1, 'IV', 'AAA', '', 'HHZ');
                INSERT INTO sources VALUES (1, 'INGV');
                """;
            cmd.ExecuteNonQuery();
            var origin = new DateTime(2026, 1, 1, 0, 0, 0, DateTimeKind.Utc);
            cmd.CommandText = "INSERT INTO events VALUES (1, 1024, 'ev1', $t, 42.7, 13.2, 8, 2.5, 0, $d)";
            cmd.Parameters.AddWithValue("$t", origin.Ticks);
            cmd.Parameters.AddWithValue("$d", System.Text.Encoding.UTF8.GetBytes("""
                {"Id":"ev1","Catalog":{"OriginTime":"2026-01-01T00:00:00Z","Latitude":42.7,"Longitude":13.2,"DepthKm":8},
                 "Relocated":{"OriginTime":"2026-01-01T00:00:00.5Z","Latitude":42.71,"Longitude":13.21,"DepthKm":9},"Magnitude":2.5}
                """));
            cmd.ExecuteNonQuery();
            cmd.Parameters.Clear();
            // Flags: phase in bits 0-3 (S = 1), origin in 8-11 (manual = 2), disabled in bit 15.
            cmd.CommandText = $$"""
                INSERT INTO picks VALUES (1, 0, 1, {{origin.AddSeconds(3).Ticks}}, 0.05, NULL, {{0 | (2 << 8)}}, 1, NULL);
                INSERT INTO picks VALUES (1, 1, 1, {{origin.AddSeconds(5).Ticks}}, 0.1, NULL, {{1 | (1 << 8)}}, 1, NULL);
                INSERT INTO picks VALUES (1, 2, 1, {{origin.AddSeconds(6).Ticks}}, 0.1, NULL, {{1 | (1 << 15)}}, 1, NULL);
                INSERT INTO tstar VALUES (1, 'ev1', 'IV.AAA', '{"EventId":"ev1","StationId":"IV.AAA","Phase":"P","TStar":0.021,"Uncertainty":0.003}');
                """;
            cmd.ExecuteNonQuery();
        }
        Assert.True(CatalogueReader.IsQuiverProject(project));
        var info = CatalogueReader.ReadQuiverProjectInfo(project);
        Assert.Equal(12, info.Grid!.Nz);
        Assert.Equal("mine", info.StartingModel!.Name);
        Assert.Equal(12, info.Tomography!.Smoothing);
        Assert.Equal(OutsideData.WhenRaysCross, info.OutsideData);

        var cat = CatalogueReader.ReadQuiverProject(project);
        var s = Assert.Single(cat.Stations);
        Assert.Equal("IV.AAA", s.Id);
        var e = Assert.Single(cat.Events);
        Assert.Equal(9, e.DepthKm); // the relocation, not the catalogue location
        Assert.Equal(3, e.Picks.Count);
        Assert.Equal(PickOrigin.Manual, e.Picks[0].Origin);
        Assert.Equal(Phase.S, e.Picks[1].Phase);
        Assert.Equal(PickOrigin.Automatic, e.Picks[1].Origin);
        Assert.True(e.Picks[2].Disabled);
        Assert.Equal(0.021, Assert.Single(cat.TStar).TStar, 9);
        Assert.Equal(8, CatalogueReader.ReadQuiverProject(project, catalogLocations: true).Events[0].DepthKm);
    }
}
