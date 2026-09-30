using System.Collections.Concurrent;
using TomoStar.Core.Forward;
using TomoStar.Core.Geo;
using TomoStar.Core.Model;

namespace TomoStar.Core.Tomography;

/// <summary>
/// One traced arrival: the path integral of the trilinear basis along the ray (the row of G for
/// slowness), the computed travel time and the travel-time gradient at the source.
/// </summary>
public sealed record RayRow(int Arrival, int[] Nodes, float[] Length, double TCalc, double Ge, double Gn, double Gz, StoredRay? Ray)
{
    /// <summary>Travel time spent outside the inversion grid, in the fixed background, s (0 when the ray stays inside).</summary>
    public double OutsideTime { get; init; }
}

/// <summary>
/// Forward step shared by every path-integral inversion (travel time and t*): traces the ray of
/// every arrival in the current model and integrates it on the inversion grid.
/// </summary>
public sealed class RayKernels(SphericalGrid grid, ForwardDomain domain, RayMethod method, string workFolder, Action<string>? log = null)
{
    public RayKernels(SphericalGrid grid, GridDefinition forwardGrid, RayMethod method, string workFolder, Action<string>? log = null)
        : this(grid, new ForwardDomain(forwardGrid), method, workFolder, log) { }

    public List<RayRow> Compute(ObservationSet data, double[] sP, double[] sS, EikonalOpenCl? gpu, int keepRayPoints,
        IProgress<(double, string)>? progress = null, double progressBase = 0, CancellationToken ct = default, double progressSpan = 0)
    {
        // The pass fills [progressBase, progressBase + progressSpan]: the tables the first half, the rays the second.
        var tableSpan = method == RayMethod.FastMarching ? 0.5 * progressSpan : 0;
        var rows = new ConcurrentBag<RayRow>();
        TravelTimeTableSet? tables = null;
        var phases = data.Arrivals.Select(a => a.Phase).Distinct().ToList();
        try
        {
            if (method == RayMethod.FastMarching)
            {
                var fm = ForwardModel.FromInversionGrid(grid, sP.Select(x => 1 / x).ToArray(), sS.Select(x => 1 / x).ToArray(), domain);
                var used = data.Arrivals.Select(a => a.Station).ToHashSet();
                var sources = data.Stations.Where((_, i) => used.Contains(i))
                    .Select(st => new TableSource(st.Id, st.Lon, st.Lat, st.DepthKm)).ToList();
                tables = TravelTimeTableSet.Build(fm, sources, phases, Path.Combine(workFolder, "tables"), gpu,
                    progress == null ? null : new Progress<(double, string)>(p => progress.Report((progressBase + tableSpan * p.Item1, p.Item2))), log, ct);
            }
            var step = 0.5 * Math.Min(new SphericalGrid(domain.Forward).MinSpacingKm(), grid.MinSpacingKm());
            var bg = domain.Extended ? domain.Background : null;
            Func<double, double>? outsideP = bg == null ? null : d => 1 / bg.Vp(d);
            Func<double, double>? outsideS = bg == null ? null : d => 1 / bg.Vs(d);
            var done = 0;
            Parallel.ForEach(Enumerable.Range(0, data.Arrivals.Count), TomoStar.Core.Compute.ComputeSettings.Options(ct), ai =>
            {
                var a = data.Arrivals[ai];
                var ev = data.Events[a.Event];
                var stn = data.Stations[a.Station];
                var slow = a.Phase == Phase.P ? sP : sS;
                RayPath? ray;
                double tcalc, ge, gn, gz;
                if (tables != null)
                {
                    var table = tables.Get(stn.Id, a.Phase);
                    if (table == null) return;
                    ray = RayTracing.Backtrack(table, ev.Lon, ev.Lat, ev.DepthKm, step);
                    if (ray == null) return;
                    tcalc = table.Time(ev.Lon, ev.Lat, ev.DepthKm);
                    (ge, gn, gz) = table.Gradient(ev.Lon, ev.Lat, ev.DepthKm);
                }
                else
                {
                    ray = RayTracing.Straight(ev.Lon, ev.Lat, ev.DepthKm, stn.Lon, stn.Lat, stn.DepthKm, step);
                    var k0 = RayTracing.Kernel(grid, ray, null, a.Phase == Phase.P ? outsideP : outsideS, out var tOut);
                    tcalc = RayTracing.Time(k0, slow) + tOut;
                    var sSrc = bg != null && !grid.Contains(ev.Lon, ev.Lat, ev.DepthKm, 1e-6)
                        ? 1 / bg.Velocity(a.Phase, ev.DepthKm)
                        : grid.Interpolate(slow, ev.Lon, ev.Lat, ev.DepthKm);
                    (ge, gn, gz) = RayTracing.StraightGradient(ev.Lon, ev.Lat, ev.DepthKm, stn.Lon, stn.Lat, stn.DepthKm, sSrc);
                }
                if (!double.IsFinite(tcalc)) return;
                var kernel = RayTracing.Kernel(grid, ray, null, a.Phase == Phase.P ? outsideP : outsideS, out var outsideTime);
                var nodes = kernel.Keys.OrderBy(x => x).ToArray();
                var len = nodes.Select(x => (float)kernel[x]).ToArray();
                StoredRay? stored = null;
                if (keepRayPoints > 1)
                {
                    var d = ray.Decimated(keepRayPoints);
                    stored = new StoredRay(ai, a.Event, a.Station, a.Phase,
                        d.Lon.Select(x => (float)x).ToArray(), d.Lat.Select(x => (float)x).ToArray(), d.Depth.Select(x => (float)x).ToArray());
                }
                rows.Add(new RayRow(ai, nodes, len, tcalc, ge, gn, gz, stored) { OutsideTime = outsideTime });
                var c = Interlocked.Increment(ref done);
                if (c % 2000 == 0) progress?.Report((progressBase + tableSpan + (progressSpan - tableSpan) * c / data.Arrivals.Count, $"Rays {c}/{data.Arrivals.Count}"));
            });
        }
        finally
        {
            tables?.Dispose();
        }
        return rows.OrderBy(r => r.Arrival).ToList();
    }
}
