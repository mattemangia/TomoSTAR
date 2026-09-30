using TomoStar.Core.Geo;

namespace TomoStar.Core.Forward;

/// <summary>A node whose travel time is imposed rather than computed.</summary>
public readonly record struct SeedNode(int Index, double Time);

/// <summary>
/// The per-node metric of a <see cref="SphericalGrid"/>, precomputed once: the physical length of
/// one grid step along longitude, latitude and depth at every node, km.
/// </summary>
public sealed class GridMetric
{
    public GridMetric(SphericalGrid grid)
    {
        Grid = grid;
        HLon = new double[grid.Ny * grid.Nz];
        HLat = new double[grid.Nz];
        for (var k = 0; k < grid.Nz; k++)
        {
            HLat[k] = grid.RadiusKm[k] * grid.DLatRad;
            for (var j = 0; j < grid.Ny; j++)
                HLon[k * grid.Ny + j] = Math.Max(1e-9, grid.RadiusKm[k] * grid.CosLat[j] * grid.DLonRad);
        }
        HR = grid.DDepthKm;
    }

    public SphericalGrid Grid { get; }
    public double[] HLon { get; }
    public double[] HLat { get; }
    public double HR { get; }
}

/// <summary>
/// Fast marching on a spherical grid: solves |∇T| = s with the metric of the sphere,
///
///     (∂T/∂r)² + (1/r)²(∂T/∂φ)² + (1/(r cosφ))²(∂T/∂λ)² = s²,
///
/// discretised with the Godunov upwind scheme (Rouy &amp; Tourin 1992, SIAM J. Numer. Anal. 29(3),
/// 867-884) and marched in order of increasing time with a min-heap (Sethian 1996, PNAS 93(4),
/// 1591-1595), second-order where two upwind nodes are available and monotone (Sethian &amp;
/// Popovici 1999, Geophysics 64(2), 516-523). The spherical metric enters through the node
/// spacings; this is the forward solver of spherical teleseismic and regional tomography
/// (Rawlinson &amp; Sambridge 2004, Geophys. J. Int. 156(3), 631-647).
///
/// Around a point source the first-order scheme is at its worst (infinite wavefront curvature), so
/// the nodes within two cells of every seed get the straight-ray time with the trapezoidal mean
/// slowness, which is exact in a homogeneous medium; without this the near-source error propagates
/// to every downstream node.
/// </summary>
public static class SphericalEikonal
{
    private const byte Far = 0, Trial = 1, Alive = 2, Seed = 3;

    /// <summary>
    /// Seeds for a point source at an arbitrary position inside the grid: every node within
    /// <paramref name="radiusCells"/> times the LARGEST node spacing of the source (a sphere in km,
    /// not a cube of cells) gets the straight-ray time with the trapezoidal mean slowness.
    ///
    /// The radius is physical because the upwind stencil fails where the wavefront radius is
    /// comparable to the spacing along any axis: on a grid with 1 km vertical and 2.8 km lateral
    /// spacing a two-cell cube reaches only 2 km down, where the front is still far too curved, and
    /// the resulting 7 % error below the source then propagates everywhere.
    /// </summary>
    public static List<SeedNode> PointSource(SphericalGrid grid, double[] slowness, double lon, double lat, double depthKm,
        double originTime = 0, double radiusCells = 3)
    {
        var dep = Math.Clamp(depthKm, grid.DepthKm[0], grid.DepthKm[^1]);
        var (fx, fy, fz) = grid.Fractional(lon, lat, dep);
        fx = Math.Clamp(fx, 0, grid.Nx - 1);
        fy = Math.Clamp(fy, 0, grid.Ny - 1);
        fz = Math.Clamp(fz, 0, grid.Nz - 1);
        var src = GeoMath.ToCartesian(lon, lat, dep);
        var sSrc = grid.Interpolate(slowness, lon, lat, dep);
        var r = GeoMath.EarthRadiusKm - dep;
        var hLon = Math.Max(1e-9, r * Math.Cos(lat * GeoMath.Deg2Rad) * grid.DLonRad);
        var hLat = r * grid.DLatRad;
        var hR = grid.DDepthKm;
        var radiusKm = radiusCells * Math.Max(hLon, Math.Max(hLat, hR));
        int ri = (int)Math.Ceiling(radiusKm / hLon), rj = (int)Math.Ceiling(radiusKm / hLat), rk = (int)Math.Ceiling(radiusKm / hR);
        var i0 = (int)Math.Round(fx);
        var j0 = (int)Math.Round(fy);
        var k0 = (int)Math.Round(fz);
        var seeds = new List<SeedNode>();
        for (var k = Math.Max(0, k0 - rk); k <= Math.Min(grid.Nz - 1, k0 + rk); k++)
        for (var j = Math.Max(0, j0 - rj); j <= Math.Min(grid.Ny - 1, j0 + rj); j++)
        for (var i = Math.Max(0, i0 - ri); i <= Math.Min(grid.Nx - 1, i0 + ri); i++)
        {
            var idx = grid.Index(i, j, k);
            var d = (grid.Cartesian(i, j, k) - src).Length;
            // Always keep the enclosing cell so the source is never isolated.
            var inside = Math.Abs(i - fx) <= 1 && Math.Abs(j - fy) <= 1 && Math.Abs(k - fz) <= 1;
            if (d > radiusKm && !inside) continue;
            seeds.Add(new SeedNode(idx, originTime + d * 0.5 * (sSrc + slowness[idx])));
        }
        return seeds;
    }

    public static double[] Solve(GridMetric metric, double[] slowness, IReadOnlyList<SeedNode> seeds, CancellationToken ct = default)
    {
        var g = metric.Grid;
        if (slowness.Length != g.Count) throw new ArgumentException("Slowness and grid differ in size.");
        if (seeds.Count == 0) throw new ArgumentException("No seed nodes.");
        var time = GC.AllocateUninitializedArray<double>(g.Count);
        Array.Fill(time, double.PositiveInfinity);
        // Seeded nodes are imposed, not computed: the upwind stencil right next to a point source
        // underestimates the time (the front is curved on the scale of one cell), so letting the
        // march "improve" an analytic seed replaces an exact value with a worse one and the error
        // propagates to every node downstream. The relaxation solvers hold seeds fixed too; a seed
        // waiting in the heap is Seed rather than Trial, so no neighbour's update touches it.
        var status = new byte[g.Count];
        var heap = new FrontHeap(g.Count, Math.Max(1024, g.Count / 64));
        foreach (var s in seeds)
        {
            if ((uint)s.Index >= (uint)g.Count || !(s.Time < time[s.Index])) continue;
            time[s.Index] = s.Time;
            status[s.Index] = Seed;
            heap.Set(s.Index, s.Time);
        }
        var march = new March(metric, slowness, time, status);
        int nx = g.Nx, ny = g.Ny, nz = g.Nz, nxy = nx * ny;
        var counter = 0;
        while (heap.TryPop(out var idx))
        {
            status[idx] = Alive;
            if ((++counter & 0xFFFF) == 0) ct.ThrowIfCancellationRequested();
            var k = idx / nxy;
            var rest = idx - k * nxy;
            var j = rest / nx;
            var i = rest - j * nx;
            if (i > 0) Update(idx - 1, i - 1, j, k);
            if (i + 1 < nx) Update(idx + 1, i + 1, j, k);
            if (j > 0) Update(idx - nx, i, j - 1, k);
            if (j + 1 < ny) Update(idx + nx, i, j + 1, k);
            if (k > 0) Update(idx - nxy, i, j, k - 1);
            if (k + 1 < nz) Update(idx + nxy, i, j, k + 1);
        }
        return time;

        void Update(int n, int i, int j, int k)
        {
            if (status[n] >= Alive) return; // Alive, or a seed
            var t = march.SolveNode(n, i, j, k);
            if (!(t < time[n])) return;
            time[n] = t;
            status[n] = Trial;
            heap.Set(n, t);
        }
    }

    /// <summary>
    /// The upwind update of one node from its Alive neighbours. The arithmetic is the reference
    /// one term for term (the same expressions in the same order), so a field is the same to the
    /// last bit; only the bookkeeping is lighter: the node's indices come from the node that was
    /// just frozen instead of being divided out again, and the (at most three) upwind terms sit in
    /// locals instead of stack spans.
    /// </summary>
    private readonly struct March(GridMetric m, double[] slowness, double[] time, byte[] status)
    {
        private readonly int _nx = m.Grid.Nx, _ny = m.Grid.Ny, _nz = m.Grid.Nz;
        private readonly double[] _hLon = m.HLon, _hLat = m.HLat;
        private readonly double _hR = m.HR;

        public double SolveNode(int index, int i, int j, int k)
        {
            var s = slowness[index];
            if (!(s > 0) || !double.IsFinite(s)) return double.PositiveInfinity;
            var nxy = _nx * _ny;
            Terms q = default;
            Axis(index, i, _nx, 1, _hLon[k * _ny + j], ref q);
            Axis(index, j, _ny, _nx, _hLat[k], ref q);
            Axis(index, k, _nz, nxy, _hR, ref q);
            if (q.Used == 0) return double.PositiveInfinity;
            q.Sort();

            for (var terms = q.Used; terms >= 1; terms--)
            {
                double sw = 0, swt = 0, swt2 = 0;
                for (var r = 0; r < terms; r++)
                {
                    var he = q.H(r);
                    var te = q.T(r);
                    var w = 1.0 / (he * he);
                    sw += w;
                    swt += w * te;
                    swt2 += w * te * te;
                }
                var disc = swt * swt - sw * (swt2 - s * s);
                if (disc < 0) continue;
                var root = (swt + Math.Sqrt(disc)) / sw;
                var ok = true;
                for (var r = 0; r < terms; r++)
                    if (root < q.N(r)) { ok = false; break; }
                if (ok) return root;
            }
            return double.PositiveInfinity;
        }

        private void Axis(int index, int pos, int n, int stride, double h, ref Terms q)
        {
            var minus = pos > 0 && status[index - stride] == Alive ? time[index - stride] : double.PositiveInfinity;
            var plus = pos + 1 < n && status[index + stride] == Alive ? time[index + stride] : double.PositiveInfinity;
            double first;
            int second;
            if (minus <= plus)
            {
                first = minus;
                second = pos > 1 ? index - 2 * stride : -1;
            }
            else
            {
                first = plus;
                second = pos + 2 < n ? index + 2 * stride : -1;
            }
            if (!double.IsFinite(first)) return;
            var t2 = second >= 0 && status[second] == Alive ? time[second] : double.PositiveInfinity;
            if (double.IsFinite(t2) && t2 <= first) q.Add((4 * first - t2) / 3, 2 * h / 3, first);
            else q.Add(first, h, first);
        }
    }

    /// <summary>Up to three upwind terms: effective time, effective spacing, neighbour time.</summary>
    private struct Terms
    {
        private double _t0, _t1, _t2, _h0, _h1, _h2, _n0, _n1, _n2;
        public int Used;

        public readonly double T(int r) => r == 0 ? _t0 : r == 1 ? _t1 : _t2;
        public readonly double H(int r) => r == 0 ? _h0 : r == 1 ? _h1 : _h2;
        public readonly double N(int r) => r == 0 ? _n0 : r == 1 ? _n1 : _n2;

        public void Add(double t, double h, double n)
        {
            switch (Used++)
            {
                case 0: _t0 = t; _h0 = h; _n0 = n; break;
                case 1: _t1 = t; _h1 = h; _n1 = n; break;
                default: _t2 = t; _h2 = h; _n2 = n; break;
            }
        }

        /// <summary>By effective time, stable: the insertion sort of <see cref="SphericalEikonal.Sort"/>.</summary>
        public void Sort()
        {
            if (Used >= 2 && _t0 > _t1) Swap01();
            if (Used == 3 && _t1 > _t2)
            {
                (_t1, _t2) = (_t2, _t1); (_h1, _h2) = (_h2, _h1); (_n1, _n2) = (_n2, _n1);
                if (_t0 > _t1) Swap01();
            }
        }

        private void Swap01()
        {
            (_t0, _t1) = (_t1, _t0); (_h0, _h1) = (_h1, _h0); (_n0, _n1) = (_n1, _n0);
        }
    }

    private static void Sort(Span<double> t, Span<double> h, Span<double> n, int count)
    {
        for (var a = 1; a < count; a++)
        {
            double tt = t[a], hh = h[a], nn = n[a];
            var b = a - 1;
            while (b >= 0 && t[b] > tt)
            {
                t[b + 1] = t[b]; h[b + 1] = h[b]; n[b + 1] = n[b];
                b--;
            }
            t[b + 1] = tt; h[b + 1] = hh; n[b + 1] = nn;
        }
    }

    /// <summary>
    /// The same upwind update applied Gauss-Seidel fashion until nothing moves: the eight sweep
    /// orderings of the fast sweeping method, each walked plane by plane along i'+j'+k' = const
    /// (first order to its fixed point, then with second-order promotion). Slower than marching on
    /// a CPU; it exists as the host reference the OpenCL kernel is checked against, and it visits
    /// the planes in exactly the device's order, so the two agree to the last bit.
    /// </summary>
    public static double[] Relax(GridMetric m, double[] slowness, IReadOnlyList<SeedNode> seeds, int maxRounds = 0)
    {
        var g = m.Grid;
        var a = new double[g.Count];
        Array.Fill(a, double.PositiveInfinity);
        var fixedNode = new bool[g.Count];
        foreach (var s in seeds)
        {
            a[s.Index] = Math.Min(a[s.Index], s.Time);
            fixedNode[s.Index] = true;
        }
        int nx = g.Nx, ny = g.Ny, nz = g.Nz;
        var cap = maxRounds > 0 ? maxRounds : SweepRoundCap(g);
        for (var phase = 0; phase < 2; phase++)
        for (var round = 0; round < cap; round++)
        {
            var changed = false;
            for (var dir = 0; dir < 8; dir++)
            for (var p = 0; p < nx + ny + nz - 2; p++)
            {
                PlaneBounds(nx, ny, nz, p, out var i0, out var i1, out var k0, out var k1);
                for (var kp = k0; kp <= k1; kp++)
                for (var ip = i0; ip <= i1; ip++)
                {
                    var jp = p - ip - kp;
                    if (jp < 0 || jp >= ny) continue;
                    var i = (dir & 1) != 0 ? nx - 1 - ip : ip;
                    var j = (dir & 2) != 0 ? ny - 1 - jp : jp;
                    var k = (dir & 4) != 0 ? nz - 1 - kp : kp;
                    var idx = g.Index(i, j, k);
                    if (fixedNode[idx]) continue;
                    var cur = a[idx];
                    var cand = RelaxNode(m, slowness, a, idx, i, j, k, phase == 1);
                    if (!(cand < cur)) continue;
                    a[idx] = cand;
                    if (!double.IsFinite(cur) || cur - cand > SettleTolerance) changed = true;
                }
            }
            if (!changed) break;
        }
        return a;
    }

    /// <summary>
    /// A round that lowers no node by more than this (s) ends a relaxation phase: well under the
    /// float32 resolution of a stored table (about 4 µs at 60 s) and of any pick, and orders of
    /// magnitude above the creep of the last digits, which would otherwise cost several more rounds.
    /// </summary>
    public const double SettleTolerance = 1e-7;

    /// <summary>Most rounds of eight sweeps a relaxation phase may take before it is declared unsettled.</summary>
    internal static int SweepRoundCap(SphericalGrid g) => Math.Max(32, (g.Nx + g.Ny + g.Nz) / 4);

    /// <summary>
    /// The i' and k' ranges of the diagonal plane i'+j'+k' = <paramref name="p"/> (j' follows and
    /// may fall off the grid inside the box).
    /// </summary>
    internal static void PlaneBounds(int nx, int ny, int nz, int p, out int i0, out int i1, out int k0, out int k1)
    {
        i0 = Math.Max(0, p - (ny - 1) - (nz - 1));
        i1 = Math.Min(nx - 1, p);
        k0 = Math.Max(0, p - (nx - 1) - (ny - 1));
        k1 = Math.Min(nz - 1, p);
    }

    internal static double RelaxNode(GridMetric m, double[] slowness, double[] field, int idx, int i, int j, int k, bool secondOrder)
    {
        var g = m.Grid;
        var s = slowness[idx];
        if (!(s > 0)) return double.PositiveInfinity;
        Span<double> te = stackalloc double[3];
        Span<double> he = stackalloc double[3];
        Span<double> tn = stackalloc double[3];
        var used = 0;
        var nxy = g.Nx * g.Ny;
        AxisR(i, g.Nx, 1, m.HLon[k * g.Ny + j], ref used, te, he, tn);
        AxisR(j, g.Ny, g.Nx, m.HLat[k], ref used, te, he, tn);
        AxisR(k, g.Nz, nxy, m.HR, ref used, te, he, tn);
        if (used == 0) return double.PositiveInfinity;
        Sort(te, he, tn, used);
        for (var terms = used; terms >= 1; terms--)
        {
            double sw = 0, swt = 0, swt2 = 0;
            for (var q = 0; q < terms; q++)
            {
                var w = 1.0 / (he[q] * he[q]);
                sw += w; swt += w * te[q]; swt2 += w * te[q] * te[q];
            }
            var disc = swt * swt - sw * (swt2 - s * s);
            if (disc < 0) continue;
            var root = (swt + Math.Sqrt(disc)) / sw;
            var ok = true;
            for (var q = 0; q < terms; q++) if (root < tn[q]) { ok = false; break; }
            if (ok) return root;
        }
        return double.PositiveInfinity;

        void AxisR(int pos, int n, int stride, double h, ref int u, Span<double> tEff, Span<double> hEff, Span<double> tNb)
        {
            var minus = pos > 0 ? field[idx - stride] : double.PositiveInfinity;
            var plus = pos + 1 < n ? field[idx + stride] : double.PositiveInfinity;
            double first;
            int second;
            if (minus <= plus) { first = minus; second = pos > 1 ? idx - 2 * stride : -1; }
            else { first = plus; second = pos + 2 < n ? idx + 2 * stride : -1; }
            if (!double.IsFinite(first)) return;
            var t2 = secondOrder && second >= 0 ? field[second] : double.PositiveInfinity;
            if (double.IsFinite(t2) && t2 <= first) { tEff[u] = (4 * first - t2) / 3; hEff[u] = 2 * h / 3; }
            else { tEff[u] = first; hEff[u] = h; }
            tNb[u] = first;
            u++;
        }
    }

    /// <summary>
    /// The Trial front: a 4-ary min-heap of (time, node) with one entry per node, lowered in place
    /// when a node's time drops. A heap that re-pushes instead holds about three entries per node,
    /// pops and discards the stale ones, and is twice as deep; the pops are most of a march's time.
    /// Four children share a cache line, and each level halves the depth of a binary heap.
    /// </summary>
    private sealed class FrontHeap
    {
        private struct Entry
        {
            public double Key;
            public int Node;
        }

        private Entry[] _e;
        private readonly int[] _at; // heap slot of each node + 1; 0 = not in the heap
        private int _count;

        public FrontHeap(int nodes, int capacity)
        {
            _e = new Entry[capacity];
            _at = new int[nodes];
        }

        /// <summary>Inserts a node, or lowers its time when it is already in the heap.</summary>
        public void Set(int node, double key)
        {
            var i = _at[node] - 1;
            if (i < 0)
            {
                if (_count == _e.Length) Array.Resize(ref _e, _count * 2);
                i = _count++;
            }
            var e = _e;
            var at = _at;
            while (i > 0)
            {
                var p = (i - 1) >> 2;
                if (e[p].Key <= key) break;
                e[i] = e[p];
                at[e[i].Node] = i + 1;
                i = p;
            }
            e[i].Key = key;
            e[i].Node = node;
            at[node] = i + 1;
        }

        public bool TryPop(out int node)
        {
            if (_count == 0) { node = -1; return false; }
            var e = _e;
            var at = _at;
            node = e[0].Node;
            at[node] = 0;
            var last = e[--_count];
            var count = _count;
            if (count == 0) return true;
            var i = 0;
            while (true)
            {
                var c = 4 * i + 1;
                if (c >= count) break;
                var best = c;
                var bestKey = e[c].Key;
                var end = Math.Min(c + 4, count);
                for (var q = c + 1; q < end; q++)
                    if (e[q].Key < bestKey) { best = q; bestKey = e[q].Key; }
                if (bestKey >= last.Key) break;
                e[i] = e[best];
                at[e[i].Node] = i + 1;
                i = best;
            }
            e[i] = last;
            at[last.Node] = i + 1;
            return true;
        }
    }
}
