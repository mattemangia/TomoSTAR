using TomoStar.Core.Geo;
using TomoStar.Core.Numerics;

namespace TomoStar.Core.Tomography;

/// <summary>
/// Where the velocity unknowns live when they are not the grid nodes themselves: a regular lattice
/// of parameter nodes, rotated about the vertical and shifted by fractions of its spacing; models on
/// several such lattices averaged remove the imprint of any one node geometry (e.g. Zaharia et al.
/// 2025, Appl. Sci. 15, 6084).
/// </summary>
public sealed class LatticeSettings
{
    public bool Enabled { get; set; }

    /// <summary>Horizontal spacing of the parameter nodes, km.</summary>
    public double SpacingKm { get; set; } = 5;

    /// <summary>Vertical spacing, km (0: the grid's own).</summary>
    public double DepthSpacingKm { get; set; }

    /// <summary>Rotation of the lattice axes from east-north, degrees anticlockwise.</summary>
    public double RotationDeg { get; set; }

    /// <summary>Shift of the lattice origin along its axes, in fractions of the spacing (0-1).</summary>
    public double OffsetU { get; set; }

    public double OffsetV { get; set; }
    public double OffsetZ { get; set; }

    public LatticeSettings Clone() => (LatticeSettings)MemberwiseClone();

    public string Describe() => Enabled
        ? string.Create(System.Globalization.CultureInfo.InvariantCulture,
            $"parameter lattice {SpacingKm:0.##} km{(DepthSpacingKm > 0 ? $" × {DepthSpacingKm:0.##} km" : "")}, rotated {RotationDeg:0.#}°, shifted ({OffsetU:0.##}, {OffsetV:0.##}, {OffsetZ:0.##})")
        : "parameters on the grid nodes";
}

/// <summary>
/// A rotated, shifted lattice of velocity parameters over a <see cref="SphericalGrid"/>: every grid
/// node takes its value by trilinear interpolation from the 8 lattice nodes around it (the operator
/// P of the system G·P), so the forward, the rays and the output stay on the grid. Only lattice
/// nodes that some grid node depends on are unknowns. The lattice lives in a local frame at the
/// centre of the grid (east, north, depth in km), rotated by the given angle.
/// </summary>
public sealed class ParameterLattice
{
    private readonly int[] _start;
    private readonly int[] _column;
    private readonly double[] _weight;
    private readonly int _nu, _nv, _nz;
    private readonly int[] _columnOfLattice; // lattice index → column, −1 inactive
    private readonly int[] _latticeOfColumn;
    private readonly double[] _mass;          // Σ weights of each column over the grid nodes

    public LatticeSettings Settings { get; }
    public double SpacingKm { get; }
    public double DepthSpacingKm { get; }

    /// <summary>Number of unknowns (active lattice nodes).</summary>
    public int Count => _latticeOfColumn.Length;

    public ParameterLattice(SphericalGrid grid, LatticeSettings s)
    {
        Settings = s;
        SpacingKm = Math.Max(1e-3, s.SpacingKm);
        DepthSpacingKm = s.DepthSpacingKm > 0 ? s.DepthSpacingKm : Math.Max(1e-3, grid.Definition.DDepth);
        double lon0 = grid.CentreLon, lat0 = 0.5 * (grid.LatDeg[0] + grid.LatDeg[^1]);
        var kx = GeoMath.Deg2Rad * GeoMath.EarthRadiusKm * Math.Cos(lat0 * GeoMath.Deg2Rad);
        var ky = GeoMath.Deg2Rad * GeoMath.EarthRadiusKm;
        double c = Math.Cos(s.RotationDeg * GeoMath.Deg2Rad), sn = Math.Sin(s.RotationDeg * GeoMath.Deg2Rad);
        var n = grid.Count;
        var u = new double[n]; var v = new double[n]; var z = new double[n];
        for (var node = 0; node < n; node++)
        {
            var (i, j, k) = grid.Decompose(node);
            double x = (GeoMath.UnwrapLon(grid.LonDeg[i], lon0) - lon0) * kx, y = (grid.LatDeg[j] - lat0) * ky;
            u[node] = x * c + y * sn;
            v[node] = -x * sn + y * c;
            z[node] = grid.DepthKm[k];
        }
        // Origin one spacing (times the shift) below the smallest coordinate, so every node is inside.
        double u0 = u.Min() - (s.OffsetU % 1) * SpacingKm, v0 = v.Min() - (s.OffsetV % 1) * SpacingKm, z0 = z.Min() - (s.OffsetZ % 1) * DepthSpacingKm;
        _nu = (int)Math.Floor((u.Max() - u0) / SpacingKm) + 2;
        _nv = (int)Math.Floor((v.Max() - v0) / SpacingKm) + 2;
        _nz = (int)Math.Floor((z.Max() - z0) / DepthSpacingKm) + 2;
        _columnOfLattice = Enumerable.Repeat(-1, _nu * _nv * _nz).ToArray();
        var entries = new List<(int Lattice, double W)>[n];
        for (var node = 0; node < n; node++)
        {
            double fu = (u[node] - u0) / SpacingKm, fv = (v[node] - v0) / SpacingKm, fz = (z[node] - z0) / DepthSpacingKm;
            int iu = Math.Min((int)Math.Floor(fu), _nu - 2), iv = Math.Min((int)Math.Floor(fv), _nv - 2), iz = Math.Min((int)Math.Floor(fz), _nz - 2);
            double tu = fu - iu, tv = fv - iv, tz = fz - iz;
            var list = new List<(int, double)>(8);
            for (var a = 0; a < 2; a++)
            for (var b = 0; b < 2; b++)
            for (var d = 0; d < 2; d++)
            {
                var w = (a == 0 ? 1 - tu : tu) * (b == 0 ? 1 - tv : tv) * (d == 0 ? 1 - tz : tz);
                if (w < 1e-9) continue;
                list.Add((Index(iu + a, iv + b, iz + d), w));
            }
            entries[node] = list;
        }
        // Active lattice nodes get columns in lattice order.
        foreach (var (l, _) in entries.SelectMany(e => e)) _columnOfLattice[l] = 0;
        var cols = new List<int>();
        for (var l = 0; l < _columnOfLattice.Length; l++)
            if (_columnOfLattice[l] == 0) { _columnOfLattice[l] = cols.Count; cols.Add(l); }
        _latticeOfColumn = cols.ToArray();
        _start = new int[n + 1];
        for (var node = 0; node < n; node++) _start[node + 1] = _start[node] + entries[node].Count;
        _column = new int[_start[n]];
        _weight = new double[_start[n]];
        _mass = new double[Count];
        for (var node = 0; node < n; node++)
            for (var q = 0; q < entries[node].Count; q++)
            {
                var col = _columnOfLattice[entries[node][q].Lattice];
                _column[_start[node] + q] = col;
                _weight[_start[node] + q] = entries[node][q].W;
                _mass[col] += entries[node][q].W;
            }
    }

    private int Index(int iu, int iv, int iz) => (iz * _nv + iv) * _nu + iu;

    /// <summary>The columns and interpolation weights of grid node <paramref name="node"/>.</summary>
    public (ReadOnlyMemory<int> Columns, ReadOnlyMemory<double> Weights) Of(int node) =>
        (_column.AsMemory(_start[node], _start[node + 1] - _start[node]), _weight.AsMemory(_start[node], _start[node + 1] - _start[node]));

    /// <summary>Adds value × weight to every column of <paramref name="node"/>.</summary>
    public void Spread(int node, double value, Action<int, double> add)
    {
        for (var q = _start[node]; q < _start[node + 1]; q++) add(_column[q], value * _weight[q]);
    }

    /// <summary>The grid-node field of a parameter vector (P·x).</summary>
    public double[] Prolong(ReadOnlySpan<double> parameters, int offset = 0)
    {
        var n = _start.Length - 1;
        var f = new double[n];
        for (var node = 0; node < n; node++)
        {
            double sum = 0;
            for (var q = _start[node]; q < _start[node + 1]; q++) sum += _weight[q] * parameters[offset + _column[q]];
            f[node] = sum;
        }
        return f;
    }

    /// <summary>Parameter values that best stand for a grid-node field: the interpolation-weighted mean of the nodes each parameter reaches.</summary>
    public double[] Restrict(ReadOnlySpan<double> nodeField)
    {
        var sum = new double[Count];
        for (var node = 0; node < nodeField.Length; node++)
            for (var q = _start[node]; q < _start[node + 1]; q++) sum[_column[q]] += _weight[q] * nodeField[node];
        for (var c = 0; c < Count; c++) sum[c] = _mass[c] > 0 ? sum[c] / _mass[c] : 0;
        return sum;
    }

    /// <summary>λ√m·x = 0 per parameter, m the number of grid nodes it stands for (Σ weights): the node damping carried to the lattice.</summary>
    public void AddDamping(CsrBuilder b, int offset, double lambda)
    {
        if (lambda <= 0) return;
        Span<int> c = stackalloc int[1];
        Span<double> v = stackalloc double[1];
        for (var p = 0; p < Count; p++)
        {
            c[0] = offset + p;
            v[0] = 1;
            b.AddRow(c, v, 0, lambda * Math.Sqrt(Math.Max(1e-6, _mass[p])));
        }
    }

    /// <summary>
    /// Second differences on the lattice (its own axes: rotated u, v and depth), of the total deviation
    /// from the reference: λ√m·Σ_axes w(2x − x₋ − x₊) = −(same on <paramref name="current"/>). The axis
    /// weights follow <see cref="Regularization.AxisWeights"/>: in km, a vertical difference costs
    /// verticalWeight·(h/dz)² of a horizontal one of the same size; between nodes, verticalWeight.
    /// </summary>
    public void AddLaplacian(CsrBuilder b, int offset, double lambda, double verticalWeight, SmoothingScale scale, double[]? current)
    {
        if (lambda <= 0) return;
        var wz = scale == SmoothingScale.Nodes ? verticalWeight : verticalWeight * Math.Pow(SpacingKm / DepthSpacingKm, 2);
        var cols = new List<int>(7);
        var vals = new List<double>(7);
        for (var p = 0; p < Count; p++)
        {
            var l = _latticeOfColumn[p];
            int iu = l % _nu, iv = l / _nu % _nv, iz = l / (_nu * _nv);
            cols.Clear();
            vals.Clear();
            double diag = 0;
            void Pair(int du, int dv, int dz, double w)
            {
                int a = iu - du, bb = iv - dv, cz = iz - dz, e = iu + du, f = iv + dv, g = iz + dz;
                if (a < 0 || bb < 0 || cz < 0 || e >= _nu || f >= _nv || g >= _nz) return;
                int c1 = _columnOfLattice[Index(a, bb, cz)], c2 = _columnOfLattice[Index(e, f, g)];
                if (c1 < 0 || c2 < 0) return;
                cols.Add(offset + c1); vals.Add(-w);
                cols.Add(offset + c2); vals.Add(-w);
                diag += 2 * w;
            }
            Pair(1, 0, 0, 1); Pair(0, 1, 0, 1); Pair(0, 0, 1, wz);
            if (diag == 0) continue;
            cols.Add(offset + p);
            vals.Add(diag);
            double rhs = 0;
            if (current != null) for (var q = 0; q < cols.Count; q++) rhs -= vals[q] * current[cols[q] - offset];
            var merged = cols.Zip(vals).GroupBy(x => x.First).Select(gp => (gp.Key, gp.Sum(x => x.Second))).OrderBy(x => x.Key).ToArray();
            b.AddRow(merged.Select(x => x.Key).ToArray(), merged.Select(x => x.Item2).ToArray(), rhs, lambda * Math.Sqrt(Math.Max(1e-6, _mass[p])));
        }
    }

    public string Summary() => string.Create(System.Globalization.CultureInfo.InvariantCulture,
        $"{Count:N0} lattice parameters ({_nu}×{_nv}×{_nz} lattice, {SpacingKm:0.##} km × {DepthSpacingKm:0.##} km, rotated {Settings.RotationDeg:0.#}°)");
}
