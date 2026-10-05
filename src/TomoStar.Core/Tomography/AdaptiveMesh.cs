// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using System.Text.Json;
using TomoStar.Core.Geo;
using TomoStar.Core.Numerics;

namespace TomoStar.Core.Tomography;

/// <summary>What the coverage field of an adaptive grid counts.</summary>
public enum AdaptiveCoverage
{
    /// <summary>Rays (or paths, or shots) through each node: a cell is refined where enough of them cross it.</summary>
    Hits,

    /// <summary>Derivative weight sum, km of ray per node (Toomey &amp; Foulger 1989).</summary>
    Dws
}

/// <summary>
/// Settings of an adaptive parameterisation. The inversion grid stays what it is (the forward and the
/// output live on it); the unknowns become the leaves of an octree over its nodes, fine where the data
/// sample the model densely and coarse where they do not.
/// </summary>
public sealed class AdaptiveGridSettings
{
    /// <summary>Invert on adaptive cells instead of on every node.</summary>
    public bool Enabled { get; set; }

    /// <summary>Largest cell, in nodes along longitude and latitude (a power of two).</summary>
    public int MaxCellNodes { get; set; } = 8;

    /// <summary>Largest cell, in nodes along depth (a power of two).</summary>
    public int MaxCellNodesDepth { get; set; } = 4;

    /// <summary>What the coverage counts (rays per node or ray length per node).</summary>
    public AdaptiveCoverage Coverage { get; set; } = AdaptiveCoverage.Hits;

    /// <summary>
    /// Coverage a child cell needs to exist: a cell is divided when at least
    /// <see cref="MinCoveredFraction"/> of its children would carry this much. In the units of
    /// <see cref="Coverage"/> (rays through the child's nodes, or km of ray).
    /// </summary>
    public double Threshold { get; set; } = 10;

    /// <summary>Fraction of the children that must reach the threshold for a cell to be divided.</summary>
    public double MinCoveredFraction { get; set; } = 0.5;

    /// <summary>
    /// Keep face-adjacent cells within a factor 2 of each other in size, so the smoothing between
    /// them compares like with like and no fine cell borders a much coarser one.
    /// </summary>
    public bool Balance { get; set; } = true;

    /// <summary>Rebuild the cells from the coverage of every iteration (rays move as the model changes).</summary>
    public bool RefineEachIteration { get; set; } = true;

    /// <summary>
    /// Allow a rebuild to make cells coarser than they were. Off (default), a zone once refined stays
    /// at least as fine, so a later step only adds detail.
    /// </summary>
    public bool AllowCoarsening { get; set; }

    public AdaptiveGridSettings Clone() => (AdaptiveGridSettings)MemberwiseClone();

    public string Describe() => Enabled
        ? $"adaptive cells up to {MaxCellNodes}×{MaxCellNodes}×{MaxCellNodesDepth} nodes, divided at {Threshold:0.###} {(Coverage == AdaptiveCoverage.Hits ? "rays" : "km of ray")} per child"
        : "every node";
}

/// <summary>
/// The index space of a regular grid of Nx × Ny × Nz nodes, x fastest (the layout of
/// <see cref="SphericalGrid"/>, of the FWI meshes with z as the slowest axis, and of 2-D maps with Nz = 1).
/// </summary>
public readonly record struct IndexSpace(int Nx, int Ny, int Nz)
{
    public int Count => Nx * Ny * Nz;
    public int Index(int i, int j, int k) => (k * Ny + j) * Nx + i;

    public (int I, int J, int K) Decompose(int index)
    {
        var k = index / (Nx * Ny);
        var rest = index - k * Nx * Ny;
        var j = rest / Nx;
        return (rest - j * Nx, j, k);
    }

    public static IndexSpace Of(SphericalGrid g) => new(g.Nx, g.Ny, g.Nz);
}

/// <summary>One leaf of an adaptive grid: the node block [I0,I1)×[J0,J1)×[K0,K1).</summary>
public readonly record struct AdaptiveCell(int I0, int I1, int J0, int J1, int K0, int K1)
{
    public int Ex => I1 - I0;
    public int Ey => J1 - J0;
    public int Ez => K1 - K0;
    public int NodeCount => Ex * Ey * Ez;
    public int MaxExtent => Math.Max(Ex, Math.Max(Ey, Ez));
    public double CentreI => 0.5 * (I0 + I1 - 1);
    public double CentreJ => 0.5 * (J0 + J1 - 1);
    public double CentreK => 0.5 * (K0 + K1 - 1);

    public int Extent(int axis) => axis switch { 0 => Ex, 1 => Ey, _ => Ez };
    public double Centre(int axis) => axis switch { 0 => CentreI, 1 => CentreJ, _ => CentreK };
}

/// <summary>
/// An adaptive parameterisation of a <see cref="SphericalGrid"/>: an octree whose leaves are blocks of
/// nodes. Every node belongs to exactly one leaf, and the leaf's single value is the unknown of all its
/// nodes, so a model on the leaves becomes a node field by the operator P (each node takes its leaf's
/// value) and a node sensitivity becomes a leaf sensitivity by Pᵀ (sums over the leaf). Nothing in the
/// forward changes: rays, eikonal tables and wave propagators keep running on the regular grid, and
/// only the unknowns of the inversion are fewer where the data cannot tell nodes apart.
///
/// This is the irregular parameterisation of seismic tomography (Abers &amp; Roecker 1991; Spakman &amp;
/// Bijwaard 2001, Pure Appl. Geophys. 158, 1401-1423; Sambridge &amp; Rawlinson 2005), built on the
/// octree so that cells nest and a later step can refine where an earlier one found resolution:
/// a cell is divided when its children would still be sampled (see
/// <see cref="AdaptiveGridSettings.Threshold"/>), adjacent cells are kept within a factor two of each
/// other, and a rebuild never coarsens a zone unless asked to.
///
/// Regularisation is written on the graph of face-adjacent cells so that it is the node Laplacian
/// when every cell is one node, and its continuum limit otherwise (second differences over the
/// distance between centres, rows weighted by √(cell volume)); damping of a cell of N nodes is
/// weighted by √N, the penalty the same update would pay on every node it moves.
/// </summary>
public sealed class AdaptiveMesh
{
    private readonly SphericalGrid? _grid;

    private AdaptiveMesh(IndexSpace space, SphericalGrid? grid, AdaptiveCell[] cells, int[] leafOfNode)
    {
        Space = space;
        _grid = grid;
        Cells = cells;
        LeafOfNode = leafOfNode;
    }

    /// <summary>The node index space the cells partition.</summary>
    public IndexSpace Space { get; }

    /// <summary>
    /// The spherical grid of the mesh; its geometry (spacings in km, axis weights, the saved definition)
    /// is what the smoothing, the cell sizes and <see cref="Save"/> need. A mesh built on a bare index
    /// space (an ambient-noise map, an FWI mesh) has none.
    /// </summary>
    public SphericalGrid Grid => _grid ?? throw new InvalidOperationException("This adaptive mesh partitions an index space, not a spherical grid.");

    public bool HasGrid => _grid != null;
    public AdaptiveCell[] Cells { get; }

    /// <summary>The leaf (column) of every grid node: the operator P as an index map.</summary>
    public int[] LeafOfNode { get; }

    public int CellCount => Cells.Length;

    /// <summary>Every node its own cell: the regular parameterisation as a mesh.</summary>
    public static AdaptiveMesh Regular(SphericalGrid grid)
    {
        var g = IndexSpace.Of(grid);
        var cells = new AdaptiveCell[g.Count];
        var map = new int[g.Count];
        for (var n = 0; n < g.Count; n++)
        {
            var (i, j, k) = g.Decompose(n);
            cells[n] = new AdaptiveCell(i, i + 1, j, j + 1, k, k + 1);
            map[n] = n;
        }
        return new AdaptiveMesh(g, grid, cells, map);
    }

    /// <summary>The mesh of the given leaves, which must partition the grid's nodes.</summary>
    public static AdaptiveMesh FromCells(SphericalGrid grid, IReadOnlyList<AdaptiveCell> cells) => FromCells(IndexSpace.Of(grid), grid, cells);

    private static AdaptiveMesh FromCells(IndexSpace g, SphericalGrid? grid, IReadOnlyList<AdaptiveCell> cells)
    {
        var map = new int[g.Count];
        Array.Fill(map, -1);
        for (var c = 0; c < cells.Count; c++)
        {
            var b = cells[c];
            if (b.I0 < 0 || b.J0 < 0 || b.K0 < 0 || b.I1 > g.Nx || b.J1 > g.Ny || b.K1 > g.Nz || b.NodeCount <= 0)
                throw new ArgumentException($"Cell {c} lies outside the grid.");
            for (var k = b.K0; k < b.K1; k++)
            for (var j = b.J0; j < b.J1; j++)
            for (var i = b.I0; i < b.I1; i++)
            {
                var n = g.Index(i, j, k);
                if (map[n] >= 0) throw new ArgumentException($"Cells {map[n]} and {c} overlap.");
                map[n] = c;
            }
        }
        if (Array.IndexOf(map, -1) >= 0) throw new ArgumentException("The cells do not cover every node.");
        return new AdaptiveMesh(g, grid, cells.ToArray(), map);
    }

    // ---- Construction ----

    /// <summary>
    /// Builds the cells from a coverage field on the nodes. Starts from the largest cells allowed and
    /// divides each one while its children would stay sampled; with a <paramref name="previous"/> mesh
    /// and coarsening off, a zone is divided at least as finely as it was.
    /// </summary>
    public static AdaptiveMesh Build(SphericalGrid grid, ReadOnlySpan<double> coverage, AdaptiveGridSettings settings, AdaptiveMesh? previous = null) =>
        Build(IndexSpace.Of(grid), grid, coverage, settings, previous);

    /// <summary>The same on a bare index space (an ambient-noise map, an FWI mesh): no geometry, no smoothing rows.</summary>
    public static AdaptiveMesh Build(IndexSpace space, ReadOnlySpan<double> coverage, AdaptiveGridSettings settings, AdaptiveMesh? previous = null) =>
        Build(space, null, coverage, settings, previous);

    private static AdaptiveMesh Build(IndexSpace g, SphericalGrid? grid, ReadOnlySpan<double> coverage, AdaptiveGridSettings settings, AdaptiveMesh? previous)
    {
        if (coverage.Length != g.Count) throw new ArgumentException("Coverage has the wrong length.");
        if (previous != null && previous.Space != g) throw new ArgumentException("The previous mesh is on another grid.");
        var sum = new PrefixSum3(g, coverage);
        // Node count of the previous leaf of each node, summed like the coverage: a box must be divided
        // when some previous leaf inside it was smaller than the box.
        int[]? prevSize = null;
        if (previous != null && !settings.AllowCoarsening)
        {
            prevSize = new int[g.Count];
            for (var n = 0; n < g.Count; n++) prevSize[n] = previous.Cells[previous.LeafOfNode[n]].NodeCount;
        }
        var minPrev = prevSize == null ? null : new PrefixMin3(g, prevSize);

        var hx = Math.Max(1, NextPow2(settings.MaxCellNodes));
        var hz = Math.Max(1, NextPow2(settings.MaxCellNodesDepth));
        var leaves = new List<AdaptiveCell>();
        var stack = new Stack<AdaptiveCell>();
        for (var k0 = 0; k0 < g.Nz; k0 += hz)
        for (var j0 = 0; j0 < g.Ny; j0 += hx)
        for (var i0 = 0; i0 < g.Nx; i0 += hx)
            stack.Push(new AdaptiveCell(i0, Math.Min(g.Nx, i0 + hx), j0, Math.Min(g.Ny, j0 + hx), k0, Math.Min(g.Nz, k0 + hz)));

        var children = new List<AdaptiveCell>(8);
        while (stack.Count > 0)
        {
            var box = stack.Pop();
            if (box.NodeCount == 1) { leaves.Add(box); continue; }
            Split(box, children);
            var divide = false;
            if (minPrev != null && minPrev.Min(box) < box.NodeCount) divide = true;
            if (!divide)
            {
                var covered = children.Count(c => sum.Sum(c) >= settings.Threshold);
                divide = covered >= Math.Max(1, (int)Math.Ceiling(settings.MinCoveredFraction * children.Count - 1e-12));
            }
            if (divide) foreach (var c in children) stack.Push(c);
            else leaves.Add(box);
        }
        var mesh = Order(g, grid, leaves);
        return settings.Balance ? mesh.Balanced() : mesh;
    }

    /// <summary>
    /// The children of a box: every axis longer than one node is halved (the lower half gets the
    /// floor), so a box of 8×8×4 nodes has eight children and one of 1×1×4 two.
    /// </summary>
    private static void Split(AdaptiveCell b, List<AdaptiveCell> children)
    {
        children.Clear();
        Span<int> xs = stackalloc int[3];
        Span<int> ys = stackalloc int[3];
        Span<int> zs = stackalloc int[3];
        var nx = Cuts(b.I0, b.I1, xs);
        var ny = Cuts(b.J0, b.J1, ys);
        var nz = Cuts(b.K0, b.K1, zs);
        for (var c = 0; c < nz; c++)
        for (var bj = 0; bj < ny; bj++)
        for (var a = 0; a < nx; a++)
            children.Add(new AdaptiveCell(xs[a], xs[a + 1], ys[bj], ys[bj + 1], zs[c], zs[c + 1]));
    }

    private static int Cuts(int lo, int hi, Span<int> cuts)
    {
        cuts[0] = lo;
        if (hi - lo <= 1) { cuts[1] = hi; return 1; }
        cuts[1] = lo + (hi - lo) / 2;
        cuts[2] = hi;
        return 2;
    }

    /// <summary>
    /// Divides cells until every face-adjacent pair is within a factor two in size (largest extent).
    /// Only divides, so it terminates, and a balanced mesh comes back unchanged.
    /// </summary>
    public AdaptiveMesh Balanced()
    {
        var mesh = this;
        var children = new List<AdaptiveCell>(8);
        for (var pass = 0; pass < 64; pass++)
        {
            var mark = new bool[mesh.CellCount];
            var any = false;
            foreach (var e in mesh.Edges)
            {
                int sa = mesh.Cells[e.A].MaxExtent, sb = mesh.Cells[e.B].MaxExtent;
                if (sa > 2 * sb) { mark[e.A] = true; any = true; }
                else if (sb > 2 * sa) { mark[e.B] = true; any = true; }
            }
            if (!any) return mesh;
            var leaves = new List<AdaptiveCell>(mesh.CellCount + 8 * mark.Count(m => m));
            for (var c = 0; c < mesh.CellCount; c++)
            {
                if (!mark[c] || mesh.Cells[c].NodeCount == 1) { leaves.Add(mesh.Cells[c]); continue; }
                Split(mesh.Cells[c], children);
                leaves.AddRange(children);
            }
            mesh = Order(Space, _grid, leaves);
        }
        return mesh;
    }

    /// <summary>Cells in grid order of their first node (depth slowest), with the node map.</summary>
    private static AdaptiveMesh Order(IndexSpace g, SphericalGrid? grid, List<AdaptiveCell> leaves)
    {
        var sorted = leaves.OrderBy(c => c.K0).ThenBy(c => c.J0).ThenBy(c => c.I0).ToArray();
        var map = new int[g.Count];
        Parallel.For(0, sorted.Length, c =>
        {
            var b = sorted[c];
            for (var k = b.K0; k < b.K1; k++)
            for (var j = b.J0; j < b.J1; j++)
            for (var i = b.I0; i < b.I1; i++)
                map[g.Index(i, j, k)] = c;
        });
        return new AdaptiveMesh(g, grid, sorted, map);
    }

    private static int NextPow2(int v)
    {
        var p = 1;
        while (p < v) p <<= 1;
        return p;
    }

    // ---- Adjacency ----

    /// <summary>
    /// A face contact between two cells across one axis (0 = longitude, 1 = latitude, 2 = depth):
    /// <see cref="Contact"/> node pairs of A and B are neighbours across it. A &lt; B.
    /// </summary>
    public readonly record struct Edge(int A, int B, int Axis, int Contact);

    private Edge[]? _edges;

    /// <summary>Every pair of face-adjacent cells, with their contact area in node pairs.</summary>
    public Edge[] Edges => _edges ??= BuildEdges();

    private Edge[] BuildEdges()
    {
        var g = Space;
        var counts = new Dictionary<(int, int, int), int>();
        for (var k = 0; k < g.Nz; k++)
        for (var j = 0; j < g.Ny; j++)
        for (var i = 0; i < g.Nx; i++)
        {
            var a = LeafOfNode[g.Index(i, j, k)];
            if (i + 1 < g.Nx) Count(a, LeafOfNode[g.Index(i + 1, j, k)], 0);
            if (j + 1 < g.Ny) Count(a, LeafOfNode[g.Index(i, j + 1, k)], 1);
            if (k + 1 < g.Nz) Count(a, LeafOfNode[g.Index(i, j, k + 1)], 2);
        }
        return counts.Select(kv => new Edge(kv.Key.Item1, kv.Key.Item2, kv.Key.Item3, kv.Value))
            .OrderBy(e => e.A).ThenBy(e => e.B).ThenBy(e => e.Axis).ToArray();

        void Count(int a, int b, int axis)
        {
            if (a == b) return;
            var key = a < b ? (a, b, axis) : (b, a, axis);
            counts[key] = counts.GetValueOrDefault(key) + 1;
        }
    }

    /// <summary>
    /// Coefficients of the cell Laplacian: for every cell a, the neighbours b with
    /// w_ab = W_axis · (contact / face area of a) / d_ab², d_ab the distance between the centres in
    /// nodes along the axis and W_axis the grid's axis weight (<see cref="Regularization.AxisWeights"/>).
    /// With one-node cells every w is W_axis: the node Laplacian.
    /// </summary>
    private List<(int B, double W)>[] LaplacianWeights(double verticalWeight, SmoothingScale scale)
    {
        var (wx, wy, wz) = Regularization.AxisWeights(Grid, verticalWeight, scale);
        var list = new List<(int, double)>[CellCount];
        for (var c = 0; c < CellCount; c++) list[c] = [];
        foreach (var e in Edges)
        {
            var wa = e.Axis switch { 0 => wx, 1 => wy, _ => wz };
            var ca = Cells[e.A];
            var cb = Cells[e.B];
            var d = Math.Abs(ca.Centre(e.Axis) - cb.Centre(e.Axis));
            if (d <= 0) continue;
            double FaceArea(AdaptiveCell c) => e.Axis switch { 0 => c.Ey * c.Ez, 1 => c.Ex * c.Ez, _ => c.Ex * c.Ey };
            list[e.A].Add((e.B, wa * e.Contact / FaceArea(ca) / (d * d)));
            list[e.B].Add((e.A, wa * e.Contact / FaceArea(cb) / (d * d)));
        }
        return list;
    }

    // ---- Regularisation rows ----

    /// <summary>
    /// Laplacian rows on the cells, the adaptive counterpart of <see cref="Regularization.AddLaplacian"/>:
    /// for each cell a, λ√Nₐ·Σ_b w_ab (x_a − x_b) = −(same on <paramref name="current"/>), so the roughness
    /// of the total deviation from the reference is penalised.
    /// </summary>
    public void AddLaplacian(CsrBuilder b, int offset, double lambda, double verticalWeight, SmoothingScale scale, double[]? current)
    {
        if (lambda <= 0) return;
        var weights = LaplacianWeights(verticalWeight, scale);
        var cols = new List<int>();
        var vals = new List<double>();
        for (var a = 0; a < CellCount; a++)
        {
            if (weights[a].Count == 0) continue;
            cols.Clear();
            vals.Clear();
            double diag = 0;
            foreach (var (nb, w) in weights[a])
            {
                cols.Add(offset + nb);
                vals.Add(-w);
                diag += w;
            }
            cols.Add(offset + a);
            vals.Add(diag);
            double rhs = 0;
            if (current != null) for (var q = 0; q < cols.Count; q++) rhs -= vals[q] * current[cols[q] - offset];
            var order = cols.Select((c, q) => (c, v: vals[q])).OrderBy(x => x.c).ToArray();
            b.AddRow(order.Select(x => x.c).ToArray(), order.Select(x => x.v).ToArray(), rhs, lambda * Math.Sqrt(Cells[a].NodeCount));
        }
    }

    /// <summary>λ√Nₐ·xₐ = 0 for every cell: damping of the update, with the weight of all the nodes it moves.</summary>
    public void AddDamping(CsrBuilder b, int offset, double lambda)
    {
        if (lambda <= 0) return;
        Span<int> c = stackalloc int[1];
        Span<double> v = stackalloc double[1];
        for (var a = 0; a < CellCount; a++)
        {
            c[0] = offset + a;
            v[0] = 1;
            b.AddRow(c, v, 0, lambda * Math.Sqrt(Cells[a].NodeCount));
        }
    }

    /// <summary>Roughness ‖L x‖ of a cell field, with the weights of <see cref="AddLaplacian"/>.</summary>
    public double Roughness(ReadOnlySpan<double> cellValues, double verticalWeight, SmoothingScale scale)
    {
        var weights = LaplacianWeights(verticalWeight, scale);
        double s = 0;
        for (var a = 0; a < CellCount; a++)
        {
            double lap = 0;
            foreach (var (nb, w) in weights[a]) lap += w * (cellValues[a] - cellValues[nb]);
            s += Cells[a].NodeCount * lap * lap;
        }
        return Math.Sqrt(s);
    }

    // ---- Operators ----

    /// <summary>P: the node field of a cell model (each node takes its cell's value).</summary>
    public double[] Prolong(ReadOnlySpan<double> cellValues)
    {
        if (cellValues.Length != CellCount) throw new ArgumentException("Cell field has the wrong length.");
        var f = new double[Space.Count];
        for (var n = 0; n < f.Length; n++) f[n] = cellValues[LeafOfNode[n]];
        return f;
    }

    /// <summary>Pᵀ: sums a node field over each cell (a node sensitivity becomes a cell sensitivity).</summary>
    public double[] SumToCells(ReadOnlySpan<double> nodeField)
    {
        if (nodeField.Length != Space.Count) throw new ArgumentException("Node field has the wrong length.");
        var s = new double[CellCount];
        for (var n = 0; n < nodeField.Length; n++) s[LeafOfNode[n]] += nodeField[n];
        return s;
    }

    /// <summary>Mean of a node field over each cell (the cell value of a node model).</summary>
    public double[] Restrict(ReadOnlySpan<double> nodeField)
    {
        var s = SumToCells(nodeField);
        for (var c = 0; c < CellCount; c++) s[c] /= Cells[c].NodeCount;
        return s;
    }

    /// <summary>
    /// Merges the node entries of a sensitivity row into cell columns starting at
    /// <paramref name="cellOffset"/> (one row of G·P); the result is sorted by column.
    /// </summary>
    public void MergeRow(ReadOnlySpan<int> nodes, ReadOnlySpan<double> values, int cellOffset, List<int> cols, List<double> vals)
    {
        cols.Clear();
        vals.Clear();
        var acc = new Dictionary<int, double>(nodes.Length);
        for (var q = 0; q < nodes.Length; q++)
        {
            var c = cellOffset + LeafOfNode[nodes[q]];
            acc[c] = acc.GetValueOrDefault(c) + values[q];
        }
        foreach (var (c, v) in acc.OrderBy(kv => kv.Key)) { cols.Add(c); vals.Add(v); }
    }

    // ---- Description ----

    /// <summary>Size of the cell of every node, km (geometric mean of its three physical extents).</summary>
    public double[] CellSizeKm()
    {
        var g = Grid;
        var cellSize = new double[CellCount];
        for (var c = 0; c < CellCount; c++)
        {
            var b = Cells[c];
            var (hx, hy, hz) = g.Spacing(Math.Clamp((int)Math.Round(b.CentreI), 0, g.Nx - 1),
                Math.Clamp((int)Math.Round(b.CentreJ), 0, g.Ny - 1), Math.Clamp((int)Math.Round(b.CentreK), 0, g.Nz - 1));
            cellSize[c] = Math.Cbrt(b.Ex * hx * b.Ey * hy * b.Ez * hz);
        }
        return Prolong(cellSize);
    }

    /// <summary>Cells by size, for a log line: "1 234 cells (1×1×1: 800, 2×2×2: 300, …)".</summary>
    public string Summary()
    {
        var groups = Cells.GroupBy(c => (c.Ex, c.Ey, c.Ez)).OrderBy(gr => gr.Key.Ex * gr.Key.Ey * gr.Key.Ez)
            .Select(gr => $"{gr.Key.Ex}×{gr.Key.Ey}×{gr.Key.Ez}: {gr.Count()}");
        return $"{CellCount:N0} cells for {Space.Count:N0} nodes ({string.Join(", ", groups)})";
    }

    private sealed record Stored(GridDefinition Grid, int[][] Cells);

    public void Save(string path) =>
        File.WriteAllText(path, JsonSerializer.Serialize(new Stored(Grid.Definition,
            Cells.Select(c => new[] { c.I0, c.I1, c.J0, c.J1, c.K0, c.K1 }).ToArray()), Model.TomoJson.Options));

    public static AdaptiveMesh Load(string path)
    {
        var s = JsonSerializer.Deserialize<Stored>(File.ReadAllText(path), Model.TomoJson.Options)
                ?? throw new InvalidDataException("Empty mesh file.");
        return FromCells(new SphericalGrid(s.Grid), s.Cells.Select(c => new AdaptiveCell(c[0], c[1], c[2], c[3], c[4], c[5])).ToList());
    }

    // ---- Box queries ----

    /// <summary>Summed-volume table: the sum of a node field over any box in O(1).</summary>
    private sealed class PrefixSum3
    {
        private readonly double[] _s;
        private readonly int _nx, _ny;

        public PrefixSum3(IndexSpace g, ReadOnlySpan<double> f)
        {
            _nx = g.Nx + 1;
            _ny = g.Ny + 1;
            _s = new double[(long)_nx * _ny * (g.Nz + 1) > int.MaxValue ? throw new ArgumentException("Grid too large.") : _nx * _ny * (g.Nz + 1)];
            for (var k = 0; k < g.Nz; k++)
            for (var j = 0; j < g.Ny; j++)
            for (var i = 0; i < g.Nx; i++)
            {
                var v = f[g.Index(i, j, k)];
                if (!double.IsFinite(v) || v < 0) v = 0;
                _s[At(i + 1, j + 1, k + 1)] = v
                    + _s[At(i, j + 1, k + 1)] + _s[At(i + 1, j, k + 1)] + _s[At(i + 1, j + 1, k)]
                    - _s[At(i, j, k + 1)] - _s[At(i, j + 1, k)] - _s[At(i + 1, j, k)]
                    + _s[At(i, j, k)];
            }
        }

        private int At(int i, int j, int k) => (k * _ny + j) * _nx + i;

        public double Sum(AdaptiveCell b) =>
            _s[At(b.I1, b.J1, b.K1)] - _s[At(b.I0, b.J1, b.K1)] - _s[At(b.I1, b.J0, b.K1)] - _s[At(b.I1, b.J1, b.K0)]
            + _s[At(b.I0, b.J0, b.K1)] + _s[At(b.I0, b.J1, b.K0)] + _s[At(b.I1, b.J0, b.K0)] - _s[At(b.I0, b.J0, b.K0)];
    }

    /// <summary>Minimum of an integer node field over a box (by scanning; boxes are small or rare).</summary>
    private sealed class PrefixMin3(IndexSpace g, int[] f)
    {
        public int Min(AdaptiveCell b)
        {
            var m = int.MaxValue;
            for (var k = b.K0; k < b.K1; k++)
            for (var j = b.J0; j < b.J1; j++)
            {
                var row = g.Index(0, j, k);
                for (var i = b.I0; i < b.I1; i++) m = Math.Min(m, f[row + i]);
            }
            return m;
        }
    }
}
