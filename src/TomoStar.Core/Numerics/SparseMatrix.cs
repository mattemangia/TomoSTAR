// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

namespace TomoStar.Core.Numerics;

/// <summary>What LSQR needs of a matrix: its shape and the products A·x and Aᵀ·x.</summary>
public interface ILinearOperator
{
    int RowCount { get; }
    int ColumnCount { get; }

    /// <summary>y = A x</summary>
    void Multiply(ReadOnlySpan<double> x, Span<double> y);

    /// <summary>y = Aᵀ x</summary>
    void TransposeMultiply(ReadOnlySpan<double> x, Span<double> y);
}

/// <summary>
/// Compressed-sparse-row matrix. The transpose is built once, lazily, so Aᵀx is also a row-parallel
/// product with no write contention and a fixed summation order (the result does not depend on
/// the thread count).
/// </summary>
public sealed class CsrMatrix : ILinearOperator
{
    private CsrMatrix? _transpose;

    public CsrMatrix(int rowCount, int columnCount, int[] rowStart, int[] columnIndex, double[] values)
    {
        if (rowStart.Length != rowCount + 1) throw new ArgumentException("rowStart must have RowCount+1 entries.");
        if (columnIndex.Length != values.Length) throw new ArgumentException("columnIndex and values differ in length.");
        RowCount = rowCount;
        ColumnCount = columnCount;
        RowStart = rowStart;
        ColumnIndex = columnIndex;
        Values = values;
    }

    public int RowCount { get; }
    public int ColumnCount { get; }
    public int[] RowStart { get; }
    public int[] ColumnIndex { get; }
    public double[] Values { get; }
    public int NonZeroCount => Values.Length;

    public void Multiply(ReadOnlySpan<double> x, Span<double> y)
    {
        if (x.Length != ColumnCount || y.Length != RowCount) throw new ArgumentException("Dimension mismatch.");
        var xs = x.ToArray();
        var ys = new double[RowCount];
        MultiplyInto(xs, ys);
        ys.CopyTo(y);
    }

    /// <summary>Array form, which avoids the span copies on the hot path of LSQR.</summary>
    public void MultiplyInto(double[] x, double[] y)
    {
        var rowStart = RowStart;
        var cols = ColumnIndex;
        var vals = Values;
        const int block = 256;
        var blocks = (RowCount + block - 1) / block;
        Parallel.For(0, blocks, b =>
        {
            var end = Math.Min(RowCount, (b + 1) * block);
            for (var i = b * block; i < end; i++)
            {
                double sum = 0;
                for (var k = rowStart[i]; k < rowStart[i + 1]; k++) sum += vals[k] * x[cols[k]];
                y[i] = sum;
            }
        });
    }

    public void TransposeMultiply(ReadOnlySpan<double> x, Span<double> y) => Transpose().Multiply(x, y);

    public void TransposeMultiplyInto(double[] x, double[] y) => Transpose().MultiplyInto(x, y);

    public CsrMatrix Transpose()
    {
        if (_transpose != null) return _transpose;
        var counts = new int[ColumnCount + 1];
        foreach (var c in ColumnIndex) counts[c + 1]++;
        for (var c = 0; c < ColumnCount; c++) counts[c + 1] += counts[c];
        var next = (int[])counts.Clone();
        var tCols = new int[ColumnIndex.Length];
        var tVals = new double[Values.Length];
        for (var i = 0; i < RowCount; i++)
        for (var k = RowStart[i]; k < RowStart[i + 1]; k++)
        {
            var p = next[ColumnIndex[k]]++;
            tCols[p] = i;
            tVals[p] = Values[k];
        }
        _transpose = new CsrMatrix(ColumnCount, RowCount, counts, tCols, tVals) { _transpose = this };
        return _transpose;
    }

    /// <summary>Sum of |a_ij| per column: the ray-coverage (hit-length) of each model parameter.</summary>
    public double[] ColumnAbsSums(int rowLimit = int.MaxValue)
    {
        var sums = new double[ColumnCount];
        var rows = Math.Min(rowLimit, RowCount);
        for (var i = 0; i < rows; i++)
        for (var k = RowStart[i]; k < RowStart[i + 1]; k++)
            sums[ColumnIndex[k]] += Math.Abs(Values[k]);
        return sums;
    }
}

/// <summary>
/// Builds a <see cref="CsrMatrix"/> and its right-hand side one row at a time. A row weight is
/// folded into both the row and its right-hand side, so the solver minimises ‖W(Ax − b)‖² without
/// knowing weights exist. Regularisation rows (damping, smoothing) are ordinary rows.
/// </summary>
public sealed class CsrBuilder
{
    private readonly List<int> _rowStart = [0];
    private readonly List<int> _columns = [];
    private readonly List<double> _values = [];
    private readonly List<double> _rhs = [];

    public CsrBuilder(int columnCount) => ColumnCount = columnCount;

    public int ColumnCount { get; }
    public int RowCount => _rowStart.Count - 1;

    public void AddRow(ReadOnlySpan<int> columns, ReadOnlySpan<double> values, double rhs, double weight = 1.0)
    {
        if (columns.Length != values.Length) throw new ArgumentException("columns and values differ in length.");
        for (var n = 0; n < columns.Length; n++)
        {
            if (values[n] == 0) continue;
            if ((uint)columns[n] >= (uint)ColumnCount) throw new ArgumentOutOfRangeException(nameof(columns));
            _columns.Add(columns[n]);
            _values.Add(values[n] * weight);
        }
        _rowStart.Add(_columns.Count);
        _rhs.Add(rhs * weight);
    }

    public void AddRow(IReadOnlyDictionary<int, double> entries, double rhs, double weight = 1.0)
    {
        foreach (var (c, v) in entries.OrderBy(e => e.Key))
        {
            if (v == 0) continue;
            if ((uint)c >= (uint)ColumnCount) throw new ArgumentOutOfRangeException(nameof(entries));
            _columns.Add(c);
            _values.Add(v * weight);
        }
        _rowStart.Add(_columns.Count);
        _rhs.Add(rhs * weight);
    }

    /// <summary>
    /// The matrix and right-hand side. The builder is emptied: each list is copied and released in
    /// turn, so a large system does not hold the growing lists and all the arrays at once.
    /// </summary>
    public (CsrMatrix Matrix, double[] Rhs) Build()
    {
        var rows = RowCount;
        static T[] Take<T>(List<T> list)
        {
            var a = list.ToArray();
            list.Clear();
            list.TrimExcess();
            return a;
        }
        var start = Take(_rowStart);
        var columns = Take(_columns);
        var values = Take(_values);
        var rhs = Take(_rhs);
        _rowStart.Add(0);
        return (new CsrMatrix(rows, ColumnCount, start, columns, values), rhs);
    }
}
