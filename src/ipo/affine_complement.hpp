#pragma once

#include <ipo/config.hpp>
#include <ipo/sparse_vector.hpp>

#include "lu.hpp"

namespace ipo
{
  template <typename T, typename IsZero>
  class AffineComplement
  {
  public:
    AffineComplement(std::size_t numVariables, const IsZero isZero = IsZero())
      : _numVariables(numVariables), _isZero(isZero), _lu(isZero),
      _columns(numVariables + 1)
    {
      for (ColumnData& columnData : _columns)
        columnData.basisIndex = std::numeric_limits<std::size_t>::max();
    }

    inline std::size_t numVariables() const
    {
      return _numVariables;
    }

    inline std::size_t rank() const
    {
      return _basisIndexToColumn.size();
    }

    void add(const sparse_vector<T>& row, const T& last, std::size_t newBasicColumn)
    {
      // Add column to basis.
      assert(newBasicColumn < _columns.size());
      assert(_columns[newBasicColumn].basisIndex == std::numeric_limits<std::size_t>::max());
      _columns[newBasicColumn].basisIndex = _basisIndexToColumn.size();
      _basisIndexToColumn.push_back(newBasicColumn);

      // Copy part of new row to basis matrix.
      std::vector<T> newRow(rank(), T(0));
      for (const auto& iter : row)
      {
        std::size_t basisIndex = _columns[iter.first].basisIndex;
        if (basisIndex < std::numeric_limits<std::size_t>::max())
          newRow[basisIndex] = iter.second;
      }
      if (last != 0 && _columns.back().basisIndex != std::numeric_limits<std::size_t>::max())
        newRow[_columns.back().basisIndex] = last;

      // Copy new basic column.
      std::vector<T> newColumn(rank() - 1, T(0));
      for (std::size_t i = 0; i < _columns[newBasicColumn].rows.size(); ++i)
        newColumn[_columns[newBasicColumn].rows[i]] = _columns[newBasicColumn].entries[i];

      // Update LU decomposition.
      _lu.extend(&newRow[0], &newColumn[0], newRow.back());

      // Add the row.
      for (const auto& iter : row)
      {
        ColumnData& colulmnData = _columns[iter.first];
        colulmnData.rows.push_back(_basisIndexToColumn.size() - 1);
        colulmnData.entries.push_back(T(0));
        colulmnData.entries.back() = iter.second;
      }
    }

    void computeKernelVector(std::size_t column, std::vector<std::size_t>& coordinates,
      std::vector<T>& entries)
    {
      std::vector<T> rhs(rank(), T(0));
      const ColumnData& columnData = _columns[column];
      for (std::size_t i = 0; i < columnData.rows.size(); ++i)
        rhs[columnData.rows[i]] = -columnData.entries[i];

      _lu.solveLeft(rhs);
      _lu.solveRight(rhs);

      coordinates.clear();
      entries.clear();
      for (std::size_t c = 0; c < numVariables(); ++c)
      {
        if (c == column)
        {
          coordinates.push_back(c);
          entries.push_back(T(1));
        }
        else if (_columns[c].basisIndex < std::numeric_limits<std::size_t>::max() && !_isZero(rhs[c]))
        {
          coordinates.push_back(c);
          entries.push_back(T(rhs[c]));
        }
      }
    }

  protected:
    struct ColumnData
    {
      std::size_t basisIndex;
      std::vector<std::size_t> rows;
      std::vector<T> entries;
    };

    std::size_t _numVariables;
    IsZero _isZero;
    IncrementalLUFactorization<T, IsZero> _lu;
    std::vector<ColumnData> _columns;
    std::vector<std::size_t> _basisIndexToColumn;
  };

} /* namespace ipo */


