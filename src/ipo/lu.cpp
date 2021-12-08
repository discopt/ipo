#include "lu.hpp"

namespace ipo
{

  IPO_NO_EXPORT
  std::size_t rowEchelon(std::size_t numColumns, std::vector<std::vector<double>>& matrix,
    std::size_t* rowPermutation, std::size_t* columnPermutation)
  {
    bool deleteRowPermutation = rowPermutation == 0;
    bool deleteColumnPermutation = columnPermutation == 0;
    if (deleteRowPermutation)
      rowPermutation = new std::size_t[matrix.size()];
    if (deleteColumnPermutation)
      columnPermutation = new std::size_t[numColumns];

    std::size_t numRows = matrix.size();
    for (std::size_t r = 0; r < numRows; ++r)
      rowPermutation[r] = r;
    for (std::size_t c = 0; c < numColumns; ++c)
      columnPermutation[c] = c;

    std::size_t rank = 0;
    while (true)
    {
      double pivotAbsolute = 0.0;
      std::size_t pivotRow = 0;
      std::size_t pivotColumn = 0;
      for (std::size_t r = rank; r < numRows; ++r)
      {
        for (std::size_t c = rank; c < numColumns; ++c)
        {
          double x = fabs(matrix[rowPermutation[r]][columnPermutation[c]]);
          if (x > pivotAbsolute)
          {
            pivotAbsolute = x;
            pivotRow = r;
            pivotColumn = c;
          }
        }
      }

      if (pivotAbsolute <= 1.0e-12)
        break;

      // Permute pivot to rank,rank.

      std::swap(rowPermutation[rank], rowPermutation[pivotRow]);
      std::swap(columnPermutation[rank], columnPermutation[pivotColumn]);

      // Perform row operations.
      for (std::size_t r = rank+1; r < numRows; ++r)
      {
        double factor = matrix[rowPermutation[r]][columnPermutation[rank]]
          / matrix[rowPermutation[rank]][columnPermutation[rank]];
        if (fabs(factor) <= 1.0e-7)
          continue;

        std::size_t pr = rowPermutation[r];
        for (std::size_t c = rank; c < numColumns; ++c)
        {
          double x = matrix[pr][columnPermutation[c]]
            -factor * matrix[rowPermutation[rank]][columnPermutation[c]];
          matrix[pr][columnPermutation[c]] = x;
        }
      }

      ++rank;
    }
    
    if (deleteColumnPermutation)
      delete[] columnPermutation;
    if (deleteRowPermutation)
      delete[] rowPermutation;

    return rank;
  }

#if defined(IPO_WITH_GMP)

  IPO_NO_EXPORT
  std::size_t rowEchelon(std::size_t numColumns, std::vector<std::vector<mpq_class>>& matrix,
    std::size_t* rowPermutation, std::size_t* columnPermutation)
  {
    bool deleteRowPermutation = rowPermutation == 0;
    bool deleteColumnPermutation = columnPermutation == 0;
    if (deleteRowPermutation)
      rowPermutation = new std::size_t[matrix.size()];
    if (deleteColumnPermutation)
      columnPermutation = new std::size_t[numColumns];

    std::size_t numRows = matrix.size();
    for (std::size_t r = 0; r < numRows; ++r)
      rowPermutation[r] = r;
    for (std::size_t c = 0; c < numColumns; ++c)
      columnPermutation[c] = c;

    std::size_t rank = 0;
    while (true)
    {
      double pivotAbsolute = 0.0;
      std::size_t pivotRow = 0;
      std::size_t pivotColumn = 0;
      for (std::size_t r = rank; r < numRows; ++r)
      {
        for (std::size_t c = rank; c < numColumns; ++c)
        {
          double x = fabs(matrix[rowPermutation[r]][columnPermutation[c]].get_d());
          if (x > pivotAbsolute)
          {
            pivotAbsolute = x;
            pivotRow = r;
            pivotColumn = c;
          }
        }
      }

      if (pivotAbsolute <= 1.0e-6)
        break;

      // Permute pivot to rank,rank.

      std::swap(rowPermutation[rank], rowPermutation[pivotRow]);
      std::swap(columnPermutation[rank], columnPermutation[pivotColumn]);

      // Perform row operations.
      for (std::size_t r = rank+1; r < numRows; ++r)
      {
        mpq_class factor = matrix[rowPermutation[r]][columnPermutation[rank]]
          / matrix[rowPermutation[rank]][columnPermutation[rank]];
        if (factor != 0)
          continue;

        std::size_t pr = rowPermutation[r];
        for (std::size_t c = rank; c < numColumns; ++c)
        {
          mpq_class x = matrix[pr][columnPermutation[c]]
            - factor * matrix[rowPermutation[rank]][columnPermutation[c]];
          matrix[pr][columnPermutation[c]] = x;
        }
      }

      ++rank;
    }
    
    if (deleteColumnPermutation)
      delete[] columnPermutation;
    if (deleteRowPermutation)
      delete[] rowPermutation;

    return rank;
  }

#endif /* IPO_WITH_GMP */

}
