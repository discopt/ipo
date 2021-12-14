#pragma once

// #define IPO_DEBUG_LU_CHECK // Uncomment to enable dense checking of LU computations.
// #define IPO_DEBUG_LU_PRINT // Uncomment to print activity.
// #define IPO_DEBUG_LU_SAGE // Uncomdment to print matrices for sage.

#if defined(IPO_DEBUG_LU_SAGE)
#define IPO_DEBUG_LU_CHECK
#endif /* IPO_DEBUG_LU_SAGE */

#include <ipo/config.hpp>
#include <ipo/sparse_vector.hpp>

#include <vector>
#include <cmath>

namespace ipo
{
  IPO_NO_EXPORT
  std::size_t rowEchelon(std::size_t numColumns, std::vector<std::vector<double>>& matrix,
    std::size_t* rowPermutation = 0, std::size_t* columnPermutation = 0);

#if defined(IPO_RATIONAL)

  IPO_NO_EXPORT
  std::size_t rowEchelon(std::size_t numColumns, std::vector<std::vector<rational>>& matrix,
    std::size_t* rowPermutation = 0, std::size_t* columnPermutation = 0);

#endif /* IPO_RATIONAL */

  /**
   * \brief Class for an LU factorization.
   */

  template <typename T>
  class IncrementalLUFactorization
  {
  public:
    IncrementalLUFactorization()
    {

    }

    inline std::size_t size() const
    {
      return _leftRows.size();
    }

    /**
     * \brief Extends the factorized matrix by one row and one column.
     *
     * Extends the factorized matrix by one row and one column. Suppose we extend the factorization
     * \f$ L \cdot U \f$ to the factorization \f$ L' \cdot U' \f$ of
     * \f$ \begin{pmatrix} {L \cdot U} & b \\ a^\intercal & \beta \end{pmatrix} \f$.
     * We compute \f$ L' \coloneqq \begin{pmatrix} L & 0 \\ (U^{-\intercal} a)^\intercal & 1 \end{pmatrix} \f$ and
     * \f$ U' \coloneqq \begin{pmatrix} U & L^{-1} b \\ 0 & \beta - a^\intercal U^{-1} L^{-1} b \end{pmatrix} \f$.
     * One verifies \f$ L' \cdot U' = \begin{pmatrix} L \cdot U & L L^{-1} b \\ a^\intercal U^{-1} U & a^\intercal U^{-1} L^{-1} b
     * + \beta - a^\intercal U^{-1} L^{-1} b \end{pmatrix} = \begin{pmatrix} {L \cdot U} & b \\ a^\intercal & \beta \end{pmatrix} \f$.
     *
     * \param newRow Dense row to be added excluding diagonal entry. Will be modified.
     * \param newColumn Dense column to be added. Will be modified.
     * \param newDiagonal New diagonal entry to be added.
     * \param zeroEpsilon Tolerance for when an entry shall be considered zero.
     * \return Whether the matrix is regular. If it is singular, the factorization is not modified.
     */

    bool extend(T* newRow, T* newColumn, T newDiagonal, double epsilonEntry)
    {
#if defined(IPO_DEBUG_LU_PRINT)
      std::cout << "IncrementalLUFactorization::extend(\nrow:\n";
      for (std::size_t r = 0; r < size(); ++r)
        std::cout << " " << newRow[r];
      std::cout << " " << newDiagonal << "\ncolumn:\n";
      for (std::size_t r = 0; r < size(); ++r)
        std::cout << " " << newColumn[r];
      std::cout << " " << newDiagonal << "\nepsilon = " << epsilonEntry << ")" << std::endl;
#endif /* IPO_DEBUG_LU_PRINT */

#if defined(IPO_DEBUG_LU_SAGE)
      std::cout << "lu_A = Matrix(QQ, [";
      for (std::size_t r = 0; r < size(); ++r)
      {
        if (r > 0)
          std::cout << ",";
        std::cout << "[";
        for (std::size_t c = 0; c < size(); ++c)
        {
          if (c > 0)
            std::cout << ",";
          std::cout << _debugDenseMatrix[r][c];
        }
        std::cout << "]";
      }
      std::cout << "]); # SAGE\n";
      std::cout << "lu_b = matrix(QQ, [";
      for (std::size_t c = 0; c < size(); ++c)
      {
        if (c > 0)
          std::cout << ",";
        std::cout << "[" << newColumn[c] << "]";
      }
      std::cout << "]); # SAGE\n";
      std::cout << "lu_a = matrix(QQ, [";
      for (std::size_t r = 0; r < size(); ++r)
      {
        if (r > 0)
          std::cout << ",";
        std::cout << "[" << newRow[r] << "]";
      }
      std::cout << "]); # SAGE\n";
      std::cout << "lu_beta = matrix(QQ, [[" << newDiagonal << "]]); # SAGE\n";
      std::cout << "lu_A2 = block_matrix([[lu_A, lu_b], [lu_a.transpose(), lu_beta]]); # SAGE\n";
      std::cout << "det(lu_A) # SAGE\n";
      std::cout << "det(lu_A2) # SAGE\n";
#endif /* IPO_DEBUG_LU_SAGE */

#if defined(IPO_DEBUG_LU_CHECK)
      for (std::size_t r = 0; r < size(); ++r)
        _debugDenseMatrix[r].push_back(newColumn[r]);
      _debugDenseMatrix.push_back(std::vector<T>());
      for (std::size_t c = 0; c < size(); ++c)
        _debugDenseMatrix.back().push_back(newRow[c]);
      _debugDenseMatrix.back().push_back(newDiagonal);
#endif /* IPO_DEBUG_LU_CHECK */

      /* Compute L^{-1} b */

      solveLeft(newColumn);

#if defined(IPO_DEBUG_LU_PRINT)
      std::cout << "L^{-1}b:\n";
      for (std::size_t i = 0; i < size(); ++i)
        std::cout << " " << newColumn[i];
      std::cout << std::endl;
#endif /* IPO_DEBUG_LU_PRINT */

      /* Compute U^{-T} a. The transpose U^T has a row-wise sparse representation.
       * We swipe over U^T from left to right. In step i, we compute x_i as a_i / U_{i,i} and
       * subtract the remaining part of the i-th column of U^T times x_i from the corresponding part
       * of a.
       */

#if defined(IPO_DEBUG_LU_PRINT)
      std::cout << "U^{-T}a:\n";
#endif /* IPO_DEBUG_LU_PRINT */

      std::vector<std::size_t> positions(size(), 0);
      for (std::size_t i = 0; i < size(); ++i)
      {
        std::size_t p_i = positions[i];
        assert(p_i == _upperRows[i].size() - 1);
        assert(_upperRows[i][p_i] == i);
        const T& x_i = newRow[i] /= _upperEntries[i][p_i];
#if defined(IPO_DEBUG_LU_PRINT)
        std::cout << " " << x_i << std::flush;
#endif /* IPO_DEBUG_LU_PRINT */
        for (std::size_t j = i+1; j < size(); ++j)
        {
          std::size_t p_j = positions[j];
          if (p_j == _upperRows[j].size())
            continue;
          else if (_upperRows[j][p_j] == i)
          {
            newRow[j] -= _upperEntries[j][p_j] * x_i;
            positions[j]++;
          }
        }
      }
#if defined(IPO_DEBUG_LU_PRINT)
      std::cout << std::endl;
#endif /* IPO_DEBUG_LU_PRINT */

      /* Cache size since we change it. */
      const std::size_t n = size();

      /* Subtract the dot product of the new row and new column from the new diagonal. */
      for (std::size_t i = 0; i < n; ++i)
        newDiagonal -= newRow[i] * newColumn[i];
#if defined(IPO_DEBUG_LU_PRINT)
      std::cout << "new diagonal of U = " << newDiagonal << std::endl;
#endif /* IPO_DEBUG_LU_PRINT */
      if (fabs(convertNumber<double>(newDiagonal)) <= epsilonEntry)
      {
#if defined(IPO_DEBUG_LU_CHECK)
        for (std::size_t r = 0; r < size(); ++r)
          _debugDenseMatrix[r].pop_back();
        _debugDenseMatrix.pop_back();
#endif /* IPO_DEBUG_LU_CHECK */
        return false;
      }

      /* Add the new row to L. */
      for (std::size_t i = 0; i < n; ++i)
      {
        if (fabs(convertNumber<double>(newRow[i])) > epsilonEntry)
        {
          _leftRows[i].push_back(n);
          _leftEntries[i].push_back(newRow[i]);
        }
      }

      /* Add a new column with a 1 to L. */
      _leftRows.push_back(std::vector<std::size_t>(1, n));
      _leftEntries.push_back(std::vector<T>(1, T(1)));

      /* Add the new column to U. */
      _upperRows.push_back(std::vector<std::size_t>());
      std::vector<std::size_t>& lastColumnRows = _upperRows.back();
      _upperEntries.push_back(std::vector<T>());
      std::vector<T>& lastColumnEntries = _upperEntries.back();
      for (std::size_t i = 0; i < n; ++i)
      {
        if (fabs(convertNumber<double>(newColumn[i])) > epsilonEntry)
        {
          lastColumnRows.push_back(i);
          lastColumnEntries.push_back(newColumn[i]);
        }
      }
      lastColumnRows.push_back(n);
      lastColumnEntries.push_back(newDiagonal);

      assert(size() == n+1);

#if defined(IPO_DEBUG_LU_CHECK)

      // Dense versions of L and U.

      _debugDenseLeft.resize(size());
      _debugDenseUpper.resize(size());
      std::vector<std::vector<T>> _debugProduct(size());
      for (std::size_t r = 0; r < size(); ++r)
      {
        _debugDenseLeft[r] = std::vector<T>(size(), 0);
        _debugDenseUpper[r] = std::vector<T>(size(), 0);
        _debugProduct[r] = std::vector<T>(size(), 0);
      }
      for (std::size_t c = 0; c < size(); ++c)
      {
        for (std::size_t p = 0; p < _leftRows[c].size(); ++p)
          _debugDenseLeft[_leftRows[c][p]][c] = _leftEntries[c][p];
      }
      for (std::size_t c = 0; c < size(); ++c)
      {
        for (std::size_t p = 0; p < _upperRows[c].size(); ++p)
          _debugDenseUpper[_upperRows[c][p]][c] = _upperEntries[c][p];
      }

      for (std::size_t r = 0; r < size(); ++r)
      {
        for (std::size_t c = 0; c < size(); ++c)
        {
          for (std::size_t k = 0; k < size(); ++k)
            _debugProduct[r][c] += _debugDenseLeft[r][k] * _debugDenseUpper[k][c];
        }
      }

      // Check if factorized matrix is equal to the product.

      bool failure = false;
      for (std::size_t r = 0; r < size(); ++r)
      {
        for (std::size_t c = 0; c < size(); ++c)
        {
          if (fabs(convertNumber<double, T>(_debugDenseMatrix[r][c] - _debugProduct[r][c])) > 1.0e-9)
            failure = true;
        }
      }
      if (failure)
      {
        std::cout << "L*U is not equal to the factorized matrix." << std::endl;
        std::cout << "Matrix:\n";
        for (std::size_t r = 0; r < size(); ++r)
        {
          for (std::size_t c = 0; c < size(); ++c)
            std::cout << " " << _debugDenseMatrix[r][c];
          std::cout << std::endl;
        }
        std::cout << "L:\n";
        for (std::size_t r = 0; r < size(); ++r)
        {
          for (std::size_t c = 0; c < size(); ++c)
            std::cout << " " << _debugDenseLeft[r][c];
          std::cout << std::endl;
        }
        std::cout << "U:\n";
        for (std::size_t r = 0; r < size(); ++r)
        {
          for (std::size_t c = 0; c < size(); ++c)
            std::cout << " " << _debugDenseUpper[r][c];
          std::cout << std::endl;
        }
        std::cout << "L*U:\n";
        for (std::size_t r = 0; r < size(); ++r)
        {
          for (std::size_t c = 0; c < size(); ++c)
            std::cout << " " << _debugProduct[r][c];
          std::cout << std::endl;
        }
        std::cout << "A-L*U:\n";
        for (std::size_t r = 0; r < size(); ++r)
        {
          for (std::size_t c = 0; c < size(); ++c)
          {
            std::cout << " " << (_debugDenseMatrix[r][c] - _debugProduct[r][c]);
          }
          std::cout << std::endl;
        }
      }

#endif /* IPO_DEBUG_LU_CHECK */
      return true;
    }

    /**
     * \brief Solves Lx = r for x.
     *
     * Solves Lx = r for x.
     *
     * \param vector input vector r and output vector x.
     */

    void solveLeft(T* vector) const
    {
      for (std::size_t r = 0; r < size(); ++r)
      {
        assert(_leftRows[r].front() == r); // First entry should be diagonal.
        vector[r] /= _leftEntries[r][0];
        for (std::size_t p = 1; p < _leftRows[r].size(); ++p)
          vector[_leftRows[r][p]] -= _leftEntries[r][p] * vector[r];
      }
    }

    /**
     * \brief Solves Ux = r for x.
     *
     * Solves Ux = r for x.
     *
     * \param vector input vector r and output vector x.
     */

    void solveUpper(T* vector) const
    {
      for (std::size_t i = size(); i > 0; --i)
      {
        const std::size_t r = i-1;
        assert(_upperRows[r].back() == r); // Last entry should be diagonal.
        vector[r] /= _upperEntries[r].back();
        for (std::size_t p = 0; p < _upperRows[r].size() - 1; ++p)
          vector[_upperRows[r][p]] -= _upperEntries[r][p] * vector[r];
      }
    }

  protected:
    std::vector<std::vector<std::size_t>> _leftRows;
    std::vector<std::vector<T>> _leftEntries;
    std::vector<std::vector<std::size_t>> _upperRows;
    std::vector<std::vector<T>> _upperEntries;
#if defined(IPO_DEBUG_LU_CHECK)
    std::vector<std::vector<T>> _debugDenseLeft;
    std::vector<std::vector<T>> _debugDenseUpper;
    std::vector<std::vector<T>> _debugDenseMatrix;
#endif /* IPO_DEBUG_LU_CHECK */
  };

} /* namespace ipo */
