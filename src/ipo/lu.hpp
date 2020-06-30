#pragma once

#include <ipo/config.hpp>

#include <vector>
#include <cmath>

#if defined(IPO_WITH_GMP)
#include <gmpxx.h>
#endif /* IPO_WITH_GMP */

namespace ipo
{

  /**
   * \brief Class for an LU factorization.
   */
  
  template <typename T, typename IsZero>
  class IncrementalLUFactorization
  {
  public:
    IncrementalLUFactorization(const IsZero isZero)
      : _isZero(isZero)
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
     *
     * 
     * \param newRow Dense row to be added excluding diagonal entry. Will be modified.
     * \param newColumn Dense column to be added. Will be modified.
     * \param newDiagonal New diagonal entry to be added.
     */

    void extend(T* newRow, T* newColumn, T newDiagonal)
    {
#ifdef IPO_DEBUG_LU
      std::cout << "IncrementalLUFactorization::extend" << std::endl;
      std::cout << "  a = [";
      for (std::size_t i = 0; i < size(); ++i)
        std::cout << (i > 0 ? " " : "") << newRow[i];
      std::cout << "]" << std::endl;
      std::cout << "  b = [";
      for (std::size_t i = 0; i < size(); ++i)
        std::cout << (i > 0 ? " " : "") << newColumn[i];
      std::cout << "]" << std::endl;
      std::cout << "  beta = " << newDiagonal << std::endl;
#endif /* IPO_DEBUG_LU */        

      /* Compute L^{-1} b */

      solveLeft(newColumn);

#ifdef IPO_DEBUG_LU
      std::cout << "L^{-1}b = [";
      for (std::size_t i = 0; i < size(); ++i)
        std::cout << (i > 0 ? " " : "") << newColumn[i];
      std::cout << "]" << std::endl;
#endif /* IPO_DEBUG_LU */        

      /* Compute U^{-T} a. The transpose U^T has a row-wise sparse representation.
       * We swipe over U^T from left to right. In step i, we compute x_i as a_i / U_{i,i} and
       * subtract the remaining part of the i-th column of U^T times x_i from the corresponding part
       * of a.
       */

#ifdef IPO_DEBUG_LU
      std::cout << "U^{-T}a = [";
#endif /* IPO_DEBUG_LU */

      std::vector<std::size_t> positions(size(), 0);
      for (std::size_t i = 0; i < size(); ++i)
      {
        std::size_t p_i = positions[i];
        assert(p_i == _upperRows[i].size() - 1);
        assert(_upperRows[i][p_i] == i);
        const T& x_i = newRow[i] /= _upperEntries[i][p_i];
#ifdef IPO_DEBUG_LU
        std::cout << (i > 0 ? " " : "") << x_i << std::flush;
#endif /* IPO_DEBUG_LU */
        for (std::size_t j = i+1; j < size(); ++j)
        {
          std::size_t p_j = positions[j];
          if (p_j == _upperRows[j].size())
            continue;
          else if (_upperRows[j][p_j] == i)
            newRow[j] -= _upperEntries[j][p_j] * x_i;
          positions[j]++;
        }
      }
#ifdef IPO_DEBUG_LU
      std::cout << "]" << std::endl;
#endif /* IPO_DEBUG_LU */

      /* Cache size since we change it. */
      const std::size_t n = size();

      /* Subtract the dot product of the new row and new column from the new diagonal. */
      for (std::size_t i = 0; i < n; ++i)
        newDiagonal -= newRow[i] * newColumn[i];
#ifdef IPO_DEBUG_LU
      std::cout << "last diagonal = " << newDiagonal << std::endl;
#endif /* IPO_DEBUG_LU */
      assert(!_isZero(newDiagonal));

      /* Add the new row to L. */
      for (std::size_t i = 0; i < n; ++i)
      {
        if (_isZero(newRow[i]))
          continue;

        _leftRows[i].push_back(n);
        _leftEntries[i].push_back(newRow[i]);
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
        if (_isZero(newColumn[i]))
          continue;

        lastColumnRows.push_back(i);
        lastColumnEntries.push_back(newColumn[i]);
      }
      if (!_isZero(newDiagonal))
      {
        lastColumnRows.push_back(n);
        lastColumnEntries.push_back(newDiagonal);
      }

      assert(size() == n+1);

#ifdef IPO_DEBUG_LU
      std::vector<std::vector<T>> dump(size());
      for (std::size_t r = 0; r < size(); ++r)
        dump[r] = std::vector<T>(size(), 0);
      for (std::size_t c = 0; c < size(); ++c)
      {
        for (std::size_t p = 0; p < _leftRows[c].size(); ++p)
          dump[_leftRows[c][p]][c] = _leftEntries[c][p];
      }
      for (std::size_t c = 0; c < size(); ++c)
      {
        for (std::size_t p = 0; p < _upperRows[c].size(); ++p)
          dump[_upperRows[c][p]][c] = _upperEntries[c][p];
      }
      std::cout << "Combined L and U matrix:\n";
      for (std::size_t r = 0; r < size(); ++r)
      {
        for (std::size_t c = 0; c < size(); ++c)
          std::cout << " " << dump[r][c];
        std::cout << std::endl;
      }
#endif /* IPO_DEBUG_LU */
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
    const IsZero _isZero;
    std::vector<std::vector<std::size_t>> _leftRows;
    std::vector<std::vector<T>> _leftEntries;
    std::vector<std::vector<std::size_t>> _upperRows;
    std::vector<std::vector<T>> _upperEntries;
  };

#if defined(IPO_WITH_GMP)
  struct RationalIsZero
  {
    bool operator()(const mpq_class& x) const
    {
      return x == 0;
    }
  };
#endif /* IPO_WITH_GMP */

  struct RealIsZero
  {
    double epsilon;

    RealIsZero(double eps)
      : epsilon(eps)
    {

    }

    bool operator()(double x) const
    {
      return fabs(x) < epsilon;
    }
  };

} /* namespace ipo */
