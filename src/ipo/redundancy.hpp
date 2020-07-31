#pragma once

// #define IPO_DEBUG_REDUNDANCY_PRINT

#include <ipo/config.hpp>
#include <ipo/constraint.hpp>
#include "lu.hpp"

namespace ipo
{
  template <typename T>
  class EquationRedundancyCheck
  {
  public:
    struct Result
    {
      double maxViolation;
      std::size_t maxCoordinate;
      T rhs;
    };

    EquationRedundancyCheck(std::size_t numVariables)
      : _numVariables(numVariables), _lu()
    {

    }

    inline std::size_t numVariables() const
    {
      return _numVariables;
    }

    inline std::size_t rank() const
    {
      return _basis.size();
    }

    inline Constraint<T>& getEquation(std::size_t e)
    {
      return _equations[e];
    }

    Result test(const Constraint<T>& constraint, double norm) const
    {
      assert(constraint.type() == EQUATION);

      return test(constraint.vector(), norm, -constraint.rhs());
    }

    bool add(const Constraint<T>& constraint, std::size_t maxCoordinate,
      double epsilonFactorization)
    {
      /* Add the constraint to the equations. */
      _equations.push_back(constraint);

      /* We enlarge the basis */
      std::vector<T> basisRow(_basis.size(), 0);
      std::vector<T> basisColumn(_basis.size() + 1, 0);
      _basis.push_back(maxCoordinate);
      for (std::size_t b = 0; b + 1 < _basis.size(); ++b)
      {
        basisRow[b] = _equations[b].vector().find(maxCoordinate, 0);
      }
      for (std::size_t b = 0; b < _basis.size(); ++b)
      {
        basisColumn[b] = constraint.vector().find(_basis[b], 0);
      }

      bool luResult = _lu.extend(&basisRow[0], &basisColumn[0], basisColumn.back(),
        epsilonFactorization);
      if (luResult)
        return true;
      else
      {
        _basis.pop_back();
        _equations.pop_back();
        return false;
      }
    }

    Result test(const sparse_vector<T>& vector, double norm, const T& rhs = 0) const
    {
#if defined(IPO_DEBUG_REDUNDANCY_PRINT)
      std::cout << "EquationRedundancy.test(" << vector << ". " << rhs << ")." << std::endl;
#endif /* IPO_DEBUG_REDUNDANCY_PRINT */

      Result result;
      result.rhs = rhs;

      /* Find multipliers such that combined equation agrees with given one on basic variables.
       * To this end, solve Ax = r, where r is the basic part of the constraint vector. */
      std::vector<T> multipliers(_basis.size());
      for (std::size_t b = 0; b < _basis.size(); ++b)
        multipliers[b] = vector.find(_basis[b], 0);
      _lu.solveLeft(&multipliers[0]);
      _lu.solveUpper(&multipliers[0]);

      /* Combine equations according to multipliers. */
      std::vector<T> dense(numVariables(), 0);
      for (std::size_t e = 0; e < _equations.size(); ++e)
      {
        for (const auto& iter : _equations[e].vector())
          dense[iter.first] += multipliers[e] * iter.second;
        result.rhs += multipliers[e] * _equations[e].rhs();
      }
      /* Subtract given equation normal. */
      for (const auto& iter : vector)
        dense[iter.first] -= iter.second;
      for (std::size_t b = 0; b < _basis.size(); ++b)
        dense[_basis[b]] = 0;

      result.maxCoordinate = std::numeric_limits<std::size_t>::max();
      result.maxViolation = 0.0;
      for (std::size_t v = 0; v < numVariables(); ++v)
      {
        double violation = fabs(convertNumber<double>(dense[v])) / norm;
        if (violation > result.maxViolation)
        {
          result.maxCoordinate = v;
          result.maxViolation = violation;
        }
      }
      return result;
    }

  protected:
    std::size_t _numVariables;
    IncrementalLUFactorization<T> _lu;
    std::vector<std::size_t> _basis;
    std::vector<Constraint<T>> _equations;
  };

} /* namespace ipo */


