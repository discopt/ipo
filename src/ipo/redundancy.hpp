#pragma once

// #define IPO_DEBUG_REDUNDANCY

#include <ipo/config.hpp>
#include <ipo/constraint.hpp>
#include "lu.hpp"

namespace ipo
{
  enum EquationRedundancy
  {
    EQUATION_INDEPENDENT,
    EQUATION_REDUNDANT,
    EQUATION_INCONSISTENT,
    EQUATION_INVALID
  };

  template <typename T>
  class EquationRedundancyCheck
  {
  public:
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

    EquationRedundancy test(const Constraint<T>& constraint, double epsilonConstraint,
      double epsilonCoefficient, double* pNorm = NULL) const
    {
      if (constraint.lhs() != constraint.rhs())
        EQUATION_INVALID;

      std::size_t newBasic;
      double norm = pNorm ? *pNorm : euclideanNorm(constraint.vector());
      T rhs = -constraint.rhs();
      EquationRedundancy result = testImplementation(constraint.vector(), rhs, newBasic,
        epsilonCoefficient);
      if (result == EQUATION_REDUNDANT && fabs(convertNumber<double>(rhs)) > epsilonConstraint * norm)
        return EQUATION_INCONSISTENT;
      else
        return result;
    }

    EquationRedundancy test(const sparse_vector<T>& vector, double epsilonCoefficient) const
    {
      T rhs(0);
      std::size_t newBasic;
      return testImplementation(vector, rhs, newBasic, epsilonCoefficient);
    }

    EquationRedundancy add(const Constraint<T>& constraint, double epsilonConstraint,
      double epsilonFactorization, double epsilonCoefficient, double* pNorm = NULL)
    {
#if defined(IPO_DEBUG_REDUNDANCY)
      std::cout << "EquationRedundancy.add(" << constraint << ")." << std::endl;
#endif /* IPO_DEBUG_REDUNDANCY */

      if (constraint.lhs() != constraint.rhs())
        return EQUATION_INVALID;

      double norm = pNorm ? *pNorm : euclideanNorm(constraint.vector());
      std::size_t newBasic;
      T rhs = -constraint.rhs();
      EquationRedundancy result = testImplementation(constraint.vector(), rhs, newBasic,
        epsilonCoefficient);
      if (result == EQUATION_REDUNDANT && fabs(convertNumber<double>(rhs)) > epsilonConstraint * norm)
        return EQUATION_INCONSISTENT;
      else if (result != EQUATION_INDEPENDENT)
        return result;

      /* Add the constraint to the equations. */
      _equations.push_back(constraint);

      /* We enlarge the basis */
      std::vector<T> basisRow(_basis.size(), 0);
      std::vector<T> basisColumn(_basis.size() + 1, 0);
      _basis.push_back(newBasic);
      for (std::size_t b = 0; b + 1 < _basis.size(); ++b)
      {
        basisRow[b] = _equations[b].vector().find(newBasic, 0);
      }
      for (std::size_t b = 0; b < _basis.size(); ++b)
      {
        basisColumn[b] = constraint.vector().find(_basis[b], 0);
      }

      _lu.extend(&basisRow[0], &basisColumn[0], basisColumn.back(), epsilonFactorization);

      return result;
    }

  protected:

    EquationRedundancy testImplementation(const sparse_vector<T>& vector, T& rhs,
      std::size_t& newBasic, double epsilonCoefficient) const
    {
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
        rhs += multipliers[e] * _equations[e].rhs();
      }
      /* Subtract given equation normal. */
      for (const auto& iter : vector)
        dense[iter.first] -= iter.second;

#if defined(IPO_DEBUG_REDUNDANCY)
      for (std::size_t v = 0; v < numVariables(); ++v)
        std::cout << " " << dense[v];
      std::cout << "; " << rhs << std::endl;
#endif /* IPO_DEBUG_REDUNDANCY */

      newBasic = std::numeric_limits<std::size_t>::max();
      double newBasicValue = 0.0;
      for (std::size_t v = 0; v < numVariables(); ++v)
      {
        if (fabs(convertNumber<double>(dense[v])) > epsilonCoefficient)
        {
          double value = fabs(convertNumber<double>(dense[v]));
          if (value > newBasicValue)
          {
            newBasic = v;
            newBasicValue = value;
          }
        }
      }
      if (newBasic != std::numeric_limits<std::size_t>::max())
        return EQUATION_INDEPENDENT;
      else
        return EQUATION_REDUNDANT;
    }

  protected:
    std::size_t _numVariables;
    IncrementalLUFactorization<T> _lu;
    std::vector<std::size_t> _basis;
    std::vector<Constraint<T>> _equations;
  };

} /* namespace ipo */


