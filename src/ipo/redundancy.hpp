#pragma once

#include <ipo/config.hpp>
#include <ipo/data.hpp>
#include "lu.hpp"

namespace ipo
{
  enum EquationRedundancy
  {
    EQUATION_INDEPENDENT,
    EQUATION_REDUNDANT,
    EQUATION_INFEASIBLE,
    EQUATION_INVALID
  };
  
  template <typename T, typename IsZero>
  class EquationRedundancyCheck
  {
  public:
    EquationRedundancyCheck(std::size_t numVariables, const IsZero isZero)
      : _numVariables(numVariables), _isZero(isZero), _lu(isZero)
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
    
    EquationRedundancy test(const Constraint& constraint) const
    {
      if (!isZero(constraint.lhs - constraint.rhs, T(0)))
        EQUATION_INVALID;

      std::size_t newBasic;
      return test(constraint, newBasic);
    }

    EquationRedundancy add(const Constraint& constraint)
    {
#if defined(IPO_DEBUG_REDUNDANCY)
      std::cout << "EquationRedundancy.add(" << constraint.lhs.real << " <= [";
      for (std::size_t v = 0; v < numVariables(); ++v)
        std::cout << (v > 0 ? " " : "") << constraint.vector.findReal(v);
      std::cout << "] <= " << constraint.rhs.real << std::endl;
#endif /* IPO_DEBUG_REDUNDANCY */

      if (!isZero(constraint.lhs - constraint.rhs, T(0)))
        EQUATION_INVALID;

      std::size_t newBasic;
      EquationRedundancy result = test(constraint, newBasic);
      if (result != EQUATION_INDEPENDENT)
        return result;

      /* Add the constraint to the equations. */
      _equations.push_back(constraint);

      /* We enlarge the basis */
      std::vector<T> basisRow(_basis.size(), 0);
      std::vector<T> basisColumn(_basis.size() + 1, 0);
      _basis.push_back(newBasic);
      for (std::size_t b = 0; b + 1 < _basis.size(); ++b)
        _equations[b].vector.findEntry(newBasic, basisRow[b]);
      for (std::size_t b = 0; b < _basis.size(); ++b)
        constraint.vector.findEntry(_basis[b], basisColumn[b]);

      _lu.extend(&basisRow[0], &basisColumn[0], basisColumn.back());

      return result;
    }

  protected:
    inline bool isZero(const Value& value, double dummy) const
    {
      return _isZero(value.real);
    }

    inline bool isZero(const Value& value, const mpq_class& dummy) const
    {
      return _isZero(value.rational);
    }

    EquationRedundancy test(const Constraint& constraint, std::size_t& newBasic) const
    {
      /* Find multipliers such that combined equation agrees with given one on basic variables.
       * To this end, solve Ax = r, where r is the basic part of the constraint vector. */
      std::vector<T> multipliers(_basis.size());
      for (std::size_t b = 0; b < _basis.size(); ++b)
        constraint.vector.get(_basis[b], multipliers[b]);
      _lu.solveLeft(&multipliers[0]);
      _lu.solveUpper(&multipliers[0]);

      /* Combine equations according to multipliers. */
      std::vector<T> dense(numVariables(), 0);
      T rhs = 0;
      T x;
      for (std::size_t e = 0; e < _equations.size(); ++e)
      {
        const Vector& vector = _equations[e].vector;
        for (std::size_t i = 0; i < vector.size(); ++i)
        {
          vector.get(i, x);
          dense[vector.coordinate(i)] += multipliers[e] * x;
        }
        _equations[e].rhs.get(x);
        rhs += multipliers[e] * x;
      }
      /* Subtract given equation normal. */
      for (std::size_t i = 0; i < constraint.vector.size(); ++i)
      {
        constraint.vector.get(i, x);
        dense[constraint.vector.coordinate(i)] -= x;
      }
      constraint.rhs.get(x);
      rhs -= x;

#if defined(IPO_DEBUG_REDUNDANCY)
      for (std::size_t v = 0; v < numVariables(); ++v)
        std::cout << " " << dense[v];
      std::cout << "; " << rhs << std::endl;
#endif /* IPO_DEBUG_REDUNDANCY */

      for (std::size_t v = 0; v < numVariables(); ++v)
      {
        if (!_isZero(dense[v]))
        {
          newBasic = v;
          return EQUATION_INDEPENDENT;
        }
      }
      newBasic = std::numeric_limits<std::size_t>::max();
      return _isZero(rhs) ? EQUATION_REDUNDANT : EQUATION_INFEASIBLE;
    }

  protected:
    std::size_t _numVariables;
    IsZero _isZero;
    IncrementalLUFactorization<T, IsZero> _lu;
    std::vector<std::size_t> _basis;
    std::vector<Constraint> _equations;
  };

} /* namespace ipo */


