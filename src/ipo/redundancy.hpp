#pragma once

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
  
  template <typename T, typename IsZero>
  class EquationRedundancyCheck
  {
  public:
    EquationRedundancyCheck(std::size_t numVariables, const IsZero isZero = IsZero())
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

    inline Constraint<T>& getEquation(std::size_t e)
    {
      return _equations[e];
    }

    EquationRedundancy test(const Constraint<T>& constraint) const
    {
      if (!_isZero(constraint.lhs() - constraint.rhs()))
        EQUATION_INVALID;

      std::size_t newBasic;
      return test(constraint, newBasic);
    }

    EquationRedundancy add(const Constraint<T>& constraint)
    {
#if defined(IPO_DEBUG_REDUNDANCY)
      std::cout << "EquationRedundancy.add(" << constraint.lhs.real << " <= [";
      for (std::size_t v = 0; v < numVariables(); ++v)
        std::cout << (v > 0 ? " " : "") << constraint.vector.findReal(v);
      std::cout << "] <= " << constraint.rhs.real << std::endl;
#endif /* IPO_DEBUG_REDUNDANCY */

      if (!_isZero(constraint.lhs() - constraint.rhs()))
        return EQUATION_INVALID;

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
      {
        basisRow[b] = _equations[b].vector().find(newBasic, 0);
      }
      for (std::size_t b = 0; b < _basis.size(); ++b)
      {
        basisColumn[b] = constraint.vector().find(_basis[b], 0);
      }

      _lu.extend(&basisRow[0], &basisColumn[0], basisColumn.back());

      return result;
    }

  protected:
//     inline bool isZero(const Value& value, double dummy) const
//     {
//       return _isZero(value.real);
//     }
// 
//     inline bool isZero(const Value& value, const mpq_class& dummy) const
//     {
//       return _isZero(value.rational);
//     }

    EquationRedundancy test(const Constraint<T>& constraint, std::size_t& newBasic) const
    {
      /* Find multipliers such that combined equation agrees with given one on basic variables.
       * To this end, solve Ax = r, where r is the basic part of the constraint vector. */
      std::vector<T> multipliers(_basis.size());
      for (std::size_t b = 0; b < _basis.size(); ++b)
        multipliers[b] = constraint.vector().find(_basis[b], 0);
      _lu.solveLeft(&multipliers[0]);
      _lu.solveUpper(&multipliers[0]);

      /* Combine equations according to multipliers. */
      std::vector<T> dense(numVariables(), 0);
      T rhs = 0;
      for (std::size_t e = 0; e < _equations.size(); ++e)
      {
        for (const auto& iter : _equations[e].vector())
          dense[iter.first] += multipliers[e] * iter.second;
        rhs += multipliers[e] * _equations[e].rhs();
      }
      /* Subtract given equation normal. */
      for (const auto& iter : constraint.vector())
        dense[iter.first] -= iter.second;
      rhs -= constraint.rhs();

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
      return _isZero(rhs) ? EQUATION_REDUNDANT : EQUATION_INCONSISTENT;
    }

  protected:
    std::size_t _numVariables;
    IsZero _isZero;
    IncrementalLUFactorization<T, IsZero> _lu;
    std::vector<std::size_t> _basis;
    std::vector<Constraint<T>> _equations;
  };

} /* namespace ipo */


