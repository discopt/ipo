#pragma once

#include <ipo/config.hpp>
#include <ipo/export.hpp>
#include <ipo/space.hpp>

#include <memory>

#ifdef IPO_WITH_SOPLEX

#define SOPLEX_WITH_GMP
#include <soplex.h>
#include <ipo/oracles.hpp>

namespace ipo
{

  class MakeRationalSolver
  {
  public:
    MakeRationalSolver(const std::vector<bool>& integrality,
      const std::vector<Constraint>& constraints);

    ~MakeRationalSolver();

    void addRow(const Constraint& constraint);

    void removeLastRow();

    void setObjective(const double* objectiveVector);
    
    void setObjective(const mpq_class* objectiveVector);

    void solve(OptimizationOracle::Result& result);

  private:
    soplex::SoPlex _spx;
    int* _indices;
    mpq_t* _coefficients;
    mpq_t* _originalLowerBounds;
    mpq_t* _originalUpperBounds;
    std::vector<bool> _integrality;
  };

} /* namespace ipo */

#endif /* IPO_WITH_SOPLEX */
