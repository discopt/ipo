#include "lp.h"

#include "reconstruct.h"
#include "scip_oracle.h"
#include "scip_exception.h"

namespace ipo {

  LinearSet::LinearSet(const std::vector<Rational>& lowerBounds, const std::vector<Rational>& upperBounds, 
    const std::vector<LinearConstraint>& rows, const std::vector<std::string>& variableNames,
    const std::vector<std::string>& rowNames)
    : _rowNames(rowNames)
  {
    std::size_t n = lowerBounds.size();
    assert(lowerBounds.size() == n);
    assert(variableNames.empty() || (variableNames.size() == n));
    assert(rowNames.empty() || (rowNames.size() == rows.size()));

    _solver.setIntParam(soplex::SoPlex::SOLVEMODE, soplex::SoPlex::SOLVEMODE_RATIONAL);
    _solver.setIntParam(soplex::SoPlex::SYNCMODE, soplex::SoPlex::SYNCMODE_AUTO);
    _solver.setRealParam(soplex::SoPlex::FEASTOL, 0.0);
    _solver.setBoolParam(soplex::SoPlex::RATREC, true);
    _solver.setBoolParam(soplex::SoPlex::RATFAC, true);
    _solver.setIntParam(soplex::SoPlex::OBJSENSE, soplex::SoPlex::OBJSENSE_MAXIMIZE);
    _solver.setIntParam(soplex::SoPlex::VERBOSITY, soplex::SoPlex::VERBOSITY_ERROR);

    // Setup columns.

    soplex::DSVectorRational zero;
    soplex::LPColSetRational cols(n);
    for (std::size_t v = 0; v < n; ++v)
      cols.add(Rational(0), lowerBounds[v], zero, upperBounds[v]);
    _solver.addColsRational(cols);

    // Setup space.
    
    SpaceData* spaceData = new SpaceData();
    for (std::size_t v = 0; v < n; ++v)
    {
      if (variableNames.empty())
      {
        char buffer[256];
        snprintf(buffer, 256, "var#%ld", v);
//         std::stringstream str;
//         str << "var#" << v << std::endl;
        spaceData->addVariable(buffer);
      }
      else
      {
        spaceData->addVariable(variableNames[v]);
      }
    }
    _space = Space(spaceData);

    // Setup rows.

    addToLP(_solver, rows);

    if (_rowNames.empty())
    {
      _rowNames.resize(_solver.numRowsRational());
      for (std::size_t r = 0; r < _solver.numRowsRational(); ++r)
      {
        char buffer[256];
        snprintf(buffer, 256, "row#%ld", r);
        _rowNames[r] = buffer;
      }
    }
  }

  LinearSet::LinearSet(const LinearSet& other)
    : _space(other.space()), _rowNames(other._rowNames)
  {
    _solver.setIntParam(soplex::SoPlex::SOLVEMODE, soplex::SoPlex::SOLVEMODE_RATIONAL);
    _solver.setIntParam(soplex::SoPlex::SYNCMODE, soplex::SoPlex::SYNCMODE_AUTO);
    _solver.setRealParam(soplex::SoPlex::FEASTOL, 0.0);
    _solver.setBoolParam(soplex::SoPlex::RATREC, true);
    _solver.setBoolParam(soplex::SoPlex::RATFAC, true);
    _solver.setIntParam(soplex::SoPlex::OBJSENSE, soplex::SoPlex::OBJSENSE_MAXIMIZE);
    _solver.setIntParam(soplex::SoPlex::VERBOSITY, soplex::SoPlex::VERBOSITY_ERROR);

    // Setup columns.
    
    std::size_t n = other._solver.numColsRational();
    soplex::LPColSetRational cols(n);
    other._solver.getColsRational(0, n - 1, cols);
    // Clear column vectors - they are added later.
    for (std::size_t c = 0; c < n; ++c)
      cols.colVector_w(c).clear();
    _solver.addColsRational(cols);

    // Setup rows.

    std::size_t m = other._solver.numRowsRational();
    soplex::LPRowSetRational rows(m);
    other._solver.getRowsRational(0, m - 1, rows);
    _solver.addRowsRational(rows);

    assert(m == _rowNames.size());
  }

#ifdef IPO_WITH_SCIP
  
  LinearSet::LinearSet(SCIP* scip)
  {
    if (scip != NULL)
      constructLinearSet(scip);
  }

  LinearSet::LinearSet(const std::string& fileName)
  {
    SCIP* scip = readSCIP(fileName);
    constructLinearSet(scip);
    freeSCIP(scip);
  }
  
  void LinearSet::constructLinearSet(SCIP* scip)
  {
    _solver.setIntParam(soplex::SoPlex::SOLVEMODE, soplex::SoPlex::SOLVEMODE_RATIONAL);
    _solver.setIntParam(soplex::SoPlex::SYNCMODE, soplex::SoPlex::SYNCMODE_AUTO);
    _solver.setRealParam(soplex::SoPlex::FEASTOL, 0.0);
    _solver.setBoolParam(soplex::SoPlex::RATREC, true);
    _solver.setBoolParam(soplex::SoPlex::RATFAC, true);
    _solver.setIntParam(soplex::SoPlex::OBJSENSE, soplex::SoPlex::OBJSENSE_MAXIMIZE);
    _solver.setIntParam(soplex::SoPlex::VERBOSITY, soplex::SoPlex::VERBOSITY_ERROR);

    std::size_t n = SCIPgetNOrigVars(scip);
    SCIP_VAR** origVars = SCIPgetOrigVars(scip);

    // Create variable mapping.

    SCIPvarToIndexMap varToIndexMap;
    getSCIPvarToIndexMap(scip, varToIndexMap);

    // Setup variables.

    soplex::LPColSetRational cols(n);
    soplex::DSVectorRational zero;
    SpaceData* spaceData = new SpaceData();
    for (std::size_t v = 0; v < n; ++v)
    {
      // Lower bound.

      double realLowerBound = SCIPvarGetLbGlobal(origVars[v]);
      Rational lowerBound;
      if (realLowerBound > -SCIPinfinity(scip))
        reconstruct(realLowerBound, lowerBound, SCIPfeastol(scip));
      else
        lowerBound = -soplex::infinity;

      // Upper bound.

      double realUpperBound = SCIPvarGetUbGlobal(origVars[v]);
      Rational upperBound;
      if (realUpperBound < SCIPinfinity(scip))
        reconstruct(realUpperBound, upperBound, SCIPfeastol(scip));
      else
        upperBound = soplex::infinity;

      cols.add(Rational(0), lowerBound, zero, upperBound);

      // Variable name.

      spaceData->addVariable(SCIPvarGetName(origVars[v]));
    }
    _solver.addColsRational(cols);
    _space = Space(spaceData);

    // Setup row constraints.

    std::size_t m = SCIPgetNConss(scip);
    soplex::LPRowSetRational rows(m);
    if (m == 0)
      return;

    SCIP_CONS** constraints = SCIPgetConss(scip);
    _rowNames.reserve(2 * m);
    soplex::DSVectorRational normal;
    for (std::size_t i = 0; i < m; ++i)
    {
      SCIP_CONS* cons = constraints[i];
      if (std::string(SCIPconshdlrGetName(SCIPconsGetHdlr(cons))) != "linear")
        continue;

      double realLhs = SCIPgetLhsLinear(scip, cons);
      Rational lhs;
      bool hasLhs = realLhs > -SCIPinfinity(scip);
      if (hasLhs)
        reconstruct(realLhs, lhs, SCIPfeastol(scip));
      else
          lhs = -soplex::infinity;
      double realRhs = SCIPgetRhsLinear(scip, cons);
      Rational rhs = 0;
      bool hasRhs = realRhs < SCIPinfinity(scip);
      if (hasRhs)
        reconstruct(realRhs, rhs, SCIPfeastol(scip));
      else rhs = soplex::infinity;
      int nvars = SCIPgetNVarsLinear(scip, cons);
      SCIP_VAR** consVars = SCIPgetVarsLinear(scip, cons);
      double* vals = SCIPgetValsLinear(scip, cons);

      // Create normal vector.

      normal.clear();
      for (int j = 0; j < nvars; ++j)
      {
        Rational x;
        reconstruct(vals[j], x, SCIPfeastol(scip));
        if (x != 0)
          normal.add(varToIndexMap[consVars[j]], x);
      }

      // Create <= inequality or equation.

      if ((hasRhs && !hasLhs) || (!hasRhs && hasLhs) || (lhs == rhs))
      {
        rows.add(lhs, normal, rhs);
        _rowNames.push_back(SCIPconsGetName(cons));
      }
      else if (!hasLhs && !hasRhs)
      {
        std::cerr << "Warning: Ignoring SCIP constraint " << SCIPconsGetName(cons) << " with two infinite sides." << std::endl;
      }
      else
      {
        rows.add(-soplex::infinity, normal, rhs);
        _rowNames.push_back(std::string(SCIPconsGetName(cons)) + "_rhs");
        rows.add(lhs, normal, soplex::infinity);
        _rowNames.push_back(std::string(SCIPconsGetName(cons)) + "_lhs");
      }
    }

    _solver.addRowsRational(rows);
  }

  SCIP* LinearSet::readSCIP(const std::string& fileName)
  {
    SCIP* scip = NULL;
    SCIP_CALL_EXC(SCIPcreate(&scip));
    SCIP_CALL_EXC(SCIPincludeDefaultPlugins(scip));
    SCIP_CALL_EXC(SCIPsetIntParam(scip, "display/verblevel", 0));
    SCIP_CALL_EXC(SCIPreadProb(scip, fileName.c_str(), NULL));
    SCIP_CALL_EXC(SCIPtransformProb(scip));

    return scip;
  }

  void LinearSet::freeSCIP(SCIP* scip)
  {
    SCIP_CALL_EXC(SCIPfree(&scip));
  }


#endif /* IPO_WITH_SCIP */

  LinearSet::~LinearSet()
  {

  }

  LinearConstraint LinearSet::lowerBoundConstraint(std::size_t variableIndex) const
  {
    VectorData* data = new VectorData(1);
    data->add(variableIndex, Rational(1));
    const Rational& lower = lowerBound(variableIndex);
    const Rational& upper = upperBound(variableIndex);
    return LinearConstraint(lower == upper ? '=' : '>', Vector(data), lower);
  }

  LinearConstraint LinearSet::upperBoundConstraint(std::size_t variableIndex) const
  {
    VectorData* data = new VectorData(1);
    data->add(variableIndex, Rational(1));
    const Rational& lower = lowerBound(variableIndex);
    const Rational& upper = upperBound(variableIndex);
    return LinearConstraint(lower == upper ? '=' : '<', Vector(data), upper);
  }

  LinearConstraint LinearSet::rowConstraint(std::size_t rowIndex) const
  {
    assert(rowIndex < _solver.numRowsRational());

    Vector normal = sparseToVector(_solver.rowVectorRational(rowIndex));
    const Rational& lhs = _solver.lhsRational(rowIndex);
    const Rational& rhs = _solver.rhsRational(rowIndex);
    soplex::LPRowBase<Rational>::Type type = _solver.rowTypeRational(rowIndex);
    if (type == soplex::LPRowBase<Rational>::EQUAL)
      return LinearConstraint('=', normal, rhs);
    else if (type == soplex::LPRowBase<Rational>::LESS_EQUAL)
      return LinearConstraint('<', normal, rhs);
    else if (type == soplex::LPRowBase<Rational>::GREATER_EQUAL)
      return LinearConstraint('>', normal, lhs);
    else
      assert(false);
  }

  void LinearSet::getConstraints(std::vector<LinearConstraint>& constraints, bool includeBounds, bool includeRows,
    bool excludeEquations, bool excludeInequalities) const
  {
    if (includeBounds)
    {
      for (std::size_t v = 0; v < numVariables(); ++v)
      {
        const Rational& lower = lowerBound(v);
        const Rational& upper = upperBound(v);
        if (excludeEquations && lower == upper)
          continue;
        if (excludeInequalities && lower < upper)
          continue;
        if (lower > -soplex::infinity)
          constraints.push_back(lowerBoundConstraint(v));
        if (upper < soplex::infinity)
          constraints.push_back(upperBoundConstraint(v));
      }
    }

    if (includeRows)
    {
      for (std::size_t r = 0; r < numRows(); ++r)
      {
        if (excludeEquations && rowConstraint(r).isEquation())
          continue;
        if (excludeInequalities && !rowConstraint(r).isEquation())
          continue;
        constraints.push_back(rowConstraint(r));
      }
    }        
  }

  void LinearSet::changeLowerBound(std::size_t variableIndex, const Rational& newLower)
  {
    _solver.changeLowerRational(variableIndex, newLower);
  }

  void LinearSet::changeUpperBound(std::size_t variableIndex, const Rational& newUpper)
  {
    _solver.changeUpperRational(variableIndex, newUpper);
  }

  void LinearSet::changeBounds(std::size_t variableIndex, const Rational& newLower, const Rational& newUpper)
  {
    _solver.changeBoundsRational(variableIndex, newLower, newUpper);
  }

  std::size_t LinearSet::addConstraint(const LinearConstraint& constraint, const std::string& rowName)
  {
    addToLP(_solver, constraint);
    std::size_t rowIndex = _solver.numRowsRational() - 1;
    if (rowName.empty())
    {
      std::stringstream str;
      str << "row#" << rowIndex;
      _rowNames.push_back(str.str());
    }
    else
      _rowNames.push_back(rowName);
  }

  void LinearSet::removeConstraint(std::size_t rowIndex)
  {
    _solver.removeRowRational(rowIndex);
    _rowNames[rowIndex] = _rowNames.back();
    _rowNames.pop_back();
  }

  void LinearSet::removeLastConstraints(std::size_t newNumRows)
  {
    _solver.removeRowRangeRational(newNumRows + 1, _solver.numRowsRational() - 1);
    _rowNames.resize(newNumRows);
  }

  MixedIntegerLinearSet::MixedIntegerLinearSet(const std::vector<bool>& integrality, const std::vector<Rational>& lowerBounds,   
    const std::vector<Rational>& upperBounds, const std::vector<LinearConstraint>& rows,
    const std::vector<std::string>& variableNames, const std::vector<std::string>& rowNames)
    : LinearSet(lowerBounds, upperBounds, rows, variableNames, rowNames), _integrality(integrality)
  {
    assert(integrality.size() == lowerBounds.size());
  }

  MixedIntegerLinearSet::MixedIntegerLinearSet(const MixedIntegerLinearSet& other)
    : LinearSet(other), _integrality(other._integrality)
  {

  }

  MixedIntegerLinearSet::MixedIntegerLinearSet(const LinearSet& linearSet)
    : LinearSet(linearSet)
  {
    _integrality.resize(linearSet.numVariables(), false);
  }
  
#ifdef IPO_WITH_SCIP
  
  MixedIntegerLinearSet::MixedIntegerLinearSet(SCIP* scip)
    : LinearSet(scip)
  {
    constructIntegrality(scip);
  }

  MixedIntegerLinearSet::MixedIntegerLinearSet(const std::string& fileName)
    : LinearSet((SCIP*) NULL)
  {
    SCIP* scip = readSCIP(fileName);
    constructLinearSet(scip);
    constructIntegrality(scip);
    freeSCIP(scip);
  }

  void MixedIntegerLinearSet::constructIntegrality(SCIP *scip)
  {
    _integrality.resize(numVariables());
    SCIP_VAR** vars = SCIPgetOrigVars(scip);
    for (std::size_t v = 0; v < numVariables(); ++v)
      _integrality[v] = SCIPvarIsIntegral(vars[v]);
  }
  
#endif /* IPO_WITH_SCIP */

  MixedIntegerLinearSet::~MixedIntegerLinearSet()
  {

  }

  LinearProgram::LinearProgram(const std::vector<Rational>& objective, const std::vector<Rational>& lowerBounds,
    const std::vector<Rational>& upperBounds,
    const std::vector<LinearConstraint>& rows, const std::vector<std::string>& variableNames,
    const std::vector<std::string>& rowNames)
    : LinearSet(lowerBounds, upperBounds, rows, variableNames, rowNames)
  {
    changeObjective(objective);
  }

  LinearProgram::LinearProgram(const LinearProgram& linearProgram)
    : LinearSet(linearProgram)
  {
    soplex::DVectorRational objective(linearProgram.numVariables());
    linearProgram._solver.getObjRational(objective);
    _solver.changeObjRational(objective);
  }

  LinearProgram::LinearProgram(const LinearSet& linearSet)
    : LinearSet(linearSet)
  {

  }

#ifdef IPO_WITH_SCIP

  LinearProgram::LinearProgram(SCIP* scip)
    : LinearSet(scip)
  {
    constructObjective(scip);
  }

  LinearProgram::LinearProgram(const std::string& fileName)
    : LinearSet((SCIP*) NULL)
  {
    SCIP* scip = readSCIP(fileName);
    constructLinearSet(scip);
    constructObjective(scip);
    freeSCIP(scip);
  }

  void LinearProgram::constructObjective(SCIP* scip)
  {
    std::size_t n = SCIPgetNOrigVars(scip);
    SCIP_VAR** vars = SCIPgetOrigVars(scip);
    bool negate = SCIPgetObjsense(scip) == SCIP_OBJSENSE_MINIMIZE;
    for (std::size_t v = 0; v < n; ++v)
    {
      double value = SCIPvarGetObj(vars[v]);
      Rational x;
      if (SCIPisZero(scip, value))
        x = 0;
      else
      {
        reconstruct(value, x, SCIPfeastol(scip));
        if (negate)
          x = -x;
      }
      _solver.changeObjRational(v, x);
    }
  }

#endif /* IPO_WITH_SCIP */

  LinearProgram::~LinearProgram()
  {

  }

  void LinearProgram::changeObjective(const std::vector<Rational>& objective)
  {
    for (std::size_t v = 0; v < numVariables(); ++v)
      _solver.changeObjRational(v, objective[v]);
  }
  
  const Rational LinearProgram::getObjective(std::size_t index) const
  {
    Rational result;
    _solver.getObjRational(index, result);
    return result;
  }

  LinearProgram::Result LinearProgram::solve(Vector& vector, Rational& objectiveValue)
  {
    soplex::SPxSolver::Status status;
    try
    {
      status = _solver.solve();    
    }
    catch (soplex::SPxException& e)
    {
      std::cerr << e.what() << std::endl;
      throw std::runtime_error("Error while solving a LinearProgram using SoPlex.");
    }
    if (status == soplex::SPxSolver::OPTIMAL)
    {
      soplex::DVectorRational primalSolution(numVariables());
      _solver.getPrimalRational(primalSolution);
      vector = denseToVector(primalSolution);
      objectiveValue = _solver.objValueRational();
      return OPTIMAL;
    }
    else if (status == soplex::SPxSolver::UNBOUNDED)
    {
      soplex::DVectorRational primalRay(numVariables());
      _solver.getPrimalRayRational(primalRay);
      vector = denseToVector(primalRay);
      objectiveValue = soplex::infinity;
      return UNBOUNDED;
    }
    else if (status == soplex::SPxSolver::INFEASIBLE)
    {
      objectiveValue = -soplex::infinity;
      return INFEASIBLE;
    }
    else
      throw std::runtime_error("Error while solving a LinearProgram using SoPlex.");
  }

  void LinearProgram::writeToFile(const std::string& fileName)
  {
    _solver.writeFileRational(fileName.c_str());
  }

} /* namespace ipo */
