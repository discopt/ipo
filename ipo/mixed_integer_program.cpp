#include "mixed_integer_program.h"

#include "spx_gmp.h"
#include "reconstruct.h"

#ifdef WITH_SCIP
#include <scip/cons_linear.h>
#endif

using namespace soplex;

namespace ipo {

#ifdef WITH_SCIP

  MixedIntegerProgram::MixedIntegerProgram(SCIP* scip) :
      _face(NULL)
  {
    std::size_t n = SCIPgetNOrigVars(scip);
    SCIP_VAR** origVars = SCIPgetOrigVars(scip);
    _worker.reDim(n, true);

    /// Create map

    SCIPvarToIndexMap varToIndexMap;
    getSCIPvarToIndexMap(scip, varToIndexMap);

    /// Setup columns.

    _columns.reMax(n);
    _variableNames.resize(n);
    _integrality.resize(n);
    Rational obj, lower, upper;
    DSVectorRational vector;
    for (std::size_t v = 0; v < n; ++v)
    {
      reconstruct(SCIPvarGetObj(origVars[v]), obj, SCIPdualfeastol(scip));
      double realLower = SCIPvarGetLbGlobal(origVars[v]);
      if (realLower > -SCIPinfinity(scip))
        reconstruct(realLower, lower, SCIPfeastol(scip));
      else
        lower = -infinity;
      double realUpper = SCIPvarGetUbGlobal(origVars[v]);
      if (realUpper < SCIPinfinity(scip))
        reconstruct(realUpper, upper, SCIPfeastol(scip));
      else
        upper = infinity;
      _columns.add(SCIPgetObjsense(scip) == SCIP_OBJSENSE_MAXIMIZE ? obj : -obj, lower, vector, upper);
      _integrality[v] = SCIPvarGetType(origVars[v]) != SCIP_VARTYPE_CONTINUOUS;
      _variableNames[v] = SCIPvarGetName(origVars[v]);
    }

    /// Setup rows.

    SCIP_CONSHDLR* linearHandler = SCIPfindConshdlr(scip, "linear");
    std::size_t numLinearConstraints = SCIPconshdlrGetNConss(linearHandler);
    if (numLinearConstraints > 0)
    {
      SCIP_CONS** linearConstraints = SCIPconshdlrGetConss(linearHandler);
      _rows.reMax(numLinearConstraints);
      _constraintNames.resize(numLinearConstraints);
      for (std::size_t i = 0; i < numLinearConstraints; ++i)
      {
        SCIP_CONS* cons = linearConstraints[i];
        double realLower = SCIPgetLhsLinear(scip, cons);
        if (realLower > -SCIPinfinity(scip))
          reconstruct(realLower, lower, SCIPfeastol(scip));
        else
          lower = -infinity;
        double realUpper = SCIPgetRhsLinear(scip, cons);
        if (realUpper < SCIPinfinity(scip))
          reconstruct(realUpper, upper, SCIPfeastol(scip));
        else
          upper = infinity;

        int nvars = SCIPgetNVarsLinear(scip, cons);
        SCIP_VAR** consVars = SCIPgetVarsLinear(scip, cons);
        double* vals = SCIPgetValsLinear(scip, cons);
        vector.clear();
        for (int j = 0; j < nvars; ++j)
        {
          double apxValue = vals[j];
          reconstruct(vals[j], obj, SCIPfeastol(scip));
          if (obj != 0)
            vector.add(varToIndexMap[consVars[j]], obj);
        }
        _rows.add(lower, vector, upper);
        _constraintNames[i] = SCIPconsGetName(cons);
      }
    }
  }

#endif

  MixedIntegerProgram::~MixedIntegerProgram()
  {

  }

  bool MixedIntegerProgram::checkPointBounds(const SVectorRational* point) const
  {
    for (int p = point->size() - 1; p >= 0; --p)
    {
      std::size_t v = point->index(p);
      const Rational& x = point->value(p);
      if (_columns.lower(v) > -infinity && x < _columns.lower(v))
        return false;
      if (_columns.upper(v) < infinity && x > _columns.upper(v))
        return false;
    }
    return true;
  }

  bool MixedIntegerProgram::checkPointConstraints(const SVectorRational* point)
  {
    _worker.clear();
    _worker.assign(*point);

    for (int r = 0; r < _rows.num(); ++r)
    {
      const Rational activity = _rows.rowVector(r) * _worker;
      if (_rows.lhs(r) > -infinity && activity < _rows.lhs(r))
        return false;
      if (_rows.rhs(r) < infinity && activity > _rows.rhs(r))
        return false;
    }
    return true;
  }

  bool MixedIntegerProgram::checkPointIntegral(const SVectorRational* point) const
  {
    for (int p = point->size() - 1; p >= 0; --p)
    {
      if (!_integrality[point->index(p)])
        continue;
      if (rational2mpzDen(point->value(p)) != 1)
        return false;
    }
    return true;
  }

  bool MixedIntegerProgram::checkPoint(const SVectorRational* point)
  {
    return checkPointIntegral(point) && checkPointBounds(point) && checkPointConstraints(point);
  }

  bool MixedIntegerProgram::checkRayBounds(const SVectorRational* ray) const
  {
    for (int p = ray->size() - 1; p >= 0; --p)
    {
      std::size_t v = ray->index(p);
      const Rational& x = ray->value(p);
      if (_columns.lower(v) > -infinity && x < 0)
        return false;
      if (_columns.upper(v) < infinity && x > 0)
        return false;
    }
    return true;
  }

  bool MixedIntegerProgram::checkRayConstraints(const SVectorRational* ray)
  {
    _worker.clear();
    _worker.assign(*ray);

    for (int r = 0; r < _rows.num(); ++r)
    {
      const Rational activity = _rows.rowVector(r) * _worker;
      if (_rows.lhs(r) > -infinity && activity < 0)
        return false;
      if (_rows.rhs(r) < infinity && activity > 0)
        return false;
    }
    return true;
  }

  bool MixedIntegerProgram::checkRay(const SVectorRational* ray)
  {
    return checkRayBounds(ray) && checkRayConstraints(ray);
  }

  void MixedIntegerProgram::faceEnabled(Face* face)
  {
    if (_face != NULL && _face != face)
    {
      throw std::runtime_error("Error while enabling two distinct faces for MixedIntegerProgram!");
    }
    else if (_face != face)
    {
      _face = face;
      _rows.add(face->rhs(), face->normal(), face->rhs());
    }
  }

  void MixedIntegerProgram::faceDisabled(Face* face)
  {
    if (_face == NULL)
      return;
    if (_face != face)
    {
      throw std::runtime_error("Error while disabling unknown face for MixedIntegerProgram!");
    }
    _face = NULL;
    _rows.remove(_rows.num() - 1);
  }

  void MixedIntegerProgram::getConstraints(LPRowSetRational& rows, bool inequalities, bool equations,
      std::vector<std::string>* names)
  {
    rows.reMax(_rows.max());
    DSVectorRational vector;
    for (int r = 0; r < _rows.num(); ++r)
    {
      const Rational& lower = _rows.lhs(r);
      const Rational& upper = _rows.rhs(r);
      if (lower == upper)
      {
        if (equations)
        {
          rows.add(lower, _rows.rowVector(r), upper);
          if (names)
            names->push_back(constraintName(r));
        }
        continue;
      }
      if (!inequalities)
        continue;

      if (upper < infinity)
      {
        rows.add(-infinity, _rows.rowVector(r), upper);
        if (names)
          names->push_back(constraintName(r) + "-rhs");
      }
      if (lower > -infinity)
      {
        vector = _rows.rowVector(r);
        vector *= -1;
        rows.add(-infinity, vector, -lower);
        if (names)
          names->push_back(constraintName(r) + "-lhs");
      }
    }
  }

  void MixedIntegerProgram::getFixedVariableEquations(LPRowSetRational& rows, std::vector<std::string>* names)
  {
    DSVectorRational vector;
    for (std::size_t c = 0; c < _columns.num(); ++c)
    {
      if (_columns.upper(c) != _columns.lower(c))
        continue;

      const Rational& value = _columns.upper(c);
      vector.clear();
      vector.add(c, Rational(1));
      rows.add(value, vector, value);
      if (names)
        names->push_back("fixed-" + variableName(c));
    }
  }

  MixedIntegerProgramCorrectorOracle::MixedIntegerProgramCorrectorOracle(const std::string& name,
      MixedIntegerProgram& mip, FaceOptimizationOracleBase* inexact, bool correctAlways) :
      FaceOptimizationOracleBase(name), _mip(mip), _inexact(inexact), _worker(mip.numVariables()), _correctAlways(
          correctAlways)
  {
    if (_mip.numVariables() != _inexact->numVariables())
      throw std::runtime_error("Dimensions of MixedIntegerProgram and inexact oracle differ!");
    initialize(_inexact);

    /// Setup LP for corrections.

    _spx.setIntParam(SoPlex::SOLVEMODE, SoPlex::SOLVEMODE_RATIONAL);
    _spx.setIntParam(SoPlex::SYNCMODE, SoPlex::SYNCMODE_AUTO);
    _spx.setRealParam(SoPlex::FEASTOL, 0.0);
    _spx.setBoolParam(SoPlex::RATREC, true);
    _spx.setBoolParam(SoPlex::RATFAC, true);
    _spx.setIntParam(SoPlex::OBJSENSE, SoPlex::OBJSENSE_MAXIMIZE);
    _spx.setIntParam(SoPlex::VERBOSITY, SoPlex::VERBOSITY_ERROR);
    _spx.setIntParam(SoPlex::SIMPLIFIER, SoPlex::SIMPLIFIER_OFF);
    _spx.addColsRational(_mip.columns());
    _spx.addRowsRational(_mip.rows());
  }

  MixedIntegerProgramCorrectorOracle::~MixedIntegerProgramCorrectorOracle()
  { 

  }

  DSVectorRational* MixedIntegerProgramCorrectorOracle::correctPoint(const SVectorRational* point,
      const VectorRational& objective)
  {
    if (!_correctAlways && _mip.checkPoint(point))
      return NULL;

    _spx.clearBasis(); // TODO: This should not be necessary, but produced a bug!
//     std::cerr << _spx.numRowsRational() << "x" << _spx.numColsRational() << std::endl;
//     _spx.writeFileRational("correct-point.lp");

    /// Fix integers to zero.

    for (std::size_t v = 0; v < numVariables(); ++v)
    {
      if (!_mip.isIntegral(v))
      {
        assert(_spx.lowerRational(v) == _mip.columns().lower(v));
        assert(_spx.upperRational(v) == _mip.columns().upper(v));
//        std::cerr << "Continuous var " << v << ": " << _spx.lowerRational(v) << " -- " << _spx.upperRational(v) << std::endl;
        continue;
      }
      _spx.changeBoundsRational(v, 0, 0);
    }

    /// Fix integers to point's values for its nonzeros.

    for (int p = point->size() - 1; p >= 0; --p)
    {
      std::size_t v = point->index(p);
      const Rational& x = point->value(p);
      if (_mip.isIntegral(v))
        _spx.changeBoundsRational(v, x, x);
    }

    /// Set objective.

    _spx.changeObjRational(objective);

    /// Solve LP.
    
//     std::cerr << "Solving" << std::endl;

    SPxSolver::Status status = _spx.solve();
    if (status != SPxSolver::OPTIMAL)
      throw std::runtime_error("MixedIntegerProgram: Unexpected LP status while correcting point!");

    /// Extract solution.

    DSVectorRational* result = new DSVectorRational;
    _spx.getPrimalRational(_worker);
    (*result) = _worker;
    return result;
  }

  DSVectorRational* MixedIntegerProgramCorrectorOracle::correctRay(const SVectorRational* ray,
      const VectorRational& objective)
  {
    if (!_correctAlways && _mip.checkRay(ray))
      return NULL;

    /// Relax integers to original bounds.

    for (std::size_t v = 0; v < numVariables(); ++v)
    {
      if (!_mip.isIntegral(v))
      {
        assert(_spx.lowerRational(v) == _mip.columns().lower(v));
        assert(_spx.upperRational(v) == _mip.columns().upper(v));
        continue;
      }
      _spx.changeBoundsRational(v, _mip.columns().lower(v), _mip.columns().upper(v));
    }

    /// Set objective.

    _spx.changeObjRational(objective);

    /// Solve LP.

    SPxSolver::Status status = _spx.solve();
    if (status != SPxSolver::UNBOUNDED)
      throw std::runtime_error("MixedIntegerProgram: Unexpected LP status while correcting ray!");

    /// Extract ray.

    DSVectorRational* result = new DSVectorRational;
    _spx.getPrimalRayRational(_worker);
    (*result) = _worker;
    return result;
  }

  void MixedIntegerProgramCorrectorOracle::run(OptimizationResult& result, const VectorRational& objective,
      const Rational* improveValue, bool forceOptimal)
  {
    /// Run inexact oracle.

    if (improveValue != NULL)
      _inexact->improve(result, objective, *improveValue, forceOptimal);
    else
      _inexact->maximize(result, objective, forceOptimal);

    bool correctedOptimal = false;
    if (!result.optimal)
    {
      /// If returned solution is not optimal, but all variables are continuous
      /// and we correct all solutions, then we are optimal.

      correctedOptimal = true;
      for (std::size_t v = 0; v < _mip.numVariables(); ++v)
      {
        if (_mip.isIntegral(v))
        {
          correctedOptimal = false;
          break;
        }
      }
    }

    if (result.isUnbounded())
    {
      /// Unbounded case.

      bool filterDuplicates = false;
      for (std::size_t r = 0; r < result.directions.size(); ++r)
      {
        DSVectorRational* ray = correctRay(result.directions[r], objective);
        if (ray == NULL)
          correctedOptimal = false;
        else
        {
          delete result.directions[r];
          result.directions[r] = ray;
          filterDuplicates = true;
        }
      }
      if (filterDuplicates)
        result.filterDuplicates();
    }
    else if (result.isFeasible())
    {
      /// Feasible case.

      bool filterDuplicates = false;
      for (std::size_t p = 0; p < result.points.size(); ++p)
      {
        DSVectorRational* point = correctPoint(result.points[p], objective);
        if (point == NULL)
          correctedOptimal = false;
        else
        {
          delete result.points[p];
          result.points[p] = point;
          filterDuplicates = true;
        }
      }

      if (filterDuplicates)
        result.filterDuplicates();

      /// Recompute objective values.

      result.setFeasible(objective);
    }

    if (correctedOptimal)
      result.optimal = true;
  }

  void MixedIntegerProgramCorrectorOracle::faceEnabled(Face* face)
  {
    _inexact->setFace(face);
    _mip.faceEnabled(face);
    _spx.addRowRational(LPRowRational(face->rhs(), face->normal(), face->rhs()));
  }

  void MixedIntegerProgramCorrectorOracle::faceDisabled(Face* face)
  {
    _inexact->setFace(NULL);
    _spx.removeRowRational(_spx.numRowsRational() - 1);
    _mip.faceDisabled(face);
  }

} /* namespace ipo */
