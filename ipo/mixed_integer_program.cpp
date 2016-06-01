#include "mixed_integer_program.h"

#include "spx_gmp.h"
#include "reconstruct.h"

#ifdef WITH_SCIP
#include <scip/cons_linear.h>
#endif

using namespace soplex;

namespace ipo {

#ifdef WITH_SCIP

  MixedIntegerProgram::MixedIntegerProgram(Space& space, SCIP* scip) : _space(space),
    _currentFace(NULL)
  {
    std::size_t n = SCIPgetNOrigVars(scip);
    SCIP_VAR** origVars = SCIPgetOrigVars(scip);
    _worker.reDim(n, true);

    /// Create map

    SCIPvarToIndexMap varToIndexMap;
    getSCIPvarToIndexMap(scip, varToIndexMap);

    /// Setup columns.

    bool buildSpace = space.dimension() == 0;
    if (!buildSpace && space.dimension() != n)
    {
      throw std::runtime_error(
        "Dimension mismatch while constructing a MixedIntegerProgram for a non-empty space.");
    }

    _columns.reMax(n);
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
      if (buildSpace)
        space.addVariable(SCIPvarGetName(origVars[v]));
    }

    /// Setup rows.

    SCIP_CONSHDLR* linearHandler = SCIPfindConshdlr(scip, "linear");
    std::size_t numLinearConstraints = SCIPconshdlrGetNConss(linearHandler);
    if (numLinearConstraints > 0)
    {
      SCIP_CONS** linearConstraints = SCIPconshdlrGetConss(linearHandler);
      _rows.reMax(numLinearConstraints);
      _rowNames.resize(numLinearConstraints);
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
        _rowNames[i] = SCIPconsGetName(cons);
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

  bool MixedIntegerProgram::checkPointRows(const SVectorRational* point)
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
    return checkPointIntegral(point) && checkPointBounds(point) && checkPointRows(point);
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

  bool MixedIntegerProgram::checkRayRows(const SVectorRational* ray)
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
    return checkRayBounds(ray) && checkRayRows(ray);
  }

  void MixedIntegerProgram::setFace(Face* newFace)
  {
    if (newFace == _currentFace)
      return;

    if (_currentFace != NULL)
      _rows.remove(_rows.num() - 1);

    _currentFace = newFace;

    if (_currentFace != NULL)
      _rows.add(_currentFace->rhs(), _currentFace->sparseNormal(), _currentFace->rhs());
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
            names->push_back(rowName(r));
        }
        continue;
      }
      if (!inequalities)
        continue;

      if (upper < infinity)
      {
        rows.add(-infinity, _rows.rowVector(r), upper);
        if (names)
          names->push_back(rowName(r) + "-rhs");
      }
      if (lower > -infinity)
      {
        vector = _rows.rowVector(r);
        vector *= -1;
        rows.add(-infinity, vector, -lower);
        if (names)
          names->push_back(rowName(r) + "-lhs");
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
        names->push_back("fixed-" + space()[c]);
    }
  }

  MixedIntegerProgramCorrectorOracle::MixedIntegerProgramCorrectorOracle(const std::string& name,
    MixedIntegerProgram& mip, OracleBase* approximateOracle)
    : OracleBase(name, mip.space()), _mip(mip), _approximateOracle(approximateOracle),
    _denseVector(mip.space().dimension())
  {
    if (_mip.space() != _approximateOracle->space())
      throw std::runtime_error("Spaces of MixedIntegerProgram and approximate oracle differ.");

    OracleBase::initializedSpace();

    initializeLP();
  }

  MixedIntegerProgramCorrectorOracle::MixedIntegerProgramCorrectorOracle(const std::string& name,
    MixedIntegerProgram& mip, OracleBase* approximateOracle,
    OracleBase* nextOracle)
    : OracleBase(name, nextOracle), _mip(mip), _approximateOracle(approximateOracle),
    _denseVector(mip.space().dimension())
  {
    if (_mip.space() != _approximateOracle->space())
      throw std::runtime_error("Spaces of MixedIntegerProgram and approximate oracle differ.");
    if (_mip.space() != _nextOracle->space())
      throw std::runtime_error("Spaces of MixedIntegerProgram and next oracle differ.");

    OracleBase::initializedSpace();

    initializeLP();
  }

  MixedIntegerProgramCorrectorOracle::~MixedIntegerProgramCorrectorOracle()
  {
    delete _spx;
  }

  void MixedIntegerProgramCorrectorOracle::initializeLP()
  {
//     std::cerr << "Creating SoPlex instance for MixedIntegerProgramCorrectorOracle." << std::endl;
    _spx = new SoPlex;
    _spx->setIntParam(SoPlex::SOLVEMODE, SoPlex::SOLVEMODE_RATIONAL);
    _spx->setIntParam(SoPlex::SYNCMODE, SoPlex::SYNCMODE_AUTO);
    _spx->setRealParam(SoPlex::FEASTOL, 0.0);
    _spx->setBoolParam(SoPlex::RATREC, true);
    _spx->setBoolParam(SoPlex::RATFAC, true);
    _spx->setIntParam(SoPlex::OBJSENSE, SoPlex::OBJSENSE_MAXIMIZE);
    _spx->setIntParam(SoPlex::VERBOSITY, SoPlex::VERBOSITY_ERROR);
    _spx->setIntParam(SoPlex::SIMPLIFIER, SoPlex::SIMPLIFIER_OFF);
    _spx->addColsRational(_mip.columns());
    _spx->addRowsRational(_mip.rows());
  }


  DSVectorRational* MixedIntegerProgramCorrectorOracle::correctPoint(const SVectorRational* point,
      const VectorRational& objective)
  {
    _spx->clearBasis(); // TODO: This should not be necessary, but produced a bug!

    /// Fix integers to zero.

    for (std::size_t v = 0; v < space().dimension(); ++v)
    {
      if (!_mip.isIntegral(v))
      {
        assert(_spx->lowerRational(v) == _mip.columns().lower(v));
        assert(_spx->upperRational(v) == _mip.columns().upper(v));
        continue;
      }

      _spx->changeBoundsRational(v, Rational(0), Rational(0));
    }

    /// Fix integers to point's values for its nonzeros.

    for (int p = point->size() - 1; p >= 0; --p)
    {
      std::size_t v = point->index(p);
      const Rational& x = point->value(p);
      if (_mip.isIntegral(v))
      {
        _spx->changeBoundsRational(v, x, x);
      }
    }

    /// Set objective.

    _spx->changeObjRational(objective);

    /// Solve LP.

    SPxSolver::Status status = _spx->solve();
    if (status != SPxSolver::OPTIMAL)
    {
      std::stringstream error;
      error << "MIP corrector: Unexpected LP status ";
      error << status << " while correcting point.";
      throw std::runtime_error(error.str());
    }

    /// Extract solution.

    DSVectorRational* result = new DSVectorRational;
    _spx->getPrimalRational(_denseVector);
    (*result) = _denseVector;
    return result;
  }

  DSVectorRational* MixedIntegerProgramCorrectorOracle::correctDirection(
    const SVectorRational* direction, const VectorRational& objective)
  {
    /// Relax integers to original bounds.

    for (std::size_t v = 0; v < space().dimension(); ++v)
    {
      if (!_mip.isIntegral(v))
      {
        assert(_spx->lowerRational(v) == _mip.columns().lower(v));
        assert(_spx->upperRational(v) == _mip.columns().upper(v));
        continue;
      }
      _spx->changeBoundsRational(v, _mip.columns().lower(v), _mip.columns().upper(v));
    }

    /// Set objective.

    _spx->changeObjRational(objective);

    /// Solve LP.

    SPxSolver::Status status = _spx->solve();
    if (status != SPxSolver::UNBOUNDED)
    {
      std::stringstream error;
      error << "MIP corrector: Unexpected LP status ";
      error << status << " while correcting direction.";
      throw std::runtime_error(error.str());
    }

    /// Extract direction.

    DSVectorRational* result = new DSVectorRational;
    _spx->getPrimalRayRational(_denseVector);
    (*result) = _denseVector;
    return result;
  }

  void MixedIntegerProgramCorrectorOracle::setFace(Face* newFace)
  {
    if (newFace == currentFace())
      return;

    if (currentFace() != NULL)
      _spx->removeRowRational(_spx->numRowsRational() - 1);

    OracleBase::setFace(newFace);
    _approximateOracle->setFace(newFace);
    _mip.setFace(newFace);

    if (currentFace() != NULL)
    {
      _spx->addRowRational(
        LPRowRational(currentFace()->rhs(), currentFace()->sparseNormal(), currentFace()->rhs()));
    }
  }

  void MixedIntegerProgramCorrectorOracle::maximize(OracleResult& result,
    const VectorRational& objective, const ObjectiveBound& objectiveBound,
    std::size_t maxHeuristic, std::size_t minHeuristic)
  {
    assert((thisHeuristic() == 0 && _nextOracle == NULL)
      || thisHeuristic() > 0 && _nextOracle != NULL);

    // Forward call if requested.

    if (thisHeuristic() > maxHeuristic)
      return _nextOracle->maximize(result, objective, objectiveBound, maxHeuristic, minHeuristic);

    // Call approximate oracle with precisely its heuristic level to avoid forwarding there.

    OracleResult approximateResult;
    _approximateOracle->maximize(approximateResult, objective, objectiveBound,
      _approximateOracle->thisHeuristic(), _approximateOracle->thisHeuristic());

    if (approximateResult.isFeasible())
    {
      result.buildStart(objective);
      for (std::size_t i = 0; i < approximateResult.points.size(); ++i)
      {
        result.buildAddPoint(correctPoint(approximateResult.points[i].point, objective));
        delete approximateResult.points[i].point;
      }
      return result.buildFinish(thisHeuristic(), true, true, true);
    }
    else if (approximateResult.isUnbounded())
    {
      result.buildStart(objective);
      for (std::size_t i = 0; i < approximateResult.directions.size(); ++i)
      {
        result.buildAddDirection(correctDirection(approximateResult.directions[i].direction,
          objective));
        delete approximateResult.directions[i].direction;
      }
      return result.buildFinish(thisHeuristic(), false, false, true);
    }
    else
    {
      assert(approximateResult.isInfeasible());
      result.buildStart(objective);
      return result.buildFinish(thisHeuristic(), false, false, false);
    }
  }


} /* namespace ipo */
