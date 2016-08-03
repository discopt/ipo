#include "mixed_integer_program.h"

#include "spx_gmp.h"
#include "reconstruct.h"

#ifdef WITH_SCIP
#include <scip/cons_linear.h>
#include "scip_oracles.h"
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
  
  MIPOracleBase::MIPOracleBase(const std::string& name, const Space& space, const MixedIntegerProgram& mip)
    : OracleBase(name, space)
  {
    OracleBase::initializedSpace();
    
    initializeLP(mip);
  }

  MIPOracleBase::MIPOracleBase(const std::string& name, OracleBase* nextOracle, const MixedIntegerProgram& mip)
    : OracleBase(name, nextOracle)
  {
    OracleBase::initializedSpace();
    
    initializeLP(mip);
  }

  MIPOracleBase::~MIPOracleBase()
  {
    delete _objective;
    delete _spx;
  }
  
  void MIPOracleBase::initializeLP(const MixedIntegerProgram& mip)
  {
    assert(space().dimension() == mip.numColumns());

    _spx = new SoPlex;
    _spx->setIntParam(SoPlex::SOLVEMODE, SoPlex::SOLVEMODE_RATIONAL);
    _spx->setIntParam(SoPlex::SYNCMODE, SoPlex::SYNCMODE_AUTO);
    _spx->setRealParam(SoPlex::FEASTOL, 0.0);
    _spx->setBoolParam(SoPlex::RATREC, true);
    _spx->setBoolParam(SoPlex::RATFAC, true);
    _spx->setIntParam(SoPlex::OBJSENSE, SoPlex::OBJSENSE_MAXIMIZE);
    _spx->setIntParam(SoPlex::VERBOSITY, SoPlex::VERBOSITY_ERROR);
    _spx->setIntParam(SoPlex::SIMPLIFIER, SoPlex::SIMPLIFIER_OFF);
    _spx->addColsRational(mip.columns());
    _spx->addRowsRational(mip.rows());
    _numRows = _spx->numRowsRational();

    // Initialize column data.

    _columns.resize(mip.numColumns());
    for (std::size_t v = 0; v < _columns.size(); ++v)
    {
      _columns[v].integral = mip.isIntegral(v);
      _columns[v].upper = mip.columns().upper(v);
      _columns[v].lower = mip.columns().lower(v);
    }
    
    _objective = new double[space().dimension()];
    _lpResult.reDim(space().dimension());
  }

  std::size_t MIPOracleBase::maximizeImplementation(OracleResult& result, const VectorRational& objective,
    const ObjectiveBound& objectiveBound, std::size_t minHeuristic, std::size_t maxHeuristic, bool& sort, bool& checkDups)
  {
    std::size_t n = space().dimension();

    // Scale objective vector.
    
    Rational largest = 0;
    for (std::size_t v = 0; v < n; ++v)
    {
      if (objective[v] > 0)
      {
        if (objective[v] > largest)
          largest = objective[v];
      }
      else
      {
        if (-objective[v] > largest)
          largest = -objective[v];
      }
    }

    // Convert exact objective to floating-point one.

    double factor = 1.0;
    if (largest > 0 && (largest > 1024 || largest < 1))
    {
      factor = 1024.0 / double(largest);
    }
    for (std::size_t v = 0; v < n; ++v)
      _objective[v] = factor * objective[v];

    // Call the solver's method.
    
    assert(_points.empty());
    assert(_rays.empty());

    solverMaximize(_objective, infinity, _points, _rays);

    if (!_points.empty())
    {
      sort = true;
      checkDups = true;
      result.points.reserve(_points.size());
      prepareSolver(objective);

      for (std::size_t i = 0; i < _points.size(); ++i)
      {
        Rational objValue;
        DSVectorRational* pt = findPoint(_points[i], objValue);
        if (pt != NULL)
          result.points.push_back(OracleResult::Point(pt, objValue));
        delete[] _points[i];
      }

      restoreSolver();
      _points.clear();
    }
    else if (!_rays.empty())
    {
      for (std::size_t i = 0; i < _rays.size(); ++i)
      {
        if (_rays[i] != NULL)
          delete[] _rays[i];
      }
      _rays.clear();

      prepareSolver(objective);
      soplex::DSVectorRational* ray = findRay();
      restoreSolver();
      
      if (ray)
        result.directions.push_back(OracleResult::Direction(ray));
      else if (heuristicLevel() == 0)
        throw std::runtime_error("MIPOracle: Ray search failed.");

      // If no ray could be found and heuristicLevel is positive, then we force forwarding by claiming infeasibility.
    }
    return heuristicLevel();
  }

  void MIPOracleBase::separatePoint(const VectorRational& point, LPRowSetRational& cuts)
  {
    
  }
  
  void MIPOracleBase::separateRay(const VectorRational& ray, LPRowSetRational& cuts)
  {
    
  }

  void MIPOracleBase::setFace(Face* newFace)
  {
    OracleBase::setFace(newFace);
    
    throw std::runtime_error("MIPOracleBase::setFace not implemented, yet!");
  }

  void MIPOracleBase::prepareSolver(const VectorRational& objective)
  {
    _spx->changeObjRational(objective);
  }

  void MIPOracleBase::restoreSolver()
  {
    // Remove all added rows.

    _lpRowPermutation.resize(_spx->numRowsRational());
    for (std::size_t i = 0; i < _numRows; ++i)
      _lpRowPermutation[i] = 0;
    for (std::size_t i = _numRows; i < _spx->numRowsRational(); ++i)
      _lpRowPermutation[i] = -1;
    _spx->removeRowsRational(&_lpRowPermutation[0]);    
  }

  DSVectorRational* MIPOracleBase::findPoint(double* approxPoint, Rational& objectiveValue)
  {
    std::size_t n = space().dimension();
    
    // Fix variable bounds for the integer variables.

    for (std::size_t v = 0; v < n; ++v)
    {
      if (_columns[v].integral)
      {
        Rational value = Rational(int(approxPoint[v] + 0.5));
        _spx->changeBoundsRational(v, value, value);
      }
      else
      {
        assert(_columns[v].lower == _spx->lowerRational(v));
        assert(_columns[v].upper == _spx->upperRational(v));
      }
    }

    while (true)
    {
      SPxSolver::Status status = _spx->solve();
      
      if (status == SPxSolver::UNBOUNDED)
      {
        _spx->getPrimalRayRational(_lpResult);
        _separateResult.clear();
        separateRay(_lpResult, _separateResult);
        if (_separateResult.num() == 0)
          throw std::runtime_error("MIPOracle: Claim is bounded, but candidate ray is not separated.");

        _spx->addRowsRational(_separateResult);  
      }
      else if (status != SPxSolver::OPTIMAL)
        return NULL;

      _spx->getPrimalRational(_lpResult);
      _separateResult.clear();
      separatePoint(_lpResult, _separateResult);
      if (_separateResult.num() == 0)
        break;

      _spx->addRowsRational(_separateResult);
    }

    // Create point as sparse vector.

    DSVectorRational* point = new DSVectorRational;
    for (std::size_t v = 0; v < n; ++v)
    {
      if (_lpResult[v] != 0)
        point->add(v, _lpResult[v]);
    }
    point->sort();
    objectiveValue = _spx->objValueRational();
    return point;
  }

  DSVectorRational* MIPOracleBase::findRay()
  {
    std::size_t n = space().dimension();

    // Set original variable bounds for the integer variables.

    for (std::size_t v = 0; v < n; ++v)
    {
      if (_columns[v].integral)
      {
        _spx->changeBoundsRational(v, _columns[v].lower, _columns[v].upper);
      }
      else
      {
        assert(_columns[v].lower == _spx->lowerRational(v));
        assert(_columns[v].upper == _spx->upperRational(v));
      }
    }

    while (true)
    {
      SPxSolver::Status status = _spx->solve();
      if (status != SPxSolver::UNBOUNDED)
        return NULL;
     
      _spx->getPrimalRayRational(_lpResult);
      _separateResult.clear();
      separateRay(_lpResult, _separateResult);
      if (_separateResult.num() == 0)
        break;

      _spx->addRowsRational(_separateResult);
    }

    // Create ray as sparse vector.

    DSVectorRational* ray = new DSVectorRational;
    for (std::size_t v = 0; v < n; ++v)
    {
      if (_lpResult[v] != 0)
        ray->add(v, _lpResult[v]);
    }
    ray->sort();
    return ray;
  }

} /* namespace ipo */
