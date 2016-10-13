#include "mip.h"

#include "rational.h"
#include "reconstruct.h"

#ifdef WITH_SCIP
#include <scip/cons_linear.h>
#include "scip_oracles.h"
#endif

namespace ipo {
  
  MixedIntegerSet::Variable::Variable()
    : lowerBound(-soplex::infinity), upperBound(soplex::infinity), integral(false)
  {

  }

  MixedIntegerSet::Variable::Variable(const Rational& lower, const Rational& upper, bool integ)
    : lowerBound(lower), upperBound(upper), integral(integ)
  {

  }


#ifdef WITH_SCIP

  MixedIntegerSet::MixedIntegerSet(SCIP* scip) : _currentFace()
  {
    std::size_t n = SCIPgetNOrigVars(scip);
    SCIP_VAR** origVars = SCIPgetOrigVars(scip);
//     _worker.reDim(n, true);

    // Create map

    SCIPvarToIndexMap varToIndexMap;
    getSCIPvarToIndexMap(scip, varToIndexMap);

    // Setup variables.

    SpaceData* spaceData = new SpaceData();
    _variables.resize(n);
    for (std::size_t v = 0; v < n; ++v)
    {
      double realLower = SCIPvarGetLbGlobal(origVars[v]);
      if (realLower > -SCIPinfinity(scip))
        reconstruct(realLower, _variables[v].lowerBound, SCIPfeastol(scip));
      else
        _variables[v].lowerBound = -soplex::infinity;
      double realUpper = SCIPvarGetUbGlobal(origVars[v]);
      if (realUpper < SCIPinfinity(scip))
        reconstruct(realUpper, _variables[v].upperBound, SCIPfeastol(scip));
      else
        _variables[v].upperBound = soplex::infinity;
      _variables[v].integral = SCIPvarGetType(origVars[v]) != SCIP_VARTYPE_CONTINUOUS;

      spaceData->addVariable(SCIPvarGetName(origVars[v]));
    }
    _space = Space(spaceData);

    // Setup row constraints.

    SCIP_CONSHDLR* linearHandler = SCIPfindConshdlr(scip, "linear");
    std::size_t numLinearConstraints = SCIPconshdlrGetNConss(linearHandler);
    if (numLinearConstraints > 0)
    {
      SCIP_CONS** linearConstraints = SCIPconshdlrGetConss(linearHandler);
      _rowConstraints.reserve(2*numLinearConstraints);
      _rowNames.reserve(numLinearConstraints);
      for (std::size_t i = 0; i < numLinearConstraints; ++i)
      {
        SCIP_CONS* cons = linearConstraints[i];

        double realLhs = SCIPgetLhsLinear(scip, cons);
        Rational lhs;
        bool hasLhs = realLhs > -SCIPinfinity(scip);
        if (hasLhs)
          reconstruct(realLhs, lhs, SCIPfeastol(scip));
        double realRhs = SCIPgetRhsLinear(scip, cons);
        Rational rhs = 0;
        bool hasRhs = realRhs < SCIPinfinity(scip);
        if (hasRhs)
          reconstruct(realRhs, rhs, SCIPfeastol(scip));
        int nvars = SCIPgetNVarsLinear(scip, cons);
        SCIP_VAR** consVars = SCIPgetVarsLinear(scip, cons);
        double* vals = SCIPgetValsLinear(scip, cons);

        // Create normal vector.

        VectorData* data = new VectorData();
        for (int j = 0; j < nvars; ++j)
        {
          Rational x;
          reconstruct(vals[j], x, SCIPfeastol(scip));
          if (x != 0)
            data->add(varToIndexMap[consVars[j]], x);
        }
        Vector normal = Vector(data);

        // Create <= inequality or equation.

        if (hasRhs && !hasLhs)
        {
          _rowConstraints.push_back(LinearConstraint('<', normal, rhs));
          _rowNames.push_back(SCIPconsGetName(cons));
        }
        else if (!hasRhs && hasLhs)
        {
          _rowConstraints.push_back(LinearConstraint('>', normal, lhs));
          _rowNames.push_back(SCIPconsGetName(cons));
        }
        else if (!hasLhs && !hasRhs)
        {
          std::cerr << "Warning: Ignoring SCIP constraint " << SCIPconsGetName(cons) << " with two infinite sides." << std::endl;
        }
        else if (lhs == rhs)
        {
          _rowConstraints.push_back(LinearConstraint('=', normal, rhs));
          _rowNames.push_back(SCIPconsGetName(cons));
        }
        else
        {
          _rowConstraints.push_back(LinearConstraint('<', normal, rhs));
          _rowNames.push_back(std::string(SCIPconsGetName(cons)) + "_rhs");
          _rowConstraints.push_back(LinearConstraint('>', normal, lhs));
          _rowNames.push_back(std::string(SCIPconsGetName(cons)) + "_lhs");
        }
      }
    }
  }

#endif

  MixedIntegerSet::~MixedIntegerSet()
  {

  }

//   bool MixedIntegerSet::checkPointBounds(const SVectorRational* point) const
//   {
//     for (int p = point->size() - 1; p >= 0; --p)
//     {
//       std::size_t v = point->index(p);
//       const Rational& x = point->value(p);
//       if (_columns.lower(v) > -infinity && x < _columns.lower(v))
//         return false;
//       if (_columns.upper(v) < infinity && x > _columns.upper(v))
//         return false;
//     }
//     return true;
//   }
// 
//   bool MIP::checkPointRows(const SVectorRational* point)
//   {
//     _worker.clear();
//     _worker.assign(*point);
// 
//     for (int r = 0; r < _rows.num(); ++r)
//     {
//       const Rational activity = _rows.rowVector(r) * _worker;
//       if (_rows.lhs(r) > -infinity && activity < _rows.lhs(r))
//         return false;
//       if (_rows.rhs(r) < infinity && activity > _rows.rhs(r))
//         return false;
//     }
//     return true;
//   }
// 
//   bool MIP::checkPointIntegral(const SVectorRational* point) const
//   {
//     for (int p = point->size() - 1; p >= 0; --p)
//     {
//       if (!_integrality[point->index(p)])
//         continue;
//       if (!ipo::isIntegral(point->value(p)))
//         return false;
//     }
//     return true;
//   }
// 
//   bool MIP::checkPoint(const SVectorRational* point)
//   {
//     return checkPointIntegral(point) && checkPointBounds(point) && checkPointRows(point);
//   }
// 
//   bool MIP::checkRayBounds(const SVectorRational* ray) const
//   {
//     for (int p = ray->size() - 1; p >= 0; --p)
//     {
//       std::size_t v = ray->index(p);
//       const Rational& x = ray->value(p);
//       if (_columns.lower(v) > -infinity && x < 0)
//         return false;
//       if (_columns.upper(v) < infinity && x > 0)
//         return false;
//     }
//     return true;
//   }
// 
//   bool MIP::checkRayRows(const SVectorRational* ray)
//   {
//     _worker.clear();
//     _worker.assign(*ray);
// 
//     for (int r = 0; r < _rows.num(); ++r)
//     {
//       const Rational activity = _rows.rowVector(r) * _worker;
//       if (_rows.lhs(r) > -infinity && activity < 0)
//         return false;
//       if (_rows.rhs(r) < infinity && activity > 0)
//         return false;
//     }
//     return true;
//   }
// 
//   bool MIP::checkRay(const SVectorRational* ray)
//   {
//     return checkRayBounds(ray) && checkRayRows(ray);
//   }

  void MixedIntegerSet::setFace(const LinearConstraint& newFace)
  {
    _currentFace = newFace;
  }

//   void MixedIntegerSet::getConstraints(LPRowSetRational& rows, bool inequalities, bool equations,
//       std::vector<std::string>* names)
//   {
//     rows.reMax(_rows.max());
//     DSVectorRational vector;
//     for (int r = 0; r < _rows.num(); ++r)
//     {
//       const Rational& lower = _rows.lhs(r);
//       const Rational& upper = _rows.rhs(r);
//       if (lower == upper)
//       {
//         if (equations)
//         {
//           rows.add(lower, _rows.rowVector(r), upper);
//           if (names)
//             names->push_back(rowName(r));
//         }
//         continue;
//       }
//       if (!inequalities)
//         continue;
// 
//       if (upper < infinity)
//       {
//         rows.add(-infinity, _rows.rowVector(r), upper);
//         if (names)
//           names->push_back(rowName(r) + "-rhs");
//       }
//       if (lower > -infinity)
//       {
//         vector = _rows.rowVector(r);
//         vector *= -1;
//         rows.add(-infinity, vector, -lower);
//         if (names)
//           names->push_back(rowName(r) + "-lhs");
//       }
//     }
//   }
// 
//   void MIP::getFixedVariableEquations(LPRowSetRational& rows, std::vector<std::string>* names)
//   {
//     DSVectorRational vector;
//     for (std::size_t c = 0; c < _columns.num(); ++c)
//     {
//       if (_columns.upper(c) != _columns.lower(c))
//         continue;
// 
//       const Rational& value = _columns.upper(c);
//       vector.clear();
//       vector.add(c, Rational(1));
//       rows.add(value, vector, value);
//       if (names)
//         names->push_back("fixed-" + space()[c]);
//     }
//   }

  MIPOracleBase::MIPOracleBase(const std::string& name, const std::shared_ptr<OracleBase>& nextOracle)
    : OracleBase(name, nextOracle)
  {
    
    
  }

  void MIPOracleBase::initialize(const std::shared_ptr<MixedIntegerSet>& mixedIntegerSet)
  {
    _mixedIntegerSet = mixedIntegerSet;

    OracleBase::initializeSpace(_mixedIntegerSet->space());

    assert(!_nextOracle || _nextOracle->space() == space());

    _spx = new soplex::SoPlex;
    _spx->setIntParam(soplex::SoPlex::SOLVEMODE, soplex::SoPlex::SOLVEMODE_RATIONAL);
    _spx->setIntParam(soplex::SoPlex::SYNCMODE, soplex::SoPlex::SYNCMODE_AUTO);
    _spx->setRealParam(soplex::SoPlex::FEASTOL, 0.0);
    _spx->setBoolParam(soplex::SoPlex::RATREC, true);
    _spx->setBoolParam(soplex::SoPlex::RATFAC, true);
    _spx->setIntParam(soplex::SoPlex::OBJSENSE, soplex::SoPlex::OBJSENSE_MAXIMIZE);
    _spx->setIntParam(soplex::SoPlex::VERBOSITY, soplex::SoPlex::VERBOSITY_ERROR);

    soplex::LPColSetRational cols(mixedIntegerSet->numVariables());
    soplex::DSVectorRational zero;
    for (std::size_t v = 0; v < mixedIntegerSet->numVariables(); ++v)
    {
      const MixedIntegerSet::Variable& variable = mixedIntegerSet->variable(v);
      cols.add(Rational(0), variable.lowerBound, zero, variable.upperBound);
    }
    _spx->addColsRational(cols);

    addToLP(*_spx, mixedIntegerSet->rowConstraints());

    _numRows = _spx->numRowsRational();

    _objective = new double[space().dimension()];
    _lpResult.reDim(space().dimension());
  }

  MIPOracleBase::~MIPOracleBase()
  {
    delete[] _objective;
    delete _spx;
  }

  HeuristicLevel MIPOracleBase::maximizeImplementation(OracleResult& result, const soplex::VectorRational& objective,
    const ObjectiveBound& objectiveBound, HeuristicLevel minHeuristic, HeuristicLevel maxHeuristic, bool& sort, bool& checkDups)
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
      _objective[v] = double(factor * objective[v]);

    // Call the solver's method.
    
    assert(_points.empty());
    assert(_rays.empty());
    
    solverMaximize(_objective, soplex::infinity, _points, _rays);

    if (!_points.empty())
    {
      sort = true;
      checkDups = true;
      result.points.reserve(_points.size());
      prepareSolver(objective);

      for (std::size_t i = 0; i < _points.size(); ++i)
      {
        Rational objValue;
        Vector vector = extendPoint(_points[i], objValue);
        result.points.push_back(OracleResult::Point(vector, objValue));
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
      Vector vector = computeRay();
      restoreSolver();

      result.rays.push_back(OracleResult::Ray(vector));
      return 0;
    }
    return heuristicLevel();
  }

  void MIPOracleBase::separatePoint(const soplex::VectorRational& point, soplex::LPRowSetRational& cuts)
  {
    
  }
  
  void MIPOracleBase::separateRay(const soplex::VectorRational& ray, soplex::LPRowSetRational& cuts)
  {
    
  }

  void MIPOracleBase::setFace(const LinearConstraint& newFace)
  {
    OracleBase::setFace(newFace);
    
    throw std::runtime_error("MIPOracleBase::setFace not implemented, yet!");
  }

  void MIPOracleBase::prepareSolver(const soplex::VectorRational& objective)
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

  Vector MIPOracleBase::extendPoint(double* approxPoint, Rational& objectiveValue)
  {
    std::size_t n = space().dimension();
    
    // Fix variable bounds for the integer variables.

    for (std::size_t v = 0; v < n; ++v)
    {
      if (_mixedIntegerSet->isIntegral(v))
      {
        Rational value = Rational(int(approxPoint[v] + 0.5));
        _spx->changeBoundsRational(v, value, value);
      }
      else
      {
        assert(_mixedIntegerSet->variable(v).lowerBound == _spx->lowerRational(v));
        assert(_mixedIntegerSet->variable(v).upperBound == _spx->upperRational(v));
      }
    }

//     for (int r = 0; r < _spx->numRowsRational(); ++r)
//     {
//       const SVectorRational& rowVector = _spx->rowVectorRational(r);
//       std::cout << "// Row " << r << std::endl;
//       std::cout << "row.clear();\n";
//       for (int p = 0; p < rowVector.size(); ++p)
//         std::cout << "row.add(" << rowVector.index(p) << ", Rational(" << rowVector.value(p) << "));" << std::endl;
//       if (_spx->lhsRational(r) <= -infinity)
//       {
//         if (_spx->rhsRational(r) >= infinity)
//         {
//           std::cout << "spx.addRowRational(LPRowRational(-infinity, row, infinity));\n";
//         }
//         else
//         {
//           std::cout << "spx.addRowRational(LPRowRational(-infinity, row, Rational(" << _spx->rhsRational(r) << ")));\n";
//         }
//       }
//       else
//       {
//         if (_spx->rhsRational(r) >= infinity)
//         {
//           std::cout << "spx.addRowRational(LPRowRational(Rational(" << _spx->lhsRational(r) << "), row, infinity));\n";
//         }
//         else
//         {
//           std::cout << "spx.addRowRational(LPRowRational(Rational(" << _spx->lhsRational(r) << "), row, Rational("
//         << _spx->rhsRational(r) << ")));\n";
//         }
//       }
//     }
//     
//     _spx->clearBasis();

    while (true)
    {
//      _spx->clearBasis();
//       _spx->writeBasisFile("extend.bas");
//      _spx->writeFileRational("extend.lp");
      
//      std::cerr << "Calling solve." << std::endl;
      soplex::SPxSolver::Status status = _spx->solve();
//      std::cerr << "Called solve: " << status << std::endl;
      
      if (status == soplex::SPxSolver::UNBOUNDED)
      {
        _spx->getPrimalRayRational(_lpResult);
        _separateResult.clear();
        separateRay(_lpResult, _separateResult);
        if (_separateResult.num() == 0)
          throw std::runtime_error("MIPOracle: Claim is bounded, but candidate ray is not separated.");

        _spx->addRowsRational(_separateResult);  
      }
      else if (status != soplex::SPxSolver::OPTIMAL)
        throw std::runtime_error("MIPOracle: Claim is bounded, point could not be extended.");

      _spx->getPrimalRational(_lpResult);
      _separateResult.clear();
      separatePoint(_lpResult, _separateResult);
      if (_separateResult.num() == 0)
        break;

      _spx->addRowsRational(_separateResult);
    }

    // Create point as sparse vector.

    std::size_t size = 0;
    for (std::size_t v = 0; v < n; ++v)
    {
      if (_lpResult[v] != 0)
        ++size;
    }
    VectorData* data = new VectorData(size);
    for (std::size_t v = 0; v < n; ++v)
    {
      if (_lpResult[v] != 0)
        data->add(v, _lpResult[v]);
    }
    assert(data->isSorted());
    objectiveValue = _spx->objValueRational();
    return Vector(data);
  }

  Vector MIPOracleBase::computeRay()
  {
    std::size_t n = space().dimension();

    // Set original variable bounds for the integer variables.

    for (std::size_t v = 0; v < n; ++v)
    {
      if (_mixedIntegerSet->isIntegral(v))
      {
        _spx->changeBoundsRational(v, _mixedIntegerSet->variable(v).lowerBound, _mixedIntegerSet->variable(v).upperBound);
      }
      else
      {
        assert(_mixedIntegerSet->variable(v).lowerBound == _spx->lowerRational(v));
        assert(_mixedIntegerSet->variable(v).upperBound == _spx->upperRational(v));
      }
    }

    while (true)
    {
      soplex::SPxSolver::Status status = _spx->solve();
      if (status != soplex::SPxSolver::UNBOUNDED)
        throw std::runtime_error("MIPOracle: Claim is unbounded, no ray found.");

      _spx->getPrimalRayRational(_lpResult);
      _separateResult.clear();
      separateRay(_lpResult, _separateResult);
      if (_separateResult.num() == 0)
        break;

      _spx->addRowsRational(_separateResult);
    }

    // Create ray as sparse vector.

    std::size_t size = 0;
    for (std::size_t v = 0; v < n; ++v)
    {
      if (_lpResult[v] != 0)
        ++size;
    }
    VectorData* data = new VectorData(size);
    for (std::size_t v = 0; v < n; ++v)
    {
      if (_lpResult[v] != 0)
        data->add(v, _lpResult[v]);
    }
    assert(data->isSorted());
    return Vector(data);
  }

} /* namespace ipo */
