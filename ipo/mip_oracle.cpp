#include "mip_oracle.h"

#include "rational.h"
#include "reconstruct.h"

#ifdef IPO_WITH_SCIP
  #ifdef NDEBUG
    #undef NDEBUG
    #include <scip/cons_linear.h>
    #include "scip_oracle.h"
    #define NDEBUG
  #else
    #include <scip/cons_linear.h>
    #include "scip_oracle.h"
  #endif
#endif

namespace ipo {

//   MixedIntegerLin::Variable::Variable()
//     : lowerBound(-soplex::infinity), upperBound(soplex::infinity), integral(false)
//   {
// 
//   }
// 
//   MixedIntegerSet::Variable::Variable(const Rational& lower, const Rational& upper, bool integ)
//     : lowerBound(lower), upperBound(upper), integral(integ)
//   {
// 
//   }
// 
// 
// #ifdef IPO_WITH_SCIP
// 
//   MixedIntegerSet::MixedIntegerSet(SCIP* scip) : _currentFace()
//   {
//     std::size_t n = SCIPgetNOrigVars(scip);
//     SCIP_VAR** origVars = SCIPgetOrigVars(scip);
// //     _worker.reDim(n, true);
// 
//     // Create map
// 
//     SCIPvarToIndexMap varToIndexMap;
//     getSCIPvarToIndexMap(scip, varToIndexMap);
// 
//     // Setup variables.
// 
//     SpaceData* spaceData = new SpaceData();
//     _variables.resize(n);
//     for (std::size_t v = 0; v < n; ++v)
//     {
//       double realLower = SCIPvarGetLbGlobal(origVars[v]);
//       if (realLower > -SCIPinfinity(scip))
//         reconstruct(realLower, _variables[v].lowerBound, SCIPfeastol(scip));
//       else
//         _variables[v].lowerBound = -soplex::infinity;
//       double realUpper = SCIPvarGetUbGlobal(origVars[v]);
//       if (realUpper < SCIPinfinity(scip))
//         reconstruct(realUpper, _variables[v].upperBound, SCIPfeastol(scip));
//       else
//         _variables[v].upperBound = soplex::infinity;
//       _variables[v].integral = SCIPvarGetType(origVars[v]) != SCIP_VARTYPE_CONTINUOUS;
// 
//       spaceData->addVariable(SCIPvarGetName(origVars[v]));
//     }
//     _space = Space(spaceData);
// 
//     // Setup row constraints.
// 
//     SCIP_CONSHDLR* linearHandler = SCIPfindConshdlr(scip, "linear");
//     std::size_t numLinearConstraints = SCIPconshdlrGetNConss(linearHandler);
//     if (numLinearConstraints > 0)
//     {
//       SCIP_CONS** linearConstraints = SCIPconshdlrGetConss(linearHandler);
//       _rowConstraints.reserve(2*numLinearConstraints);
//       _rowNames.reserve(numLinearConstraints);
//       for (std::size_t i = 0; i < numLinearConstraints; ++i)
//       {
//         SCIP_CONS* cons = linearConstraints[i];
// 
//         double realLhs = SCIPgetLhsLinear(scip, cons);
//         Rational lhs;
//         bool hasLhs = realLhs > -SCIPinfinity(scip);
//         if (hasLhs)
//           reconstruct(realLhs, lhs, SCIPfeastol(scip));
//         double realRhs = SCIPgetRhsLinear(scip, cons);
//         Rational rhs = 0;
//         bool hasRhs = realRhs < SCIPinfinity(scip);
//         if (hasRhs)
//           reconstruct(realRhs, rhs, SCIPfeastol(scip));
//         int nvars = SCIPgetNVarsLinear(scip, cons);
//         SCIP_VAR** consVars = SCIPgetVarsLinear(scip, cons);
//         double* vals = SCIPgetValsLinear(scip, cons);
// 
//         // Create normal vector.
// 
//         VectorData* data = new VectorData();
//         for (int j = 0; j < nvars; ++j)
//         {
//           Rational x;
//           reconstruct(vals[j], x, SCIPfeastol(scip));
//           if (x != 0)
//             data->add(varToIndexMap[consVars[j]], x);
//         }
//         Vector normal = Vector(data);
// 
//         // Create <= inequality or equation.
// 
//         if (hasRhs && !hasLhs)
//         {
//           _rowConstraints.push_back(LinearConstraint('<', normal, rhs));
//           _rowNames.push_back(SCIPconsGetName(cons));
//         }
//         else if (!hasRhs && hasLhs)
//         {
//           _rowConstraints.push_back(LinearConstraint('>', normal, lhs));
//           _rowNames.push_back(SCIPconsGetName(cons));
//         }
//         else if (!hasLhs && !hasRhs)
//         {
//           std::cerr << "Warning: Ignoring SCIP constraint " << SCIPconsGetName(cons) << " with two infinite sides." << std::endl;
//         }
//         else if (lhs == rhs)
//         {
//           _rowConstraints.push_back(LinearConstraint('=', normal, rhs));
//           _rowNames.push_back(SCIPconsGetName(cons));
//         }
//         else
//         {
//           _rowConstraints.push_back(LinearConstraint('<', normal, rhs));
//           _rowNames.push_back(std::string(SCIPconsGetName(cons)) + "_rhs");
//           _rowConstraints.push_back(LinearConstraint('>', normal, lhs));
//           _rowNames.push_back(std::string(SCIPconsGetName(cons)) + "_lhs");
//         }
//       }
//     }
//   }
// 
// #endif /* IPO_WITH_SCIP */
// 
//   MixedIntegerSet::~MixedIntegerSet()
//   {
// 
//   }
// 
//   LinearConstraint MixedIntegerSet::upperBoundConstraint(std::size_t variable) const
//   {
//     VectorData* data = new VectorData(1);
//     data->add(variable, Rational(1));
//     return LinearConstraint(_variables[variable].upperBound == _variables[variable].lowerBound ? '=' : '<', Vector(data),
//       _variables[variable].upperBound);
//   }
// 
//   LinearConstraint MixedIntegerSet::lowerBoundConstraint(std::size_t variable) const
//   {
//     VectorData* data = new VectorData(1);
//     data->add(variable, Rational(1));
//     return LinearConstraint(_variables[variable].upperBound == _variables[variable].lowerBound ? '=' : '>', Vector(data),
//       _variables[variable].lowerBound);
//   }
// 
// 
// //   bool MixedIntegerSet::checkPointBounds(const SVectorRational* point) const
// //   {
// //     for (int p = point->size() - 1; p >= 0; --p)
// //     {
// //       std::size_t v = point->index(p);
// //       const Rational& x = point->value(p);
// //       if (_columns.lower(v) > -infinity && x < _columns.lower(v))
// //         return false;
// //       if (_columns.upper(v) < infinity && x > _columns.upper(v))
// //         return false;
// //     }
// //     return true;
// //   }
// //
// //   bool MIP::checkPointRows(const SVectorRational* point)
// //   {
// //     _worker.clear();
// //     _worker.assign(*point);
// //
// //     for (int r = 0; r < _rows.num(); ++r)
// //     {
// //       const Rational activity = _rows.rowVector(r) * _worker;
// //       if (_rows.lhs(r) > -infinity && activity < _rows.lhs(r))
// //         return false;
// //       if (_rows.rhs(r) < infinity && activity > _rows.rhs(r))
// //         return false;
// //     }
// //     return true;
// //   }
// //
// //   bool MIP::checkPointIntegral(const SVectorRational* point) const
// //   {
// //     for (int p = point->size() - 1; p >= 0; --p)
// //     {
// //       if (!_integrality[point->index(p)])
// //         continue;
// //       if (!ipo::isIntegral(point->value(p)))
// //         return false;
// //     }
// //     return true;
// //   }
// //
// //   bool MIP::checkPoint(const SVectorRational* point)
// //   {
// //     return checkPointIntegral(point) && checkPointBounds(point) && checkPointRows(point);
// //   }
// //
// //   bool MIP::checkRayBounds(const SVectorRational* ray) const
// //   {
// //     for (int p = ray->size() - 1; p >= 0; --p)
// //     {
// //       std::size_t v = ray->index(p);
// //       const Rational& x = ray->value(p);
// //       if (_columns.lower(v) > -infinity && x < 0)
// //         return false;
// //       if (_columns.upper(v) < infinity && x > 0)
// //         return false;
// //     }
// //     return true;
// //   }
// //
// //   bool MIP::checkRayRows(const SVectorRational* ray)
// //   {
// //     _worker.clear();
// //     _worker.assign(*ray);
// //
// //     for (int r = 0; r < _rows.num(); ++r)
// //     {
// //       const Rational activity = _rows.rowVector(r) * _worker;
// //       if (_rows.lhs(r) > -infinity && activity < 0)
// //         return false;
// //       if (_rows.rhs(r) < infinity && activity > 0)
// //         return false;
// //     }
// //     return true;
// //   }
// //
// //   bool MIP::checkRay(const SVectorRational* ray)
// //   {
// //     return checkRayBounds(ray) && checkRayRows(ray);
// //   }
// 
//   void MixedIntegerSet::setFace(const LinearConstraint& newFace)
//   {
//     _currentFace = newFace;
//   }
// 
// //   void MixedIntegerSet::getConstraints(LPRowSetRational& rows, bool inequalities, bool equations,
// //       std::vector<std::string>* names)
// //   {
// //     rows.reMax(_rows.max());
// //     DSVectorRational vector;
// //     for (int r = 0; r < _rows.num(); ++r)
// //     {
// //       const Rational& lower = _rows.lhs(r);
// //       const Rational& upper = _rows.rhs(r);
// //       if (lower == upper)
// //       {
// //         if (equations)
// //         {
// //           rows.add(lower, _rows.rowVector(r), upper);
// //           if (names)
// //             names->push_back(rowName(r));
// //         }
// //         continue;
// //       }
// //       if (!inequalities)
// //         continue;
// //
// //       if (upper < infinity)
// //       {
// //         rows.add(-infinity, _rows.rowVector(r), upper);
// //         if (names)
// //           names->push_back(rowName(r) + "-rhs");
// //       }
// //       if (lower > -infinity)
// //       {
// //         vector = _rows.rowVector(r);
// //         vector *= -1;
// //         rows.add(-infinity, vector, -lower);
// //         if (names)
// //           names->push_back(rowName(r) + "-lhs");
// //       }
// //     }
// //   }
// //
// //   void MIP::getFixedVariableEquations(LPRowSetRational& rows, std::vector<std::string>* names)
// //   {
// //     DSVectorRational vector;
// //     for (std::size_t c = 0; c < _columns.num(); ++c)
// //     {
// //       if (_columns.upper(c) != _columns.lower(c))
// //         continue;
// //
// //       const Rational& value = _columns.upper(c);
// //       vector.clear();
// //       vector.add(c, Rational(1));
// //       rows.add(value, vector, value);
// //       if (names)
// //         names->push_back("fixed-" + space()[c]);
// //     }
// //   }

  MIPOracleBase::MIPOracleBase(const std::string& name, const std::shared_ptr<OracleBase>& nextOracle)
    : OracleBase(name, nextOracle)
  {


  }

  void MIPOracleBase::initialize(const std::shared_ptr<MixedIntegerLinearSet>& mixedIntegerLinearSet)
  {
    _mixedIntegerLinearSet = mixedIntegerLinearSet;
    _correctionLP = std::make_shared<LinearProgram>(mixedIntegerLinearSet);

    OracleBase::initializeSpace(_mixedIntegerLinearSet->space());

    assert(!_nextOracle || _nextOracle->space() == space());

//     _spx = new soplex::SoPlex;
//     _spx->setIntParam(soplex::SoPlex::SOLVEMODE, soplex::SoPlex::SOLVEMODE_RATIONAL);
//     _spx->setIntParam(soplex::SoPlex::SYNCMODE, soplex::SoPlex::SYNCMODE_AUTO);
//     _spx->setRealParam(soplex::SoPlex::FEASTOL, 0.0);
//     _spx->setBoolParam(soplex::SoPlex::RATREC, true);
//     _spx->setBoolParam(soplex::SoPlex::RATFAC, true);
//     _spx->setIntParam(soplex::SoPlex::OBJSENSE, soplex::SoPlex::OBJSENSE_MAXIMIZE);
// //     _spx->setIntParam(soplex::SoPlex::VERBOSITY, soplex::SoPlex::VERBOSITY_ERROR);
// 
//     soplex::LPColSetRational cols(mixedIntegerLinearSet->numVariables());
//     soplex::DSVectorRational zero;
//     for (std::size_t v = 0; v < mixedIntegerSet->numVariables(); ++v)
//     {
//       const MixedIntegerSet::Variable& variable = mixedIntegerSet->variable(v);
//       cols.add(Rational(0), variable.lowerBound, zero, variable.upperBound);
//     }
//     _spx->addColsRational(cols);
// 
//     addToLP(*_spx, mixedIntegerSet->rowConstraints());
// 
//     _numRows = _spx->numRowsRational();

    _objective = new double[space().dimension()];
  }

  MIPOracleBase::~MIPOracleBase()
  {
    delete[] _objective;
//     delete _spx;
  }

  HeuristicLevel MIPOracleBase::maximizeImplementation(OracleResult& result, const soplex::VectorRational& objective,
    const ObjectiveBound& objectiveBound, HeuristicLevel minHeuristic, HeuristicLevel maxHeuristic, bool& sort, bool& checkDups)
  {
    bool hitLimit = false;
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

    solverMaximize(_objective, soplex::infinity, _points, _rays, hitLimit);

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

    if (heuristicLevel() == 0 && hitLimit)
    {
      std::stringstream str;
      str << "Oracle \"" << name() << "\" reached a limit.";
      throw std::runtime_error(str.str());
    }

    return  heuristicLevel();
  }

  void MIPOracleBase::separatePoint(const Vector& point, std::vector< LinearConstraint >& cuts)
  {

  }

  void MIPOracleBase::separateRay(const Vector& point, std::vector< LinearConstraint >& cuts)
  {

  }

  void MIPOracleBase::setFace(const LinearConstraint& newFace)
  {
    OracleBase::setFace(newFace);

    throw std::runtime_error("MIPOracleBase::setFace not implemented, yet!");
  }

  void MIPOracleBase::prepareSolver(const soplex::VectorRational& objective)
  {
    std::vector<Rational> obj(objective.dim());
    for (std::size_t v = 0; v < obj.size(); ++v)
      obj[v] = objective[v];
    _correctionLP->changeObjective(obj);
  }

  void MIPOracleBase::restoreSolver()
  {
    // Remove all added rows.
    
    _correctionLP->removeLastConstraints(_mixedIntegerLinearSet->numRows());
  }

  Vector MIPOracleBase::extendPoint(double* approxPoint, Rational& objectiveValue)
  {
    std::size_t n = space().dimension();

    // Fix variable bounds for the integer variables.

    for (std::size_t v = 0; v < n; ++v)
    {
      if (_mixedIntegerLinearSet->isIntegral(v))
      {
        Rational value = Rational(int(approxPoint[v] + 0.5));
        _correctionLP->changeBounds(v, value, value);
      }
      else
      {
        assert(_mixedIntegerLinearSet->lowerBound(v) == _correctionLP->lowerBound(v));
        assert(_mixedIntegerLinearSet->upperBound(v) == _correctionLP->upperBound(v));
      }
    }

    Vector vector;
    while (true)
    {
      LinearProgram::Result result = _correctionLP->solve(vector, objectiveValue);
      
      if (result == LinearProgram::UNBOUNDED)
      {
        separateRay(vector, _cuts);
        if (_cuts.size() == 0)
          throw std::runtime_error("MIPOracle: Claim is bounded, but candidate ray is not separated.");

        for (std::size_t i = 0; i < _cuts.size(); ++i)
          _correctionLP->addConstraint(_cuts[i]);
        _cuts.clear();
      }
      else if (result != LinearProgram::OPTIMAL)
        throw std::runtime_error("MIPOracle: Claim is bounded, point could not be extended.");

      separatePoint(vector, _cuts);
      if (_cuts.size() == 0)
        break;

      for (std::size_t i = 0; i < _cuts.size(); ++i)
        _correctionLP->addConstraint(_cuts[i]);
      _cuts.clear();
    }

    return vector;
  }

  Vector MIPOracleBase::computeRay()
  {
    std::size_t n = space().dimension();

    // Set original variable bounds for the integer variables.

    for (std::size_t v = 0; v < n; ++v)
    {
      if (_mixedIntegerLinearSet->isIntegral(v))
      {
        _correctionLP->changeBounds(v, _mixedIntegerLinearSet->lowerBound(v), _mixedIntegerLinearSet->upperBound(v));
      }
      else
      {
        assert(_mixedIntegerLinearSet->lowerBound(v) == _correctionLP->lowerBound(v));
        assert(_mixedIntegerLinearSet->upperBound(v) == _correctionLP->upperBound(v));
      }
    }

    Vector vector;
    while (true)
    {
      Rational objectiveValue;
      LinearProgram::Result result = _correctionLP->solve(vector, objectiveValue);

      if (result != LinearProgram::UNBOUNDED)
        throw std::runtime_error("MIPOracle: Claim is unbounded, no ray found.");

      separateRay(vector, _cuts);
      if (_cuts.empty())
        break;

      for (std::size_t i = 0; i < _cuts.size(); ++i)
        _correctionLP->addConstraint(_cuts[i]);
      _cuts.clear();
    }

    return vector;
  }

} /* namespace ipo */
