#include <ipo/oracles_scip.hpp>

#include <cassert>
#include <functional>
#include <chrono>

#ifdef NDEBUG
  #undef NDEBUG
  #include <scip/scipdefplugins.h>
  #include <scip/cons_linear.h>
  #include <scip/pub_misc.h>
  #include <scip/pub_var.h>
  #define NDEBUG
#else
  #include <scip/scipdefplugins.h>
  #include <scip/cons_linear.h>
  #include <scip/pub_misc.h>
  #include <scip/pub_var.h>
#endif


/**
 * \brief Macro to raise a SCIPException in case of a SCIP error.
 **/

#define SCIP_CALL_EXC(x) \
{ \
  SCIP_RETCODE _retcode; \
  if ((_retcode = (x)) != SCIP_OKAY) \
    throw ipo::SCIPException(_retcode); \
}

static const int SCIP_SEPARATION_CHECK_TIMELIMIT_FREQUENCY = 1000;

namespace ipo
{
  /**
   * \brief Exception handling class for SCIP.
   *
   * Represents a SCIP error in C++.
   **/

  class SCIPException: public std::exception
  {
  public:

    /**
     * \brief Constructs a SCIPException from an error code.
     **/

    SCIPException(SCIP_RETCODE retcode)
      : _retcode(retcode)
    {
      switch (retcode)
      {
      case SCIP_OKAY:
        SCIPsnprintf(_message, 256, "normal termination");
      break;
      case SCIP_ERROR:
        SCIPsnprintf(_message, 256, "unspecified error");
      break;
      case SCIP_NOMEMORY:
        SCIPsnprintf(_message, 256, "insufficient memory error");
      break;
      case SCIP_READERROR:
        SCIPsnprintf(_message, 256, "read error");
      break;
      case SCIP_WRITEERROR:
        SCIPsnprintf(_message, 256, "write error");
      break;
      case SCIP_NOFILE:
        SCIPsnprintf(_message, 256, "file not found error");
      break;
      case SCIP_FILECREATEERROR:
        SCIPsnprintf(_message, 256, "cannot create file");
      break;
      case SCIP_LPERROR:
        SCIPsnprintf(_message, 256, "error in LP solver");
      break;
      case SCIP_NOPROBLEM:
        SCIPsnprintf(_message, 256, "no problem exists");
      break;
      case SCIP_INVALIDCALL:
        SCIPsnprintf(_message, 256, "method cannot be called at this time in solution process");
      break;
      case SCIP_INVALIDDATA:
        SCIPsnprintf(_message, 256, "method cannot be called with this type of data");
      break;
      case SCIP_INVALIDRESULT:
        SCIPsnprintf(_message, 256, "method returned an invalid result code");
      break;
      case SCIP_PLUGINNOTFOUND:
        SCIPsnprintf(_message, 256, "a required plugin was not found");
      break;
      case SCIP_PARAMETERUNKNOWN:
        SCIPsnprintf(_message, 256, "the parameter with the given name was not found");
      break;
      case SCIP_PARAMETERWRONGTYPE:
        SCIPsnprintf(_message, 256, "the parameter is not of the expected type");
      break;
      case SCIP_PARAMETERWRONGVAL:
        SCIPsnprintf(_message, 256, "the value is invalid for the given parameter");
      break;
      case SCIP_KEYALREADYEXISTING:
        SCIPsnprintf(_message, 256, "the given key is already existing in table");
      break;
      case SCIP_MAXDEPTHLEVEL:
        SCIPsnprintf(_message, 256, "maximal branching depth level exceeded");
      break;
      case SCIP_BRANCHERROR:
        SCIPsnprintf(_message, 256, "branching could not be performed (e.g. too large values in variable domain)");
      break;
      default:
        SCIPsnprintf(_message, 256, "unknown error code %d", retcode);
      break;
      }
    }

    /**
     * \brief Destructor.
     **/

    ~SCIPException(void) throw ()
    {

    }

    /**
     * \brief Returns the error message
     **/

    const char* what(void) const throw ()
    {
      return _message;
    }

  private:
    /// Buffer for the error message.
    char _message[256];
    /// SCIP error code.
    SCIP_RETCODE _retcode;
  };

  /**
   * \brief Auxiliary method that iterates over the bounds as constraints.
   *
   * Auxiliary method that iterates over the bound constraints of all variables calling \p visitor
   * with each of them as a constraint.
   *
   * \p ranged Whether to returned ranged rows.
   */

  void SCIPiterateBounds(SCIP* scip,
    const std::unordered_map<SCIP_VAR*, std::size_t>& variablesToCoordinates,
    std::function<void(Constraint<double>&& constraint)> visitor, bool ranged)
  {
    SCIP_VAR** vars = SCIPgetOrigVars(scip);
    int n = SCIPgetNOrigVars(scip);
    sparse_vector<double> entries;
    for (int v = 0; v < n; ++v)
    {
      double lhs = SCIPvarGetLbOriginal(vars[v]);
      if (SCIPisInfinity(scip, -lhs))
        lhs = -std::numeric_limits<double>::infinity();
      double rhs = SCIPvarGetUbOriginal(vars[v]);
      if (SCIPisInfinity(scip, rhs))
        rhs = std::numeric_limits<double>::infinity();
      if (ranged || SCIPisEQ(scip, lhs, rhs))
      {
        auto vector = std::make_shared<sparse_vector<double>>();
        vector->push_back(variablesToCoordinates.at(vars[v]), 1.0);
        visitor(Constraint<double>(lhs, vector, rhs));
      }
      else
      {
        if (!ipo::isPlusInfinity(rhs))
        {
          auto vector = std::make_shared<sparse_vector<double>>();
          vector->push_back(variablesToCoordinates.at(vars[v]), 1.0);
          visitor(Constraint<double>(-std::numeric_limits<double>::infinity(), vector, rhs));
        }
        if (!ipo::isMinusInfinity(lhs))
        {
          auto vector = std::make_shared<sparse_vector<double>>();
          vector->push_back(variablesToCoordinates.at(vars[v]), -1.0);
          visitor(Constraint<double>(-std::numeric_limits<double>::infinity(), vector, -lhs));
        }
      }
    }
  }

  /**
   * \brief Auxiliary method that iterates over the rows.
   *
   * Auxiliary method that iterates over the rows of linear and setppc constraints, calling
   * \p visitor with each of them. It passes the lhs, number of nonzeros, variables and
   * coefficients of the nonzeros, and the rhs of each row in that order.
   *
   * \p ranged Whether to returned ranged rows.
   */

  void SCIPiterateRows(SCIP* scip,
    const std::unordered_map<SCIP_VAR*, std::size_t>& variablesToCoordinates,
    std::function<void(Constraint<double>&& constraint)> visitor, bool ranged)
  {
    SCIP_CONS** conss = SCIPgetConss(scip);
    std::size_t m = SCIPgetNConss(scip);
    std::vector<double> ones;
    ones.resize(SCIPgetNOrigVars(scip), 1.0);
    std::vector<std::pair<std::size_t, double>> unsortedEntries;
    sparse_vector<double> entries;
    for (std::size_t c = 0; c < m; ++c)
    {
      SCIP_CONS* cons = conss[c];
      SCIP_CONSHDLR* conshdlr = SCIPconsGetHdlr(cons);
      std::size_t k = 0;
      SCIP_VAR** vars = nullptr;
      double* vals = nullptr;
      double lhs = -std::numeric_limits<double>::infinity();
      double rhs = std::numeric_limits<double>::infinity();
      const std::string name = SCIPconshdlrGetName(conshdlr);
      if (name == "linear")
      {
        k = SCIPgetNVarsLinear(scip, cons);
        vars = SCIPgetVarsLinear(scip, cons);
        vals = SCIPgetValsLinear(scip, cons);
        lhs = SCIPgetLhsLinear(scip, cons);
        rhs = SCIPgetRhsLinear(scip, cons);
      }
      else if (name == "setppc")
      {
        k = SCIPgetNVarsSetppc(scip, cons);
        vars = SCIPgetVarsSetppc(scip, cons);
        vals = &ones[0];
        lhs = SCIPgetTypeSetppc(scip, cons) == SCIP_SETPPCTYPE_PACKING ? 0.0 : 1.0;
        rhs = SCIPgetTypeSetppc(scip, cons) == SCIP_SETPPCTYPE_COVERING ? SCIPinfinity(scip) : 1.0;
      }
      else
      {
        throw std::runtime_error("SCIP oracles do not implement constraint type <" + name + ">.");
      }

      if (SCIPisInfinity(scip, rhs))
        rhs = std::numeric_limits<double>::infinity();
      if (SCIPisInfinity(scip, -lhs))
        lhs = -std::numeric_limits<double>::infinity();

      if (ranged || SCIPisEQ(scip, lhs, rhs))
      {
        unsortedEntries.clear();
        for (std::size_t i = 0; i < k; ++i)
          unsortedEntries.push_back(std::make_pair(variablesToCoordinates.at(vars[i]), vals[i]));
        auto vector = std::make_shared<sparse_vector<double>>(std::move(unsortedEntries));
        visitor(Constraint<double>(lhs, vector, rhs));
      }
      else
      {
        if (!ipo::isPlusInfinity(rhs))
        {
          unsortedEntries.clear();
          for (std::size_t i = 0; i < k; ++i)
            unsortedEntries.push_back(std::make_pair(variablesToCoordinates.at(vars[i]), vals[i]));
          auto vector = std::make_shared<sparse_vector<double>>(std::move(unsortedEntries));
          visitor(Constraint<double>(-std::numeric_limits<double>::infinity(), vector, rhs));
        }
        if (!ipo::isMinusInfinity(lhs))
        {
          unsortedEntries.clear();
          for (std::size_t i = 0; i < k; ++i)
            unsortedEntries.push_back(std::make_pair(variablesToCoordinates.at(vars[i]), -vals[i]));
          auto vector = std::make_shared<sparse_vector<double>>(std::move(unsortedEntries));
          visitor(Constraint<double>(-std::numeric_limits<double>::infinity(), vector, -lhs));
        }
      }
    }
  }

  SCIPSolver::SCIPSolver(SCIP* scip)
    : _scip(scip), _currentFace(std::make_shared<Constraint<double>>(alwaysSatisfiedConstraint<double>()))
  {
    SCIP_CALL_EXC( SCIPsetIntParam(_scip, "display/verblevel", 0) );
    initialize();
  }

  SCIPSolver::SCIPSolver(const std::string& fileName)
    : _currentFace(std::make_shared<Constraint<double>>(alwaysSatisfiedConstraint<double>()))
  {
    SCIP_CALL_EXC( SCIPcreate(&_scip) );
    SCIP_CALL_EXC( SCIPincludeDefaultPlugins(_scip) );
    SCIP_CALL_EXC( SCIPsetIntParam(_scip, "display/verblevel", 0) );
    SCIP_CALL_EXC( SCIPreadProb(_scip, fileName.c_str(), NULL) );

    initialize();
  }

  void SCIPSolver::initialize()
  {
    std::size_t n = SCIPgetNOrigVars(_scip);
    _variables.resize(n);
    std::vector<std::string> variableNames;
    variableNames.resize(n);
    _instanceObjective = new double[n];
    SCIP_VAR** vars = SCIPgetOrigVars(_scip);
    double scale = SCIPgetObjsense(_scip) == SCIP_OBJSENSE_MAXIMIZE ? 1.0 : -1.0;
    for (std::size_t i = 0; i < n; ++i)
    {
      _variables[i] = vars[i];
      variableNames[i] = SCIPvarGetName(vars[i]);
      _variablesToCoordinates[vars[i]] = i;
      _instanceObjective[i] = scale * SCIPvarGetObj(vars[i]);
    }
    SCIP_CALL_EXC( SCIPsetObjsense(_scip, SCIP_OBJSENSE_MAXIMIZE) );

    _name = SCIPgetProbName(_scip);
    _space = std::make_shared<Space>(variableNames);

    std::vector<bool> integrality(n);
    std::vector<std::pair<double, double>> bounds(n);
    for (std::size_t i = 0; i < n; ++i)
    {
      integrality[i] = SCIPvarGetType(vars[i]) != SCIP_VARTYPE_CONTINUOUS;
      bounds[i].first = SCIPvarGetLbOriginal(vars[i]);
      if (SCIPisInfinity(_scip, -bounds[i].first))
        bounds[i].first = -std::numeric_limits<double>::infinity();
      bounds[i].second = SCIPvarGetUbOriginal(vars[i]);
      if (SCIPisInfinity(_scip, bounds[i].second))
        bounds[i].second = std::numeric_limits<double>::infinity();
    }

    _extender = std::make_shared<RationalMIPExtender>(integrality, bounds);

    struct Visitor
    {
      std::shared_ptr<RationalMIPExtender> extender;

      void operator()(const Constraint<double>& constraint)
      {
        extender->addConstraint(constraint);
      }
    };

    Visitor visitor = { _extender };
    SCIPiterateRows(_scip, _variablesToCoordinates, visitor, true);
  }

  SCIPSolver::~SCIPSolver()
  {
    for (auto& iter : _faceConstraints)
    {
      if (iter.second != nullptr)
        SCIPreleaseCons(_scip, &iter.second);
    }

    SCIPfree(&_scip);
  }

  void SCIPSolver::setFace(std::shared_ptr<Constraint<double>> face)
  {
    if (face == _currentFace)
      return;

    // Remove from SCIP.
    SCIP_CONS* currentCons = _currentFace->isAlwaysSatisfied() ? nullptr
      : _faceConstraints.at(_currentFace);
    if (currentCons != nullptr)
      SCIP_CALL_EXC( SCIPdelCons(_scip, currentCons) );

    if (!face->isAlwaysSatisfied())
    {
      // Create face constraint if it does not exist.

      auto found = _faceConstraints.find(face);
      SCIP_CONS* cons = nullptr;
      if (found == _faceConstraints.end())
      {
        // Add it to SCIP.
        char consName[16];
        SCIPsnprintf(consName, 16, "face#%d", _faceConstraints.size());
        std::vector<SCIP_VAR*> vars;
        std::vector<double> coefficients;
        vars.resize(face->vector().size());
        coefficients.resize(face->vector().size());
        for (const auto& iter : face->vector())
        {
          vars.push_back(_variables[iter.first]);
          coefficients.push_back(iter.second);
        }
        double lhs = face->lhs() == -std::numeric_limits<double>::infinity() ? -SCIPinfinity(_scip) : face->lhs();
        double rhs = face->rhs() == std::numeric_limits<double>::infinity() ? SCIPinfinity(_scip) : face->rhs();
        SCIP_CALL_EXC( SCIPcreateConsBasicLinear(_scip, &cons, consName, vars.size(), &vars[0],
          &coefficients[0], lhs, rhs));
        _faceConstraints.insert(std::make_pair(face, cons));
      }
      else
        cons = found->second;

      // Add to SCIP.
      if (cons != nullptr)
        SCIP_CALL_EXC( SCIPaddCons(_scip, cons) );
    }

    _currentFace = face;
  }

// #if defined(IPO_WITH_GMP) && defined(IPO_WITH_SOPLEX)
//
//   void SCIPSolver::makeRational(OptimizationOracle::Result& result,
//     const mpq_class* objectiveVector)
//   {
//     _makeRationalSolver->setObjective(objectiveVector);
//     _makeRationalSolver->solve(result);
//   }
//
//   void SCIPSolver::makeRational(OptimizationOracle::Result& result, const double* objectiveVector)
//   {
//     _makeRationalSolver->setObjective(objectiveVector);
//     _makeRationalSolver->solve(result);
//   }
//
// #endif

  SCIPOptimizationOracleDouble::SCIPOptimizationOracleDouble(std::shared_ptr<SCIPSolver> solver,
    std::shared_ptr<Constraint<double>> face)
    : OptimizationOracle<double>(solver->name()), _solver(solver), _face(face)
  {
    _space = solver->space();
  }

  SCIPOptimizationOracleDouble::~SCIPOptimizationOracleDouble()
  {

  }

  OptimizationOracle<double>::Result SCIPOptimizationOracleDouble::maximize(const double* objectiveVector,
      const OptimizationOracle<double>::Query& query)
  {
    OptimizationOracle<double>::Result result;

    _solver->setFace(_face);

    SCIP_CALL_EXC( SCIPsetRealParam(_solver->_scip, "limits/time",
      query.timeLimit == std::numeric_limits<double>::infinity() ? SCIPinfinity(_solver->_scip)
      : query.timeLimit) );

    SCIP_CALL_EXC( SCIPsetIntParam(_solver->_scip, "limits/solutions", query.maxNumSolutions) );
    SCIP_CALL_EXC( SCIPsetObjlimit(_solver->_scip,
      isMinusInfinity(query.minObjectiveValue) ? query.minObjectiveValue
      : -SCIPinfinity(_solver->_scip) ) );

    std::size_t n = space()->dimension();

    for (std::size_t i = 0; i < n; ++i)
    {
      SCIP_CALL_EXC( SCIPchgVarObj(_solver->_scip, _solver->_variables[i], objectiveVector[i]) );
    }

    int oldMaxRounds;
    SCIP_CALL_EXC( SCIPgetIntParam(_solver->_scip, "presolving/maxrounds", &oldMaxRounds) );

    for (int attempt = 1; attempt <= 2; ++attempt)
    {
      SCIP_CALL_EXC( SCIPsolve(_solver->_scip) );

      bool hasRay = SCIPhasPrimalRay(_solver->_scip);
      if (hasRay)
      {
        auto vector = std::make_shared<sparse_vector<double>>();
        result.dualBound = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < n; ++i)
        {
          double y = SCIPgetPrimalRayVal(_solver->_scip, _solver->_variables[i]);
          if (!SCIPisZero(_solver->_scip, y))
            vector->push_back(i, y);
        }
        result.rays.push_back(OptimizationOracle::Result::Ray(vector));
        break;
      }

      SCIP_STATUS status = SCIPgetStatus(_solver->_scip);
      if (status == SCIP_STATUS_UNBOUNDED && attempt == 2 && !hasRay)
      {
        throw std::runtime_error("SCIPOptimizationOracle: SCIP reports unboundedness without ray!");
      }

      std::size_t numSolutions = SCIPgetNSols(_solver->_scip);
      if (status != SCIP_STATUS_INFEASIBLE && status != SCIP_STATUS_INFORUNBD
        && status != SCIP_STATUS_UNBOUNDED && numSolutions > 0)
      {
        result.dualBound = SCIPgetDualbound(_solver->_scip);
        SCIP_SOL** solutions = SCIPgetSols(_solver->_scip);
        for (std::size_t solIndex = 0; solIndex < numSolutions; ++solIndex)
        {
          SCIP_SOL* sol = solutions[solIndex];
          double objectiveValue = SCIPgetSolOrigObj(_solver->_scip, sol);
          auto vector = std::make_shared<sparse_vector<double>>();
          if (objectiveValue > query.minObjectiveValue)
          {
            for (std::size_t i = 0; i < n; ++i)
            {
              double x = SCIPgetSolVal(_solver->_scip, sol, _solver->_variables[i]);
              if (!SCIPisZero(_solver->_scip, x))
                vector->push_back(i, x);
            }
            result.points.push_back(OptimizationOracle<double>::Result::Point(vector,
              objectiveValue));
          }
        }
        result.primalBound = SCIPgetPrimalbound(_solver->_scip);
        break;
      }
      else if (status == SCIP_STATUS_INFEASIBLE)
      {
        assert(numSolutions == 0);
        result.dualBound = -std::numeric_limits<double>::infinity();
        result.primalBound = -std::numeric_limits<double>::infinity();
        break;
      }
      else if (status == SCIP_STATUS_UNBOUNDED)
      {
        result.dualBound = std::numeric_limits<double>::infinity();
        if (SCIPhasPrimalRay(_solver->_scip))
        {
          auto vector = std::make_shared<sparse_vector<double>>();
          for (std::size_t i = 0; i < n; ++i)
          {
            double y = SCIPgetPrimalRayVal(_solver->_scip, _solver->_variables[i]);
            if (!SCIPisZero(_solver->_scip, y))
              vector->push_back(i, y);
          }
          result.rays.push_back(OptimizationOracle<double>::Result::Ray(vector));
        }
      }
      else if (attempt == 2)
      {
        std::stringstream ss;
        ss << "SCIPOptimizationOracle: unhandled SCIP status code " << status << ".";
        throw std::runtime_error(ss.str());
      }

      // Disable presolving for the second round.

      SCIP_CALL_EXC( SCIPsetIntParam(_solver->_scip, "presolving/maxrounds", 0) );
      SCIP_CALL_EXC( SCIPfreeSolve(_solver->_scip, true) );
      SCIP_CALL_EXC( SCIPfreeTransform(_solver->_scip) );
    }

    SCIP_CALL_EXC( SCIPsetIntParam(_solver->_scip, "presolving/maxrounds", oldMaxRounds) );
    SCIP_CALL_EXC( SCIPfreeSolve(_solver->_scip, true) );
    SCIP_CALL_EXC( SCIPfreeTransform(_solver->_scip) );

    if (!result.rays.empty() && result.points.empty())
    {
      result.primalBound = std::numeric_limits<double>::infinity();
      if (query.minObjectiveValue == -std::numeric_limits<double>::infinity())
      {
        // We have to check feasibility.
        
        for (std::size_t i = 0; i < n; ++i)
          SCIP_CALL_EXC( SCIPchgVarObj(_solver->_scip, _solver->_variables[i], 0.0) );

        SCIP_CALL_EXC( SCIPsolve(_solver->_scip) );

        if (SCIPgetStatus(_solver->_scip) == SCIP_STATUS_OPTIMAL)
        {
          SCIP_SOL* sol = SCIPgetBestSol(_solver->_scip);
          auto vector = std::make_shared<sparse_vector<double>>();
          double objectiveValue = 0.0;
          for (std::size_t i = 0; i < n; ++i)
          {
            double x = SCIPgetSolVal(_solver->_scip, sol, _solver->_variables[i]);
            if (!SCIPisZero(_solver->_scip, x))
            {
              vector->push_back(i, x);
              objectiveValue =+ objectiveVector[i] * x;
            }
          }
          result.points.push_back(OptimizationOracle<double>::Result::Point(vector, objectiveValue)); 
        }
        else
        {
          result.rays.clear();
          result.primalBound = -std::numeric_limits<double>::infinity();
        }
        SCIP_CALL_EXC( SCIPfreeSolve(_solver->_scip, true) );
        SCIP_CALL_EXC( SCIPfreeTransform(_solver->_scip) );
      }
    }

    result.sortPoints();

    return result;
  }

  SCIPSeparationOracleDouble::SCIPSeparationOracleDouble(std::shared_ptr<SCIPSolver> solver,
    std::shared_ptr<Constraint<double>> face)
    : SeparationOracle(solver->name()), _solver(solver), _face(face)
  {
    _space = solver->space();
  }

  SCIPSeparationOracleDouble::~SCIPSeparationOracleDouble()
  {

  }

  SeparationOracle<double>::Result SCIPSeparationOracleDouble::getInitial(
    const SeparationOracle<double>::Query& query)
  {
    SeparationOracle<double>::Result result;

    if (query.maxNumInequalities > 0 && !_face->isAlwaysSatisfied())
      result.constraints.push_back(*_face);

    struct Visitor
    {
      std::shared_ptr<SCIPSolver> solver;
      const SCIPSeparationOracleDouble::Query& query;
      SCIPSeparationOracleDouble::Result& result;
      std::size_t iteration;
      std::chrono::time_point<std::chrono::system_clock> started;

      void operator()(const Constraint<double>& constraint)
      {
        // Did we reach a limit?
        if (result.hitTimeLimit || result.constraints.size() == query.maxNumInequalities)
          return;

        iteration = (iteration + 1) % SCIP_SEPARATION_CHECK_TIMELIMIT_FREQUENCY;
        if (iteration == 0)
        {
          std::chrono::duration<double> duration = std::chrono::system_clock::now() - started;
          if (duration.count() > query.timeLimit)
          {
            result.hitTimeLimit = true;
            return;
          }
        }

        result.constraints.push_back(constraint);
      }
    };

    Visitor visitor = { _solver, query, result, 0, std::chrono::system_clock::now() };

    SCIPiterateBounds(_solver->_scip, _solver->_variablesToCoordinates, visitor, true);
    SCIPiterateRows(_solver->_scip, _solver->_variablesToCoordinates, visitor, true);

    return result;
  }

  SeparationOracle<double>::Result SCIPSeparationOracleDouble::separate(const double* vector, bool isPoint,
      const SeparationOracle<double>::Query& query)
  {
    SeparationOracle<double>::Result result;

    _solver->setFace(_face);

    struct Visitor
    {
      std::shared_ptr<SCIPSolver> solver;
      const double* vector;
      bool isPoint;
      const SCIPSeparationOracleDouble::Query& query;
      SCIPSeparationOracleDouble::Result& result;
      std::size_t iteration;
      std::chrono::time_point<std::chrono::system_clock> started;

      void operator()(const Constraint<double>& constraint)
      {
        // Did we reach a limit?
        if (result.hitTimeLimit || result.constraints.size() == query.maxNumInequalities)
          return;

        iteration = (iteration + 1) % SCIP_SEPARATION_CHECK_TIMELIMIT_FREQUENCY;
        if (iteration == 0)
        {
          std::chrono::duration<double> duration = std::chrono::system_clock::now() - started;
          if (duration.count() > query.timeLimit)
          {
            result.hitTimeLimit = true;
            return;
          }
        }

        // Set lhs/rhs to 0 if we are separating a ray.
        double rhs = constraint.rhs();
        if (!isPoint && rhs < std::numeric_limits<double>::infinity())
          rhs = 0.0;
        double lhs = constraint.lhs();
        if (!isPoint && lhs > -std::numeric_limits<double>::infinity())
          lhs = 0.0;

        double activity = 0.0;
        for (const auto& iter : constraint.vector())
          activity += vector[iter.first] * iter.second;

        if (!SCIPisFeasLE(solver->_scip, activity, rhs)
          || !SCIPisFeasGE(solver->_scip, activity, lhs))
        {
          result.constraints.push_back(constraint);
        }
      }
    };

    Visitor visitor = { _solver, vector, isPoint, query, result, 0,
        std::chrono::system_clock::now() };

    SCIPiterateBounds(_solver->_scip, _solver->_variablesToCoordinates, visitor, false);
    SCIPiterateRows(_solver->_scip, _solver->_variablesToCoordinates, visitor, false);

    return result;
  }

} /* namespace ipo */
