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
   * \brief Auxiliary method that iterates over the rows.
   *
   * Auxiliary method that iterates over the rows of linear and setppc constraints, calling
   * \p visitor with each of them. It passes the lower bound, the value 1, a pointer to the
   * variable, a pointer to the coefficient 1, and the upper bound in that order.
   */

  void SCIPiterateBoundsReal(SCIP* scip,
    std::function<void(const Constraint& constraint)> visitor)
  {
    SCIP_VAR** vars = SCIPgetOrigVars(scip);
    std::size_t n = SCIPgetNOrigVars(scip);
    std::vector<std::pair<std::size_t, double>> entries;
    entries.resize(1);
    entries.front().second = 1.0;
    for (std::size_t v = 0; v < n; ++v)
    {
      SCIP_VAR* var = vars[v];
      entries.front().first = v;
      visitor(Constraint(Value(SCIPvarGetLbGlobal(var)), Vector(entries),
        Value(SCIPvarGetUbGlobal(var))));
    }
  }

  void SCIPiterateBoundsRational(SCIP* scip,
    std::function<void(const Constraint& constraint)> visitor)
  {
    SCIP_VAR** vars = SCIPgetOrigVars(scip);
    std::size_t n = SCIPgetNOrigVars(scip);
    std::vector<std::pair<std::size_t, mpq_class>> entries;
    entries.resize(1);
    entries.front().second = 1;
    for (std::size_t v = 0; v < n; ++v)
    {
      SCIP_VAR* var = vars[v];
      entries.front().first = v;
      visitor(Constraint(Value(mpq_class(SCIPvarGetLbGlobal(var))), Vector(entries),
        Value(mpq_class(SCIPvarGetUbGlobal(var)))));
    }
  }

  /**
   * \brief Auxiliary method that iterates over the rows.
   *
   * Auxiliary method that iterates over the rows of linear and setppc constraints, calling
   * \p visitor with each of them. It passes the lhs, number of nonzeros, variables and
   * coefficients of the nonzeros, and the rhs of each row in that order.
   */

  void SCIPiterateRowsReal(SCIP* scip,
    const std::unordered_map<SCIP_VAR*, std::size_t>& variablesToCoordinates,
    std::function<void(const Constraint& constraint)> visitor)
  {
    SCIP_CONS** conss = SCIPgetConss(scip);
    std::size_t m = SCIPgetNConss(scip);
    std::vector<double> ones;
    ones.resize(SCIPgetNOrigVars(scip), 1.0);
    std::vector<std::pair<std::size_t, double>> entries;
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

      entries.resize(k);
      for (std::size_t i = 0; i < k; ++i)
        entries[i] = std::pair<std::size_t, double>(variablesToCoordinates.at(vars[i]), vals[i]);
      visitor(Constraint(Value(lhs), Vector(entries), Value(rhs)));
    }
  }

  void SCIPiterateRowsRational(SCIP* scip,
    const std::unordered_map<SCIP_VAR*, std::size_t>& variablesToCoordinates,
    std::function<void(const Constraint& constraint)> visitor)
  {
    SCIP_CONS** conss = SCIPgetConss(scip);
    std::size_t m = SCIPgetNConss(scip);
    std::vector<double> ones;
    ones.resize(SCIPgetNOrigVars(scip), 1.0);
    std::vector<std::pair<std::size_t, mpq_class>> entries;
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

      Value valLhs = SCIPisInfinity(scip, rhs) ? plusInfinity() : Value(mpq_class(lhs));
      Value valRhs = SCIPisInfinity(scip, -lhs) ? minusInfinity() : Value(mpq_class(rhs));

      entries.resize(k);
      for (std::size_t i = 0; i < k; ++i)
        entries[i] = std::pair<std::size_t, mpq_class>(variablesToCoordinates.at(vars[i]), mpq_class(vals[i]));
      std::sort(entries.begin(), entries.end());
      Vector vector(entries, true);
      visitor(Constraint(valLhs, vector, valRhs));
    }
  }

  SCIPSolver::SCIPSolver(SCIP* scip)
    : _scip(scip), _currentFace(Value(), Vector(), Value())
  {
    initialize();
  }

  SCIPSolver::SCIPSolver(const std::string& fileName)
    : _currentFace(Value(), Vector(), Value())
  {
    SCIP_CALL_EXC( SCIPcreate(&_scip) );
    SCIP_CALL_EXC( SCIPincludeDefaultPlugins(_scip) );
//     SCIP_CALL_EXC(SCIPsetIntParam(_scip, "display/verblevel", 0) );
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
    
#if defined(IPO_WITH_GMP) && defined(IPO_WITH_SOPLEX)

    std::vector<bool> integrality(n);
    for (std::size_t i = 0; i < n; ++i)
      integrality[i] = SCIPvarGetType(vars[i]) != SCIP_VARTYPE_CONTINUOUS;

    // Extract all constraints from the MIP.
    std::vector<Constraint> constraints;

    struct Visitor
    {
      std::vector<Constraint>* constraints;

      void operator()(const Constraint& constraint)
      {
        constraints->push_back(constraint);
      }
    };

    Visitor visitor = { &constraints };

    SCIPiterateBoundsRational(_scip, visitor);
    SCIPiterateRowsRational(_scip, _variablesToCoordinates, visitor);
    
    _makeRationalSolver = new MakeRationalSolver(integrality, constraints);

#endif /* IPO_WITH_GMP && IPO_WITH_SOPLEX */
  }

  SCIPSolver::~SCIPSolver()
  {
#if defined(IPO_WITH_GMP) && defined(IPO_WITH_SOPLEX)
    delete _makeRationalSolver;
#endif /* IPO_WITH_GMP && IPO_WITH_SOPLEX */

    for (auto& iter : _faceConstraints)
    {
      if (iter.second != nullptr)
        SCIPreleaseCons(_scip, &iter.second);
    }

    SCIPfree(&_scip);
  }

//   SCIPSolver::Face SCIPSolver::addFace(std::size_t numNonzeros, std::size_t* coordinates,
//     double* coefficients, double rhs)
//   {
//     _faces.push_back(FaceData());
//     FaceData& faceData = _faces.back();
// 
//     // Add it to SCIP.
//     char consName[16];
//     SCIPsnprintf(consName, 16, "face%d", _faces.size() - 1);
//     std::vector<SCIP_VAR*> vars;
//     vars.resize(numNonzeros);
//     for (std::size_t i = 0; i < numNonzeros; ++i)
//       vars[i] = _variables[coordinates[i]];
//     SCIP_CALL_EXC( SCIPcreateConsBasicLinear(_scip, &faceData.constraint, consName, numNonzeros,
//       &vars[0], &coefficients[0], -SCIPinfinity(_scip), rhs));
// 
//     faceData.numNonzeros = numNonzeros;
//     faceData.coordinates = new std::size_t[numNonzeros];
//     faceData.coefficients = new double[numNonzeros];
//     faceData.rhs = rhs;
// #if defined(IPO_WITH_GMP)
//     faceData.rationalCoefficients = new mpq_class[numNonzeros];
//     faceData.rationalRhs = rhs;
// #endif /* IPO_WITH_GMP */
//     for (std::size_t i = 0; i < numNonzeros; ++i)
//     {
//       faceData.coordinates[i] = coordinates[i];
//       faceData.coefficients[i] = coefficients[i];
// #if defined(IPO_WITH_GMP)
//       faceData.rationalCoefficients[i] = coefficients[i];
// #endif /* IPO_WITH_GMP */
//     }
// 
//     return _faces.size() - 1;
//   }
// 
// #if defined(IPO_WITH_GMP)
// 
//   SCIPSolver::Face SCIPSolver::addFace(std::size_t numNonzeros, std::size_t* coordinates,
//     mpq_class* coefficients, const mpq_class& rhs)
//   {
//     double* approxCoefficients = new double[numNonzeros];
//     double approxRhs = rhs.get_d();
// 
//     for (std::size_t i = 0; i < numNonzeros; ++i)
//       approxCoefficients[i] = coefficients[i].get_d();
// 
//     Face face = addFace(numNonzeros, coordinates, approxCoefficients, approxRhs);
//     _faces.back().rationalRhs = rhs;
//     for (std::size_t i = 0; i < numNonzeros; ++i)
//       _faces.back().rationalCoefficients[i] = coefficients[i];
// 
//     delete[] approxCoefficients;
// 
//     return face;
//   }
// 
// #endif

  void SCIPSolver::setFace(const Constraint& face)
  {
    if (face == _currentFace)
      return;
    

    // Remove from SCIP.
    SCIP_CONS* currentCons = _currentFace.isAlwaysSatisfied() ? nullptr 
      : _faceConstraints.at(_currentFace);
    if (currentCons != nullptr)
      SCIP_CALL_EXC( SCIPdelCons(_scip, currentCons) );

#if defined(IPO_WITH_GMP) && defined (IPO_WITH_SOPLEX)
    // Remove from makeRational.
    if (currentCons != nullptr)
    {
      _makeRationalSolver->removeLastRow();
    }
#endif /* IPO_WITH_GMP && IPO_WITH_SOPLEX */

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
      vars.resize(face.vector.size());
      coefficients.resize(face.vector.size());
      for (std::size_t i = 0; i < vars.size(); ++i)
      {
        vars[i] = _variables[face.vector.coordinate(i)];
        coefficients[i] = face.vector.real(i);
      }
      double lhs = face.lhs.isMinusInfinity() ? -SCIPinfinity(_scip) : face.lhs.real;
      double rhs = face.rhs.isPlusInfinity() ? SCIPinfinity(_scip) : face.rhs.real;
      SCIP_CALL_EXC( SCIPcreateConsBasicLinear(_scip, &cons, consName, face.vector.size(), &vars[0],
        &coefficients[0], lhs, rhs));
      _faceConstraints.insert(std::make_pair(face, cons));
    }
    else
      cons = found->second;

    // Add to SCIP.
    if (cons != nullptr)
      SCIP_CALL_EXC( SCIPaddCons(_scip, cons) );

#if defined(IPO_WITH_GMP) && defined(IPO_WITH_SOPLEX)
    // Add to makeRational.
    if (cons != nullptr)
      _makeRationalSolver->addRow(face);
#endif /* IPO_WITH_GMP && IPO_WITH_SOPLEX */

    _currentFace = face;
  }

#if defined(IPO_WITH_GMP) && defined(IPO_WITH_SOPLEX)

  void SCIPSolver::makeRational(OptimizationOracle::Result& result,
    const mpq_class* objectiveVector)
  {
    _makeRationalSolver->setObjective(objectiveVector);
    _makeRationalSolver->solve(result);
  }

  void SCIPSolver::makeRational(OptimizationOracle::Result& result, const double* objectiveVector)
  {
    _makeRationalSolver->setObjective(objectiveVector);
    _makeRationalSolver->solve(result);
  }

#endif

  SCIPOptimizationOracle::SCIPOptimizationOracle(std::shared_ptr<SCIPSolver> solver,
    Constraint face)
    : OptimizationOracle(solver->name()), _solver(solver), _face(face)
  {
    _space = solver->space();
  }

  SCIPOptimizationOracle::~SCIPOptimizationOracle()
  {

  }

  bool SCIPOptimizationOracle::isExact() const
  {
#if defined(IPO_WITH_SOPLEX) && defined(IPO_WITH_GMP)
    return true;
#else
    return false;
#endif /* IPO_WITH_SOPLEX */
  }

  void SCIPOptimizationOracle::maximize(const double* objectiveVector, const Query& query,
    Result& result)
  {
    _solver->setFace(_face);

    SCIP_CALL_EXC( SCIPsetRealParam(_solver->_scip, "limits/time",
      query.timeLimit == std::numeric_limits<double>::infinity() ? SCIPinfinity(_solver->_scip)
      : query.timeLimit) );

    double oldObjectiveLimit = SCIPgetObjlimit(_solver->_scip);
    int oldSolutionLimit;
    SCIP_CALL_EXC( SCIPgetIntParam(_solver->_scip, "limits/solutions", &oldSolutionLimit) );
    SCIP_CALL_EXC( SCIPsetObjlimit(_solver->_scip,
      query.minObjectiveValue > -std::numeric_limits<double>::infinity() ? query.minObjectiveValue
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
        std::vector<std::pair<std::size_t, double>> entries;
        result.dualBound = plusInfinity();
        for (std::size_t i = 0; i < n; ++i)
        {
          double y = SCIPgetPrimalRayVal(_solver->_scip, _solver->_variables[i]);
          if (!SCIPisZero(_solver->_scip, y))
            entries.push_back(std::make_pair(i, y));
        }
        result.rays.push_back(OptimizationOracle::Result::Ray(
          Vector(entries)));
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
        std::vector<std::pair<std::size_t, double>> entries;
        for (std::size_t solIndex = 0; solIndex < numSolutions; ++solIndex)
        {
          SCIP_SOL* sol = solutions[solIndex];
          entries.clear();
          for (std::size_t i = 0; i < n; ++i)
          {
            double x = SCIPgetSolVal(_solver->_scip, sol, _solver->_variables[i]);
            if (!SCIPisZero(_solver->_scip, x))
              entries.push_back(std::make_pair(i, x));
          }
          result.points.push_back(OptimizationOracle::Result::Point(Vector(entries),
            Value(SCIPgetSolOrigObj(_solver->_scip, sol))));
        }
        break;
      }
      else if (status == SCIP_STATUS_INFEASIBLE)
      {
        assert(numSolutions == 0);
        result.dualBound = minusInfinity();
        break;
      }
      else if (status == SCIP_STATUS_UNBOUNDED)
      {
        result.dualBound = plusInfinity();
        if (SCIPhasPrimalRay(_solver->_scip))
        {
          std::vector<std::pair<std::size_t, double>> entries;
          for (std::size_t i = 0; i < n; ++i)
          {
            double y = SCIPgetPrimalRayVal(_solver->_scip, _solver->_variables[i]);
            if (!SCIPisZero(_solver->_scip, y))
              entries.push_back(std::make_pair(i, y));
          }
          result.rays.push_back(Vector(entries));
          break;
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
    SCIP_CALL_EXC( SCIPsetObjlimit(_solver->_scip, oldObjectiveLimit) );
    SCIP_CALL_EXC( SCIPsetIntParam(_solver->_scip, "limits/solutions", oldSolutionLimit) );
    
#if defined(IPO_WITH_GMP) && defined(IPO_WITH_SOPLEX)
    if (query.rational)
      _solver->makeRational(result, objectiveVector);
#endif /* IPO_WITH_GMP && IPO_WITH_SOPLEX */
    
    result.sortPoints();
  }

#ifdef IPO_WITH_GMP

  void SCIPOptimizationOracle::maximize(const mpq_class* objectiveVector,
    const OptimizationOracle::Query& query, OptimizationOracle::Result& result)
  {
    // Create floating-point approximation of objective vector.
    double* approximateObjectiveVector = new double[space()->dimension()];
    for (std::size_t i = 0; i < space()->dimension(); ++i)
      approximateObjectiveVector[i] = objectiveVector[i].get_d();

#if defined(IPO_WITH_SOPLEX)
    if (query.rational)
    {
      OptimizationOracle::Query noRationalQuery = query;
      noRationalQuery.rational = false;
      this->maximize(approximateObjectiveVector, noRationalQuery, result);
      this->_solver->makeRational(result, objectiveVector);
    }
    else
#endif /* IPO_WITH_SOPLEX */
    {
      this->maximize(approximateObjectiveVector, query, result);
    }

    delete[] approximateObjectiveVector;

    result.sortPoints();
  }

#endif /* IPO_WITH_GMP */

  SCIPSeparationOracle::SCIPSeparationOracle(std::shared_ptr<SCIPSolver> solver,
    Constraint face)
    : SeparationOracle(solver->name()), _solver(solver), _face(face)
  {
    _space = solver->space();
  }

  SCIPSeparationOracle::~SCIPSeparationOracle()
  {

  }

  bool SCIPSeparationOracle::separate(const double* vector, bool isPoint, const Query& query,
    Result& result)
  {
    _solver->setFace(_face);

    struct Visitor
    {
      std::shared_ptr<SCIPSolver> solver;
      const double* vector;
      bool isPoint;
      const SCIPSeparationOracle::Query& query;
      SCIPSeparationOracle::Result& result;
      std::size_t iteration;
      std::chrono::time_point<std::chrono::system_clock> started;

      void operator()(const Constraint& constraint)
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
        double rhs = constraint.rhs.real;
        if (!isPoint && rhs < std::numeric_limits<double>::infinity())
          rhs = 0.0;
        double lhs = constraint.lhs.real;
        if (!isPoint && lhs > -std::numeric_limits<double>::infinity())
          lhs = 0.0;

        double activity = 0.0;
        for (std::size_t i = 0; i < constraint.vector.size(); ++i)
        {
          activity += vector[constraint.vector.coordinate(i)] * constraint.vector.real(i);
        }

        if (!SCIPisFeasLE(solver->_scip, activity, rhs)
          || !SCIPisFeasGE(solver->_scip, activity, lhs))
        {
          result.constraints.push_back(constraint);
        }
      }
    };

    Visitor visitor = { _solver, vector, isPoint, query, result, 0,
        std::chrono::system_clock::now() };

#if defined(IPO_WITH_GMP)
    if (query.rational)
    {
      SCIPiterateBoundsRational(_solver->_scip, visitor);
      SCIPiterateRowsRational(_solver->_scip, _solver->_variablesToCoordinates, visitor);
    }
    else
    {
      SCIPiterateBoundsReal(_solver->_scip, visitor);
      SCIPiterateRowsReal(_solver->_scip, _solver->_variablesToCoordinates, visitor);
    }
#else
    SCIPiterateBoundsReal(_solver->_scip, visitor);
    SCIPiterateRowsReal(_solver->_scip, _solver->_variablesToCoordinates, visitor);
#endif /* IPO_WITH_GMP */

    return !result.constraints.empty();
  }

#if defined(IPO_WITH_GMP)

  bool SCIPSeparationOracle::separate(const mpq_class* vector, bool isPoint, const Query& query,
    Result& result)
  {
    _solver->setFace(_face);

    struct Visitor
    {
      std::shared_ptr<SCIPSolver> solver;
      const mpq_class* vector;
      bool isPoint;
      const SCIPSeparationOracle::Query& query;
      SCIPSeparationOracle::Result& result;
      std::size_t iteration;
      std::chrono::time_point<std::chrono::system_clock> started;

      void operator()(const Constraint& constraint)
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
        double rhs = constraint.rhs.real;
        if (!isPoint && rhs < std::numeric_limits<double>::infinity())
          rhs = 0.0;
        double lhs = constraint.lhs.real;
        if (!isPoint && lhs > -std::numeric_limits<double>::infinity())
          lhs = 0.0;

        mpq_class activity = 0;
        for (std::size_t i = 0; i < constraint.vector.size(); ++i)
        {
          activity += vector[constraint.vector.coordinate(i)] * constraint.vector.real(i);
        }

        if (activity < lhs || activity > rhs)
          result.constraints.push_back(constraint);
      }
    };

    Visitor visitor = { _solver, vector, isPoint, query, result, 0,
        std::chrono::system_clock::now() };

#if defined(IPO_WITH_GMP)
    if (query.rational)
    {
      SCIPiterateBoundsRational(_solver->_scip, visitor);
      SCIPiterateRowsRational(_solver->_scip, _solver->_variablesToCoordinates, visitor);
    }
    else
    {
      SCIPiterateBoundsReal(_solver->_scip, visitor);
      SCIPiterateRowsReal(_solver->_scip, _solver->_variablesToCoordinates, visitor);
    }
#else
    SCIPiterateBoundsReal(_solver->_scip, visitor);
    SCIPiterateRowsReal(_solver->_scip, _solver->_variablesToCoordinates, visitor);
#endif /* IPO_WITH_GMP */

    return !result.constraints.empty();
  }

#endif /* IPO_WITH_GMP */

} /* namespace ipo */
