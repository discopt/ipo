// #define IPO_DEBUG // Uncomment to debug this file.

#include <ipo/constraint.hpp>
#include <ipo/oracles_scip.hpp>

// TODO: Remove objective limit and add event handler for primal and dual bound.

#include <cassert>
#include <functional>
#include <chrono>
#include <sstream>
#include <iostream>
#include <unordered_map>

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

#define EVENTHDLR_NAME "boundchange"
#define EVENTHDLR_DESC "event handler for primal or dual bound change"

extern "C"
{
  struct SCIP_EventhdlrData
  {
    double minPrimalBound;
    double maxDualBound;
  };
  
  /** copy method for event handler plugins (called when SCIP copies plugins) */
  static SCIP_DECL_EVENTCOPY(eventCopyGlobalBoundChange) 
  {
    return SCIP_OKAY;
  }

  /** initialization method of event handler (called after problem was transformed) */
  static SCIP_DECL_EVENTINIT(eventInitGlobalBoundChange)
  {
    assert(scip != NULL);
    assert(eventhdlr != NULL);
    assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

    /* notify SCIP that your event handler wants to react on the event type best solution found */
    SCIP_CALL( SCIPcatchEvent( scip,
      SCIP_EVENTTYPE_BESTSOLFOUND | SCIP_EVENTTYPE_NODESOLVED | SCIP_EVENTTYPE_LPSOLVED, eventhdlr,
      NULL, NULL) );

    return SCIP_OKAY;
  }

  /** deinitialization method of event handler (called before transformed problem is freed) */
  static SCIP_DECL_EVENTEXIT(eventExitGlobalBoundChange)
  {
    assert(scip != NULL);
    assert(eventhdlr != NULL);
    assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

    /* notify SCIP that your event handler wants to drop the event type best solution found */
    SCIP_CALL( SCIPdropEvent( scip,
      SCIP_EVENTTYPE_BESTSOLFOUND | SCIP_EVENTTYPE_NODESOLVED | SCIP_EVENTTYPE_LPSOLVED, eventhdlr,
      NULL, -1) );

    return SCIP_OKAY;
  }

  /** execution method of event handler */
  static SCIP_DECL_EVENTEXEC(eventExecGlobalBoundChange)
  {
    assert(scip);

    SCIP_EVENTHDLRDATA* eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
    if ((eventhdlrdata->minPrimalBound > -std::numeric_limits<double>::infinity()
      && SCIPgetPrimalbound(scip) > eventhdlrdata->minPrimalBound
      && SCIPisLT(scip, SCIPgetDualbound(scip), SCIPinfinity(scip)))
      || (eventhdlrdata->maxDualBound < std::numeric_limits<double>::infinity()
      && SCIPisLE(scip, SCIPgetDualbound(scip), eventhdlrdata->maxDualBound)))
    {
#if defined(IPO_DEBUG)
      std::cout << "Setting node limit to current number of nodes to stop optimization." << std::endl;
#endif /* IPO_DEBUG */
      SCIP_CALL( SCIPsetLongintParam(scip, "limits/totalnodes", SCIPgetNTotalNodes(scip)) );
    }

    return SCIP_OKAY;
  }

}

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

  static
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
      double rhs = SCIPvarGetUbOriginal(vars[v]);
      std::size_t coordinate = variablesToCoordinates.at(vars[v]);
      if (SCIPisEQ(scip, lhs, rhs))
      {
        double average = (lhs + rhs) / 2.0;
        auto vector = std::make_shared<sparse_vector<double>>();
        vector->push_back(coordinate, 1.0);
        visitor(Constraint<double>(average, vector, average, ConstraintType::EQUATION));
        continue;
      }
      if (SCIPisInfinity(scip, -lhs) && SCIPisInfinity(scip, rhs))
        continue;
      if (SCIPisInfinity(scip, -lhs) || (!SCIPisInfinity(scip, rhs) && !ranged))
      {
        auto vector = std::make_shared<sparse_vector<double>>();
        vector->push_back(coordinate, 1.0);
        visitor(Constraint<double>(vector, rhs));
      }
      if (SCIPisInfinity(scip, rhs) || (!SCIPisInfinity(scip, -lhs) && !ranged))
      {
        auto vector = std::make_shared<sparse_vector<double>>();
        vector->push_back(coordinate, 1.0);
        visitor(Constraint<double>(lhs, vector));
      }
      if (ranged && !SCIPisInfinity(scip, -lhs) && !SCIPisInfinity(scip, rhs))
      {
        auto vector = std::make_shared<sparse_vector<double>>();
        vector->push_back(coordinate, 1.0);
        visitor(Constraint<double>(lhs, vector, rhs, ConstraintType::RANGED));
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

  static
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

      if (SCIPisEQ(scip, lhs, rhs))
      {
        double average = (lhs + rhs) / 2.0;
        unsortedEntries.clear();
        for (std::size_t i = 0; i < k; ++i)
          unsortedEntries.push_back(std::make_pair(variablesToCoordinates.at(vars[i]), vals[i]));
        auto vector = std::make_shared<sparse_vector<double>>(std::move(unsortedEntries), true);
        visitor(Constraint<double>(average, vector, average, ConstraintType::EQUATION));
        continue;
      }
      if (SCIPisInfinity(scip, -lhs) && SCIPisInfinity(scip, rhs))
        continue;
      if (SCIPisInfinity(scip, -lhs) || (!SCIPisInfinity(scip, rhs) && !ranged))
      {
        unsortedEntries.clear();
        for (std::size_t i = 0; i < k; ++i)
          unsortedEntries.push_back(std::make_pair(variablesToCoordinates.at(vars[i]), vals[i]));
        auto vector = std::make_shared<sparse_vector<double>>(std::move(unsortedEntries), true);
        visitor(Constraint<double>(vector, rhs));
      }
      if (SCIPisInfinity(scip, rhs) || (!SCIPisInfinity(scip, -lhs) && !ranged))
      {
        unsortedEntries.clear();
        for (std::size_t i = 0; i < k; ++i)
          unsortedEntries.push_back(std::make_pair(variablesToCoordinates.at(vars[i]), vals[i]));
        auto vector = std::make_shared<sparse_vector<double>>(std::move(unsortedEntries), true);
        visitor(Constraint<double>(lhs, vector));
      }
      if (ranged && !SCIPisInfinity(scip, -lhs) && !SCIPisInfinity(scip, rhs))
      {
        unsortedEntries.clear();
        for (std::size_t i = 0; i < k; ++i)
          unsortedEntries.push_back(std::make_pair(variablesToCoordinates.at(vars[i]), vals[i]));
        auto vector = std::make_shared<sparse_vector<double>>(std::move(unsortedEntries), true);
        visitor(Constraint<double>(lhs, vector, rhs, ConstraintType::RANGED));
      }
    }
  }

  SCIPSolver::SCIPSolver(SCIP*&& scip)
    : _scip(scip), _currentFace(NULL)
  {
    scip = NULL;
#if !defined(IPO_DEBUG)
    SCIP_CALL_EXC( SCIPsetIntParam(_scip, "display/verblevel", 0) );
#endif /* !IPO_DEBUG */
    initialize();
  }

  SCIPSolver::SCIPSolver(const std::string& fileName)
    : _currentFace(NULL)
  {
    SCIP_CALL_EXC( SCIPcreate(&_scip) );
    SCIP_CALL_EXC( SCIPincludeDefaultPlugins(_scip) );
#if !defined(IPO_DEBUG)
    SCIP_CALL_EXC( SCIPsetIntParam(_scip, "display/verblevel", 0) );
#endif /* !IPO_DEBUG */
    SCIP_CALL_EXC( SCIPreadProb(_scip, fileName.c_str(), NULL) );

    initialize();
  }

  void SCIPSolver::initialize()
  {
    SCIP_CALL_EXC( SCIPsetBoolParam(_scip, "misc/catchctrlc", false) );
    std::size_t n = SCIPgetNOrigVars(_scip);
    _variables.resize(n);
    std::vector<std::string> variableNames;
    variableNames.resize(n);
    _instanceObjective = new double[n+1];
    SCIP_VAR** vars = SCIPgetOrigVars(_scip);
    double scale = SCIPgetObjsense(_scip) == SCIP_OBJSENSE_MAXIMIZE ? 1.0 : -1.0;
    for (std::size_t i = 0; i < n; ++i)
    {
      _variables[i] = vars[i];
      variableNames[i] = SCIPvarGetName(vars[i]);
      _variablesToCoordinates[vars[i]] = i;
      _instanceObjective[i] = scale * SCIPvarGetObj(vars[i]);
    }
    _instanceObjective[n] = scale * SCIPgetOrigObjoffset(_scip);
    SCIP_CALL_EXC( SCIPsetObjsense(_scip, SCIP_OBJSENSE_MAXIMIZE) );
    SCIP_CALL_EXC( SCIPaddOrigObjoffset(_scip, -SCIPgetOrigObjoffset(_scip)) );

    _name = SCIPgetProbName(_scip);
    _space = std::make_shared<Space>(std::move(variableNames));

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

#if defined(IPO_RATIONAL_LP)
    _extender = new RationalMIPExtender(integrality, bounds);

    struct Visitor
    {
      RationalMIPExtender* extender;

      void operator()(const Constraint<double>& constraint)
      {
        extender->addConstraint(constraint);
      }
    };

    Visitor visitor = { _extender };
    SCIPiterateRows(_scip, _variablesToCoordinates, visitor, true);
#endif /* IPO_RATIONAL_LP */
    
    SCIP_EVENTHDLR* eventhdlr = NULL;
    SCIP_CALL_EXC( SCIPincludeEventhdlrBasic(_scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
      eventExecGlobalBoundChange, (SCIP_EVENTHDLRDATA*)(&_boundLimits)) );
    assert(eventhdlr != NULL);

    SCIP_CALL_EXC( SCIPsetEventhdlrCopy(_scip, eventhdlr, eventCopyGlobalBoundChange) );
    SCIP_CALL_EXC( SCIPsetEventhdlrInit(_scip, eventhdlr, eventInitGlobalBoundChange) );
    SCIP_CALL_EXC( SCIPsetEventhdlrExit(_scip, eventhdlr, eventExitGlobalBoundChange) );
  }

  SCIPSolver::~SCIPSolver()
  {
    delete[] _instanceObjective;

#if defined(IPO_RATIONAL_LP)
    delete _extender;
#endif /* IPO_RATIONAL_LP */

    for (auto& iter : _faceConstraints)
    {
      if (iter.second != nullptr)
        SCIPreleaseCons(_scip, &iter.second);
    }

    SCIPfree(&_scip);
  }

  void SCIPSolver::addFace(Constraint<double>* face)
  {
    assert(face);
    if (face->isAlwaysSatisfied())
      return;

    if (_faceConstraints.find(face) != _faceConstraints.end())
      throw std::runtime_error("SCIPSolver::addFace() failed: the face is already known.");

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
    double lhs, rhs;
    if (face->type() == ConstraintType::LESS_OR_EQUAL)
    {
      lhs = face->rhs();
      rhs = face->rhs();
    }
    else if (face->type() == ConstraintType::GREATER_OR_EQUAL)
    {
      lhs = face->lhs();
      rhs = face->lhs();
    }
    else if (face->type() == ConstraintType::EQUATION)
    {
      lhs = face->lhs();
      rhs = face->rhs();
      assert(SCIPisEQ(_scip, lhs, rhs));
    }
    else
    {
      throw std::runtime_error("Cannot use a ranged constraint or equation to define a face.");
    }
    SCIP_CONS* cons = NULL;
    SCIP_CALL_EXC( SCIPcreateConsBasicLinear(_scip, &cons, consName, vars.size(), &vars[0],
      &coefficients[0], lhs, rhs));
    _faceConstraints.insert(std::make_pair(face, cons));

#if defined(IPO_DEBUG)
    printf("Created new face constraint for face %p with lhs %f and rhs %f:\n", face, face->lhs(), face->rhs());
    SCIPprintCons(_scip, cons, stdout);
    fflush(stdout);
#endif /* IPO_DEBUG */
  }

  void SCIPSolver::deleteFace(Constraint<double>* face)
  {
    assert(face);
    if (face->isAlwaysSatisfied())
      return;

    if (face == _currentFace)
      selectFace(NULL);

    auto iter = _faceConstraints.find(face);
    if (iter == _faceConstraints.end())
      throw std::runtime_error("SCIPSolver::deleteFace() failed: the face is not known.");

    SCIP_CALL_EXC( SCIPreleaseCons(_scip, &iter->second) );
    _faceConstraints.erase(iter);
  }

  void SCIPSolver::selectFace(Constraint<double>* face)
  {
    if (face && face->isAlwaysSatisfied())
      face = NULL;

    if (face == _currentFace)
    {
#if defined(IPO_DEBUG)
      std::cout << "Face is already selected." << std::endl;
#endif /* IPO_DEBUG */
      return;
    }

    // Remove from SCIP.
    if (_currentFace)
    {
#if defined(IPO_DEBUG)
      std::cout << "Disabling old face." << std::endl;
#endif /* IPO_DEBUG */

      SCIP_CALL_EXC( SCIPdelCons(_scip, _faceConstraints.at(_currentFace)) );
    }
    if (face)
    {
#if defined(IPO_DEBUG)
      printf("Enabling new face %p with lhs %f and rhs %f:\n", face, face->lhs(), face->rhs());
      SCIPprintCons(_scip, _faceConstraints.at(face), stdout);
      fflush(stdout);
#endif /* IPO_DEBUG */
      SCIP_CALL_EXC( SCIPaddCons(_scip, _faceConstraints.at(face)) );
    }

    _currentFace = face;
  }

  /**
   * \brief Returns a double-arithmetic optimization oracle for the requested \p face.
   */

  template <>
  IPO_EXPORT
  std::shared_ptr<SCIPOptimizationOracle<double>> SCIPSolver::getOptimizationOracle<double>(
    const Constraint<double>& face)
  {
    return std::make_shared<SCIPOptimizationOracle<double>>(shared_from_this(), face);
  }

#if defined(IPO_RATIONAL_LP)

  /**
   * \brief Returns a rational optimization oracle for the requested \p face.
   **/

  template <>
  IPO_EXPORT
  std::shared_ptr<SCIPOptimizationOracle<rational>> SCIPSolver::getOptimizationOracle<rational>(
    const Constraint<rational>& face)
  {
    auto approximateFace = convertConstraint<double>(face);
    auto approximateOracle = getOptimizationOracle<double>(approximateFace);
    return std::make_shared<SCIPOptimizationOracle<rational>>(_extender, approximateOracle, face);
  }

#endif /* IPO_RATIONAL_LP */

  SCIPOptimizationOracle<double>::SCIPOptimizationOracle(std::shared_ptr<SCIPSolver> solver,
    const Constraint<double>& face)
    : OptimizationOracle<double>(solver->name()), _solver(solver), _face(face)
  {
    _space = solver->space();
    _solver->addFace(&_face);
  }

  SCIPOptimizationOracle<double>::~SCIPOptimizationOracle()
  {
    _solver->deleteFace(&_face);
  }

  /**
   * \brief Returns a double arithmetic separation oracle for the \p face.
   **/

  template <>
  IPO_EXPORT
  std::shared_ptr<SCIPSeparationOracle<double>> SCIPSolver::getSeparationOracle<double>(
    const Constraint<double>& face)
  {
    return std::make_shared<SCIPSeparationOracle<double>>(shared_from_this(), face);
  }

#if defined(IPO_RATIONAL_LP)

  /**
   * \brief Returns a rational arithmetic separation oracle for the \p face.
   **/

  template <>
  IPO_EXPORT
  std::shared_ptr<SCIPSeparationOracle<rational>> SCIPSolver::getSeparationOracle<rational>(
    const Constraint<rational>& face)
  {
    return std::make_shared<SCIPSeparationOracle<rational>>(shared_from_this(), face);
  }

#endif /* IPO_RATIONAL_LP */

  static
  void extractPoints(SCIP* scip, std::vector<SCIP_VAR*>& variables, const double* objectiveVector,
    OptimizationOracle<double>::Response& response)
  {
    std::size_t numSolutions = SCIPgetNSols(scip);
    SCIP_SOL** solutions = SCIPgetSols(scip);
    for (std::size_t solIndex = 0; solIndex < numSolutions; ++solIndex)
    {
      SCIP_SOL* sol = solutions[solIndex];
      auto vector = std::make_shared<sparse_vector<double>>();
      bool isFinite = true;
#if defined(IPO_DEBUG)
      double maxAbsValue = 0.0;
#endif /* IPO_DEBUG */
      for (std::size_t i = 0; i < variables.size(); ++i)
      {
        double x = SCIPgetSolVal(scip, sol, variables[i]);
        if (fabs(x) > 1.0e15)
        {
#if defined(IPO_DEBUG)
          std::cout << "A solution had an extremely large entry: " << SCIPvarGetName(variables[i]) << "=" << x
            << "." << std::endl;
#endif /* IPO_DEBUG */
          isFinite = false;
          break;
        }
#if defined(IPO_DEBUG)
        maxAbsValue = std::max(maxAbsValue, fabs(x));
#endif /* IPO_DEBUG */
        if (!SCIPisZero(scip, x))
          vector->push_back(i, x);
      }
      if (isFinite)
      {
#if defined(IPO_DEBUG)
        std::cout << "Extracted a solution with objective value " << (objectiveVector * *vector) << " and "
          "absolute maximum entry " << maxAbsValue << "." << std::endl;
#endif /* IPO_DEBUG */
        response.points.push_back(OptimizationOracle<double>::Response::Point(vector, objectiveVector * *vector));
      }
    }
  }

  static
  void extractRays(SCIP* scip, std::vector<SCIP_VAR*>& variables, OptimizationOracle<double>::Response& response)
  {
    if (SCIPhasPrimalRay(scip))
    {
      auto vector = std::make_shared<sparse_vector<double>>();
      for (std::size_t i = 0; i < variables.size(); ++i)
      {
        double y = SCIPgetPrimalRayVal(scip, variables[i]);
        if (!SCIPisZero(scip, y))
          vector->push_back(i, y);
      }
      response.rays.push_back(OptimizationOracle<double>::Response::Ray(vector));
    }
    else
    {
      throw std::runtime_error("SCIPOptimizationOracle<double>: SCIP reports unboundedness without ray!");
    }
  }

  OptimizationOracle<double>::Response SCIPOptimizationOracle<double>::maximize(const double* objectiveVector,
    const OptimizationOracle<double>::Query& query)
  {
#if defined(IPO_DEBUG)
      std::cout << "SCIPOptimizationOracle<double>::maximize() called." << std::endl;
#endif // IPO_DEBUG

    OptimizationOracle<double>::Response response;

    std::size_t n = space()->dimension();

#if defined(IPO_DEBUG)
    std::cout << "Maximizing objective";
    for (std::size_t v = 0; v < n; ++v)
    {
      if (objectiveVector[v])
        std::cout << " + " << objectiveVector[v] << "*" << space()->variable(v);
    }
    std::cout << ".\n";
    std::cout << "Setting SCIP face to " << _face.vector() << " with rhs " << _face.rhs() << std::endl;
    std::cout << "Setting Time limit to " << query.timeLimit << "." << std::endl;
#endif /* IPO_DEBUG */
    _solver->selectFace(&_face);

    // SCIP settings.
    int oldPresolvingMaxRounds;
    double remainingTime = (query.timeLimit == std::numeric_limits<double>::infinity())
      ? (SCIPinfinity(_solver->_scip)) : query.timeLimit;
    SCIP_CALL_EXC( SCIPgetIntParam(_solver->_scip, "presolving/maxrounds", &oldPresolvingMaxRounds) );
    SCIP_CALL_EXC( SCIPsetLongintParam(_solver->_scip, "limits/totalnodes", -1L) );
    SCIP_CALL_EXC( SCIPsetRealParam(_solver->_scip, "limits/time", remainingTime) );
    SCIP_CALL_EXC( SCIPsetIntParam(_solver->_scip, "limits/solutions", query.maxNumSolutions) );

    // Scale the objective vector down before passing it to SCIP.
    double objectiveScalingFactor = 1.0 / std::max(1.0, maxAbsoluteValue(objectiveVector, n) / 1.0e5 );
#if defined(IPO_DEBUG)
    std::cout << "Scaling objective vector by " << objectiveScalingFactor << "." << std::endl;
#endif /* IPO_DEBUG */
    for (std::size_t i = 0; i < n; ++i)
    {
      SCIP_CALL_EXC( SCIPchgVarObj(_solver->_scip, _solver->_variables[i],
        objectiveVector[i] * objectiveScalingFactor) );
    }

    // Set bound limits of SCIP.
#if defined(IPO_DEBUG)
    if (query.hasMinPrimalBound())
      std::cout << "minPrimalBound = " << query.minPrimalBound() << std::endl;
    if (query.hasMaxDualBound())
      std::cout << "maxDualBound = " << query.maxDualBound() << std::endl;
#endif /* IPO_DEBUG */

    _solver->_boundLimits.minPrimalBound = query.hasMinPrimalBound()
      ? (query.minPrimalBound() * objectiveScalingFactor + 1.0e-6) : -std::numeric_limits<double>::infinity();
    _solver->_boundLimits.maxDualBound = query.hasMaxDualBound()
      ? query.maxDualBound() * objectiveScalingFactor : std::numeric_limits<double>::infinity();

    // First call to SCIP with presolving.
#if defined(IPO_DEBUG)
    std::cout << "Calling SCIP with presolving";
    double timeLimit;
    SCIP_CALL_EXC( SCIPgetRealParam(_solver->_scip, "limits/time", &timeLimit) );
    if (SCIPisFinite(timeLimit))
      std::cout << " using a time limit of " << timeLimit << "s";
    std::cout << "." << std::endl;
    SCIP_CALL_EXC( SCIPwriteOrigProblem(_solver->_scip, "SCIPOptimizationOracle-with-presolve.lp", 0, FALSE) );
    SCIP_CALL_EXC( SCIPwriteParams(_solver->_scip, "SCIPOptimizationOracle-with-presolve.set", FALSE, TRUE) );
#endif /* IPO_DEBUG */

    remainingTime += SCIPgetTotalTime(_solver->_scip);
    SCIP_RETCODE retcode = SCIPsolve(_solver->_scip);
    if (retcode != SCIP_OKAY)
    {
      std::stringstream message;
      message << "SCIPOptimizationOracle<double>: received return code " << retcode << " from SCIPsolve() call.";
      throw std::runtime_error(message.str());
    }
    remainingTime -= SCIPgetTotalTime(_solver->_scip);
    SCIP_STATUS status = SCIPgetStatus(_solver->_scip);

#if defined(IPO_DEBUG)
    std::cout << "SCIP with presolving returned with status " << status << ", primal bound "
      << SCIPgetPrimalbound(_solver->_scip) << " and dual bound " << SCIPgetDualbound(_solver->_scip)
      << "." << std::endl;
#endif /* IPO_DEBUG */

    // Case distinction depending on status.
    if (status == SCIP_STATUS_OPTIMAL || status == SCIP_STATUS_TIMELIMIT)
    {
      assert(!SCIPhasPrimalRay(_solver->_scip));
      assert(SCIPgetNSols(_solver->_scip) > 0);
      extractPoints(_solver->_scip, _solver->_variables, objectiveVector, response);
      response.setPrimalBound(SCIPgetPrimalbound(_solver->_scip) / objectiveScalingFactor);
      response.dualBound = SCIPgetDualbound(_solver->_scip) / objectiveScalingFactor;
      response.hasDualBound = true;
      response.outcome = OptimizationOutcome::FEASIBLE;
      response.hitTimeLimit = status == SCIP_STATUS_TIMELIMIT;
      if (response.points.empty())
      {
        SCIP_CALL_EXC( SCIPwriteOrigProblem(_solver->_scip, "SCIPOptimizationOracle-all-sols-infinite.lp", 0,
          FALSE) );
        throw std::runtime_error("SCIPOptimizationOracle<double>: All SCIP solutions had an infinite entry."
          " Wrote instance to file <SCIPOptimizationOracle-all-sols-infinite.lp>.");
      }
    }
    else if (status == SCIP_STATUS_INFEASIBLE)
    {
      response.outcome = OptimizationOutcome::INFEASIBLE;
    }
    else
    {
      if (status == SCIP_STATUS_INFORUNBD)
      {
        // For the second call we disable presolving.
        SCIP_CALL_EXC( SCIPfreeSolve(_solver->_scip, true) );
        SCIP_CALL_EXC( SCIPfreeTransform(_solver->_scip) );
        SCIP_CALL_EXC( SCIPsetIntParam(_solver->_scip, "presolving/maxrounds", 0) );
        SCIP_CALL_EXC( SCIPsetRealParam(_solver->_scip, "limits/time", remainingTime) );

        // Second call to SCIP, without presolving.
#if defined(IPO_DEBUG)
        std::cout << "Calling SCIP without presolving.";
        double timeLimit;
        SCIP_CALL_EXC( SCIPgetRealParam(_solver->_scip, "limits/time", &timeLimit) );
        if (SCIPisFinite(timeLimit))
          std::cout << " with time limit " << timeLimit << "s";
        std::cout << "." << std::endl;
#endif /* IPO_DEBUG */

        remainingTime += SCIPgetTotalTime(_solver->_scip);
        SCIP_RETCODE retcode = SCIPsolve(_solver->_scip);
        if (retcode != SCIP_OKAY)
        {
          std::stringstream message;
          message << "SCIPOptimizationOracle<double>: received return code " << retcode << " from SCIPsolve() call.";
          throw std::runtime_error(message.str());
        }
        remainingTime -= SCIPgetTotalTime(_solver->_scip);
        status = SCIPgetStatus(_solver->_scip);

#if defined(IPO_DEBUG)
        std::cout << "SCIP without presolving returned with status " << status << ", primal bound "
          << SCIPgetPrimalbound(_solver->_scip) << " and dual bound " << SCIPgetDualbound(_solver->_scip)
          << "." << std::endl;
#endif /* IPO_DEBUG */
      }

      // Case distinction of status if first call was neither optimal, time limit nor infeasible.
      if (status == SCIP_STATUS_INFEASIBLE)
      {
        response.outcome = OptimizationOutcome::INFEASIBLE;
      }
      else if (status == SCIP_STATUS_UNBOUNDED)
      {
        assert(SCIPhasPrimalRay(_solver->_scip));
        extractPoints(_solver->_scip, _solver->_variables, objectiveVector, response);
        extractRays(_solver->_scip, _solver->_variables, response);
        response.outcome = OptimizationOutcome::UNBOUNDED;

        // We're unbounded but but don't have a point, so we solve a feasibility problem.
        if (response.points.empty())
        {
          SCIP_CALL_EXC( SCIPfreeSolve(_solver->_scip, true) );
          SCIP_CALL_EXC( SCIPfreeTransform(_solver->_scip) );
          SCIP_CALL_EXC( SCIPsetIntParam(_solver->_scip, "presolving/maxrounds", oldPresolvingMaxRounds) );
          SCIP_CALL_EXC( SCIPsetRealParam(_solver->_scip, "limits/time", remainingTime) );

          // Update objective vector.
          for (std::size_t i = 0; i < n; ++i)
            SCIP_CALL_EXC( SCIPchgVarObj(_solver->_scip, _solver->_variables[i], 0.0) );

          // Remove limits on primal and dual bounds.
          _solver->_boundLimits.maxDualBound = std::numeric_limits<double>::infinity();
          _solver->_boundLimits.minPrimalBound = -std::numeric_limits<double>::infinity();

          // Third call to SCIP, for feasibility.
#if defined(IPO_DEBUG)
          std::cout << "Calling SCIP for feasibility.";
          double timeLimit;
          SCIP_CALL_EXC( SCIPgetRealParam(_solver->_scip, "limits/time", &timeLimit) );
          if (SCIPisFinite(timeLimit))
            std::cout << " with time limit " << timeLimit << "s";
          std::cout << "." << std::endl;
#endif /* IPO_DEBUG */

          remainingTime += SCIPgetTotalTime(_solver->_scip);
          SCIP_RETCODE retcode = SCIPsolve(_solver->_scip);
          if (retcode != SCIP_OKAY)
          {
            std::stringstream message;
            message << "SCIPOptimizationOracle<double>: received return code " << retcode << " from SCIPsolve() call.";
            throw std::runtime_error(message.str());
          }
          remainingTime -= SCIPgetTotalTime(_solver->_scip);
          status = SCIPgetStatus(_solver->_scip);

#if defined(IPO_DEBUG)
          std::cout << "SCIP for feasibility returned with status " << status << ", primal bound "
            << SCIPgetPrimalbound(_solver->_scip) << " and dual bound " << SCIPgetDualbound(_solver->_scip)
            << "." << std::endl;
#endif /* IPO_DEBUG */

          if (status == SCIP_STATUS_INFEASIBLE)
          {
            response.rays.clear();
            response.outcome = OptimizationOutcome::INFEASIBLE;
          }
          else if (status == SCIP_STATUS_OPTIMAL || status == SCIP_STATUS_TIMELIMIT)
          {
            extractPoints(_solver->_scip, _solver->_variables, objectiveVector, response);
            response.hitTimeLimit = status == SCIP_STATUS_TIMELIMIT;
            if (response.points.empty())
            {
              std::stringstream message;
              message << "SCIPOptimizationOracle<double>: "
                "Feasbility SCIP has found points all of which have extremely large entries.";
              throw std::runtime_error(message.str());
            }
          }
          else
          {
            std::stringstream message;
            message << "SCIPOptimizationOracle<double>: Unhandled SCIP status " << status << " in feasibility call.";
            throw std::runtime_error(message.str());
          }
        }
      }
      else
      {
        std::stringstream message;
        message << "SCIPOptimizationOracle<double>: Unhandled SCIP status " << status;
        if (status == SCIP_STATUS_OPTIMAL)
          message << " (optimal) after a previous status " << SCIP_STATUS_INFORUNBD << " (infeasible/unbounded).";
        else if (status == SCIP_STATUS_INFORUNBD)
          message << " (infeasible/unbounded) even after disabling presolving.";
        else
          message << ".";
        throw std::runtime_error(message.str());
      }
    }

    SCIP_CALL_EXC( SCIPfreeSolve(_solver->_scip, true) );
    SCIP_CALL_EXC( SCIPfreeTransform(_solver->_scip) );
    SCIP_CALL_EXC( SCIPsetIntParam(_solver->_scip, "presolving/maxrounds", oldPresolvingMaxRounds) );

#if !defined(NDEBUG)
    for (auto& point : response.points)
    {
      assert(fabs(*point.vector * objectiveVector - point.objectiveValue) < 1.0e-9);
    }
#endif /* !NDEBUG */

#if defined(IPO_DEBUG)
    std::cout << "Sorting " << response.points.size() << " points." << std::endl;
#endif /* IPO_DEBUG */

    std::sort(response.points.begin(), response.points.end());

#if defined(IPO_DEBUG)
    std::cout << "Returning response " << response << " with primal bound ";
    if (response.hasPrimalBound())
      std::cout << response.primalBound() << std::endl;
    else
      std::cout << "+/-inf" << std::endl;
#endif /* IPO_DEBUG */

    return response;
  }

  
  template <>
  SCIPSeparationOracle<double>::SCIPSeparationOracle(std::shared_ptr<SCIPSolver> solver,
    const Constraint<double>& face)
    : SeparationOracle<double>(solver->name()), _solver(solver), _face(face),
    _approximateFace(face)
  {
    _space = solver->space();
    _solver->addFace(&_approximateFace);
  }

  template <>
  SCIPSeparationOracle<double>::~SCIPSeparationOracle()
  {
    _solver->deleteFace(&_approximateFace);
  }

  template <>
  SeparationResponse<double> SCIPSeparationOracle<double>::getInitial(
    const SeparationQuery& query)
  {
    SeparationResponse<double> result;

    if (query.maxNumInequalities > 0 && !_face.isAlwaysSatisfied())
      result.constraints.push_back(_face);

    struct Visitor
    {
      std::shared_ptr<SCIPSolver> solver;
      const SeparationQuery& query;
      SeparationResponse<double>& response;
      std::size_t iteration;
      std::chrono::time_point<std::chrono::system_clock> started;

      void operator()(const Constraint<double>& constraint)
      {
        // Did we reach a limit?
        if (response.hitTimeLimit || response.constraints.size() == query.maxNumInequalities)
          return;

        iteration = (iteration + 1) % SCIP_SEPARATION_CHECK_TIMELIMIT_FREQUENCY;
        if (iteration == 0)
        {
          std::chrono::duration<double> duration = std::chrono::system_clock::now() - started;
          if (duration.count() > query.timeLimit)
          {
            response.hitTimeLimit = true;
            return;
          }
        }

        response.constraints.push_back(constraint);
      }
    };

    Visitor visitor = { _solver, query, result, 0, std::chrono::system_clock::now() };

    SCIPiterateBounds(_solver->_scip, _solver->_variablesToCoordinates, visitor, true);
    SCIPiterateRows(_solver->_scip, _solver->_variablesToCoordinates, visitor, true);

    return result;
  }

  template <>
  SeparationResponse<double> SCIPSeparationOracle<double>::separate(const double* vector,
    bool isPoint, const SeparationQuery& query)
  {
    SeparationResponse<double> result;

    _solver->selectFace(&_approximateFace);

    struct Visitor
    {
      std::shared_ptr<SCIPSolver> solver;
      const double* vector;
      bool isPoint;
      const SeparationQuery& query;
      SeparationResponse<double>& response;
      std::size_t iteration;
      std::chrono::time_point<std::chrono::system_clock> started;

      void operator()(const Constraint<double>& constraint)
      {
        // Did we reach a limit?
        if (response.hitTimeLimit || response.constraints.size() == query.maxNumInequalities)
          return;

        iteration = (iteration + 1) % SCIP_SEPARATION_CHECK_TIMELIMIT_FREQUENCY;
        if (iteration == 0)
        {
          std::chrono::duration<double> duration = std::chrono::system_clock::now() - started;
          if (duration.count() > query.timeLimit)
          {
            response.hitTimeLimit = true;
            return;
          }
        }

        // Set lhs/rhs to 0 if we are separating a ray.
        double lhs, rhs;
        if (!constraint.hasLhs())
          lhs = -std::numeric_limits<double>::infinity();
        else if (isPoint)
          lhs = constraint.lhs();
        else
          lhs = 0.0;
        if (!constraint.hasRhs())
          rhs = std::numeric_limits<double>::infinity();
        else if (isPoint)
          rhs = constraint.rhs();
        else
          rhs = 0.0;

        double activity = 0.0;
        for (const auto& iter : constraint.vector())
          activity += vector[iter.first] * iter.second;

        if (!SCIPisFeasLE(solver->_scip, activity, rhs)
          || !SCIPisFeasGE(solver->_scip, activity, lhs))
        {
          response.constraints.push_back(constraint);
        }
      }
    };

    Visitor visitor = { _solver, vector, isPoint, query, result, 0,
        std::chrono::system_clock::now() };

    SCIPiterateBounds(_solver->_scip, _solver->_variablesToCoordinates, visitor, false);
    SCIPiterateRows(_solver->_scip, _solver->_variablesToCoordinates, visitor, false);

    return result;
  }

  template <>
  SeparationResponse<double> SCIPSeparationOracle<double>::separateDouble(const double* vector,
    bool isPoint, const SeparationQuery& query)
  {
    return separate(vector, isPoint, query);
  }

#if defined(IPO_RATIONAL_LP)

  template <>
  SCIPSeparationOracle<rational>::SCIPSeparationOracle(std::shared_ptr<SCIPSolver> solver,
    const Constraint<rational>& face)
    : SeparationOracle<rational>(solver->name()), _solver(solver), _face(face),
    _approximateFace(convertConstraint<double>(face))
  {
    _space = solver->space();
    _solver->addFace(&_approximateFace);
  }

  template <>
  SCIPSeparationOracle<rational>::~SCIPSeparationOracle()
  {
    _solver->deleteFace(&_approximateFace);
  }

  template <>
  SeparationResponse<rational> SCIPSeparationOracle<rational>::getInitial(
    const SeparationQuery& query)
  {
    SeparationResponse<rational> result;

    if (query.maxNumInequalities > 0 && !_face.isAlwaysSatisfied())
      result.constraints.push_back(_face);

    struct Visitor
    {
      std::shared_ptr<SCIPSolver> solver;
      const SeparationQuery& query;
      SeparationResponse<rational>& response;
      std::size_t iteration;
      std::chrono::time_point<std::chrono::system_clock> started;

      void operator()(const Constraint<double>& constraint)
      {
        // Did we reach a limit?
        if (response.hitTimeLimit || response.constraints.size() == query.maxNumInequalities)
          return;

        iteration = (iteration + 1) % SCIP_SEPARATION_CHECK_TIMELIMIT_FREQUENCY;
        if (iteration == 0)
        {
          std::chrono::duration<double> duration = std::chrono::system_clock::now() - started;
          if (duration.count() > query.timeLimit)
          {
            response.hitTimeLimit = true;
            return;
          }
        }

        response.constraints.push_back(convertConstraint<rational>(constraint));
      }
    };

    Visitor visitor = { _solver, query, result, 0, std::chrono::system_clock::now() };

    SCIPiterateBounds(_solver->_scip, _solver->_variablesToCoordinates, visitor, true);
    SCIPiterateRows(_solver->_scip, _solver->_variablesToCoordinates, visitor, true);

    return result;
  }

  template <>
  SeparationResponse<rational> SCIPSeparationOracle<rational>::separate(const rational* vector,
    bool isPoint, const SeparationQuery& query)
  {
    SeparationResponse<rational> result;

    _solver->selectFace(&_approximateFace);

    struct Visitor
    {
      std::shared_ptr<SCIPSolver> solver;
      const rational* vector;
      bool isPoint;
      const SeparationQuery& query;
      SeparationResponse<rational>& response;
      std::size_t iteration;
      std::chrono::time_point<std::chrono::system_clock> started;

      void operator()(const Constraint<double>& constraint)
      {
        // Did we reach a limit?
        if (response.hitTimeLimit || response.constraints.size() == query.maxNumInequalities)
          return;

        iteration = (iteration + 1) % SCIP_SEPARATION_CHECK_TIMELIMIT_FREQUENCY;
        if (iteration == 0)
        {
          std::chrono::duration<double> duration = std::chrono::system_clock::now() - started;
          if (duration.count() > query.timeLimit)
          {
            response.hitTimeLimit = true;
            return;
          }
        }

        // Set lhs/rhs to 0 if we are separating a ray.
        double lhs, rhs;
        if (!constraint.hasLhs())
          lhs = -std::numeric_limits<double>::infinity();
        else if (isPoint)
          lhs = constraint.lhs();
        else
          lhs = 0.0;
        if (!constraint.hasRhs())
          rhs = std::numeric_limits<double>::infinity();
        else if (isPoint)
          rhs = constraint.rhs();
        else
          rhs = 0.0;

        double activity = 0.0;
        for (const auto& iter : constraint.vector())
          activity += convertNumber<double>(vector[iter.first]) * iter.second;

        if (!SCIPisFeasLE(solver->_scip, activity, rhs)
          || !SCIPisFeasGE(solver->_scip, activity, lhs))
        {
          response.constraints.push_back(convertConstraint<rational>(constraint));
        }
      }
    };

    Visitor visitor = { _solver, vector, isPoint, query, result, 0,
        std::chrono::system_clock::now() };

    SCIPiterateBounds(_solver->_scip, _solver->_variablesToCoordinates, visitor, false);
    SCIPiterateRows(_solver->_scip, _solver->_variablesToCoordinates, visitor, false);

    return result;
  }

  template <>
  SeparationResponse<rational> SCIPSeparationOracle<rational>::separateDouble(const double* vector,
    bool isPoint, const SeparationQuery& query)
  {
    SeparationResponse<rational> result;

    _solver->selectFace(&_approximateFace);

    struct Visitor
    {
      std::shared_ptr<SCIPSolver> solver;
      const double* vector;
      bool isPoint;
      const SeparationQuery& query;
      SeparationResponse<rational>& response;
      std::size_t iteration;
      std::chrono::time_point<std::chrono::system_clock> started;

      void operator()(const Constraint<double>& constraint)
      {
        // Did we reach a limit?
        if (response.hitTimeLimit || response.constraints.size() == query.maxNumInequalities)
          return;

        iteration = (iteration + 1) % SCIP_SEPARATION_CHECK_TIMELIMIT_FREQUENCY;
        if (iteration == 0)
        {
          std::chrono::duration<double> duration = std::chrono::system_clock::now() - started;
          if (duration.count() > query.timeLimit)
          {
            response.hitTimeLimit = true;
            return;
          }
        }

        // Set lhs/rhs to 0 if we are separating a ray.
        double lhs, rhs;
        if (!constraint.hasLhs())
          lhs = -std::numeric_limits<double>::infinity();
        else if (isPoint)
          lhs = constraint.lhs();
        else
          lhs = 0.0;
        if (!constraint.hasRhs())
          rhs = std::numeric_limits<double>::infinity();
        else if (isPoint)
          rhs = constraint.rhs();
        else
          rhs = 0.0;

        double activity = 0.0;
        for (const auto& iter : constraint.vector())
          activity += convertNumber<double>(vector[iter.first]) * iter.second;

        if (!SCIPisFeasLE(solver->_scip, activity, rhs)
          || !SCIPisFeasGE(solver->_scip, activity, lhs))
        {
          response.constraints.push_back(convertConstraint<rational>(constraint));
        }
      }
    };

    Visitor visitor = { _solver, vector, isPoint, query, result, 0,
        std::chrono::system_clock::now() };

    SCIPiterateBounds(_solver->_scip, _solver->_variablesToCoordinates, visitor, false);
    SCIPiterateRows(_solver->_scip, _solver->_variablesToCoordinates, visitor, false);

    return result;
  }

#endif /* IPO_RATIONAL_LP */
  
} /* namespace ipo */
