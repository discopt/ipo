// #define IPO_DEBUG // Uncomment to debug this file.

#include <ipo/oracles_scip.hpp>

// TODO: Remove objective limit and add event handler for primal and dual bound.

#include <cassert>
#include <functional>
#include <chrono>
#include <sstream>
#include <iostream>

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
      && SCIPgetPrimalbound(scip) >= eventhdlrdata->minPrimalBound)
      || (eventhdlrdata->maxDualBound < std::numeric_limits<double>::infinity()
      && SCIPgetDualbound(scip) <= eventhdlrdata->maxDualBound))
    {
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

#if defined(IPO_WITH_GMP) && defined(IPO_WITH_SOPLEX)
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
#endif /* IPO_WITH_GMP && IPO_WITH_SOPLEX */
    
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

#if defined(IPO_WITH_GMP) && defined(IPO_WITH_SOPLEX)
    delete _extender;
#endif /* IPO_WITH_GMP && IPO_WITH_SOPLEX */

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

  SCIPRealOptimizationOracle::SCIPRealOptimizationOracle(std::shared_ptr<SCIPSolver> solver,
    const Constraint<double>& face)
    : OptimizationOracle<double>(solver->name()), _solver(solver), _face(face)
  {
    _space = solver->space();
    _solver->addFace(&_face);
  }

  SCIPRealOptimizationOracle::~SCIPRealOptimizationOracle()
  {
    _solver->deleteFace(&_face);
  }

  OptimizationOracle<double>::Response SCIPRealOptimizationOracle::maximize(const double* objectiveVector,
      const OptimizationOracle<double>::Query& query)
  {
    OptimizationOracle<double>::Response response;

#if defined(IPO_DEBUG)
    std::cout << "Setting SCIP face to " << _face.vector() << " with rhs " << _face.rhs() << std::endl;
    std::cout << "Setting Time limit to " << query.timeLimit << "." << std::endl;
#endif /* IPO_DEBUG */
    _solver->selectFace(&_face);

    SCIP_CALL_EXC( SCIPsetLongintParam(_solver->_scip, "limits/totalnodes", -1L) );
    SCIP_CALL_EXC( SCIPsetRealParam(_solver->_scip, "limits/time",
      query.timeLimit == std::numeric_limits<double>::infinity() ? SCIPinfinity(_solver->_scip)
      : query.timeLimit) );

    SCIP_CALL_EXC( SCIPsetIntParam(_solver->_scip, "limits/solutions", query.maxNumSolutions) );

    std::size_t n = space()->dimension();

    // Compute factor for scaling the objective vector down.

    double maxEntry = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
      double ci = fabs(objectiveVector[i]);
      if (ci > maxEntry)
        maxEntry = ci;
    }

    const double maxAllowedEntry = 10e5;
    double scalingFactor = 1.0 / std::max(1.0, maxEntry / maxAllowedEntry);

    for (std::size_t i = 0; i < n; ++i)
    {
      SCIP_CALL_EXC( SCIPchgVarObj(_solver->_scip, _solver->_variables[i],
        objectiveVector[i] * scalingFactor) );
    }

    // Bound limits.

#if defined(IPO_DEBUG)
    if (query.hasMinPrimalBound())
      std::cout << "minPrimalBound = " << query.minPrimalBound() << std::endl;
    if (query.hasMaxDualBound())
      std::cout << "maxDualBound = " << query.maxDualBound() << std::endl;
#endif /* IPO_DEBUG */

    _solver->_boundLimits.minPrimalBound = query.hasMinPrimalBound()
      ? query.minPrimalBound() * scalingFactor : -std::numeric_limits<double>::infinity();
    _solver->_boundLimits.maxDualBound = query.hasMaxDualBound()
      ? query.maxDualBound() * scalingFactor : std::numeric_limits<double>::infinity();

    // Save presolve settings.

    int oldMaxRounds;
    SCIP_CALL_EXC( SCIPgetIntParam(_solver->_scip, "presolving/maxrounds", &oldMaxRounds) );

    for (int attempt = 1; attempt <= 2; ++attempt)
    {
#if defined(IPO_DEBUG)
      std::cout << "Solving optimization problem";
      double timeLimit;
      SCIP_CALL_EXC( SCIPgetRealParam(_solver->_scip, "limits/time", &timeLimit) );
      if (SCIPisFinite(timeLimit))
        std::cout << " with time limit " << timeLimit << "s";
      std::cout << "." << std::endl;
#endif /* IPO_DEBUG */

      double solutionTime = -SCIPgetTotalTime(_solver->_scip);
      SCIP_RETCODE retcode = SCIPsolve(_solver->_scip);
      if (retcode != SCIP_OKAY)
      {
        std::cerr << "SCIPRealOptimizationOracle received return code " << retcode
          << " from SCIPsolve() call." << std::endl;
      }
      solutionTime += SCIPgetTotalTime(_solver->_scip);

#if defined(IPO_DEBUG)
      std::cout << "SCIP returned with return code " << retcode << ", status "
        << SCIPgetStatus(_solver->_scip) << ", primal bound " << SCIPgetPrimalbound(_solver->_scip)
        << " and dual bound " << SCIPgetDualbound(_solver->_scip) << "." << std::endl;
#endif /* IPO_DEBUG */

      bool hasRay = SCIPhasPrimalRay(_solver->_scip);
      if (hasRay)
      {
        auto vector = std::make_shared<sparse_vector<double>>();
        for (std::size_t i = 0; i < n; ++i)
        {
          double y = SCIPgetPrimalRayVal(_solver->_scip, _solver->_variables[i]);
          if (!SCIPisZero(_solver->_scip, y))
            vector->push_back(i, y);
        }
        response.rays.push_back(OptimizationOracle<double>::Response::Ray(vector));
        response.outcome = OptimizationOutcome::UNBOUNDED;
        break;
      }

      SCIP_STATUS status = SCIPgetStatus(_solver->_scip);
      if (status == SCIP_STATUS_UNBOUNDED && attempt == 2 && !hasRay)
      {
        throw std::runtime_error("SCIPOptimizationOracle: SCIP reports unboundedness without ray!");
      }

      std::size_t numSolutions = SCIPgetNSols(_solver->_scip);
      if (status != SCIP_STATUS_INFORUNBD && status != SCIP_STATUS_UNBOUNDED && numSolutions > 0)
      {
        // Note that SCIP_STATUS_INFEASIBLE is fine since this is issued if the objective limit was
        // reached. However, if solutions were found, we want to process them.

        SCIP_SOL** solutions = SCIPgetSols(_solver->_scip);
        for (std::size_t solIndex = 0; solIndex < numSolutions; ++solIndex)
        {
          SCIP_SOL* sol = solutions[solIndex];
          auto vector = std::make_shared<sparse_vector<double>>();
          for (std::size_t i = 0; i < n; ++i)
          {
            double x = SCIPgetSolVal(_solver->_scip, sol, _solver->_variables[i]);
            if (!SCIPisZero(_solver->_scip, x))
              vector->push_back(i, x);
          }
          response.points.push_back(OptimizationOracle<double>::Response::Point(vector,
            objectiveVector * *vector));
        }
        if (status == SCIP_STATUS_TIMELIMIT)
          response.hitTimeLimit = true;
        response.primalBound = SCIPgetPrimalbound(_solver->_scip) / scalingFactor;
        response.dualBound = SCIPgetDualbound(_solver->_scip) / scalingFactor;
        response.hasDualBound = true;
        response.outcome = OptimizationOutcome::FEASIBLE;
#if defined(IPO_DEBUG)
        std::cout << "Created response." << std::endl;
#endif /* IPO_DEBUG */
        break;
      }
      else if (status == SCIP_STATUS_INFEASIBLE)
      {
        assert(numSolutions == 0);
        response.outcome = OptimizationOutcome::INFEASIBLE;
        break;
      }
      else if (status == SCIP_STATUS_UNBOUNDED)
      {
        if (SCIPhasPrimalRay(_solver->_scip))
        {
          auto vector = std::make_shared<sparse_vector<double>>();
          for (std::size_t i = 0; i < n; ++i)
          {
            double y = SCIPgetPrimalRayVal(_solver->_scip, _solver->_variables[i]);
            if (!SCIPisZero(_solver->_scip, y))
              vector->push_back(i, y);
          }
          response.rays.push_back(OptimizationOracle<double>::Response::Ray(vector));
          response.outcome = OptimizationOutcome::UNBOUNDED;
        }
      }
      else if (status == SCIP_STATUS_TIMELIMIT)
      {
        response.hitTimeLimit = true;
        response.outcome = OptimizationOutcome::TIMEOUT;
        if (SCIPisFinite(SCIPgetPrimalbound(_solver->_scip)))
          response.primalBound = SCIPgetPrimalbound(_solver->_scip) / scalingFactor;
        if (SCIPisFinite(SCIPgetDualbound(_solver->_scip)))
        {
          response.dualBound = SCIPgetDualbound(_solver->_scip) / scalingFactor;
          response.hasDualBound = true;
        }
      }
      else if (attempt == 2)
      {
        std::ostringstream ss;
        ss << "SCIPOptimizationOracle: unhandled SCIP status code " << status << ".";
        throw std::runtime_error(ss.str());
      }

      // Disable presolving for the second round.

      SCIP_CALL_EXC( SCIPsetIntParam(_solver->_scip, "presolving/maxrounds", 0) );
      SCIP_CALL_EXC( SCIPfreeSolve(_solver->_scip, true) );
      SCIP_CALL_EXC( SCIPfreeTransform(_solver->_scip) );

      // Adapt time limit.

#if defined(IPO_DEBUG)
      std::cout << "Time limit is " << query.timeLimit << "." << std::endl;
      std::cout << "SCIP used " << solutionTime << "s already." << std::endl;
#endif
      if (query.timeLimit < std::numeric_limits<double>::infinity())
      {
        SCIP_CALL_EXC( SCIPsetRealParam(_solver->_scip, "limits/time",
          std::max(0.0, query.timeLimit - solutionTime)) );
      }

#if defined(IPO_DEBUG)
      std::cout << "Disabling presolve for a second optimization attempt." << std::endl;
#endif /* IPO_DEBUG */
    }

    SCIP_CALL_EXC( SCIPsetIntParam(_solver->_scip, "presolving/maxrounds", oldMaxRounds) );
    SCIP_CALL_EXC( SCIPfreeSolve(_solver->_scip, true) );
    SCIP_CALL_EXC( SCIPfreeTransform(_solver->_scip) );

    if (!response.rays.empty() && response.points.empty())
    {
      response.primalBound = std::numeric_limits<double>::infinity();

      // We have to check feasibility.

      for (std::size_t i = 0; i < n; ++i)
        SCIP_CALL_EXC( SCIPchgVarObj(_solver->_scip, _solver->_variables[i], 0.0) );

      // Remove limits on primal and dual bounds.
      _solver->_boundLimits.maxDualBound = std::numeric_limits<double>::infinity();
      _solver->_boundLimits.minPrimalBound = -std::numeric_limits<double>::infinity();

#if defined(IPO_DEBUG)
      std::cout << "Solving feasibility problem." << std::endl;
#endif /* IPO_DEBUG */

      SCIP_CALL_EXC( SCIPsolve(_solver->_scip) );

#if defined(IPO_DEBUG)
      std::cout << "SCIP returned with status " << SCIPgetStatus(_solver->_scip) << "."
        << std::endl;
      SCIPprintTimingStatistics(_solver->_scip, stdout);
      fflush(stdout);
#endif /* IPO_DEBUG */

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
        response.points.push_back(OptimizationOracle<double>::Response::Point(vector, objectiveValue)); 
      }
      else
      {
        response.rays.clear();
        response.primalBound = -std::numeric_limits<double>::infinity();
      }
      SCIP_CALL_EXC( SCIPfreeSolve(_solver->_scip, true) );
      SCIP_CALL_EXC( SCIPfreeTransform(_solver->_scip) );
    }

#if !defined(NDEBUG)
    for (auto& point : response.points)
    {
      assert(fabs(*point.vector * objectiveVector - point.objectiveValue) < 1.0e-9);
    }
#endif /* !NDEBUG */

#if defined(IPO_DEBUG)
    std::cout << "Sorting points." << std::endl;
#endif /* IPO_DEBUG */

    std::sort(response.points.begin(), response.points.end());

    return response;
  }

  SCIPRealSeparationOracle::SCIPRealSeparationOracle(std::shared_ptr<SCIPSolver> solver,
    const Constraint<double>& face)
    : RealSeparationOracle(solver->name()), _solver(solver), _face(face)
  {
    _space = solver->space();
    _solver->addFace(&_face);
  }

  SCIPRealSeparationOracle::~SCIPRealSeparationOracle()
  {
    _solver->deleteFace(&_face);
  }

  RealSeparationOracle::Response SCIPRealSeparationOracle::getInitial(
    const RealSeparationOracle::Query& query)
  {
    RealSeparationOracle::Response result;

    if (query.maxNumInequalities > 0 && !_face.isAlwaysSatisfied())
      result.constraints.push_back(_face);

    struct Visitor
    {
      std::shared_ptr<SCIPSolver> solver;
      const SCIPRealSeparationOracle::Query& query;
      SCIPRealSeparationOracle::Response& response;
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

  RealSeparationOracle::Response SCIPRealSeparationOracle::separate(const double* vector,
    bool isPoint, const RealSeparationOracle::Query& query)
  {
    RealSeparationOracle::Response result;

    _solver->selectFace(&_face);

    struct Visitor
    {
      std::shared_ptr<SCIPSolver> solver;
      const double* vector;
      bool isPoint;
      const SCIPRealSeparationOracle::Query& query;
      SCIPRealSeparationOracle::Response& response;
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

} /* namespace ipo */
