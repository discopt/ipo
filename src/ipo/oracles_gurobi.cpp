// #define IPO_DEBUG // Uncomment to debug this file.

#include <ipo/constraint.hpp>
#include <ipo/oracles_gurobi.hpp>

// TODO: Remove objective limit and add event handler for primal and dual bound.

#include <cassert>
#include <functional>
#include <chrono>
#include <sstream>
#include <iostream>
#include <unordered_map>

/**
 * \brief Macro to raise a SCIPException in case of a SCIP error.
 **/

#define GUROBI_CALL_EXC(x) \
{ \
  int _retcode; \
  if ((_retcode = (x)) != 0) \
    throw ipo::GurobiException(_retcode, __FILE__, __LINE__); \
}

static const int GUROBI_SEPARATION_CHECK_TIMELIMIT_FREQUENCY = 1000;

// #define EVENTHDLR_NAME "boundchange"
// #define EVENTHDLR_DESC "event handler for primal or dual bound change"

namespace ipo
{
  /**
   * \brief Exception handling class for SCIP.
   *
   * Represents a SCIP error in C++.
   **/

  class GurobiException: public std::exception
  {
  public:

    /**
     * \brief Constructs a SCIPException from an error code.
     **/

    GurobiException(int retcode, const char* fileName, int line)
      : _retcode(retcode)
    {
      switch (retcode)
      {
      case 0:
        snprintf(_message, 256, "Normal termination");
      break;
      case GRB_ERROR_OUT_OF_MEMORY:
        snprintf(_message, 256, "Available memory was exhausted");
      break;
      case GRB_ERROR_NO_LICENSE:
        snprintf(_message, 256, "Insufficient memory error");
      break;
      case GRB_ERROR_SIZE_LIMIT_EXCEEDED:
        snprintf(_message, 256, "Attempted to solve a model that is larger than the limit for a demo license");
      break;
      case GRB_ERROR_FILE_READ:
        snprintf(_message, 256, "Failed to read the requested file");
      break;
      case GRB_ERROR_FILE_WRITE:
        snprintf(_message, 256, "Failed to write the requested file");
      break;
      case GRB_ERROR_NUMERIC:
        snprintf(_message, 256, "Numerical error during requested operation");
      break;
      case GRB_ERROR_NODEFILE:
        snprintf(_message, 256, "Error in reading or writing a node file during MIP optimization");
      break;
      break;
      default:
        snprintf(_message, 256, "Unknown error code %d", retcode);
      break;
      }
      int len = strlen(_message);
      snprintf(&_message[len], 256 - len, " at %s:%d", fileName, line);
    }

    /**
     * \brief Destructor.
     **/

    ~GurobiException(void) throw ()
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
    /// Gurobi error code.
    int _retcode;
  };

  /**
   * \brief Auxiliary method that iterates over the bounds as constraints.
   *
   * Auxiliary method that iterates over the bound constraints of all variables calling \p visitor
   * with each of them as a constraint.
   *
   * \p ranged Whether to return ranged rows.
   */

  static
  void GRBiterateBounds(GRBmodel* model, std::shared_ptr<Space> space,
    std::function<void(Constraint<double>&& constraint)> visitor, bool ranged)
  {
    sparse_vector<double> entries;
    for (std::size_t v = 0; v < space->dimension(); ++v)
    {
      double lhs, rhs;
      GUROBI_CALL_EXC( GRBgetdblattrelement(model, GRB_DBL_ATTR_LB, v, &lhs) );
      GUROBI_CALL_EXC( GRBgetdblattrelement(model, GRB_DBL_ATTR_UB, v, &rhs) );
      if (fabs(lhs - rhs) < 1.0e-9)
      {
        double average = (lhs + rhs) / 2.0;
        auto vector = std::make_shared<sparse_vector<double>>();
        vector->push_back(v, 1.0);
        visitor(Constraint<double>(average, vector, average, ConstraintType::EQUATION));
        continue;
      }
      if (lhs <= -GRB_INFINITY && rhs >= GRB_INFINITY)
        continue;
      if (lhs <= -GRB_INFINITY || (rhs < GRB_INFINITY && !ranged))
      {
        auto vector = std::make_shared<sparse_vector<double>>();
        vector->push_back(v, 1.0);
        visitor(Constraint<double>(vector, rhs));
      }
      if (rhs >= GRB_INFINITY || (lhs > -GRB_INFINITY && !ranged))
      {
        auto vector = std::make_shared<sparse_vector<double>>();
        vector->push_back(v, 1.0);
        visitor(Constraint<double>(lhs, vector));
      }
      if (ranged && lhs > -GRB_INFINITY && rhs < GRB_INFINITY)
      {
        auto vector = std::make_shared<sparse_vector<double>>();
        vector->push_back(v, 1.0);
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
  void GRBiterateRows(GRBmodel* model, std::shared_ptr<Space> space,
    std::function<void(Constraint<double>&& constraint)> visitor, bool ranged)
  {
    int numConstraints;
    GUROBI_CALL_EXC( GRBgetintattr(model, GRB_INT_ATTR_NUMCONSTRS, &numConstraints) );

    int numNonzeros;
    GUROBI_CALL_EXC( GRBgetconstrs(model, &numNonzeros, nullptr, nullptr, nullptr, 0, numConstraints) );

    std::vector<char> senses(numConstraints);
    std::vector<double> rhs(numConstraints);
    std::vector<int> slices(numConstraints + 1);
    std::vector<int> columns(numNonzeros);
    std::vector<double> coefficients(numNonzeros);

    GUROBI_CALL_EXC( GRBgetcharattrarray(model, GRB_CHAR_ATTR_SENSE, 0, numConstraints, &senses[0]) );
    GUROBI_CALL_EXC( GRBgetdblattrarray(model, GRB_DBL_ATTR_RHS, 0, numConstraints, &rhs[0]) );
    GUROBI_CALL_EXC( GRBgetconstrs(model, &numNonzeros, &slices[0], &columns[0], &coefficients[0], 0, numConstraints) );
    slices.back() = numNonzeros;

    std::vector<std::pair<std::size_t, double>> unsortedEntries;
    sparse_vector<double> entries;
    for (int c = 0; c < numConstraints; ++c)
    {
      std::size_t beyond = slices[c + 1];
      unsortedEntries.clear();
      for (std::size_t i = slices[c]; i < beyond; ++i)
        unsortedEntries.push_back(std::make_pair(static_cast<std::size_t>(columns[i]), coefficients[i]));

      auto vector = std::make_shared<sparse_vector<double>>(std::move(unsortedEntries), true);
      if (senses[c] == '<')
        visitor(Constraint<double>(vector, rhs[c]));
      else if (senses[c] == '>')
        visitor(Constraint<double>(rhs[c], vector));
      else if (senses[c] == '=')
        visitor(Constraint<double>(rhs[c], vector, rhs[c]));
    }
  }

  GurobiSolver::GurobiSolver(GRBmodel*&& model)
    : _env(nullptr), _model(model), _currentFace(nullptr)
  {
    model = nullptr;
#if !defined(IPO_DEBUG)
    GRBsetintparam(GRBgetenv(_model), "OutputFlag", 0);
#endif /* !IPO_DEBUG */
    initialize();
  }

  GurobiSolver::GurobiSolver(const std::string& fileName)
    : _currentFace(nullptr)
  {
    GUROBI_CALL_EXC( GRBemptyenv(&_env) );
#if !defined(IPO_DEBUG)
    GUROBI_CALL_EXC( GRBsetintparam(_env, "OutputFlag", 0) );
#endif /* !IPO_DEBUG */
    GUROBI_CALL_EXC( GRBstartenv(_env) );
    GUROBI_CALL_EXC( GRBreadmodel(_env, fileName.c_str(), &_model) );

    initialize();
  }

  void GurobiSolver::initialize()
  {
    GUROBI_CALL_EXC( GRBupdatemodel(_model) );

    int n;
    GUROBI_CALL_EXC( GRBgetintattr(_model, GRB_INT_ATTR_NUMVARS, &n) );

    _instanceObjective = new double[n+1];
    GUROBI_CALL_EXC( GRBgetdblattrarray(_model, GRB_DBL_ATTR_OBJ, 0, n, &_instanceObjective[0]) );
    GUROBI_CALL_EXC( GRBgetdblattr(_model, GRB_DBL_ATTR_OBJCON, &_instanceObjective[n]) );
    int modelSense;
    GUROBI_CALL_EXC( GRBgetintattr(_model, GRB_INT_ATTR_MODELSENSE, &modelSense) );
    std::vector<std::string> variableNames;
    variableNames.resize(n);
    char* name;
    GRBgetstrattr(_model, GRB_STR_ATTR_MODELNAME, &name);
    _name = name;
    for (int i = 0; i < n; ++i)
    {
      GUROBI_CALL_EXC( GRBgetstrattrelement(_model, GRB_STR_ATTR_VARNAME, i, &name) );
      variableNames[i] = name;
      _instanceObjective[i] *= -modelSense;
    }
    _instanceObjective[n] *= -modelSense;
    _space = std::make_shared<Space>(std::move(variableNames));

    // Make it a maximization problem.
    GUROBI_CALL_EXC( GRBsetintattr(_model, GRB_INT_ATTR_MODELSENSE, -1) );

    // Extract bounds and integrality for initializing the MIP extended.
    std::vector<bool> integrality(n);
    std::vector<std::pair<double, double>> bounds(n);
    for (int i = 0; i < n; ++i)
    {
      char value;
      GUROBI_CALL_EXC( GRBgetcharattrelement(_model, GRB_CHAR_ATTR_VTYPE, i, &value) );
      integrality[i] = (value == 'B' || value == 'I' || value == 'N');
      GUROBI_CALL_EXC( GRBgetdblattrelement(_model, GRB_DBL_ATTR_LB, i, &bounds[i].first) );
      if (bounds[i].first <= -1e20)
        bounds[i].first = -std::numeric_limits<double>::infinity();
      GUROBI_CALL_EXC( GRBgetdblattrelement(_model, GRB_DBL_ATTR_UB, i, &bounds[i].second) );
      if (bounds[i].second >= 1e20)
        bounds[i].second = std::numeric_limits<double>::infinity();
      if (value == 'S' || value == 'N')
      {
        bounds[i].first = std::min(bounds[i].first, 0.0);
        bounds[i].second = std::max(bounds[i].second, 0.0);
      }

#if defined(IPO_DEBUG)
      std::cout << "Variable bounds are " << bounds[i].first << " <= x_" << i << " <= " << bounds[i].second
        << std::endl;
#endif /* IPO_DEBUG */
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
    GRBiterateRows(_model, space(), visitor, true);
#endif /* IPO_RATIONAL_LP */

    // SCIP_EVENTHDLR* eventhdlr = NULL;
    // SCIP_CALL_EXC( SCIPincludeEventhdlrBasic(_scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
      // eventExecGlobalBoundChange, (SCIP_EVENTHDLRDATA*)(&_boundLimits)) );
    // assert(eventhdlr != NULL);

    // SCIP_CALL_EXC( SCIPsetEventhdlrCopy(_scip, eventhdlr, eventCopyGlobalBoundChange) );
    // SCIP_CALL_EXC( SCIPsetEventhdlrInit(_scip, eventhdlr, eventInitGlobalBoundChange) );
    // SCIP_CALL_EXC( SCIPsetEventhdlrExit(_scip, eventhdlr, eventExitGlobalBoundChange) );
  }

  GurobiSolver::~GurobiSolver()
  {
    delete[] _instanceObjective;

#if defined(IPO_RATIONAL_LP)
    delete _extender;
#endif /* IPO_RATIONAL_LP */

    GRBfreemodel(_model);
    if (_env)
      GRBfreeenv(_env);
  }

  void GurobiSolver::addFace(Constraint<double>* face)
  {
    assert(face);
  }

  void GurobiSolver::deleteFace(Constraint<double>* face)
  {
    assert(face);

    if (face->isAlwaysSatisfied())
      return;

    if (face == _currentFace)
      selectFace(nullptr);
  }

  void GurobiSolver::selectFace(Constraint<double>* face)
  {
    if (face && face->isAlwaysSatisfied())
      face = nullptr;

    if (face == _currentFace)
    {
#if defined(IPO_DEBUG)
      std::cout << "Face is already selected." << std::endl;
#endif /* IPO_DEBUG */
      return;
    }

    int numConstraints;
    GUROBI_CALL_EXC( GRBgetintattr(_model, GRB_INT_ATTR_NUMCONSTRS, &numConstraints) );

    // Remove from Gurobi model.
    if (_currentFace)
    {
#if defined(IPO_DEBUG)
      std::cout << "Disabling old face." << std::endl;
#endif /* IPO_DEBUG */


      int faceConstraint = numConstraints - 1;
      GRBdelconstrs(_model, 1, &faceConstraint);
    }

    if (face)
    {
#if defined(IPO_DEBUG)
      printf("Enabling new face %p with lhs %f and rhs %f.\n", face, face->lhs(), face->rhs());
#endif /* IPO_DEBUG */
//       SCIP_CALL_EXC( SCIPaddCons(_scip, _faceConstraints.at(face)) );

      std::vector<int> indices(face->vector().size());
      std::vector<double> coefficients(face->vector().size());
      size_t i = 0;
      for (auto iter : face->vector())
      {
        indices[i] = iter.first;
        coefficients[i] = iter.second;
        ++i;
      }
      double rhs;
      if (face->type() == ConstraintType::LESS_OR_EQUAL)
        rhs = face->rhs();
      else if (face->type() == ConstraintType::GREATER_OR_EQUAL)
        rhs = face->lhs();
      else if (face->type() == ConstraintType::EQUATION)
      {
        rhs = face->rhs();
        assert(abs(face->lhs() - rhs) < 1.0e-9);
      }
      else
        throw std::runtime_error("Cannot use a ranged constraint or equation to define a face.");

      GUROBI_CALL_EXC( GRBaddconstr(_model, face->vector().size(), &indices[0], &coefficients[0], '=', rhs, "face") );
    }

    _currentFace = face;
  }

  /**
   * \brief Returns a double-arithmetic optimization oracle for the requested \p face.
   */

  template <>
  IPO_EXPORT
  std::shared_ptr<GurobiOptimizationOracle<double>> GurobiSolver::getOptimizationOracle<double>(
    const Constraint<double>& face)
  {
    return std::make_shared<GurobiOptimizationOracle<double>>(shared_from_this(), face);
  }

#if defined(IPO_RATIONAL_LP)

  /**
   * \brief Returns a rational optimization oracle for the requested \p face.
   **/

  template <>
  IPO_EXPORT
  std::shared_ptr<GurobiOptimizationOracle<rational>> GurobiSolver::getOptimizationOracle<rational>(
    const Constraint<rational>& face)
  {
    auto approximateFace = convertConstraint<double>(face);
    auto approximateOracle = getOptimizationOracle<double>(approximateFace);
    return std::make_shared<GurobiOptimizationOracle<rational>>(_extender, approximateOracle, face);
  }

#endif /* IPO_RATIONAL_LP */

  GurobiOptimizationOracle<double>::GurobiOptimizationOracle(std::shared_ptr<GurobiSolver> solver,
    const Constraint<double>& face)
    : OptimizationOracle<double>(solver->name()), _solver(solver), _face(face)
  {
    _space = solver->space();
    _solver->addFace(&_face);
  }

  GurobiOptimizationOracle<double>::~GurobiOptimizationOracle()
  {
    _solver->deleteFace(&_face);
  }

  /**
   * \brief Returns a double arithmetic separation oracle for the \p face.
   **/

  template <>
  IPO_EXPORT
  std::shared_ptr<GurobiSeparationOracle<double>> GurobiSolver::getSeparationOracle<double>(
    const Constraint<double>& face)
  {
    return std::make_shared<GurobiSeparationOracle<double>>(shared_from_this(), face);
  }

#if defined(IPO_RATIONAL_LP)

  /**
   * \brief Returns a rational arithmetic separation oracle for the \p face.
   **/

  template <>
  IPO_EXPORT
  std::shared_ptr<GurobiSeparationOracle<rational>> GurobiSolver::getSeparationOracle<rational>(
    const Constraint<rational>& face)
  {
    return std::make_shared<GurobiSeparationOracle<rational>>(shared_from_this(), face);
  }

#endif /* IPO_RATIONAL_LP */

  static
  void extractPoints(GRBmodel* model, std::shared_ptr<Space> space, const double* objectiveVector, OptimizationOracle<double>::Response& response)
  {
    int numSolutions;
    GUROBI_CALL_EXC( GRBgetintattr(model, GRB_INT_ATTR_SOLCOUNT, &numSolutions) );
    std::vector<double> dense(space->dimension());
    for (int solIndex = 0; solIndex < numSolutions; ++solIndex)
    {
      GUROBI_CALL_EXC( GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_SOLUTIONNUMBER, solIndex) );
      GUROBI_CALL_EXC( GRBgetdblattrarray(model, GRB_DBL_ATTR_Xn, 0, space->dimension(), &dense[0]) );
      auto vector = std::make_shared<sparse_vector<double>>();
#if defined(IPO_DEBUG)
      double maxAbsValue = 0.0;
#endif /* IPO_DEBUG */
      bool isFinite = true;
      for (std::size_t i = 0; i < space->dimension(); ++i)
      {
        double x = dense[i];
        if (fabs(x) > 1.0e15)
        {
#if defined(IPO_DEBUG)
          std::cout << "A solution had an extremely large entry: " << space->variable(i) << "=" << x
            << "." << std::endl;
#endif /* IPO_DEBUG */
          isFinite = false;
          break;
        }
#if defined(IPO_DEBUG)
        maxAbsValue = std::max(maxAbsValue, fabs(x));
#endif /* IPO_DEBUG */
        if (fabs(x) > 1.0e-12)
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
  void extractRays(GRBmodel* model, std::shared_ptr<Space> space, OptimizationOracle<double>::Response& response)
  {
    auto vector = std::make_shared<sparse_vector<double>>();
    std::vector<double> dense(space->dimension());
    GUROBI_CALL_EXC( GRBgetdblattrarray(model, GRB_DBL_ATTR_UNBDRAY, 0, space->dimension(), &dense[0]) );

    for (std::size_t i = 0; i < space->dimension(); ++i)
    {
      double y = dense[i];
      if (fabs(y) > 1.0e-12)
        vector->push_back(i, y);
    }
    response.rays.push_back(OptimizationOracle<double>::Response::Ray(vector));
  }

  OptimizationOracle<double>::Response GurobiOptimizationOracle<double>::maximize(const double* objectiveVector,
    const OptimizationOracle<double>::Query& query)
  {
#if defined(IPO_DEBUG)
      std::cout << "GurobiOptimizationOracle<double>::maximize() called." << std::endl;
#endif // IPO_DEBUG

    OptimizationOracle<double>::Response response;

    std::size_t n = space()->dimension();
    double remainingTime = query.timeLimit;

#if defined(IPO_DEBUG)
    std::cout << "Maximizing objective";
    for (std::size_t v = 0; v < n; ++v)
    {
      if (objectiveVector[v])
        std::cout << " + " << objectiveVector[v] << "*" << space()->variable(v);
    }
    std::cout << ".\n";
    std::cout << "Setting Gurobi face to " << _face.vector() << " with rhs " << _face.rhs() << std::endl;
    std::cout << "Setting time limit to " << remainingTime << "." << std::endl;
#endif /* IPO_DEBUG */

    _solver->selectFace(&_face);

    // Gurobi settings.
    int oldDualReductions;
    GUROBI_CALL_EXC( GRBgetintparam(GRBgetenv(_solver->_model), GRB_INT_PAR_DUALREDUCTIONS, &oldDualReductions) );
    GUROBI_CALL_EXC( GRBsetdblparam(GRBgetenv(_solver->_model), GRB_DBL_PAR_TIMELIMIT, remainingTime) );
    GUROBI_CALL_EXC( GRBsetintparam(GRBgetenv(_solver->_model), GRB_INT_PAR_POOLSOLUTIONS,
      std::min(query.maxNumSolutions, 2000000000UL)) );
    GUROBI_CALL_EXC( GRBsetintparam(GRBgetenv(_solver->_model), GRB_INT_PAR_INFUNBDINFO, 1) );

    // Scale the objective vector down before passing it to SCIP.
    double objectiveScalingFactor = 1.0 / std::max(1.0, maxAbsoluteValue(objectiveVector, n) / 1.0e5 );
#if defined(IPO_DEBUG)
    std::cout << "Scaling objective vector by " << objectiveScalingFactor << "." << std::endl;
#endif /* IPO_DEBUG */
    for (std::size_t i = 0; i < n; ++i)
    {
      GUROBI_CALL_EXC( GRBsetdblattrelement(_solver->_model, GRB_DBL_ATTR_OBJ, i,
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

    // First call to Gurobi with dual reductions.
#if defined(IPO_DEBUG)
    std::cout << "Calling Gurobi";
    double timeLimit;
    GUROBI_CALL_EXC( GRBgetdblparam(GRBgetenv(_solver->_model), GRB_DBL_PAR_TIMELIMIT, &timeLimit) );
    if (timeLimit < std::numeric_limits<double>::infinity())
      std::cout << " using a time limit of " << timeLimit << "s";
    std::cout << "." << std::endl;
    GRBwrite(_solver->_model, "GurobiOptimizationOracle.lp");
    GRBwriteparams(GRBgetenv(_solver->_model), "GurobiOptimizationOracle.prm");
#endif /* IPO_DEBUG */

    int retcode = GRBoptimize(_solver->_model);
    if (retcode)
    {
      std::stringstream message;
      message << "GurobiOptimizationOracle<double>: received return code " << retcode << " from GRBoptimize() call.";
      throw std::runtime_error(message.str());
    }
    double runtime, primalBound, dualBound;
    int status;
    GUROBI_CALL_EXC( GRBgetdblattr(_solver->_model, GRB_DBL_ATTR_RUNTIME, &runtime) );
    remainingTime -= runtime;
    GUROBI_CALL_EXC( GRBgetintattr(_solver->_model, GRB_INT_ATTR_STATUS, &status) );
    GUROBI_CALL_EXC( GRBgetdblattr(_solver->_model, GRB_DBL_ATTR_OBJBOUND, &dualBound) );
#if defined(IPO_DEBUG)
    std::cout << "Gurobi returned with status " << status << " and dual bound " << dualBound << "." << std::endl;
#endif /* IPO_DEBUG */

    // Case distinction depending on status.
    if (status == GRB_OPTIMAL || status == GRB_TIME_LIMIT || status == GRB_NODE_LIMIT)
    {
      GUROBI_CALL_EXC( GRBgetdblattr(_solver->_model, GRB_DBL_ATTR_OBJVAL, &primalBound) );
      extractPoints(_solver->_model, _space, objectiveVector, response);
      response.setPrimalBound(primalBound / objectiveScalingFactor);
      response.dualBound = dualBound / objectiveScalingFactor;
      response.hasDualBound = true;
      response.outcome = OptimizationOutcome::FEASIBLE;
      response.hitTimeLimit = status == GRB_TIME_LIMIT;
      if (response.points.empty())
      {
        GUROBI_CALL_EXC( GRBwrite(_solver->_model, "GurobiOptimizationOracle-all-sols-infinite.lp") );
        throw std::runtime_error("GurobiOptimizationOracle<double>: All Gurobi solutions had an infinite entry."
          " Wrote instance to file <GurobiOptimizationOracle-all-sols-infinite.lp>.");
      }
    }
    else if (status == GRB_INFEASIBLE)
    {
      response.outcome = OptimizationOutcome::INFEASIBLE;
    }
    else
    {
      if (status == GRB_INF_OR_UNBD)
      {
        // For the second call we disable dual reductions.
        GUROBI_CALL_EXC( GRBsetintparam(GRBgetenv(_solver->_model), GRB_INT_PAR_DUALREDUCTIONS, 0) );
        GUROBI_CALL_EXC( GRBsetdblparam(GRBgetenv(_solver->_model), GRB_DBL_PAR_TIMELIMIT, remainingTime) );

        // Second call to Gurobi, without dual reductions.
#if defined(IPO_DEBUG)
        std::cout << "Calling Gurobi without dual reductions";
        if (remainingTime < std::numeric_limits<double>::infinity())
          std::cout << " with time limit " << remainingTime << "s";
        std::cout << "." << std::endl;
#endif /* IPO_DEBUG */

        int retcode = GRBoptimize(_solver->_model);
        if (retcode)
        {
          std::stringstream message;
          message << "GurobiOptimizationOracle<double>: received return code " << retcode << " from GRBoptimize() call.";
          throw std::runtime_error(message.str());
        }

        GUROBI_CALL_EXC( GRBgetdblattr(_solver->_model, GRB_DBL_ATTR_RUNTIME, &runtime) );
        remainingTime -= runtime;
        GUROBI_CALL_EXC( GRBgetintattr(_solver->_model, GRB_INT_ATTR_STATUS, &status) );
        GUROBI_CALL_EXC( GRBgetdblattr(_solver->_model, GRB_DBL_ATTR_OBJBOUND, &dualBound) );
#if defined(IPO_DEBUG)
        std::cout << "Gurobi without dual reductions returned with status " << status << " and dual bound "
          << dualBound << "." << std::endl;
#endif /* IPO_DEBUG */
      }

      // Case distinction of status if first call was neither optimal, time limit nor infeasible.
      if (status == GRB_INFEASIBLE)
      {
        response.outcome = OptimizationOutcome::INFEASIBLE;
      }
      else if (status == GRB_UNBOUNDED)
      {
        extractPoints(_solver->_model, _space, objectiveVector, response);
        response.outcome = OptimizationOutcome::UNBOUNDED;

        // We're unbounded but but don't have a point, so we solve a feasibility problem.
        if (response.points.empty())
        {
          GUROBI_CALL_EXC( GRBsetintparam(GRBgetenv(_solver->_model), GRB_INT_PAR_DUALREDUCTIONS, oldDualReductions) );
          GUROBI_CALL_EXC( GRBsetdblparam(GRBgetenv(_solver->_model), GRB_DBL_PAR_TIMELIMIT, remainingTime) );

          // Update objective vector.
          for (std::size_t i = 0; i < n; ++i)
            GUROBI_CALL_EXC( GRBsetdblattrelement(_solver->_model, GRB_DBL_ATTR_OBJ, i, 0.0) );

          // Remove limits on primal and dual bounds.
          // _solver->_boundLimits.maxDualBound = std::numeric_limits<double>::infinity();
          // _solver->_boundLimits.minPrimalBound = -std::numeric_limits<double>::infinity();

          // Third call to Gurobi, for feasibility.
#if defined(IPO_DEBUG)
          std::cout << "Calling Gurobi for feasibility.";
          if (remainingTime < std::numeric_limits<double>::infinity())
            std::cout << " with time limit " << timeLimit << "s";
          std::cout << "." << std::endl;
#endif /* IPO_DEBUG */

          int retcode = GRBoptimize(_solver->_model);
          if (retcode)
          {
            std::stringstream message;
            message << "GurobiOptimizationOracle<double>: received return code " << retcode << " from GRBoptimize() call.";
            throw std::runtime_error(message.str());
          }

          GUROBI_CALL_EXC( GRBgetdblattr(_solver->_model, GRB_DBL_ATTR_RUNTIME, &runtime) );
          remainingTime -= runtime;
          GUROBI_CALL_EXC( GRBgetintattr(_solver->_model, GRB_INT_ATTR_STATUS, &status) );
          GUROBI_CALL_EXC( GRBgetdblattr(_solver->_model, GRB_DBL_ATTR_OBJBOUND, &dualBound) );
#if defined(IPO_DEBUG)
          std::cout << "Gurobi for feasibility returned with status " << status << " and dual bound " << dualBound
            << "." << std::endl;
#endif /* IPO_DEBUG */

          if (status == GRB_INFEASIBLE)
          {
            response.rays.clear();
            response.outcome = OptimizationOutcome::INFEASIBLE;
          }
          else if (status == GRB_OPTIMAL || status == GRB_TIME_LIMIT)
          {
            extractPoints(_solver->_model, _space, objectiveVector, response);
            response.hitTimeLimit = status == GRB_TIME_LIMIT;
            if (response.points.empty())
            {
              std::stringstream message;
              message << "GurobiOptimizationOracle<double>: "
                "Feasbility model has found points all of which have extremely large entries.";
              throw std::runtime_error(message.str());
            }
          }
          else
          {
            std::stringstream message;
            message << "GurobiOptimizationOracle<double>: Unhandled Gurobi status " << status
              << " in feasibility call.";
            throw std::runtime_error(message.str());
          }
        }

        // Solve the LP relaxation to get an unbounded ray.
        {
          GUROBI_CALL_EXC( GRBsetintparam(GRBgetenv(_solver->_model), GRB_INT_PAR_DUALREDUCTIONS, 0) );
          GUROBI_CALL_EXC( GRBsetdblparam(GRBgetenv(_solver->_model), GRB_DBL_PAR_TIMELIMIT, remainingTime) );

          // Update objective vector.
          for (std::size_t i = 0; i < n; ++i)
          {
            GUROBI_CALL_EXC( GRBsetdblattrelement(_solver->_model, GRB_DBL_ATTR_OBJ, i,
              objectiveVector[i] * objectiveScalingFactor) );
          }

          // Remove limits on primal and dual bounds.
          // _solver->_boundLimits.maxDualBound = std::numeric_limits<double>::infinity();
          // _solver->_boundLimits.minPrimalBound = -std::numeric_limits<double>::infinity();

          std::vector<char> vtypesOriginal(n);
          GUROBI_CALL_EXC( GRBgetcharattrarray(_solver->_model, GRB_CHAR_ATTR_VTYPE, 0, n, &vtypesOriginal[0]) );
          std::vector<char> vtypesRelaxed(n, 'C');
          GUROBI_CALL_EXC( GRBsetcharattrarray(_solver->_model, GRB_CHAR_ATTR_VTYPE, 0, n, &vtypesRelaxed[0]) );

          // Fourth call to Gurobi, for unbounded ray.
#if defined(IPO_DEBUG)
          std::cout << "Calling Gurobi for unbounded ray.";
          if (remainingTime < std::numeric_limits<double>::infinity())
            std::cout << " with time limit " << timeLimit << "s";
          std::cout << "." << std::endl;
#endif /* IPO_DEBUG */

          int retcode = GRBoptimize(_solver->_model);
          if (retcode)
          {
            std::stringstream message;
            message << "GurobiOptimizationOracle<double>: received return code " << retcode << " from GRBoptimize() call.";
            throw std::runtime_error(message.str());
          }

          GUROBI_CALL_EXC( GRBgetdblattr(_solver->_model, GRB_DBL_ATTR_RUNTIME, &runtime) );
          remainingTime -= runtime;
          GUROBI_CALL_EXC( GRBgetintattr(_solver->_model, GRB_INT_ATTR_STATUS, &status) );
          GUROBI_CALL_EXC( GRBgetdblattr(_solver->_model, GRB_DBL_ATTR_OBJBOUND, &dualBound) );
#if defined(IPO_DEBUG)
          std::cout << "Gurobi for unbounded ray returned with status " << status << " and dual bound " << dualBound
            << "." << std::endl;
#endif /* IPO_DEBUG */

          GUROBI_CALL_EXC( GRBsetcharattrarray(_solver->_model, GRB_CHAR_ATTR_VTYPE, 0, n, &vtypesOriginal[0]) );
          GUROBI_CALL_EXC( GRBsetintparam(GRBgetenv(_solver->_model), GRB_INT_PAR_DUALREDUCTIONS, oldDualReductions) );

          if (status == GRB_UNBOUNDED)
          {
            extractRays(_solver->_model, _space, response);
          }
          else
          {
            std::stringstream message;
            message << "GurobiOptimizationOracle<double>: Unhandled Gurobi status " << status
              << " in unbounded ray call.";
            throw std::runtime_error(message.str());
          }
        }
      }
      else
      {
        std::stringstream message;
        message << "GurobiOptimizationOracle<double>: Unhandled Gurobi status " << status;
        if (status == GRB_OPTIMAL)
          message << " (optimal) after a previous status " << GRB_INF_OR_UNBD << " (infeasible/unbounded).";
        else if (status == GRB_INF_OR_UNBD)
          message << " (infeasible/unbounded) even after disabling dual reductions.";
        else if (status == GRB_NODE_LIMIT)
        {
          double nodeLimit;
          GUROBI_CALL_EXC( GRBgetdblparam(GRBgetenv(_solver->_model), GRB_DBL_PAR_NODELIMIT, &nodeLimit) );
          message << " (total node limit was " << nodeLimit << ").";
        }
        else
          message << ".";
        throw std::runtime_error(message.str());
      }
    }

    GUROBI_CALL_EXC( GRBsetintparam(GRBgetenv(_solver->_model), GRB_INT_PAR_DUALREDUCTIONS, oldDualReductions) );

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
  GurobiSeparationOracle<double>::GurobiSeparationOracle(std::shared_ptr<GurobiSolver> solver,
    const Constraint<double>& face)
    : SeparationOracle<double>(solver->name()), _solver(solver), _face(face),
    _approximateFace(face)
  {
    _space = solver->space();
    _solver->addFace(&_approximateFace);
  }

  template <>
  GurobiSeparationOracle<double>::~GurobiSeparationOracle()
  {
    _solver->deleteFace(&_approximateFace);
  }

  template <>
  SeparationResponse<double> GurobiSeparationOracle<double>::getInitial(
    const SeparationQuery& query)
  {
    SeparationResponse<double> result;

    if (query.maxNumInequalities > 0 && !_face.isAlwaysSatisfied())
      result.constraints.push_back(_face);

    struct Visitor
    {
      std::shared_ptr<GurobiSolver> solver;
      const SeparationQuery& query;
      SeparationResponse<double>& response;
      std::size_t iteration;
      std::chrono::time_point<std::chrono::system_clock> started;

      void operator()(const Constraint<double>& constraint)
      {
        // Did we reach a limit?
        if (response.hitTimeLimit || response.constraints.size() == query.maxNumInequalities)
          return;

        iteration = (iteration + 1) % GUROBI_SEPARATION_CHECK_TIMELIMIT_FREQUENCY;
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

    GRBiterateBounds(_solver->_model, space(), visitor, true);
    GRBiterateRows(_solver->_model, space(), visitor, true);

    return result;
  }

  template <>
  SeparationResponse<double> GurobiSeparationOracle<double>::separate(const double* vector,
    bool isPoint, const SeparationQuery& query)
  {
    SeparationResponse<double> result;

    _solver->selectFace(&_approximateFace);

    struct Visitor
    {
      std::shared_ptr<GurobiSolver> solver;
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

        iteration = (iteration + 1) % GUROBI_SEPARATION_CHECK_TIMELIMIT_FREQUENCY;
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

        if (activity < lhs - 1.0e-6 || activity > rhs + 1.0e-6)
          response.constraints.push_back(constraint);
      }
    };

    Visitor visitor = { _solver, vector, isPoint, query, result, 0,
      std::chrono::system_clock::now() };

    GRBiterateBounds(_solver->_model, space(), visitor, false);
    GRBiterateRows(_solver->_model, space(), visitor, false);

    return result;
  }

  template <>
  SeparationResponse<double> GurobiSeparationOracle<double>::separateDouble(const double* vector,
    bool isPoint, const SeparationQuery& query)
  {
    return separate(vector, isPoint, query);
  }

#if defined(IPO_RATIONAL_LP)

  template <>
  GurobiSeparationOracle<rational>::GurobiSeparationOracle(std::shared_ptr<GurobiSolver> solver,
    const Constraint<rational>& face)
    : SeparationOracle<rational>(solver->name()), _solver(solver), _face(face),
    _approximateFace(convertConstraint<double>(face))
  {
    _space = solver->space();
    _solver->addFace(&_approximateFace);
  }

  template <>
  GurobiSeparationOracle<rational>::~GurobiSeparationOracle()
  {
    _solver->deleteFace(&_approximateFace);
  }

  template <>
  SeparationResponse<rational> GurobiSeparationOracle<rational>::getInitial(
    const SeparationQuery& query)
  {
    SeparationResponse<rational> result;

    if (query.maxNumInequalities > 0 && !_face.isAlwaysSatisfied())
      result.constraints.push_back(_face);

    struct Visitor
    {
      std::shared_ptr<GurobiSolver> solver;
      const SeparationQuery& query;
      SeparationResponse<rational>& response;
      std::size_t iteration;
      std::chrono::time_point<std::chrono::system_clock> started;

      void operator()(const Constraint<double>& constraint)
      {
        // Did we reach a limit?
        if (response.hitTimeLimit || response.constraints.size() == query.maxNumInequalities)
          return;

        iteration = (iteration + 1) % GUROBI_SEPARATION_CHECK_TIMELIMIT_FREQUENCY;
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

    GRBiterateBounds(_solver->_model, space(), visitor, true);
    GRBiterateRows(_solver->_model, space(), visitor, true);

    return result;
  }

  template <>
  SeparationResponse<rational> GurobiSeparationOracle<rational>::separate(const rational* vector,
    bool isPoint, const SeparationQuery& query)
  {
    SeparationResponse<rational> result;

    _solver->selectFace(&_approximateFace);

    struct Visitor
    {
      std::shared_ptr<GurobiSolver> solver;
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

        iteration = (iteration + 1) % GUROBI_SEPARATION_CHECK_TIMELIMIT_FREQUENCY;
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

        if (activity < lhs - 1.0e-12 || activity > rhs + 1.0e-12)
        {
          response.constraints.push_back(convertConstraint<rational>(constraint));
        }
      }
    };

    Visitor visitor = { _solver, vector, isPoint, query, result, 0,
        std::chrono::system_clock::now() };

    GRBiterateBounds(_solver->_model, space(), visitor, false);
    GRBiterateRows(_solver->_model, space(), visitor, false);

    return result;
  }

  template <>
  SeparationResponse<rational> GurobiSeparationOracle<rational>::separateDouble(const double* vector,
    bool isPoint, const SeparationQuery& query)
  {
    SeparationResponse<rational> result;

    _solver->selectFace(&_approximateFace);

    struct Visitor
    {
      std::shared_ptr<GurobiSolver> solver;
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

        iteration = (iteration + 1) % GUROBI_SEPARATION_CHECK_TIMELIMIT_FREQUENCY;
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

        if (activity < lhs - 1.0e-6 || activity > rhs + 1.0e-6)
        {
          response.constraints.push_back(convertConstraint<rational>(constraint));
        }
      }
    };

    Visitor visitor = { _solver, vector, isPoint, query, result, 0,
        std::chrono::system_clock::now() };

    GRBiterateBounds(_solver->_model, space(), visitor, false);
    GRBiterateRows(_solver->_model, space(), visitor, false);

    return result;
  }

#endif /* IPO_RATIONAL_LP */
  
} /* namespace ipo */
