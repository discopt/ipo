#include "smallest_face.h"

#include <vector>
#include <algorithm>
#include <map>

#include "polar_lp.h"

using namespace soplex;

namespace ipo {

  namespace SmallestFace {

    class NormalConeOracle: public OracleBase, public PolarLP
    {
    public:
      NormalConeOracle(Space& space, UniqueRationalVectorsBase& originalPoints,
        UniqueRationalVectorsBase& originalRays, OracleBase* originalOracle,
        const SVectorRational* target, OutputBase* output)
        : OracleBase("normal cone of " + originalOracle->name(), space),
        PolarLP(originalPoints, originalRays, originalOracle, 1024.0, 1000), _output(output)
      {
        if (space.dimension() == 0)
        {
          for (std::size_t v = 0; v < n(); ++v)
            space.addVariable("a#" + originalOracle->space()[v]);
          space.addVariable("beta");
        }
        else if (space.dimension() != originalOracle->space().dimension() + 1)
        {
          throw std::runtime_error(
            "Spaces differ while constructing NormalConeOracle.");
        }

        OracleBase::initializedSpace();

        /// Initialize Polar LP.

        DSVectorRational vector;
        vector = *target;
        vector.add(n(), Rational(-1));
        _pointConstraint = addConstraint(Rational(0), vector, Rational(0));
        for (std::size_t v = 0; v < n(); ++v)
          setBounds(v, Rational(-1), Rational(1));
      }

      virtual ~NormalConeOracle()
      {

      }

      /**
      * \brief Restricts the oracle to the face defined by \p newFace.
      *
      * Restricts the optimization oracle to the face \f$ F \f$ of \f$ P \f$ defined by \p newFace.
      * For \p newFace equal to \c NULL we define \f$ F := P \f$.
      *
      * This implementation throws an exception since we should not use this oracle to
      * optimize over the faces of a normal cone.
      */

      virtual void setFace(Face* newFace = NULL)
      {
        throw std::runtime_error("NormalConeOracle does not support faces.");
      }

      /**
      * \brief Runs the oracle to maximize the dense rational \p objective.
      *
      * Runs the optimization oracle to maximize the given dense rational \p objective
      * over the current face \f$ F \f$ (see setFace()) and returns \p result.
      * If \p maxHeuristic is less than thisHeuristic() or if the objective value
      * requested by \p objectiveBound is not exceeded, then the call must be forwarded to the
      * next oracle.
      *
      * \param result         After the call, contains the oracle's answer.
      * \param objective      Objective vector \f$ c \in \mathbb{Q}^n \f$ to be maximized.
      * \param objectiveBound Objective value \f$ \gamma \f$ that should be exceeded.
      * \param maxHeuristic   Requested maximum heuristic level.
      * \param minHeuristic   Requested minimum heuristic level.
      *
      * This implementation
      */

      virtual void maximize(OracleResult& result, const soplex::VectorRational& objective,
        const ObjectiveBound& objectiveBound = ObjectiveBound(),
        std::size_t maxHeuristic = std::numeric_limits<std::size_t>::max(),
        std::size_t minHeuristic = 0)
      {
        updateObjective(objective);

        /// Perform stabilization.

        stabilizedPresolve();

        /// Optimize

        optimize(false);

        DSVectorRational* solution = new DSVectorRational;
        getPrimalSolution(*solution);

        result.buildStart(objective);
        if (*solution * objective == 0)
          result.buildAddPoint(solution);
        else
          result.buildAddDirection(solution);
        result.buildFinish(0, true, false, false);
        return result.buildAddPoint(solution);
      }

      friend class InformationBase;


    protected:

      virtual void onBeforeSolve(bool stabilizing)
      {
        _output->onConeBeforeSolve(stabilizing);
      }

      virtual void onAfterSolve(bool stabilizing)
      {
        _output->onConeAfterSolve(stabilizing);
      }

      virtual void onPenaltyDecrease()
      {
        _output->onConePenaltyDecrease();
      }

      virtual void onBeforeCache()
      {
        _output->onConeBeforeCache();
      }

      virtual void onAfterCache(std::size_t numPoints, std::size_t numRays)
      {
         _output->onConeAfterCache(numPoints, numRays);
      }

      virtual void onBeforeOracleCall(bool forceOptimal)
      {
        _output->onConeBeforeOracleCall(forceOptimal);
      }

      virtual void onAfterOracleCall(bool forceOptimal, bool feasible, std::size_t numPoints,
std::size_t numRays,
          bool lastIteration)
      {
        _output->onConeAfterOracleCall(forceOptimal, feasible, numPoints, numRays);
      }

      virtual void onBeforeAddPoint()
      {
        _output->onConeBeforeAddPoint();
      }

      virtual void onAfterAddPoint()
      {
        _output->onConeAfterAddPoint();
      }

      virtual void onBeforeAddRay()
      {
        _output->onConeBeforeAddRay();
      }

      virtual void onAfterAddRay()
      {
        _output->onConeAfterAddRay();
      }

    protected:
      std::size_t _pointConstraint;
      OutputBase* _output;
    };

    class Implementation
    {
    public:
      Implementation(UniqueRationalVectorsBase& points, UniqueRationalVectorsBase& directions,
OracleBase* oracle) :
          _originalPoints(points), _originalDirections(directions), _originalOracle(oracle),
_output(NULL), _normalConeOracle(
          NULL), _dimension(-2)
      {
        _maximizingObjective.reDim(_originalOracle->space().dimension() + 1);
      }

      virtual ~Implementation()
      {

      }

      bool verifyNormalConeElements(OracleResult& result, const SVectorRational& objective,
        const VectorRational& target)
      {
        Rational product = objective * target;
        _originalOracle->maximize(result, objective, ObjectiveBound());
        result.addToContainers(_originalPoints, _originalDirections);
        if (result.isFeasible() && result.points.front().objectiveValue == product)
          return true;
        else
          return false;
      }

      int run(const Point* target, OutputBase& output, LPRowSetRational& coneEquations,
          const std::vector<DSVectorRational*>& freeNormalConeElements,
          const std::vector<const SVectorRational*>& copyNormalConeElements, bool verifyElements)
      {
        std::size_t n = _originalOracle->space().dimension();
        _output = &output;
        _output->_implementation = this;
        _dimension = -2;

        _output->onStart();

        Space normaleConeSpace;
        NormalConeOracle normalConeOracle(normaleConeSpace, _originalPoints,
          _originalDirections, _originalOracle, target, _output);
        _normalConeOracle = &normalConeOracle;
        UniqueRationalVectors conePoints(n + 1);
        UniqueRationalVectors coneRays(n + 1);

        /// Initialize points and directions.

        conePoints.insertFree(new Point);

        DVectorRational denseTarget;
        OracleResult verifyResult;
        std::size_t numElements = freeNormalConeElements.size() + copyNormalConeElements.size();
        if (verifyElements)
        {
          denseTarget.reDim(n, true);
          denseTarget.assign(*target);
          if (verifyElements)
            _output->onBeforeVerifyElements(numElements);
        }
        for (std::size_t i = 0; i < numElements; ++i)
        {
          bool free = i < freeNormalConeElements.size();
          std::size_t index = free ? i : i - freeNormalConeElements.size();
          const SVectorRational& element = free ? *freeNormalConeElements[index] :
*copyNormalConeElements[index];

          if (verifyElements)
          {
            _output->onBeforeOracleVerifyElement(i);
            bool verified = verifyNormalConeElements(verifyResult, element, denseTarget);
            _output->onAfterOracleVerifyElement(verifyResult.points.size(),
verifyResult.directions.size(), verified);
            if (!verified)
              continue;
          }

          Direction* ray;
          if (free)
            ray = freeNormalConeElements[i];
          else
          {
            ray = new Direction;
            *ray = *copyNormalConeElements[i];
          }
          Rational rhs = *ray * *target;
          if (rhs != 0)
            ray->add(n, rhs);

          coneRays.insertFree(ray);
        }
        if (verifyElements)
          _output->onAfterVerifyElements(numElements);

        /// Initialize equations.

        DSVectorRational vector;
        vector = *target;
        vector.add(n, Rational(-1));
        coneEquations.add(Rational(0), vector, Rational(0));

        _output->onAddedInitials(conePoints.size(), coneRays.size(), coneEquations.num());

        /// Start affine hull algorithm.

        AffineHull::OutputBase& hullOutput = output.normalConeHullOutput();
        AffineHull::Result hull;
        hull.run(conePoints, coneRays, coneEquations, &normalConeOracle, hullOutput);
        _dimension = n - hull.dimension();

        /// Extract ray.

        _maximizingObjective.clear();
        for (std::size_t i = 0; i < hull.numSpanningPoints(); ++i)
          _maximizingObjective += *conePoints[hull.spanningPoints()[i]];
        for (std::size_t i = 0; i < hull.numSpanningDirections(); ++i)
          _maximizingObjective += *coneRays[hull.spanningDirections()[i]];

        _output->onEnd();
        _output->_implementation = NULL;
        _normalConeOracle = NULL;

        return _dimension;
      }

      int runAdjacency(const Point* first, const Point* second, OutputBase& output,
          const std::vector<DSVectorRational*>& freeNormalConeElements,
          const std::vector<const SVectorRational*>& copyNormalConeElements, bool verifyElements)
      {
        DVectorRational dense;
        dense.reDim(_originalOracle->space().dimension(), true);
        Point sparse;

        /// Compute difference and create equation valid for normal cone.

        LPRowSetRational coneEquations;
        dense.assign(*first);
        dense -= *second;
        sparse = dense;
        coneEquations.add(Rational(0), sparse, Rational(0));

        /// Compute barycenter.

        dense.clear();
        dense.assign(*first);
        dense += *second;
        sparse.clear();
        for (std::size_t v = 0; v < _originalOracle->space().dimension(); ++v)
        {
          if (dense[v] != 0)
            sparse.add(v, dense[v] / 2);
        }

        return run(&sparse, output, coneEquations, freeNormalConeElements, copyNormalConeElements,
verifyElements);
      }

      friend class InformationBase;
      friend class Result;

    protected:
      UniqueRationalVectorsBase& _originalPoints;
      UniqueRationalVectorsBase& _originalDirections;
      OracleBase* _originalOracle;
      OutputBase* _output;
      NormalConeOracle* _normalConeOracle;

      int _dimension;
      DVectorRational _maximizingObjective;
    };

    Result::Result(UniqueRationalVectorsBase& points, UniqueRationalVectorsBase& directions,
OracleBase* oracle)
    {
      _implementation = new Implementation(points, directions, oracle);
    }

    Result::~Result()
    {
      assert(hasImplementation());
      delete _implementation;
    }

    int Result::dimension() const
    {
      assert(hasImplementation());
      if (_implementation->_dimension < -1)
        throw std::runtime_error("SmallestFace: Dimension not yet determined.");
      return _implementation->_dimension;
    }

    void Result::getMaximizingObjective(soplex::DSVectorRational& maximizingObjective) const
    {
      assert(hasImplementation());
      if (_implementation->_dimension < -1)
        throw std::runtime_error("SmallestFace: Maximizing ray not yet determined.");

      /// Copy maximizing objective by hand since the implementation variable has dimension n+1.

      maximizingObjective.clear();
      for (std::size_t v = 0; v < numVariables(); ++v)
      {
        const Rational& x = _implementation->_maximizingObjective[v];
        if (x != 0)
          maximizingObjective.add(v, x);
      }
    }

    int Result::run(const Point* target, OutputBase& output,
        const std::vector<const soplex::SVectorRational*>& normalConeElements, bool verifyElements)
    {
      LPRowSetRational coneEquations;
      std::vector<DSVectorRational*> freeNormalConeElements;
      return _implementation->run(target, output, coneEquations, freeNormalConeElements,
normalConeElements,
          verifyElements);
    }

    int Result::run(const Point* target, OutputBase& output, const LPRowSetRational&
affineHullEquations)
    {
      LPRowSetRational coneEquations;
      std::vector<DSVectorRational*> freeNormalConeElements;
      std::vector<const SVectorRational*> copyNormalConeElements;
      for (int i = affineHullEquations.num() - 1; i >= 0; --i)
      {
        DSVectorRational* positiveNormal = new DSVectorRational;
        *positiveNormal = affineHullEquations.rowVector(i);
        DSVectorRational* negativeNormal = new DSVectorRational;
        *negativeNormal = affineHullEquations.rowVector(i);
        *negativeNormal *= -1;
        freeNormalConeElements.push_back(positiveNormal);
        freeNormalConeElements.push_back(negativeNormal);
      }
      return _implementation->run(target, output, coneEquations, freeNormalConeElements,
copyNormalConeElements, false);
    }

    int Result::run(const Point* target, OutputBase& output)
    {
      LPRowSetRational coneEquations;
      std::vector<DSVectorRational*> freeNormalConeElements;
      std::vector<const SVectorRational*> copyNormalConeElements;
      return _implementation->run(target, output, coneEquations, freeNormalConeElements,
copyNormalConeElements, false);
    }

    bool Result::isVertex(const Point* target, OutputBase& output,
        const std::vector<const SVectorRational*>& normalConeElements, bool verifyElements)
    {
      return run(target, output, normalConeElements, verifyElements) == 0;
    }

    bool Result::isVertex(const Point* target, OutputBase& output, const LPRowSetRational&
affineHullEquations)
    {
      return run(target, output, affineHullEquations) == 0;
    }

    bool Result::isVertex(const Point* target, OutputBase& output)
    {
      return run(target, output) == 0;
    }

    bool Result::areAdjacent(const Point* first, const Point* second, OutputBase& output,
        const std::vector<const SVectorRational*>& copyNormalConeElements, bool verifyElements)
    {
      std::vector<DSVectorRational*> freeNormalConeElements;
      return _implementation->runAdjacency(first, second, output, freeNormalConeElements,
copyNormalConeElements,
          verifyElements);
    }

    bool Result::areAdjacent(const Point* first, const Point* second, OutputBase& output,
        const LPRowSetRational& affineHullEquations)
    {
      std::vector<DSVectorRational*> freeNormalConeElements;
      std::vector<const SVectorRational*> copyNormalConeElements;
      for (int i = affineHullEquations.num() - 1; i >= 0; --i)
      {
        DSVectorRational* positiveNormal = new DSVectorRational;
        *positiveNormal = affineHullEquations.rowVector(i);
        DSVectorRational* negativeNormal = new DSVectorRational;
        *negativeNormal = affineHullEquations.rowVector(i);
        *negativeNormal *= -1;
        freeNormalConeElements.push_back(positiveNormal);
        freeNormalConeElements.push_back(negativeNormal);
      }
      return _implementation->runAdjacency(first, second, output, freeNormalConeElements,
copyNormalConeElements, false);
    }

    bool Result::areAdjacent(const Point* first, const Point* second, OutputBase& output)
    {
      std::vector<const SVectorRational*> copyNormalConeElements;
      return areAdjacent(first, second, output, copyNormalConeElements, false);
    }

    InformationBase::InformationBase() :
        _implementation(NULL)
    {

    }

    InformationBase::~InformationBase()
    {

    }

    const std::string& InformationBase::oracleName() const
    {
      ensureImplementation();
      return _implementation->_originalOracle->name();
    }

    std::size_t InformationBase::numVariables() const
    {
      ensureImplementation();
      return _implementation->_originalOracle->space().dimension();
    }

    std::size_t InformationBase::numRowsLP() const
    {
      ensureImplementation();
      return _implementation->_normalConeOracle->numRows();
    }

    std::size_t InformationBase::numColumnsLP() const
    {
      ensureImplementation();
      return _implementation->_normalConeOracle->numColumns();
    }

    std::size_t InformationBase::numNonzerosLP() const
    {
      ensureImplementation();
      return _implementation->_normalConeOracle->numNonzeros();
    }

    bool InformationBase::hasImplementation() const
    {
      return _implementation != NULL;
    }

    void InformationBase::ensureImplementation() const
    {
      if (!hasImplementation())
      {
        throw std::runtime_error(
          "Facets.Separation: Output is not associated to an implementation.");
      }
    }

    OutputBase::OutputBase()
    {

    }
    OutputBase::~OutputBase()
    {

    }

    bool OutputBase::isRunning() const
    {
      return hasImplementation();
    }

    void OutputBase::onStart()
    {

    }

    void OutputBase::onBeforeVerifyElements(std::size_t numVerifications)
    {

    }

    void OutputBase::onAfterVerifyElements(std::size_t numVerifications)
    {

    }

    void OutputBase::onBeforeOracleVerifyElement(std::size_t verifyIndex)
    {

    }

    void OutputBase::onAfterOracleVerifyElement(std::size_t numPoints, std::size_t numRays,
      bool verified)
    {

    }

    void OutputBase::onAddedInitials(std::size_t numPoints, std::size_t numRays,
      std::size_t numEquations)
    {

    }

    void OutputBase::onConeBeforeSolve(bool stabilizing)
    {

    }

    void OutputBase::onConeAfterSolve(bool stabilizing)
    {

    }

    void OutputBase::onConePenaltyDecrease()
    {

    }

    void OutputBase::onConeBeforeCache()
    {

    }

    void OutputBase::onConeAfterCache(std::size_t numPoints, std::size_t numRays)
    {

    }

    void OutputBase::onConeBeforeOracleCall(bool forceOptimal)
    {

    }

    void OutputBase::onConeAfterOracleCall(bool forceOptimal, bool feasible, std::size_t numPoints,
      std::size_t numRays)
    {

    }

    void OutputBase::onConeBeforeAddPoint()
    {

    }

    void OutputBase::onConeAfterAddPoint()
    {

    }

    void OutputBase::onConeBeforeAddRay()
    {

    }

    void OutputBase::onConeAfterAddRay()
    {

    }

    void OutputBase::onEnd()
    {

    }

    QuietOutput::QuietOutput()
    {

    }

    QuietOutput::~QuietOutput()
    {

    }

    AffineHull::OutputBase& QuietOutput::normalConeHullOutput()
    {
      return _normalConeHullOutput;
    }

    DebugOutput::DebugOutput()
    {

    }

    DebugOutput::~DebugOutput()
    {

    }

    void DebugOutput::onStart()
    {
      std::cout << "Determining smallest face via the affine hull of the normal cone.\n" <<
std::flush;
    }

    void DebugOutput::onBeforeVerifyElements(std::size_t numVerifications)
    {

    }

    void DebugOutput::onAfterVerifyElements(std::size_t numVerifications)
    {

    }

    void DebugOutput::onBeforeOracleVerifyElement(std::size_t verifyIndex)
    {

    }

    void DebugOutput::onAfterOracleVerifyElement(std::size_t numPoints, std::size_t numRays, bool
verified)
    {

    }

    void DebugOutput::onAddedInitials(std::size_t numPoints, std::size_t numRays, std::size_t
numEquations)
    {

    }

    void DebugOutput::onConeBeforeSolve(bool stabilizing)
    {
      std::cout << "Solving " << (stabilizing ? "Stab-LP" : "Main-LP") << "..." << std::flush;
    }

    void DebugOutput::onConeAfterSolve(bool stabilizing)
    {
      std::cout << " done.\n" << std::flush;
    }

    void DebugOutput::onConeBeforeCache()
    {
      std::cout << "Searching known points and directions:" << std::flush;
    }

    void DebugOutput::onConeAfterCache(std::size_t numPoints, std::size_t numRays)
    {
      std::cout << " (" << numPoints << " known points and " << numRays << " directions).\n" <<
std::flush;
    }

    void DebugOutput::onConeBeforeOracleCall(bool forceOptimal)
    {
      std::cout << "Calling oracle" << (forceOptimal ? "" : " (w/ non-optimal sols)") << ":" <<
std::flush;
    }

    void DebugOutput::onConeAfterOracleCall(bool forceOptimal, bool feasible, std::size_t numPoints,
        std::size_t numRays)
    {
      if (feasible)
      {
        if (numPoints > 0)
          std::cout << " (" << numPoints << " points).\n" << std::flush;
        else
          std::cout << " (feasible).\n" << std::flush;
      }
      else if (numRays > 0)
        std::cout << " (" << numRays << " directions).\n" << std::flush;
      else
        std::cout << " (infeasible).\n" << std::flush;
    }

    void DebugOutput::onEnd()
    {

    }

    AffineHull::OutputBase& DebugOutput::normalConeHullOutput()
    {
      return _normalConeHullOutput;
    }

    ProgressOutput::AffineHullOutput::AffineHullOutput(std::size_t indent)
    {
      _indent.resize(indent, ' ');
      _numVerificationCalls = std::numeric_limits<std::size_t>::max();
      _timeApproxDirections = 0;
      _timeExactDirections = 0;
      _timeCache = 0;
      _timeFactorization = 0;
      _timeStarted = 0;
      _timeLP = 0;
      _timeMainCache = 0;
      _timeMainHeuristics = 0;
      _timeMainOracles = 0;
      _lastTime = 0;
      _lastApprox = false;
      _numMainCache = 0;
      _numMainHeuristics = 0;
      _numMainOracles = 0;
    }

    ProgressOutput::AffineHullOutput::~AffineHullOutput()
    {

    }

    double ProgressOutput::AffineHullOutput::timeStamp()
    {
      double time = _timer.time();
      double elapsed = time - _lastTime;
      _lastTime = time;
      return elapsed;
    }

    void ProgressOutput::AffineHullOutput::onProgress()
    {
      assert(dimensionSafeUpperBound() == dimensionUnsafeUpperBound());
      std::cout << _indent << "Points: " << numSpanningPoints() << ", Rays: " <<
numSpanningDirections() << ",  "
          << dimensionLowerBound() << " <= dim <= " << dimensionSafeUpperBound();
      std::cout << ",  Norm. cone opts: " << numHeuristicCalls() << ".\n" << std::flush;
    }

    void ProgressOutput::AffineHullOutput::onStart()
    {
      _timer.reset();
      _timeApproxDirections = 0;
      _timeExactDirections = 0;
      _timeCache = 0;
      _timeFactorization = 0;
      _timeLP = 0;
      _timeMainCache = 0;
      _timeMainHeuristics = 0;
      _timeMainOracles = 0;
      _timer.start();
      _lastTime = _timer.time();
      _timeStarted = _lastTime;
    }

    void ProgressOutput::AffineHullOutput::onAddedInitialEquations(std::size_t numAllEquations)
    {
      onProgress();
    }

    void ProgressOutput::AffineHullOutput::onBeforeApproximateDirections()
    {
      timeStamp();
    }

    void ProgressOutput::AffineHullOutput::onAfterApproximateDirections(std::size_t numComputed)
    {
      _timeApproxDirections += timeStamp();
    }

    void ProgressOutput::AffineHullOutput::onBeforeExactDirections()
    {
      timeStamp();
    }

    void ProgressOutput::AffineHullOutput::onAfterExactDirections(std::size_t numComputed)
    {
      _timeExactDirections += timeStamp();
    }

    void ProgressOutput::AffineHullOutput::onBeforeCache()
    {
      timeStamp();
    }

    void ProgressOutput::AffineHullOutput::onAfterCache(std::size_t numPoints, std::size_t numRays)
    {
      _timeCache += timeStamp();
    }

    void ProgressOutput::AffineHullOutput::onBeforeOracleZero(bool forceOptimal)
    {

    }

    void ProgressOutput::AffineHullOutput::onAfterOracleZero(std::size_t numPoints)
    {

    }

    void ProgressOutput::AffineHullOutput::onBeforeOracleMaximize(bool forceOptimal)
    {
//      _lastApprox = !forceOptimal;
//      timeStamp();
    }

    void ProgressOutput::AffineHullOutput::onAfterOracleMaximize(std::size_t numPoints, std::size_t
numRays)
    {
//      if (_lastApprox)
//        _timeNormalizationConeOracle += timeStamp();
//      else
//        _timeExactOracles += timeStamp();
    }

    void ProgressOutput::AffineHullOutput::onBeforeOracleMinimize(bool forceOptimal)
    {
//      _lastApprox = !forceOptimal;
//      timeStamp();
    }

    void ProgressOutput::AffineHullOutput::onAfterOracleMinimize(std::size_t numPoints, std::size_t
numRays)
    {
//      if (_lastApprox)
//        _timeNormalizationConeOracle += timeStamp();
//      else
//        _timeExactOracles += timeStamp();
    }

    void ProgressOutput::AffineHullOutput::onBeforeOracleVerify(std::size_t verifyIndex)
    {
//      std::cout << _indent << "(calling exact oracle " << (verifyIndex + 1) << "/" <<
// 7_numVerificationCalls << "..."
//          << std::flush;
    }

    void ProgressOutput::AffineHullOutput::onAfterOracleVerify(std::size_t numPoints, std::size_t
numRays)
    {
//      _timeExactOracles += timeStamp();
//      std::cout << " done)" << std::endl;
    }

    void ProgressOutput::AffineHullOutput::onBeforePoint(bool twoPoints)
    {
      timeStamp();
    }

    void ProgressOutput::AffineHullOutput::onAfterPoint(bool twoPoints)
    {
      onProgress();
      _timeFactorization += timeStamp();
    }

    void ProgressOutput::AffineHullOutput::onBeforeDirection()
    {
      timeStamp();
    }

    void ProgressOutput::AffineHullOutput::onAfterDirection()
    {
      onProgress();
      _timeFactorization += timeStamp();
    }

    void ProgressOutput::AffineHullOutput::onPotentialEquation()
    {
      onProgress();
    }

    void ProgressOutput::AffineHullOutput::onEquation()
    {
      onProgress();
    }

    void ProgressOutput::AffineHullOutput::onBeforeVerifyImmediate()
    {
      _numVerificationCalls = 2;
    }

    void ProgressOutput::AffineHullOutput::onBeforeVerifyDelayed(std::size_t numVerifications)
    {
      assert(false);
      _numVerificationCalls = numVerifications + 1;
    }

    void ProgressOutput::AffineHullOutput::onAfterVerifyDelayed(std::size_t numVerifications)
    {
      onProgress();
    }

    void ProgressOutput::AffineHullOutput::onEnd()
    {
      std::cout << _indent << "Max. direction bitsize: " << maxDirectionBitsize() << "\n";
      std::cout << _indent << "Number of invocations:            LP: " << std::setw(6) <<
numHeuristicCalls()
          << ", Main heuristics: " << std::setw(6) << _numMainHeuristics << ", Main oracles: " <<
std::setw(6)
          << _numMainOracles << ", Main cache: " << std::setw(6) << _numMainCache << ", Cache: " <<
std::setw(6)
          << numCacheQueries() << ", Directions: " << std::setw(6) <<
numApproximateDirectionSolves() << " (approx.), "
          << std::setw(6) << numExactDirectionSolves() << " (exact), Factorizations: " <<
std::setw(6)
          << (dimensionLowerBound() + 1) << "\n";
      std::cout << _indent << "Timings (s): Overall:   " << std::setw(8) << (_timer.time() -
_timeStarted) << ", LP: "
          << std::setw(6) << _timeLP << ", Main heuristics: " << std::setw(6) << _timeMainHeuristics
          << ", Main oracles: " << std::setw(6) << _timeMainOracles << ", Main cache: " <<
std::setw(6)
          << _timeMainCache << ", Cache: " << std::setw(6) << _timeCache << ", Directions: " <<
std::setw(6)
          << _timeApproxDirections << " (approx.), " << std::setw(6) << _timeExactDirections
          << " (exact), Factorizations: " << std::setw(6) << _timeFactorization << std::endl;
      _timer.stop();
    }

    ProgressOutput::ProgressOutput(std::size_t indent) :
        _normalConeHullOutput(indent), _numVerifications(0)
    {
      _indent.resize(indent, ' ');
    }

    ProgressOutput::~ProgressOutput()
    {

    }

    void ProgressOutput::onBeforeVerifyElements(std::size_t numVerifications)
    {
      _numVerifications = numVerifications;
    }

    void ProgressOutput::onBeforeOracleVerifyElement(std::size_t verifyIndex)
    {
      std::cout << _indent << "Verifying ray " << (verifyIndex + 1) << "/" << _numVerifications
          << " with enforcing oracle..." << std::flush;
    }

    void ProgressOutput::onAfterOracleVerifyElement(std::size_t numPoints, std::size_t numRays, bool
verified)
    {
      std::cout << (verified ? " verified.\n" : " failed!\n") << std::flush;
    }

    void ProgressOutput::onAddedInitials(std::size_t numPoints, std::size_t numRays, std::size_t
numEquations)
    {
      std::cout << _indent << "Added " << numPoints << " initial points, " << numRays <<
        " directions, and " << numEquations
          << " equations.\n" << std::flush;
    }

    void ProgressOutput::onConeBeforeSolve(bool stabilizing)
    {
      _normalConeHullOutput.timeStamp();
      std::cout << _indent << (stabilizing ?  "  Stab-LP " : "  Main-LP ") << numRowsLP() << "x" <<
numColumnsLP() << ": " << std::flush;
    }

    void ProgressOutput::onConeAfterSolve(bool stabilizing)
    {
      _normalConeHullOutput._timeLP += _normalConeHullOutput.timeStamp();
      std::cout << "solved." << std::flush;
    }

    void ProgressOutput::onConePenaltyDecrease()
    {
      std::cout << " Decreasing penalty.\n" << std::flush;
    }

    void ProgressOutput::onConeBeforeCache()
    {
      _normalConeHullOutput.timeStamp();
      ++_normalConeHullOutput._numMainCache;
      std::cout << " Cache: " << std::flush;
    }

    void ProgressOutput::onConeAfterCache(std::size_t numPoints, std::size_t numRays)
    {
      _normalConeHullOutput._timeMainCache += _normalConeHullOutput.timeStamp();
      if (numPoints > 0)
      {
        if (numRays > 0)
          std::cout << numPoints << " points, " << numRays << " directions.\n" << std::flush;
        else
          std::cout << numPoints << " points.\n" << std::flush;
      }
      else if (numRays > 0)
        std::cout << numRays << " directions.\n" << std::flush;
      else
        std::cout << "nothing." << std::flush;
    }

    void ProgressOutput::onConeBeforeOracleCall(bool forceOptimal)
    {
      if (forceOptimal)
        ++_normalConeHullOutput._numMainOracles;
      else
        ++_normalConeHullOutput._numMainHeuristics;
      _normalConeHullOutput.timeStamp();

      if (forceOptimal)
        std::cout << " Oracle: " << std::flush;
      else
        std::cout << " Heuristic: " << std::flush;
    }

    void ProgressOutput::onConeAfterOracleCall(bool forceOptimal, bool feasible, std::size_t
      numPoints, std::size_t numRays)
    {
      double time = _normalConeHullOutput.timeStamp();
      if (forceOptimal)
        _normalConeHullOutput._timeMainOracles += time;
      else
        _normalConeHullOutput._timeMainHeuristics += time;

      if (feasible)
      {
        if (numPoints > 0)
          std::cout << numPoints << " points.\n" << std::flush;
        else if (numRays > 0)
          std::cout << numRays << " directions.\n" << std::flush;
        else
        {
          std::cout << "feasible." << (forceOptimal ? "\n" : "") << std::flush;
        }
      }
      else
        std::cout << "infeasible!" << std::flush;
    }

    AffineHull::OutputBase& ProgressOutput::normalConeHullOutput()
    {
      return _normalConeHullOutput;
    }

  }

}
