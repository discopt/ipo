#include "smallest_face.h"

#include <vector>
#include <algorithm>
#include <map>

#include "polar_lp.h"

using namespace soplex;

namespace ipo {

  namespace SmallestFace {

    class NormalConeOptimizationOracle: public OptimizationOracleBase, public PolarLP
    {
    public:
      NormalConeOptimizationOracle(UniqueRationalVectorsBase& originalPoints, UniqueRationalVectorsBase& originalRays,
          OptimizationOracleBase* originalOracle, const SVectorRational* target, OutputBase* output) :
          OptimizationOracleBase("normal.cone(" + originalOracle->name() + ")"), PolarLP(originalPoints, originalRays,
              originalOracle), _output(output)
      {
        /// Initialize oracle.

        std::vector<std::string> varNames(n() + 1);
        for (std::size_t v = 0; v < n(); ++v)
          varNames[v] = "a#" + originalOracle->variableName(v);
        varNames[n()] = "beta";
        initialize(varNames);

        /// Initialize Polar LP.

        DSVectorRational vector;
        vector = *target;
        vector.add(n(), Rational(-1));
        _pointConstraint = addConstraint(Rational(0), vector, Rational(0));
        for (std::size_t v = 0; v < n(); ++v)
          setBounds(v, Rational(-1), Rational(1));
      }

      virtual ~NormalConeOptimizationOracle()
      {

      }

      friend class InformationBase;

    protected:

      virtual void run(OptimizationResult& result, const VectorRational& objective, const Rational* improveValue,
          bool approximate)
      {
        updateObjective(objective);

        optimize();

        DSVectorRational* solution = new DSVectorRational;
        getPrimalSolution(*solution);

        result.points.push_back(solution);
        result.setFeasible(objective);
        if (result.bestValue > 0)
        {
          result.directions.push_back(result.points.front());
          result.points.pop_back();
          result.setUnbounded();
        }
        result.optimal = true;
      }

      virtual void onBeforeSolve()
      {
        _output->onConeBeforeSolve();
      }

      virtual void onAfterSolve()
      {
        _output->onConeAfterSolve();
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

      virtual void onAfterOracleCall(bool forceOptimal, bool feasible, std::size_t numPoints, std::size_t numRays)
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
      Implementation(UniqueRationalVectorsBase& points, UniqueRationalVectorsBase& rays, OptimizationOracleBase* oracle) :
          _originalPoints(points), _originalRays(rays), _originalOracle(oracle), _output(NULL), _normalConeOptimizationOracle(
          NULL), _dimension(-2)
      {
        _maximizingObjective.reDim(_originalOracle->numVariables() + 1);
      }

      virtual ~Implementation()
      {

      }

      bool verifyNormalConeElements(OptimizationResult& result, const SVectorRational& objective,
          const VectorRational& target)
      {
        Rational product = objective * target;
        _originalOracle->maximize(result, objective, true);
        for (std::size_t i = 0; i < result.points.size(); ++i)
          _originalPoints.insertFree(result.points[i]);
        for (std::size_t i = 0; i < result.directions.size(); ++i)
          _originalRays.insertFree(result.directions[i]);
        if (result.isFeasible() && result.bestValue == product)
          return true;
        else
          return false;
      }

      int run(const Point* target, OutputBase& output, LPRowSetRational& coneEquations,
          const std::vector<DSVectorRational*>& freeNormalConeElements,
          const std::vector<const SVectorRational*>& copyNormalConeElements, bool verifyElements)
      {
        std::size_t n = _originalOracle->numVariables();
        _output = &output;
        _output->_implementation = this;
        _dimension = -2;

        _output->onStart();

        NormalConeOptimizationOracle normalConeOracle(_originalPoints, _originalRays, _originalOracle, target, _output);
        _normalConeOptimizationOracle = &normalConeOracle;
        UniqueRationalVectors conePoints(n + 1);
        UniqueRationalVectors coneRays(n + 1);

        /// Initialize points and rays.

        conePoints.insertFree(new Point);

        DVectorRational denseTarget;
        OptimizationResult verifyResult;
        std::size_t numElements = freeNormalConeElements.size() + copyNormalConeElements.size();
        if (verifyElements)
        {
          verifyResult.reset(n);
          denseTarget.reDim(n, true);
          denseTarget.assign(*target);
          if (verifyElements)
            _output->onBeforeVerifyElements(numElements);
        }
        for (std::size_t i = 0; i < numElements; ++i)
        {
          bool free = i < freeNormalConeElements.size();
          std::size_t index = free ? i : i - freeNormalConeElements.size();
          const SVectorRational& element = free ? *freeNormalConeElements[index] : *copyNormalConeElements[index];

          if (verifyElements)
          {
            _output->onBeforeOracleVerifyElement(i);
            bool verified = verifyNormalConeElements(verifyResult, element, denseTarget);
            _output->onAfterOracleVerifyElement(verifyResult.points.size(), verifyResult.directions.size(), verified);
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
        _normalConeOptimizationOracle = NULL;

        return _dimension;
      }

      int runAdjacency(const Point* first, const Point* second, OutputBase& output,
          const std::vector<DSVectorRational*>& freeNormalConeElements,
          const std::vector<const SVectorRational*>& copyNormalConeElements, bool verifyElements)
      {
        DVectorRational dense;
        dense.reDim(_originalOracle->numVariables(), true);
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
        for (std::size_t v = 0; v < _originalOracle->numVariables(); ++v)
        {
          if (dense[v] != 0)
            sparse.add(v, dense[v] / 2);
        }

        return run(&sparse, output, coneEquations, freeNormalConeElements, copyNormalConeElements, verifyElements);
      }

      friend class InformationBase;
      friend class Result;

    protected:
      UniqueRationalVectorsBase& _originalPoints;
      UniqueRationalVectorsBase& _originalRays;
      OptimizationOracleBase* _originalOracle;
      OutputBase* _output;
      NormalConeOptimizationOracle* _normalConeOptimizationOracle;

      int _dimension;
      DVectorRational _maximizingObjective;
    };

    Result::Result(UniqueRationalVectorsBase& points, UniqueRationalVectorsBase& rays, OptimizationOracleBase* oracle)
    {
      _implementation = new Implementation(points, rays, oracle);
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
      return _implementation->run(target, output, coneEquations, freeNormalConeElements, normalConeElements,
          verifyElements);
    }

    int Result::run(const Point* target, OutputBase& output, const LPRowSetRational& affineHullEquations)
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
      return _implementation->run(target, output, coneEquations, freeNormalConeElements, copyNormalConeElements, false);
    }

    int Result::run(const Point* target, OutputBase& output)
    {
      LPRowSetRational coneEquations;
      std::vector<DSVectorRational*> freeNormalConeElements;
      std::vector<const SVectorRational*> copyNormalConeElements;
      return _implementation->run(target, output, coneEquations, freeNormalConeElements, copyNormalConeElements, false);
    }

    bool Result::isVertex(const Point* target, OutputBase& output,
        const std::vector<const SVectorRational*>& normalConeElements, bool verifyElements)
    {
      return run(target, output, normalConeElements, verifyElements) == 0;
    }

    bool Result::isVertex(const Point* target, OutputBase& output, const LPRowSetRational& affineHullEquations)
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
      return _implementation->runAdjacency(first, second, output, freeNormalConeElements, copyNormalConeElements,
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
      return _implementation->runAdjacency(first, second, output, freeNormalConeElements, copyNormalConeElements, false);
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
      return _implementation->_originalOracle->numVariables();
    }

    std::size_t InformationBase::numRowsLP() const
    {
      ensureImplementation();
      return _implementation->_normalConeOptimizationOracle->numRows();
    }

    std::size_t InformationBase::numColumnsLP() const
    {
      ensureImplementation();
      return _implementation->_normalConeOptimizationOracle->numColumns();
    }

    std::size_t InformationBase::numNonzerosLP() const
    {
      ensureImplementation();
      return _implementation->_normalConeOptimizationOracle->numNonzeros();
    }

    bool InformationBase::hasImplementation() const
    {
      return _implementation != NULL;
    }

    void InformationBase::ensureImplementation() const
    {
      if (!hasImplementation())
        throw std::runtime_error("Facets.Separation: Output is not associated to an implementation.");
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

    void OutputBase::onAfterOracleVerifyElement(std::size_t numPoints, std::size_t numRays, bool verified)
    {

    }

    void OutputBase::onAddedInitials(std::size_t numPoints, std::size_t numRays, std::size_t numEquations)
    {

    }

    void OutputBase::onConeBeforeSolve()
    {

    }

    void OutputBase::onConeAfterSolve()
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

    void OutputBase::onConeAfterOracleCall(bool forceOptimal, bool feasible, std::size_t numPoints, std::size_t numRays)
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
      std::cout << "Determining smallest face via the affine hull of the normal cone.\n" << std::flush;
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

    void DebugOutput::onAfterOracleVerifyElement(std::size_t numPoints, std::size_t numRays, bool verified)
    {

    }

    void DebugOutput::onAddedInitials(std::size_t numPoints, std::size_t numRays, std::size_t numEquations)
    {

    }

    void DebugOutput::onConeBeforeSolve()
    {
      std::cout << "Solving target cut LP..." << std::flush;
    }

    void DebugOutput::onConeAfterSolve()
    {
      std::cout << " done.\n" << std::flush;
    }

    void DebugOutput::onConeBeforeCache()
    {
      std::cout << "Searching known points and rays:" << std::flush;
    }

    void DebugOutput::onConeAfterCache(std::size_t numPoints, std::size_t numRays)
    {
      std::cout << " (" << numPoints << " known points and " << numRays << " rays).\n" << std::flush;
    }

    void DebugOutput::onConeBeforeOracleCall(bool forceOptimal)
    {
      std::cout << "Calling oracle" << (forceOptimal ? "" : " (w/ non-optimal sols)") << ":" << std::flush;
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
        std::cout << " (" << numRays << " rays).\n" << std::flush;
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
      std::cout << _indent << "Points: " << numSpanningPoints() << ", Rays: " << numSpanningDirections() << ",  "
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

    void ProgressOutput::AffineHullOutput::onAfterOracleMaximize(std::size_t numPoints, std::size_t numRays)
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

    void ProgressOutput::AffineHullOutput::onAfterOracleMinimize(std::size_t numPoints, std::size_t numRays)
    {
//      if (_lastApprox)
//        _timeNormalizationConeOracle += timeStamp();
//      else
//        _timeExactOracles += timeStamp();
    }

    void ProgressOutput::AffineHullOutput::onBeforeOracleVerify(std::size_t verifyIndex)
    {
//      std::cout << _indent << "(calling exact oracle " << (verifyIndex + 1) << "/" << _numVerificationCalls << "..."
//          << std::flush;
    }

    void ProgressOutput::AffineHullOutput::onAfterOracleVerify(std::size_t numPoints, std::size_t numRays)
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

    void ProgressOutput::AffineHullOutput::onBeforeRay()
    {
      timeStamp();
    }

    void ProgressOutput::AffineHullOutput::onAfterRay()
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
      std::cout << _indent << "Number of invocations:            LP: " << std::setw(6) << numHeuristicCalls()
          << ", Main heuristics: " << std::setw(6) << _numMainHeuristics << ", Main oracles: " << std::setw(6)
          << _numMainOracles << ", Main cache: " << std::setw(6) << _numMainCache << ", Cache: " << std::setw(6)
          << numCacheQueries() << ", Directions: " << std::setw(6) << numApproximateDirectionSolves() << " (approx.), "
          << std::setw(6) << numExactDirectionSolves() << " (exact), Factorizations: " << std::setw(6)
          << (dimensionLowerBound() + 1) << "\n";
      std::cout << _indent << "Timings (s): Overall:   " << std::setw(8) << (_timer.time() - _timeStarted) << ", LP: "
          << std::setw(6) << _timeLP << ", Main heuristics: " << std::setw(6) << _timeMainHeuristics
          << ", Main oracles: " << std::setw(6) << _timeMainOracles << ", Main cache: " << std::setw(6)
          << _timeMainCache << ", Cache: " << std::setw(6) << _timeCache << ", Directions: " << std::setw(6)
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

    void ProgressOutput::onAfterOracleVerifyElement(std::size_t numPoints, std::size_t numRays, bool verified)
    {
      std::cout << (verified ? " verified.\n" : " failed!\n") << std::flush;
    }

    void ProgressOutput::onAddedInitials(std::size_t numPoints, std::size_t numRays, std::size_t numEquations)
    {
      std::cout << _indent << "Added " << numPoints << " initial points, " << numRays << " rays, and " << numEquations
          << " equations.\n" << std::flush;
    }

    void ProgressOutput::onConeBeforeSolve()
    {
      _normalConeHullOutput.timeStamp();
      std::cout << _indent << "  LP " << numRowsLP() << "x" << numColumnsLP() << ": " << std::flush;
    }

    void ProgressOutput::onConeAfterSolve()
    {
      _normalConeHullOutput._timeLP += _normalConeHullOutput.timeStamp();
      std::cout << "solved." << std::flush;
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
          std::cout << numPoints << " points, " << numRays << " rays.\n" << std::flush;
        else
          std::cout << numPoints << " points.\n" << std::flush;
      }
      else if (numRays > 0)
        std::cout << numRays << " rays.\n" << std::flush;
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

    void ProgressOutput::onConeAfterOracleCall(bool forceOptimal, bool feasible, std::size_t numPoints,
        std::size_t numRays)
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
          std::cout << numRays << " rays.\n" << std::flush;
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
