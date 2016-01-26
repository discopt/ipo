#include "facets.h"

#include <vector>
#include <algorithm>
#include <map>

#include "polar_lp.h"

using namespace soplex;

namespace ipo {

  namespace Separation {

    class Implementation: public PolarLP
    {
    public:
      Implementation(UniqueRationalVectorsBase& points, UniqueRationalVectorsBase& rays,
          const VectorSubset& spanningPoints, const VectorSubset& spanningRays, const VectorSubset& columnBasis,
          OptimizationOracleBase* oracle) :
          PolarLP(points, rays, oracle, 16.0, 30), _separatingPoint(false), _output(NULL), _separatedEquation(false), _separatedFacet(
              false)
      {
        SVectorRational vector;
        _normalizationConstraint = addConstraint(-infinity, vector, Rational(0));
        _dimension = int(columnBasis.size()) - 1;

        /// Go through spanning points.

        _interiorPoint.reDim(n(), true);
        for (std::size_t i = 0; i < spanningPoints.size(); ++i)
        {
          std::size_t index = spanningPoints[i];
          _basisConstraints.push_back(addPointContraint(index));
          _basisConstraintTypes.push_back('p');
          _basisConstraintIndices.push_back(index);
          _interiorPoint += *_points[index];
        }
        for (std::size_t v = 0; v < n(); ++v)
          _interiorPoint[v] /= int(spanningPoints.size());

        /// Go through spanning rays.

        _interiorRay.reDim(n(), true);
        for (std::size_t i = 0; i < spanningRays.size(); ++i)
        {
          std::size_t index = spanningRays[i];
          _basisConstraints.push_back(addRayContraint(index));
          _basisConstraintTypes.push_back('r');
          _basisConstraintIndices.push_back(index);
          _interiorPoint += *_rays[index];
          _interiorRay += *_rays[index];
        }

        /// Set initial basis.

        _basis.columnStatus.resize(n() + 1, SPxSolver::ZERO);
        for (std::size_t i = 0; i < columnBasis.size(); ++i)
          _basis.columnStatus[columnBasis[i]] = SPxSolver::BASIC;
        _basis.constraintStatus.resize(numConstraints(), SPxSolver::ON_UPPER);
        _basis.constraintStatus[_normalizationConstraint] = SPxSolver::BASIC;
      }

      virtual ~Implementation()
      {

      }

      virtual void onBeforeSolve(bool stabilizing)
      {
        _output->onBeforeSolve(stabilizing, currentStabilizationPenalty());
      }

      virtual void onAfterSolve(bool stabilizing)
      {
        _output->onAfterSolve(stabilizing, currentStabilizationPenalty(), lastMainObjective(),
            lastStabilizationPenaltyCosts());
      }

      virtual void onPenaltyDecrease()
      {
        _output->onPenaltyDecrease(currentStabilizationPenalty());
      }

      virtual void onBeforeCache()
      {
        _output->onBeforeCache();
      }

      virtual void onAfterCache(std::size_t numPoints, std::size_t numRays)
      {
        _output->onAfterCache(numPoints, numRays);
      }

      virtual void onBeforeOracleCall(bool forceOptimal)
      {
        _output->onBeforeOracleCall(forceOptimal);
      }

      virtual void onAfterOracleCall(bool forceOptimal, bool feasible, std::size_t numPoints, std::size_t numRays,
          bool lastIteration)
      {
        _output->onAfterOracleCall(forceOptimal, feasible, numPoints, numRays, lastIteration);
      }

      virtual void onBeforeAddPoint()
      {
        _output->onBeforeAddPoint();
      }

      virtual void onAfterAddPoint()
      {
        _output->onAfterAddPoint();
      }

      virtual void onBeforeAddRay()
      {
        _output->onBeforeAddRay();
      }

      virtual void onAfterAddRay()
      {
        _output->onAfterAddRay();
      }

      bool run(const SVectorRational* target, bool separatingPoint, OutputBase& output)
      {
        _output = &output;
        _output->_implementation = this;
        _separatedFacet = false;
        _separatedEquation = false;
        _separatingPoint = separatingPoint;

        _output->onStart();

        /// Set objective.

        DVectorRational dense;
        dense.reDim(n() + 1, true);
        for (int p = target->size() - 1; p >= 0; --p)
          dense[target->index(p)] = target->value(p);
        if (separatingPoint)
          dense[n()] = -1;
        updateObjective(dense);

        /// Set normalization

        dense.reDim(n());
        if (separatingPoint)
          dense -= _interiorPoint;
        else
          dense -= _interiorRay;
        updateConstraint(_normalizationConstraint, -infinity, dense, 1);

        /// Perform stabilization.
        
        stabilizedPresolve();

        /// Reset basis.

        setBasis(_basis);

        /// Optimize
        
        optimize();
        
        /// Extract solution

        dense.reDim(n() + 1);
        getPrimalSolution(dense);
        _inequality.setRhs(dense[n()]);
        DSVectorRational sparse;
        for (std::size_t v = 0; v < n(); ++v)
        {
          if (dense[v] != 0)
            sparse.add(v, dense[v]);
        }
        _inequality.setRowVector(sparse);
        _violation = getObjectiveValue();

        /// Use basis information for the certificate.

        Basis basis;
        getBasis(basis);
        _certificate.pointIndices.clear();
        std::copy(basis.tightPoints.begin(), basis.tightPoints.end(), std::back_inserter(_certificate.pointIndices));
        _certificate.directionIndices.clear();
        std::copy(basis.tightRays.begin(), basis.tightRays.end(), std::back_inserter(_certificate.directionIndices));
        for (std::size_t i = 0; i < _basisConstraints.size(); ++i)
        {
          std::size_t constraint = _basisConstraints[i];
          if (basis.constraintStatus[constraint] != SPxSolver::ON_UPPER)
            continue;
          if (_basisConstraintTypes[i] == 'p')
            _certificate.pointIndices.push_back(_basisConstraintIndices[i]);
          else
            _certificate.directionIndices.push_back(_basisConstraintIndices[i]);
        }

        /// Based on certificate we know what we separated.

        if (_certificate.pointIndices.size() + _certificate.directionIndices.size() == _dimension + 1)
        {
          _separatedEquation = true;
          _inequality.setLhs(_inequality.rhs());
        }
        else
        {
          _inequality.setLhs(-infinity);
          if (_certificate.pointIndices.size() + _certificate.directionIndices.size() == _dimension)
            _separatedFacet = true;
          else
          {
            std::cerr << "\n!!! IPO computed lower-dimensional face. This is a BUG, probably caused by a SoPlex restart !!!" << std::endl;
            exit(1);
          }
        }

        _output->onEnd(_violation > 0);
        _output->_implementation = NULL;

        return _violation > 0;
      }

      friend class InformationBase;
      friend class Result;

    protected:
      std::size_t n() const
      {
        return _oracle->numVariables();
      }

    protected:
      bool _separatingPoint;
      OutputBase* _output;
      int _dimension;
      DVectorRational _interiorPoint;
      DVectorRational _interiorRay;
      std::size_t _normalizationConstraint;
      std::vector<std::size_t> _basisConstraints;
      std::vector<char> _basisConstraintTypes;
      std::vector<std::size_t> _basisConstraintIndices;
      Basis _basis;

      /// Result.

      LPRowRational _inequality;
      Certificate _certificate;
      Rational _violation;
      bool _separatedFacet;
      bool _separatedEquation;
    };

    Result::Result(UniqueRationalVectorsBase& points, UniqueRationalVectorsBase& rays,
        const VectorSubset& spanningPoints, const VectorSubset& spanningRays, const VectorSubset& columnBasis,
        OptimizationOracleBase* oracle)
    {
      _implementation = new Implementation(points, rays, spanningPoints, spanningRays, columnBasis, oracle);
    }

    Result::~Result()
    {
      assert(hasImplementation());
      delete _implementation;
    }

    void Result::inequality(soplex::LPRowRational& inequality) const
    {
      inequality = _implementation->_inequality;
    }

    void Result::inequality(soplex::LPRowSetRational& inequalities) const
    {
      inequalities.add(-infinity, _implementation->_inequality.rowVector(), _implementation->_inequality.rhs());
    }

    void Result::certificate(Certificate& certificate) const
    {
      certificate = _implementation->_certificate;
    }

    const Rational& Result::violation() const
    {
      return _implementation->_violation;
    }

    bool Result::separatePoint(const Point* targetPoint, OutputBase& output)
    {
      return _implementation->run(targetPoint, true, output);
    }

    bool Result::separateRay(const Point* targetRay, OutputBase& output)
    {
      return _implementation->run(targetRay, false, output);
    }

    InformationBase::InformationBase() :
        _implementation(NULL)
    {

    }
    InformationBase::~InformationBase()
    {

    }

    bool InformationBase::separatingPoint() const
    {
      ensureImplementation();
      return _implementation->_separatingPoint;
    }

    bool InformationBase::separatedFacet() const
    {
      ensureImplementation();
      return _implementation->_separatedFacet;
    }

    bool InformationBase::separatedEquation() const
    {
      ensureImplementation();
      return _implementation->_separatedEquation;
    }

    const std::string& InformationBase::oracleName() const
    {
      ensureImplementation();
      return _implementation->_oracle->name();
    }

    std::size_t InformationBase::numVariables() const
    {
      ensureImplementation();
      return _implementation->n();
    }

    std::size_t InformationBase::numRowsLP() const
    {
      ensureImplementation();
      return _implementation->numRows();
    }

    std::size_t InformationBase::numColumnsLP() const
    {
      ensureImplementation();
      return _implementation->numColumns();
    }

    std::size_t InformationBase::numNonzerosLP() const
    {
      ensureImplementation();
      return _implementation->numNonzeros();
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

    void OutputBase::onEnd(bool separated)
    {

    }

    void OutputBase::onBeforeSolve(bool stabilizing, double penalty)
    {

    }

    void OutputBase::onAfterSolve(bool stabilizing, double penalty, double mainObjective, double penaltyCosts)
    {

    }

    void OutputBase::onPenaltyDecrease(double penalty)
    {

    }

    void OutputBase::onBeforeCache()
    {

    }

    void OutputBase::onAfterCache(std::size_t numPoints, std::size_t numRays)
    {

    }

    void OutputBase::onBeforeOracleCall(bool forceOptimal)
    {

    }

    void OutputBase::onAfterOracleCall(bool forceOptimal, bool feasible, std::size_t numPoints, std::size_t numRays,
        bool lastIteration)
    {

    }

    void OutputBase::onBeforeAddPoint()
    {

    }

    void OutputBase::onAfterAddPoint()
    {

    }

    void OutputBase::onBeforeAddRay()
    {

    }

    void OutputBase::onAfterAddRay()
    {

    }

    QuietOutput::QuietOutput()
    {

    }

    QuietOutput::~QuietOutput()
    {

    }

    DebugOutput::DebugOutput()
    {

    }

    DebugOutput::~DebugOutput()
    {

    }

    void DebugOutput::onStart()
    {
      std::cout << "Starting to separate a " << (separatingPoint() ? "point" : "ray") << ".\n" << std::flush;
    }

    void DebugOutput::onEnd(bool separated)
    {
      if (separated)
      {
        std::cout << "Separated a " << (separatingPoint() ? "point" : "ray") << " by "
            << (separatedEquation() ? "an equation" : (separatedFacet() ? "a facet" : "a face")) << ".\n" << std::flush;
      }
      else
      {
        std::cout << "Given " << (separatingPoint() ? "point" : "ray") << " is feasible.\n" << std::flush;
      }
    }

    void DebugOutput::onBeforeSolve(bool stabilizing, double penalty)
    {
      std::cout << "Solving " << (stabilizing ? "stabilization" : "main") << " LP ";
      if (penalty > 0.0)
        std::cout << " with penalty " << penalty << " ";
      std::cout << std::flush;
    }

    void DebugOutput::onAfterSolve(bool stabilizing, double penalty, double mainObjective, double penaltyCosts)
    {
      std::cout << " done (obj" << mainObjective << " - " << penaltyCosts << ").\n" << std::flush;
    }

    void DebugOutput::onPenaltyDecrease(double penalty)
    {
      std::cout << "Decreasing stabilization penalty to " << penalty << ".\n" << std::flush;
    }

    void DebugOutput::onBeforeCache()
    {
      std::cout << "Searching known points and rays:" << std::flush;
    }

    void DebugOutput::onAfterCache(std::size_t numPoints, std::size_t numRays)
    {
      std::cout << " (" << numPoints << " known points and " << numRays << " rays).\n" << std::flush;
    }

    void DebugOutput::onBeforeOracleCall(bool forceOptimal)
    {
      std::cout << "Calling oracle" << (forceOptimal ? "" : " (w/ non-optimal sols)") << ":" << std::flush;
    }

    void DebugOutput::onAfterOracleCall(bool forceOptimal, bool feasible, std::size_t numPoints, std::size_t numRays,
        bool lastIteration)
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

    void DebugOutput::onBeforeAddPoint()
    {

    }

    void DebugOutput::onAfterAddPoint()
    {

    }

    void DebugOutput::onBeforeAddRay()
    {

    }

    void DebugOutput::onAfterAddRay()
    {

    }

    ProgressOutput::ProgressOutput(std::size_t indent)
    {
      _indent.resize(indent, ' ');
      _numStabilizationLP = 0;
      _numMainLP = 0;
      _numCache = 0;
      _numHeuristics = 0;
      _numOracles = 0;
      _lastTime = 0.0;
      _timeOverall = 0.0;
      _timeStarted = 0.0;
      _timeStabilizationLP = 0.0;
      _timeMainLP = 0.0;
      _timeCache = 0.0;
      _timeHeuristics = 0.0;
      _timeOracles = 0.0;
    }

    ProgressOutput::~ProgressOutput()
    {

    }

    void ProgressOutput::printStatistics()
    {
      std::cout << std::fixed;
      std::cout << _indent << "Number of invocations:            Heuristics: " << std::setw(6) << _numHeuristics
          << ", Oracles: " << std::setw(6) << _numOracles << ", Cache: " << std::setw(6) << _numCache << ", Stab-LP: "
          << std::setw(6) << _numStabilizationLP << ", Main-LP: " << std::setw(6) << _numMainLP << "\n";
      std::cout << _indent << "Timings (s): Overall:   " << std::setw(8) << (_timeOverall) << ", Heuristics: "
          << std::setw(6) << _timeHeuristics << ", Oracles: " << std::setw(6) << _timeOracles << ", Cache: "
          << std::setw(6) << _timeCache << ", Stab-LP: " << std::setw(6) << _timeStabilizationLP << ", Main-LP: "
          << std::setw(6) << _timeMainLP << std::endl;
    }

    double ProgressOutput::timeStamp()
    {
      double time = _timer.time();
      double elapsed = time - _lastTime;
      _lastTime = time;
      return elapsed;
    }

    void ProgressOutput::onStart()
    {
      _timer.start();
      _lastTime = _timer.time();
      _timeStarted = _lastTime;
      std::cout << _indent << "Starting to separate a " << (separatingPoint() ? "point" : "ray") << ".\n" << std::flush;
    }

    void ProgressOutput::onEnd(bool separated)
    {
      _timeOverall += _timer.time() - _timeStarted;
      _timer.stop();
      if (separated)
      {
        std::cout << _indent << "Separated a " << (separatingPoint() ? "point" : "ray") << " by "
            << (separatedEquation() ? "an equation" : (separatedFacet() ? "a facet" : "a face")) << ".\n" << std::flush;
      }
      else
      {
        std::cout << _indent << "Given " << (separatingPoint() ? "point" : "ray") << " is feasible.\n" << std::flush;
      }
    }

    void ProgressOutput::onBeforeSolve(bool stabilizing, double penalty)
    {
      std::cout << _indent << (stabilizing ? "Stab-LP " : "Main-LP ") << numRowsLP() << "x" << numColumnsLP() << ": "
          << std::flush;
      timeStamp();
    }

    void ProgressOutput::onAfterSolve(bool stabilizing, double penalty, double mainObjective, double penaltyCosts)
    {
      double elapsed = timeStamp();
      if (stabilizing)
      {
        _numStabilizationLP++;
        _timeStabilizationLP += elapsed;
      }
      else
      {
        _numMainLP++;
        _timeMainLP += elapsed;
      }

      std::cout << "solved, viol.: " << std::fixed << std::setprecision(1) << (100.0 * mainObjective) << "%";
      if (penaltyCosts > 0.0)
        std::cout << ", pen.costs: " << std::fixed << std::setprecision(1) << (100.0 * penaltyCosts) << "%).";
      else
        std::cout << ".";
      std::cout << std::scientific << std::flush;
    }

    void ProgressOutput::onPenaltyDecrease(double penalty)
    {
      std::cout << " Decreasing stabilization penalty to " << std::scientific << std::setprecision(0) << penalty
          << ".\n" << std::flush;
    }

    void ProgressOutput::onBeforeCache()
    {
      std::cout << " Cache: " << std::flush;
      timeStamp();
    }

    void ProgressOutput::onAfterCache(std::size_t numPoints, std::size_t numRays)
    {
      _numCache++;
      _timeCache += timeStamp();

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

    void ProgressOutput::onBeforeOracleCall(bool forceOptimal)
    {
      if (forceOptimal)
        std::cout << " Oracle: " << std::flush;
      else
        std::cout << " Heuristic: " << std::flush;
      timeStamp();
    }

    void ProgressOutput::onAfterOracleCall(bool forceOptimal, bool feasible, std::size_t numPoints, std::size_t numRays,
        bool lastIteration)
    {
      double elapsed = timeStamp();
      if (forceOptimal)
      {
        _numOracles++;
        _timeOracles += elapsed;
      }
      else
      {
        _numHeuristics++;
        _timeHeuristics += elapsed;
      }

      if (feasible)
      {
        if (numPoints > 0)
          std::cout << numPoints << " points.\n" << std::flush;
        else if (numRays > 0)
          std::cout << numRays << " rays.\n" << std::flush;
        else
        {
          std::cout << "feasible." << (forceOptimal && lastIteration ? "\n" : "") << std::flush;
        }
      }
      else
        std::cout << "infeasible!\n" << std::flush;
    }

  }

}

