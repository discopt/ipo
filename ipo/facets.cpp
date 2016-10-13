#include "facets.h"

#include <vector>
#include <algorithm>
#include <map>

#include "polar_lp.h"

using namespace soplex;

namespace ipo {

  FacetSeparationState::FacetSeparationState()
  {

  }

  FacetSeparationState::~FacetSeparationState()
  {

  }

  class FacetSeparation: public FacetSeparationState, PolarLPHandler
  {
  public:
    FacetSeparation(std::vector<FacetSeparationHandler*>& handlers, const std::shared_ptr<OracleBase>& oracle,
      const Vector& targetVector, bool separatingRay)
      : _handlers(handlers), _oracle(oracle), _targetVector(targetVector), _separatingRay(separatingRay),
      _approximateLP(oracle, *this, true), _exactLP(oracle, *this, false), _approximateSolve(false), _exactSolve(false),
      _separatedByFacet(false), _separatedByEquation(false)
    {

    }

    virtual ~FacetSeparation()
    {

    }

    virtual const Space& oracleSpace() const
    {
      return _oracle->space();
    }

    virtual const Space& polarSpace() const
    {
      return _approximateLP.polarSpace();
    }

    virtual double polarTolerance() const
    {
      return _approximateSolve ? _approximateLP.getTolerance() : _exactLP.getTolerance();
    }

    virtual std::size_t polarNumPoints() const
    {
      return _approximateSolve ? _approximateLP.numPointsLP() : _exactLP.numPointsLP();
    }

    virtual std::size_t polarNumRays() const
    {
      return _approximateSolve ? _approximateLP.numRaysLP() : _exactLP.numRaysLP();
    }

    virtual std::size_t polarNumPointsAdded() const
    {
      return _approximateSolve ? _approximateLP.numPointsAdded() : _exactLP.numPointsAdded();
    }

    virtual std::size_t polarNumRaysAdded() const
    {
      return _approximateSolve ? _approximateLP.numRaysAdded() : _exactLP.numRaysAdded();
    }

    virtual std::size_t polarNumRowsLP() const
    {
      return _approximateSolve ? _approximateLP.numRowsLP() : _exactLP.numRowsLP();
    }

    virtual std::size_t polarNumColumnsLP() const
    {
      return _approximateSolve ? _approximateLP.numColumnsLP() : _exactLP.numColumnsLP();
    }

    virtual std::size_t polarNumNonzerosLP() const
    {
      return _approximateSolve ? _approximateLP.numNonzerosLP() : _exactLP.numNonzerosLP();
    }

    virtual LinearConstraint currentInequality() const
    {
      return _approximateSolve ? _approximateLP.currentInequality() : _exactLP.currentInequality();
    }

    virtual soplex::Rational oracleObjectiveValue() const
    {
      return _approximateSolve ? _approximateLP.oracleObjectiveValue() : _exactLP.oracleObjectiveValue();
    }

    virtual bool approximateSolve() const
    {
      return _approximateSolve;
    }

    virtual bool exactSolve() const
    {
      return _exactSolve;
    }

    virtual HeuristicLevel oracleMaxHeuristicLevel() const
    {
      return _approximateSolve ? _approximateLP.oracleMaxHeuristicLevel() : _exactLP.oracleMaxHeuristicLevel();
    }

    virtual HeuristicLevel oracleMinHeuristicLevel() const
    {
      return _approximateSolve ? _approximateLP.oracleMinHeuristicLevel() : _exactLP.oracleMinHeuristicLevel();
    }

    virtual HeuristicLevel oracleResultHeuristicLevel() const
    {
      return _approximateSolve ? _approximateLP.oracleResultHeuristicLevel() : _exactLP.oracleResultHeuristicLevel();
    }

    virtual std::size_t oracleNumPoints() const
    {
      return _approximateSolve ? _approximateLP.oracleNumPoints() : _exactLP.oracleNumPoints();
    }

    virtual std::size_t oracleNumRays() const
    {
      return _approximateSolve ? _approximateLP.oracleNumRays() : _exactLP.oracleNumRays();
    }

    virtual Vector separationTarget() const
    {
      return _targetVector;
    }

    virtual bool separatingRay() const
    {
      return _separatingRay;
    }

    virtual bool separatedByFacet() const
    {
      return _separatedByFacet;
    }

    virtual bool separatedByEquation() const
    {
      return _separatedByEquation;
    }

    virtual void notify(PolarLPHandler::Event event, XPolarLP& polarLP)
    {
      assert(static_cast<FacetSeparationHandler::Event>(PolarLPHandler::LP_BEGIN) == FacetSeparationHandler::LP_BEGIN);
      assert(static_cast<FacetSeparationHandler::Event>(PolarLPHandler::LP_END) == FacetSeparationHandler::LP_END);
      assert(static_cast<FacetSeparationHandler::Event>(PolarLPHandler::ORACLE_BEGIN) == FacetSeparationHandler::ORACLE_BEGIN);
      assert(static_cast<FacetSeparationHandler::Event>(PolarLPHandler::ORACLE_END) == FacetSeparationHandler::ORACLE_END);
      assert(static_cast<FacetSeparationHandler::Event>(PolarLPHandler::POINTS_BEGIN) == FacetSeparationHandler::POINTS_BEGIN);
      assert(static_cast<FacetSeparationHandler::Event>(PolarLPHandler::POINTS_END) == FacetSeparationHandler::POINTS_END);
      assert(static_cast<FacetSeparationHandler::Event>(PolarLPHandler::RAYS_BEGIN) == FacetSeparationHandler::RAYS_BEGIN);
      assert(static_cast<FacetSeparationHandler::Event>(PolarLPHandler::RAYS_END) == FacetSeparationHandler::RAYS_END);
      switch (event)
      {
        case LP_BEGIN:
        case LP_END:
        case ORACLE_BEGIN:
        case ORACLE_END:
        case POINTS_BEGIN:
        case POINTS_END:
        case RAYS_BEGIN:
        case RAYS_END:
          notify(static_cast<FacetSeparationHandler::Event>(event));
        break;
      }
    }

  protected:

    void notify(FacetSeparationHandler::Event event)
    {
      for (std::size_t i = 0; i < _handlers.size(); ++i)
        _handlers[i]->notify(event, *this);
    }

  public:

    bool run(const InnerDescription& spanning, LinearConstraint& constraint, InnerDescription* certificate)
    {
      _separatedByEquation = false;
      _separatedByFacet = false;
      notify(FacetSeparationHandler::BEGIN);

      DVectorRational dense;
      dense.reDim(polarSpace().dimension());
      _approximateLP.clear();
      _exactLP.clear();

      if (separatingRay())
      {
        assert(false);
      }
      else
      {
        // Normalization constraint.

        for (std::size_t i = 0; i < spanning.points.size(); ++i)
        {
          dense -= spanning.points[i];
          _approximateLP.addPointRow(spanning.points[i], false);
          _exactLP.addPointRow(spanning.points[i], false);
        }
        for (std::size_t i = 0; i < spanning.rays.size(); ++i)
        {
          dense -= spanning.rays[i];
          _approximateLP.addPointRow(spanning.rays[i], false);
          _exactLP.addPointRow(spanning.rays[i], false);
        }
        dense /= int(spanning.points.size() + spanning.rays.size());

        for (std::size_t p = 0; p < _targetVector.size(); ++p)
          dense[_targetVector.index(p)] += _targetVector.value(p);

        DSVectorRational sparse(dense);
        _approximateLP.addRow(-soplex::infinity, sparse, Rational(1), false);
        _exactLP.addRow(-soplex::infinity, sparse, Rational(1), false);

        // Objective.

        dense.clear();
        for (std::size_t p = 0; p < _targetVector.size(); ++p)
          dense[_targetVector.index(p)] = _targetVector.value(p);
        dense[polarSpace().dimension() - 1] = -1;
        _approximateLP.setObjective(dense);
        _exactLP.setObjective(dense);
      }

      notify(FacetSeparationHandler::INITIALIZED);

      // Solve approximate polar LP.

      _approximateSolve = true;
      notify(FacetSeparationHandler::APPROXIMATE_SOLVE_BEGIN);
      _approximateLP.solve();
      notify(FacetSeparationHandler::APPROXIMATE_SOLVE_END);
      _approximateSolve = false;

      // Copy relevant rows from approximate to exact LP.

      InnerDescription pointsRays;
      _approximateLP.getTightPointsRays(pointsRays, true);
      for (std::size_t i = 0; i < pointsRays.points.size(); ++i)
        _exactLP.addPointRow(pointsRays.points[i], true);
      for (std::size_t i = 0; i < pointsRays.rays.size(); ++i)
        _exactLP.addRayRow(pointsRays.rays[i], true);

      // Solve exact polar LP.

      _exactSolve = true;
      notify(FacetSeparationHandler::EXACT_SOLVE_BEGIN);
      _exactLP.solve();
      notify(FacetSeparationHandler::EXACT_SOLVE_END);
      _exactSolve = false;

      if (_exactLP.getObjectiveValue() == 0)
      {
        notify(FacetSeparationHandler::END);
        return false;
      }

      // Extract tight points and rays.

      _exactLP.getTightPointsRays(pointsRays);
      if (certificate != NULL)
        *certificate = pointsRays;

      constraint = _exactLP.currentInequality();
      std::size_t resultDimMinus1 = pointsRays.points.size() + pointsRays.rays.size();
      std::size_t givenDimMinus1 = spanning.points.size() + spanning.rays.size();
      if (resultDimMinus1 == givenDimMinus1)
      {
        _separatedByEquation = true;
        constraint = LinearConstraint('=', constraint.normal(), constraint.rhs());
      }
      else
      {
        _separatedByFacet = true;
      }

      notify(FacetSeparationHandler::END);

      return true;
    }

  public:
    std::size_t paramMaxAge;

  protected:
    std::vector<FacetSeparationHandler*>& _handlers;
    std::shared_ptr<OracleBase> _oracle;
    const Vector& _targetVector;
    bool _separatingRay;
    XPolarLP _approximateLP;
    XPolarLP _exactLP;

    bool _approximateSolve;
    bool _exactSolve;
    bool _separatedByFacet;
    bool _separatedByEquation;
  };

  FacetSeparationHandler::FacetSeparationHandler()
  {

  }

  FacetSeparationHandler::~FacetSeparationHandler()
  {

  }

  DebugFacetSeparationHandler::DebugFacetSeparationHandler(std::ostream& stream, bool printPointsAndRays, bool printInequalities)
    : _stream(stream), _printPointsAndRays(printPointsAndRays), _printInequalities(printInequalities)
  {

  }

  DebugFacetSeparationHandler::~DebugFacetSeparationHandler()
  {

  }

  void DebugFacetSeparationHandler::notify(FacetSeparationHandler::Event event, FacetSeparationState& state)
  {
    switch (event)
    {
      case  BEGIN:
        _stream << "FS: Separating " << (state.separatingRay() ? "ray" : "point");
        if (_printPointsAndRays)
        {
          _stream << " ";
          state.oracleSpace().printVector(_stream, state.separationTarget());
        }
        _stream << ".\n";
      break;
      case INITIALIZED:
        _stream << "FS: Initialized LPs.\n";
      break;
      case APPROXIMATE_SOLVE_BEGIN:
        _stream << "FS: Entering approximate mode.\n";
      break;
      case APPROXIMATE_SOLVE_END:
        _stream << "FS: Leaving approximate mode.\n";
      break;
      case EXACT_SOLVE_BEGIN:
        _stream << "FS: Entering exact mode.\n";
      break;
      case EXACT_SOLVE_END:
        _stream << "FS: Leaving exact mode.\n";
      break;
      case END:
        if (state.separatedByFacet())
          _stream << "FS: Separated " << (state.separatingRay() ? "ray" : "point") << " by a facet.\n";
        else if (state.separatedByEquation())
          _stream << "FS: Separated " << (state.separatingRay() ? "ray" : "point") << " by an equation.\n";
        else if (state.separatingRay())
          _stream << "FS: Given ray lies in recession cone.\n";
        else
          _stream << "FS: Given point lies in polyhedron.\n";
      break;
      case LP_BEGIN:
        _stream << "FS: Solving polar LP of size " << state.polarNumRowsLP() << "x" << state.polarNumColumnsLP()
          << ", #nonzeros: " << state.polarNumNonzerosLP() << ", #points: " << state.polarNumPoints() << ", #rays: "
          << state.polarNumRays() << ", tolerance: " << state.polarTolerance() << ".\n";
      break;
      case LP_END:
        _stream << "FS: Solved polar LP.\n";
        if (_printInequalities)
        {
          _stream << "FS: Current inequality: ";
          state.oracleSpace().printLinearConstraint(_stream, state.currentInequality());
          _stream << "\n";
        }
      break;
      case ORACLE_BEGIN:
        _stream << "FS: Maximizing current normal vector";
        if (state.oracleMaxHeuristicLevel() == std::numeric_limits<HeuristicLevel>::max())
          _stream << " (heurLevel >= " << state.oracleMinHeuristicLevel() << ")\n";
        else
          _stream << " (" << state.oracleMaxHeuristicLevel() << " >= heurLevel >= " << state.oracleMinHeuristicLevel() << ")\n";
      break;
      case ORACLE_END:
        if (state.oracleNumPoints() > 0)
        {
          _stream << "FS: Oracle returned " << state.oracleNumPoints() << " points";
          if (_printInequalities)
            _stream << ", the best one having rhs " << state.oracleObjectiveValue();
        }
        else if (state.oracleNumRays() > 0)
          _stream << "FS: Oracle returned " << state.oracleNumRays() << " rays";
        else
          _stream << "FS: Oracle claimed infeasible";
        _stream << ", heurLevel = " << state.oracleResultHeuristicLevel() << "\n";
      break;
      case POINTS_BEGIN:
        _stream << "FS: Adding points.\n";
      break;
      case POINTS_END:
        _stream << "FS: Added " << state.polarNumPointsAdded() << " points.\n";
      break;
      case RAYS_BEGIN:
        _stream << "FS: Adding rays.\n";
        if (_printPointsAndRays)
        {
//           _stream << "FS: ";
//           state.space().printVector(_stream, state.innerDescription().rays.back()); // TODO: ???
//           _stream << "\n";
        }
      break;
      case RAYS_END:
        _stream << "FS: Added " << state.polarNumRaysAdded() << " rays.\n";
      break;
      default:
        _stream << "FS: Unhandled event " << event << ".\n" << std::endl;
      break;
    }
    _stream << std::flush;
  }

  StatisticsFacetSeparationHandler::StatisticsFacetSeparationHandler()
  {

  }

  StatisticsFacetSeparationHandler::~StatisticsFacetSeparationHandler()
  {

  }

  void StatisticsFacetSeparationHandler::reset()
  {

  }

  void StatisticsFacetSeparationHandler::notify(FacetSeparationHandler::Event event, FacetSeparationState& state)
  {

  }

  bool separatePoint(const std::shared_ptr<OracleBase>& oracle, const Vector& point, const InnerDescription& spanning,
    std::vector<FacetSeparationHandler*>& handlers, LinearConstraint& constraint, InnerDescription* certificate)
  {
    FacetSeparation algorithm(handlers, oracle, point, false);
    algorithm.paramMaxAge = 10;
    return algorithm.run(spanning, constraint, certificate);
  }

  bool separateRay(const std::shared_ptr<OracleBase>& oracle, const Vector& ray, const InnerDescription& spanning,
    std::vector<FacetSeparationHandler*>& handlers, LinearConstraint& constraint, InnerDescription* certificate)
  {
    FacetSeparation algorithm(handlers, oracle, ray, true);
    algorithm.paramMaxAge = 10;
    return algorithm.run(spanning, constraint, certificate);
  }


//   LinearConstraint separatePoint(const std::shared_ptr<OracleBase>& oracle, const Vector& point,
//     const InnerDescription& spanning, std::vector<FacetSeparationHandler*>& handlers, InnerDescription* certificate)
//   {
//     FacetSeparation algorithm(handlers, oracle, point, false);
//     algorithm.paramMaxAge = 10;
//     return algorithm.run(spanning, certificate);
//   }
//
//   LinearConstraint separateRay(const std::shared_ptr<OracleBase>& oracle, const Vector& ray,
//     const InnerDescription& spanning, std::vector<FacetSeparationHandler*>& handlers, InnerDescription* certificate)
//   {
//     FacetSeparation algorithm(handlers, oracle, ray, true);
//     algorithm.paramMaxAge = 10;
//     return algorithm.run(spanning, certificate);
//   }




  namespace Separation {

    class Implementation: public PolarLP
    {
    public:
      Implementation(const std::vector<Vector>& spanningPoints, const std::vector<Vector>& spanningRays,
        const std::vector<std::size_t>& columnBasis, const std::shared_ptr<OracleBase>& oracle)
        : PolarLP(oracle, 16.0, 30), _separatingPoint(false), _output(NULL), _separatedEquation(false), _separatedFacet(false)
      {
        SVectorRational vector;
        _normalizationConstraint = addConstraint(-infinity, vector, Rational(0));
        _dimension = int(columnBasis.size()) - 1;

        /// Go through spanning points.

        _interiorPoint.reDim(n(), true);
        for (std::size_t i = 0; i < spanningPoints.size(); ++i)
        {
          addPointConstraint(spanningPoints[i]);
          _interiorPoint += spanningPoints[i];
        }
        for (std::size_t v = 0; v < n(); ++v)
          _interiorPoint[v] /= int(spanningPoints.size());

        /// Go through spanning rays.

        _interiorRay.reDim(n(), true);
        for (std::size_t i = 0; i < spanningRays.size(); ++i)
        {
          addRayConstraint(spanningRays[i]);
          _interiorRay += _interiorRay, spanningRays[i];
        }
        _interiorPoint += _interiorRay;

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

      bool run(const Vector& target, bool separatingPoint, OutputBase& output)
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
        for (std::size_t p = 0; p < target.size(); ++p)
          dense[target.index(p)] = target.value(p);
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

        // Perform stabilization.

        stabilizedPresolve();

        // Reset basis to the one with maximum size basis matrix.

        setBasis(_basis);

        // Optimize

        optimize(true);

        /// Extract solution

        dense.reDim(n() + 1);
        getPrimalSolution(dense);
        Rational rhs = 0;
        std::swap(rhs, dense[n()]);
        _violation = getObjectiveValue();

        /// Use basis information for the certificate.

        Basis basis;
        getBasis(basis);
        _certificate.points.clear();
        std::copy(basis.tightPoints.begin(), basis.tightPoints.end(), std::back_inserter(_certificate.points));
        _certificate.rays.clear();
        std::copy(basis.tightRays.begin(), basis.tightRays.end(), std::back_inserter(_certificate.rays));

        /// Based on certificate we know what we separated.

        if (_certificate.points.size() + _certificate.rays.size() == _dimension + 1)
        {
          _separatedEquation = true;
          _inequality = LinearConstraint('=', denseToVector(dense) , rhs);
        }
        else
        {
          _inequality = LinearConstraint('<', denseToVector(dense) , rhs);
          if (_certificate.points.size() + _certificate.rays.size() == _dimension)
            _separatedFacet = true;
          else
          {
            std::cerr << "\n!!! IPO computed lower-dimensional face.";
            std::cerr << "This is a BUG, probably caused by a SoPlex restart !!!" << std::endl;
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
        return _oracle->space().dimension();
      }

    protected:
      bool _separatingPoint;
      OutputBase* _output;
      int _dimension;
      soplex::DVectorRational _interiorPoint;
      soplex::DVectorRational _interiorRay;
      std::size_t _normalizationConstraint;
      Basis _basis;

      /// Result.

      LinearConstraint _inequality;
      Certificate _certificate;
      Rational _violation;
      bool _separatedFacet;
      bool _separatedEquation;
    };

    Result::Result(const std::vector<Vector>& spanningPoints, const std::vector<Vector>& spanningRays,
      const std::vector<std::size_t>& columnBasis, const std::shared_ptr<OracleBase>& oracle)
    {
      _implementation = new Implementation(spanningPoints, spanningRays, columnBasis, oracle);
    }

    Result::~Result()
    {
      assert(hasImplementation());
      delete _implementation;
    }

    LinearConstraint Result::inequality() const
    {
      return _implementation->_inequality;
    }

    void Result::certificate(Certificate& certificate) const
    {
      certificate = _implementation->_certificate;
    }

    const Rational& Result::violation() const
    {
      return _implementation->_violation;
    }

    bool Result::separatePoint(const Vector& targetPoint, OutputBase& output)
    {
      return _implementation->run(targetPoint, true, output);
    }

    bool Result::separateRay(const Vector& targetRay, OutputBase& output)
    {
      return _implementation->run(targetRay, false, output);
    }

    InformationBase::InformationBase()
      : _implementation(NULL)
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

