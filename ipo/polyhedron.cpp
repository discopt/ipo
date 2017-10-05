#include "polyhedron.h"

#include "affine_hull.h"

// Uncomment the following line for debugging.
// #define DEBUG

using namespace soplex;

namespace ipo {

  Polyhedron::VectorInfo::VectorInfo(Vector& vector, bool isPoint)
    : _vector(vector), _isPoint(isPoint)
  {

  }

  Polyhedron::VectorInfo::~VectorInfo()
  {

  }

  Polyhedron::Face::Face(const LinearConstraint& inequality)
    : _inequality(inequality), _hasDimension(false)
  {

  }

  Polyhedron::Face::~Face()
  {

  }

  Polyhedron::CollectOracle::CollectOracle(const std::shared_ptr<OracleBase>& nextOracle)
    : OracleBase("CollectOracle(" + nextOracle->name() + ")", nextOracle), _points(nextOracle->space().dimension()),
    _rays(nextOracle->space().dimension()), _inequalities(nextOracle->space().dimension())
  {
    LinearConstraint constraint = completeFaceConstraint();
    Vector normal = constraint.normal();
    _inequalities.insert(normal, std::make_shared<Face>(constraint));

    _heuristicLevel--; // Effectively set heuristicLevel to that of nextOracle.

    initializeSpace(nextOracle->space());
  }

  Polyhedron::CollectOracle::~CollectOracle()
  {

  }

  void Polyhedron::CollectOracle::setFace(const LinearConstraint& newFace)
  {
    OracleBase::setFace(newFace);
  }

  HeuristicLevel Polyhedron::CollectOracle::maximizeController(OracleResult& result, const VectorRational& objective,
    const ObjectiveBound& objectiveBound, HeuristicLevel minHeuristic, HeuristicLevel maxHeuristic, bool& sort, bool& checkDups)
  {
    HeuristicLevel level = OracleBase::maximizeController(result, objective, objectiveBound, minHeuristic, maxHeuristic, sort,
      checkDups);

    if (level == 0 && result.isFeasible())
    {
      if (sort)
      {
        result.computeMissingObjectiveValues();
        result.sortPoints();
        sort = false;
      }

      Vector normal = denseToVector(objective);
      const Rational& rhs = result.points.front().objectiveValue;
      Rational factor;
      scaleIntegral(normal, &factor);
      _inequalities.insert(normal, std::make_shared<Face>(LinearConstraint('<', normal, rhs * factor)));
    }

    assert(level <= heuristicLevel());

    for (std::size_t i = 0; i < result.points.size(); ++i)
      _points.insert(result.points[i].vector, VectorInfo(result.points[i].vector, true));
    for (std::size_t i = 0; i < result.rays.size(); ++i)
    {
      Vector ray = integralScaled(result.rays[i].vector);
      _rays.insert(ray, VectorInfo(result.rays[i].vector, false));
    }

    return level;
  }

  HeuristicLevel Polyhedron::CollectOracle::maximizeImplementation(OracleResult& result, const VectorRational& objective,
    const ObjectiveBound& objectiveBound, HeuristicLevel minHeuristic, HeuristicLevel maxHeuristic, bool& sort, bool& checkDups)
  {
    return heuristicLevel();
  }

  Polyhedron::Polyhedron(const std::shared_ptr<OracleBase>& oracle)
    : _collectOracle(std::make_shared<CollectOracle>(oracle)), _completeFace(_collectOracle->_inequalities[zeroVector()]),
    _affineHullLastCheapHeuristic(1), _affineHullLastModerateHeuristic(0), _affineHullApproximateDirections(true), 
    _affineHullExactDirectionTimeLimit(std::numeric_limits<double>::max())
  {

  }
  
  Polyhedron::Polyhedron(const std::shared_ptr<DefaultOracleWrapper>& oracleWrapper)
    : _collectOracle(std::make_shared<CollectOracle>(oracleWrapper->queryOracle())),
    _completeFace(_collectOracle->_inequalities[zeroVector()]),
    _affineHullLastCheapHeuristic(1), _affineHullLastModerateHeuristic(0), _affineHullApproximateDirections(true),
    _affineHullExactDirectionTimeLimit(std::numeric_limits<double>::max())
  {

  }


  Polyhedron::~Polyhedron()
  {

  }

  void Polyhedron::affineHull(std::shared_ptr<Face>& face, std::vector<AffineHullHandler*>& extraHandlers,
    const std::vector<LinearConstraint>& givenEquations)
  {
    if (face->hasDimension())
      return;

    LinearConstraint faceConstraint = LinearConstraint('=', face->inequality().normal(), face->inequality().rhs());
    _collectOracle->setFace(faceConstraint);

    std::vector<AffineHullHandler*> handlers;
    std::copy(extraHandlers.begin(), extraHandlers.end(), std::back_inserter(handlers));
#ifdef DEBUG
    DebugAffineHullHandler debugHandler(std::cout);
    handlers.push_back(&debugHandler);
#endif

    std::vector<LinearConstraint> initialEquations;
    std::copy(givenEquations.begin(), givenEquations.end(), std::back_inserter(initialEquations));
    if (!faceConstraint.definesCompleteFace())
    {
      std::copy(_completeFace->_outerDescription.begin(), _completeFace->_outerDescription.end(),
        std::back_inserter(initialEquations));
    }

    ipo::affineHull(_collectOracle, face->_innerDescription, face->_outerDescription, handlers, _affineHullLastCheapHeuristic,
      _affineHullLastModerateHeuristic, initialEquations, _affineHullApproximateDirections, _affineHullExactDirectionTimeLimit);
    face->_hasDimension = true;
  }

  void Polyhedron::addConstraint(const LinearConstraint& constraint, bool doNormalize)
  {
    LinearConstraint normalized = normalize(constraint);
    Vector normal = normalized.normal();
    _collectOracle->_inequalities.insert(normal, std::make_shared<Face>(doNormalize ? normalized : constraint));
  }

  void Polyhedron::getFaces(std::vector<std::shared_ptr<Face> >& constraints, bool onlyInequalities, bool onlyWithDimension)
  {
    constraints.clear();
    for (VectorMap<std::shared_ptr<Face> >::Iterator iter = _collectOracle->_inequalities.begin();
      iter != _collectOracle->_inequalities.end(); ++iter)
    {
      std::shared_ptr<Face> face = *iter;
      if (onlyInequalities && face->inequality().type() == '=')
        continue;
      if (onlyWithDimension && !face->hasDimension())
        continue;
      constraints.push_back(face);
    }
  }

  bool Polyhedron::separatePoint(const Vector& point, std::vector<FacetSeparationHandler*>& extraHandlers, LinearConstraint& constraint,
    InnerDescription* certificate)
  {
    // We need a set of spanning points and rays.

    affineHull();
    
    std::vector<FacetSeparationHandler*> handlers;
    std::copy(extraHandlers.begin(), extraHandlers.end(), std::back_inserter(handlers));
#ifdef DEBUG
    DebugFacetSeparationHandler debugHandler(std::cout);
    handlers.push_back(&debugHandler);
#endif

    return ipo::separatePoint(_collectOracle, point, _completeFace->_innerDescription, handlers, constraint, certificate);

    // TODO: Add facet with certificates to list of faces.
  }

  bool Polyhedron::separateRay(const Vector& ray, std::vector<FacetSeparationHandler*>& extraHandlers, LinearConstraint& constraint,
    InnerDescription* certificate)
  {
    // We need a set of spanning points and rays.

    affineHull();
    
    std::vector<FacetSeparationHandler*> handlers;
    std::copy(extraHandlers.begin(), extraHandlers.end(), std::back_inserter(handlers));
#ifdef DEBUG
    DebugFacetSeparationHandler debugHandler(std::cout);
    handlers.push_back(&debugHandler);
#endif

    return ipo::separateRay(_collectOracle, ray, _completeFace->_innerDescription, handlers, constraint, certificate);

    // TODO: Add facet with certificates to list of faces.
  }

}
