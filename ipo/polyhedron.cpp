#include "polyhedron.h"

#include "affine_hull.h"

using namespace soplex;

namespace ipo {

  Polyhedron::CollectOracle::CollectOracle(const std::shared_ptr<OracleBase>& nextOracle)
    : OracleBase("CollectOracle(" + nextOracle->name() + ")", nextOracle), _points(nextOracle->space().dimension()), 
    _rays(nextOracle->space().dimension()), _normals(nextOracle->space().dimension())
  {
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
      
      if (_normals.insert(normal))
        _inequalities.push_back(LinearConstraint('<', normal, rhs * factor));
    }

    assert(level <= heuristicLevel());

    for (std::size_t i = 0; i < result.points.size(); ++i)
      _points.insert(result.points[i].vector);
    for (std::size_t i = 0; i < result.rays.size(); ++i)
    {
      Vector ray = integralScaled(result.rays[i].vector);
      _rays.insert(ray);
    }

    return level;
  }

  HeuristicLevel Polyhedron::CollectOracle::maximizeImplementation(OracleResult& result, const VectorRational& objective,
    const ObjectiveBound& objectiveBound, HeuristicLevel minHeuristic, HeuristicLevel maxHeuristic, bool& sort, bool& checkDups)
  {
    return heuristicLevel();
  }

  Polyhedron::Polyhedron(const std::shared_ptr<OracleBase>& oracle)
    : _collectOracle(std::make_shared<CollectOracle>(oracle)), _affineHullComputed(false), _affineHullLastCheapHeuristic(1),
    _affineHullLastModerateHeuristic(0), _affineHullApproximateDirections(true)
  {

  }

  Polyhedron::~Polyhedron()
  {

  }

  void Polyhedron::affineHull()
  {
    if (_affineHullComputed)
      return;

    std::vector<AffineHullHandler*> handlers;
    DebugAffineHullHandler debugHandler(std::cout);
    handlers.push_back(&debugHandler);
    std::vector<LinearConstraint> givenEquations;
    ipo::affineHull(_collectOracle, _affineHullInner, _affineHullOuter, handlers, _affineHullLastModerateHeuristic,
      _affineHullLastCheapHeuristic, givenEquations, _affineHullApproximateDirections);

    _affineHullComputed = true;
  }

}
