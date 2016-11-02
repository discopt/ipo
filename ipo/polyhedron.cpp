#include "polyhedron.h"

#include "affine_hull.h"

using namespace soplex;

namespace ipo {

  Polyhedron::CollectOracle::CollectOracle(const std::shared_ptr<OracleBase>& nextOracle)
    : OracleBase("CollectOracle(" + nextOracle->name() + ")", nextOracle)
  {
    
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

      const Rational& optimum = result.points.front().objectiveValue;

      // TODO: new valid inequality.
    }

    // TODO: new points or rays.

    assert(level < heuristicLevel());

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
    std::vector<LinearConstraint> givenEquations;
    ipo::affineHull(_collectOracle, _affineHullInner, _affineHullOuter, handlers, _affineHullLastModerateHeuristic,
      _affineHullLastCheapHeuristic, givenEquations, _affineHullApproximateDirections);

    _affineHullComputed = true;
  }

  
}
