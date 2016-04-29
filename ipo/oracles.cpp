#include "oracles.h"

#include <cassert>

#include "unique_rational_vectors.h"

using namespace soplex;

namespace ipo {

  SparseVectorDuplicateChecker::SparseVectorDuplicateChecker(std::size_t dimension)
  {
    _firstValues.resize(dimension, NULL);
  }


  SparseVectorDuplicateChecker::~SparseVectorDuplicateChecker()
  {

  }

  bool SparseVectorDuplicateChecker::operator()(const SVectorRational* first,
    const SVectorRational* second)
  {
    if (first->size() != second->size())
      return false;

    bool result = true;

    for (int p = first->size() - 1; p >= 0; --p)
      _firstValues[first->index(p)] = &(first->value(p));

    for (int p = second->size() - 1; p >= 0; --p)
    {
      int index = second->index(p);
      if ((_firstValues[index] == NULL) || (*_firstValues[index] != second->value(p)))
      {
        result = false;
        break;
      }
    }

    for (int p = first->size() - 1; p >= 0; --p)
      _firstValues[first->index(p)] = NULL;

    return result;
  }

  Face::Face(const Space& space) : _synced(false), _denseNormal(space.dimension())
  {

  }

  Face::Face(const Space& space, const LPRowRational& inequality) : _synced(false),
    _denseNormal(space.dimension())
  {
    add(inequality);
  }

  Face::Face(const Space& space, const LPRowSetRational& inequalities) : _synced(false),
    _denseNormal(space.dimension())
  {
    add(inequalities);
  }

  Face::~Face()
  {

  }

  void Face::add(const LPRowRational& inequality)
  {
    if (inequality.lhs() > -infinity && inequality.rhs() < infinity && inequality.lhs() < inequality.rhs())
    {
      throw std::runtime_error("Cannot add a proper ranged row as face-defining inequality.");
    }
    _inequalities.add(inequality);
    _synced = false;
  }

  void Face::add(const LPRowSetRational& inequalities)
  {
    for (std::size_t i = 0; i < inequalities.num(); ++i)
    {
      if (inequalities.lhs(i) > -infinity && inequalities.rhs(i) < infinity
          && inequalities.lhs(i) < inequalities.rhs(i))
      {
        throw std::runtime_error("Cannot add a proper ranged row as face-defining inequality.");
      }
    }

    if (_inequalities.max() < _inequalities.num() + inequalities.num())
      _inequalities.reMax(_inequalities.num() + inequalities.num());
    _inequalities.add(inequalities);
    _synced = false;
  }

  bool Face::containsPoint(const SVectorRational& point)
  {
    ensureSync();
    return _denseNormal * point == _rhs;
  }

  bool Face::containsDirection(const SVectorRational& direction)
  {
    ensureSync();
    return _denseNormal * direction == 0;
  }

  void Face::ensureSync()
  {
    if (_synced)
      return;

    _denseNormal.clear();
    _rhs = 0;

    for (std::size_t i = 0; i < _inequalities.num(); ++i)
    {
      const Rational& lhs = _inequalities.lhs(i);
      const Rational& rhs = _inequalities.rhs(i);
      if (lhs <= -infinity)
      {
        _rhs += rhs;
        const SVectorRational& row = _inequalities.rowVector(i);
        for (int p = row.size() - 1; p >= 0; --p)
          _denseNormal[row.index(p)] += row.value(p);
      }
      else
      {
        assert(rhs >= -infinity);
        _rhs -= lhs;
        const SVectorRational& row = _inequalities.rowVector(i);
        for (int p = row.size() - 1; p >= 0; --p)
          _denseNormal[row.index(p)] -= row.value(p);
      }
    }

    _sparseNormal.clear();
    _sparseNormal = _denseNormal;
    _maximumNorm = 0;
    for (int p = _sparseNormal.size() - 1; p >= 0; --p)
    {
      const Rational& x = _sparseNormal.value(p);
      if (x > 0 && x > _maximumNorm)
        _maximumNorm = x;
      if (x < 0 && x < -_maximumNorm)
        _maximumNorm = -x;
    }
  }

  OracleResult::OracleResult() : _heuristic(0)
  {

  }

  OracleResult::~OracleResult()
  {

  }

  void OracleResult::buildStart(const VectorRational& objective)
  {
    _objective = &objective;
    points.clear();
    directions.clear();
  }

  void OracleResult::buildAddPoint(DSVectorRational const* point)
  {
    Point pt = { 0, point, std::numeric_limits<std::size_t>::max() };
    points.push_back(pt);
  }

  void OracleResult::buildAddDirection(DSVectorRational const* direction)
  {
    Direction dir = { direction, std::numeric_limits<std::size_t>::max() };
    directions.push_back(dir);
  }

  void OracleResult::buildFinish(std::size_t heuristic, bool computeObjectiveValues, bool sort,
    bool removeDups)
  {
    _heuristic = heuristic;

    /// If objective is passed, recompute the objective values of all points.

    if (computeObjectiveValues)
    {
      for (std::size_t i = 0; i < points.size(); ++i)
        points[i].objectiveValue = *_objective * *(points[i].point);
    }

    /// If requested, sort all points by objective value.

    if (sort)
      std::sort(points.begin(), points.end());

    /// If requested, remove duplicates.

    if (removeDups)
      removeDuplicates(false);
  }

  void OracleResult::addToContainers(UniqueRationalVectorsBase& pointContainer,
    UniqueRationalVectorsBase& directionContainer)
  {
    for (std::size_t i = 0; i < points.size(); ++i)
    {
      if (points[i].index == std::numeric_limits<std::size_t>::max())
      {
        points[i].index = pointContainer.insertFree(points[i].point);
      }
      else
      {
        assert(pointContainer.vector(points[i].index) == points[i].point);
      }
    }
    for (std::size_t i = 0; i < directions.size(); ++i)
    {
      if (directions[i].index == std::numeric_limits<std::size_t>::max())
      {
        directions[i].index = directionContainer.insertFree(directions[i].direction);
      }
      else
      {
        assert(directionContainer.vector(directions[i].index) == directions[i].direction);
      }
    }
  }


  void OracleResult::checkConsistent()
  {
    for (std::size_t i = 1; i < points.size(); ++i)
    {
      if (points[i-1].objectiveValue < points[i].objectiveValue)
        throw std::runtime_error("Inconsistent OracleResult: Points are not sorted.");
    }

    if (removeDuplicates(true))
      throw std::runtime_error("Inconsistent OracleResult: Duplicate point or direction found.");
  }

  bool OracleResult::removeDuplicates(bool abort)
  {
    SparseVectorDuplicateChecker dup(_objective->dim());

    /// Check points.

    bool found = false;
    for (std::size_t i = 0; i < points.size(); ++i)
    {
      const Rational& firstObjective = points[i].objectiveValue;
      for (std::size_t j = i + 1; j < points.size(); ++j)
      {
        const Rational& secondObjective = points[j].objectiveValue;
        if (secondObjective != firstObjective)
          break;

        if (!dup(points[i].point, points[j].point))
          continue;

        if (abort)
          return true;

        found = true;
        delete points[j].point;
        for (std::size_t k = j + 1; k < points.size(); ++k)
          points[k-1] = points[k];
        points.pop_back();
        --j;
      }
    }

    /// Check directions.

    for (std::size_t i = 0; i < directions.size(); ++i)
    {
      for (std::size_t j = i + 1; j < directions.size(); ++j)
      {
        if (!dup(directions[i].direction, directions[j].direction))
          continue;

        if (abort)
          return true;

        found = true;
        delete directions[j].direction;
        for (std::size_t k = j + 1; k < directions.size(); ++k)
          directions[k-1] = directions[k];
        directions.pop_back();
        --j;
      }
    }

    return found;
  }

  OracleBase::OracleBase(const std::string& name, const Space& space) :
      _name(name), _space(space), _nextOracle(NULL), _thisHeuristic(0), _currentFace(NULL)
  {

  }

  OracleBase::OracleBase(const std::string& name,
    OracleBase* nextOracle) : _name(name), _space(nextOracle->space()),
    _nextOracle(nextOracle), _thisHeuristic(nextOracle->thisHeuristic() + 1), _currentFace(NULL)
  {

  }

  void OracleBase::initializedSpace()
  {
    _tempObjective.reDim(_space.dimension(), false);
  }

  OracleBase::~OracleBase()
  {

  }

  void OracleBase::setFace(Face* newFace)
  {
    if (newFace == currentFace())
      return;
    _currentFace = newFace;
    if (_nextOracle != NULL)
      _nextOracle->setFace(newFace);
  }

  void OracleBase::maximize(OracleResult& result,
    const VectorReal& objective, const ObjectiveBound& objectiveBound, std::size_t maxHeuristic,
    std::size_t minHeuristic)
  {
    if (maxHeuristic < thisHeuristic())
      return _nextOracle->maximize(result, objective, objectiveBound, maxHeuristic, minHeuristic);
    _tempObjective = objective;
    return maximize(result, _tempObjective, objectiveBound, maxHeuristic, minHeuristic);
  }

  void OracleBase::maximize(OracleResult& result,
    const SVectorRational& objective, const ObjectiveBound& objectiveBound,
    std::size_t maxHeuristic, std::size_t minHeuristic)
  {
    if (maxHeuristic < thisHeuristic())
      return _nextOracle->maximize(result, objective, objectiveBound, maxHeuristic, minHeuristic);
    _tempObjective.clear();
    _tempObjective.assign(objective);
    return maximize(result, _tempObjective, objectiveBound, maxHeuristic, minHeuristic);
  }

  void OracleBase::maximize(OracleResult& result,
    const SVectorReal& objective, const ObjectiveBound& objectiveBound, std::size_t maxHeuristic,
    std::size_t minHeuristic)
  {
    if (maxHeuristic < thisHeuristic())
      return _nextOracle->maximize(result, objective, objectiveBound, maxHeuristic, minHeuristic);
    _tempObjective.clear();
    _tempObjective.assign(objective);
    return maximize(result, _tempObjective, objectiveBound, maxHeuristic, minHeuristic);
  }

  FaceOracleBase::FaceOracleBase(const std::string& name,
    OracleBase* nextOracle, std::size_t numBlindIterations, double initialM)
    : OracleBase(name, nextOracle), _numBlindIterations(numBlindIterations),
    _M(initialM)
  {

  }

  FaceOracleBase::FaceOracleBase(const std::string& name,
    const Space& space, std::size_t numBlindIterations, double initialM)
    : OracleBase(name, space), _numBlindIterations(numBlindIterations), _M(initialM)
  {

  }

  FaceOracleBase::~FaceOracleBase()
  {

  }

  void FaceOracleBase::initializedSpace()
  {
    OracleBase::initializedSpace();

    _modifiedObjective.reDim(_space.dimension(), false);
  }

  void FaceOracleBase::setFace(Face* newFace)
  {
    if (newFace == currentFace())
      return;

    OracleBase::setFace(newFace);
    assert(newFace == currentFace());

    if (newFace != NULL)
      _factor = _M / newFace->maxNorm();
  }

  void FaceOracleBase::maximize(OracleResult& result,
    const VectorRational& objective, const ObjectiveBound& objectiveBound, std::size_t maxHeuristic,
    std::size_t minHeuristic)
  {
    assert((thisHeuristic() == 0 && _nextOracle == NULL)
      || thisHeuristic() > 0 && _nextOracle != NULL);

    if (thisHeuristic() > maxHeuristic)
      return _nextOracle->maximize(result, objective, objectiveBound, maxHeuristic, minHeuristic);

    if (currentFace() == NULL)
      return unrestrictedMaximize(result, objective, objectiveBound, objective, objectiveBound,
        maxHeuristic, minHeuristic);

    // Compute max-norm of c.

    Rational maxObjectiveCoefficient = 0;
    for (int v = 0; v < space().dimension(); ++v)
    {
      if (objective[v] >= 0)
      {
        if (objective[v] > maxObjectiveCoefficient)
          maxObjectiveCoefficient = objective[v];
      }
      else
      {
        if (-objective[v] > maxObjectiveCoefficient)
          maxObjectiveCoefficient = -objective[v];
      }
    }

    // Zero vector does not need change.

    if (maxObjectiveCoefficient == 0)
    {
      return unrestrictedMaximize(result, objective, objectiveBound, objective, objectiveBound,
        maxHeuristic, minHeuristic);
    }

    // Copy objective since modification might only affect a subset of the coordinates.

    for (std::size_t v = 0; v < space().dimension(); ++v)
    {
      _modifiedObjective[v] = objective[v];
    }
    const VectorRational& denseNormal = currentFace()->denseNormal();
    const SVectorRational& sparseNormal = currentFace()->sparseNormal();

    // Loop that increases _M as long as results are invalid for current face.
    OracleResult unrestrictedResult;
    bool verifiedProperFace = false;
    ObjectiveBound modifiedObjectiveBound;
    Rational modifiedObjectiveShift;
    for (std::size_t iteration = 0; true; ++iteration)
    {
      /// So far no point in the face was found, so we check whether it is empty.

      if (iteration == _numBlindIterations)
      {
        /// Setup objective via face's normal vector.

        for (std::size_t v = 0; v < space().dimension(); ++v)
          _modifiedObjective[v] = 0;
        for (int p = sparseNormal.size() - 1; p >= 0; --p)
          _modifiedObjective[sparseNormal.index(p)] = sparseNormal.value(p);
        modifiedObjectiveBound.value = currentFace()->rhs();
        modifiedObjectiveBound.strict = false;

        unrestrictedMaximize(unrestrictedResult, _modifiedObjective, modifiedObjectiveBound,
          _modifiedObjective, modifiedObjectiveBound, maxHeuristic, minHeuristic);

        /*
          * The following happens after maximizing F's normal:
          *
          * infeasible: Return infeasible which (heuristically) proves that F is empty.
          * unbounded: Raise an error.
          * feasible: Continue as usual.
          */

        if (unrestrictedResult.isInfeasible())
        {
          result.buildStart(objective);
          return result.buildFinish(unrestrictedResult.heuristic(), false, false, false);
        }
        else if (unrestrictedResult.isUnbounded())
        {
          for (std::size_t i = 0; i < unrestrictedResult.directions.size(); ++i)
            delete unrestrictedResult.directions[i].direction;
          throw std::runtime_error(
            "Error in FaceOracleBase::maximize: Maximizing F's normal is unbounded.");
        }

        continue;
      }

      Rational scaling = _factor / maxObjectiveCoefficient;
      modifiedObjectiveShift = scaling * currentFace()->rhs();
      modifiedObjectiveBound.value = objectiveBound.value + modifiedObjectiveShift;
      modifiedObjectiveBound.strict = objectiveBound.strict;
      for (int p = sparseNormal.size() - 1; p >= 0; --p)
      {
        int v = sparseNormal.index(p);
        _modifiedObjective[v] = objective[v] + scaling * sparseNormal.value(p);
      }

      unrestrictedMaximize(unrestrictedResult, _modifiedObjective, modifiedObjectiveBound,
        objective, objectiveBound, maxHeuristic, minHeuristic);

      if (unrestrictedResult.isFeasible())
      {
        // Remove points that do not lie in F.

        result.buildStart(objective);
        for (std::size_t i = 0; i < unrestrictedResult.points.size(); ++i)
        {
          if (currentFace()->containsPoint(*unrestrictedResult.points[i].point))
            result.buildAddPoint(unrestrictedResult.points[i].point);
          else
            delete unrestrictedResult.points[i].point;
        }

        if (!result.points.empty())
        {
          // If best point satisfied objective bound, return everything, otherwise increase M.

          Rational bestValue = result.points.front().objectiveValue - modifiedObjectiveShift;
          if (objectiveBound.satisfiedBy(bestValue))
            return result.buildFinish(unrestrictedResult.heuristic(), true, false, false);
        }
      }
      else if (unrestrictedResult.isUnbounded())
      {
        // Remove directions that do not lie in F.

        result.buildStart(objective);
        for (std::size_t i = 0; i < unrestrictedResult.directions.size(); ++i)
        {
          if (currentFace()->containsDirection(*unrestrictedResult.directions[i].direction))
            result.buildAddDirection(unrestrictedResult.directions[i].direction);
          else
            delete unrestrictedResult.directions[i].direction;
        }

        // If a direction remains, return everything, otherwise increase M.

        if (!result.directions.empty())
          return result.buildFinish(unrestrictedResult.heuristic(), false, false, false);
      }
      else
      {
        assert(unrestrictedResult.isInfeasible());
        result.buildStart(objective);
        return result.buildFinish(unrestrictedResult.heuristic(), false, false, false);
      }

      _factor *= 2; // Increase current scaling.
      _M *= 2; // Increase scaling for future calls.
    }
  }

} /* namespace ipo */
