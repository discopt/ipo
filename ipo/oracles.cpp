#include "oracles.h"

#include <cassert>

#include "unique_vectors.h"

using namespace soplex;

namespace ipo {

  OracleResult::Point::Point(Vector& vec)
    : vector(vec)
  {
    this->objectiveValue = -infinity;
  }

  OracleResult::Point::Point(Vector& vec, const Rational& value)
    : vector(vec)
  {
    this->objectiveValue = value;
  }

  OracleResult::Ray::Ray(Vector& vec)
    : vector(vec)
  {

  }

  OracleResult::OracleResult() : _heuristicLevel(0)
  {

  }

  OracleResult::~OracleResult()
  {

  }

  void OracleResult::checkConsistent()
  {
    for (std::size_t i = 1; i < points.size(); ++i)
    {
      if (points[i-1].objectiveValue < points[i].objectiveValue)
        throw std::runtime_error("Inconsistent OracleResult: Points are not sorted.");
    }
  }

  void OracleResult::computeMissingObjectiveValues()
  {
    for (std::size_t i = 0; i < points.size(); ++i)
    {
      if (points[i].objectiveValue == -infinity)
      {
        points[i].objectiveValue = *_objective * points[i].vector;
      }
    }
  }

  void OracleResult::sortPoints()
  {
    std::sort(points.begin(), points.end());
  }

  void OracleResult::removeDuplicates()
  {
    // Check points.

    UniqueVectors usv(_objective->dim());
    std::size_t write = 0;
    for (std::size_t read = 0; read < points.size(); ++read)
    {
      if (usv.insert(points[read].vector))
      {
        if (write != read)
          points[write] = points[read];
        ++write;
      }
    }
    while (points.size() > write)
      points.pop_back();

    // Check rays.

    usv.clear();
    write = 0;
    for (std::size_t read = 0; read < rays.size(); ++read)
    {
      if (usv.insert(rays[read].vector))
      {
        if (write != read)
          rays[write] = rays[read];
        ++write;
      }
    }
    while (rays.size() > write)
      rays.pop_back();
  }



  OracleBase::OracleBase(const std::string& name, const std::shared_ptr<OracleBase>& nextOracle)
    : _name(name), _nextOracle(nextOracle), _heuristicLevel(nextOracle ? nextOracle->heuristicLevel() + 1 : 0)
  {

  }

  void OracleBase::initializeSpace(const Space& space)
  {
    _space = space;

    _tempObjective.reDim(_space.dimension(), false);
  }

  OracleBase::~OracleBase()
  {

  }

  void OracleBase::setFace(const LinearConstraint& newFace)
  {
    if (newFace == currentFace())
      return;
    _currentFace = newFace;
    if (_nextOracle != NULL)
      _nextOracle->setFace(newFace);
  }

  void OracleBase::maximize(OracleResult& result, const Vector& objective, const ObjectiveBound& objectiveBound,
    HeuristicLevel minHeuristic, HeuristicLevel maxHeuristic)
  {
    vectorToDense(objective, _tempObjective);
    return maximize(result, _tempObjective, objectiveBound, minHeuristic, maxHeuristic);
  }

  void OracleBase::maximize(OracleResult& result, const soplex::VectorRational& objective, const ObjectiveBound& objectiveBound,
    HeuristicLevel minHeuristic, HeuristicLevel maxHeuristic)
  {
    assert(objective.dim() == space().dimension());

    // Initialize result.

    result._objective = &objective;
    result.points.clear();
    result.rays.clear();

    bool sort = false;
    bool checkDuplicates = false;

    result._heuristicLevel = maximizeController(result, objective, objectiveBound, minHeuristic, maxHeuristic, sort,
      checkDuplicates);

    // Compute missing objective values.
    result.computeMissingObjectiveValues();

    // If requested, sort all points by objective value.

    if (sort)
      std::sort(result.points.begin(), result.points.end());

    // If requested, remove duplicates.

    if (checkDuplicates)
      result.removeDuplicates();
  }

  HeuristicLevel OracleBase::maximizeController(OracleResult& result, const soplex::VectorRational& objective,
    const ObjectiveBound& objectiveBound, HeuristicLevel minHeuristic, HeuristicLevel maxHeuristic, bool& sort, bool& checkDups)
  {
    assert((heuristicLevel() == 0 && _nextOracle == NULL)
      || heuristicLevel() > 0 && _nextOracle != NULL);
    assert(minHeuristic <= maxHeuristic);

    // If requested, forward to next oracle.

    if (heuristicLevel() > maxHeuristic)
    {
      return _nextOracle->maximizeController(result, objective, objectiveBound, minHeuristic, maxHeuristic, sort, checkDups);
    }

    // Call implementation and check whether results are satisfactory.

    HeuristicLevel implHeuristicLevel = maximizeImplementation(result, objective, objectiveBound, minHeuristic, maxHeuristic,
      sort, checkDups);
    if (implHeuristicLevel == 0 || !result.rays.empty())
      return implHeuristicLevel;

    if (!result.points.empty())
    {
      result.computeMissingObjectiveValues();

      // Check if one of them satisfies the objective bound.

      bool satisfied = false;
      for (std::size_t i = 0; i < result.points.size(); ++i)
      {
        if (objectiveBound.satisfiedBy(result.points[i].objectiveValue))
        {
          satisfied = true;
          break;
        }
      }

      if (satisfied)
        return heuristicLevel();
    }

    // Answer not satisfactory.

    if (heuristicLevel() > minHeuristic || (result.points.empty() && result.rays.empty()))
    {
      return _nextOracle->maximizeController(result, objective, objectiveBound, minHeuristic, maxHeuristic, sort, checkDups);
    }
    else
    {
      return heuristicLevel();
    }
  }

  const LinearConstraint& OracleBase::currentFace()
  {
    return _currentFace;
  }

  FaceOracleBase::FaceOracleBase(const std::string& name, const std::shared_ptr<OracleBase>& nextOracle,
    std::size_t maxInfeasibleIterations, double initialM)
    : OracleBase(name, nextOracle), _maxInfeasibleIterations(maxInfeasibleIterations), _M(initialM), _completeFace()
  {

  }

  FaceOracleBase::~FaceOracleBase()
  {

  }

  void FaceOracleBase::initializeSpace(const Space& space)
  {
    OracleBase::initializeSpace(space);

    _modifiedObjective.reDim(_space.dimension(), false);
  }

  void FaceOracleBase::setFace(const LinearConstraint& newFace)
  {
    if (newFace == currentFace())
      return;

    OracleBase::setFace(newFace);
    assert(newFace == OracleBase::currentFace());

    if (newFace.definesTrivialFace())
      _factor = _M / currentFace().getMaximumNorm();
  }

  const LinearConstraint& FaceOracleBase::currentFace()
  {
    return _completeFace;
  }

  HeuristicLevel FaceOracleBase::maximizeController(OracleResult& result, const soplex::VectorRational& objective,
    const ObjectiveBound& objectiveBound, HeuristicLevel minHeuristic, HeuristicLevel maxHeuristic, bool& sort, bool& checkDups)
  {
    assert((heuristicLevel() == 0 && _nextOracle == NULL)
      || heuristicLevel() > 0 && _nextOracle != NULL);

    // If requested, forward to next oracle.

    if (heuristicLevel() > maxHeuristic)
    {
      return _nextOracle->maximizeController(result, objective, objectiveBound, minHeuristic, maxHeuristic, sort, checkDups);
    }

    // Get face (note that currentFace() always yields NULL.

    const LinearConstraint& face = OracleBase::currentFace();
    if (face.definesCompleteFace())
    {
      maximizeImplementation(result, objective, objectiveBound, minHeuristic, maxHeuristic, sort, checkDups);
    }
    else if (face.definesEmptyFace())
    {
      return 0;
    }
    else // We have a face and must tilt the objective vector.
    {
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

      // Check if c is the zero vector.

      if (maxObjectiveCoefficient == 0)
      {
        maximizeImplementation(result, objective, objectiveBound, minHeuristic, maxHeuristic, sort, checkDups);
      }
      else
      {
        // Copy objective since modification might only affect a subset of the coordinates.

        for (std::size_t v = 0; v < space().dimension(); ++v)
        {
          _modifiedObjective[v] = objective[v];
        }
        const Vector& faceNormal = face.normal();

        // Loop that increases M as long as results are not feasible.

        OracleResult unrestrictedResult;
        bool verifiedNonEmpty = false;
        ObjectiveBound modifiedObjectiveBound;
        Rational modifiedObjectiveShift;
        for (std::size_t iteration = 0; iteration <= _maxInfeasibleIterations; ++iteration)
        {
          if (iteration == _maxInfeasibleIterations)
          {
            // Setup objective via face's normal vector.

            for (std::size_t v = 0; v < space().dimension(); ++v)
              _modifiedObjective[v] = 0;
            for (int p = faceNormal.size() - 1; p >= 0; --p)
              _modifiedObjective[faceNormal.index(p)] = faceNormal.value(p);
            modifiedObjectiveBound.value = face.rhs();
            modifiedObjectiveBound.strict = false;
          }
          else
          {
            // Setup objective via tilting.

            Rational scaling = _factor / maxObjectiveCoefficient;
            modifiedObjectiveShift = scaling * face.rhs();
            modifiedObjectiveBound.value = objectiveBound.value + modifiedObjectiveShift;
            modifiedObjectiveBound.strict = objectiveBound.strict;
            for (int p = faceNormal.size() - 1; p >= 0; --p)
            {
              int v = faceNormal.index(p);
              _modifiedObjective[v] = objective[v] + scaling * faceNormal.value(p);
            }
          }

          // Carry out the optimization over P.

          maximizeImplementation(unrestrictedResult, _modifiedObjective, modifiedObjectiveBound, minHeuristic, maxHeuristic,
            sort, checkDups);

          // Check results.

          if (!unrestrictedResult.points.empty())
          {
            bool satisfied = false;
            std::size_t bestIndex = std::numeric_limits<std::size_t>::max();
            Rational bestValue = -infinity;
            for (std::size_t i = 0; i < unrestrictedResult.points.size(); ++i)
            {
              Vector& vector = unrestrictedResult.points[i].vector;
              if (face.evaluatePoint(vector) != 0)
                continue;

              Rational objectiveValue = objective * vector;
              if (objectiveBound.satisfiedBy(objectiveValue))
              {
                result.points.push_back(OracleResult::Point(vector, objectiveValue));
                satisfied = true;
              }
              else if (objectiveValue > bestValue)
              {
                bestIndex = i;
                bestValue = objectiveValue;
              }
            }

            // If we found a good feasible point, we can return.

            if (satisfied)
              return heuristicLevel();

            // If feasible points were found, but they were too bad, we still add the best for feasibility, aborting the loop.

            if (bestIndex < std::numeric_limits<std::size_t>::max())
            {
              result.points.push_back(OracleResult::Point(unrestrictedResult.points[bestIndex].vector, bestValue));
              break;
            }
          }
          else if (!unrestrictedResult.rays.empty())
          {
            for (std::size_t i = 0; i < unrestrictedResult.rays.size(); ++i)
            {
              Vector& vector = unrestrictedResult.rays[i].vector;
              if (face.evaluateRay(vector) != 0)
                continue;

              result.rays.push_back(OracleResult::Ray(vector));
            }

            if (!result.rays.empty())
              return heuristicLevel();
          }
          else
          {
            // Infeasible: Forward to next oracle.

            break;
          }

          _factor *= 2; // Increase current scaling.
          _M *= 2; // Increase scaling for future calls.
        }
      }
    }

    // Answer not satisfactory.

    if (heuristicLevel() > minHeuristic || (result.points.empty() && result.rays.empty()))
    {
      return _nextOracle->maximizeController(result, objective, objectiveBound, minHeuristic, maxHeuristic, sort, checkDups);
    }
    else
    {
      return heuristicLevel();
    }
  }

} /* namespace ipo */
