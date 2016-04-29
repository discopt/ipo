#include "cache_oracle.h"

using namespace soplex;

namespace ipo {

  CacheOracle::CacheOracle(const std::string& name, const Space& space,
    const UniqueRationalVectorsBase& points, const UniqueRationalVectorsBase& directions)
    : OracleBase(name, space), _points(points), _directions(directions),
    _endFacePoints(0), _endFaceDirections(0)
  {
    OracleBase::initializedSpace();
  }

  CacheOracle::CacheOracle(const std::string& name, OracleBase* nextOracle,
    const UniqueRationalVectorsBase& points, const UniqueRationalVectorsBase& directions)
    : OracleBase(name, nextOracle), _points(points), _directions(directions),
    _endFacePoints(0), _endFaceDirections(0)
  {
    OracleBase::initializedSpace();
  }
  
  CacheOracle::~CacheOracle()
  {

  }

  void CacheOracle::setFace(Face* newFace)
  {
    if (newFace == currentFace())
      return;

    OracleBase::setFace(newFace);

    _facePoints.clear();
    _endFacePoints = 0;
    updateFaceIndices(_points, _facePoints, _endFacePoints, true);

    _faceDirections.clear();
    _endFaceDirections = 0;
    updateFaceIndices(_directions, _faceDirections, _endFaceDirections, false);
  }

  void CacheOracle::maximize(OracleResult& result, const VectorRational& objective,
    const ObjectiveBound& objectiveBound, std::size_t maxHeuristic, std::size_t minHeuristic)
  {
    assert((thisHeuristic() == 0 && _nextOracle == NULL)
      || thisHeuristic() > 0 && _nextOracle != NULL);

    // Forward call if requested.

    if (thisHeuristic() > maxHeuristic)
      return _nextOracle->maximize(result, objective, objectiveBound, maxHeuristic, minHeuristic);

    result.buildStart(objective);
    DVectorReal approxObjective(objective.dim());
    approxObjective = objective;
    std::vector<std::size_t> searchResult;

    // Search directions.

    search(_directions, _faceDirections, approxObjective, 0.0, false, searchResult);
    for (std::size_t i = 0; i < searchResult.size(); ++i)
    {
      std::size_t d = _faceDirections[searchResult[i]];
      Rational activity = *_directions.vector(d) * objective;
      if (activity <= 0)
        continue;
      result.buildAddDirection(_directions.vector(d));
      result.directions.back().index = d;
    }
    if (!result.directions.empty())
      return result.buildFinish(thisHeuristic(), false, false, false);

    // Search points.

    search(_points, _facePoints, approxObjective,
      objectiveBound.value > -infinity ? double(objectiveBound.value) : 0.0, true, searchResult);

    bool foundSatisfying = false;
    for (std::size_t i = 0; i < searchResult.size(); ++i)
    {
      std::size_t p = _facePoints[searchResult[i]];
      Rational activity = *_points.vector(p) * objective;
      bool satisfied = objectiveBound.satisfiedBy(activity);
      if (i == 0 || satisfied)
      {
        result.buildAddPoint(_points.vector(p));
        result.points.back().index = p;
        result.points.back().objectiveValue = activity;
        if (satisfied)
          foundSatisfying = true;
      }
    }

    if (foundSatisfying || thisHeuristic() <= minHeuristic)
      return result.buildFinish(thisHeuristic(), false, true, true);

    // Forward call if no sufficiently good points were found.

    return _nextOracle->maximize(result, objective, objectiveBound, maxHeuristic, minHeuristic);
  }

  void CacheOracle::updateFaceIndices(const UniqueRationalVectorsBase& vectors,
    CacheOracle::FaceIndices& faceIndices, std::size_t& end, bool handlingPoints)
  {
    for (; end < vectors.size(); ++end)
    {
      if (currentFace() != NULL)
      {
        if ((handlingPoints &&
          (*vectors[end] * currentFace()->denseNormal() != currentFace()->rhs()))
          || (!handlingPoints &&
          (*vectors[end] * currentFace()->denseNormal() != 0)))
        {
          continue;
        }
      }
      faceIndices.push_back(end);
    }
  }

  void CacheOracle::search(const UniqueRationalVectorsBase& vectors,
    const FaceIndices& faceIndices, const VectorReal& approxObjective, double approxObjectiveBound,
    bool handlingPoints, std::vector<std::size_t>& result)
  {
    // Fill stats vector.

    _vectorStats.resize(faceIndices.size());
    for (std::size_t i = 0; i < faceIndices.size(); ++i)
    {
      const SVectorReal& approxVector = *vectors.approximation(faceIndices[i]);
      double approxActivity = approxVector * approxObjective;
      _vectorStats[i] = VectorStats(approxActivity - approxObjectiveBound,
        approxVector.size(), faceIndices[i]);
    }

    // Sort it.

    std::sort(_vectorStats.begin(), _vectorStats.end());

    // Extract the top element from every exponent group.

    std::cerr << "CacheOracle::search()" << std::endl;

    int lastSign = 2;
    int lastExponent = std::numeric_limits<int>::min();
    for (std::size_t i = 0; i < _vectorStats.size(); ++i)
    {
      const VectorStats& stats = _vectorStats[i];
      int sign = (stats.valueMantissa > 0 ? 1 : 0) - (stats.valueMantissa < 0 ? 1 : 0);
      if (sign != lastSign || stats.valueExponent != lastExponent)
      {
        lastSign = sign;
        lastExponent = stats.valueExponent;
        std::cerr << "! ";
        result.push_back(i);
      }
      std::cerr << stats.value << ", sparsity = " << stats.sparsity << std::endl;
    }
  }

}
