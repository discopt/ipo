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
    _faceDirections.clear();
    _endFaceDirections = 0;
  }
  
  std::size_t CacheOracle::maximizeImplementation(OracleResult& result, const VectorRational& objective,
    const ObjectiveBound& objectiveBound, std::size_t minHeuristic, std::size_t maxHeuristic, bool& sort, bool& checkDups)
  {
    // Update face indices since in the meantime points / directions could have been added.

    updateFaceIndices(_points, _facePoints, _endFacePoints, true);
    updateFaceIndices(_directions, _faceDirections, _endFaceDirections, false);

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
      result.directions.push_back(OracleResult::Direction(_directions.vector(d)));
      result.directions.back().index = d;
    }
    if (!result.directions.empty())
      return heuristicLevel();

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
        result.points.push_back(OracleResult::Point(_points.vector(p)));
        result.points.back().index = p;
        result.points.back().objectiveValue = activity;
      }
    }

    return heuristicLevel();
  }

  CacheOracle::VectorStats::VectorStats()
  {

  }

  CacheOracle::VectorStats::VectorStats(double theObjectiveValue, std::size_t theSparsity,
    std::size_t theIndex) : sparsity(theSparsity), index(theIndex)
  {
    valueMantissa = frexp(theObjectiveValue, &valueExponent);
    value = theObjectiveValue;
  }

  CacheOracle::VectorStats& CacheOracle::VectorStats::operator=(
    const CacheOracle::VectorStats& other)
  {
    valueExponent = other.valueExponent;
    valueMantissa = other.valueMantissa;
    sparsity = other.sparsity;
    index = other.index;
    return *this;
  }

  bool CacheOracle::VectorStats::operator<(const CacheOracle::VectorStats& other) const
  {
    if (valueMantissa > 0 && other.valueMantissa <= 0)
      return true;
    if (valueMantissa >= 0 && other.valueMantissa > 0)
      return false;
    if (valueExponent > other.valueExponent)
      return true;
    if (valueExponent < other.valueExponent)
      return false;
    if (sparsity < other.sparsity)
      return true;
    if (sparsity > other.sparsity)
      return false;
    return valueMantissa > other.valueMantissa;
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

    for (std::size_t i = 0; i < _vectorStats.size(); ++i)
    {
      std::size_t bestIndex = i;
      for (std::size_t j = i + 1; j < _vectorStats.size(); ++j)
      {
        if (_vectorStats[j] < _vectorStats[bestIndex])
          bestIndex = j;
      }
      if (bestIndex != i)
        std::swap(_vectorStats[i], _vectorStats[bestIndex]);
    }
//    std::sort(_vectorStats.begin(), _vectorStats.end()); TODO: Why does std::sort fail?

    // Extract the top element from every exponent group.

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
        result.push_back(i);
      }
    }
  }

}
