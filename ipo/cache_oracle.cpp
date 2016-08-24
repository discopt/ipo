#include "cache_oracle.h"

using namespace soplex;

namespace ipo {
  
  CacheOracle::Data::Data(Vector& vec)
    : vector(vec)
  {

  }

  CacheOracle::Data::~Data()
  {

  }

  void CacheOracle::Data::updateObjective(const soplex::VectorReal& approximateObjective, double approximateObjectiveBound)
  {
    double approximateObjectiveValue = 0.0;
    for (std::size_t p = 0; p < vector.size(); ++p)
      approximateObjectiveValue += approximateObjective[vector.index(p)] * vector.approximation(p);

    valueMantissa = frexp(approximateObjectiveValue, &valueExponent);
  }

  bool CacheOracle::Data::operator<(const CacheOracle::Data& other) const
  {
    if (valueMantissa > 0 && other.valueMantissa <= 0)
      return true;
    if (valueMantissa >= 0 && other.valueMantissa > 0)
      return false;
    if (valueExponent > other.valueExponent)
      return true;
    if (valueExponent < other.valueExponent)
      return false;
    if (vector.size()  < other.vector.size())
      return true;
    if (vector.size() > other.vector.size())
      return false;
    return valueMantissa > other.valueMantissa;
  }

  CacheOracle::CacheOracle(const std::string& name, const std::shared_ptr<OracleBase>& nextOracle)
    : OracleBase(name, nextOracle), _uniquePoints(nextOracle->space().dimension()), _uniqueRays(nextOracle->space().dimension())
  {
    assert(nextOracle);

    OracleBase::initializeSpace(nextOracle->space());
  }

  CacheOracle::~CacheOracle()
  {

  }

  void CacheOracle::setFace(const LinearConstraint& newFace)
  {
    if (newFace == currentFace())
      return;

    OracleBase::setFace(newFace);

    _facePoints.clear();
    _faceRays.clear();
    if (!currentFace().definesCompleteFace())
    {
      for (UniqueVectors::Iterator iter = _uniquePoints.begin(); iter != _uniquePoints.end(); ++iter)
      {
        if (currentFace().evaluatePoint(*iter) == 0)
          _facePoints.push_back(Data(*iter));
      }
      for (UniqueVectors::Iterator iter = _uniqueRays.begin(); iter != _uniqueRays.end(); ++iter)
      {
        if (currentFace().evaluateRay(*iter) == 0)
          _faceRays.push_back(Data(*iter));
      }
    }
  }

  std::size_t CacheOracle::maximizeController(OracleResult& result, const soplex::VectorRational& objective,
    const ObjectiveBound& objectiveBound, std::size_t minHeuristic, std::size_t maxHeuristic, bool& sort, bool& checkDups)
  {
    std::size_t level = OracleBase::maximizeController(result, objective, objectiveBound, maxHeuristic, minHeuristic, sort, 
      checkDups);

    if (level < heuristicLevel())
    {
      // Result was not produced by the cache.

      std::size_t numAddedPoints = 0;
      for (std::size_t i = 0; i < result.points.size(); ++i)
      {
        if (addPoint(result.points[i].vector))
          ++numAddedPoints;
      }
      std::size_t numAddedRays = 0;
      for (std::size_t i = 0; i < result.rays.size(); ++i)
      {
        if (addRay(result.rays[i].vector))
          ++numAddedRays;
      }
    }

    return level;
  }

  std::size_t CacheOracle::maximizeImplementation(OracleResult& result, const soplex::VectorRational& objective,
    const ObjectiveBound& objectiveBound, std::size_t minHeuristic, std::size_t maxHeuristic, bool& sort, bool& checkDups)
  {
    soplex::DVectorReal approximateObjective(objective.dim());
    for (std::size_t i = 0; i < approximateObjective.dim(); ++i)
      approximateObjective[i] = double(objective[i]);
    std::vector<Vector> searchResult;

    // Search rays.

    search(_faceRays, approximateObjective, 0.0, false, searchResult);
    for (std::size_t i = 0; i < searchResult.size(); ++i)
    {
      Rational activity = objective * searchResult[i];
      if (activity <= 0)
        continue;
      result.rays.push_back(OracleResult::Ray(searchResult[i]));
    }
    if (!result.rays.empty())
      return heuristicLevel();

    // Search points.

    search(_facePoints, approximateObjective, objectiveBound.value > minusInfinity ? double(objectiveBound.value) : 0.0, true, 
      searchResult);

    bool foundSatisfying = false;
    for (std::size_t i = 0; i < searchResult.size(); ++i)
    {
      Rational activity = objective * searchResult[i];

      bool satisfied = objectiveBound.satisfiedBy(activity);
      if (i == 0 || satisfied)
      {
        result.points.push_back(OracleResult::Point(searchResult[i]));
        result.points.back().objectiveValue = activity;
      }
    }

    return heuristicLevel();
  }

  bool CacheOracle::addPoint(const Vector& point)
  {
    Vector p = point;
    if (!_uniquePoints.insert(p))
      return false;

    if (currentFace().definesCompleteFace() || currentFace().evaluatePoint(p) == 0)
      _facePoints.push_back(Data(p));
    return true;
  }

  bool CacheOracle::addRay(const Vector& ray)
  {
    Vector r = ray;
    if (!_uniqueRays.insert(r))
      return false;

    if (currentFace().definesCompleteFace() || currentFace().evaluateRay(r) == 0)
      _faceRays.push_back(Data(r));
    return true;
  }

  void CacheOracle::search(std::vector<Data>& vectors, const soplex::VectorReal& approximateObjective,
    double approximateObjectiveBound, bool handlingPoints, std::vector<Vector>& result)
  {
    // Fill stats vector.

    for (std::size_t i = 0; i < vectors.size(); ++i)
    {
      vectors[i].updateObjective(approximateObjective, approximateObjectiveBound);
    }

    // Sort it.

    for (std::size_t i = 0; i < vectors.size(); ++i)
    {
      std::size_t bestIndex = i;
      for (std::size_t j = i + 1; j < vectors.size(); ++j)
      {
        if (vectors[j] < vectors[bestIndex])
          bestIndex = j;
      }
      if (bestIndex != i)
        std::swap(vectors[i], vectors[bestIndex]);
    }
//    std::sort(_vectorStats.begin(), _vectorStats.end()); TODO: Why does std::sort fail?

    // Extract the top element from every exponent group.

    int lastSign = 2;
    int lastExponent = std::numeric_limits<int>::min();
    for (std::size_t i = 0; i < vectors.size(); ++i)
    {
      const Data& data = vectors[i];
      int sign = (data.valueMantissa > 0 ? 1 : 0) - (data.valueMantissa < 0 ? 1 : 0);
      if (sign != lastSign || data.valueExponent != lastExponent)
      {
        lastSign = sign;
        lastExponent = data.valueExponent;
        result.push_back(vectors[i].vector);
      }
    }
  }

}
