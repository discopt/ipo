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

  CacheOracle::CacheOracle(const std::string& name, const Space& space)
    : OracleBase(name, space), _uniquePoints(space.dimension()), _uniqueDirections(space.dimension())
  {
    OracleBase::initializedSpace();
  }

  CacheOracle::CacheOracle(const std::string& name, OracleBase* nextOracle)
    : OracleBase(name, nextOracle), _uniquePoints(nextOracle->space().dimension()), 
    _uniqueDirections(nextOracle->space().dimension())
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
    _faceDirections.clear();
    if (currentFace() != NULL)
    {
      for (UniqueVectors::Iterator iter = _uniquePoints.begin(); iter != _uniquePoints.end(); ++iter)
      {
        if (currentFace()->containsPoint(*iter))
          _facePoints.push_back(Data(*iter));
      }
      for (UniqueVectors::Iterator iter = _uniqueDirections.begin(); iter != _uniqueDirections.end(); ++iter)
      {
        if (currentFace()->containsDirection(*iter))
          _faceDirections.push_back(Data(*iter));
      }
    }
  }

  std::size_t CacheOracle::maximizeController(OracleResult& result, const soplex::VectorRational& objective,
    const ObjectiveBound& objectiveBound, std::size_t maxHeuristic, std::size_t minHeuristic, bool& sort, bool& checkDups)
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
      for (std::size_t i = 0; i < result.directions.size(); ++i)
      {
        if (addRay(result.directions[i].vector))
          ++numAddedRays;
      }
//       std::cerr << "CacheOracle: added " << numAddedPoints << " of " << result.points.size() << " points and "
//         << numAddedRays << " of " << result.directions.size() << " rays to cache." << std::endl;
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

    // Search directions.

    search(_faceDirections, approximateObjective, 0.0, false, searchResult);
    for (std::size_t i = 0; i < searchResult.size(); ++i)
    {
      Rational activity = objective * searchResult[i];
      if (activity <= 0)
        continue;
      result.directions.push_back(OracleResult::Direction(searchResult[i]));
    }
    if (!result.directions.empty())
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

    if (currentFace() != NULL)
    {
      if (currentFace()->containsPoint(p))
        _facePoints.push_back(Data(p));
    }
    return true;
  }

  bool CacheOracle::addRay(const Vector& ray)
  {
    Vector r = ray;
    if (!_uniqueDirections.insert(r))
      return false;

    if (currentFace() != NULL)
    {
      if (currentFace()->containsDirection(r))
        _faceDirections.push_back(Data(r));
    }
    return true;
  }

  
/*
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
  }*/

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
