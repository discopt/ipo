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

  CacheOracle::CacheOracle(const std::shared_ptr<OracleBase>& nextOracle, Behavior outerBehavior, Behavior innerBehavior)
    : OracleBase("CacheOracle", nextOracle), _uniquePoints(nextOracle->space().dimension()),
    _uniqueRays(nextOracle->space().dimension()), _outerBehavior(outerBehavior), _innerBehavior(innerBehavior)
  {
    assert(nextOracle);

    OracleBase::initializeSpace(nextOracle->space());

    _inequalities.setIntParam(SoPlex::OBJSENSE, SoPlex::OBJSENSE_MAXIMIZE);
    _inequalities.setIntParam(SoPlex::SOLVEMODE, SoPlex::SOLVEMODE_RATIONAL);
    _inequalities.setIntParam(SoPlex::SYNCMODE, SoPlex::SYNCMODE_AUTO);
    _inequalities.setRealParam(SoPlex::FEASTOL, 0.0);
    _inequalities.setBoolParam(SoPlex::RATREC, true);
    _inequalities.setBoolParam(SoPlex::RATFAC, true);
    _inequalities.setIntParam(SoPlex::VERBOSITY, SoPlex::VERBOSITY_ERROR);

    LPColSetRational cols;
    DSVectorRational zeroVector;
    for (std::size_t v = 0; v < space().dimension(); ++v)
    {
      cols.add(Rational(0), -infinity, zeroVector, infinity);
    }
    _inequalities.addColsRational(cols);

  }

  CacheOracle::~CacheOracle()
  {

  }

  void CacheOracle::setOuterBehavior(CacheOracle::Behavior outerBehavior)
  {
    _outerBehavior = outerBehavior;
    if (outerBehavior == DISABLED)
    {
      _uniquePoints.clear();
      _uniqueRays.clear();
      _facePoints.clear();
      _faceRays.clear();
    }
  }

  void CacheOracle::setInnerBehavior(CacheOracle::Behavior innerBehavior)
  {
    _innerBehavior = innerBehavior;

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

  HeuristicLevel CacheOracle::maximizeController(OracleResult& result, const soplex::VectorRational& objective,
    const ObjectiveBound& objectiveBound, HeuristicLevel minHeuristic, HeuristicLevel maxHeuristic, bool& sort, bool& checkDups)
  {
    HeuristicLevel level = OracleBase::maximizeController(result, objective, objectiveBound, minHeuristic, maxHeuristic, sort,
      checkDups);

    if (level > 0 && _outerBehavior == CACHE_AND_SEARCH && result.isFeasible())
    {
      if (sort)
      {
        result.computeMissingObjectiveValues();
        result.sortPoints();
        sort = false;
      }

      _inequalities.changeObjRational(objective);
      SPxSolver::Status status = _inequalities.solve();

      if (status == SPxSolver::OPTIMAL && _inequalities.objValueRational() == result.points.front().objectiveValue)
      {
        level = 0;
      }
    }
    else if (level == 0 && _outerBehavior != DISABLED && result.isFeasible())
    {
      if (sort)
      {
        result.computeMissingObjectiveValues();
        result.sortPoints();
        sort = false;
      }

      const Rational& optimum = result.points.front().objectiveValue;
      std::size_t numNonzeros = 0;
      std::size_t nonzero = 0;
      for (std::size_t v = 0; v < space().dimension(); ++v)
      {
        if (objective[v] != 0)
        {
          numNonzeros++;
          nonzero = v;
          if (numNonzeros >= 2)
            break;
        }
      }
      if (numNonzeros == 1)
      {
        // Update bound constraint.

        if (objective[nonzero] > 0)
        {
          Rational newBound = optimum / objective[nonzero];
          if (newBound < _inequalities.upperRational(nonzero))
            _inequalities.changeBoundsRational(nonzero, _inequalities.lowerRational(nonzero), newBound);
        }
        else
        {
          assert(objective[nonzero] < 0);
          Rational newBound = optimum / objective[nonzero];
          if (newBound > _inequalities.lowerRational(nonzero))
            _inequalities.changeBoundsRational(nonzero, newBound, _inequalities.upperRational(nonzero));
        }
      }
      else if (numNonzeros > 1)
      {
        DSVectorRational normal(objective);
        _inequalities.addRowRational(LPRowRational(-infinity, normal, optimum));
      }
    }

    if (_innerBehavior != DISABLED && level < heuristicLevel())
    {
      // Result was not produced by the cache. so we have to add it.

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

  HeuristicLevel CacheOracle::maximizeImplementation(OracleResult& result, const soplex::VectorRational& objective,
    const ObjectiveBound& objectiveBound, HeuristicLevel minHeuristic, HeuristicLevel maxHeuristic, bool& sort, bool& checkDups)
  {
    if (_innerBehavior != CACHE_AND_SEARCH)
      return heuristicLevel();

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
      return 0;
    else
      searchResult.clear();

    // Search points.

    search(_facePoints, approximateObjective, objectiveBound.value > -infinity ? double(objectiveBound.value) : 0.0, true,
      searchResult);

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
    sort = true;

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
    Vector r = integralScaled(ray);
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
