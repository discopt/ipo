#include <ipo/polyhedron.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>

// TODO: Create shared pointer to the CacheOptimizationOracle that is separate from _optimization
// array.

namespace ipo
{
  class CacheOptimizationOracle : public OptimizationOracle
  {
  public:
    CacheOptimizationOracle(std::shared_ptr<Space> space)
      : OptimizationOracle("CacheOptimizationOracle"), _normalizedRayEpsilon(1.0e-6),
      _normalizedPointEpsilon(1.0e-6), _queryCount(0)
    {
      _space = space;
    }

    void setEpsilon(double normalizedRayEpsilon, double normalizedPointEpsilon)
    {
      _normalizedRayEpsilon = normalizedRayEpsilon;
      _normalizedPointEpsilon = normalizedPointEpsilon;
    }

    bool isExact() const override
    {
      return true;
    }

    virtual void maximize(const double* objectiveVector, const OptimizationOracle::Query& query,
      OptimizationOracle::Result& result) override
    {
      double objectiveNorm = 0.0;
      for (std::size_t v = 0; v < space()->dimension(); ++v)
      {
        double c = objectiveVector[v];
        objectiveNorm += c * c;
      }
      objectiveNorm = (objectiveNorm == 0.0) ? 1.0 : sqrt(objectiveNorm);
      ++_queryCount;

      // TODO: Respect time limit.
      
      // Compute scalar product with each ray.

      for (auto& rayData : _rays)
      {
#if defined(IPO_WITH_GMP)
        rayData.product = (rayData.vector * objectiveVector) * rayData.normalization;
        rayData.rationalProduct = 0;
#endif /* IPO_WITH_GMP */
      }

      // Sort rays according to product.

      std::sort(_rays.begin(), _rays.end());

      // Add best rays as long as they have sufficiently positive product.

      assert(result.rays.size() + result.points.size() <= query.maxNumSolutions);
      std::size_t maxNum = std::min(_rays.size(),
        query.maxNumSolutions - result.rays.size() - result.points.size());
      double epsilon = _normalizedRayEpsilon * objectiveNorm;
      for (std::size_t i = 0; i < maxNum; ++i)
      {
        if (_rays[i].product > epsilon)
        {
          result.rays.push_back(OptimizationOracle::Result::Ray(_rays[i].vector));
          _rays[i].lastSuccess = _queryCount;
        }
        else
          break;
      }

      // Skip further computations if we are done already.

      if (query.maxNumSolutions == result.rays.size() + result.points.size())
        return;

      // Compute scalar product with each point.

      for (auto& pointData : _points)
      {
#if defined(IPO_WITH_GMP)
        pointData.product = pointData.vector * objectiveVector;
        pointData.rationalProduct = 0;
#endif /* IPO_WITH_GMP */
      }

      // Sort points according to product.

      std::sort(_points.begin(), _points.end());

      // Add best points as longs as they have sufficiently positive product.

      assert(result.rays.size() + result.points.size() <= query.maxNumSolutions);
      maxNum = std::min(_points.size(),
        query.maxNumSolutions - result.rays.size() - result.points.size());
      double threshold = query.minObjectiveValue + _normalizedPointEpsilon * objectiveNorm;
      for (std::size_t i = 0; i < maxNum; ++i)
      {
        if (_points[i].product < threshold)
          break;
        result.points.push_back(OptimizationOracle::Result::Point(_points[i].vector,
          _points[i].product));
        _points[i].lastSuccess = _queryCount;
      }
    }

#if defined(IPO_WITH_GMP)

    virtual void maximize(const mpq_class* objectiveVector, const OptimizationOracle::Query& query,
      OptimizationOracle::Result& result) override
    {
      ++_queryCount;

      // TODO: Respect time limit.
      
      // Compute scalar product with each ray.

      for (auto& rayData : _rays)
      {
#if defined(IPO_WITH_GMP)
        rayData.rationalProduct = rayData.vector * objectiveVector;
        if (rayData.rationalProduct > 0)
          rayData.product = rayData.rationalProduct.get_d() * rayData.normalization;
#endif /* IPO_WITH_GMP */
      }

      // Sort rays according to product.

      std::sort(_rays.begin(), _rays.end());

      // Add best rays as long as they have sufficiently positive product.

      assert(result.rays.size() + result.points.size() <= query.maxNumSolutions);
      std::size_t maxNum = std::min(_rays.size(),
        query.maxNumSolutions - result.rays.size() - result.points.size());
      for (std::size_t i = 0; i < maxNum; ++i)
      {
        if (_rays[i].product > 0)
        {
          result.rays.push_back(OptimizationOracle::Result::Ray(_rays[i].vector));
          _rays[i].lastSuccess = _queryCount;
        }
        else
          break;
      }

      // Skip further computations if we are done already.

      if (query.maxNumSolutions == result.rays.size() + result.points.size())
        return;

      // Compute scalar product with each point.

      for (auto& pointData : _points)
      {
#if defined(IPO_WITH_GMP)
        pointData.rationalProduct = pointData.vector * objectiveVector;
        pointData.product = pointData.rationalProduct.get_d();
#endif /* IPO_WITH_GMP */
      }

      // Sort points according to product.

      std::sort(_points.begin(), _points.end());

      // Add best points as longs as they have sufficiently positive product.
      
      assert(result.rays.size() + result.points.size() <= query.maxNumSolutions);
      maxNum = std::min(_points.size(),
        query.maxNumSolutions - result.rays.size() - result.points.size());
      double threshold = query.minObjectiveValue;
      for (std::size_t i = 0; i < maxNum; ++i)
      {
        if (_points[i].product <= threshold)
          break;
        result.points.push_back(OptimizationOracle::Result::Point(_points[i].vector,
          _points[i].product));
        _points[i].lastSuccess = _queryCount;
      }
    }

#endif /* IPO_WITH_GMP */

    void reduceCacheSize()
    {

    }

  protected:
    struct Data
    {
      double product;
#if defined(IPO_WITH_GMP)
      mpq_class rationalProduct;
#endif /* IPO_WITH_GMP */
      std::size_t lastSuccess; /// Largest query number at which this vector was returned.
      Vector vector; /// The cached point/ray.
      double normalization; /// Inverse of Euclidean norm of vector.

      Data(Vector v, std::size_t ls)
        : lastSuccess(ls), vector(v)
      {
        normalization = vector.squaredRealNorm();
        normalization = (normalization == 0.0) ? 1.0 : 1.0 / sqrt(normalization);
      }

      inline bool operator<(const Data& other) const
      {
        if (product > other.product)
          return true;
        if (product < other.product)
          return false;

#if defined(IPO_WITH_GMP)
        int result = cmp(rationalProduct, other.rationalProduct);
        if (result != 0)
          return result > 0;          
#endif /* IPO_WITH_GMP */

        return lastSuccess > other.lastSuccess;
      }
    };

    std::vector<Data> _points;
    std::vector<Data> _rays;
    double _normalizedRayEpsilon;
    double _normalizedPointEpsilon;
    std::size_t _queryCount;
  };
  
  template <>
  Polyhedron::Data<OptimizationOracle>::Data(std::shared_ptr<OptimizationOracle> o, bool cache)
    : sumRunningTime(0.0), sumSuccess(0), oracle(o), isCache(cache)
  {

  }

  template <>
  Polyhedron::Data<SeparationOracle>::Data(std::shared_ptr<SeparationOracle> o, bool cache)
    : sumRunningTime(0.0), sumSuccess(0), oracle(o), isCache(cache)
  {

  }

  template <>
  bool Polyhedron::Data<OptimizationOracle>::operator<(const Data<OptimizationOracle>& other) const
  {
    return priority < other.priority;
  }

  template <>
  void Polyhedron::Data<OptimizationOracle>::swap(Data<OptimizationOracle>& other)
  {
    std::swap(sumRunningTime, other.sumRunningTime);
    std::swap(sumSuccess, other.sumSuccess);
    std::swap(priority, other.priority);
    oracle.swap(other.oracle);
    history.swap(other.history);
    std::swap(isCache, other.isCache);
  }

  template <>
  void Polyhedron::Data<OptimizationOracle>::updateHistory(double runningTime, bool success,
    std::size_t historySize)
  {
    history.push_back(QueryStatistics(runningTime, success));
    sumRunningTime += runningTime;
    if (success)
      ++sumSuccess;

    // If history is too long, remove the first entry.

    if (history.size() > historySize)
    {
      sumRunningTime -= history.front().runningTime;
      if (history.front().success)
        --sumSuccess;
      history.pop_front();
    }

    // We always add a bonus of +1 to the number of successful runs. This avoids division by 0
    // and gives some advantage to unsuccessful oracles that do not require too much time.

    priority = sumRunningTime / (sumSuccess + 1);  
  }

  Polyhedron::QueryStatistics::QueryStatistics(double rt, bool succ)
    : runningTime(rt), success(succ)
  {

  }

  Polyhedron::Polyhedron(std::shared_ptr<OptimizationOracle> optimizationOracle)
    : _historySize(16)
  {
    auto cache = std::make_shared<CacheOptimizationOracle>(optimizationOracle->space());
    _optimization.push_back(Data<OptimizationOracle>(cache, true));
    _optimization.push_back(Data<OptimizationOracle>(optimizationOracle));
  }

  Polyhedron::Polyhedron(std::shared_ptr<SeparationOracle> o)
    : _historySize(16)
  {
    _separation.push_back(Data<SeparationOracle>(o));
  }

  Polyhedron::~Polyhedron()
  {

  }

  std::shared_ptr<Space> Polyhedron::space() const
  {
    if (!_optimization.empty())
      return _optimization.front().oracle->space();
    if (!_separation.empty())
      return _separation.front().oracle->space();
    return std::shared_ptr<Space>();
  }

  void Polyhedron::maximize(const double* objectiveVector, const OptimizationOracle::Query& query,
    OptimizationOracle::Result& result)
  {
    tuneOracles();

    result.reset();
    std::size_t lastOracle = 0;
    for (; lastOracle < _optimization.size(); ++lastOracle)
    {
      Data<OptimizationOracle>& data = _optimization[lastOracle];
      // Run current oracle.

      std::chrono::time_point<std::chrono::system_clock> started = std::chrono::system_clock::now();
      data.oracle->maximize(objectiveVector, query, result);
      std::chrono::duration<double> duration = std::chrono::system_clock::now() - started;

      data.updateHistory(duration.count(), !result.isInfeasible(), _historySize);

      // If the current oracle found a ray or a point, we stop.
      if (!result.isInfeasible())
        return;
    }
  }

  void Polyhedron::maximize(const mpq_class* objectiveVector,
    const OptimizationOracle::Query& query, OptimizationOracle::Result& result)
  {
    tuneOracles();

    result.reset();
    std::size_t lastOracle = 0;
    for (; lastOracle < _optimization.size(); ++lastOracle)
    {
      Data<OptimizationOracle>& data = _optimization[lastOracle];
      // Run current oracle.

      std::chrono::time_point<std::chrono::system_clock> started = std::chrono::system_clock::now();
      data.oracle->maximize(objectiveVector, query, result);
      std::chrono::duration<double> duration = std::chrono::system_clock::now() - started;

      data.updateHistory(duration.count(), !result.isInfeasible(), _historySize);

      // If the current oracle found a ray or a point, we stop.
      if (!result.isInfeasible())
        return;
    }
  }

  void Polyhedron::tuneOracles()
  {
    // In case the expected successful running times have changed, we reorder the oracles such that
    // the one with the smallest value is queried first.

    if (!std::is_sorted(_optimization.begin(), _optimization.end()))
      std::sort(_optimization.begin(), _optimization.end());

    if (_optimization.size() > 1)
    {
      // Maybe adjust allowed size of cache oracle.
      std::size_t cacheOracle = std::numeric_limits<std::size_t>::max();
      double maxRunningTime = 0.0;
      for (std::size_t o = 0; o < _optimization.size(); ++o)
      {
        if (_optimization[o].isCache)
        {
          cacheOracle = o;
          continue;
        }

        if (!_optimization[o].history.empty())
        {
          double runningTime = _optimization[o].sumRunningTime / _optimization[o].history.size();
          maxRunningTime = std::max(maxRunningTime, runningTime);
        }
      }

      if (cacheOracle < std::numeric_limits<std::size_t>::max())
      {
        // If the expected running time of cache oracle is at least that of the slowest oracle, then
        // we have to reduce the cache size.

        double cacheRunningTime = _optimization[cacheOracle].sumRunningTime
          / _optimization[cacheOracle].history.size();
        if (cacheRunningTime >= maxRunningTime)
        {
          auto oracle = dynamic_cast<CacheOptimizationOracle&>(*_optimization[cacheOracle].oracle);
          oracle.reduceCacheSize();
        }
      }
    }
  }

} /* namespace ipo */
