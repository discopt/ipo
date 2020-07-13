#pragma once

#include <iostream>

#include <ipo/config.hpp>
#include <ipo/export.hpp>
#include <ipo/space.hpp>
#include <ipo/oracles.hpp>

#include <chrono>
#include <memory>
#include <deque>

namespace ipo
{
  template <typename T>
  class CacheOptimizationOracle: public OptimizationOracle<T>
  {
  public:
    CacheOptimizationOracle(std::shared_ptr<Space> space)
      : OptimizationOracle<T>("CacheOptimizationOracle"), _normalizedRayEpsilon(1.0e-6),
      _normalizedPointEpsilon(1.0e-6), _queryCount(0)
    {
      this->_space = space;
    }

    void setEpsilon(double normalizedRayEpsilon, double normalizedPointEpsilon)
    {
      _normalizedRayEpsilon = normalizedRayEpsilon;
      _normalizedPointEpsilon = normalizedPointEpsilon;
    }

    virtual typename OptimizationOracle<T>::Result maximizeDouble(const double* objectiveVector,
      const typename OptimizationOracle<T>::Query& query) override
    {
      typename OptimizationOracle<T>::Result result;
      ++_queryCount;

      // Compute norm of objective vector.
      double objectiveNormalization = 0.0;
      for (std::size_t v = 0; v < this->space()->dimension(); ++v)
      {
        double x = double(objectiveVector[v]);
        objectiveNormalization += x*x;
      }
      objectiveNormalization = sqrt(objectiveNormalization);

      // TODO: Respect time limit.

      // Compute scalar product with each normalized ray.
      for (auto& rayData : _rays)
      {
        rayData.product = 0.0;
        for (const auto& iter : rayData.vector)
          rayData.product += double(objectiveVector[iter.first]) * double(iter.second);
        rayData.product *= rayData.normalization;
      }

      // Sort rays according to product.

      std::sort(_rays.begin(), _rays.end());

      // Add best rays as long as they have sufficiently positive product.
      double epsilon = objectiveNormalization * _normalizedRayEpsilon;
      for (auto& rayData : _rays)
      {
        T product = 0.0;
        for (const auto& iter : rayData.vector)
          product += double(objectiveVector[iter.first]) * double(iter.second);
        product *= rayData.normalization;

        if (product > epsilon)
        {
          result.rays.push_back(typename OptimizationOracle<T>::Result::Ray(rayData.vector));
          rayData.lastSuccess = _queryCount;
          if (result.rays.size() >= query.maxNumSolutions)
          {
            assert(result.checkPointsSorted());
            return result;
          }
        }
        else
          break;
      }

      // Compute scalar product with each point.

      for (auto& pointData : _points)
      {
        pointData.product = 0.0;
        for (const auto& iter : pointData.vector)
          pointData.product += double(objectiveVector[iter.first]) * double(iter.second);
      }

      // Sort points according to product.

      std::sort(_points.begin(), _points.end());

      // Add best points as longs as they have sufficiently positive product.

      double threshold = query.minObjectiveValue;
      epsilon = objectiveNormalization * _normalizedPointEpsilon;
      for (auto& pointData : _points)
      {
        T product = 0.0;
        for (const auto& iter : pointData.vector)
          product += double(objectiveVector[iter.first]) * double(iter.second);

        if (product <= threshold)
          break;
        result.points.push_back(typename OptimizationOracle<T>::Result::Point(pointData.vector,
          pointData.product));
        pointData.lastSuccess = _queryCount;
      }

      assert(result.checkPointsSorted());
      return result;
    }

    virtual typename OptimizationOracle<T>::Result maximize(const T* objectiveVector,
      const typename OptimizationOracle<T>::Query& query) override
    {
      typename OptimizationOracle<T>::Result result;
      ++_queryCount;

      // Compute norm of objective vector.
      double objectiveNormalization = 0.0;
      for (std::size_t v = 0; v < this->space()->dimension(); ++v)
      {
        double x = double(objectiveVector[v]);
        objectiveNormalization += x*x;
      }
      objectiveNormalization = sqrt(objectiveNormalization);

      // TODO: Respect time limit.

      // Compute scalar product with each normalized ray.
      for (auto& rayData : _rays)
      {
        rayData.product = 0.0;
        for (const auto& iter : rayData.vector)
          rayData.product += double(objectiveVector[iter.first]) * double(iter.second);
        rayData.product *= rayData.normalization;
      }

      // Sort rays according to product.

      std::sort(_rays.begin(), _rays.end());

      // Add best rays as long as they have sufficiently positive product.
      double epsilon = objectiveNormalization * _normalizedRayEpsilon;
      for (auto& rayData : _rays)
      {
        T product = 0.0;
        for (const auto& iter : rayData.vector)
          product += double(objectiveVector[iter.first]) * double(iter.second);
        product *= rayData.normalization;

        if (product > epsilon)
        {
          result.rays.push_back(typename OptimizationOracle<T>::Result::Ray(rayData.vector));
          rayData.lastSuccess = _queryCount;
          if (result.rays.size() >= query.maxNumSolutions)
          {
            assert(result.checkPointsSorted());
            return result;
          }
        }
        else
          break;
      }

      // Compute scalar product with each point.

      for (auto& pointData : _points)
      {
        pointData.product = 0.0;
        for (const auto& iter : pointData.vector)
          pointData.product += double(objectiveVector[iter.first]) * double(iter.second);
      }

      // Sort points according to product.

      std::sort(_points.begin(), _points.end());

      // Add best points as longs as they have sufficiently positive product.

      double threshold = query.minObjectiveValue;
      epsilon = objectiveNormalization * _normalizedPointEpsilon;
      for (auto& pointData : _points)
      {
        T product = 0.0;
        for (const auto& iter : pointData.vector)
          product += double(objectiveVector[iter.first]) * double(iter.second);

        if (product <= threshold)
          break;
        result.points.push_back(typename OptimizationOracle<T>::Result::Point(pointData.vector,
          pointData.product));
        pointData.lastSuccess = _queryCount;
      }

      assert(result.checkPointsSorted());
      return result;
    }

    void reduceCacheSize()
    {

    }

  protected:

    struct Data
    {
      double product;
      double hash;
      sparse_vector<T> vector; /// The cached point/ray.
      std::size_t lastSuccess; /// Largest query number at which this vector was returned.
      double normalization; /// Inverse of Euclidean norm of vector.

      Data(sparse_vector<T>&& v, std::size_t ls)
        : vector(std::move(v)), lastSuccess(ls)
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

        return lastSuccess > other.lastSuccess;
      }
    };

    std::vector<Data> _points;
    std::vector<Data> _rays;
    double _normalizedRayEpsilon;
    double _normalizedPointEpsilon;
    std::size_t _queryCount;
  };

  template <typename T, typename IsZero>
  class Polyhedron: public std::enable_shared_from_this<Polyhedron<T, IsZero>>
  {
  public:
    Polyhedron() = delete;

    Polyhedron(std::shared_ptr<Space> space, IsZero isZero)
      : _space(space), _historySize(16)
    {

    }

    Polyhedron(std::shared_ptr<OptimizationOracle<T>> optimizationOracle, IsZero isZero)
      : _historySize(16)
    {
      _space = optimizationOracle->space();
      _cacheOptimization = std::make_shared<CacheOptimizationOracle<T>>(_space);
      _optimization.push_back(Data<OptimizationOracle<T>>(_cacheOptimization, true));
      _optimization.push_back(Data<OptimizationOracle<T>>(optimizationOracle));
    }

    Polyhedron(std::shared_ptr<SeparationOracle<T>> separationOracle)
      : _historySize(16)
    {
      _separation.push_back(Data<SeparationOracle<T>>(separationOracle));
    }

    ~Polyhedron()
    {

    }

    /**
     * \brief Returns the polyhedron's ambient space.
     */

    IPO_EXPORT
    std::shared_ptr<Space> space() const
    {
      return _space;
    }

    /**
     * \brief Maximize a floating-point objective vector.
     *
     * \param objectiveVector Array that maps coordinates to objective value coefficients.
     * \param query Structure for query.
     **/

    typename OptimizationOracle<T>::Result maximizeDouble(const double* objectiveVector,
      const typename OptimizationOracle<T>::Query& query)
    {
      tuneOracles();

      typename OptimizationOracle<T>::Result result;
      for (std::size_t o = 0; o < _optimization.size(); ++o)
      {
        Data<OptimizationOracle<T>>& data = _optimization[o];

        // Run current oracle.

        std::chrono::time_point<std::chrono::system_clock> timeStarted =
          std::chrono::system_clock::now();
        result = data.oracle->maximizeDouble(objectiveVector, query);
        double elapsed = (std::chrono::system_clock::now() - timeStarted).count();

        data.updateHistory(elapsed, !result.isInfeasible(), _historySize);

        // If the current oracle found a ray or a point, we stop.
        if (!result.isInfeasible())
          return result;
      }

      return result;
    }

    /**
     * \brief Maximize an objective vector of the corresponding type.
     *
     * \param objectiveVector Array that maps coordinates to objective value coefficients.
     * \param query Structure for query.
     **/

    IPO_EXPORT
    typename OptimizationOracle<T>::Result maximize(const T* objectiveVector,
      const typename OptimizationOracle<T>::Query& query)
    {
      tuneOracles();

      typename OptimizationOracle<T>::Result result;
      std::size_t lastOracle = 0;
      for (; lastOracle < _optimization.size(); ++lastOracle)
      {
        Data<OptimizationOracle<T>>& data = _optimization[lastOracle];

        // Run current oracle.
        
        std::cout << "Calling oracle <" << data.oracle->name() << ">." << std::endl;

        std::chrono::time_point<std::chrono::system_clock> started = std::chrono::system_clock::now();
        result = data.oracle->maximize(objectiveVector, query);
        std::chrono::duration<double> duration = std::chrono::system_clock::now() - started;

        data.updateHistory(duration.count(), !result.isInfeasible(), _historySize);

        // If the current oracle found a ray or a point, we stop.
        if (!result.isInfeasible())
          return result;
      }

      return result;
    }

  protected:

    IPO_EXPORT
    void tuneOracles()
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
            auto oracle = dynamic_cast<CacheOptimizationOracle<T>&>(
              *_optimization[cacheOracle].oracle);
            oracle.reduceCacheSize();
          }
        }
      }
    }

  protected:
    struct QueryStatistics
    {
      double runningTime;
      bool success;

      QueryStatistics(double rt, bool succ)
        : runningTime(rt), success(succ)
      {

      }
    };

    template <typename O>
    struct Data
    {
      double sumRunningTime; /// Sum of running times of history.
      unsigned int sumSuccess; /// Number of successful runs.
      double priority; /// \ref sumRunningTime / (\ref sumSuccess + 1)
      std::shared_ptr<O> oracle; /// Pointer to oracle.
      std::deque<QueryStatistics> history;
      bool isCache; /// \c true if this is a cache oracle.

      Data(std::shared_ptr<O> o, bool cache = false)
        : sumRunningTime(0.0), sumSuccess(0), priority(0.0), oracle(o), isCache(cache)
      {

      }

      bool operator<(const Data<O>& other) const
      {
        return priority < other.priority;
      }

      void swap(Data<O>& other)
      {
        std::swap(sumRunningTime, other.sumRunningTime);
        std::swap(sumSuccess, other.sumSuccess);
        std::swap(priority, other.priority);
        oracle.swap(other.oracle);
        history.swap(other.history);
        std::swap(isCache, other.isCache);
      }

      void updateHistory(double runningTime, bool success, std::size_t historySize)
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
    };

    std::shared_ptr<Space> _space;
    std::shared_ptr<CacheOptimizationOracle<T>> _cacheOptimization;
    std::vector<Data<OptimizationOracle<T>>> _optimization;
    std::vector<Data<SeparationOracle<T>>> _separation;
    double _normalizedRayEpsilon;
    double _normalizedPointEpsilon;
    std::size_t _historySize;
  };

} /* namespace ipo */

