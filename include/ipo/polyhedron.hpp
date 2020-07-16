#pragma once

#include <iostream>

#include <ipo/config.hpp>
#include <ipo/export.hpp>
#include <ipo/space.hpp>
#include <ipo/oracles.hpp>

#include <map>
#include <chrono>
#include <memory>
#include <deque>
#include <random>

namespace ipo
{
  template <typename T>
  class Polyhedron;

  template <typename T>
  class CacheOptimizationOracle: public OptimizationOracle<T>
  {
  public:
    CacheOptimizationOracle(Polyhedron<T>* polyhedron)
      : OptimizationOracle<T>("CacheOptimizationOracle"), _polyhedron(polyhedron), _queryCount(0)
    {
      this->_space = polyhedron->space();
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
      for (auto& rayProduct : _polyhedron->_rayProducts)
      {
        rayProduct.product = 0.0;
        for (const auto& iter : *_polyhedron->_rays[rayProduct.vectorIndex].vector)
          rayProduct.product += objectiveVector[iter.first] * double(iter.second);
        rayProduct.product *= _polyhedron->_rays[rayProduct.vectorIndex].inverseNorm;
      }

      // Sort rays according to product.

      if (query.maxNumSolutions < _polyhedron->_rayProducts.size())
      {
        std::nth_element(_polyhedron->_rayProducts.begin(),
          _polyhedron->_rayProducts.begin() + query.maxNumSolutions - 1,
          _polyhedron->_rayProducts.end());
        std::sort(_polyhedron->_rayProducts.begin(),
          _polyhedron->_rayProducts.begin() + query.maxNumSolutions - 1);
      }
      else
        std::sort(_polyhedron->_rayProducts.begin(), _polyhedron->_rayProducts.end());

      // Add best rays as long as they have sufficiently positive product.
      double epsilon = objectiveNormalization * _polyhedron->_normalizedRayEpsilon;
      for (auto& rayProduct : _polyhedron->_rayProducts)
      {
        // Exactly compute the product.
        T product = 0.0;
        for (const auto& iter : *_polyhedron->_rays[rayProduct.vectorIndex].vector)
          product += objectiveVector[iter.first] * iter.second;
        product *= _polyhedron->_rays[rayProduct.vectorIndex].inverseNorm;

        if (product > epsilon)
        {
          result.rays.push_back(typename OptimizationOracle<T>::Result::Ray(
            _polyhedron->_rays[rayProduct.vectorIndex].vector));
          _polyhedron->_rays[rayProduct.vectorIndex].lastSuccess = _queryCount;
          if (result.rays.size() >= query.maxNumSolutions)
            return result;
        }
        else
        {
          // Since we sorted, all subsequent rays will have a smaller scalar product.
          break;
        }
      }

      if (!result.rays.empty())
      {
        result.primalBound = std::numeric_limits<double>::infinity();
        result.dualBound = std::numeric_limits<double>::signaling_NaN();
        return result;
      }

      // Compute scalar product with each point.

      for (auto& pointProduct : _polyhedron->_pointProducts)
      {
        pointProduct.product = 0.0;
        for (const auto& iter : *_polyhedron->_points[pointProduct.vectorIndex].vector)
          pointProduct.product += objectiveVector[iter.first] * double(iter.second);
        pointProduct.product *= _polyhedron->_points[pointProduct.vectorIndex].inverseNorm;
      }

      // Sort points according to product.

      std::size_t maxNumPoints = query.maxNumSolutions - result.rays.size();
      if (maxNumPoints < _polyhedron->_pointProducts.size())
      {
        std::nth_element(_polyhedron->_pointProducts.begin(),
          _polyhedron->_pointProducts.begin() + maxNumPoints - 1,
          _polyhedron->_pointProducts.end());
        std::sort(_polyhedron->_pointProducts.begin(),
          _polyhedron->_pointProducts.begin() + maxNumPoints  - 1);
      }
      else
        std::sort(_polyhedron->_pointProducts.begin(), _polyhedron->_pointProducts.end());

      // Add best points as longs as they have sufficiently positive product.

      double threshold = query.minObjectiveValue;
      if (std::isfinite(threshold))
        threshold += objectiveNormalization * _polyhedron->_normalizedPointEpsilon;
      for (auto& pointProduct : _polyhedron->_pointProducts)
      {
        T product = 0.0;
        for (const auto& iter : *_polyhedron->_points[pointProduct.vectorIndex].vector)
          product += objectiveVector[iter.first] * iter.second;

        if (product <= threshold)
          break;

        result.points.push_back(typename OptimizationOracle<T>::Result::Point(
          _polyhedron->_points[pointProduct.vectorIndex].vector, pointProduct.product));
        _polyhedron->_points[pointProduct.vectorIndex].lastSuccess = _queryCount;
        if (result.points.size() >= maxNumPoints)
          break;
      }

      if (!result.points.empty())
      {
        result.sortPoints();
        result.primalBound = result.points.front().objectiveValue;
      }

      return result;
    }

    virtual typename OptimizationOracle<T>::Result maximize(const T* objectiveVector,
      const typename OptimizationOracle<T>::Query& query) override
    {
      typename OptimizationOracle<T>::Result result;
      return result;

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
      for (auto& rayProduct : _polyhedron->_rayProducts)
      {
        rayProduct.product = 0.0;
        for (const auto& iter : *_polyhedron->_rays[rayProduct.vectorIndex].vector)
          rayProduct.product += double(objectiveVector[iter.first]) * double(iter.second);
        rayProduct.product *= _polyhedron->_rays[rayProduct.vectorIndex].inverseNorm;
      }

      // Sort rays according to product.

      if (query.maxNumSolutions < _polyhedron->_rayProducts.size())
      {
        std::nth_element(_polyhedron->_rayProducts.begin(),
          _polyhedron->_rayProducts.begin() + query.maxNumSolutions - 1,
          _polyhedron->_rayProducts.end());
        std::sort(_polyhedron->_rayProducts.begin(),
          _polyhedron->_rayProducts.begin() + query.maxNumSolutions - 1);
      }
      else
        std::sort(_polyhedron->_rayProducts.begin(), _polyhedron->_rayProducts.end());

      // Add best rays as long as they have sufficiently positive product.
      double epsilon = objectiveNormalization * _polyhedron->_normalizedRayEpsilon;
      for (auto& rayProduct : _polyhedron->_rayProducts)
      {
        // Compute the product exactly.
        T product = 0.0;
        for (const auto& iter : *_polyhedron->_rays[rayProduct.vectorIndex].vector)
          product += objectiveVector[iter.first] * iter.second;
        product *= _polyhedron->_rays[rayProduct.vectorIndex].inverseNorm;

        if (product > epsilon)
        {
          result.rays.push_back(typename OptimizationOracle<T>::Result::Ray(
            _polyhedron->_rays[rayProduct.vectorIndex].vector));
          _polyhedron->_rays[rayProduct.vectorIndex].lastSuccess = _queryCount;
          if (result.rays.size() >= query.maxNumSolutions)
            return result;
        }
        else
        {
          // Since we sorted, all subsequent rays will have a smaller scalar product.
          break;
        }
      }

      if (!result.rays.empty())
      {
        result.primalBound = std::numeric_limits<double>::infinity();
        result.dualBound = std::numeric_limits<double>::signaling_NaN();
        return result;
      }

      // Compute scalar product with each point.

      for (auto& pointProduct : _polyhedron->_pointProducts)
      {
        pointProduct.product = 0.0;
        for (const auto& iter : *_polyhedron->_points[pointProduct.vectorIndex].vector)
          pointProduct.product += double(objectiveVector[iter.first]) * double(iter.second);
        pointProduct.product *= _polyhedron->_points[pointProduct.vectorIndex].inverseNorm;
      }

      // Sort points according to product.

      std::size_t maxNumPoints = query.maxNumSolutions - result.rays.size();
      if (maxNumPoints < _polyhedron->_pointProducts.size())
      {
        std::nth_element(_polyhedron->_pointProducts.begin(),
          _polyhedron->_pointProducts.begin() + maxNumPoints - 1,
          _polyhedron->_pointProducts.end());
        std::sort(_polyhedron->_pointProducts.begin(),
          _polyhedron->_pointProducts.begin() + maxNumPoints  - 1);
      }
      else
        std::sort(_polyhedron->_pointProducts.begin(), _polyhedron->_pointProducts.end());

      // Add best points as longs as they have sufficiently positive product.

      double threshold = double(query.minObjectiveValue);
      if (std::isfinite(threshold))
        threshold += objectiveNormalization * _polyhedron->_normalizedPointEpsilon;
      for (auto& pointProduct : _polyhedron->_pointProducts)
      {
        T product = 0.0;
        for (const auto& iter : *_polyhedron->_points[pointProduct.vectorIndex].vector)
          product += objectiveVector[iter.first] * iter.second;

        if (product <= threshold || result.points.size() >= maxNumPoints)
          break;
        result.points.push_back(typename OptimizationOracle<T>::Result::Point(
          _polyhedron->_points[pointProduct.vectorIndex].vector, pointProduct.product));
        _polyhedron->_points[pointProduct.vectorIndex].lastSuccess = _queryCount;
      }

      assert(result.checkPointsSorted());
      return result;
    }


  protected:
    Polyhedron<T>* _polyhedron;
    std::size_t _queryCount;
  };

  template <typename T>
  class Polyhedron: public std::enable_shared_from_this<Polyhedron<T>>
  {
  public:
    Polyhedron() = delete;

    Polyhedron(std::shared_ptr<Space> space)
      : _space(space), _historySize(16)
    {
      initialize();
    }

    Polyhedron(std::shared_ptr<OptimizationOracle<T>> optimizationOracle)
      : _historySize(16)
    {
      _space = optimizationOracle->space();
      _cacheOptimization = std::make_shared<CacheOptimizationOracle<T>>(this);
      _optimization.push_back(Data<OptimizationOracle<T>>(_cacheOptimization, true));
      _optimization.push_back(Data<OptimizationOracle<T>>(optimizationOracle));
      initialize();
    }

    Polyhedron(std::shared_ptr<SeparationOracle<T>> separationOracle)
      : _historySize(16)
    {
      _separation.push_back(Data<SeparationOracle<T>>(separationOracle));
      initialize();
    }

    void initialize()
    {
      // To find duplicates of points or rays, we need a random vector.

      std::default_random_engine generator;
      std::normal_distribution<double> distribution;
      _hashVector = new double[_space->dimension()];
      double squaredNorm = 0.0;
      while (squaredNorm < 1.0e-3)
      {
        for (std::size_t v = 0; v < _space->dimension(); ++v)
        {
          double x = distribution(generator);
          _hashVector[v] = x;
          squaredNorm += x*x;
        }
      }
      for (std::size_t v = 0; v < _space->dimension(); ++v)
        _hashVector[v] /= squaredNorm;
      _hashEpsilon = 1.0e-9;
    }

    ~Polyhedron()
    {
      delete[] _hashVector;
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
      sortOptimizationOracles();

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
      sortOptimizationOracles();

      typename OptimizationOracle<T>::Result result;
      std::size_t lastOracle = 0;
      for (; lastOracle < _optimization.size(); ++lastOracle)
      {
        Data<OptimizationOracle<T>>& data = _optimization[lastOracle];

        // Run current oracle.

        std::cout << "\nCalling oracle <" << data.oracle->name() << ">." << std::endl;

        std::chrono::time_point<std::chrono::system_clock> started = std::chrono::system_clock::now();
        result = data.oracle->maximize(objectiveVector, query);
        std::chrono::duration<double> duration = std::chrono::system_clock::now() - started;

        data.updateHistory(duration.count(), !result.isInfeasible(), _historySize);

        // If the current oracle found a ray or a point, we stop.
        if (result.isUnbounded() || result.isFeasible() || result.isInfeasible())
        {
          if (!data.isCache)
          {
            // We add all points to the cache.

            for (auto& point : result.points)
            {
              std::shared_ptr<sparse_vector<T>> vector = point.vector;
              
              std::cout << "Checking " << *vector << " for duplicate." << std::endl;
              
              double inverseNorm = 0.0;
              for (auto& iter : *vector)
              {
                double x = iter.second;
                inverseNorm += x*x;
              }
              if (inverseNorm > 0.0)
                inverseNorm = 1.0 / sqrt(inverseNorm);
              else
                inverseNorm = 1.0;
              double hash = *vector * _hashVector;
              auto afterIter = _hashToPointIndex.lower_bound(hash);
              double afterDist = std::numeric_limits<double>::infinity();
              if (afterIter != _hashToPointIndex.end())
                afterDist = afterIter->first - hash;
              auto beforeIter = afterIter;
              double beforeDist = std::numeric_limits<double>::infinity();
              if (afterIter != _hashToPointIndex.begin())
              {
                --beforeIter;
                beforeDist = hash - beforeIter->first;
              }
              if (beforeDist < afterDist && beforeDist * inverseNorm <= _hashEpsilon)
              {
                std::cout << "The before iterator is close enough with distance " << beforeDist << "." << std::endl;
                std::cout << "Replacing " << *vector << " by " << *_points[beforeIter->second].vector << std::endl;
                point.vector = _points[beforeIter->second].vector;
                // TODO: Recompute objective value.
              }
              else if (afterDist <= beforeDist && afterDist * inverseNorm <= _hashEpsilon)
              {
                std::cout << "The after iterator is close enough with distance " << afterDist << "." << std::endl;
                std::cout << "Replacing " << *vector << " by " << *_points[afterIter->second].vector << std::endl;
                point.vector = _points[afterIter->second].vector;
                // TODO: Recompute objective value.
              }
              else
              {
                std::cout << "No iterator is close enough." << std::endl;
                _hashToPointIndex.insert(afterIter, std::make_pair(hash, _points.size()));
                _pointProducts.push_back(ProductVector(_points.size()));
                _points.push_back(CachedVector(point.vector, hash, inverseNorm));
              }
            }
          }

          return result;
        }
      }

      return result;
    }

  protected:

    IPO_EXPORT
    void sortOptimizationOracles()
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
//             oracle.reduceCacheSize();
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

    struct CachedVector
    {
      std::shared_ptr<sparse_vector<T>> vector; /// The cached point/ray.
      double hash; /// The scalar product with the random vector used for hashing.
      std::size_t lastSuccess; /// Largest query number at which this vector was returned.
      double inverseNorm; /// Inverse of Euclidean norm of vector.

      CachedVector(std::shared_ptr<sparse_vector<T>> vec, double hsh, double inorm)
        : vector(vec), hash(hsh), lastSuccess(0), inverseNorm(inorm)
      {

      }
    };

    struct ProductVector
    {
      std::size_t vectorIndex;
      double product;

      ProductVector(std::size_t index)
        : vectorIndex(index), product(0)
      {

      }

      bool operator<(const ProductVector& other) const
      {
        return product > other.product;
      }
    };

    friend CacheOptimizationOracle<T>;

    std::shared_ptr<Space> _space;
    std::shared_ptr<CacheOptimizationOracle<T>> _cacheOptimization;
    std::vector<Data<OptimizationOracle<T>>> _optimization;
    std::vector<Data<SeparationOracle<T>>> _separation;
    double* _hashVector;
    std::vector<CachedVector> _points;
    std::map<double, std::size_t> _hashToPointIndex;
    std::vector<ProductVector> _pointProducts;
    std::vector<CachedVector> _rays;
    std::map<double, std::size_t> _hashToRayIndex;
    std::vector<ProductVector> _rayProducts;
    double _hashEpsilon;
    double _normalizedRayEpsilon;
    double _normalizedPointEpsilon;
    std::size_t _historySize;
  };

} /* namespace ipo */

