#include <ipo/polyhedron.hpp>

// #define IPO_DEBUG_POLYHEDRON_PRINT // Uncomment to debug activity.

#include <map>
#include <deque>
#include <random>
#include <chrono>

namespace ipo
{
  template <typename R, typename OptOracle, typename SepaOracle, typename CacheOptOracle>
  struct PolyhedronImplementation;

  class RealCacheOptimizationOracle: public RealOptimizationOracle
  {
  public:
    RealCacheOptimizationOracle(
      PolyhedronImplementation<double, RealOptimizationOracle, RealSeparationOracle, RealCacheOptimizationOracle>* implementation);

    RealOptimizationOracle::Response maximize(const double* objectiveVector,
      const RealOptimizationOracle::Query& query) override;

  protected:
    PolyhedronImplementation<double, RealOptimizationOracle, RealSeparationOracle, RealCacheOptimizationOracle>* _implementation;
    std::size_t _queryCount;
  };

#if defined(IPO_WITH_GMP)

  class RationalCacheOptimizationOracle: public RationalOptimizationOracle
  {
  public:
    RationalCacheOptimizationOracle(
      PolyhedronImplementation<mpq_class, RationalOptimizationOracle, RationalSeparationOracle,
      RationalCacheOptimizationOracle>* implementation);

    RationalOptimizationOracle::Response maximize(const mpq_class* objectiveVector,
      const RationalOptimizationOracle::Query& query) override;

    RationalOptimizationOracle::Response maximize(const double* objectiveVector,
      const RationalOptimizationOracle::Query& query) override;

  protected:
    PolyhedronImplementation<mpq_class, RationalOptimizationOracle, RationalSeparationOracle,
      RationalCacheOptimizationOracle>* _implementation;
    std::size_t _queryCount;
  };
#endif /* IPO_WITH_GMP */

  template <typename R, typename OptOracle, typename SepaOracle, typename CacheOptOracle>
  struct PolyhedronImplementation
  {
    PolyhedronImplementation(std::shared_ptr<Space> space,
      const std::vector<std::shared_ptr<OptOracle>>& optOracles,
      const std::vector<std::shared_ptr<SepaOracle>>& sepaOracles)
      : _space(space)
    {
      for (auto& oracle : optOracles)
        _optimization.push_back(Data<OptOracle>(oracle));
      for (auto& oracle : sepaOracles)
        _separation.push_back(Data<SepaOracle>(oracle));

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
      _squaredDuplicateEpsilon = 1.0e-12;

      _historySize = 256;
      _maxCacheSize = std::numeric_limits<std::size_t>::max();
      _normalizedPointEpsilon = 1.0e-9;
      _normalizedRayEpsilon = 1.0e-9;
    }

    IPO_EXPORT
    typename OptOracle::Response maximize(const R* objectiveVector,
      const typename OptOracle::Query& query)
    {
      sortOptimizationOracles();

      typename OptOracle::Response response;
      std::size_t lastOracle = 0;
      for (; lastOracle < _optimization.size(); ++lastOracle)
      {
        Data<OptOracle>& data = _optimization[lastOracle];

        // Run current oracle.

#if defined (IPO_DEBUG_POLYHEDRON_PRINT)
        std::cout << "\nCalling oracle <" << data.oracle->name() << ">." << std::endl;
#endif /* IPO_DEBUG_POLYHEDRON_PRINT */

        std::chrono::time_point<std::chrono::system_clock> timeStarted = std::chrono::system_clock::now();
        response = data.oracle->maximize(objectiveVector, query);
        double elapsed =  std::chrono::duration<double>(std::chrono::system_clock::now() - timeStarted).count();

        data.updateHistory(elapsed, response.wasSuccessful(), _historySize);

#if defined (IPO_DEBUG_POLYHEDRON_PRINT)
        std::cout << "Oracle " << data.oracle->name() << " took "
          << (data.sumRunningTime / data.history.size()) << "s averaged over " << data.history.size()
          << " queries " << data.sumSuccess << " of which were successful." << std::endl;
#endif /* IPO_DEBUG_POLYHEDRON_PRINT */

        // If the current oracle found a ray or a point, we stop.
        if (response.wasSuccessful())
        {
          if (!data.isCache)
          {
            // We add points / rays to the cache.

            for (auto& ray : response.rays)
              addToCache(_rays, _hashToRayIndex, _rayProducts, ray.vector);

            for (auto& point : response.points)
              addToCache(_points, _hashToPointIndex, _pointProducts, point.vector);
          }

          return response;
        }
      }

      return response;
    }

    /**
     * \brief Maximize a floating-point objective vector.
     *
     * \param objectiveVector Array that maps coordinates to objective value coefficients.
     * \param query Structure for query.
     **/

    typename OptOracle::Response maximizeDouble(const double* objectiveVector,
      const typename OptOracle::Query& query)
    {
      sortOptimizationOracles();

      typename OptOracle::Response response;
      for (std::size_t o = 0; o < _optimization.size(); ++o)
      {
        Data<OptOracle>& data = _optimization[o];

        // Run current oracle.

#if defined (IPO_DEBUG_POLYHEDRON_PRINT)
        std::cout << "\nCalling oracle <" << data.oracle->name() << ">." << std::endl;
#endif /* IPO_DEBUG_POLYHEDRON_PRINT */

        std::chrono::time_point<std::chrono::system_clock> timeStarted =
          std::chrono::system_clock::now();
        response = data.oracle->maximize(objectiveVector, query);
        double elapsed = std::chrono::duration<double>(std::chrono::system_clock::now() - timeStarted).count();

        data.updateHistory(elapsed, response.wasSuccessful(), _historySize);
#if defined (IPO_DEBUG_POLYHEDRON_PRINT)
        std::cout << "Oracle " << data.oracle->name() << " took "
          << (data.sumRunningTime / data.history.size()) << "s averaged over " << data.history.size()
          << " queries " << data.sumSuccess << " of which were successful." << std::endl;
#endif /* IPO_DEBUG_POLYHEDRON_PRINT */

        // If the current oracle found a ray or a point, we stop.
        if (response.wasSuccessful())
          return response;
      }

      return response;
    }

    struct CachedVector
    {
      std::shared_ptr<sparse_vector<R>> vector; /// The cached point/ray.
      double hash; /// The scalar product with the random vector used for hashing.
      std::size_t lastSuccess; /// Largest query number at which this vector was returned.
      double inverseNorm; /// Inverse of Euclidean norm of vector.

      CachedVector(std::shared_ptr<sparse_vector<R>> vec, double hsh)
        : vector(vec), hash(hsh), lastSuccess(0)
      {
        inverseNorm = euclideanNorm(*vec);
        if (inverseNorm > 0.0)
          inverseNorm = 1.0 / inverseNorm;
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

    void addToCache(std::vector<CachedVector>& vectors, std::map<double, std::size_t>& hashToIndex,
      std::vector<ProductVector>& products, std::shared_ptr<sparse_vector<R>>& vector)
    {
      double hash = _hashVector * *vector;
      auto afterIter = hashToIndex.lower_bound(hash);
      double afterDist = std::numeric_limits<double>::infinity();
      if (afterIter != hashToIndex.end())
        afterDist = afterIter->first - hash;
      auto beforeIter = afterIter;
      double beforeDist = std::numeric_limits<double>::infinity();
      if (afterIter != hashToIndex.begin())
      {
        --beforeIter;
        beforeDist = hash - beforeIter->first;
      }
      if (beforeDist < afterDist && squaredEuclideanDistance(
        *vector, *vectors[beforeIter->second].vector) < _squaredDuplicateEpsilon)
      {
        vector = vectors[beforeIter->second].vector;
        // TODO: Recompute objective value.
      }
      else if (afterDist <= beforeDist && std::isfinite(afterDist) && squaredEuclideanDistance(
        *vector, *vectors[afterIter->second].vector) < _squaredDuplicateEpsilon)
      {
        vector = vectors[afterIter->second].vector;
        // TODO: Recompute objective value.
      }
      else
      {
        hashToIndex.insert(afterIter, std::make_pair(hash, vectors.size()));
        products.push_back(ProductVector(vectors.size()));
        vectors.push_back(CachedVector(vector, hash));
      }
    }

    void sortOptimizationOracles()
    {
      // In case the expected successful running times have changed, we reorder the oracles such that
      // the one with the smallest value is queried first.

      std::sort(_optimization.begin(), _optimization.end());

      // The cache oracle should not be the last one.
      if (_optimization.size() > 1 && _optimization.back().isCache)
        std::swap(_optimization[_optimization.size() - 2], _optimization[_optimization.size() - 1]);

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
          // If the expected running time of cache oracle is at least half of the slowest oracle,
          // then we reduce the cache size.

          double cacheRunningTime = _optimization[cacheOracle].sumRunningTime
            / _optimization[cacheOracle].history.size();
          if (cacheRunningTime >= 0.5 * maxRunningTime)
            reduceCacheSize();
        }
      }
    }

  protected:

    void reduceCacheSize()
    {
#if defined(IPO_DEBUG_POLYHEDRON_PRINT)
      std::cout << "We have to reduce the cache size. Currently " << (_points.size() + _rays.size())
        << " of " << _maxCacheSize << std::endl;
#endif /* IPO_DEBUG_POLYHEDRON_PRINT */
    }

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

  public:
    std::shared_ptr<Space> _space;
    std::shared_ptr<CacheOptOracle> _cacheOptimization;
    std::vector<Data<OptOracle>> _optimization;
    std::vector<Data<SepaOracle>> _separation;
    double* _hashVector;
    std::vector<CachedVector> _points;
    std::map<double, std::size_t> _hashToPointIndex;
    std::vector<ProductVector> _pointProducts;
    std::vector<CachedVector> _rays;
    std::map<double, std::size_t> _hashToRayIndex;
    std::vector<ProductVector> _rayProducts;
    /**
     * \brief Square of Euclidean distance above which points/rays are considered as different.
     */
    double _squaredDuplicateEpsilon;
    double _normalizedRayEpsilon;
    double _normalizedPointEpsilon;
    std::size_t _historySize;
    std::size_t _maxCacheSize;
  };
  
  RealCacheOptimizationOracle::RealCacheOptimizationOracle(
    PolyhedronImplementation<double, RealOptimizationOracle, RealSeparationOracle,
    RealCacheOptimizationOracle>* implementation)
    : RealOptimizationOracle("cache"), _implementation(implementation), _queryCount(0)
  {
    this->_space = implementation->_space;
  }

  RealOptimizationOracle::Response RealCacheOptimizationOracle::maximize(
    const double* objectiveVector, const RealOptimizationOracle::Query& query)
  {
    RealOptimizationOracle::Response response;
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
    for (auto& rayProduct : _implementation->_rayProducts)
    {
      rayProduct.product = objectiveVector * *_implementation->_rays[rayProduct.vectorIndex].vector;
      rayProduct.product *= _implementation->_rays[rayProduct.vectorIndex].inverseNorm;
    }

    // Sort rays according to product.

    if (query.maxNumSolutions < _implementation->_rayProducts.size())
    {
      std::nth_element(_implementation->_rayProducts.begin(),
        _implementation->_rayProducts.begin() + query.maxNumSolutions - 1,
        _implementation->_rayProducts.end());
      std::sort(_implementation->_rayProducts.begin(),
        _implementation->_rayProducts.begin() + query.maxNumSolutions - 1);
    }
    else
      std::sort(_implementation->_rayProducts.begin(), _implementation->_rayProducts.end());

    // Add best rays as long as they have sufficiently positive product.
    double epsilon = objectiveNormalization * _implementation->_normalizedRayEpsilon;
    for (auto& rayProduct : _implementation->_rayProducts)
    {
      // Exactly compute the product.
      double product = 0.0;
      for (const auto& iter : *_implementation->_rays[rayProduct.vectorIndex].vector)
        product += objectiveVector[iter.first] * iter.second;
      product *= _implementation->_rays[rayProduct.vectorIndex].inverseNorm;

      if (product > epsilon)
      {
        response.outcome = OptimizationOutcome::UNBOUNDED;
        response.rays.push_back(RealOptimizationOracle::Response::Ray(
          _implementation->_rays[rayProduct.vectorIndex].vector));
        _implementation->_rays[rayProduct.vectorIndex].lastSuccess = _queryCount;
        if (response.rays.size() >= query.maxNumSolutions)
          return response;
      }
      else
      {
        // Since we sorted, all subsequent rays will have a smaller scalar product.
        break;
      }
    }

    if (!response.rays.empty())
      return response;

    // Compute scalar product with each point.

    for (auto& pointProduct : _implementation->_pointProducts)
    {
      pointProduct.product = objectiveVector * *_implementation->_points[pointProduct.vectorIndex].vector;
      pointProduct.product *= _implementation->_points[pointProduct.vectorIndex].inverseNorm;
    }

    // Sort points according to product.

    std::size_t maxNumPoints = query.maxNumSolutions - response.rays.size();
    if (maxNumPoints < _implementation->_pointProducts.size())
    {
      std::nth_element(_implementation->_pointProducts.begin(),
        _implementation->_pointProducts.begin() + maxNumPoints - 1,
        _implementation->_pointProducts.end());
      std::sort(_implementation->_pointProducts.begin(),
        _implementation->_pointProducts.begin() + maxNumPoints  - 1);
    }
    else
      std::sort(_implementation->_pointProducts.begin(), _implementation->_pointProducts.end());

    // Add best points as longs as they have sufficiently positive product.

    double threshold = convertNumber<double>(query.minObjectiveValue);
    if (std::isfinite(threshold))
      threshold += objectiveNormalization * _implementation->_normalizedPointEpsilon;
    for (auto& pointProduct : _implementation->_pointProducts)
    {
      double product = 0.0;
      for (const auto& iter : *_implementation->_points[pointProduct.vectorIndex].vector)
        product += objectiveVector[iter.first] * iter.second;

      if (product <= threshold)
        break;

      response.points.push_back(RealOptimizationOracle::Response::Point(
        _implementation->_points[pointProduct.vectorIndex].vector, pointProduct.product));
      _implementation->_points[pointProduct.vectorIndex].lastSuccess = _queryCount;
      if (response.points.size() >= maxNumPoints)
        break;
    }

    if (!response.points.empty())
    {
      std::sort(response.points.begin(), response.points.end());
      response.primalBound = response.points.front().objectiveValue;
      response.outcome = OptimizationOutcome::FEASIBLE;
    }

    return response;
  }

  typedef PolyhedronImplementation<double, RealOptimizationOracle, RealSeparationOracle,
    RealCacheOptimizationOracle> RealPolyhedronImplementation;

  RealPolyhedron::RealPolyhedron(std::shared_ptr<Space> space)
  {
    _space = space;
    std::vector<std::shared_ptr<RealOptimizationOracle>> optOracles;
    std::vector<std::shared_ptr<RealSeparationOracle>> sepaOracles;
    _implementation = new RealPolyhedronImplementation(_space, optOracles, sepaOracles);
  }

  RealPolyhedron::RealPolyhedron(std::shared_ptr<RealOptimizationOracle> optOracle)
  {
    _space = optOracle->space();
    std::vector<std::shared_ptr<RealOptimizationOracle>> optOracles;
    std::vector<std::shared_ptr<RealSeparationOracle>> sepaOracles;
    optOracles.push_back(optOracle);
    _implementation = new RealPolyhedronImplementation(_space, optOracles, sepaOracles);
  }

  RealPolyhedron::RealPolyhedron(std::shared_ptr<RealSeparationOracle> sepaOracle)
  {
    _space = sepaOracle->space();
    std::vector<std::shared_ptr<RealOptimizationOracle>> optOracles;
    std::vector<std::shared_ptr<RealSeparationOracle>> sepaOracles;
    sepaOracles.push_back(sepaOracle);
    _implementation = new RealPolyhedronImplementation(_space, optOracles, sepaOracles);
  }

  RealPolyhedron::~RealPolyhedron()
  {
    delete static_cast<RealPolyhedronImplementation*>(_implementation);
  }

  RealOptimizationOracle::Response RealPolyhedron::maximize(
    const double* objectiveVector, const RealOptimizationOracle::Query& query)
  {
    return static_cast<RealPolyhedronImplementation*>(_implementation)->maximize(
      objectiveVector, query);
  }

#if defined(IPO_WITH_GMP)

  RationalCacheOptimizationOracle::RationalCacheOptimizationOracle(
    PolyhedronImplementation<mpq_class, RationalOptimizationOracle, RationalSeparationOracle,
    RationalCacheOptimizationOracle>* implementation)
    : RationalOptimizationOracle("cache"), _implementation(implementation), _queryCount(0)
  {
    this->_space = implementation->_space;
  }

  RationalOptimizationOracle::Response RationalCacheOptimizationOracle::maximize(
    const mpq_class* objectiveVector, const RationalOptimizationOracle::Query& query)
  {
    RationalOptimizationOracle::Response response;
    ++_queryCount;

    // Compute norm of objective vector.
    double objectiveNormalization = 0.0;
    for (std::size_t v = 0; v < this->space()->dimension(); ++v)
    {
      double x = convertNumber<double>(objectiveVector[v]);
      objectiveNormalization += x*x;
    }
    objectiveNormalization = sqrt(objectiveNormalization);

    // TODO: Respect time limit.

    // Compute scalar product with each normalized ray.
    for (auto& rayProduct : _implementation->_rayProducts)
    {
      rayProduct.product = 0.0;
      for (const auto& iter : *_implementation->_rays[rayProduct.vectorIndex].vector)
      {
        rayProduct.product += convertNumber<double>(objectiveVector[iter.first])
          * convertNumber<double>(iter.second);
      }
      rayProduct.product *= _implementation->_rays[rayProduct.vectorIndex].inverseNorm;
    }

    // Sort rays according to product.

    if (query.maxNumSolutions < _implementation->_rayProducts.size())
    {
      std::nth_element(_implementation->_rayProducts.begin(),
        _implementation->_rayProducts.begin() + query.maxNumSolutions - 1,
        _implementation->_rayProducts.end());
      std::sort(_implementation->_rayProducts.begin(),
        _implementation->_rayProducts.begin() + query.maxNumSolutions - 1);
    }
    else
      std::sort(_implementation->_rayProducts.begin(), _implementation->_rayProducts.end());

    // Add best rays as long as they have sufficiently positive product.
    double epsilon = objectiveNormalization * _implementation->_normalizedRayEpsilon;
    for (auto& rayProduct : _implementation->_rayProducts)
    {
      // Compute the product exactly.
      mpq_class product = 0;
      for (const auto& iter : *_implementation->_rays[rayProduct.vectorIndex].vector)
        product += objectiveVector[iter.first] * iter.second;
      product *= _implementation->_rays[rayProduct.vectorIndex].inverseNorm;

      if (product > epsilon)
      {
        response.outcome = OptimizationOutcome::UNBOUNDED;
        response.rays.push_back(RationalOptimizationOracle::Response::Ray(
          _implementation->_rays[rayProduct.vectorIndex].vector));
        _implementation->_rays[rayProduct.vectorIndex].lastSuccess = _queryCount;
        if (response.rays.size() >= query.maxNumSolutions)
          return response;
      }
      else
      {
        // Since we sorted, all subsequent rays will have a smaller scalar product.
        break;
      }
    }

    if (!response.rays.empty())
      return response;

    // Compute scalar product with each point.

    for (auto& pointProduct : _implementation->_pointProducts)
    {
      pointProduct.product = 0.0;
      for (const auto& iter : *_implementation->_points[pointProduct.vectorIndex].vector)
      {
        pointProduct.product += convertNumber<double>(objectiveVector[iter.first])
          * convertNumber<double>(iter.second);
      }
      pointProduct.product *= _implementation->_points[pointProduct.vectorIndex].inverseNorm;
    }

    // Sort points according to product.

    std::size_t maxNumPoints = query.maxNumSolutions - response.rays.size();
    if (maxNumPoints < _implementation->_pointProducts.size())
    {
      std::nth_element(_implementation->_pointProducts.begin(),
        _implementation->_pointProducts.begin() + maxNumPoints - 1,
        _implementation->_pointProducts.end());
      std::sort(_implementation->_pointProducts.begin(),
        _implementation->_pointProducts.begin() + maxNumPoints  - 1);
    }
    else
      std::sort(_implementation->_pointProducts.begin(), _implementation->_pointProducts.end());

    // Add best points as longs as they have sufficiently positive product.

    double threshold = convertNumber<double>(query.minObjectiveValue);
    if (std::isfinite(threshold))
      threshold += objectiveNormalization * _implementation->_normalizedPointEpsilon;
    for (auto& pointProduct : _implementation->_pointProducts)
    {
      if (response.points.size() >= maxNumPoints)
        break;
      if (std::isfinite(threshold))
      {
        mpq_class product = 0;
        for (const auto& iter : *_implementation->_points[pointProduct.vectorIndex].vector)
          product += objectiveVector[iter.first] * iter.second;
        if (convertNumber<double>(product) <= threshold)
          break;
      }
      response.points.push_back(RationalOptimizationOracle::Response::Point(
        _implementation->_points[pointProduct.vectorIndex].vector, pointProduct.product));
      _implementation->_points[pointProduct.vectorIndex].lastSuccess = _queryCount;
    }

    if (!response.points.empty())
    {
      response.primalBound = response.points.front().objectiveValue;
      response.outcome = OptimizationOutcome::FEASIBLE;
    }

    return response;
  }

  RationalOptimizationOracle::Response RationalCacheOptimizationOracle::maximize(
    const double* objectiveVector, const RationalOptimizationOracle::Query& query)
  {
    RationalOptimizationOracle::Response response;
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
    for (auto& rayProduct : _implementation->_rayProducts)
    {
      rayProduct.product = objectiveVector * *_implementation->_rays[rayProduct.vectorIndex].vector;
      rayProduct.product *= _implementation->_rays[rayProduct.vectorIndex].inverseNorm;
    }

    // Sort rays according to product.

    if (query.maxNumSolutions < _implementation->_rayProducts.size())
    {
      std::nth_element(_implementation->_rayProducts.begin(),
        _implementation->_rayProducts.begin() + query.maxNumSolutions - 1,
        _implementation->_rayProducts.end());
      std::sort(_implementation->_rayProducts.begin(),
        _implementation->_rayProducts.begin() + query.maxNumSolutions - 1);
    }
    else
      std::sort(_implementation->_rayProducts.begin(), _implementation->_rayProducts.end());

    // Add best rays as long as they have sufficiently positive product.
    double epsilon = objectiveNormalization * _implementation->_normalizedRayEpsilon;
    for (auto& rayProduct : _implementation->_rayProducts)
    {
      // Exactly compute the product.
      mpq_class product = 0;
      for (const auto& iter : *_implementation->_rays[rayProduct.vectorIndex].vector)
        product += objectiveVector[iter.first] * iter.second;
      product *= _implementation->_rays[rayProduct.vectorIndex].inverseNorm;

      if (product > epsilon)
      {
        response.outcome = OptimizationOutcome::UNBOUNDED;
        response.rays.push_back(RationalOptimizationOracle::Response::Ray(
          _implementation->_rays[rayProduct.vectorIndex].vector));
        _implementation->_rays[rayProduct.vectorIndex].lastSuccess = _queryCount;
        if (response.rays.size() >= query.maxNumSolutions)
          return response;
      }
      else
      {
        // Since we sorted, all subsequent rays will have a smaller scalar product.
        break;
      }
    }

    if (!response.rays.empty())
      return response;

    // Compute scalar product with each point.

    for (auto& pointProduct : _implementation->_pointProducts)
    {
      pointProduct.product = objectiveVector * *_implementation->_points[pointProduct.vectorIndex].vector;
      pointProduct.product *= _implementation->_points[pointProduct.vectorIndex].inverseNorm;
    }

    // Sort points according to product.

    std::size_t maxNumPoints = query.maxNumSolutions - response.rays.size();
    if (maxNumPoints < _implementation->_pointProducts.size())
    {
      std::nth_element(_implementation->_pointProducts.begin(),
        _implementation->_pointProducts.begin() + maxNumPoints - 1,
        _implementation->_pointProducts.end());
      std::sort(_implementation->_pointProducts.begin(),
        _implementation->_pointProducts.begin() + maxNumPoints  - 1);
    }
    else
      std::sort(_implementation->_pointProducts.begin(), _implementation->_pointProducts.end());

    // Add best points as longs as they have sufficiently positive product.

    double threshold = convertNumber<double>(query.minObjectiveValue);
    if (std::isfinite(threshold))
      threshold += objectiveNormalization * _implementation->_normalizedPointEpsilon;
    for (auto& pointProduct : _implementation->_pointProducts)
    {
      mpq_class product = 0;
      for (const auto& iter : *_implementation->_points[pointProduct.vectorIndex].vector)
        product += objectiveVector[iter.first] * iter.second;

      if (product <= threshold)
        break;

      response.points.push_back(RationalOptimizationOracle::Response::Point(
        _implementation->_points[pointProduct.vectorIndex].vector, pointProduct.product));
      _implementation->_points[pointProduct.vectorIndex].lastSuccess = _queryCount;
      if (response.points.size() >= maxNumPoints)
        break;
    }

    if (!response.points.empty())
    {
      response.primalBound = response.points.front().objectiveValue;
      response.outcome = OptimizationOutcome::FEASIBLE;
    }

    return response;
  }

  typedef PolyhedronImplementation<mpq_class, RationalOptimizationOracle, RationalSeparationOracle,
    RationalCacheOptimizationOracle> RationalPolyhedronImplementation;

  RationalPolyhedron::RationalPolyhedron(std::shared_ptr<Space> space)
  {
    _space = space;
    std::vector<std::shared_ptr<RationalOptimizationOracle>> optOracles;
    std::vector<std::shared_ptr<RationalSeparationOracle>> sepaOracles;
    _implementation = new RationalPolyhedronImplementation(_space, optOracles, sepaOracles);
  }

  RationalPolyhedron::RationalPolyhedron(std::shared_ptr<RationalOptimizationOracle> optOracle)
  {
    _space = optOracle->space();
    std::vector<std::shared_ptr<RationalOptimizationOracle>> optOracles;
    std::vector<std::shared_ptr<RationalSeparationOracle>> sepaOracles;
    optOracles.push_back(optOracle);
    _implementation = new RationalPolyhedronImplementation(_space, optOracles, sepaOracles);
  }

  RationalPolyhedron::RationalPolyhedron(std::shared_ptr<RationalSeparationOracle> sepaOracle)
  {
    _space = sepaOracle->space();
    std::vector<std::shared_ptr<RationalOptimizationOracle>> optOracles;
    std::vector<std::shared_ptr<RationalSeparationOracle>> sepaOracles;
    sepaOracles.push_back(sepaOracle);
    _implementation = new RationalPolyhedronImplementation(_space, optOracles, sepaOracles);
  }

  RationalPolyhedron::~RationalPolyhedron()
  {
    delete static_cast<RationalPolyhedronImplementation*>(_implementation);
  }

  RationalOptimizationOracle::Response RationalPolyhedron::maximize(
    const mpq_class* objectiveVector, const RationalOptimizationOracle::Query& query)
  {
    return static_cast<RationalPolyhedronImplementation*>(_implementation)->maximize(
      objectiveVector, query);
  }

  RationalOptimizationOracle::Response RationalPolyhedron::maximize(const double* objectiveVector,
    const RationalOptimizationOracle::Query& query)
  {
    return static_cast<RationalPolyhedronImplementation*>(_implementation)->maximizeDouble(
      objectiveVector, query);
  }

#endif /* IPO_WITH_GMP */


//   template <typename T>
//   class Polyhedron: public std::enable_shared_from_this<Polyhedron<T>>
//   {
//   public:
//     Polyhedron() = delete;
// 
//     Polyhedron(std::shared_ptr<Space> space)
//       : _space(space)
//     {
//       initialize();
//     }
// 
//     Polyhedron(std::shared_ptr<OptimizationOracle<T>> optimizationOracle)
//     {
//       _space = optimizationOracle->space();
//       _cacheOptimization = std::make_shared<CacheOptimizationOracle<T>>(this);
//       _optimization.push_back(Data<OptimizationOracle<T>>(_cacheOptimization, true));
//       _optimization.push_back(Data<OptimizationOracle<T>>(optimizationOracle));
//       initialize();
//     }
// 
//     Polyhedron(std::shared_ptr<SeparationOracle<T>> separationOracle)
//     {
//       _separation.push_back(Data<SeparationOracle<T>>(separationOracle));
//       _space = separationOracle->space();
//       initialize();
//     }
// 
//     void initialize()
//     {
//       // To find duplicates of points or rays, we need a random vector.
// 
//       std::default_random_engine generator;
//       std::normal_distribution<double> distribution;
//       _hashVector = new double[_space->dimension()];
//       double squaredNorm = 0.0;
//       while (squaredNorm < 1.0e-3)
//       {
//         for (std::size_t v = 0; v < _space->dimension(); ++v)
//         {
//           double x = distribution(generator);
//           _hashVector[v] = x;
//           squaredNorm += x*x;
//         }
//       }
//       for (std::size_t v = 0; v < _space->dimension(); ++v)
//         _hashVector[v] /= squaredNorm;
//       _squaredDuplicateEpsilon = 1.0e-12;
// 
//       _historySize = 256;
//       _maxCacheSize = std::numeric_limits<std::size_t>::max();
//       _normalizedPointEpsilon = 1.0e-9;
//       _normalizedRayEpsilon = 1.0e-9;
//     }
// 
//     ~Polyhedron()
//     {
//       delete[] _hashVector;
//     }
// 
//     /**
//      * \brief Returns the polyhedron's ambient space.
//      */
// 
//     IPO_EXPORT
//     std::shared_ptr<Space> space() const
//     {
//       return _space;
//     }
// 
//     /**
//      * \brief Maximize a floating-point objective vector.
//      *
//      * \param objectiveVector Array that maps coordinates to objective value coefficients.
//      * \param query Structure for query.
//      **/
// 
//     typename OptimizationOracle<T>::Response maximizeDouble(const double* objectiveVector,
//       const typename OptimizationOracle<T>::Query& query)
//     {
//       sortOptimizationOracles();
// 
//       typename OptimizationOracle<T>::Response response;
//       for (std::size_t o = 0; o < _optimization.size(); ++o)
//       {
//         Data<OptimizationOracle<T>>& data = _optimization[o];
// 
//         // Run current oracle.
// 
// #if defined (IPO_DEBUG_POLYHEDRON_PRINT)
//         std::cout << "\nCalling oracle <" << data.oracle->name() << ">." << std::endl;
// #endif /* IPO_DEBUG_POLYHEDRON_PRINT */
// 
//         std::chrono::time_point<std::chrono::system_clock> timeStarted =
//           std::chrono::system_clock::now();
//         response = data.oracle->maximizeDouble(objectiveVector, query);
//         double elapsed = std::chrono::duration<double>(std::chrono::system_clock::now() - timeStarted).count();
// 
//         data.updateHistory(elapsed, response.wasSuccessful(), _historySize);
// #if defined (IPO_DEBUG_POLYHEDRON_PRINT)
//         std::cout << "Oracle " << data.oracle->name() << " took "
//           << (data.sumRunningTime / data.history.size()) << "s averaged over " << data.history.size()
//           << " queries " << data.sumSuccess << " of which were successful." << std::endl;
// #endif /* IPO_DEBUG_POLYHEDRON_PRINT */
// 
//         // If the current oracle found a ray or a point, we stop.
//         if (response.wasSuccessful())
//           return response;
//       }
// 
//       return response;
//     }
// 
//     /**
//      * \brief Maximize an objective vector of the corresponding type.
//      *
//      * \param objectiveVector Array that maps coordinates to objective value coefficients.
//      * \param query Structure for query.
//      **/
// 
//     IPO_EXPORT
//     typename OptimizationOracle<T>::Response maximize(const T* objectiveVector,
//       const typename OptimizationOracle<T>::Query& query)
//     {
//       sortOptimizationOracles();
// 
//       typename OptimizationOracle<T>::Response response;
//       std::size_t lastOracle = 0;
//       for (; lastOracle < _optimization.size(); ++lastOracle)
//       {
//         Data<OptimizationOracle<T>>& data = _optimization[lastOracle];
// 
//         // Run current oracle.
// 
// #if defined (IPO_DEBUG_POLYHEDRON_PRINT)
//         std::cout << "\nCalling oracle <" << data.oracle->name() << ">." << std::endl;
// #endif /* IPO_DEBUG_POLYHEDRON_PRINT */
// 
//         std::chrono::time_point<std::chrono::system_clock> timeStarted = std::chrono::system_clock::now();
//         response = data.oracle->maximize(objectiveVector, query);
//         double elapsed =  std::chrono::duration<double>(std::chrono::system_clock::now() - timeStarted).count();
// 
//         data.updateHistory(elapsed, response.wasSuccessful(), _historySize);
// 
// #if defined (IPO_DEBUG_POLYHEDRON_PRINT)
//         std::cout << "Oracle " << data.oracle->name() << " took "
//           << (data.sumRunningTime / data.history.size()) << "s averaged over " << data.history.size()
//           << " queries " << data.sumSuccess << " of which were successful." << std::endl;
// #endif /* IPO_DEBUG_POLYHEDRON_PRINT */
// 
//         // If the current oracle found a ray or a point, we stop.
//         if (response.wasSuccessful())
//         {
//           if (!data.isCache)
//           {
//             // We add points / rays to the cache.
// 
//             for (auto& ray : response.rays)
//               addToCache(_rays, _hashToRayIndex, _rayProducts, ray.vector);
// 
//             for (auto& point : response.points)
//               addToCache(_points, _hashToPointIndex, _pointProducts, point.vector);
//           }
// 
//           return response;
//         }
//       }
// 
//       return response;
//     }
// 
//   protected:
// 
//     struct CachedVector
//     {
//       std::shared_ptr<sparse_vector<T>> vector; /// The cached point/ray.
//       double hash; /// The scalar product with the random vector used for hashing.
//       std::size_t lastSuccess; /// Largest query number at which this vector was returned.
//       double inverseNorm; /// Inverse of Euclidean norm of vector.
// 
//       CachedVector(std::shared_ptr<sparse_vector<T>> vec, double hsh)
//         : vector(vec), hash(hsh), lastSuccess(0)
//       {
//         inverseNorm = euclideanNorm(*vec);
//         if (inverseNorm > 0.0)
//           inverseNorm = 1.0 / inverseNorm;
//       }
//     };
// 
//     struct ProductVector
//     {
//       std::size_t vectorIndex;
//       double product;
// 
//       ProductVector(std::size_t index)
//         : vectorIndex(index), product(0)
//       {
// 
//       }
// 
//       bool operator<(const ProductVector& other) const
//       {
//         return product > other.product;
//       }
//     };
// 
//     void addToCache(std::vector<CachedVector>& vectors, std::map<double, std::size_t>& hashToIndex,
//       std::vector<ProductVector>& products, std::shared_ptr<sparse_vector<T>>& vector)
//     {
//       double hash = _hashVector * *vector;
//       auto afterIter = hashToIndex.lower_bound(hash);
//       double afterDist = std::numeric_limits<double>::infinity();
//       if (afterIter != hashToIndex.end())
//         afterDist = afterIter->first - hash;
//       auto beforeIter = afterIter;
//       double beforeDist = std::numeric_limits<double>::infinity();
//       if (afterIter != hashToIndex.begin())
//       {
//         --beforeIter;
//         beforeDist = hash - beforeIter->first;
//       }
//       if (beforeDist < afterDist && squaredEuclideanDistance(
//         *vector, *vectors[beforeIter->second].vector) < _squaredDuplicateEpsilon)
//       {
//         vector = vectors[beforeIter->second].vector;
//         // TODO: Recompute objective value.
//       }
//       else if (afterDist <= beforeDist && std::isfinite(afterDist) && squaredEuclideanDistance(
//         *vector, *vectors[afterIter->second].vector) < _squaredDuplicateEpsilon)
//       {
//         vector = vectors[afterIter->second].vector;
//         // TODO: Recompute objective value.
//       }
//       else
//       {
//         hashToIndex.insert(afterIter, std::make_pair(hash, vectors.size()));
//         products.push_back(ProductVector(vectors.size()));
//         vectors.push_back(CachedVector(vector, hash));
//       }
//     }
// 
//     void sortOptimizationOracles()
//     {
//       // In case the expected successful running times have changed, we reorder the oracles such that
//       // the one with the smallest value is queried first.
// 
//       std::sort(_optimization.begin(), _optimization.end());
// 
//       // The cache oracle should not be the last one.
//       if (_optimization.size() > 1 && _optimization.back().isCache)
//         std::swap(_optimization[_optimization.size() - 2], _optimization[_optimization.size() - 1]);
// 
//       if (_optimization.size() > 1)
//       {
//         // Maybe adjust allowed size of cache oracle.
//         std::size_t cacheOracle = std::numeric_limits<std::size_t>::max();
//         double maxRunningTime = 0.0;
//         for (std::size_t o = 0; o < _optimization.size(); ++o)
//         {
//           if (_optimization[o].isCache)
//           {
//             cacheOracle = o;
//             continue;
//           }
// 
//           if (!_optimization[o].history.empty())
//           {
//             double runningTime = _optimization[o].sumRunningTime / _optimization[o].history.size();
//             maxRunningTime = std::max(maxRunningTime, runningTime);
//           }
//         }
// 
//         if (cacheOracle < std::numeric_limits<std::size_t>::max())
//         {
//           // If the expected running time of cache oracle is at least half of the slowest oracle,
//           // then we reduce the cache size.
// 
//           double cacheRunningTime = _optimization[cacheOracle].sumRunningTime
//             / _optimization[cacheOracle].history.size();
//           if (cacheRunningTime >= 0.5 * maxRunningTime)
//             reduceCacheSize();
//         }
//       }
//     }
// 
//   protected:
// 
//     void reduceCacheSize()
//     {
// #if defined(IPO_DEBUG_POLYHEDRON_PRINT)
//       std::cout << "We have to reduce the cache size. Currently " << (_points.size() + _rays.size())
//         << " of " << _maxCacheSize << std::endl;
// #endif /* IPO_DEBUG_POLYHEDRON_PRINT */
//     }
// 
//     struct QueryStatistics
//     {
//       double runningTime;
//       bool success;
// 
//       QueryStatistics(double rt, bool succ)
//         : runningTime(rt), success(succ)
//       {
// 
//       }
//     };
// 
//     template <typename O>
//     struct Data
//     {
//       double sumRunningTime; /// Sum of running times of history.
//       unsigned int sumSuccess; /// Number of successful runs.
//       double priority; /// \ref sumRunningTime / (\ref sumSuccess + 1)
//       std::shared_ptr<O> oracle; /// Pointer to oracle.
//       std::deque<QueryStatistics> history;
//       bool isCache; /// \c true if this is a cache oracle.
// 
//       Data(std::shared_ptr<O> o, bool cache = false)
//         : sumRunningTime(0.0), sumSuccess(0), priority(0.0), oracle(o), isCache(cache)
//       {
// 
//       }
// 
//       bool operator<(const Data<O>& other) const
//       {
//         return priority < other.priority;
//       }
// 
//       void swap(Data<O>& other)
//       {
//         std::swap(sumRunningTime, other.sumRunningTime);
//         std::swap(sumSuccess, other.sumSuccess);
//         std::swap(priority, other.priority);
//         oracle.swap(other.oracle);
//         history.swap(other.history);
//         std::swap(isCache, other.isCache);
//       }
// 
//       void updateHistory(double runningTime, bool success, std::size_t historySize)
//       {
//         history.push_back(QueryStatistics(runningTime, success));
//         sumRunningTime += runningTime;
//         if (success)
//           ++sumSuccess;
// 
//         // If history is too long, remove the first entry.
// 
//         if (history.size() > historySize)
//         {
//           sumRunningTime -= history.front().runningTime;
//           if (history.front().success)
//             --sumSuccess;
//           history.pop_front();
//         }
// 
//         // We always add a bonus of +1 to the number of successful runs. This avoids division by 0
//         // and gives some advantage to unsuccessful oracles that do not require too much time.
// 
//         priority = sumRunningTime / (sumSuccess + 1);
//       }
//     };
// 
//     friend CacheOptimizationOracle<T>;
// 
//     std::shared_ptr<Space> _space;
//     std::shared_ptr<CacheOptimizationOracle<T>> _cacheOptimization;
//     std::vector<Data<OptimizationOracle<T>>> _optimization;
//     std::vector<Data<SeparationOracle<T>>> _separation;
//     double* _hashVector;
//     std::vector<CachedVector> _points;
//     std::map<double, std::size_t> _hashToPointIndex;
//     std::vector<ProductVector> _pointProducts;
//     std::vector<CachedVector> _rays;
//     std::map<double, std::size_t> _hashToRayIndex;
//     std::vector<ProductVector> _rayProducts;
//     /**
//      * \brief Square of Euclidean distance above which points/rays are considered as different.
//      */
//     double _squaredDuplicateEpsilon;
//     double _normalizedRayEpsilon;
//     double _normalizedPointEpsilon;
//     std::size_t _historySize;
//     std::size_t _maxCacheSize;
//   };

} /* namespace ipo */
