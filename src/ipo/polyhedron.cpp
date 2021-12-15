// #define IPO_DEBUG /* Uncomment to debug this file. */

#include <ipo/polyhedron.hpp>

#include <map>
#include <deque>
#include <random>
#include <chrono>

namespace ipo
{
  template <typename NumberType>
  struct PolyhedronImplementation;

  template <typename NumberType>
  class CacheOptimizationOracle: public OptimizationOracle<NumberType>
  {
  public:
    typedef NumberType Number;

    CacheOptimizationOracle(PolyhedronImplementation<Number>* implementation)
      : OptimizationOracle<NumberType>("cache"), _implementation(implementation), _queryCount(0)
    {
      this->_space = implementation->space();
    }

    virtual OptimizationResponse<Number> maximize(const Number* objectiveVector,
      const OptimizationQuery<Number>& query) override
    {
      typename OptimizationOracle<Number>::Response response;
      ++_queryCount;

#if defined(IPO_DEBUG)
      std::cout << "CacheOptimizationOracle::maximize() call #" << _queryCount
        << ", searching among " << _implementation->_pointProducts.size() << " points and "
        << _implementation->_rayProducts.size() << " rays." << std::endl;
#endif /* IPO_DEBUG */

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
        rayProduct.product = convertNumber<double>(objectiveVector * *_implementation->_rays[rayProduct.vectorIndex].vector);
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
        Number product = 0;
        for (const auto& iter : *_implementation->_rays[rayProduct.vectorIndex].vector)
          product += objectiveVector[iter.first] * iter.second;
        product *= _implementation->_rays[rayProduct.vectorIndex].inverseNorm;

        if (product > epsilon)
        {
          response.outcome = OptimizationOutcome::UNBOUNDED;
          response.rays.push_back(typename OptimizationResponse<Number>::Ray(
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
        pointProduct.product = convertNumber<double>(objectiveVector * *_implementation->_points[pointProduct.vectorIndex].vector);
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

      for (auto& pointProduct : _implementation->_pointProducts)
      {
        Number product = 0;
        for (const auto& iter : *_implementation->_points[pointProduct.vectorIndex].vector)
          product += objectiveVector[iter.first] * iter.second;

#if defined(IPO_DEBUG)
        std::cout << "Product with cached point is " << convertNumber<double>(product) << std::endl;
#endif /* IPO_DEBUG */

        if (query.hasMinPrimalBound() && (product <= query.minPrimalBound()))
          break;

        response.points.push_back(typename OptimizationResponse<Number>::Point(
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

    virtual OptimizationResponse<Number> maximizeDouble(
      const double* objectiveVector, const OptimizationQuery<Number>& query) override
    {
      typename OptimizationOracle<Number>::Response response;
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
        Number product = 0;
        for (const auto& iter : *_implementation->_rays[rayProduct.vectorIndex].vector)
          product += objectiveVector[iter.first] * iter.second;
        product *= _implementation->_rays[rayProduct.vectorIndex].inverseNorm;

        if (product > epsilon)
        {
          response.outcome = OptimizationOutcome::UNBOUNDED;
          response.rays.push_back(typename OptimizationResponse<Number>::Ray(
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

      double threshold;
      if (query.hasMinPrimalBound())
      {
        threshold = convertNumber<double>(query.minPrimalBound())
          + objectiveNormalization * _implementation->_normalizedPointEpsilon;
      }
      else
        threshold = -std::numeric_limits<double>::infinity();
      for (auto& pointProduct : _implementation->_pointProducts)
      {
        Number product = 0;
        for (const auto& iter : *_implementation->_points[pointProduct.vectorIndex].vector)
          product += objectiveVector[iter.first] * iter.second;

        if (product <= threshold)
          break;

        response.points.push_back(typename OptimizationResponse<Number>::Point(
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

  protected:
    PolyhedronImplementation<Number>* _implementation;
    std::size_t _queryCount;
  };

  template <typename NumberType>
  struct PolyhedronImplementation
  {
    typedef NumberType Number;
    typedef OptimizationOracle<Number> OptOracle;
    typedef SeparationOracle<Number> SepaOracle;
    typedef CacheOptimizationOracle<Number> CacheOptOracle;

    PolyhedronImplementation(std::shared_ptr<Space> space,
      const std::vector<std::shared_ptr<OptOracle>>& optOracles,
      const std::vector<std::shared_ptr<SepaOracle>>& sepaOracles)
      : _space(space)
    {
      if (space->dimension() == 0)
        throw std::runtime_error("Tried to initialize polyhedron in 0-dimensional space.");

      _cacheOptimization = std::make_shared<CacheOptOracle>(this); 
      _optimization.push_back(Data<OptOracle>(_cacheOptimization, true));
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

    ~PolyhedronImplementation()
    {
      delete[] _hashVector;
    }

    typename OptOracle::Response maximize(const Number* objectiveVector,
      const typename OptOracle::Query& query)
    {
      // TODO: collect returned dual bounds and consider minimum among them.

      sortOptimizationOracles();

      typename OptOracle::Response bestResponse;
      std::size_t lastOracle = 0;
      for (; lastOracle < _optimization.size(); ++lastOracle)
      {
        Data<OptOracle>& data = _optimization[lastOracle];

        // Run current oracle.

#if defined (IPO_DEBUG)
        std::cout << "\nCalling oracle <" << data.oracle->name() << ">." << std::endl;
#endif /* IPO_DEBUG */

        std::chrono::time_point<std::chrono::system_clock> timeStarted = std::chrono::system_clock::now();
        typename OptOracle::Response response = data.oracle->maximize(objectiveVector, query);
        double elapsed =  std::chrono::duration<double>(std::chrono::system_clock::now() - timeStarted).count();

        data.updateHistory(elapsed, response.wasSuccessful(), _historySize);

#if defined (IPO_DEBUG)
        std::cout << "    Oracle " << data.oracle->name() << " took "
          << (data.sumRunningTime / data.history.size()) << "s averaged over " << data.history.size()
          << " queries " << data.sumSuccess << " of which were successful.\n   ";
        if (query.hasMinPrimalBound())
          std::cout << " queried points better than " << query.minPrimalBound() << ";";
        std::cout << " response is " << response << ".\n" << std::flush;
#endif /* IPO_DEBUG */

        // If the current oracle found a ray or a point, we stop.
        if (!response.wasSuccessful())
          continue;

        // Add points/rays to cache.
        if (!data.isCache)
        {
          // We add points / rays to the cache.

          for (auto& ray : response.rays)
            addToCache(_rays, _hashToRayIndex, _rayProducts, ray.vector);

          for (auto& point : response.points)
            addToCache(_points, _hashToPointIndex, _pointProducts, point.vector);
        }

        if (response.outcome == OptimizationOutcome::INFEASIBLE
          || response.outcome == OptimizationOutcome::UNBOUNDED
          || (response.hasDualBound && response.primalBound == response.dualBound)
          || !query.hasMinPrimalBound() || response.primalBound >= query.minPrimalBound())
        {
          return response;
        }

        assert(response.outcome == OptimizationOutcome::FEASIBLE);

        bestResponse = response;

#if defined (IPO_DEBUG)
        std::cout << "    Did not find a solution of desired quality. Trying next oracle." << std::endl;
#endif /* IPO_DEBUG */
      }

#if defined (IPO_DEBUG)
        std::cout << "    Using last mentioned solution." << std::endl;
#endif /* IPO_DEBUG */

      return bestResponse;
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

      typename OptOracle::Response bestResponse;
      for (std::size_t o = 0; o < _optimization.size(); ++o)
      {
        Data<OptOracle>& data = _optimization[o];

        // Run current oracle.

#if defined (IPO_DEBUG)
        std::cout << "\nCalling oracle <" << data.oracle->name() << ">." << std::endl;
#endif /* IPO_DEBUG */

        std::chrono::time_point<std::chrono::system_clock> timeStarted =
          std::chrono::system_clock::now();
        typename OptOracle::Response response = data.oracle->maximizeDouble(objectiveVector, query);
        double elapsed = std::chrono::duration<double>(std::chrono::system_clock::now() - timeStarted).count();

        data.updateHistory(elapsed, response.wasSuccessful(), _historySize);
#if defined (IPO_DEBUG)
        std::cout << "Oracle " << data.oracle->name() << " took "
          << (data.sumRunningTime / data.history.size()) << "s averaged over " << data.history.size()
          << " queries " << data.sumSuccess << " of which were successful." << std::endl;
#endif /* IPO_DEBUG */

        // If the current oracle found a ray or a point, we stop.
        if (!response.wasSuccessful())
          break;

        // Add points/rays to cache.
        if (!data.isCache)
        {
          // We add points / rays to the cache.

          for (auto& ray : response.rays)
            addToCache(_rays, _hashToRayIndex, _rayProducts, ray.vector);

          for (auto& point : response.points)
            addToCache(_points, _hashToPointIndex, _pointProducts, point.vector);
        }

        if (response.outcome == OptimizationOutcome::INFEASIBLE
          || response.outcome == OptimizationOutcome::UNBOUNDED)
        {
          return response;
        }

        assert(response.outcome == OptimizationOutcome::FEASIBLE);

        bestResponse = response;
      }

      return bestResponse;
    }

    /**
     * \brief Adds given point to cache.
     *
     * Checks for duplicates.
     * 
     * \returns \c true if point was not cached before.
     */

    bool cachePoint(std::shared_ptr<sparse_vector<Number>> point)
    {
      return addToCache(_points, _hashToPointIndex, _pointProducts, point);
    }

    /**
     * \brief Adds given point to cache.
     *
     * Checks for duplicates.
     * 
     * \returns \c true if ray was not cached before.
     */

    bool cacheRay(std::shared_ptr<sparse_vector<Number>> ray)
    {
      return addToCache(_rays, _hashToRayIndex, _rayProducts, ray);
    }

    std::size_t numCachedSolutions() const
    {
      return _points.size() + _rays.size();
    }
    
    inline std::shared_ptr<Space> space() const
    {
      return _space;
    }

    std::pair<bool, std::shared_ptr<sparse_vector<Number>>> getCachedSolution(
      std::size_t index)
    {
      if (index < _rays.size())
        return std::pair<bool, std::shared_ptr<sparse_vector<Number>>>(false, _rays[index].vector);
      return std::pair<bool, std::shared_ptr<sparse_vector<Number>>>(true,
        _points[index - _rays.size()].vector);
    }
    
    typename SepaOracle::Response separate(const Number* vector, bool isPoint, const typename SepaOracle::Query& query)
    {
      typename SepaOracle::Response response;

      for (auto& sepa : _separation)
      {
        typename SepaOracle::Query currentQuery;
        currentQuery.maxNumInequalities = query.maxNumInequalities - response.constraints.size();
        auto currentResponse = sepa.oracle->separate(vector, isPoint, currentQuery);
        response.constraints.insert(response.constraints.end(),
          std::make_move_iterator(currentResponse.constraints.begin()),
          std::make_move_iterator(currentResponse.constraints.end()));
      }

      return response;
    }

    struct CachedVector
    {
      std::shared_ptr<sparse_vector<Number>> vector; /// The cached point/ray.
      double hash; /// The scalar product with the random vector used for hashing.
      std::size_t lastSuccess; /// Largest query number at which this vector was returned.
      double inverseNorm; /// Inverse of Euclidean norm of vector.

      CachedVector(std::shared_ptr<sparse_vector<Number>> vec, double hsh)
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

    /**
     * \brief Adds given vector to cache.
     * 
     * \returns \c true if the vector was not cached before.
     */

    bool addToCache(std::vector<CachedVector>& vectors, std::map<double, std::size_t>& hashToIndex,
      std::vector<ProductVector>& products, std::shared_ptr<sparse_vector<Number>>& vector)
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
        // TODO: Maybe Recompute objective value.
        return false;
      }
      else if (afterDist <= beforeDist && std::isfinite(afterDist) && squaredEuclideanDistance(
        *vector, *vectors[afterIter->second].vector) < _squaredDuplicateEpsilon)
      {
        vector = vectors[afterIter->second].vector;
        // TODO: Maybe Recompute objective value.
        return false;
      }
      else
      {
        hashToIndex.insert(afterIter, std::make_pair(hash, vectors.size()));
        products.push_back(ProductVector(vectors.size()));
        vectors.push_back(CachedVector(vector, hash));
        return true;
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
#if defined(IPO_DEBUG)
      std::cout << "We have to reduce the cache size. Currently " << (_points.size() + _rays.size()) << " of "
        << _maxCacheSize << std::endl;
#endif /* IPO_DEBUG */
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

  template<typename Number>
  Polyhedron<Number>::Polyhedron(std::shared_ptr<Space> space, const std::string& name)
    : OptimizationOracle<Number>(name), SeparationOracle<Number>(name)
  {
    Oracle<Number>::_space = space;
    std::vector<std::shared_ptr<OptimizationOracle<Number>>> optOracles;
    std::vector<std::shared_ptr<SeparationOracle<Number>>> sepaOracles;
    _implementation = new PolyhedronImplementation<Number>(Oracle<Number>::_space, optOracles, sepaOracles);
  }
  
  template <typename Number>
  Polyhedron<Number>::Polyhedron(std::shared_ptr<OptimizationOracle<Number>> optOracle)
    : OptimizationOracle<Number>("Polyhedron for " + optOracle->name()),
    SeparationOracle<Number>("Polyhedron for " + optOracle->name())
  {
    Oracle<Number>::_space = optOracle->space();
    std::vector<std::shared_ptr<OptimizationOracle<Number>>> optOracles;
    std::vector<std::shared_ptr<SeparationOracle<Number>>> sepaOracles;
    optOracles.push_back(optOracle);
    _implementation = new PolyhedronImplementation<Number>(Oracle<Number>::_space, optOracles, sepaOracles);
  }

  template <typename Number>
  Polyhedron<Number>::~Polyhedron()
  {
    delete static_cast<PolyhedronImplementation<Number>*>(_implementation);
  }

  template <typename Number>
  OptimizationResponse<Number> Polyhedron<Number>::maximize(
    const Number* objectiveVector, const OptimizationQuery<Number>& query)
  {
    return static_cast<PolyhedronImplementation<Number>*>(_implementation)->maximize(objectiveVector, query);
  }

  template <typename NumberType>
  SeparationResponse<NumberType> Polyhedron<NumberType>::separate(const NumberType* vector, bool isPoint,
      const SeparationQuery& query)
  {
    return static_cast<PolyhedronImplementation<NumberType>*>(_implementation)->separate(vector, isPoint,
        query);
  }

#if defined(IPO_DOUBLE)
  
  template class Polyhedron<double>;

#endif /* IPO_DOUBLE */

#if defined(IPO_RATIONAL)

  template class Polyhedron<rational>;

#endif /* IPO_RATIONAL */

} /* namespace ipo */
