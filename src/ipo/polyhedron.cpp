#include <ipo/polyhedron.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>

namespace ipo
{
  class CacheOptimizationOracle : public OptimizationOracle
  {
  public:
    CacheOptimizationOracle()
      : OptimizationOracle("CacheOptimizationOracle"), _normalizedRayEpsilon(1.0e-6),
      _normalizedPointEpsilon(1.0e-6), _queryCount(0)
    {
      
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
      for (std::size_t i = 0; i < maxNum; ++i)
      {
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
      for (std::size_t i = 0; i < maxNum; ++i)
      {
        result.points.push_back(OptimizationOracle::Result::Point(_points[i].vector,
          _points[i].product));
        _points[i].lastSuccess = _queryCount;
      }
    }

#endif /* IPO_WITH_GMP */

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
  Polyhedron::Data<OptimizationOracle>::Data(std::shared_ptr<OptimizationOracle> o)
    : expectedRunningTime(0), expectedSuccess(0), oracle(o)
  {

  }

  template <>
  Polyhedron::Data<SeparationOracle>::Data(std::shared_ptr<SeparationOracle> o)
    : expectedRunningTime(0), expectedSuccess(0), oracle(o)
  {

  }

  Polyhedron::Polyhedron(std::shared_ptr<OptimizationOracle> optimizationOracle)
  {
    _optimization.push_back(Data<OptimizationOracle>(optimizationOracle));
  }

  Polyhedron::Polyhedron(std::shared_ptr<SeparationOracle> o)
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

} /* namespace ipo */
