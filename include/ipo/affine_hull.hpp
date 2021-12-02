#pragma once

#include <ipo/config.hpp>
#include <ipo/export.hpp>
#include <ipo/polyhedron.hpp>
#include <ipo/constraint.hpp>

#include <vector>
#include <ostream>

namespace ipo
{
  struct AffineHullQuery
  {
    /**
     * \brief If a normalized vector cannot be combined from others up to this accuracy,
     *        then we are sure that it is linearly independent from these.
     */
    double epsilonLinearDependence;
    /**
     * \brief Decides when a constraint is satisfied with equality for double-precision
     * computations.
     */
    double epsilonConstraints;
    /**
     * \brief Decides whether a kernel vector should be invalidated.
     * 
     * Decides whether a kernel vector should be invalidated. Should be at most
     * \ref epsilonConstraints. When debugging numerically difficult instances, it can be set to
     * 0 to invalidate almost all vectors.
     */
    double epsilonInvalidate;
    /**
     * \brief If maximization yields a points with objective value within that bound, numerical
     * issues are likely and thus we try whether minimization gives a better result.
     */
    double epsilonSafety;
    /**
     * \brief Smallest value for the LU factorization.
     */
    double epsilonFactorization;
    /**
     * \brief Vector entries smaller than this are considered to be zero.
     */
    double epsilonCoefficient;
    double timeLimit; /// Time limit in seconds.

    IPO_EXPORT
    AffineHullQuery();

    IPO_EXPORT
    AffineHullQuery& operator=(const AffineHullQuery& other);
  };

  const static int AFFINEHULL_ERROR_RUNNING = -2;
  const static int AFFINEHULL_ERROR_TIMEOUT = -3;
  const static int AFFINEHULL_ERROR_NUMERICS = -4;


  template <typename R>
  struct CommonAffineHullResult
  {
    int dimension; /// Dimension or an error code.
    int lowerBound; /// Lower bound on dimension.
    int upperBound; /// Upper bound on dimension.
    double timeTotal; /// Total time spent.
    double timeOracles; /// Time spent for oracle calls.
    std::size_t numKernel; // Number of computed kernel vectors.
    double timeKernel; /// Time spent for computing kernel vectors.
    double timePointsRays; /// Time spent for the LU factorization of points and rays.
    double timeEquations; /// Time spent for the LU factorization of equations.
    std::vector<std::shared_ptr<sparse_vector<R>>> points; /// Points of inner description.
    std::vector<std::shared_ptr<sparse_vector<R>>> rays; /// Rays of inner description.
    std::vector<Constraint<R>> equations; /// Linearly independent equations.

    IPO_EXPORT
    CommonAffineHullResult()
    {
      dimension = AFFINEHULL_ERROR_RUNNING;
      lowerBound = -1;
      upperBound = std::numeric_limits<int>::max();
      timeTotal = 0.0;
      timeOracles = 0.0;
      numKernel = 0;
      timeKernel = 0.0;
      timePointsRays = 0.0;
      timeEquations = 0.0;
    }

    IPO_EXPORT
    CommonAffineHullResult(CommonAffineHullResult<R>&& other)
      : points(std::move(other.points)), rays(std::move(other.rays)),
      equations(std::move(other.equations))
    {
      dimension = other.dimension;
      lowerBound = other.lowerBound;
      upperBound = other.upperBound;
      timeTotal = other.timeTotal;
      timeOracles = other.timeOracles;
      numKernel = other.numKernel;
      timeKernel = other.timeKernel;
      timePointsRays = other.timePointsRays;
      timeEquations = other.timeEquations;
    }

    IPO_EXPORT
    CommonAffineHullResult<R>& operator=(const CommonAffineHullResult<R>& other)
    {
      points = other.points;
      rays = other.rays;
      equations = other.equations;
      dimension = other.dimension;
      lowerBound = other.lowerBound;
      upperBound = other.upperBound;
      timeTotal = other.timeTotal;
      timeOracles = other.timeOracles;
      numKernel = other.numKernel;
      timeKernel = other.timeKernel;
      timePointsRays = other.timePointsRays;
      timeEquations = other.timeEquations;
      return *this;
    }

    IPO_EXPORT
    CommonAffineHullResult<R>& operator=(CommonAffineHullResult<R>&& other)
    {
      points = std::move(other.points);
      rays = std::move(other.rays);
      equations = std::move(other.equations);
      dimension = other.dimension;
      lowerBound = other.lowerBound;
      upperBound = other.upperBound;
      timeTotal = other.timeTotal;
      timeOracles = other.timeOracles;
      numKernel = other.numKernel;
      timeKernel = other.timeKernel;
      timePointsRays = other.timePointsRays;
      timeEquations = other.timeEquations;
      return *this;
    }

    inline bool success() const
    {
      return dimension >= -1;
    }
  };

  template <typename R>
  std::ostream& operator<<(std::ostream& stream, const CommonAffineHullResult<R>& result)
  {
    if (result.dimension >= -1)
      stream << "Dimension: " << result.dimension;
    else
    {
      if (result.dimension == AFFINEHULL_ERROR_RUNNING)
        stream << "Running. ";
      else if (result.dimension == AFFINEHULL_ERROR_TIMEOUT)
        stream << "Timeout. ";
      else if (result.dimension == AFFINEHULL_ERROR_NUMERICS)
        stream << "Numerical difficulties. ";
      else
        stream << "Unknown error. ";
      stream << "Dimension in [" << result.lowerBound << "," << result.upperBound << "]";
    }
    stream << " (stats -- kernel vectors: " << result.numKernel << ")";
    stream << " (running times -- total: " << result.timeTotal << ", oracles: " << result.timeOracles
      << ", kernel: " << result.timeKernel << ", points/rays: " << result.timePointsRays
      << ", equations: " << result.timeEquations << ")";
    return stream;
  }

  typedef CommonAffineHullResult<double> RealAffineHullResult;

  IPO_EXPORT
  RealAffineHullResult affineHull(
    std::shared_ptr<RealPolyhedron> polyhedron,
    const AffineHullQuery& query = AffineHullQuery(),
    const std::vector<Constraint<double>>& knownEquations = std::vector<Constraint<double>>());

  IPO_EXPORT
  bool verifyAffineHullResult(std::shared_ptr<RealPolyhedron> polyhedron,
    const RealAffineHullResult& result);

#if defined(IPO_WITH_GMP)

  typedef CommonAffineHullResult<mpq_class> RationalAffineHullResult;

  IPO_EXPORT
  RationalAffineHullResult affineHull(
    std::shared_ptr<RationalPolyhedron> polyhedron,
    const AffineHullQuery& query = AffineHullQuery(),
    const std::vector<Constraint<mpq_class>>& knownEquations = std::vector<Constraint<mpq_class>>());

#endif /* IPO_WITH_GMP */

} /* namespace ipo */
