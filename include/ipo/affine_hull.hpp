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
     * \brief Decides when a constraint is satisfied with equality for double-precision
     * computations.
     */
    double epsilonConstraints;
    /**
     * \brief If maximization yields a points with objective value within that bound, numerical
     * issues are likely and thus we try whether minimization gives a better result.
     */
    double epsilonSafety;
    double timeLimit; /// Time limit in seconds.

    IPO_EXPORT
    AffineHullQuery();
  };

  const static int AFFINEHULL_ERROR_RUNNING = -2;
  const static int AFFINEHULL_ERROR_TIMEOUT = -3;

  struct AffineHullResultCommon
  {
    int dimension; /// Dimension or an error code.
    int lowerBound; /// Lower bound on dimension.
    int upperBound; /// Upper bound on dimension.
    double timeTotal; /// Total time spent.
    double timeOracles; /// Time spent for oracle calls.
    double timeKernel; /// Time spent for computing kernel vectors.
    double timePointsRays; /// Time spent for the LU factorization of points and rays.
    double timeEquations; /// Time spent for the LU factorization of equations.
  };

  template <typename T>
  struct AffineHullResult: public AffineHullResultCommon
  {
    std::vector<std::shared_ptr<sparse_vector<T>>> points; /// Points of inner description.
    std::vector<std::shared_ptr<sparse_vector<T>>> rays; /// Rays of inner description.
    std::vector<Constraint<T>> equations; /// Linearly independent equations.

    IPO_EXPORT
    AffineHullResult()
    {
      dimension = AFFINEHULL_ERROR_RUNNING;
      lowerBound = -1;
      upperBound = std::numeric_limits<int>::max();
      timeTotal = 0.0;
      timeOracles = 0.0;
      timeKernel = 0.0;
      timePointsRays = 0.0;
      timeEquations = 0.0;
    }

    IPO_EXPORT
    AffineHullResult(AffineHullResult<T>&& other)
      : points(std::move(other.points)), rays(std::move(other.rays)),
      equations(std::move(other.equations))
    {
      dimension = other.dimension;
      lowerBound = other.lowerBound;
      upperBound = other.upperBound;
      timeTotal = other.timeTotal;
      timeOracles = other.timeOracles;
      timeKernel = other.timeKernel;
      timePointsRays = other.timePointsRays;
      timeEquations = other.timeEquations;
    }
  };

  template <typename T>
  std::ostream& operator<<(std::ostream& stream, const AffineHullResult<T>& result)
  {
    if (result.dimension >= -1)
      stream << "Dimension: " << result.dimension;
    else
    {
      if (result.dimension == AFFINEHULL_ERROR_RUNNING)
        stream << "Running. ";
      else if (result.dimension == AFFINEHULL_ERROR_TIMEOUT)
        stream << "Timeout. ";
      else
        stream << "Unknown error. ";
      stream << "Dimension in [" << result.lowerBound << "," << result.upperBound << "]";
    }
    stream << " (running times -- total: " << result.timeTotal << ", oracles: " << result.timeOracles
      << ", kernel: " << result.timeKernel << ", points/rays: " << result.timePointsRays
      << ", equations: " << result.timeEquations << ")";
    return stream;
  }

  IPO_EXPORT
  AffineHullResult<double> affineHull(
    std::shared_ptr<Polyhedron<double>> polyhedron,
    const AffineHullQuery& query,
    const std::vector<Constraint<double>>& knownEquations = std::vector<Constraint<double>>());

#if defined(IPO_WITH_GMP)

  IPO_EXPORT
  AffineHullResult<mpq_class> affineHull(
    std::shared_ptr<Polyhedron<mpq_class>> polyhedron,
    const AffineHullQuery& query,
    const std::vector<Constraint<mpq_class>>& knownEquations = std::vector<Constraint<mpq_class>>());

#endif /* IPO_WITH_GMP */

} /* namespace ipo */
