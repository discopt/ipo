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

    AffineHullQuery();

    AffineHullQuery& operator=(const AffineHullQuery& other);
  };

  const static int AFFINEHULL_ERROR_RUNNING = -2;
  const static int AFFINEHULL_ERROR_TIMEOUT = -3;
  const static int AFFINEHULL_ERROR_NUMERICS = -4;

  template <typename NumberType>
  struct AffineHullResult
  {
    typedef NumberType Number;

    int dimension; /// Dimension or an error code.
    int lowerBound; /// Lower bound on dimension.
    int upperBound; /// Upper bound on dimension.
    double timeTotal; /// Total time spent.
    double timeOracles; /// Time spent for oracle calls.
    std::size_t numKernel; // Number of computed kernel vectors.
    double timeKernel; /// Time spent for computing kernel vectors.
    double timePointsRays; /// Time spent for the LU factorization of points and rays.
    double timeEquations; /// Time spent for the LU factorization of equations.
    std::vector<std::shared_ptr<sparse_vector<NumberType>>> points; /// Points of inner description.
    std::vector<std::shared_ptr<sparse_vector<NumberType>>> rays; /// Rays of inner description.
    std::vector<Constraint<NumberType>> equations; /// Linearly independent equations.

    AffineHullResult();

    AffineHullResult(AffineHullResult<Number>&& other);

    AffineHullResult<Number>& operator=(const AffineHullResult<Number>& other);

    AffineHullResult<Number>& operator=(AffineHullResult<Number>&& other);

    inline bool success() const
    {
      return dimension >= -1;
    }
  };

  template <typename Number>
  std::ostream& operator<<(std::ostream& stream, const AffineHullResult<Number>& result);

  template <typename Number>
  AffineHullResult<Number> affineHull(
    std::shared_ptr<Polyhedron<Number>> polyhedron, const AffineHullQuery& query = AffineHullQuery(),
    const std::vector<Constraint<Number>>& knownEquations = std::vector<Constraint<Number>>());

  bool verifyAffineHullResult(std::shared_ptr<Polyhedron<double>> polyhedron, const AffineHullResult<double>& result);

} /* namespace ipo */
