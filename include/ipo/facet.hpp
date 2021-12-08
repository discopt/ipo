#pragma once

#include <ipo/config.hpp>
#include <ipo/export.hpp>
#include <ipo/polyhedron.hpp>
#include <ipo/constraint.hpp>
#include <ipo/affine_hull.hpp>

#include <vector>
#include <ostream>

namespace ipo
{
  struct FacetQuery
  {
    /**
     * \brief Decides when a constraint is satisfied with equality for double-precision
     * computations.
     */
    double epsilonConstraints;
    double timeLimit; /// Time limit in seconds.

    FacetQuery();

    FacetQuery& operator=(const FacetQuery& other);
  };

  const static int FACET_ERROR_RUNNING = -2;
  const static int FACET_ERROR_TIMEOUT = -3;
  const static int FACET_ERROR_NUMERICS = -4;

  template <typename NumberType>
  struct FacetResult
  {
    typedef NumberType Number;

    Constraint<Number> constraint; /// The computed facet or equation.
    double timeTotal; /// Total time spent.
    std::vector<std::shared_ptr<sparse_vector<Number>>> points; /// Points of inner description.
    std::vector<std::shared_ptr<sparse_vector<Number>>> rays; /// Rays of inner description.

    FacetResult();

    FacetResult(FacetResult<Number>&& other);

    FacetResult<Number>& operator=(const FacetResult<Number>& other);

    FacetResult<Number>& operator=(FacetResult<Number>&& other);
  };

  template <typename Number>
  FacetResult<Number> separatePoint(std::shared_ptr<Polyhedron<Number>> polyhedron,
    const AffineHullResult<Number>& affineHullResult, const Number* point, const FacetQuery& query = FacetQuery());

  template <typename Number>
  FacetResult<Number> separateRay(std::shared_ptr<Polyhedron<Number>> polyhedron,
    const AffineHullResult<Number>& affineHullResult, const Number* ray, const FacetQuery& query = FacetQuery());

} /* namespace ipo */
