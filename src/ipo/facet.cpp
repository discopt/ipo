#include <ipo/facet.hpp>

#include <ipo/lp.hpp>

namespace ipo
{
  FacetQuery::FacetQuery()
    : epsilonConstraints(1.0e-6), timeLimit(std::numeric_limits<double>::infinity())
  {

  }

  FacetQuery& FacetQuery::operator=(const FacetQuery& other)
  {
    epsilonConstraints = other.epsilonConstraints;
    timeLimit = other.timeLimit;
    return *this;
  }

  template <typename Number>
  FacetResult<Number>::FacetResult()
    : constraint(alwaysSatisfiedConstraint<Number>()), timeTotal(0.0)
  {

  }

  template <typename Number>
  FacetResult<Number>::FacetResult(FacetResult<Number>&& other)
    : constraint(std::move(other.constraint)), timeTotal(other.timeTotal), points(std::move(other.points)),
    rays(std::move(other.rays))
  {

  }

  template <typename Number>
  FacetResult<Number>& FacetResult<Number>::operator=(const FacetResult<Number>& other)
  {
    constraint = other.constraint;
    timeTotal = other.timeTotal;
    points = other.points;
    rays = other.rays;
    return *this;
  }

  template <typename Number>
  FacetResult<Number>& FacetResult<Number>::operator=(FacetResult<Number>&& other)
  {
    constraint = std::move(other.constraint);
    timeTotal = other.timeTotal;
    points = std::move(other.points);
    rays = std::move(other.rays);
    return *this;
  }

  template class FacetResult<double>;

#if defined(IPO_WITH_GMP) && defined(IPO_WITH_SOPLEX)

  template class FacetResult<mpq_class>;

#endif /* IPO_WITH_GMP && IPO_WITH_SOPLEX */

  template <>
  FacetResult<double> separatePoint(std::shared_ptr<Polyhedron<double>> polyhedron,
    const AffineHullResult<double>& affineHullResult, const double* point, const FacetQuery& query)
  {
    LP<double, true, false> lp;
    LPColumn x = lp.addColumn(0.0, 1.0, 0.0, "x");
    LPColumn y = lp.addColumn(0.0, 1.0, 0.0, "y");
    LPColumn z = lp.addColumn(0.0, 1.0, 0.0, "z");
    lp.update();
    LPColumn vars[3] = {x, y, z};
    double vals[3] = {1.0, 2.0, 3.0 };
    lp.addRow(-std::numeric_limits<double>::infinity(), 3, vars, vals, 2.0, "cons1");
    lp.solve();
  }

} /* namespace ipo */
