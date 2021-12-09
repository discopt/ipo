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

  template <typename Number>
  FacetResult<Number> separatePoint(std::shared_ptr<Polyhedron<Number>> polyhedron,
    const AffineHullResult<Number>& affineHull, const Number* point, const FacetQuery& query)
  {
    FacetResult<Number> result;

    LP<Number> lp;
    lp.setSense(ipo::LPSense::MAXIMIZE);
    size_t n = polyhedron->space()->dimension();
    for (size_t v = 0; v < n; ++v)
      lp.addColumn(lp.minusInfinity(), lp.plusInfinity(), point[v], "alpha_" + polyhedron->space()->variable(v));
    ipo::LPColumn beta = lp.addColumn(lp.minusInfinity(), lp.plusInfinity(), -1, "beta");

    std::vector<Number> interiorPoint(n,0);
    std::vector<ipo::LPColumn> nonzeroColumns;
    std::vector<Number> nonzeroCoefficients;

    for (const auto& point : affineHull.points)
    {
      std::cout << polyhedron->space()->printVector(*point) << std::endl;
      nonzeroColumns.clear();
      nonzeroCoefficients.clear();
      for (const auto& iter : *point)
      {
        nonzeroColumns.push_back(iter.first);
        nonzeroCoefficients.push_back(iter.second);
        interiorPoint[iter.first] += iter.second;
      }
      nonzeroColumns.push_back(beta);
      nonzeroCoefficients.push_back(-1);
      lp.addRow(lp.minusInfinity(), nonzeroColumns.size(), &nonzeroColumns[0], &nonzeroCoefficients[0], 0);
    }
    for (const auto& ray : affineHull.rays)
    {
      std::cout << polyhedron->space()->printVector(*ray) << std::endl;
      nonzeroColumns.clear();
      nonzeroCoefficients.clear();
      for (const auto& iter : *ray)
      {
        nonzeroColumns.push_back(iter.first);
        nonzeroCoefficients.push_back(iter.second);
        interiorPoint[iter.first] += iter.second;
      }
      lp.addRow(lp.minusInfinity(), nonzeroColumns.size(), &nonzeroColumns[0], &nonzeroCoefficients[0], 0);
    }

    for (size_t v = 0; v < n; ++v)
    {
      interiorPoint[v] = point[v] - interiorPoint[v] / affineHull.points.size();
    }

    nonzeroColumns.clear();
    nonzeroCoefficients.clear();
    for (size_t v = 0; v < n; ++v)
    {
      if (interiorPoint[v])
      {
        nonzeroColumns.push_back(v);
        nonzeroCoefficients.push_back(interiorPoint[v]);
      }
    }
    lp.addRow(lp.minusInfinity(), nonzeroColumns.size(), &nonzeroColumns[0], &nonzeroCoefficients[0], 1);
 
    lp.write("separation.lp");
    auto status = lp.solve();

    std::cout << "LP status: " << status << "." << std::endl;
    if (status == ipo::LPStatus::OPTIMAL)
    {
      std::cout << "The constraint violation is " << ipo::convertNumber<double>(lp.getObjectiveValue()) << "." << std::endl;
      assert(lp.hasPrimalSolution());
      std::vector<Number> solution = lp.getPrimalSolution();
      for (std::size_t c = 0; c < n; ++c)
      {
        std::cout << polyhedron->space()->variable(c) << " = " << solution[c] << std::endl;
      }
      std::cout << "beta = " << solution[n] << std::endl;
    }

    return result;
  }

  template FacetResult<double> separatePoint(std::shared_ptr<Polyhedron<double>> polyhedron,
    const AffineHullResult<double>& affineHullResult, const double* point, const FacetQuery& query);

#if defined(IPO_WITH_GMP) && defined(IPO_WITH_SOPLEX)

  template FacetResult<mpq_class> separatePoint(std::shared_ptr<Polyhedron<mpq_class>> polyhedron,
    const AffineHullResult<mpq_class>& affineHullResult, const mpq_class* point, const FacetQuery& query);

#endif /* IPO_WITH_GMP && IPO_WITH_SOPLEX */

} /* namespace ipo */
