#include <ipo/oracles_polar.hpp>

#include <ipo/lp.hpp>

namespace ipo
{

#if defined(IPO_DOUBLE_LP) || defined(IPO_RATIONAL_LP)
  
  template <typename NumberType>
  class PolarSeparationOracleImplementation
  {
  public:
    typedef NumberType Number;

    PolarSeparationOracleImplementation(std::shared_ptr<OptimizationOracle<Number>> optOracle,
      std::shared_ptr<SeparationOracle<Number>> sepaRelaxtionOracle)
      : _optOracle(optOracle), _sepaRelaxationOracle(sepaRelaxtionOracle), _interior(optOracle->space()->dimension(), 0)
    {
      _lp.setSense(ipo::LPSense::MAXIMIZE);
      for (size_t v = 0; v < _interior.size(); ++v)
        _lp.addColumn(_lp.minusInfinity(), _lp.plusInfinity(), 0, "alpha_" + optOracle->space()->variable(v));
      _lp.addColumn(_lp.minusInfinity(), _lp.plusInfinity(), -1, "beta");
      _lp.update();
      _normalizationRow = _lp.addRow(_lp.minusInfinity(), 0, NULL, NULL, _lp.plusInfinity());
    }

    void setAffineHull(const AffineHull<Number>& affineHull)
    {
      std::vector<int> nonzeroColumns;
      std::vector<Number> nonzeroCoefficients;

      for (const auto& point : affineHull.points)
      {
        std::cout << _optOracle->space()->printVector(*point) << std::endl;
        nonzeroColumns.clear();
        nonzeroCoefficients.clear();
        for (const auto& iter : *point)
        {
          nonzeroColumns.push_back(iter.first);
          nonzeroCoefficients.push_back(iter.second);
          _interior[iter.first] += iter.second;
        }
        nonzeroColumns.push_back(_interior.size());
        nonzeroCoefficients.push_back(-1);
        _lp.addRow(_lp.minusInfinity(), nonzeroColumns.size(), &nonzeroColumns[0], &nonzeroCoefficients[0], 0);
      }
      for (const auto& ray : affineHull.rays)
      {
        std::cout << _optOracle->space()->printVector(*ray) << std::endl;
        nonzeroColumns.clear();
        nonzeroCoefficients.clear();
        for (const auto& iter : *ray)
        {
          nonzeroColumns.push_back(iter.first);
          nonzeroCoefficients.push_back(iter.second);
          _interior[iter.first] += iter.second;
        }
        _lp.addRow(_lp.minusInfinity(), nonzeroColumns.size(), &nonzeroColumns[0], &nonzeroCoefficients[0], 0);
      }

      for (size_t v = 0; v < _interior.size(); ++v)
        _interior[v] /= affineHull.points.size();
    }

    /**
     * \brief Returns initially known inequalities.
     * 
     * The default implementation returns nothing.
     *
     * \param query Structure for query.
     * \return Separation result.
     **/

    SeparationResponse<Number> getInitial(const SeparationQuery& query)
    {
      return _sepaRelaxationOracle->getInitial(query);
    }

    /**
     * \brief Separates a point/ray.
     *
     * \param vector Array that maps coordinates to point/ray coordinates.
     * \param isPoint Whether a point shall be separated.
     * \param query Additional query information.
     *
     * \returns Reponse structure.
     **/

    SeparationResponse<Number> separate(const Number* vector, bool isPoint, const SeparationQuery& query)
    {
      SeparationResponse<Number> response;

      std::size_t n = _interior.size();
      for (std::size_t v = 0; v < n; ++v)
        _lp.changeObjective(v, vector[v]);

      std::vector<int> nonzeroColumns;
      std::vector<Number> nonzeroCoefficients;
      for (std::size_t v = 0; v < n; ++v)
      {
        Number coefficient = vector[v] - _interior[v];
        if (coefficient)
        {
          nonzeroColumns.push_back(v);
          nonzeroCoefficients.push_back(coefficient);
        }
      }
      _lp.changeRow(_normalizationRow, _lp.minusInfinity(), nonzeroColumns.size(), &nonzeroColumns[0],
        &nonzeroCoefficients[0], 1);

      _lp.write("separate.lp");

      while (true)
      {
        auto status = _lp.solve();

        std::cout << "LP status: " << status << "." << std::endl;
        if (status == LPStatus::OPTIMAL)
        {
          std::cout << "The constraint violation is " << convertNumber<double>(_lp.getObjectiveValue()) << "." << std::endl;
          assert(_lp.hasPrimalSolution());
          std::vector<Number> solution = _lp.getPrimalSolution();
          for (std::size_t c = 0; c < n; ++c)
          {
            std::cout << _optOracle->space()->variable(c) << " = " << solution[c] << std::endl;
          }
          std::cout << "beta = " << solution[n] << std::endl;
        }

        break;
      }

      return response;
    }

  protected:
    std::shared_ptr<OptimizationOracle<Number>> _optOracle;
    std::shared_ptr<SeparationOracle<Number>> _sepaRelaxationOracle;
    std::vector<Number> _interior;
    LP<Number> _lp;
    LPKey _normalizationRow;
  };

  template <typename Number>
  PolarSeparationOracle<Number>::PolarSeparationOracle(std::shared_ptr<OptimizationOracle<Number>> optOracle,
    std::shared_ptr<SeparationOracle<Number>> sepaRelaxtionOracle, const std::string& name)
    : SeparationOracle<Number>(name.empty() ? ("Polar separation for " + name) : name)
  {
    _implementation = new PolarSeparationOracleImplementation<Number>(optOracle, sepaRelaxtionOracle);
  }

  template <typename Number>
  void PolarSeparationOracle<Number>::setAffineHull(const AffineHull<Number>& affineHull)
  {
    static_cast<PolarSeparationOracleImplementation<Number>*>(_implementation)->setAffineHull(affineHull);
  }

  template <typename Number>
  SeparationResponse<Number> PolarSeparationOracle<Number>::getInitial(const SeparationQuery& query)
  {
    return static_cast<PolarSeparationOracleImplementation<Number>*>(_implementation)->getInitial(query);
  }

  template <typename Number>
  SeparationResponse<Number> PolarSeparationOracle<Number>::separate(const Number* vector, bool isPoint,
    const SeparationQuery& query)
  {
    return static_cast<PolarSeparationOracleImplementation<Number>*>(_implementation)->separate(vector, isPoint, query);
  }

#endif /* IPO_DOUBLE_LP || IPO_RATIONAL_LP */

#if defined(IPO_DOUBLE_LP)

  template class PolarSeparationOracleImplementation<double>;
  template class PolarSeparationOracle<double>;

#endif /* IPO_DOUBLE_LP */

#if defined(IPO_RATIONAL_LP)

  template class PolarSeparationOracleImplementation<rational>;
  template class PolarSeparationOracle<rational>;

#endif /* IPO_RATIONAL_LP */
  
//   SeparationQuery::SeparationQuery()
//     : epsilonConstraints(1.0e-6), timeLimit(std::numeric_limits<double>::infinity())
//   {
// 
//   }
// 
//   SeparationQuery& SeparationQuery::operator=(const SeparationQuery& other)
//   {
//     epsilonConstraints = other.epsilonConstraints;
//     timeLimit = other.timeLimit;
//     return *this;
//   }

//   template <typename Number>
//   Separated<Number>::Separated()
//     : constraint(alwaysSatisfiedConstraint<Number>()), timeTotal(0.0)
//   {
// 
//   }
// 
//   template <typename Number>
//   Separated<Number>::Separated(Separated<Number>&& other)
//     : constraint(std::move(other.constraint)), timeTotal(other.timeTotal), points(std::move(other.points)),
//     rays(std::move(other.rays))
//   {
// 
//   }
// 
//   template <typename Number>
//   Separated<Number>& Separated<Number>::operator=(const Separated<Number>& other)
//   {
//     constraint = other.constraint;
//     timeTotal = other.timeTotal;
//     points = other.points;
//     rays = other.rays;
//     return *this;
//   }
// 
//   template <typename Number>
//   Separated<Number>& Separated<Number>::operator=(Separated<Number>&& other)
//   {
//     constraint = std::move(other.constraint);
//     timeTotal = other.timeTotal;
//     points = std::move(other.points);
//     rays = std::move(other.rays);
//     return *this;
//   }
// 
//   template class Separated<double>;
// 
// #if defined(IPO_WITH_GMP) && defined(IPO_WITH_SOPLEX)
// 
//   template class Separated<mpq_class>;
// 
// #endif /* IPO_WITH_GMP && IPO_WITH_SOPLEX */
// 
//   template <typename Number>
//   Separated<Number> separatePoint(std::shared_ptr<Polyhedron<Number>> polyhedron,
//     const AffineHullResult<Number>& affineHull, const Number* point, const SeparationQuery& query)
//   {
//     Separated<Number> result;
// 
//     LP<Number> lp;
//     lp.setSense(ipo::LPSense::MAXIMIZE);
//     size_t n = polyhedron->space()->dimension();
//     for (size_t v = 0; v < n; ++v)
//       lp.addColumn(lp.minusInfinity(), lp.plusInfinity(), point[v], "alpha_" + polyhedron->space()->variable(v));
//     ipo::LPColumn beta = lp.addColumn(lp.minusInfinity(), lp.plusInfinity(), -1, "beta");
// 
//     std::vector<Number> interiorPoint(n,0);
//     std::vector<ipo::LPColumn> nonzeroColumns;
//     std::vector<Number> nonzeroCoefficients;
// 
//     for (const auto& point : affineHull.points)
//     {
//       std::cout << polyhedron->space()->printVector(*point) << std::endl;
//       nonzeroColumns.clear();
//       nonzeroCoefficients.clear();
//       for (const auto& iter : *point)
//       {
//         nonzeroColumns.push_back(iter.first);
//         nonzeroCoefficients.push_back(iter.second);
//         interiorPoint[iter.first] += iter.second;
//       }
//       nonzeroColumns.push_back(beta);
//       nonzeroCoefficients.push_back(-1);
//       lp.addRow(lp.minusInfinity(), nonzeroColumns.size(), &nonzeroColumns[0], &nonzeroCoefficients[0], 0);
//     }
//     for (const auto& ray : affineHull.rays)
//     {
//       std::cout << polyhedron->space()->printVector(*ray) << std::endl;
//       nonzeroColumns.clear();
//       nonzeroCoefficients.clear();
//       for (const auto& iter : *ray)
//       {
//         nonzeroColumns.push_back(iter.first);
//         nonzeroCoefficients.push_back(iter.second);
//         interiorPoint[iter.first] += iter.second;
//       }
//       lp.addRow(lp.minusInfinity(), nonzeroColumns.size(), &nonzeroColumns[0], &nonzeroCoefficients[0], 0);
//     }
// 
//     for (size_t v = 0; v < n; ++v)
//     {
//       interiorPoint[v] = point[v] - interiorPoint[v] / affineHull.points.size();
//     }
// 
//     nonzeroColumns.clear();
//     nonzeroCoefficients.clear();
//     for (size_t v = 0; v < n; ++v)
//     {
//       if (interiorPoint[v])
//       {
//         nonzeroColumns.push_back(v);
//         nonzeroCoefficients.push_back(interiorPoint[v]);
//       }
//     }
//     lp.addRow(lp.minusInfinity(), nonzeroColumns.size(), &nonzeroColumns[0], &nonzeroCoefficients[0], 1);
//  
//     lp.write("separation.lp");
//     auto status = lp.solve();
// 
//     std::cout << "LP status: " << status << "." << std::endl;
//     if (status == ipo::LPStatus::OPTIMAL)
//     {
//       std::cout << "The constraint violation is " << ipo::convertNumber<double>(lp.getObjectiveValue()) << "." << std::endl;
//       assert(lp.hasPrimalSolution());
//       std::vector<Number> solution = lp.getPrimalSolution();
//       for (std::size_t c = 0; c < n; ++c)
//       {
//         std::cout << polyhedron->space()->variable(c) << " = " << solution[c] << std::endl;
//       }
//       std::cout << "beta = " << solution[n] << std::endl;
//     }
// 
//     return result;
//   }
// 
//   template Separated<double> separatePoint(std::shared_ptr<Polyhedron<double>> polyhedron,
//     const AffineHullResult<double>& affineHullResult, const double* point, const SeparationQuery& query);
// 
// #if defined(IPO_WITH_GMP) && defined(IPO_WITH_SOPLEX)
// 
//   template Separated<mpq_class> separatePoint(std::shared_ptr<Polyhedron<mpq_class>> polyhedron,
//     const AffineHullResult<mpq_class>& affineHullResult, const mpq_class* point, const SeparationQuery& query);
// 
// #endif /* IPO_WITH_GMP && IPO_WITH_SOPLEX */

} /* namespace ipo */
