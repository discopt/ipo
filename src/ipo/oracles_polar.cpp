//#define IPO_DEBUG /* Uncomment to debug this file. */

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
      _lp.addColumn(_lp.minusInfinity(), _lp.plusInfinity(), 0, "beta");
      _lp.update();
      _normalizationRow = _lp.addRow(_lp.minusInfinity(), 0, NULL, NULL, _lp.plusInfinity());
    }

    void setAffineHull(const AffineHull<Number>& affineHull)
    {
      if (affineHull.dimension < 0)
        throw std::runtime_error("Trying to initialize PolarSeparationOracle for an infeasible polyhedron.");

#if defined(IPO_DEBUG)
      std::cout << "Adding " << affineHull.points.size() << " points and " << affineHull.rays.size()
        << " rays from affine hull." << std::endl;
#endif /* IPO_DEBUG */

      std::vector<int> nonzeroColumns;
      std::vector<Number> nonzeroCoefficients;
      for (const auto& point : affineHull.points)
      {
#if defined(IPO_DEBUG)
        std::cout << "Point " << _optOracle->space()->printVector(*point) << " contributes to relative interior."
          << std::endl;
#endif /* IPO_DEBUG */
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
#if defined(IPO_DEBUG)
        std::cout << "Ray " << _optOracle->space()->printVector(*ray) << " contributes to relative interior."
          << std::endl;
#endif /* IPO_DEBUG */
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
      {
#if defined(IPO_DEBUG)
        std::cout << "Scaling relative interior point for var#" << v << " to " << _interior[v] << " / "
          << affineHull.points.size() << " = ";
#endif /* IPO_DEBUG */
        _interior[v] /= affineHull.points.size();
#if defined(IPO_DEBUG)
        std::cout << _interior[v] << "." << std::endl;
#endif /* IPO_DEBUG */
      }
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
      _lp.changeObjective(n, isPoint ? -1 : 0);

      std::vector<int> nonzeroColumns;
      std::vector<Number> nonzeroCoefficients;
      for (std::size_t v = 0; v < n; ++v)
      {
        Number coefficient = vector[v] - (isPoint ? _interior[v] : 0);
#if defined(IPO_DEBUG)
        std::cout << "Coefficient of var#" << v << " is " << vector[v] << " - " << (isPoint ? _interior[v] : 0)
          << " = " << coefficient << std::endl;
#endif /* IPO_DEBUG */
        if (coefficient)
        {
          nonzeroColumns.push_back(v);
          nonzeroCoefficients.push_back(coefficient);
        }
      }
      _lp.changeRow(_normalizationRow, _lp.minusInfinity(), nonzeroColumns.size(), &nonzeroColumns[0],
        &nonzeroCoefficients[0], 1);

      while (true)
      {
#if defined(IPO_DEBUG)
        std::cout << "Solving Polar LP with " << _lp.numRows() << " constraints and " << _lp.numColumns()
          << " columns..." << std::flush;
        _lp.write("PolarSeparationOracle.lp");
#endif /* IPO_DEBUG */
        
        auto status = _lp.solve();

#if defined(IPO_DEBUG)
        std::cout << " done. Status: " << status << "." << std::endl;
#endif /* IPO_DEBUG */
        if (status == LPStatus::OPTIMAL)
        {
#if defined(IPO_DEBUG)
          std::cout << "The (potentially invalid) inequality has violation " << convertNumber<double>(_lp.getObjectiveValue()) << "." << std::endl;
#endif /* IPO_DEBUG */
          assert(_lp.hasPrimalSolution());
          std::vector<Number> solution = _lp.getPrimalSolution();

#if defined(IPO_DEBUG)
          bool first = true;
          for (std::size_t c = 0; c < n; ++c)
          {
            if (solution[c])
            {
              std::cout << (first ? "Coefficients: " : " + ") << solution[c] << "*" << _optOracle->space()->variable(c);
              first = false;
            }
          }
          std::cout << " <= " << convertNumber<double>(solution[n]) << " = " << solution[n] << std::endl;
#endif /* IPO_DEBUG */

          OptimizationQuery<Number> optQuery;
          optQuery.setMinPrimalBound(solution[n]);

#if defined(IPO_DEBUG)
          std::cout << "Calling optimization oracle with minimum primal bound " << solution[n] << "." << std::endl;
#endif /* IPO_DEBUG */
          OptimizationResponse<Number> optResponse = _optOracle->maximize(&solution[0], optQuery);

          auto violation = optResponse.primalBound() - solution[n];
          if (optResponse.outcome == OptimizationOutcome::UNBOUNDED ||
            (optResponse.outcome == OptimizationOutcome::FEASIBLE && violation > 1.0e-9))
          {
#if defined(IPO_DEBUG)
            std::cout << "Optimization oracle found a solution that violates the inequality by "
              << convertNumber<double>(violation) << " = "
              << (violation) << "." << std::endl;
#endif /* IPO_DEBUG */
            std::vector<int> nonzeroColumns;
            std::vector<Number> nonzeroCoefficients;
            for (const auto& point : optResponse.points)
            {
#if defined(IPO_DEBUG)
//               std::cout << _optOracle->space()->printVector(point.vector) << " with value "
//                 << convertNumber<double>(point.objectiveValue) << std::endl;
#endif /* IPO_DEBUG */
              nonzeroColumns.clear();
              nonzeroCoefficients.clear();
              for (const auto& iter : *point.vector)
              {
                nonzeroColumns.push_back(iter.first);
                nonzeroCoefficients.push_back(iter.second);
              }
              nonzeroColumns.push_back(_interior.size());
              nonzeroCoefficients.push_back(-1);
              _lp.addRow(_lp.minusInfinity(), nonzeroColumns.size(), &nonzeroColumns[0], &nonzeroCoefficients[0], 0);
            }
            for (const auto& ray : optResponse.rays)
            {
#if defined(IPO_DEBUG)
              std::cout << _optOracle->space()->printVector(ray.vector) << std::endl;
#endif /* IPO_DEBUG */
              nonzeroColumns.clear();
              nonzeroCoefficients.clear();
              for (const auto& iter : *ray.vector)
              {
                nonzeroColumns.push_back(iter.first);
                nonzeroCoefficients.push_back(iter.second);
              }
              _lp.addRow(_lp.minusInfinity(), nonzeroColumns.size(), &nonzeroColumns[0], &nonzeroCoefficients[0], 0);
            }
          }
          else
          {
#if defined(IPO_DEBUG)
            std::cout << "Primal bound is " << optResponse.primalBound() << " and dual bound is ";
            if (optResponse.hasDualBound)
              std::cout << optResponse.dualBound << "." << std::endl;
            else
              std::cout << "infinity." << std::endl;
            std::cout << "Optimization oracle proved validity of inequality." << std::endl;
#endif /* IPO_DEBUG */
            if (_lp.getObjectiveValue() > 0)
            {
              auto vector = std::make_shared<sparse_vector<Number>>();
              for (std::size_t v = 0; v < n; ++v)
              {
                if (solution[v])
                  vector->push_back(v, solution[v]);
              }
              response.constraints.emplace_back(Constraint<Number>(vector, solution[n]));
            }
            return response;
          }
        }
        else
        {
          assert(false);
        }
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
    std::shared_ptr<SeparationOracle<Number>> sepaRelaxationOracle, const std::string& name)
    : SeparationOracle<Number>(name.empty() ? ("Polar separation for " + name) : name)
  {
    _implementation = new PolarSeparationOracleImplementation<Number>(optOracle, sepaRelaxationOracle);
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

} /* namespace ipo */
