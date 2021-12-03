#pragma once

#include <ipo/config.hpp>
#include <ipo/export.hpp>
#include <ipo/oracles.hpp>

namespace ipo
{

  class ForestRealOptimizationOracle: public OptimizationOracle<double>
  {
  public:
    IPO_EXPORT
    ForestRealOptimizationOracle(std::size_t numNodes, std::pair<std::size_t, std::size_t>* edgesFirst,
      std::pair<std::size_t, std::size_t>* edgesBeyond, bool spanning, const std::string& name = "forest");

    IPO_EXPORT
    ~ForestRealOptimizationOracle();

    /**
     * \brief Maximize an objective vector of type double.
     *
     * \param objectiveVector Array that maps coordinates to objective value coefficients.
     * \param query Parameters of query.
     * \return Optimization result.
     **/

    IPO_EXPORT
    virtual OptimizationOracle<double>::Response maximize(const double* objectiveVector,
      const OptimizationOracle<double>::Query& query);

  protected:
    bool _spanning;
    std::size_t _numNodes;
    std::vector<std::pair<std::size_t, std::size_t>> _edges;
  };
  
#if defined(IPO_WITH_GMP)

  class ForestRationalOptimizationOracle: public OptimizationOracle<mpq_class>
  {
  public:
    IPO_EXPORT
    ForestRationalOptimizationOracle(std::size_t numNodes, std::pair<std::size_t, std::size_t>* edgesFirst,
      std::pair<std::size_t, std::size_t>* edgesBeyond, bool spanning, const std::string& name = "forest");

    IPO_EXPORT
    ~ForestRationalOptimizationOracle();

    /**
     * \brief Maximize an objective vector of type mpq_class.
     *
     * \param objectiveVector Array that maps coordinates to objective value coefficients.
     * \param query Parameters of query.
     * \return Optimization result.
     **/

    IPO_EXPORT
    virtual OptimizationOracle<mpq_class>::Response maximize(const mpq_class* objectiveVector,
      const OptimizationOracle<mpq_class>::Query& query);

  protected:
    bool _spanning;
    std::size_t _numNodes;
    std::vector<std::pair<std::size_t, std::size_t>> _edges;
  };

#endif /* IPO_WITH_GMP */

} /* namespace ipo */
