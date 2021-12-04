#pragma once

#include <ipo/config.hpp>
#include <ipo/export.hpp>
#include <ipo/oracles.hpp>

namespace ipo
{

  template <typename NumberType>
  class ForestOptimizationOracle: public OptimizationOracle<NumberType>
  {
  public:
    IPO_EXPORT
    ForestOptimizationOracle(std::size_t numNodes, std::pair<std::size_t, std::size_t>* edgesFirst,
      std::pair<std::size_t, std::size_t>* edgesBeyond, bool spanning, const std::string& name = "forest");

    IPO_EXPORT
    ~ForestOptimizationOracle();

    /**
     * \brief Maximize an objective vector of type double.
     *
     * \param objectiveVector Array that maps coordinates to objective value coefficients.
     * \param query Parameters of query.
     * \return Optimization result.
     **/

    IPO_EXPORT
    virtual OptimizationResponse<NumberType> maximize(const NumberType* objectiveVector,
      const OptimizationQuery<NumberType>& query);

  protected:
    bool _spanning;
    std::size_t _numNodes;
    std::vector<std::pair<std::size_t, std::size_t>> _edges;
  };

} /* namespace ipo */
