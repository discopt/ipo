// #define IPO_DEBUG /* Uncomment to debug this file. */

#include <regex>

#include <ipo/dominant.hpp>

namespace ipo
{

  template <typename NumberType>
  DominantOptimizationOracle<NumberType>::DominantOptimizationOracle(
    std::shared_ptr<OptimizationOracle<NumberType>> sourceOracle,
    const std::string& name)
    : OptimizationOracle<NumberType>(name.empty() ? ("Dominant(" + sourceOracle->name() + ")") : name),
    _sourceOracle(sourceOracle)
  {
    this->_space = sourceOracle->space();
  }

  template <typename NumberType>
  DominantOptimizationOracle<NumberType>::~DominantOptimizationOracle()
  {

  }

  template <typename NumberType>
  OptimizationResponse<NumberType> DominantOptimizationOracle<NumberType>::maximize(const NumberType* objectiveVector,
    const OptimizationQuery<NumberType>& query)
  {
    OptimizationResponse<NumberType> response;
    for (size_t v = 0; v < this->space()->dimension(); ++v)
    {
      if (objectiveVector[v] > 0)
      {
        std::vector<std::pair<std::size_t, NumberType>> vectorData(1, std::make_pair(v, (NumberType)(1)));
        response.rays.push_back(typename OptimizationResponse<NumberType>::Ray(
          std::make_shared<sparse_vector<NumberType>>(vectorData, false)));
        response.outcome = OptimizationOutcome::UNBOUNDED;
        return response;
      }
    }

    return _sourceOracle->maximize(objectiveVector, query);
  }
  
  template class DominantOptimizationOracle<double>;
#if defined(IPO_WITH_GMP)
  template class DominantOptimizationOracle<mpq_class>;
#endif /* IPO_WITH_GMP */

} /* namespace ipo */
