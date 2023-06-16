//#define IPO_DEBUG /* Uncomment to debug this file. */

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
#if defined(IPO_DEBUG)
      std::cout << "DominantOptimizationOracle::maximize() called." << std::endl;
#endif // IPO_DEBUG

    for (size_t v = 0; v < this->space()->dimension(); ++v)
    {
      if (objectiveVector[v] > 0)
      {
        OptimizationResponse<NumberType> response;
        std::vector<NumberType> zeroVector(this->space()->dimension(), 0);
        std::vector<std::pair<std::size_t, NumberType>> vectorData(1, std::make_pair(v, (NumberType)(1)));
        response.rays.push_back(typename OptimizationResponse<NumberType>::Ray(
          std::make_shared<sparse_vector<NumberType>>(vectorData, false)));
        response.outcome = OptimizationOutcome::UNBOUNDED;
        const OptimizationQuery<NumberType> sourceQuery;
        auto subResponse = _sourceOracle->maximize(&zeroVector[0], sourceQuery);

#if defined(IPO_DEBUG)
        std::cout << "Response of sub-oracle of DominantOptimizationOracle: " << subResponse << std::endl;
#endif // IPO_DEBUG

        assert(subResponse.outcome == OptimizationOutcome::FEASIBLE);
        response.points = std::move(subResponse.points);
        for (auto& point : response.points)
        {
          assert(point.objectiveValue == 0);
          point.objectiveValue = *point.vector * objectiveVector;
        }
        return response;
      }
    }

    return _sourceOracle->maximize(objectiveVector, query);
  }
  
#if defined(IPO_DOUBLE)

  template class DominantOptimizationOracle<double>;

#endif /* IPO_DOUBLE */

#if defined(IPO_RATIONAL)

  template class DominantOptimizationOracle<rational>;

#endif /* IPO_RATIONAL */

} /* namespace ipo */
