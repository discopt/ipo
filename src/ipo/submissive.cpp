// #define IPO_DEBUG /* Uncomment to debug this file. */

#include <regex>

#include <ipo/submissive.hpp>

namespace ipo
{

  template <typename Number>
  SubmissiveOptimizationOracle<Number>::SubmissiveOptimizationOracle(
    std::shared_ptr<OptimizationOracle<Number>> sourceOracle, const std::string& name)
    : OptimizationOracle<Number>(name.empty() ? ("Submissive(" + sourceOracle->name() + ")") : name),
    _sourceOracle(sourceOracle)
  {
    this->_space = sourceOracle->space();
  }

  template <typename NumberType>
  SubmissiveOptimizationOracle<NumberType>::~SubmissiveOptimizationOracle()
  {

  }

  template <typename NumberType>
  OptimizationResponse<NumberType> SubmissiveOptimizationOracle<NumberType>::maximize(const NumberType* objectiveVector,
    const OptimizationQuery<NumberType>& query)
  {
    std::vector<bool> changed(this->space()->dimension());
    std::vector<NumberType> modifiedObjective(this->space()->dimension());
    for (size_t v = 0; v < this->space()->dimension(); ++v)
    {
      if ((changed[v] = (objectiveVector[v] < 0)))
        modifiedObjective[v] = 0;
      else
        modifiedObjective[v] = objectiveVector[v];
    }

    auto response = _sourceOracle->maximize(&modifiedObjective[0], query);
    for (auto& point : response.points)
    {
      bool change = false;
      for (auto& coef : *point.vector)
      {
        if (changed[coef.first])
        {
          change = true;
          break;
        }
      }

      if (change)
      {
        auto newVector = std::make_shared<sparse_vector<NumberType>>();
        for (const auto& coef : *point.vector)
        {
          if (!changed[coef.first])
            newVector->push_back(coef.first, coef.second);
        }
        point.vector = newVector;
      }
    }
    return response;
  }

  template class SubmissiveOptimizationOracle<double>;

#if defined(IPO_WITH_GMP)

  template class SubmissiveOptimizationOracle<mpq_class>;

#endif /* IPO_WITH_GMP */

} /* namespace ipo */
