// #define IPO_DEBUG /* Uncomment to debug this file. */

#include <regex>

#include <ipo/projection.hpp>

namespace ipo
{
  template <typename NumberType>
  Projection<NumberType>::Projection()
    : _space(std::make_shared<Space>())
  {

  }

  template <typename NumberType>
  Projection<NumberType>::Projection(std::shared_ptr<Space> space, const std::vector<std::size_t>& variableIndices)
    : _space(std::make_shared<Space>())
  {
    addVariables(space, variableIndices);
  }

  template <typename NumberType>
  Projection<NumberType>::Projection(std::shared_ptr<Space> space, const std::string& regex)
    : _space(std::make_shared<Space>())
  {
    addVariables(space, regex);
  }

#if defined(IPO_WITH_GMP)

  template <>
  IPO_EXPORT
  Projection<mpq_class>::Projection(std::shared_ptr<Space> space, const std::string& regex)
    : _space(std::make_shared<Space>())
  {
    addVariables(space, regex);
  }
  
#endif /* IPO_WITH_GMP */

  template <typename NumberType>
  void Projection<NumberType>::addVariable(const std::string& name, const sparse_vector<Number>& linear,
    const Number& constant)
  {
    _space->addVariable(name);
    _mapLinear.push_back(sparse_vector<Number>(linear));
    _mapConstant.push_back(constant);
  }

  template <typename NumberType>
  Projection<NumberType>::~Projection()
  {

  }

#if defined(IPO_WITH_GMP)

  template <>
  IPO_EXPORT
  Projection<mpq_class>::~Projection()
  {

  }

#endif /* IPO_WITH_GMP */

  template <typename NumberType>
  void Projection<NumberType>::addVariable(std::shared_ptr<Space> space, std::size_t variableIndex)
  {
    std::vector<std::pair<std::size_t, NumberType>> linear(1, std::make_pair(variableIndex, NumberType(1)));
    _space->addVariable(space->variable(variableIndex));
    _mapLinear.push_back(sparse_vector<Number>(std::move(linear), false));
    _mapConstant.push_back(0);
  }

  template <typename NumberType>
  std::size_t Projection<NumberType>::addVariables(std::shared_ptr<Space> space,
    const std::vector<std::size_t>& variableIndices)
  {
    _mapLinear.reserve(_mapLinear.size() + variableIndices.size());
    for (auto variableIndex : variableIndices)
      addVariable(space, variableIndex);
    return variableIndices.size();
  }

  template <typename NumberType>
  std::size_t Projection<NumberType>::addVariables(std::shared_ptr<Space> space, const std::string& regex)
  {
    std::regex re(regex);
    std::size_t countMatches = 0;
    for (std::size_t v = 0; v < space->dimension(); ++v)
    {
      if (std::regex_match(space->variable(v), re))
      {
        addVariable(space, v);
        ++countMatches;
      }
    }
    return countMatches;
  }

  template <typename NumberType>
  std::shared_ptr<sparse_vector<NumberType>> Projection<NumberType>::projectPoint(
    std::shared_ptr<sparse_vector<NumberType>> point)
  {
    auto result = std::make_shared<sparse_vector<NumberType>>();
    for (size_t v = 0; v < space()->dimension(); ++v)
    {
      NumberType x = _mapConstant[v] + _mapLinear[v] * *point;
      if (x)
        result->push_back(v, x);
    }
    return result;
  }

  template <typename NumberType>
  std::shared_ptr<sparse_vector<NumberType>> Projection<NumberType>::projectRay(
    std::shared_ptr<sparse_vector<NumberType>> ray)
  {
    auto result = std::make_shared<sparse_vector<NumberType>>();
    for (size_t v = 0; v < space()->dimension(); ++v)
    {
      NumberType x = _mapLinear[v] * *ray;
      if (x)
        result->push_back(v, x);
    }
    return result;
  }

  template <typename NumberType>
  NumberType Projection<NumberType>::liftObjective(const NumberType* objective,
    std::vector<NumberType>& liftedObjective)
  {
    for (auto& coef : liftedObjective)
      coef = 0;
    NumberType result = 0;
    for (size_t v = 0; v < space()->dimension(); ++v)
    {
      if (objective[v])
      {
        for (auto coef : _mapLinear[v])
          liftedObjective[coef.first] += objective[v] * coef.second;
        result += objective[v] * _mapConstant[v];
      }
    }

    return result;
  }

  template class Projection<double>;
#if defined(IPO_WITH_GMP)
  template class Projection<mpq_class>;
#endif /* IPO_WITH_GMP */

  template <typename NumberType>
  ProjectionOptimizationOracle<NumberType>::ProjectionOptimizationOracle(
    std::shared_ptr<OptimizationOracle<NumberType>> sourceOracle,
    std::shared_ptr<Projection<NumberType>> projection,
    const std::string& name)
    : OptimizationOracle<NumberType>(name.empty() ? ("Projection(" + sourceOracle->name() + ")") : name),
    _projection(projection), _sourceOracle(sourceOracle)
  {
    this->_space = _projection->space();
  }

#if defined(IPO_WITH_GMP)

  template <>
  IPO_EXPORT
  ProjectionOptimizationOracle<mpq_class>::ProjectionOptimizationOracle(
    std::shared_ptr<OptimizationOracle<mpq_class>> sourceOracle,
    std::shared_ptr<Projection<mpq_class>> projection,
    const std::string& name)
    : OptimizationOracle<mpq_class>(name.empty() ? ("Projection(" + sourceOracle->name() + ")") : name),
    _projection(projection), _sourceOracle(sourceOracle)
  {
    this->_space = _projection->space();
  }

#endif /* IPO_WITH_GMP */

  template <typename NumberType>
  ProjectionOptimizationOracle<NumberType>::~ProjectionOptimizationOracle()
  {

  }

  template <typename NumberType>
  OptimizationResponse<NumberType> ProjectionOptimizationOracle<NumberType>::maximize(const NumberType* objectiveVector,
    const OptimizationQuery<NumberType>& query)
  {
    // Compute objective vector in space of source oracle.
    std::vector<NumberType> liftedObjective(_sourceOracle->space()->dimension());
    NumberType offset = _projection->liftObjective(objectiveVector, liftedObjective);

    OptimizationQuery<NumberType> liftedQuery;
    liftedQuery.maxNumSolutions = query.maxNumSolutions;
    liftedQuery.timeLimit = query.timeLimit;
    if (query.hasMinPrimalBound())
      liftedQuery.setMinPrimalBound(query.minPrimalBound() + offset);
    if (query.hasMaxDualBound())
      liftedQuery.setMaxDualBound(query.maxDualBound() + offset);

    auto liftedResponse = _sourceOracle->maximize(&liftedObjective[0], liftedQuery);

    OptimizationResponse<NumberType> response;
    response.outcome = liftedResponse.outcome;
    if (!liftedResponse.points.empty())
    {
      for (auto point : liftedResponse.points)
      {
        response.points.push_back(typename OptimizationResponse<NumberType>::Point(
          _projection->projectPoint(point.vector), point.objectiveValue + offset));
      }
    }
    else if (!liftedResponse.rays.empty())
    {
      for (auto ray : liftedResponse.rays)
      {
        response.rays.push_back(typename OptimizationResponse<NumberType>::Ray(
          _projection->projectRay(ray.vector)));
      }
    }

    return response;
  }

  
  template class ProjectionOptimizationOracle<double>;

#if defined(IPO_WITH_GMP)

  template class ProjectionOptimizationOracle<mpq_class>;

#endif /* IPO_WITH_GMP */

} /* namespace ipo */
