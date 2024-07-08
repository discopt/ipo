// #define IPO_DEBUG /* Uncomment to debug this file. */

#include <regex>
#include <unordered_map>

#include <ipo/projection.hpp>

namespace ipo
{
  template <typename Number>
  Projection<Number>::Projection()
    : _space(std::make_shared<Space>())
  {

  }

  template <typename Number>
  Projection<Number>::Projection(std::shared_ptr<Space> space, const std::vector<std::size_t>& variableIndices)
    : _space(std::make_shared<Space>())
  {
    addVariables(space, variableIndices);
  }

  template <typename Number>
  Projection<Number>::Projection(std::shared_ptr<Space> space, const std::string& regex)
    : _space(std::make_shared<Space>())
  {
    addVariables(space, regex);
  }

  template <typename Number>
  void Projection<Number>::addVariable(const std::string& name, const sparse_vector<Number>& linear,
    const Number& constant)
  {
    _space->addVariable(name);
    _mapLinear.push_back(sparse_vector<Number>(linear));
    _mapConstant.push_back(constant);
  }

  template <typename Number>
  Projection<Number>::~Projection()
  {

  }

  template <typename Number>
  void Projection<Number>::addVariable(std::shared_ptr<Space> space, std::size_t variableIndex)
  {
    std::vector<std::pair<std::size_t, Number>> linear(1, std::make_pair(variableIndex, Number(1)));
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

  template <typename Number>
  std::shared_ptr<sparse_vector<Number>> Projection<Number>::projectPoint(std::shared_ptr<sparse_vector<Number>> point)
  {
    auto result = std::make_shared<sparse_vector<Number>>();
    for (size_t v = 0; v < space()->dimension(); ++v)
    {
      Number x = _mapConstant[v] + _mapLinear[v] * *point;
      if (x)
        result->push_back(v, x);
    }
    return result;
  }

  template <typename Number>
  std::shared_ptr<sparse_vector<Number>> Projection<Number>::projectRay(std::shared_ptr<sparse_vector<Number>> ray)
  {
    auto result = std::make_shared<sparse_vector<Number>>();
    for (size_t v = 0; v < space()->dimension(); ++v)
    {
      Number x = _mapLinear[v] * *ray;
      if (x)
        result->push_back(v, x);
    }
    return result;
  }

  template <typename Number>
  Number Projection<Number>::liftObjective(const Number* objective, std::vector<Number>& liftedObjective)
  {
    for (auto& coef : liftedObjective)
      coef = 0;
    Number result = 0;
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

#if defined(IPO_DOUBLE)

  template class Projection<double>;

#endif /* IPO_DOUBLE */

#if defined(IPO_RATIONAL)

  template class Projection<rational>;

#endif /* IPO_RATIONAL */

  template <typename Number>
  ProjectionOptimizationOracle<Number>::ProjectionOptimizationOracle(
    std::shared_ptr<OptimizationOracle<Number>> sourceOracle,
    std::shared_ptr<Projection<Number>> projection,
    const std::string& name)
    : OptimizationOracle<Number>(name.empty() ? ("Projection(" + sourceOracle->name() + ")") : name),
    _projection(projection), _sourceOracle(sourceOracle)
  {
    this->_space = _projection->space();
  }

  template <typename Number>
  std::vector<Constraint<Number>> projectionEquationsImplementation(std::shared_ptr<Projection<Number>> projection,
    const std::vector<Constraint<Number>>& equations)
  {
    std::vector<Constraint<Number>> result;

#if defined(IPO_DEBUG)
    // TODO: Implement.
#endif /* IPO_DEBUG */

    return result;
  }

  template <typename Number>
  std::vector<Constraint<Number>> projectionEquations(std::shared_ptr<Projection<Number>> projection,
    const std::vector<Constraint<Number>>& equations)
  {
    return projectionEquationsImplementation(projection, equations);
  }

  /* Explicit template instantiation. */

  template std::vector<Constraint<double>> projectionEquations(std::shared_ptr<Projection<double>> projection,
    const std::vector<Constraint<double>>& equations);

#if defined(IPO_RATIONAL)

  template std::vector<Constraint<rational>> projectionEquations(std::shared_ptr<Projection<rational>> projection,
    const std::vector<Constraint<rational>>& equations);

#endif /* IPO_RATIONAL */



  template <typename Number>
  std::vector<Constraint<Number>> projectionCompactImplementation(std::shared_ptr<Projection<Number>> projection,
    const std::vector<Constraint<Number>>& constraints, bool onlyEquations)
  {
    std::vector<Constraint<Number>> result;

    // TODO: Implement completely. A temporary workaround is to determine the set of source variables that are
    //       projected orthogonally.

    // Mapping from projectable variables to their image variables.
    std::unordered_map<std::size_t, std::size_t> orthogonallyProjectedVariables;
    for (size_t v = 0; v < projection->space()->dimension(); ++v)
    {
      const sparse_vector<Number>& linear = projection->mapLinear(v);
      if (linear.size() == 1)
        orthogonallyProjectedVariables.insert(std::make_pair(linear.begin()->first, v));
    }

    for (const Constraint<Number>& constraint : constraints)
    {
      if ( onlyEquations && constraint.type() != ConstraintType::EQUATION)
        continue;

      bool projectable = true;
      std::vector<std::pair<std::size_t, Number>> nonzeros;
      Number lhs = constraint.hasLhs() ? constraint.lhs() : Number(0);
      Number rhs = constraint.hasRhs() ? constraint.rhs() : Number(0);
      for (auto& term : constraint.vector())
      {
        auto iter = orthogonallyProjectedVariables.find(term.first);
        if (iter == orthogonallyProjectedVariables.end())
        {
          projectable = false;
          break;
        }
        else
        {
          const sparse_vector<Number>& linear = projection->mapLinear(iter->second);
          const Number& factor = linear.begin()->second;
          const Number& constant = projection->mapConstant(iter->second);
          nonzeros.push_back(std::make_pair(iter->second, term.second / factor));
          if (constraint.hasLhs())
            lhs -= term.second * constant / factor;
          if (constraint.hasRhs())
            rhs -= term.second * constant / factor;
        }
      }
      if (projectable)
      {
        result.emplace_back(Constraint<Number>(std::move(lhs),
          std::make_shared<sparse_vector<Number>>(nonzeros, true), std::move(rhs), constraint.type()));
      }
    }

#if defined(IPO_DEBUG)
    // TODO: Implement.
#endif /* IPO_DEBUG */

    return result;
  }

  template <typename Number>
  std::vector<Constraint<Number>> projectionCompact(std::shared_ptr<Projection<Number>> projection,
    const std::vector<Constraint<Number>>& constraints, bool onlyEquations)
  {
    return projectionCompactImplementation(projection, constraints, onlyEquations);
  }

  /* Explicit template instantiation. */

  template std::vector<Constraint<double>> projectionCompact(std::shared_ptr<Projection<double>> projection,
    const std::vector<Constraint<double>>& constraints, bool onlyEquations);

#if defined(IPO_RATIONAL)

  template std::vector<Constraint<rational>> projectionCompact(std::shared_ptr<Projection<rational>> projection,
    const std::vector<Constraint<rational>>& constraints, bool onlyEquations);

#endif /* IPO_RATIONAL */





  template <typename Number>
  ProjectionOptimizationOracle<Number>::~ProjectionOptimizationOracle()
  {

  }

  template <typename Number>
  OptimizationResponse<Number> ProjectionOptimizationOracle<Number>::maximize(const Number* objectiveVector,
    const OptimizationQuery<Number>& query)
  {
#if defined(IPO_DEBUG)
      std::cout << "ProjectionOptimizationOracle::maximize() called." << std::endl;
#endif // IPO_DEBUG

    // Compute objective vector in space of source oracle.
    std::vector<Number> liftedObjective(_sourceOracle->space()->dimension());
    Number offset = _projection->liftObjective(objectiveVector, liftedObjective);

    OptimizationQuery<Number> liftedQuery;
    liftedQuery.maxNumSolutions = query.maxNumSolutions;
    liftedQuery.timeLimit = query.timeLimit;
    if (query.hasMinPrimalBound())
      liftedQuery.setMinPrimalBound(query.minPrimalBound() + offset);
    if (query.hasMaxDualBound())
      liftedQuery.setMaxDualBound(query.maxDualBound() + offset);

    auto liftedResponse = _sourceOracle->maximize(&liftedObjective[0], liftedQuery);
#if defined(IPO_DEBUG)
    std::cout << "Lifted response = " << liftedResponse << std::endl;
#endif /* IPO_DEBUG */

    OptimizationResponse<Number> response;
    response.outcome = liftedResponse.outcome;
    if (liftedResponse.hasDualBound)
    {
      response.hasDualBound = true;
      response.dualBound = liftedResponse.dualBound - offset;
    }
    if (!liftedResponse.points.empty())
    {
      for (auto point : liftedResponse.points)
      {
        response.points.push_back(typename OptimizationResponse<Number>::Point(
          _projection->projectPoint(point.vector), point.objectiveValue + offset));
      }
      response.setPrimalBound(response.points.front().objectiveValue);
    }
    if (!liftedResponse.rays.empty())
    {
      for (auto ray : liftedResponse.rays)
      {
        response.rays.push_back(typename OptimizationResponse<Number>::Ray(
          _projection->projectRay(ray.vector)));
      }
    }

#if defined(IPO_DEBUG)
    std::cout << "Projected response = " << response << std::endl;
#endif /* IPO_DEBUG */

    return response;
  }

  template class ProjectionOptimizationOracle<double>;

#if defined(IPO_RATIONAL)

  template class ProjectionOptimizationOracle<rational>;

#endif /* IPO_RATIONAL */

} /* namespace ipo */
