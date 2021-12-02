// #define IPO_DEBUG /* Uncomment to debug this file. */

#include <regex>

#include <ipo/oracles_projection.hpp>

namespace ipo
{
  template <typename NumberType>
  Projection<NumberType>::Projection()
    : _space(std::make_shared<Space>())
  {

  }

  template <typename NumberType>
  Projection<NumberType>::~Projection()
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

  template <typename NumberType>
  void Projection<NumberType>::addVariable(const std::string& name, const sparse_vector<Number>& linear,
    const Number& constant)
  {
    _space->addVariable(name);
    _mapLinear.push_back(sparse_vector<Number>(linear));
    _mapConstant.push_back(constant);
  }

  template <typename NumberType>
  void Projection<NumberType>::addVariable(std::shared_ptr<Space> space, std::size_t variableIndex)
  {
    std::vector<std::size_t, NumberType> linear(1, std::make_pair(variableIndex, NumberType(1)));
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


  typedef Projection<double> Projection_double;

#if defined(IPO_WITH_GMP)
  typedef Projection<mpq_class> Projection_mpq_class;
#endif /* IPO_WITH_GMP */
  

  ProjectionRealOptimizationOracle::ProjectionRealOptimizationOracle(std::shared_ptr<OptimizationOracle<double>> sourceOracle,
    const std::string& name)
    : OptimizationOracle<double>(name.empty() ? ("projection of " + sourceOracle->name()) : sourceOracle->name()),
    _sourceOracle(sourceOracle)
  {
    _space = std::make_shared<Space>();
  }

  ProjectionRealOptimizationOracle::~ProjectionRealOptimizationOracle()
  {

  }

  void ProjectionRealOptimizationOracle::addVariable(const size_t sourceVariableIndex)
  {
    sparse_vector<double> vector;
    vector.push_back(sourceVariableIndex, 1.0);
    _projectionLinear.push_back(vector);
    _projectionConstant.push_back(0.0);
    space()->addVariable(_sourceOracle->space()->variable(sourceVariableIndex));
  }

  std::size_t ProjectionRealOptimizationOracle::addVariables(const std::string& regex)
  {
    std::regex re(regex);
    Space& sourceSpace = *_sourceOracle->space();
    std::size_t countMatches = 0;
    for (std::size_t v = 0; v < sourceSpace.dimension(); ++v)
    {
#if defined(IPO_DEBUG)
      std::cout << "Checking whether '" << sourceSpace[v] << "' matches '" << regex << "'." << std::endl;
#endif /* IPO_DEBUG */
      if (std::regex_match(sourceSpace.variable(v), re))
      {
        addVariable(v);
        ++countMatches;
      }
    }
    return countMatches;
  }

  OptimizationOracle<double>::Response ProjectionRealOptimizationOracle::maximize(const double* objectiveVector,
    const Query& query)
  {
    size_t projDimension = space()->dimension();
    size_t sourceDimension = _sourceOracle->space()->dimension();

    // Compute objective vector in space of source oracle.
    auto liftedObjective = new double[sourceDimension];
    for (size_t v = 0; v < sourceDimension; ++v)
      liftedObjective[v] = 0.0;
    double offset = 0; // Product of objective with constant part of projection map.
    for (size_t v = 0; v < projDimension; ++v)
    {
      if (objectiveVector[v] != 0.0)
      {
        for (auto coef : _projectionLinear[v])
          liftedObjective[coef.first] += coef.second * objectiveVector[v];
        offset += objectiveVector[v] * _projectionConstant[v];
      }
    }

    OptimizationOracle<double>::Query liftedQuery;
    liftedQuery.maxNumSolutions = query.maxNumSolutions;
    liftedQuery.timeLimit = query.timeLimit;
    if (query.hasMinPrimalBound())
      liftedQuery.setMinPrimalBound(query.minPrimalBound() + offset);
    if (query.hasMaxDualBound())
      liftedQuery.setMaxDualBound(query.maxDualBound() + offset);

    auto liftedResponse = _sourceOracle->maximize(liftedObjective, liftedQuery);
    delete[] liftedObjective;

    OptimizationOracle<double>::Response response;
    response.outcome = liftedResponse.outcome;
    if (!liftedResponse.points.empty())
    {
      for (auto point : liftedResponse.points)
      {
        auto vector = std::make_shared<sparse_vector<double>>();
        for (size_t v = 0; v < projDimension; ++v)
        {
          double x = _projectionConstant[v] + _projectionLinear[v] * *point.vector;
          if (x != 0.0)
            vector->push_back(v, x);
        }
        response.points.push_back(OptimizationOracle<double>::Response::Point(vector, point.objectiveValue + offset));
      }
    }
    else if (!liftedResponse.rays.empty())
    {
      for (auto ray : liftedResponse.rays)
      {
        auto vector = std::make_shared<sparse_vector<double>>();
        for (size_t v = 0; v < projDimension; ++v)
        {
          double x = _projectionLinear[v] * *ray.vector;
          if (x != 0.0)
            vector->push_back(v, x);
        }
        response.rays.push_back(OptimizationOracle<double>::Response::Ray(vector));
      }
    }

    return response;
  }

#if defined(IPO_WITH_GMP)

  ProjectionRationalOptimizationOracle::ProjectionRationalOptimizationOracle(
    std::shared_ptr<OptimizationOracle<mpq_class>> sourceOracle, const std::string& name)
    : OptimizationOracle<mpq_class>(name.empty() ? ("projection of " + sourceOracle->name()) : sourceOracle->name()),
    _sourceOracle(sourceOracle)
  {
    _space = std::make_shared<Space>();
  }

  ProjectionRationalOptimizationOracle::~ProjectionRationalOptimizationOracle()
  {
    
  }

  void ProjectionRationalOptimizationOracle::addVariable(const size_t sourceVariableIndex)
  {
    sparse_vector<mpq_class> vector;
    vector.push_back(sourceVariableIndex, 1.0);
    _projectionLinear.push_back(vector);
    _projectionConstant.push_back(0.0);
    space()->addVariable(_sourceOracle->space()->variable(sourceVariableIndex));
  }

  std::size_t ProjectionRationalOptimizationOracle::addVariables(const std::string& regex)
  {
    std::regex re(regex);
    Space& sourceSpace = *_sourceOracle->space();
    std::size_t countMatches = 0;
    for (std::size_t v = 0; v < sourceSpace.dimension(); ++v)
    {
      if (std::regex_match(sourceSpace.variable(v), re))
      {
        addVariable(v);
        ++countMatches;
      }
    }
    return countMatches;
  }

  OptimizationOracle<mpq_class>::Response ProjectionRationalOptimizationOracle::maximize(const mpq_class* objectiveVector,
    const Query& query)
  {
    size_t projDimension = space()->dimension();
    size_t sourceDimension = _sourceOracle->space()->dimension();

#if defined(IPO_DEBUG)
    std::cout << "RationalOptimizationOracle::maximize() called with objective vector (";
    for (size_t v = 0; v < projDimension; ++v)
      std::cout << (v == 0 ? "" : ",") << objectiveVector[v];
    std::cout << ")" << std::endl;
#endif /* IPO_DEBUG */
    
    // Compute objective vector in space of source oracle.
    auto liftedObjective = new mpq_class[sourceDimension];
    for (size_t v = 0; v < sourceDimension; ++v)
      liftedObjective[v] = 0.0;
    mpq_class offset = 0; // Product of objective with constant part of projection map.
    for (size_t v = 0; v < projDimension; ++v)
    {
      if (objectiveVector[v] != 0.0)
      {
        for (auto coef : _projectionLinear[v])
          liftedObjective[coef.first] += objectiveVector[v] * coef.second;
        offset += objectiveVector[v] * _projectionConstant[v];
      }
    }

    OptimizationOracle<mpq_class>::Query liftedQuery;
    liftedQuery.maxNumSolutions = query.maxNumSolutions;
    liftedQuery.timeLimit = query.timeLimit;
    if (query.hasMinPrimalBound())
      liftedQuery.setMinPrimalBound(query.minPrimalBound() + offset);
    if (query.hasMaxDualBound())
      liftedQuery.setMaxDualBound(query.maxDualBound() + offset);

    auto liftedResponse = _sourceOracle->maximize(liftedObjective, liftedQuery);
    delete[] liftedObjective;

    OptimizationOracle<mpq_class>::Response response;
    response.outcome = liftedResponse.outcome;
    if (!liftedResponse.points.empty())
    {
      for (auto point : liftedResponse.points)
      {
        auto vector = std::make_shared<sparse_vector<mpq_class>>();
        for (size_t v = 0; v < projDimension; ++v)
        {
          mpq_class x = _projectionConstant[v] + _projectionLinear[v] * *point.vector;
          if (x != 0.0)
            vector->push_back(v, x);
        }
        response.points.push_back(OptimizationOracle<mpq_class>::Response::Point(vector, point.objectiveValue + offset));
      }
    }
    else if (!liftedResponse.rays.empty())
    {
      for (auto ray : liftedResponse.rays)
      {
        auto vector = std::make_shared<sparse_vector<mpq_class>>();
        for (size_t v = 0; v < projDimension; ++v)
        {
          mpq_class x = _projectionLinear[v] * *ray.vector;
          if (x != 0.0)
            vector->push_back(v, x);
        }
        response.rays.push_back(OptimizationOracle<mpq_class>::Response::Ray(vector));
      }
    }

    return response;
  }

#endif /* IPO_WITH_GMP */

} /* namespace ipo */
