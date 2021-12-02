// #define IPO_DEBUG // Uncomment to debug this file.

#include <regex>

#include <ipo/oracles_projection.hpp>

namespace ipo
{

  ProjectionRealOptimizationOracle::ProjectionRealOptimizationOracle(std::shared_ptr<RealOptimizationOracle> sourceOracle,
    const std::string& name)
    : RealOptimizationOracle(name.empty() ? ("projection of " + sourceOracle->name()) : sourceOracle->name()),
    _sourceOracle(sourceOracle)
  {

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
    space()->addVariable((*_sourceOracle->space())[sourceVariableIndex]);
  }

  void ProjectionRealOptimizationOracle::addVariables(const std::string& regex)
  {
    std::regex re(regex);
    Space& sourceSpace = *_sourceOracle->space();
    for (std::size_t v = 0; v < sourceSpace.dimension(); ++v)
    {
      if (std::regex_match(sourceSpace[v], re))
        addVariable(v);
    }
  }

  RealOptimizationOracle::Response ProjectionRealOptimizationOracle::maximize(const double* objectiveVector,
    const Query& query)
  {
    size_t projDimension = space()->dimension();
    size_t sourceDimension = _sourceOracle->space()->dimension();

    // Compute objective vector in space of source oracle.
    auto liftedObjective = new double[sourceDimension];
    for (size_t v = 0; v < sourceDimension; ++v)
      liftedObjective[v] = 0.0;
    for (size_t v = 0; v < projDimension; ++v)
    {
      if (objectiveVector[v] != 0.0)
      {
        for (auto coef : _projectionLinear[v])
          liftedObjective[coef.first] += coef.second * objectiveVector[v];
      }
    }

    RealOptimizationOracle::Query liftedQuery;
    liftedQuery.maxNumSolutions = query.maxNumSolutions;
    liftedQuery.timeLimit = query.timeLimit;
    auto liftedResponse = _sourceOracle->maximize(liftedObjective, liftedQuery);
    delete[] liftedObjective;

    RealOptimizationOracle::Response response;
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
        response.points.push_back(RealOptimizationOracle::Response::Point(vector));
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
        response.rays.push_back(RealOptimizationOracle::Response::Ray(vector));
      }
    }

    return response;
  }

#if defined(IPO_WITH_GMP)

  ProjectionRationalOptimizationOracle::ProjectionRationalOptimizationOracle(
    std::shared_ptr<RationalOptimizationOracle> sourceOracle, const std::string& name)
    : RationalOptimizationOracle(name.empty() ? ("projection of " + sourceOracle->name()) : sourceOracle->name()),
    _sourceOracle(sourceOracle)
  {

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
    space()->addVariable((*_sourceOracle->space())[sourceVariableIndex]);
  }

  void ProjectionRationalOptimizationOracle::addVariables(const std::string& regex)
  {
    std::regex re(regex);
    Space& sourceSpace = *_sourceOracle->space();
    for (std::size_t v = 0; v < sourceSpace.dimension(); ++v)
    {
      if (std::regex_match(sourceSpace[v], re))
        addVariable(v);
    }
  }

  RationalOptimizationOracle::Response ProjectionRationalOptimizationOracle::maximize(const mpq_class* objectiveVector,
    const Query& query)
  {
    size_t projDimension = space()->dimension();
    size_t sourceDimension = _sourceOracle->space()->dimension();

    // Compute objective vector in space of source oracle.
    auto liftedObjective = new mpq_class[sourceDimension];
    for (size_t v = 0; v < sourceDimension; ++v)
      liftedObjective[v] = 0.0;
    for (size_t v = 0; v < projDimension; ++v)
    {
      if (objectiveVector[v] != 0.0)
      {
        for (auto coef : _projectionLinear[v])
          liftedObjective[coef.first] += coef.second * objectiveVector[v];
      }
    }

    RationalOptimizationOracle::Query liftedQuery;
    liftedQuery.maxNumSolutions = query.maxNumSolutions;
    liftedQuery.timeLimit = query.timeLimit;
    auto liftedResponse = _sourceOracle->maximize(liftedObjective, liftedQuery);
    delete[] liftedObjective;

    RationalOptimizationOracle::Response response;
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
        response.points.push_back(RationalOptimizationOracle::Response::Point(vector));
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
        response.rays.push_back(RationalOptimizationOracle::Response::Ray(vector));
      }
    }

    return response;
  }

#endif /* IPO_WITH_GMP */

} /* namespace ipo */
