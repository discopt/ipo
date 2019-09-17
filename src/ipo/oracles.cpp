#include <ipo/oracles.hpp>

#include <algorithm>
#include <cassert>

#include <iostream>

namespace ipo
{

  Oracle::Oracle(const std::string& name)
    : _name(name), _space(nullptr)
  {

  }

  OptimizationOracle::Query::Query()
  {
    reset();
  }

  void OptimizationOracle::Query::reset()
  {
#if defined(IPO_WITH_GMP)
    rational = false;
#endif
    minNumSolutions = 1;
    maxNumSolutions = 10;
    minObjectiveValue = std::numeric_limits<double>::infinity();
    timeLimit = std::numeric_limits<double>::infinity();
  }

  OptimizationOracle::Result::Point::Point(const Vector& vec)
    : objectiveValue(-std::numeric_limits<double>::signaling_NaN()), vector(vec)
  {

  }

  OptimizationOracle::Result::Point::Point(const Vector& vec, const Value& value)
    : objectiveValue(value), vector(vec)
  {
    
  }

  OptimizationOracle::Result::Ray::Ray(const Vector& vec)
    : vector(vec)
  {

  }

  OptimizationOracle::Result::Result()
  {
    reset();
  }

  void OptimizationOracle::Result::reset()
  {
    hitTimeLimit = false;
    points.clear();
    rays.clear();
  }

  void OptimizationOracle::Result::sortPoints()
  {
    std::sort(points.begin(), points.end());
  }
  
  bool OptimizationOracle::Result::checkPointsSorted() const
  {
    return std::is_sorted(points.begin(), points.end());
  }

  OptimizationOracle::OptimizationOracle(const std::string& name)
    : Oracle(name)
  {

  }

#if defined(IPO_WITH_GMP)

  void OptimizationOracle::maximize(const mpq_class* objectiveVector,
    const OptimizationOracle::Query& query, OptimizationOracle::Result& result)
  {
    // Create floating-point approximation of objective vector.
    double* approximateObjectiveVector = new double[space()->dimension()];
    for (std::size_t i = 0; i < space()->dimension(); ++i)
      approximateObjectiveVector[i] = objectiveVector[i].get_d();

    this->maximize(approximateObjectiveVector, query, result);

    delete[] approximateObjectiveVector;
  }

#endif /* IPO_WITH_GMP */

  SeparationOracle::Query::Query()
  {
    reset();
  }

  void SeparationOracle::Query::reset()
  {
#if defined(IPO_WITH_GMP)
    rational = false;
#endif /* IPO_WITH_GMP */
    maxNumInequalities = 50;
    timeLimit = std::numeric_limits<double>::infinity();
  }

  SeparationOracle::Result::Result()
  {
    reset();
  }

  void SeparationOracle::Result::reset()
  {
    hitTimeLimit = false;
    constraints.clear();
  }

  SeparationOracle::SeparationOracle(const std::string& name)
    : Oracle(name)
  {

  }

  void SeparationOracle::getInitial(const SeparationOracle::Query& query,
    SeparationOracle::Result& result)
  {
    result.reset();
  }

#if defined(IPO_WITH_GMP)

  bool SeparationOracle::separate(const mpq_class* vector, bool isPoint,
    const SeparationOracle::Query& query, SeparationOracle::Result& result)
  {
    // Create floating-point approximation of point.
    double* approximateVector = new double[space()->dimension()];
    for (std::size_t i = 0; i < space()->dimension(); ++i)
      approximateVector[i] = vector[i].get_d();

    bool separated = this->separate(approximateVector, isPoint, query, result);

    delete[] approximateVector;

    return separated;
  }

#endif /* IPO_WITH_GMP */

} /* namespace ipo */
