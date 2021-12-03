#include <ipo/oracles.hpp>

namespace ipo
{
  
  template <>
  OptimizationOracle<double>::OptimizationOracle(const std::string& name)
    : Oracle<double>(name)
  {

  }

  template <>
  OptimizationResponse<double> OptimizationOracle<double>::maximizeDouble(
    const double* objectiveVector, const OptimizationQuery<double>& query)
  {
    return maximize(objectiveVector, query);
  }

  template <typename R>
  static std::ostream& printOptimizationRepsonse(std::ostream& stream,
    const OptimizationResponse<R>& response)
  {
    switch(response.outcome)
    {
    case OptimizationOutcome::TIMEOUT:
      return stream << "{timeout error}";
    case OptimizationOutcome::FOUND_NOTHING:
      return stream << "{found nothing}";
    case OptimizationOutcome::UNBOUNDED:
      return stream << "{unbounded, " << response.rays.size() << " rays, " << response.points.size()
        << " points}";
    case OptimizationOutcome::INFEASIBLE:
      return stream << "{infeasible}";
    case OptimizationOutcome::FEASIBLE:
      assert(response.rays.empty());
      stream << "{feasible, " << response.rays.size() << " rays, " << response.points.size()
        << " points, opt in [" << convertNumber<double>(response.primalBound) << ",";
      if (response.hasDualBound)
        stream << convertNumber<double>(response.dualBound) << "]";
      else
        stream << "inf)";
      return stream << "}";
    default:
      return stream << "{unknown error}";
    }
  }

  std::ostream& operator<<(std::ostream& stream, const OptimizationResponse<double>& response)
  {
    return printOptimizationRepsonse(stream, response);
  }

  template <>
  SeparationOracle<double>::SeparationOracle(const std::string& name)
    : Oracle<double>(name)
  {

  }

  template <>
  SeparationResponse<double> SeparationOracle<double>::getInitial(
    const SeparationQuery& query)
  {
    return Response();
  }

  template <>
  SeparationResponse<double> SeparationOracle<double>::separateDouble(const double* vector,
    bool isPoint, const Query& query)
  {
    return separate(vector, isPoint, query);
  }

#if defined(IPO_WITH_GMP)

  template <>
  OptimizationOracle<mpq_class>::OptimizationOracle(const std::string& name)
    : Oracle<mpq_class>(name)
  {

  }

  template <>
  OptimizationResponse<mpq_class> OptimizationOracle<mpq_class>::maximizeDouble(
    const double* objectiveVector, const OptimizationQuery<mpq_class>& query)
  {
    std::vector<mpq_class> convertedObjectiveVector(this->_space->dimension());
    for (std::size_t v = 0; v < this->_space->dimension(); ++v)
      convertedObjectiveVector[v] = objectiveVector[v];
    return maximize(&convertedObjectiveVector[0], query);
  }

  std::ostream& operator<<(std::ostream& stream, const OptimizationResponse<mpq_class>& response)
  {
    return printOptimizationRepsonse(stream, response);
  }
  
  template<>
  SeparationOracle<mpq_class>::SeparationOracle(const std::string& name)
    : Oracle<mpq_class>(name)
  {

  }

  template <>
  SeparationResponse<mpq_class> SeparationOracle<mpq_class>::getInitial(
    const SeparationQuery& query)
  {
    return SeparationResponse<mpq_class>();
  }

  template <>
  SeparationResponse<mpq_class> SeparationOracle<mpq_class>::separateDouble(const double* vector,
    bool isPoint, const Query& query)
  {
    std::vector<mpq_class> rationalVector(_space->dimension());
    for (std::size_t v = 0; v < _space->dimension(); ++v)
      rationalVector[v] = reconstructRational(vector[v]);
    return separate(&rationalVector[0], isPoint, query);
  }

#endif /* IPO_WITH_GMP */

  template class OptimizationOracle<double>;
  template class SeparationOracle<double>;

#if defined(IPO_WITH_GMP)
  template class OptimizationOracle<mpq_class>;
  template class SeparationOracle<mpq_class>;
#endif /* IPO_WITH_GMP */

} /* namespace ipo */
