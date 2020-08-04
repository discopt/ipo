#include <ipo/oracles.hpp>

namespace ipo
{

  RealOptimizationOracle::RealOptimizationOracle(const std::string& name)
    : CommonOracle<double>(name)
  {

  }

  template <typename R>
  static std::ostream& printOptimizationRepsonse(std::ostream& stream,
    const CommonOptimizationReponse<R>& response)
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
        << " points, " << convertNumber<double>(response.primalBound) << " <= OPT";
      if (response.hasDualBound)
        stream << " <= " << convertNumber<double>(response.dualBound);
      return stream << "}";
    default:
      return stream << "{unknown error}";
    }
  }

  std::ostream& operator<<(std::ostream& stream, const CommonOptimizationReponse<double>& response)
  {
    return printOptimizationRepsonse(stream, response);
  }

  RealSeparationOracle::RealSeparationOracle(const std::string& name)
    : CommonOracle<double>(name)
  {

  }

  RealSeparationOracle::Response RealSeparationOracle::getInitial(
    const RealSeparationOracle::Query& query)
  {
    return Response();
  }

#if defined(IPO_WITH_GMP)

  RationalOptimizationOracle::RationalOptimizationOracle(const std::string& name)
    : CommonOracle<mpq_class>(name)
  {

  }

  RationalOptimizationOracle::Response RationalOptimizationOracle::maximize(
    const double* objectiveVector, const RationalOptimizationOracle::Query& query)
  {
    std::vector<mpq_class> rationalObjectiveVector(_space->dimension());
    for (std::size_t v = 0; v < _space->dimension(); ++v)
      rationalObjectiveVector[v] = objectiveVector[v];
    return maximize(&rationalObjectiveVector[0], query);
  }

  std::ostream& operator<<(std::ostream& stream,
    const CommonOptimizationReponse<mpq_class>& response)
  {
    return printOptimizationRepsonse(stream, response);
  }
  
  RationalSeparationOracle::RationalSeparationOracle(const std::string& name)
    : CommonOracle<mpq_class>(name)
  {

  }

  RationalSeparationOracle::Response RationalSeparationOracle::getInitial(
    const RationalSeparationOracle::Query& query)
  {
    return Response();
  }

  RationalSeparationOracle::Response RationalSeparationOracle::separate(const double* vector,
    bool isPoint, const Query& query)
  {
    std::vector<mpq_class> rationalVector(_space->dimension());
    for (std::size_t v = 0; v < _space->dimension(); ++v)
      rationalVector[v] = reconstructRational(vector[v]);
    return separate(&rationalVector[0], isPoint, query);
  }

#endif /* IPO_WITH_GMP */

}
