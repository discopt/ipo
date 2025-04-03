#include <ipo/oracles.hpp>

namespace ipo
{

  template <typename Number>
  Oracle<Number>::Oracle(const std::string& name)
    : _name(name), _space(nullptr)
  {

  }

  template <typename Number>
  Oracle<Number>::~Oracle()
  {

  }

  template <typename R>
  static std::ostream& printOptimizationQuery(std::ostream& stream,
    const OptimizationQuery<R>& query)
  {
    stream << "{";
    bool isFirst = true;
    if (query.hasMinPrimalBound())
    {
      if (!isFirst)
        stream << ", ";
      stream << "only sol.values > " << formatNumberApprox(query.minPrimalBound());
      isFirst = false;
    }
    if (query.hasMaxDualBound())
    {
      if (!isFirst)
        stream << ", ";
      stream << "stops if opt <= " << formatNumberApprox(query.maxDualBound());
      isFirst = false;
    }
    if (query.timeLimit < std::numeric_limits<double>::infinity())
    {
      if (!isFirst)
        stream << ", ";
      stream << "timout at " << query.timeLimit << "s";
      isFirst = false;
    }
    return stream << "}";
  }

  std::ostream& operator<<(std::ostream& stream, const OptimizationQuery<double>& query)
  {
    return printOptimizationQuery(stream, query);
  }

#if defined(IPO_WITH_RATIONAL)

  std::ostream& operator<<(std::ostream& stream, const OptimizationQuery<rational>& query)
  {
    return printOptimizationQuery(stream, query);
  }

#endif /* IPO_WITH_RATIONAL */
  
  template <typename Number>
  OptimizationOracle<Number>::OptimizationOracle()
    : Oracle<Number>("abstract OptimizationOracle")
  {

  }

  template <typename Number>
  OptimizationOracle<Number>::~OptimizationOracle()
  {

  }

  template <>
  OptimizationResponse<double> OptimizationOracle<double>::maximizeDouble(
    const double* objectiveVector, const OptimizationQuery<double>& query)
  {
    return maximize(objectiveVector, query);
  }

  template <typename Number>
  TrustRegionOptimizationOracle<Number>::TrustRegionOptimizationOracle()
    : Oracle<Number>("abstract TrustRegionOptimizationOracle")
  {

  }

  template <typename Number>
  TrustRegionOptimizationOracle<Number>::~TrustRegionOptimizationOracle()
  {

  }

  template <>
  OptimizationResponse<double> TrustRegionOptimizationOracle<double>::maximizeTrustRegionDouble(
    const ManhattanTrustRegion<double>& trustRegion, const double* objectiveVector, const OptimizationQuery<double>& query)
  {
    return maximizeTrustRegion(trustRegion, objectiveVector, query);
  }
  
#if defined(IPO_WITH_RATIONAL)

  template <>
  OptimizationResponse<rational> OptimizationOracle<rational>::maximizeDouble(
    const double* objectiveVector, const OptimizationQuery<rational>& query)
  {
    std::vector<rational> convertedObjectiveVector(this->_space->dimension());
    for (std::size_t v = 0; v < this->_space->dimension(); ++v)
      convertedObjectiveVector[v] = objectiveVector[v];
    return maximize(&convertedObjectiveVector[0], query);
  }

  template <>
  OptimizationResponse<rational> TrustRegionOptimizationOracle<rational>::maximizeTrustRegionDouble(
    const ManhattanTrustRegion<rational>& trustRegion, const double* objectiveVector,
    const OptimizationQuery<rational>& query)
  {
    std::vector<rational> convertedObjectiveVector(this->_space->dimension());
    for (std::size_t v = 0; v < this->_space->dimension(); ++v)
      convertedObjectiveVector[v] = objectiveVector[v];
    return maximizeTrustRegion(trustRegion, &convertedObjectiveVector[0], query);
  }

#endif /* IPO_WITH_RATIONAL */

  template <typename R>
  static std::ostream& printOptimizationResponse(std::ostream& stream,
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
        << " points, opt in [";
      if (response.hasPrimalBound())
        stream << convertNumber<double>(response.primalBound()) << ",";
      else
        stream << "-inf,";
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
    return printOptimizationResponse(stream, response);
  }

#if defined(IPO_WITH_RATIONAL)

  std::ostream& operator<<(std::ostream& stream, const OptimizationResponse<rational>& response)
  {
    return printOptimizationResponse(stream, response);
  }

#endif /* IPO_WITH_RATIONAL */

  template <typename Number>
  SeparationOracle<Number>::SeparationOracle()
    : Oracle<Number>("abstract SeparationOracle")
  {

  }

  template <typename Number>
  SeparationOracle<Number>::~SeparationOracle()
  {

  }

  template <typename Number>
  SeparationResponse<Number> SeparationOracle<Number>::getInitial(
    const SeparationQuery& query)
  {
    return SeparationResponse<Number>();
  }

  template <>
  SeparationResponse<double> SeparationOracle<double>::separateDouble(const double* vector,
    bool isPoint, const Query& query)
  {
    return separate(vector, isPoint, query);
  }

#if defined(IPO_WITH_RATIONAL)

  template <>
  SeparationResponse<rational> SeparationOracle<rational>::separateDouble(const double* vector,
    bool isPoint, const Query& query)
  {
    std::vector<rational> rationalVector(_space->dimension());
    for (std::size_t v = 0; v < _space->dimension(); ++v)
      rationalVector[v] = reconstructRational(vector[v]);
    return separate(&rationalVector[0], isPoint, query);
  }

#endif /* IPO_WITH_RATIONAL */

  template class Oracle<double>;
  template class OptimizationOracle<double>;
  template class TrustRegionOptimizationOracle<double>;
  template class SeparationOracle<double>;

#if defined(IPO_WITH_RATIONAL)

  template class Oracle<rational>;
  template class OptimizationOracle<rational>;
  template class TrustRegionOptimizationOracle<rational>;
  template class SeparationOracle<rational>;

#endif /* IPO_WITH_RATIONAL */

} /* namespace ipo */
