#define IPO_DEBUG /* Uncomment to debug this file. */

#include <ipo/oracles_radialcone.hpp>

#include <chrono>

namespace ipo
{

  template <typename Number>
  RadialConeOptimizationOracle<Number>::RadialConeOptimizationOracle(
    std::shared_ptr<OptimizationOracle<Number>> sourceOracle, std::shared_ptr<sparse_vector<Number>> apex,
    const std::string& name)
    : Oracle<Number>(name.empty() ? ("RadialCone(" + sourceOracle->name() + ")") : name),
    _sourceOracle(sourceOracle),
    _trustRegionOracle(dynamic_cast<TrustRegionOptimizationOracle<Number>*>(sourceOracle.get())), _apex(apex),
    _trustRegion(), _iteration(0), _trustRegionResetIterations(0), _trustRegionInitialSize(16),
    _trustRegionMaximumSize(0), _trustRegionGrowthRate(2)
  {
    this->_space = sourceOracle->space();
  }

  template <typename Number>
  RadialConeOptimizationOracle<Number>::~RadialConeOptimizationOracle()
  {

  }

  template <typename Number>
  bool RadialConeOptimizationOracle<Number>::isTrustRegionCapable() const
  {
    return _trustRegionOracle != nullptr;
  }

  template <typename Number>
  void RadialConeOptimizationOracle<Number>::enableManhattanTrustRegion(const Number& initialSize,
    const Number& maximumSize, const Number& growthRate, std::size_t resetIterations)
  {
    if (!isTrustRegionCapable())
      throw std::runtime_error("Source oracle of RadialConeOptimizationOracle is not a TrustRegionOptimizationOracle.");
    if (growthRate <= 1)
      throw std::runtime_error("Growth rate of trust region must be greater than 1.");

    _iteration = 0;
    _trustRegionInitialSize = initialSize;
    _trustRegionMaximumSize = maximumSize;
    _trustRegionGrowthRate = growthRate;
    _trustRegionResetIterations = resetIterations;

    if ((_trustRegionEnabled = (initialSize <= maximumSize)))
      _trustRegion = ManhattanTrustRegion<Number>(_apex, initialSize);
    else
      _trustRegion = ManhattanTrustRegion<Number>();
  }

  template <typename Number>
  void RadialConeOptimizationOracle<Number>::disableManhattanTrustRegion()
  {
    _trustRegion = ManhattanTrustRegion<Number>();
  }

  template <typename Number>
  OptimizationResponse<Number> RadialConeOptimizationOracle<Number>::maximize(const Number* objectiveVector,
      const OptimizationQuery<Number>& query)
  {
    OptimizationResponse<Number> response;
    Number apexObjective = objectiveVector * *_apex;
    std::chrono::time_point<std::chrono::system_clock> started = std::chrono::high_resolution_clock::now();
    if (_trustRegionResetIterations && (((_iteration + 1) % _trustRegionResetIterations) == 0))
    {
#if defined(IPO_DEBUG)
      std::cout << this->name() << ": resetting trust region size to " << _trustRegionInitialSize << std::endl;
#endif /* IPO_DEBUG */
      Number newSize = _trustRegionInitialSize;
      if (newSize > _trustRegionMaximumSize)
        _trustRegion = ManhattanTrustRegion<Number>();
      else
        _trustRegion.updateSize(newSize);
    }

    while (true)
    {
      OptimizationResponse<Number> sourceResponse;
      OptimizationQuery<Number> sourceQuery = query;

      sourceQuery.timeLimit = std::max(0.0, query.timeLimit
        - std::chrono::duration<double>(std::chrono::system_clock::now() - started).count());

      if (_trustRegionEnabled)
      {
#if defined(IPO_DEBUG)
        std::cout << this->name() << ": calling "<< _sourceOracle->name() << " with trust region size "
          << _trustRegion.size() << "." << std::endl;
#endif /* IPO_DEBUG */

        sourceResponse = _trustRegionOracle->maximizeTrustRegion(_trustRegion, objectiveVector, sourceQuery);
      }
      else
      {
#if defined(IPO_DEBUG)
        std::cout << this->name() << ": calling " << _sourceOracle->name() << "." << std::endl;
#endif /* IPO_DEBUG */

        sourceResponse = _sourceOracle->maximize(objectiveVector, sourceQuery);
      }

#if defined(IPO_DEBUG)
      std::cout << this->name() << ": received response " << sourceResponse
        << std::endl;
#endif /* IPO_DEBUG */


      response.rays = sourceResponse.rays;
      for (auto point : sourceResponse.points)
      {
        auto ray = *point.vector - *_apex;
        if ((objectiveVector * *ray) > 0)
          response.rays.push_back(typename OptimizationResponse<Number>::Ray(ray));
      }

      if (!_trustRegionEnabled || sourceResponse.hitTimeLimit || !response.rays.empty())
      {
        // We found some rays.
        response.points.push_back(typename OptimizationResponse<Number>::Point(_apex, apexObjective));
        response.outcome = sourceResponse.outcome;
        if (!response.rays.empty())
          response.outcome = OptimizationOutcome::UNBOUNDED;
        response.hitTimeLimit = sourceResponse.hitTimeLimit;
        if (sourceResponse.hasDualBound && sourceResponse.dualBound <= apexObjective)
        {
          response.hasDualBound = true;
          response.dualBound = apexObjective;
        }
        else
          response.hasDualBound = false;

        _iteration++;
        return response;
      }
      else
      {
        // We did not hit the time limit, but found no rays, so we update the trust region.

        Number newSize = _trustRegion.size() * _trustRegionGrowthRate;
        if (newSize > _trustRegionMaximumSize)
          _trustRegion = ManhattanTrustRegion<Number>();
        else
          _trustRegion.updateSize(newSize);

#if defined(IPO_DEBUG)
        std::cout << this->name() << ": updating trust region";
        if (_trustRegion.isBounded())
          std::cout << " size to " << _trustRegion.size() << "." << std::endl;
        else
          std::cout << " to complete space." << std::endl;
#endif /* IPO_DEBUG */
      }
    }
  }

  template class RadialConeOptimizationOracle<double>;

#if defined(IPO_WITH_RATIONAL)

  template class RadialConeOptimizationOracle<rational>;

#endif /* IPO_WITH_RATIONAL */

} /* namespace ipo */
