#include "oracle_wrapper.h"

// Uncomment the following line for debugging.
// #define DEBUG

#include "statistics_oracle.h"
#include "cache_oracle.h"

using namespace soplex;

namespace ipo {

  StatisticsCacheOracleWrapper::StatisticsCacheOracleWrapper(std::shared_ptr<OracleBase>& mainOracle, Behavior outerBehavior,
    Behavior innerBehavior)
    : _mainOracle(mainOracle)
  {
    _mainStatisticsOracle = std::make_shared<StatisticsOracle>(_mainOracle);
    _cacheOracle = std::make_shared<CacheOracle>(_mainStatisticsOracle, outerBehavior, innerBehavior);
    _cacheStatisticsOracle = std::make_shared<StatisticsOracle>(_cacheOracle);
  }

  StatisticsCacheOracleWrapper::~StatisticsCacheOracleWrapper()
  {

  }

}
