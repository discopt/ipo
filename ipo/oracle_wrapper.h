#ifndef IPO_ORACLE_WRAPPER_H_
#define IPO_ORACLE_WRAPPER_H_

#include "common.h"

#include "oracles.h"
#include "cache_oracle.h"

namespace ipo {

  class StatisticsCacheOracleWrapper
  {
    typedef CacheOracle::Behavior Behavior;

  public:
    StatisticsCacheOracleWrapper(std::shared_ptr<OracleBase>& mainOracle, Behavior outerBehavior = Behavior::CACHE_ONLY,
      Behavior innerBehavior = Behavior::CACHE_AND_SEARCH);

    ~StatisticsCacheOracleWrapper();

    inline const Space& space() const
    {
      return _mainOracle->space();
    }

    inline std::shared_ptr<OracleBase> linkOracle()
    {
      return _mainStatisticsOracle;
    }

    inline std::shared_ptr<OracleBase> queryOracle()
    {
      return _cacheStatisticsOracle;
    }

    inline std::shared_ptr<OracleBase> getMainOracle()
    {
      return _mainOracle;
    }

    inline std::shared_ptr<OracleBase> getCacheOracle()
    {
      return _cacheOracle;
    }

  protected:

    std::shared_ptr<OracleBase> _mainOracle;
    std::shared_ptr<OracleBase> _mainStatisticsOracle;
    std::shared_ptr<OracleBase> _cacheOracle;
    std::shared_ptr<OracleBase> _cacheStatisticsOracle;
  };

  typedef StatisticsCacheOracleWrapper DefaultOracleWrapper;

}

#endif /* IPO_ORACLE_WRAPPER_H_ */
