#include <ipo/common.h>

#include "ipo/scip_exception.hpp"
#include "ipo/affine_hull.h"
#include "ipo/scip_oracles.h"
#include "ipo/cache_oracle.h"
#include "ipo/statistics_oracle.h"
#include <ipo/vectors.h>
using namespace ipo;

//Output Modes for DebugAffineHullHandler
#define STDOUT 1

class ScipOracleController
{
  public:
  //Constructors
  //creates ScipOracle with Statistics Wrapper and Cache
  ScipOracleController(std::string name);
  ScipOracleController(std::string name, ScipOracleController prev);

  ~ScipOracleController();

  //Functions
  std::string name();
  int heuristicLevel_ScipOracle();
  int heuristicLevel_CacheOracle();
  Space space();
  InnerDescription affineHullInner(int outputMode);
  AffineOuterDescription affineHullOuter(int outputMode);
  std::shared_ptr<StatisticsOracle> getConnectionOracle();
  
  

  private:
  std::shared_ptr<StatisticsOracle> scipOracle;
  std::shared_ptr<SCIPOracle> scipOracleImpl;
  std::shared_ptr<CacheOracle> cacheOracleImpl;
  std::shared_ptr<StatisticsOracle> cacheOracle;
};
