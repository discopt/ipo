#include <ipo/common.h>

#ifdef WITH_SCIP
#ifdef NDEBUG
  #undef NDEBUG
  #include <scip/scip.h>
  #include <scip/scipdefplugins.h>
  #include <scip/cons_linear.h>
  #include <ipo/scip_exception.hpp>
  #include <ipo/scip_oracles.h>
  #define NDEBUG
#else
  #include <scip/scip.h>
  #include <scip/scipdefplugins.h>
  #include <scip/cons_linear.h>
  #include <ipo/scip_exception.hpp>
  #include <ipo/scip_oracles.h>
#endif
#endif

#include "ipo/scip_exception.hpp"
#include "ipo/affine_hull.h"
#include "ipo/scip_oracles.h"
#include "ipo/cache_oracle.h"
#include "ipo/statistics_oracle.h"
#include "ipo/vectors-pub.h"
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
  InnerDescription affineHullInner();
  AffineOuterDescription affineHullOuter();
  std::shared_ptr<StatisticsOracle> getConnectionOracle();
  
  

  private:
  std::shared_ptr<StatisticsOracle> scipOracle;
  std::shared_ptr<SCIPOracle> scipOracleImpl;
  std::shared_ptr<CacheOracle> cacheOracleImpl;
  std::shared_ptr<StatisticsOracle> cacheOracle;
};
