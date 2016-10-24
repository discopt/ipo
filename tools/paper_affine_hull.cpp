
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

using namespace ipo;

int main(int argc, char** argv)
{
  // Read instance and create MixedIntegerSet.

  SCIP* scip = NULL;
  SCIP_CALL_EXC(SCIPcreate(&scip));
  SCIP_CALL_EXC(SCIPincludeDefaultPlugins(scip));
  SCIP_CALL_EXC(SCIPsetIntParam(scip, "display/verblevel", 0));
  SCIP_CALL_EXC(SCIPreadProb(scip, argv[1], NULL));
  SCIP_CALL_EXC(SCIPtransformProb(scip));

  std::shared_ptr<MixedIntegerSet> mixedIntegerSet= std::make_shared<MixedIntegerSet>(scip);

  SCIP_CALL_EXC(SCIPfree(&scip));

  // Initialize oracles.

//   std::shared_ptr<ExactSCIPOracle> exactSCIPOracle = std::make_shared<ExactSCIPOracle>(
//     "ExactSCIPOracle(" + std::string(argv[1]) + ")", mixedIntegerSet);
//   exactSCIPOracle->setBinaryPath("/home/matthias/software/exactscip/scip-3.0.0-ex/bin/scip");
//   std::shared_ptr<StatisticsOracle> exactScipOracleStats = std::make_shared<StatisticsOracle>(exactSCIPOracle);

  std::shared_ptr<SCIPOracle> scipOracle = std::make_shared<SCIPOracle>("SCIPOracle(" + std::string(argv[1]) + ")",
    mixedIntegerSet/*, exactScipOracleStats*/);
  std::shared_ptr<StatisticsOracle> scipOracleStats = std::make_shared<StatisticsOracle>(scipOracle);

  std::shared_ptr<CacheOracle> cacheOracle = std::make_shared<CacheOracle>(scipOracleStats);
  std::shared_ptr<StatisticsOracle> cacheOracleStats = std::make_shared<StatisticsOracle>(cacheOracle);

  std::shared_ptr<OracleBase> oracle = cacheOracleStats;

  std::vector<AffineHullHandler*> handlers;
  DebugAffineHullHandler debugHandler(std::cout);
  StatisticsAffineHullHandler statsHandler;
  handlers.push_back(&debugHandler);
  handlers.push_back(&statsHandler);

  InnerDescription inner;
  AffineOuterDescription outer;
  affineHull(oracle, inner, outer, handlers, 2, 1);

  std::cout << "\n";
  std::cout << "Algorithm statistics:\n";
  std::cout << "\n";
  std::cout << "Overall time: " << statsHandler.timeAll() << "  =  main loop time: " << statsHandler.timeMainLoop()
    << "  +  verification time: " << statsHandler.timeVerification() << "\n";
  std::cout << "Approximate directions: " << statsHandler.numDirectionApproximateSolves() << " in " <<
    statsHandler.timeApproximateDirections() << " seconds.\n";
  std::cout << "Exact directions: " << statsHandler.numDirectionExactSolves() << " in " <<
    statsHandler.timeExactDirections() << " seconds.\n";
  std::cout << "Factorizations: " << statsHandler.numFactorizations() << " in " << statsHandler.timeFactorizations()
    << " seconds.\n";
  std::cout << "Oracle queries: " << statsHandler.numOracleQueries() << " in " << statsHandler.timeOracles()
    << " seconds.\n";
  std::cout << "\n";
  std::cout << "Oracle statistics:\n";
  std::cout << "\n";
  for (std::shared_ptr<OracleBase> o = oracle; o != NULL; o = o->nextOracle())
  {
    std::size_t h = o->heuristicLevel();
    std::cout << "Oracle " << h << ": " << o->name() << "\n";
    std::shared_ptr<StatisticsOracle> s = std::dynamic_pointer_cast<StatisticsOracle>(o);
    std::cout << "  #calls:   " << s->numCalls() << "\n";
    std::cout << "  #success: " << s->numSuccess() << "\n";
    std::cout << "  time:     " << s->time() << "\n";
  }
  std::cout << std::endl;

  return 0;
}
