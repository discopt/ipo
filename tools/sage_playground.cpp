
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
  //////// Oracle with heuristicLevel 0 ////////

  // scipOracle = ipo.SCIPOracle(sys.argv[1])

  SCIP* scip = NULL;
  SCIP_CALL_EXC(SCIPcreate(&scip));
  SCIP_CALL_EXC(SCIPincludeDefaultPlugins(scip));
  SCIP_CALL_EXC(SCIPsetIntParam(scip, "display/verblevel", 0));
  SCIP_CALL_EXC(SCIPreadProb(scip, argv[1], NULL));
  SCIP_CALL_EXC(SCIPtransformProb(scip));
  std::shared_ptr<SCIPOracle> scipOracleImpl = std::make_shared<SCIPOracle>("SCIPOracle(" + std::string(argv[1]) + ")", scip);
  SCIP_CALL_EXC(SCIPfree(&scip));
  std::shared_ptr<StatisticsOracle> scipOracle = std::make_shared<StatisticsOracle>(scipOracleImpl);

  // print("%s has heuristicLevel %d" % (scipOracle.name(), scipOracle.heuristicLevel()))

  std::cout << scipOracle->name() << " has heuristicLevel " << scipOracle->heuristicLevel() << std::endl;

  //////// Oracle with heuristicLevel 1 ////////

  // heuristicSCIPOracle = ipo.SCIPOracle(sys.argv[1], scipOracle) <-- Letzter Konstruktor-Parameter = nächstes Oracle

  SCIP_CALL_EXC(SCIPcreate(&scip));
  SCIP_CALL_EXC(SCIPincludeDefaultPlugins(scip));
  SCIP_CALL_EXC(SCIPsetIntParam(scip, "display/verblevel", 0));
  SCIP_CALL_EXC(SCIPsetLongintParam(scip, "limits/nodes", 1));
  SCIP_CALL_EXC(SCIPreadProb(scip, argv[1], NULL));
  SCIP_CALL_EXC(SCIPtransformProb(scip));
  std::shared_ptr<SCIPOracle> heuristicSCIPOracleImpl = std::make_shared<SCIPOracle>("HeuristicSCIPOracle(" +
std::string(argv[1]) + ")", scip, scipOracle);
  SCIP_CALL_EXC(SCIPfree(&scip));
  std::shared_ptr<StatisticsOracle> heuristicSCIPOracle = std::make_shared<StatisticsOracle>(heuristicSCIPOracleImpl);

  // print("%s has heuristicLevel %d" % (heuristicSCIPOracle.name(), heuristicSCIPOracle.heuristicLevel()))

  std::cout << heuristicSCIPOracle->name() << " has heuristicLevel " << heuristicSCIPOracle->heuristicLevel() << std::endl;

  //////// Oracle with heuristicLevel 2 ////////

  // cacheOracle = ipo.CacheOracle(scipOracle)

  std::shared_ptr<CacheOracle> cacheOracleImpl = std::make_shared<CacheOracle>(heuristicSCIPOracle);
  std::shared_ptr<StatisticsOracle> cacheOracle = std::make_shared<StatisticsOracle>(cacheOracleImpl);

  // print("%s has heuristicLevel %d" % (cacheOracle.name(), cacheOracle.heuristicLevel()))

  std::cout << cacheOracle->name() << " has heuristicLevel " << cacheOracle->heuristicLevel() << std::endl;

  // space = scipOracle.space()

  Space space = cacheOracle->space();

  // mixedIntegerSet = scipOracle.mixedIntegerSet()

//   std::shared_ptr<MixedIntegerSet> mixedIntegerSet = scipOracleImpl->mixedIntegerSet();

  // print('MIP is %d x %d.' % (mixedIntegerSet.numRows(), mixedIntegerSet.numVariables()))

//   std::cout << "MIP is " << mixedIntegerSet->numRows() << " x " << mixedIntegerSet->numVariables() << "." << std::endl;

  // (inner, outer) = ipo.affineHull(cacheOracle, 2, 1)

  std::vector<AffineHullHandler*> handlers; // Soll pro Aufruf lokal erzeugt werden.

  DebugAffineHullHandler debugHandler(std::cout); // Später soll Ausgabedetail konfigurierbar sein.
  StatisticsAffineHullHandler statsHandler;
  handlers.push_back(&debugHandler);
  handlers.push_back(&statsHandler);
  InnerDescription inner; // In Python: Tupel (points, rays), beides jeweils Listen von Vector-Wrappern.
  AffineOuterDescription outer; // In Python: Array von LinearConstraint-Wrappern.
  affineHull(cacheOracle, inner, outer, handlers, 2, 1);

  // print('Dimension = %d' % (len(inner[0]) + len(inner[1]) - 1))

  std::cout << "Dimension = " << (long(inner.points.size() + inner.rays.size())-1) << std::endl;

  // for p in inner[0]: print("Point %s" % (p))
  
  for (std::size_t i = 0; i < inner.points.size(); ++i)
  {
    std::cout << "Point ";
    space.printVector(std::cout, inner.points[i]);
    std::cout << std::endl;
  }

  // for p in inner[0]: print("Point %s" % (p.toSage()))
  // Soll Vektor wie sage ausgeben.

  return 0;
}
