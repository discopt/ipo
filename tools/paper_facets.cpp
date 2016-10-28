
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
#include "ipo/facets.h"
#include "ipo/scip_oracles.h"
#include "ipo/cache_oracle.h"
#include "ipo/statistics_oracle.h"

using namespace ipo;
using namespace soplex;

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

  std::shared_ptr<CacheOracle> cacheOracle = std::make_shared<CacheOracle>(scipOracleStats, CacheOracle::CACHE_AND_SEARCH);
  std::shared_ptr<StatisticsOracle> cacheOracleStats = std::make_shared<StatisticsOracle>(cacheOracle);

  std::shared_ptr<OracleBase> oracle = cacheOracleStats;

  std::vector<AffineHullHandler*> affineHullHandlers;
  InnerDescription inner;
  AffineOuterDescription outer;
  affineHull(oracle, inner, outer, affineHullHandlers, 2, 0);
  std::cout << "Dimension: " << (long(inner.points.size() + inner.rays.size()) - 1) << std::endl;

//   exactScipOracleStats->reset();
  scipOracleStats->reset();
  cacheOracleStats->reset();

  std::vector<FacetSeparationHandler*> facetSeparationHandlers;
  DebugFacetSeparationHandler debugHandler(std::cout, true, true);
  StatisticsFacetSeparationHandler statsHandler;
  facetSeparationHandlers.push_back(&debugHandler);
  facetSeparationHandlers.push_back(&statsHandler);

  SoPlex spx;
  spx.setIntParam(SoPlex::SOLVEMODE, SoPlex::SOLVEMODE_RATIONAL);
  spx.setIntParam(SoPlex::SYNCMODE, SoPlex::SYNCMODE_AUTO);
  spx.setRealParam(SoPlex::FEASTOL, 0.0);
  spx.setBoolParam(SoPlex::RATREC, true);
  spx.setBoolParam(SoPlex::RATFAC, true);
  spx.setIntParam(SoPlex::OBJSENSE, SoPlex::OBJSENSE_MAXIMIZE);
  spx.setIntParam(SoPlex::VERBOSITY, SoPlex::VERBOSITY_ERROR);

  const MixedIntegerSet& mis = *scipOracle->mixedIntegerSet();
  LPColSetRational cols(mis.numVariables());
  DSVectorRational zero;
  for (std::size_t v = 0; v < oracle->space().dimension(); ++v)
  {
    cols.add(Rational(1), mis.variable(v).lowerBound, zero, mis.variable(v).upperBound);
  }
  spx.addColsRational(cols);
  addToLP(spx, mis.rowConstraints());

  DVectorRational solution(mis.numVariables());
  while (true)
  {
    std::cout << "Solving relaxation LP with " << spx.numRowsRational() << " rows." << std::endl;

    SPxSolver::Status status = spx.solve();
    if (status == SPxSolver::OPTIMAL)
    {
      spx.getPrimalRational(solution);
      ipo::Vector point = denseToVector(solution, false);
      InnerDescription certificate;
      LinearConstraint constraint;
      if (separatePoint(oracle, point, inner, facetSeparationHandlers, constraint, &certificate))
      {
        scaleIntegral(constraint);

        std::cout << "Separated point with ";
        oracle->space().printLinearConstraint(std::cout, constraint);
        std::cout << std::endl;

        addToLP(spx, constraint);
      }
      else
        break;

    }
    else
      assert(false);
  }

  std::cout << "\n";
  std::cout << "Algorithm statistics (without affine hull computation):\n";
  std::cout << "\n";
  std::cout << "Overall time: " << statsHandler.timeAll() << "\n";
  std::cout << "Approximate LPs: " << statsHandler.numApproximateLPs() << " in " << statsHandler.timeApproximateLPs() 
    << " seconds.\n";
  std::cout << "Exact LPs: " << statsHandler.numExactLPs() << " in " << statsHandler.timeExactLPs()
    << " seconds.\n";
  std::cout << "Oracle queries: " << statsHandler.numOracleQueries() << " in " << statsHandler.timeOracles() 
    << " seconds.\n";
  std::cout << "\n";
  std::cout << "Oracle statistics (without affine hull computation):\n";
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
