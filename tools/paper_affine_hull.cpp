
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
#include "ipo/polyhedron.h"

using namespace ipo;

int main(int argc, char** argv)
{
  // Parse arguments.
  
  // Read instance and create MixedIntegerSet.

  SCIP* scip = NULL;
  SCIP_CALL_EXC(SCIPcreate(&scip));
  SCIP_CALL_EXC(SCIPincludeDefaultPlugins(scip));
  SCIP_CALL_EXC(SCIPsetIntParam(scip, "display/verblevel", 0));
  SCIP_CALL_EXC(SCIPreadProb(scip, argv[1], NULL));
  SCIP_CALL_EXC(SCIPtransformProb(scip));

  std::shared_ptr<MixedIntegerSet> mixedIntegerSet = std::make_shared<MixedIntegerSet>(scip);

  SCIP_CALL_EXC(SCIPfree(&scip));

  // Initialize oracles.

  std::shared_ptr<ExactSCIPOracle> exactSCIPOracle = std::make_shared<ExactSCIPOracle>(
    "ExactSCIPOracle(" + std::string(argv[1]) + ")", mixedIntegerSet);
//****************************************************??
  exactSCIPOracle->setBinaryPath("/home/matthias/software/exactscip/scip-3.0.0-ex/bin/scip");
//****************************************************
  std::shared_ptr<StatisticsOracle> exactScipOracleStats = std::make_shared<StatisticsOracle>(exactSCIPOracle);

  std::shared_ptr<SCIPOracle> scipOracle = std::make_shared<SCIPOracle>("SCIPOracle(" + std::string(argv[1]) + ")",
    mixedIntegerSet, exactScipOracleStats);
  std::shared_ptr<StatisticsOracle> scipOracleStats = std::make_shared<StatisticsOracle>(scipOracle);

  std::shared_ptr<CacheOracle> cacheOracle = std::make_shared<CacheOracle>(scipOracleStats);
  std::shared_ptr<StatisticsOracle> cacheOracleStats = std::make_shared<StatisticsOracle>(cacheOracle);

  std::shared_ptr<OracleBase> oracle = cacheOracleStats;

  Polyhedron poly(cacheOracleStats);
  poly.setAffineHullLastModerateHeuristicLevel(1); // Use oracle 0 only for verification.
  poly.setAffineHullLastCheapHeuristicLevel(2); // Alternate min/max for oracle 2 before using oracle 1.

  for (std::size_t r = 0; r < mixedIntegerSet->numRows(); ++r)
  {
    const LinearConstraint& row = mixedIntegerSet->rowConstraint(r);
    poly.addConstraint(row);
  }

  StatisticsAffineHullHandler statsHandler;
  std::vector<AffineHullHandler*> handlers;
  handlers.push_back(&statsHandler);

  std::cout << "Dimension:\n" << std::flush;
  poly.affineHull(handlers);
  int dim = poly.dimension();
  std::cout << dim << "\n\n" << std::flush;

  // Compute constraint dimensions.

  for (std::size_t r = 0; r < mixedIntegerSet->numRows(); ++r)
  {
    const LinearConstraint& row = mixedIntegerSet->rowConstraint(r);
    std::shared_ptr<Polyhedron::Face> face = poly.inequalityToFace(row);

    if (!row.isEquation())
    {
      std::cout << "Computing dimension of face defined by ";
      poly.space().printLinearConstraint(std::cout, row);
      std::cout << ": " << std::flush;
      poly.affineHull(face, handlers);
      std::cout << face->dimension() << std::endl;
    }
  }

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

  // Output constraint dimensions.
  
  std::vector<std::shared_ptr<Polyhedron::Face> > faces;
  poly.getFaces(faces, true, true);
  for (std::size_t i = 0; i < faces.size(); ++i)
  {
    std::cout << "Constraint #" << i << ": ";
    poly.space().printLinearConstraint(std::cout, faces[i]->inequality());
    std::cout << " has dimension " << faces[i]->dimension() << "\n" << std::flush;
  }

  return 0;
}
