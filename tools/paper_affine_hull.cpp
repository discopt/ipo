#include <ipo/affine_hull.h>
#include <ipo/scip_oracle.h>
#include <ipo/scip_exception.h>
#include <ipo/exactscip_oracle.h>
#include <ipo/cache_oracle.h>
#include <ipo/statistics_oracle.h>
#include <ipo/polyhedron.h>

#ifdef IPO_WITH_SCIP
#ifdef NDEBUG
  #undef NDEBUG
  #include <scip/scip.h>
  #include <scip/scipdefplugins.h>
  #include <scip/cons_linear.h>
  #define NDEBUG
#else
  #include <scip/scip.h>
  #include <scip/scipdefplugins.h>
  #include <scip/cons_linear.h>
#endif
#endif

using namespace ipo;

int printUsage(const std::string& program)
{
  std::cerr << "Usage: " << program << " [-o|-c] [-t ORACLE-TIMELIMIT] [-d DIRECTION-TIMELIMIT] INSTANCE\n" << std::flush;
  return EXIT_FAILURE;
}

int main(int argc, char** argv)
{
  // Parse arguments.

  bool constraintDimensions = false;
  bool optimalFaceDimension = false;
  std::string instanceFile = "";
  double oracleTimeLimit = -1;
  double directionTimeLimit = std::numeric_limits<double>::max();
  for (int i = 1; i < argc; ++i)
  {
    if (std::string(argv[i]) == "-h")
      return printUsage(argv[0]);
    if (std::string(argv[i]) == "-c")
    {
      constraintDimensions = true;
      continue;
    }
    if (std::string(argv[i]) == "-o")
    {
      optimalFaceDimension = true;
      continue;
    }
    if (std::string(argv[i]) == "-t" && i+1 < argc)
    {
      std::stringstream str(argv[i+1]);
      str >> oracleTimeLimit;
      ++i;
      continue;
    }
    if (std::string(argv[i]) == "-d" && i+1 < argc)
    {
      std::stringstream str(argv[i+1]);
      str >> directionTimeLimit;
      ++i;
      continue;
    }
    if (instanceFile.empty())
      instanceFile = argv[i];
    else
      return printUsage(argv[0]);
  }
  if (instanceFile.empty() || (optimalFaceDimension && constraintDimensions))
    return printUsage(argv[0]);

  // Read instance and create MixedIntegerSet.

  std::cout << "Creating oracles for mixed-integer set defined by " << instanceFile << ".\n" << std::flush;

  SCIP* scip = NULL;
  SCIP_CALL_EXC(SCIPcreate(&scip));
  SCIP_CALL_EXC(SCIPincludeDefaultPlugins(scip));
  SCIP_CALL_EXC(SCIPsetIntParam(scip, "display/verblevel", 0));
  SCIP_CALL_EXC(SCIPreadProb(scip, instanceFile.c_str(), NULL));
  SCIP_CALL_EXC(SCIPtransformProb(scip));

  std::shared_ptr<MixedIntegerLinearSet> mixedIntegerLinearSet = std::make_shared<MixedIntegerLinearSet>(scip);

  Vector instanceObjective = getSCIPObjective(scip);

  SCIP_CALL_EXC(SCIPfree(&scip));

  // Initialize oracles.

#ifdef IPO_WITH_EXACT_SCIP
  std::shared_ptr<ExactSCIPOracle> exactSCIPOracle = std::make_shared<ExactSCIPOracle>(
    "ExactSCIPOracle(" + instanceFile + ")", mixedIntegerLinearSet);
  std::shared_ptr<StatisticsOracle> exactSCIPOracleStats = std::make_shared<StatisticsOracle>(exactSCIPOracle);

  std::shared_ptr<SCIPOracle> scipOracle = std::make_shared<SCIPOracle>("SCIPOracle(" + instanceFile + ")",
    mixedIntegerLinearSet, exactSCIPOracleStats);
#else
  std::shared_ptr<SCIPOracle> scipOracle = std::make_shared<SCIPOracle>("SCIPOracle(" + instanceFile + ")",
    mixedIntegerLinearSet);
#endif

  if (oracleTimeLimit > 0)
  {
    scipOracle->setTimeLimit(oracleTimeLimit);
#ifdef IPO_WITH_EXACT_SCIP
    exactSCIPOracle->setTimeLimit(oracleTimeLimit);
#endif
  }

  std::shared_ptr<StatisticsOracle> scipOracleStats = std::make_shared<StatisticsOracle>(scipOracle);

  std::shared_ptr<CacheOracle> cacheOracle = std::make_shared<CacheOracle>(scipOracleStats);
  std::shared_ptr<StatisticsOracle> cacheOracleStats = std::make_shared<StatisticsOracle>(cacheOracle);

  std::shared_ptr<OracleBase> oracle = cacheOracleStats;

  LinearConstraint faceConstraint = completeFaceConstraint();
  if (optimalFaceDimension)
  {
    OracleResult result;
    cacheOracleStats->maximize(result, instanceObjective);
    if (result.isInfeasible() || result.isUnbounded())
    {
      std::cout << "Dimension of optimal face: -1\n" << std::flush;
      return EXIT_SUCCESS;
    }

    faceConstraint = LinearConstraint('<', instanceObjective, result.objectiveValue());
  }
  
  Polyhedron poly(cacheOracleStats);
  poly.addConstraint(faceConstraint);

  // Set affine-hull parameters.

#ifdef IPO_WITH_EXACT_SCIP
  poly.setAffineHullLastModerateHeuristicLevel(1); // Use oracle 0 only for verification.
  poly.setAffineHullLastCheapHeuristicLevel(2); // Caching is cheap.
#else
  poly.setAffineHullLastModerateHeuristicLevel(0); // No verification.
  poly.setAffineHullLastCheapHeuristicLevel(1); // Caching is cheap.
#endif
  poly.setAffineHullExactDirectionTimeLimit(directionTimeLimit);

  for (std::size_t r = 0; r < mixedIntegerLinearSet->numRows(); ++r)
  {
    const LinearConstraint& row = mixedIntegerLinearSet->rowConstraint(r);
    poly.addConstraint(row);
  }
  for (std::size_t v = 0; v < mixedIntegerLinearSet->numVariables(); ++v)
  {
    LinearConstraint constraint;
    constraint = mixedIntegerLinearSet->lowerBoundConstraint(v);
    if (!constraint.isEquation() && constraint.rhs() > -soplex::infinity)
      poly.addConstraint(constraint);
    constraint = mixedIntegerLinearSet->upperBoundConstraint(v);
    if (!constraint.isEquation() && constraint.rhs() < soplex::infinity)
      poly.addConstraint(constraint);
  }

  StatisticsAffineHullHandler statsHandler;
  DebugAffineHullHandler debugHandler(std::cout);
  
  std::vector<LinearConstraint> givenEquations;
  mixedIntegerLinearSet->getConstraints(givenEquations, false, true, false, true);

  std::vector<AffineHullHandler*> handlers;
  handlers.push_back(&statsHandler);
  if (constraintDimensions)
    std::cout << "Dimension: " << std::flush;
  else
    handlers.push_back(&debugHandler);
  
  try
  {
    std::shared_ptr<Polyhedron::Face> face = poly.constraintToFace(faceConstraint);
    poly.affineHull(face, handlers, givenEquations);
  }
  catch(std::exception& e)
  {
    std::cout << "Error: " << e.what() << "\n" << std::flush;
    return EXIT_FAILURE;
  }
  int dim = poly.dimension();
  if (constraintDimensions)
  {
    std::cout << dim << "\n\n" << std::flush;

    for (std::size_t v = 0; v < mixedIntegerLinearSet->numVariables(); ++v)
    {
      const LinearConstraint& constraint = mixedIntegerLinearSet->lowerBoundConstraint(v);
      if (constraint.isEquation() || constraint.rhs() == -soplex::infinity)
        continue;
      std::shared_ptr<Polyhedron::Face> face = poly.constraintToFace(constraint);
      std::cout << "Computing dimension of face defined by lower bound ";
      poly.space().printLinearConstraint(std::cout, constraint);
      std::cout << ": " << std::flush;
      poly.affineHull(face, handlers);
      std::cout << face->dimension() << std::endl;
    }

    for (std::size_t v = 0; v < mixedIntegerLinearSet->numVariables(); ++v)
    {
      const LinearConstraint& constraint = mixedIntegerLinearSet->upperBoundConstraint(v);
      if (constraint.isEquation() || constraint.rhs() == soplex::infinity)
        continue;
      std::shared_ptr<Polyhedron::Face> face = poly.constraintToFace(constraint);
      std::cout << "Computing dimension of face defined by upper bound ";
      poly.space().printLinearConstraint(std::cout, constraint);
      std::cout << ": " << std::flush;
      poly.affineHull(face, handlers);
      std::cout << face->dimension() << std::endl;
    }

    for (std::size_t r = 0; r < mixedIntegerLinearSet->numRows(); ++r)
    {
      const LinearConstraint& row = mixedIntegerLinearSet->rowConstraint(r);
      std::shared_ptr<Polyhedron::Face> face = poly.constraintToFace(row);
      if (row.isEquation())
        continue;

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

  soplex::Rational::freeListMem();
  ipo::Space::freeStaticMem();
  ipo::Vector::freeStaticMem();

  return 0;
}
