#include <iostream>
#include <sstream>

#include <ipo/oracles_scip.hpp>
#include <ipo/affine_hull.hpp>

int printUsage(const std::string& program)
{
  std::cout << program << " [OPTIONS] MODEL\n";
  std::cout << "Options:\n";
#if defined(IPO_WITH_GMP)
  std::cout << " --gmp|-g        Use exact arithmetic (GMP).\n";
#endif /* IPO_WITH_GMP */
  std::cout << " --new-equations Print new equations.\n";
  std::cout << " --help|-h       Show this help and exit.\n";
  std::cout << std::flush;
  return EXIT_FAILURE;
}

int main(int argc, char** argv)
{
#if defined(IPO_WITH_GMP)
  bool gmp = false;
#endif /* IPO_WITH_GMP */
  double timeLimit = std::numeric_limits<double>::infinity();
  std::string fileName;
  bool printEquations = false;
  for (int a = 1; a < argc; ++a)
  {
    const std::string arg = argv[a];
    if (arg == "-H" || arg == "-h" || arg == "--help")
    {
      printUsage(argv[0]);
      return EXIT_SUCCESS;
    }
#if defined(IPO_WITH_GMP)
    else if (arg == "-G" || arg == "-g" || arg == "--gmp")
      gmp = true;
#endif /* IPO_WITH_GMP */
    else if ((arg == "-T" || arg == "-t" || arg == "--time") && a+1 < argc)
    {
      std::stringstream ss(argv[a+1]);
      ss >> timeLimit;
      ++a;
    }
    else if (arg == "--print-equations")
      printEquations = true;
    else if (fileName.empty())
      fileName = arg;
    else
    {
      std::cerr << "Invalid parameters. Interpreted <" << fileName << "> and <" << arg
        << "> as instance files." << std::endl;
      return EXIT_FAILURE;
    }
  }

  auto scip = std::make_shared<ipo::SCIPSolver>(fileName);
#if defined(IPO_WITH_GMP)
  if (gmp)
  {
    auto opt = scip->getOptimizationOracleRational();
    auto poly = std::make_shared<ipo::Polyhedron<mpq_class>>(opt);

    // Extract known equations from a separation oracle.
    std::vector<ipo::Constraint<mpq_class>> knownEquations;
    ipo::SeparationOracle<mpq_class>::Query sepaQuery;
    auto sepaResult = scip->getSeparationOracleRational()->getInitial(sepaQuery);
    for (const auto& constraint : sepaResult.constraints)
    {
      if (constraint.type() == ipo::EQUATION)
        knownEquations.push_back(constraint);
    }

    std::cout << "Starting affine hull computation in dimension " << poly->space()->dimension() << std::endl;
    ipo::AffineHullQuery affQuery(false);
    affQuery.timeLimit = timeLimit;
    auto affResult = ipo::affineHull(poly, affQuery, knownEquations);
    std::cout << affResult << "\nAmbient dimension: " << poly->space()->dimension() << std::endl;
    if (printEquations)
    {
      for (auto& equation : affResult.equations)
      {
        bool isKnown = false;
        for (const auto& knownEquation : knownEquations)
        {
          if (knownEquation == equation)
          {
            isKnown = true;
            break;
          }
        }
        if (!isKnown)
          scaleIntegral(equation);
        std::cout << (isKnown ? "Known " : "New ") << "equation "
          << poly->space()->printConstraint(equation, true) << std::endl;
      }
    }
  }
  else
#endif /* IPO_WITH_GMP */
  {
    auto opt = scip->getOptimizationOracleDouble();
    auto poly = std::make_shared<ipo::Polyhedron<double>>(opt);

    // Extract known equations from a separation oracle.
    std::vector<ipo::Constraint<double>> knownEquations;
    ipo::SeparationOracle<double>::Query sepaQuery;
    auto sepaResult = scip->getSeparationOracleDouble()->getInitial(sepaQuery);
    for (const auto& constraint : sepaResult.constraints)
    {
      if (constraint.type() == ipo::EQUATION)
        knownEquations.push_back(constraint);
    }

    std::cout << "Starting affine hull computation in dimension " << poly->space()->dimension() << std::endl;
    ipo::AffineHullQuery affQuery(false);
    affQuery.timeLimit = timeLimit;
    auto affResult = ipo::affineHull(poly, affQuery, knownEquations);
    std::cout << affResult << "\nAmbient dimension: " << poly->space()->dimension() << std::endl;
    if (printEquations)
    {
      for (auto& equation : affResult.equations)
      {
        bool isKnown = false;
        for (const auto& knownEquation : knownEquations)
        {
          if (knownEquation == equation)
          {
            isKnown = true;
            break;
          }
        }
        if (!isKnown)
          scaleIntegral(equation);
        std::cout << (isKnown ? "Known " : "New ") << "equation "
          << poly->space()->printConstraint(equation, true) << std::endl;
      }
    }
  }
 

//   // Parameters
// 
// #ifdef IPO_WITH_EXACT_SCIP
//   bool exactUse = false;
// #endif /* IPO_WITH_EXACT_SCIP */
//   bool affineHullDebug = false;
//   bool affineHullStats = false;
//   bool printEquations = false;
//   bool separationDebug = false;
//   bool separationStats = false;
//   std::size_t numIterations = 100;
//   std::string fileName = "";
// 
//   for (int i = 1; i < argc; ++i)
//   {
//     std::string arg = argv[i];
//     if (arg == "-ad" || arg == "--affinehull-debug")
//       affineHullDebug = true;
// #ifdef IPO_WITH_EXACT_SCIP
//     else if (arg == "-x" || arg == "--exact")
//       exactUse = true;
// #endif /* IPO_WITH_EXACT_SCIP */
//     else if (arg == "-e" || arg == "--equations")
//       printEquations = true;
//     else if (arg == "-as" || arg == "--affinehull-stats")
//       affineHullStats = true;
//     else if (arg == "-d" || arg == "--debug")
//       separationDebug = true;
//     else if (arg == "-s" || arg == "--stats")
//       separationStats = true;
//     else if ((arg == "-i" || arg == "--iterations") && (i + 1 < argc))
//     {
//       std::stringstream str(argv[i+1]);
//       str >> numIterations;
//       ++i;
//     }
//     else if (arg == "-h" || arg == "--help")
//       return printUsage(argv[0]);
//     else if (fileName.empty())
//       fileName = arg;
//     else
//     {
//       std::cout << "Two non-option arguments \"" << fileName << "\" and \"" << arg << "\".\n\n";
//       return printUsage(argv[0]);
//     }
//   }
//   if (fileName.empty())
//   {
//     std::cout << "Missing non-option arguments.\n\n";
//     return printUsage(argv[0]);
//   }

  // Read instance and create MixedIntegerSet.

//   SCIP* scip = NULL;
//   SCIP_CALL_EXC(SCIPcreate(&scip));
//   SCIP_CALL_EXC(SCIPincludeDefaultPlugins(scip));
//   SCIP_CALL_EXC(SCIPsetIntParam(scip, "display/verblevel", 0));
//   SCIP_CALL_EXC(SCIPreadProb(scip, fileName.c_str(), NULL));
//   SCIP_CALL_EXC(SCIPtransformProb(scip));
// 
//   ipo::Vector originalObjective = getSCIPObjective(scip, true);
// 
//   std::shared_ptr<MixedIntegerLinearSet> mixedIntegerSet = std::make_shared<MixedIntegerLinearSet>(scip);
// 
//   SCIP_CALL_EXC(SCIPfree(&scip));
// 
//   // Initialize oracles.
// 
//   std::shared_ptr<SCIPOracle> scipOracle;
// #ifdef IPO_WITH_EXACT_SCIP
//   std::shared_ptr<ExactSCIPOracle> exactScipOracle;
//   std::shared_ptr<StatisticsOracle> exactScipOracleStats;
//   if (exactUse)
//   {
//     exactScipOracle = std::make_shared<ExactSCIPOracle>("ExactSCIPOracle(" + fileName + ")", mixedIntegerSet);
//     exactScipOracleStats = std::make_shared<StatisticsOracle>(exactScipOracle);
//     scipOracle = std::make_shared<SCIPOracle>("SCIPOracle(" + fileName + ")",  mixedIntegerSet, exactScipOracleStats);
//   }
//   else
// #endif /* IPO_WITH_EXACT_SCIP */
//   {
//     scipOracle = std::make_shared<SCIPOracle>("SCIPOracle(" + fileName + ")", mixedIntegerSet);
//   }
//   std::shared_ptr<StatisticsOracle> scipOracleStats = std::make_shared<StatisticsOracle>(scipOracle);
//   std::shared_ptr<CacheOracle> cacheOracle = std::make_shared<CacheOracle>(scipOracleStats, CacheOracle::CACHE_AND_SEARCH);
//   std::shared_ptr<StatisticsOracle> cacheOracleStats = std::make_shared<StatisticsOracle>(cacheOracle);
//   std::shared_ptr<OracleBase> oracle = cacheOracleStats;
// 
//   std::vector<AffineHullHandler*> affineHullHandlers;
//   DebugAffineHullHandler debugAffineHull(std::cout);
//   if (affineHullDebug)
//     affineHullHandlers.push_back(&debugAffineHull);
//   StatisticsAffineHullHandler statsAffineHull;
//   if (affineHullStats)
//     affineHullHandlers.push_back(&statsAffineHull);
// 
//   InnerDescription inner;
//   AffineOuterDescription outer;
//   affineHull(oracle, inner, outer, affineHullHandlers, cacheOracle->heuristicLevel(), cacheOracle->heuristicLevel());
//   std::cout << "Dimension: " << (long(inner.points.size() + inner.rays.size()) - 1) << "\n" << std::flush;
// 
//   if (printEquations)
//   {
//     std::cout << "\nEquations:\n";
//     for (std::size_t i = 0; i < outer.size(); ++i)
//     {
//       oracle->space().printLinearConstraint(std::cout, outer[i]);
//       std::cout << "\n";
//     }
//     std::cout << std::flush;
//   }
// 
//   if (affineHullStats)
//   {
//     std::cout << "\n";
//     std::cout << "Algorithm statistics for affine-hull computation:\n";
//     std::cout << "\n";
//     std::cout << "Overall time: " << statsAffineHull.timeAll() << "  =  main loop time: " << statsAffineHull.timeMainLoop()
//       << "  +  verification time: " << statsAffineHull.timeVerification() << "\n";
//     std::cout << "Approximate directions: " << statsAffineHull.numDirectionApproximateSolves() << " in " <<
//       statsAffineHull.timeApproximateDirections() << " seconds.\n";
//     std::cout << "Exact directions: " << statsAffineHull.numDirectionExactSolves() << " in " <<
//       statsAffineHull.timeExactDirections() << " seconds.\n";
//     std::cout << "Factorizations: " << statsAffineHull.numFactorizations() << " in " << statsAffineHull.timeFactorizations()
//       << " seconds.\n";
//     std::cout << "Oracle queries: " << statsAffineHull.numOracleQueries() << " in " << statsAffineHull.timeOracles()
//       << " seconds.\n";
//     std::cout << "\n";
//     std::cout << "Oracle statistics:\n";
//     std::cout << "\n";
//     for (std::shared_ptr<OracleBase> o = oracle; o != NULL; o = o->nextOracle())
//     {
//       std::size_t h = o->heuristicLevel();
//       std::cout << "Oracle " << h << ": " << o->name() << "\n";
//       std::shared_ptr<StatisticsOracle> s = std::dynamic_pointer_cast<StatisticsOracle>(o);
//       std::cout << "  #calls:   " << s->numCalls() << "\n";
//       std::cout << "  #success: " << s->numSuccess() << "\n";
//       std::cout << "  time:     " << s->time() << "\n";
//     }
//     std::cout << std::endl;
//   }
// 
// #ifdef IPO_WITH_EXACT_SCIP
//   if (exactUse)
//     exactScipOracleStats->reset();
// #endif /* IPO_WITH_EXACT_SCIP */
//   scipOracleStats->reset();
//   cacheOracleStats->reset();
// 
//   std::vector<FacetSeparationHandler*> facetSeparationHandlers;
//   DebugFacetSeparationHandler debugSeparation(std::cout, true, true);
//   if (separationDebug)
//     facetSeparationHandlers.push_back(&debugSeparation);
//   StatisticsFacetSeparationHandler statsSeparation;
//   if (separationStats)
//     facetSeparationHandlers.push_back(&statsSeparation);
// 
//   SoPlex spx;
//   spx.setIntParam(SoPlex::SOLVEMODE, SoPlex::SOLVEMODE_RATIONAL);
//   spx.setIntParam(SoPlex::SYNCMODE, SoPlex::SYNCMODE_AUTO);
//   spx.setRealParam(SoPlex::FEASTOL, 0.0);
//   spx.setBoolParam(SoPlex::RATREC, true);
//   spx.setBoolParam(SoPlex::RATFAC, true);
//   spx.setIntParam(SoPlex::OBJSENSE, SoPlex::OBJSENSE_MAXIMIZE);
//   spx.setIntParam(SoPlex::VERBOSITY, SoPlex::VERBOSITY_ERROR);
// 
//   std::shared_ptr<MixedIntegerLinearSet> mis = scipOracle->mixedIntegerLinearSet();
//   LPColSetRational cols(mis->numVariables());
//   DSVectorRational zero;
//   DVectorRational denseOriginalObjective(mis->numVariables());
//   vectorToDense(originalObjective, denseOriginalObjective);
//   for (std::size_t v = 0; v < oracle->space().dimension(); ++v)
//   {
//     cols.add(denseOriginalObjective[v], mis->lowerBound(v), zero, mis->upperBound(v));
//   }
//   spx.addColsRational(cols);
//   std::vector<LinearConstraint> rowConstraints;
//   mis->getConstraints(rowConstraints, true, true);
//   addToLP(spx, rowConstraints);
//   addToLP(spx, outer);
// 
//   //spx.writeFileRational("init.lp");
// 
//   DVectorRational solution(mis->numVariables());
//   std::default_random_engine generator(0);
//   for (std::size_t i = 0; i < numIterations; ++i)
//   {
//     if (i > 0)
//     {
//       std::cout << "Using random objective. " << std::flush;
//       std::normal_distribution<double> distribution;
//       DVectorReal randomVector(oracle->space().dimension());
//       double norm = 0;
//       for (std::size_t c = 0; c < oracle->space().dimension(); ++c)
//       {
//         double x = distribution(generator);
//         randomVector[c] = x;
//         norm += x*x;
//       }
//       if (norm > 0)
//       {
//         norm = std::sqrt(norm);
//         for (std::size_t v = 0; v < oracle->space().dimension(); ++v)
//         {
//           double x = randomVector[v] / norm;
// //           std::cerr << "Obj#" << v << " = " << x << std::endl;
//           spx.changeObjRational(v, Rational(x));
//         }
//       }
//     }
// 
//     std::cout << "Solving relaxation LP. " << std::flush;
// 
//     while (true)
//     {
//       SPxSolver::Status status = spx.solve();
//       if (status == SPxSolver::OPTIMAL)
//       {
//         std::cout << "\nSeparating LP optimum..." << std::flush;
//         spx.getPrimalRational(solution);
//         ipo::Vector point = denseToVector(solution, false);
//         oracle->space().printVector(std::cout, point);
//         std::cout << "\n";
//         InnerDescription certificate;
//         LinearConstraint constraint;
//         if (separatePoint(oracle, point, inner, facetSeparationHandlers, constraint, &certificate))
//         {
//           scaleIntegral(constraint);
//           manhattanNormImproveInequality(mis->numVariables(), constraint, outer);
// 
//           std::cout << "\n\n separated with facet or equation ";
//           oracle->space().printLinearConstraint(std::cout, constraint);
//           std::cout << std::endl;
// 
//           addToLP(spx, constraint);
//         }
//         else
//         {
//           std::cout << " feasible." << std::endl;
//           break;
//         }
//       }
//       else if (status == SPxSolver::UNBOUNDED)
//       {
//         std::cout << "Separating LP extreme ray..." << std::flush;
//         spx.getPrimalRayRational(solution);
//         ipo::Vector ray = denseToVector(solution, false);
//         InnerDescription certificate;
//         LinearConstraint constraint;
//         if (separateRay(oracle, ray, inner, facetSeparationHandlers, constraint, &certificate))
//         {
//           scaleIntegral(constraint);
//           manhattanNormImproveInequality(mis->numVariables(), constraint, outer);
// 
//           std::cout << "\n\n separated with facet or equation ";
//           oracle->space().printLinearConstraint(std::cout, constraint);
//           std::cout << std::endl;
// 
//           addToLP(spx, constraint);
//         }
//         else
//         {
//           std::cout << " feasible." << std::endl;
//           break;
//         }
//       }
//       else if (status == SPxSolver::INFEASIBLE)
//       {
//         std::cout << "LP infeasible." << std::endl;
//         spx.writeFileRational("infeasible.lp");
//         numIterations = i;
//         break;
//       }
//       else
//       {
//         std::cout << "LP has invalid status: " << status << std::endl;
//         break;
//       }
//     }
//   }
// 
//   if (separationStats)
//   {
//     std::cout << "\n";
//     std::cout << "Algorithm statistics (without affine-hull computation):\n";
//     std::cout << "\n";
//     std::cout << "Overall time: " << statsSeparation.timeAll() << "\n";
//     std::cout << "Approximate LPs: " << statsSeparation.numApproximateLPs() << " in " << statsSeparation.timeApproximateLPs()
//       << " seconds.\n";
//     std::cout << "Exact LPs: " << statsSeparation.numExactLPs() << " in " << statsSeparation.timeExactLPs()
//       << " seconds.\n";
//     std::cout << "Oracle queries: " << statsSeparation.numOracleQueries() << " in " << statsSeparation.timeOracles()
//       << " seconds.\n";
//     std::cout << "\n";
//     std::cout << "Oracle statistics (without affine hull computation):\n";
//     std::cout << "\n";
//     for (std::shared_ptr<OracleBase> o = oracle; o != NULL; o = o->nextOracle())
//     {
//       std::size_t h = o->heuristicLevel();
//       std::cout << "Oracle " << h << ": " << o->name() << "\n";
//       std::shared_ptr<StatisticsOracle> s = std::dynamic_pointer_cast<StatisticsOracle>(o);
//       std::cout << "  #calls:   " << s->numCalls() << "\n";
//       std::cout << "  #success: " << s->numSuccess() << "\n";
//       std::cout << "  time:     " << s->time() << "\n";
//     }
//     std::cout << std::endl;
//   }
// 
//   soplex::Rational::freeListMem();
//   ipo::Space::freeStaticMem();
//   ipo::Vector::freeStaticMem();

  return 0;
}
