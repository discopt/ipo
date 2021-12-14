#include <iostream>
#include <sstream>

#include <ipo/oracles_scip.hpp>
#include <ipo/oracles_polar.hpp>
#include <ipo/projection.hpp>
#include <ipo/dominant.hpp>
#include <ipo/submissive.hpp>
#include <ipo/affine_hull.hpp>
#include <ipo/lp.hpp>

int printUsage(const std::string& program)
{
  std::cout << program << " [OPTIONS] FILE TASK...\n";
  std::cout << "Performs different computations on a polyhedron defined by FILE.\n";
  std::cout << "General options:\n";
  std::cout << " -h       Show this help and exit.\n";
  std::cout << " -t TIME  Abort computations after TIME seconds.\n";
  std::cout << "Oracle/polyhedron options:\n";
#if defined(IPO_RATIONAL)
  std::cout << " -x       Use exact arithmetic oracles instead of double precision.\n";
#endif /* IPO_RATIONAL */
  std::cout << " -p REGEX Project the polyhedron on all variables matching REGEX.\n";
  std::cout << " -d       Consider the dominant polyhedron.\n";
  std::cout << " -s       Consider the submissive polyhedron.\n";
  std::cout << "Tasks:\n";
  std::cout << " dimension       Output the dimension\n";
  std::cout << " equations       Output a complete system of valid equations.\n";
  std::cout << " interior        Output a point in the relative interior.\n";
  std::cout << " instance-facets Outputs facets encountered in a cutting plane method for the instance's objective.\n";
  std::cout << std::flush;

  return EXIT_FAILURE;
}

template <typename Number>
void run(std::shared_ptr<ipo::SCIPSolver> scip, std::shared_ptr<ipo::OptimizationOracle<Number>> baseOracle,
  std::shared_ptr<ipo::SeparationOracle<Number>> sepa, bool gmp, double timeLimit, const std::string& projectionRegex,
  bool useDominant, bool useSubmissive, bool outputDimension, bool outputEquations, bool outputInterior,
  bool outputInstanceFacets)
{
  std::shared_ptr<ipo::OptimizationOracle<Number>> projectionOracle;
  std::shared_ptr<ipo::Projection<Number>> projectionMap;
  if (projectionRegex.empty())
    projectionOracle = baseOracle;
  else
  {
    projectionMap = std::make_shared<ipo::Projection<Number>>(baseOracle->space(), projectionRegex);
    projectionOracle = std::make_shared<ipo::ProjectionOptimizationOracle<Number>>(baseOracle, projectionMap);
  }
  std::shared_ptr<ipo::OptimizationOracle<Number>> dominantOracle;
  if (useDominant)
    dominantOracle = std::make_shared<ipo::DominantOptimizationOracle<Number>>(projectionOracle);
  else
    dominantOracle = projectionOracle;

  std::shared_ptr<ipo::OptimizationOracle<Number>> submissiveOracle;
  if (useSubmissive)
    submissiveOracle = std::make_shared<ipo::SubmissiveOptimizationOracle<Number>>(dominantOracle);
  else
    submissiveOracle = dominantOracle;

  auto poly = std::make_shared<ipo::Polyhedron<Number>>(dominantOracle);

  std::cerr << "Initialized oracle with ambient dimension " << poly->space()->dimension() << std::endl;
  
  bool needAffineHull = outputDimension || outputEquations || outputInterior || outputInstanceFacets;

  ipo::SeparationResponse<Number> sepaResponse;
  if (needAffineHull)
    sepaResponse = sepa->getInitial();
  
  ipo::AffineHull<Number> affineHull;
  if (needAffineHull)
  {
    // Extract known equations from a separation oracle.
    std::vector<ipo::Constraint<Number>> knownEquations;
    for (const auto& cons : sepaResponse.constraints)
    {
      if (cons.type() == ipo::ConstraintType::EQUATION)
        knownEquations.push_back(cons);
    }
    std::vector<ipo::Constraint<Number>> projectedEquations;
    if (projectionRegex.empty())
      projectedEquations = knownEquations;
    else
      projectedEquations = projectionEquations(projectionMap, knownEquations);    

    std::cerr << "Starting affine hull computation.\n" << std::flush;
    ipo::AffineHullQuery affQuery;
    affQuery.timeLimit = timeLimit;
    affineHull = ipo::affineHull(poly, affQuery, knownEquations);
    if (outputDimension)
      std::cout << "Dimension: " << affineHull.dimension << std::endl;
    if (outputEquations)
    {
      for (auto& equation : affineHull.equations)
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
  
  if (outputInstanceFacets)
  {
    ipo::LP<Number> lp;
    lp.setSense(ipo::LPSense::MAXIMIZE);
    std::size_t scipVar = 0;
    for (std::size_t c = 0; c < poly->space()->dimension(); ++c)
    {
      const std::string& name = poly->space()->variable(c);
      while (baseOracle->space()->variable(scipVar) != name)
      {
        ++scipVar;
        assert(scipVar < baseOracle->space()->dimension());
      }
      lp.addColumn(lp.minusInfinity(), lp.plusInfinity(), Number(scip->instanceObjective()[scipVar]), name);
    }
    lp.update();
    
    std::vector<int> nonzeroColumns;
    std::vector<Number> nonzeroCoefficients;
    for (const ipo::Constraint<Number>& cons : sepaResponse.constraints)
    {
      if (cons.vector().size() == 1)
      {
        auto& coefficient = *cons.vector().begin();
        if (coefficient.second > 0)
        {
          if (cons.type() == ipo::ConstraintType::LESS_OR_EQUAL)
            lp.changeUpper(coefficient.first, cons.rhs() / coefficient.second);
          else if (cons.type() == ipo::ConstraintType::GREATER_OR_EQUAL)
            lp.changeLower(coefficient.first, cons.lhs() / coefficient.second);
          else
            lp.changeBounds(coefficient.first, cons.lhs() / coefficient.second, cons.rhs() / coefficient.second);
        }
        else
        {
          if (cons.type() == ipo::ConstraintType::LESS_OR_EQUAL)
            lp.changeLower(coefficient.first, cons.rhs() / coefficient.second);
          else if (cons.type() == ipo::ConstraintType::GREATER_OR_EQUAL)
            lp.changeUpper(coefficient.first, cons.lhs() / coefficient.second);
          else
            lp.changeBounds(coefficient.first, cons.rhs() / coefficient.second, cons.lhs() / coefficient.second); 
        }
      }
      else
      {
        Number lhs = cons.hasLhs() ? cons.lhs() : lp.minusInfinity();
        Number rhs = cons.hasRhs() ? cons.rhs() : lp.plusInfinity();
        nonzeroColumns.clear();
        nonzeroCoefficients.clear();
        for (const auto& iter : cons.vector())
        {
          nonzeroColumns.push_back(iter.first);
          nonzeroCoefficients.push_back(iter.second);
        }
        lp.addRow(lhs, nonzeroColumns.size(), &nonzeroColumns[0], &nonzeroCoefficients[0], rhs);
      }
    }

    for (std::size_t v = 0; v < poly->space()->dimension(); ++v)
      std::cout << poly->space()->variable(v) << std::endl;

    auto polarOracle = std::make_shared<ipo::PolarSeparationOracle<Number>>(poly);
    polarOracle->setAffineHull(affineHull);

    while (true)
    {
      auto status = lp.solve();
      std::cout << "Main LP status: " << status << "." << std::endl;
      if (status == ipo::LPStatus::OPTIMAL)
      {
        assert(lp.hasPrimalSolution());
        std::vector<Number> solution = lp.getPrimalSolution();
        for (std::size_t c = 0; c < poly->space()->dimension(); ++c)
        {
          std::cout << poly->space()->variable(c) << " = " << solution[c] << std::endl;
        }

        sepaResponse = polarOracle->separate(&solution[0], true);
        for (const auto& cons : sepaResponse.constraints)
          poly->space()->printConstraint(std::cout, cons);
      }

      break;
    }

//     lp.write("test.lp");
  }

}


class Bar
{
public:
  Bar()
  {
    std::cout << "Foo" << std::endl;
  }
};

int main(int argc, char** argv)
{
#if defined(IPO_DOUBLE) && defined(IPO_RATIONAL)
  bool exact = false;
#endif /* IPO_DOUBLE && IPO_RATIONAL */
  double timeLimit = std::numeric_limits<double>::infinity();
  std::string fileName;
  std::string projectionRegex = "";
  bool useDominant = false;
  bool useSubmissive = false;
  bool outputDimension = false;
  bool outputEquations = false;
  bool outputInterior = false;
  bool outputInstanceFacets = false;
  for (int a = 1; a < argc; ++a)
  {
    const std::string arg = argv[a];
    if (arg == "-h")
    {
      printUsage(argv[0]);
      return EXIT_SUCCESS;
    }
#if defined(IPO_DOUBLE) && defined(IPO_RATIONAL)
    else if (arg == "-x")
      exact = true;
#endif /* IPO_DOUBLE && IPO_RATIONAL */
    else if (arg == "-t" && a+1 < argc)
    {
      std::stringstream ss(argv[a+1]);
      ss >> timeLimit;
      ++a;
    }
    else if (arg == "-p" && a+1 < argc)
    {
      projectionRegex = argv[a+1];
      ++a;
    }
    else if (arg == "-d")
      useDominant = true;
    else if (arg == "-s")
      useSubmissive = true;
    else if (arg == "dimension")
      outputDimension = true;
    else if (arg == "equations")
      outputEquations = true;
    else if (arg == "interior")
      outputInterior = true;
    else if (arg == "instance-facets")
      outputInstanceFacets = true;
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
#if defined(IPO_DOUBLE) && defined(IPO_RATIONAL)
  if (exact)
#endif /* IPO_DOUBLE && IPO_RATIONAL */
  {
#if defined(IPO_RATIONAL)
    run<ipo::rational>(scip, scip->getOptimizationOracle<ipo::rational>(), scip->getSeparationOracle<ipo::rational>(), true,
      timeLimit, projectionRegex, useDominant, useSubmissive, outputDimension, outputEquations, outputInterior, outputInstanceFacets);
#endif /* IPO_RATIONAL */
  }
#if defined(IPO_DOUBLE) && defined(IPO_RATIONAL)
  else
#endif /* IPO_DOUBLE && IPO_RATIONAL */
  {
#if defined(IPO_DOUBLE)
    run<double>(scip, scip->getOptimizationOracle<double>(), scip->getSeparationOracle<double>(), false,
      timeLimit, projectionRegex, useDominant, useSubmissive, outputDimension, outputEquations, outputInterior, outputInstanceFacets);
#endif /* IPO_DOUBLE */
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
