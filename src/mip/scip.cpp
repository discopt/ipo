#define IPO_DEBUG // Comment out to add debug information.

#include <iostream>
#include <sstream>
#include <random>

#include <ipo/oracles_scip.hpp>
#include <ipo/oracles_polar.hpp>
#include <ipo/projection.hpp>
#include <ipo/dominant.hpp>
#include <ipo/submissive.hpp>
#include <ipo/affine_hull.hpp>
#include <ipo/lp.hpp>

template <typename Number>
void run(std::shared_ptr<ipo::SCIPSolver> scip, std::shared_ptr<ipo::OptimizationOracle<Number>> baseOracle,
  std::shared_ptr<ipo::SeparationOracle<Number>> sepa, bool gmp, double timeLimit, int randomSeed,
  const std::string& projectionRegex, bool useDominant, bool useSubmissive, bool outputDimension, bool outputEquations,
  bool outputInterior, bool outputInstanceFacets, int outputRandomFacets)
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

  auto poly = std::make_shared<ipo::Polyhedron<Number>>(submissiveOracle);

  std::cerr << "Initialized oracle with ambient dimension " << poly->space()-> dimension() << std::endl;

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

    // Project the equations.
    std::vector<ipo::Constraint<Number>> projectedEquations;
    if (useDominant || useSubmissive)
    {
    }
    else if (projectionRegex.empty())
      projectedEquations = knownEquations;
    else
      projectedEquations = projectionEquations(projectionMap, knownEquations);

    // Execute actual affine hull algorithm.

    std::cerr << "Starting affine hull computation.\n" << std::flush;
    ipo::AffineHullQuery affQuery;
    affQuery.timeLimit = timeLimit;
    affineHull = ipo::affineHull(poly, affQuery, projectedEquations);
    if (outputDimension)
      std::cout << "Dimension: " << affineHull.dimension << std::endl;
    if (outputEquations)
    {
      auto projectedEquationSet = ipo::ConstraintSet<Number>(projectionOracle->space()->dimension(),
        projectedEquations);
      for (auto& equation : affineHull.equations)
      {
        bool isKnown = projectedEquationSet.exists(equation);
        if (!isKnown)
          scaleIntegral(equation);
        std::cout << (isKnown ? "Known " : "New ") << "equation "
          << poly->space()->printConstraint(equation, true) << std::endl;
      }
    }
  }

  if (outputInterior)
  {
    std::cout << "Output of interior point is not implemented, yet." << std::endl;
  }

  std::default_random_engine generator;
  generator.seed(randomSeed);
  std::normal_distribution<double> normal_distribution;
  if (outputInstanceFacets || outputRandomFacets)
  {
#if defined(IPO_DOUBLE_LP) || defined(IPO_RATIONAL_LP)

    ipo::LP<Number> lp;
    lp.setSense(ipo::LPSense::MAXIMIZE);
    std::size_t scipVar = 0;
    bool first = true;
    for (std::size_t c = 0; c < poly->space()->dimension(); ++c)
    {
      const std::string& name = poly->space()->variable(c);
      while (baseOracle->space()->variable(scipVar) != name)
      {
        ++scipVar;
        assert(scipVar < baseOracle->space()->dimension());
      }
      lp.addColumn(lp.minusInfinity(), lp.plusInfinity(), Number(scip->instanceObjective()[scipVar]), name);
      if (scip->instanceObjective()[scipVar])
      {
        std::cerr << (first ? "Objective: " : " + ") << scip->instanceObjective()[scipVar] << "*"
          << baseOracle->space()->variable(scipVar);
        first = false;
      }
    }
    std::cerr << std::endl;
    lp.update();

    std::vector<int> nonzeroColumns;
    std::vector<Number> nonzeroCoefficients;
    if (!useSubmissive && !useDominant)
    {
      std::vector<ipo::Constraint<Number>> constraints;
      if (!projectionRegex.empty())
      {
        constraints = ipo::projectionCompact(projectionMap, sepaResponse.constraints, false);
      }
      else
        constraints = sepaResponse.constraints;

      for (const ipo::Constraint<Number>& cons : constraints)
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
    }

    // Add new equations from affine hull computation.
    for (const auto& equation : affineHull.equations)
    {
      Number rhs = equation.rhs();
      nonzeroColumns.clear();
      nonzeroCoefficients.clear();
      for (const auto& iter : equation.vector())
      {
        nonzeroColumns.push_back(iter.first);
        nonzeroCoefficients.push_back(iter.second);
      }
      lp.addRow(rhs, nonzeroColumns.size(), &nonzeroColumns[0], &nonzeroCoefficients[0], rhs);
    }

    auto polarOracle = std::make_shared<ipo::PolarSeparationOracle<Number>>(poly);
    polarOracle->setAffineHull(affineHull);

    int randomIteration = outputInstanceFacets ? -1 : 0;
    while (randomIteration < outputRandomFacets)
    {
      if (randomIteration >= 0)
      {
        // Generate a vector uniformly at random from the sphere.
        std::vector<double> samples(poly->space()->dimension());
        double squaredLength = 0.0;
        while (squaredLength < 1.0e-12)
        {
          for (std::size_t c = 0; c < poly->space()->dimension(); ++c)
          {
            double x = normal_distribution(generator);
            samples[c] = x;
            squaredLength += x*x;
          }
        }
        double length = sqrt(squaredLength);
        for (std::size_t c = 0; c < poly->space()->dimension(); ++c)
        {
          lp.changeObjective(c, Number(samples[c] / length));
        }
      }

      auto status = lp.solve();
      if (status == ipo::LPStatus::OPTIMAL)
      {
        std::cerr << "\nThe LP optimum is " << ipo::convertNumber<double>(lp.getObjectiveValue()) << "="
          << lp.getObjectiveValue() << "." << std::endl;
        assert(lp.hasPrimalSolution());
        std::vector<Number> solution = lp.getPrimalSolution();
        for (std::size_t c = 0; c < poly->space()->dimension(); ++c)
          std::cerr << (c == 0 ? "Solution vector: " : ", ") << poly->space()->variable(c) << "=" << solution[c];
        std::cerr << std::endl;

        sepaResponse = polarOracle->separate(&solution[0], true);
        std::cerr << "Found " << sepaResponse.constraints.size() << " undominated constraints." << std::endl;
        for (ipo::Constraint<Number>& cons : sepaResponse.constraints)
        {
          scaleIntegral(cons);
          std::cout << "  ";
          poly->space()->printConstraint(std::cout, cons);
          std::cout << std::endl;

          nonzeroColumns.clear();
          nonzeroCoefficients.clear();
          for (const auto& iter : cons.vector())
          {
            nonzeroColumns.push_back(iter.first);
            nonzeroCoefficients.push_back(iter.second);
          }
          lp.addRow(cons.hasLhs() ? cons.lhs() : lp.minusInfinity(), nonzeroColumns.size(), &nonzeroColumns[0],
            &nonzeroCoefficients[0], cons.hasRhs() ? cons.rhs() : lp.plusInfinity());
        }
        if (sepaResponse.constraints.empty())
        {
          ++randomIteration;
        }
      }
      else if (status == ipo::LPStatus::UNBOUNDED)
      {
        std::cerr << "\nThe LP is unbounded." << std::endl;
        assert(lp.hasPrimalRay());
        std::vector<Number> ray = lp.getPrimalRay();
        for (std::size_t c = 0; c < poly->space()->dimension(); ++c)
          std::cerr << (c == 0 ? "Unbounded ray: " : ", ") << poly->space()->variable(c) << "=" << ray[c];
        std::cerr << std::endl;

        sepaResponse = polarOracle->separate(&ray[0], false);
        std::cerr << "Found " << sepaResponse.constraints.size() << " undominated constraints." << std::endl;
        for (auto& cons : sepaResponse.constraints)
        {
          scaleIntegral(cons);
          std::cout << "  ";
          poly->space()->printConstraint(std::cout, cons);
          std::cout << std::endl;

          nonzeroColumns.clear();
          nonzeroCoefficients.clear();
          for (const auto& iter : cons.vector())
          {
            nonzeroColumns.push_back(iter.first);
            nonzeroCoefficients.push_back(iter.second);
          }
          lp.addRow(cons.hasLhs() ? cons.lhs() : lp.minusInfinity(), nonzeroColumns.size(), &nonzeroColumns[0],
            &nonzeroCoefficients[0], cons.hasRhs() ? cons.rhs() : lp.plusInfinity());
        }
        if (sepaResponse.constraints.empty())
        {
          ++randomIteration;
        }
      }
      else
      {
        std::cerr << "\nMain LP has unhandled status " << status << "." << std::endl;
        assert(false);
        break;
      }
    }

#else /* of (IPO_DOUBLE_LP || IPO_RATIONAL_LP) */

    std::cerr << "Option instance-facets requires an LP solver such as SoPlex." << std::endl;

#endif /* IPO_DOUBLE_LP || IPO_RATIONAL_LP */
  }
}


int printUsage(const std::string& program)
{
  std::cout << program << " [OPTIONS] FILE TASK...\n";
  std::cout << "Performs different computations on a polyhedron defined by FILE.\n";
  std::cout << "General options:\n";
  std::cout << " -h       Show this help and exit.\n";
  std::cout << " -t TIME  Abort computations after TIME seconds.\n";
  std::cout << " -S SEED  Use SEED to initialize the random number generator.\n";
  std::cout << "Oracle/polyhedron options:\n";
#if defined(IPO_RATIONAL)
  std::cout << " -x       Use exact arithmetic oracles instead of double precision.\n";
#endif /* IPO_RATIONAL */
  std::cout << " -p REGEX Project the polyhedron on all variables matching REGEX.\n";
  std::cout << " -d       Consider the dominant polyhedron.\n";
  std::cout << " -s       Consider the submissive polyhedron.\n";
  std::cout << "Tasks:\n";
  std::cout << " dimension         Output the dimension\n";
  std::cout << " equations         Output a complete system of valid equations.\n";
  std::cout << " interior          Output a point in the relative interior.\n";
  std::cout << " instance-facets   Outputs facets of a cutting plane method for the instance's objective.\n";
  std::cout << " random-facets NUM Outputs facets of a cutting plane method for NUM random objectives.\n";
  std::cout << std::flush;

  return EXIT_FAILURE;
}

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
  int outputRandomFacets = 0;
  int randomSeed = 0;
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
    else if (arg == "-S" && a+1 < argc)
    {
      std::stringstream ss(argv[a+1]);
      ss >> randomSeed;
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
    else if (arg == "random-facets" && a+1 < argc)
    {
      std::stringstream ss(argv[a+1]);
      ss >> outputRandomFacets;
      ++a;
    }
    else if (fileName.empty())
      fileName = arg;
    else
    {
      std::cerr << "Invalid parameters. Interpreted <" << fileName << "> and <" << arg
        << "> as instance files." << std::endl;
      return EXIT_FAILURE;
    }
  }

  if (fileName.empty())
  {
    std::cerr << "Invalid parameters. No filename provided." << std::endl;
    return EXIT_FAILURE;
  }
  if (outputRandomFacets < 0)
  {
    std::cerr << "Invalid parameters. Negative value " << outputRandomFacets << " for number of random facets."
      << std::endl;
    return EXIT_FAILURE;
  }

  auto scip = std::make_shared<ipo::SCIPSolver>(fileName);
#if defined(IPO_DOUBLE) && defined(IPO_RATIONAL)
  if (exact)
#endif /* IPO_DOUBLE && IPO_RATIONAL */
  {
#if defined(IPO_RATIONAL) && defined(IPO_RATIONAL_MIP_SCIP)
    run<ipo::rational>(scip, scip->getOptimizationOracle<ipo::rational>(), scip->getSeparationOracle<ipo::rational>(),
      true, timeLimit, randomSeed, projectionRegex, useDominant, useSubmissive, outputDimension, outputEquations, outputInterior,
      outputInstanceFacets, outputRandomFacets);
#endif /* IPO_RATIONAL && IPO_RATIONAL_MIP_SCIP */
  }
#if defined(IPO_DOUBLE) && defined(IPO_RATIONAL)
  else
#endif /* IPO_DOUBLE && IPO_RATIONAL */
  {
#if defined(IPO_DOUBLE)
    run<double>(scip, scip->getOptimizationOracle<double>(), scip->getSeparationOracle<double>(), false,
      timeLimit, randomSeed, projectionRegex, useDominant, useSubmissive, outputDimension, outputEquations, outputInterior,
      outputInstanceFacets, outputRandomFacets);
#endif /* IPO_DOUBLE */
  }

  return 0;
}
