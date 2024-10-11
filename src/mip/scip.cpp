// #define IPO_DEBUG // Comment out to add debug information.

#include <ipo/oracles_scip.hpp>

#include "ipo.hpp"

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
      ipo::printUsage(argv[0]);
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
    ipo::run<ipo::rational>(scip->getOptimizationOracle<ipo::rational>(), scip->getSeparationOracle<ipo::rational>(),
      scip->instanceObjective(), true, timeLimit, randomSeed, projectionRegex, useDominant, useSubmissive,
      outputDimension, outputEquations, outputInterior, outputInstanceFacets, outputRandomFacets);
#endif /* IPO_RATIONAL && IPO_RATIONAL_MIP_SCIP */
  }
#if defined(IPO_DOUBLE) && defined(IPO_RATIONAL)
  else
#endif /* IPO_DOUBLE && IPO_RATIONAL */
  {
#if defined(IPO_DOUBLE)
    ipo::run<double>(scip->getOptimizationOracle<double>(), scip->getSeparationOracle<double>(),
      scip->instanceObjective(), false, timeLimit, randomSeed, projectionRegex, useDominant, useSubmissive,
      outputDimension, outputEquations, outputInterior, outputInstanceFacets, outputRandomFacets);
#endif /* IPO_DOUBLE */
  }

  return 0;
}
