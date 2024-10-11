// #define IPO_DEBUG // Comment out to add debug information.

#include <ipo/oracles_scip.hpp>

#include "inverse.hpp"

int main(int argc, char** argv)
{
#if defined(IPO_DOUBLE) && defined(IPO_RATIONAL)
  bool exact = false;
#endif /* IPO_DOUBLE && IPO_RATIONAL */
  std::string fileName;
  for (int a = 1; a < argc; ++a)
  {
    const std::string arg = argv[a];
    if (arg == "-h")
    {
      inverse::printUsage(argv[0]);
      return EXIT_SUCCESS;
    }
#if defined(IPO_DOUBLE) && defined(IPO_RATIONAL)
    else if (arg == "-x")
      exact = true;
#endif /* IPO_DOUBLE && IPO_RATIONAL */
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

  auto scip = std::make_shared<ipo::SCIPSolver>(fileName);
#if defined(IPO_DOUBLE) && defined(IPO_RATIONAL)
  if (exact)
#endif /* IPO_DOUBLE && IPO_RATIONAL */
  {
#if defined(IPO_RATIONAL) && defined(IPO_RATIONAL_MIP_GUROBI)
    inverse::run<ipo::rational>(scip->getOptimizationOracle<ipo::rational>());
#endif /* IPO_RATIONAL && IPO_RATIONAL_MIP_GUROBI */
  }
#if defined(IPO_DOUBLE) && defined(IPO_RATIONAL)
  else
#endif /* IPO_DOUBLE && IPO_RATIONAL */
  {
#if defined(IPO_DOUBLE)
    inverse::run<double>(scip->getOptimizationOracle<double>());
#endif /* IPO_DOUBLE */
  }

  return 0;
}

