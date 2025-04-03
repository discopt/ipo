// #define IPO_DEBUG // Comment out to add debug information.

#include <ipo/oracles_scip.hpp>

#include "inverse.hpp"

int main(int argc, char** argv)
{
  bool exact = false;
  std::string fileName;
  for (int a = 1; a < argc; ++a)
  {
    const std::string arg = argv[a];
    if (arg == "-h")
    {
      inverse::printUsage(argv[0]);
      return EXIT_SUCCESS;
    }
#if defined(IPO_WITH_RATIONAL)
    else if (arg == "-x")
      exact = true;
#endif /* IPO_WITH_RATIONAL */
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

#if defined(IPO_WITH_RATIONAL)
  if (exact)
  {
    auto instance = inverse::readInverseProblem<ipo::SCIPSolver, ipo::rational>(fileName);

    inverse::solve(instance);
  }
#endif /* IPO_WITH_RATIONAL */

  if (!exact)
  {
    auto instance = inverse::readInverseProblem<ipo::SCIPSolver, double>(fileName);

    inverse::solve(instance);
  }

  return 0;
}

