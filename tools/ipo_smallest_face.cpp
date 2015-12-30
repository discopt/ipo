#include <iostream>
#include <cmath>
#include <sstream>
#include <map>

//#include <scip/scip.h>
//#include <scip/scipdefplugins.h>
//#include <scip/cons_linear.h>
//
//#include "ipo/scip_exception.hpp"
//#include "ipo/scip_oracles.h"
//#include "ipo/wrapper_oracle.h"
#include "ipo/affine_hull.h"
#include "ipo/smallest_face.h"

#include "common.h"

using namespace ipo;

void printUsageAndExit(const char* program)
{
  std::cerr << "Usage: " << program << " [OPTIONS] ORACLE\n";
  std::cerr << "Options:\n";
  std::cerr << "  -p POINT Specifies a point with a comma-separated list of VARNAME=VALUE, e.g., -p \"x#1=1,y#2=3\".\n";
  std::cerr << "  -r RAY Specifies a ray with a comma-separated list of VARNAME=VALUE, e.g., -p \"x#1=1,y#2=3\".\n";
  printOracleArguments(std::cerr);
  std::cerr << "\n";
  std::cerr << "Computes the dimension of the smallest face containing the given points and rays.\n";
  std::cerr << "For a single point x, the dimension is equal to 0 if and only if x is a vertex.\n";
  std::cerr << "For a two vertices x,y, the dimension is equal to 1 if and only if x and y are adjacent.\n";
  std::cerr << std::flush;
  exit(EXIT_FAILURE);
}

void parseElement(soplex::DVectorRational& element, const char* program,
    const std::map<std::string, std::size_t>& varMap, const std::string& input)
{
  std::size_t first = 0;
  std::size_t beyond = 0;
  std::size_t length = 0;
  while (beyond != std::string::npos)
  {
    beyond = input.find('=', first);
    if (beyond == std::string::npos)
    {
      if (first == input.size())
        break;
      std::cerr << "Error while parsing \"" << input << "\": " << input.substr(first, beyond) << " is invalid.\n\n";
      printUsageAndExit(program);
    }
    length = beyond - first;
    std::map<std::string, std::size_t>::const_iterator iter = varMap.find(input.substr(first, length));
    if (iter == varMap.end())
    {
      std::cerr << "Error while parsing \"" << input << "\": Variable " << input.substr(first, length)
          << " not found.\n\n";
      std::cerr << "List of all variables:\n";
      for (iter = varMap.begin(); iter != varMap.end(); ++iter)
        std::cerr << "  " << iter->first << "\n";
      std::cerr << "\n";
      printUsageAndExit(program);
    }
    first = beyond + 1;
    beyond = input.find(',', first);
    length = (beyond != std::string::npos) ? (beyond - first) : std::string::npos;
    soplex::Rational x;
    if (!soplex::readStringRational(input.substr(first, length).c_str(), x))
    {
      std::cerr << "Error while parsing \"" << input << "\": Value " << input.substr(first, length)
          << " is invalid.\n\n";
      printUsageAndExit(program);
    }
    element[iter->second] = x;
    first = beyond + 1;
  }
}

int main(int argc, char** argv)
{
  {
    std::vector<std::string> targetPoints, targetRays;
    int arg = 1;
    while (arg < argc)
    {
      if (argv[arg] == std::string("-p") && arg + 1 < argc)
      {
        targetPoints.push_back(argv[arg + 1]);
        arg += 2;
        continue;
      }
      else if (argv[arg] == std::string("-r") && arg + 1 < argc)
      {
        targetRays.push_back(argv[arg + 1]);
        arg += 2;
        continue;
      }
      else
        break;
    }
    MIPOracleInfo mipOracleInfo;
    OptimizationOracleBase* oracle;
    if (!parseOracleArguments(argc - arg, &argv[arg], oracle, mipOracleInfo))
      printUsageAndExit(argv[0]);

    if (targetPoints.empty())
    {
      std::cerr << "Error: No point specified.\n\n";
      printUsageAndExit(argv[0]);
    }

    std::size_t n = oracle->numVariables();

    /// Parse input points and rays.

    std::map<std::string, std::size_t> varMap;
    for (std::size_t v = 0; v < n; ++v)
      varMap[oracle->variableName(v)] = v;
    soplex::DVectorRational element;
    element.reDim(n);
    soplex::DVectorRational denseTarget;
    denseTarget.reDim(n, true);
    for (std::size_t i = 0; i < targetPoints.size(); ++i)
    {
      element.clear();
      parseElement(element, argv[0], varMap, targetPoints[i]);
      denseTarget += element;
    }
    denseTarget *= soplex::Rational(1) / soplex::Rational(int(targetPoints.size()));
    for (std::size_t i = 0; i < targetRays.size(); ++i)
    {
      element.clear();
      parseElement(element, argv[0], varMap, targetRays[i]);
      denseTarget += element;
    }
    soplex::DSVectorRational sparseTarget;
    sparseTarget = denseTarget;

    std::cout << "Computing dimension of smallest face of \"" << oracle->name() << "\" containing " << std::flush;
    oracle->printVector(std::cout, &sparseTarget);
    std::cout << ".\n" << std::flush;

    /// Compute affine hull.

    std::cout << "\nStarting affine hull algorithm for main polyhedron living in a " << oracle->numVariables()
        << "-dimensional space.\n" << std::flush;
    UniqueRationalVectors points(n);
    UniqueRationalVectors rays(n);
    soplex::LPRowSetRational equations;
    AffineHull::ProgressOutput hullOutput;
    AffineHull::Result hull;
    hull.run(points, rays, equations, oracle, hullOutput, AffineHull::REDUNDANT_EQUATIONS_REMOVE); // TODO: exactOracle!
    std::cout << "Dimension of main polyhedron: " << hull.dimension() << "\n" << std::flush;

    oracle->printRows(std::cout, equations);

    /// Compute smallest face.

    std::cout << "\nStarting to compute smallest face.\n" << std::flush;
    SmallestFace::ProgressOutput smallestFaceOutput;
    SmallestFace::Result smallestFace(points, rays, oracle); // TODO: exactOracle!
    smallestFace.run(&sparseTarget, smallestFaceOutput);
    std::cout << "Dimension of smallest face: " << smallestFace.dimension() << std::endl;
    Point maximizingObjective;
    smallestFace.getMaximizingObjective(maximizingObjective);
    std::cout << "It is the maximum face of ";
    oracle->printVector(std::cout, &maximizingObjective);
    std::cout << std::endl;

    delete oracle;
  }

  soplex::Rational::freeListMem();

  return EXIT_SUCCESS;
}

