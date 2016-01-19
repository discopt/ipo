#include <iostream>
#include <cmath>
#include <sstream>
#include <map>

#ifdef WITH_SCIP
#include <ipo/scip_exception.hpp>
#include <scip/scip.h>
#include <scip/scipdefplugins.h>
#include <scip/cons_linear.h>
#include "ipo/scip_oracles.h"
#endif

#include "ipo/affine_hull.h"
#include "ipo/min_norm_2d.h"

using namespace ipo;

void printUsageAndExit(const char* program)
{
  std::cerr << "Usage: " << program << " [OPTIONS] MIP-MODEL\n";
  std::cerr << "Options:\n";
  std::cerr << "  -R   Consider relaxation instead of mixed-integer hull.\n";
  std::cerr << "  -O   Consider optimal face.\n";
  std::cerr << "  -C   Consider faces corresponding to inequality constraints.\n";
  std::cerr << "  -nc  Apply algorithm with cache turned off.\n";
  std::cerr << "  -a   Use approximate oracle only.\n";
  std::cerr << "  -d   Use approximate oracle and verify equations at the end (default).\n";
  std::cerr << "  -i   Use approximate oracle and verify equations immediately.\n";
  std::cerr << "  -x   Use exact oracle only.\n";
  std::cerr << "  -nie Don't improve the equations.\n";
  std::cerr << "  -nse Don't show equations.\n";
  std::cerr << std::flush;
  exit(EXIT_FAILURE);
}

int main(int argc, char** argv)
{
#ifdef WITH_SCIP

  int fileArg = -1;
  bool relaxation = false;
  bool optimalFace = false;
  bool constraintFaces = false;
  bool enableCache = true;
  bool exactNever = false;
  bool exactDelayed = true;
  bool exactImmediate = false;
  bool exactOnly = false;
  bool showEquations = true;
  bool improveEquations = true;
  for (int a = 1; a < argc; ++a)
  {
    if (argv[a] == std::string("-R"))
      relaxation = true;
    else if (argv[a] == std::string("-O"))
      optimalFace = true;
    else if (argv[a] == std::string("-C"))
      constraintFaces = true;
    else if (argv[a] == std::string("-nc"))
      enableCache = false;
    else if (argv[a] == std::string("-nie"))
      improveEquations = false;
    else if (argv[a] == std::string("-nse"))
      showEquations = false;
    else if (argv[a] == std::string("-a"))
    {
      exactDelayed = false;
      exactNever = true;
    }
    else if (argv[a] == std::string("-x"))
    {
      exactDelayed = false;
      exactOnly = true;
    }
    else if (argv[a] == std::string("-d"))
    {
      exactDelayed = true;
    }
    else if (argv[a] == std::string("-i"))
    {
      exactDelayed = false;
      exactImmediate = true;
    }
    else if (fileArg < 0)
      fileArg = a;
    else
      printUsageAndExit(argv[0]);
  }
  if (optimalFace && constraintFaces)
  {
    std::cerr << "Invalid arguments: Cannot pass -O and -Cat the same time!\n\n";
    printUsageAndExit(argv[0]);
  }
  if ((exactDelayed ? 1 : 0) + (exactImmediate ? 1 : 0) + (exactNever ? 1 : 0) + (exactOnly ? 1 : 0) >= 2)
  {
    std::cerr << "Invalid arguments: Cannot set multiple exact/approximate options!\n\n";
    printUsageAndExit(argv[0]);
  }

  if (fileArg < 0)
    printUsageAndExit(argv[0]);

//  std::cout << "P is the "
//      << (relaxation ? "LP relaxation" : "mixed-integer hull") << " of \""
//      << argv[fileArg] << "\".\n" << std::flush;

  SCIP* scip = NULL;
  SCIP_CALL_EXC(SCIPcreate(&scip));
  SCIP_CALL_EXC(SCIPincludeDefaultPlugins(scip));
  SCIP_CALL_EXC(SCIPsetIntParam(scip, "display/verblevel", 0));
  SCIP_CALL_EXC(SCIPreadProb(scip, argv[fileArg], NULL));

  if (relaxation)
  {
    int nvars = SCIPgetNOrigVars(scip);
    bool changed = true;
    while (changed)
    {
      changed = false;
      SCIP_VAR** vars = SCIPgetOrigVars(scip);
      unsigned int infeasible;
      for (int i = 0; i < nvars; ++i)
      {
        if (SCIPvarGetType(vars[i]) != SCIP_VARTYPE_CONTINUOUS)
        {
          SCIP_CALL_EXC(SCIPchgVarType(scip, vars[i], SCIP_VARTYPE_CONTINUOUS, &infeasible));
          changed = true;
        }
      }
    }
  }

  SCIP_CALL_EXC(SCIPtransformProb(scip));

  std::string name = argv[fileArg];
  MixedIntegerProgram mip(scip);
  SCIPOptimizationOracle* scipOracle = new SCIPOptimizationOracle(name, mip);
  MixedIntegerProgramCorrectorOracle* correctorOracle = new MixedIntegerProgramCorrectorOracle(name + "-corr", mip,
      scipOracle);
  ExactSCIPOptimizationOracle* exactOracle = new ExactSCIPOptimizationOracle(name + "-exact",
      "/home/xammy/software/exact-scip/scip-3.0.0-ex/bin/scip", mip, correctorOracle, 3600.0);
  FaceOptimizationOracleBase* oracle = correctorOracle;
  std::size_t n = oracle->numVariables();

  int options = (enableCache ? AffineHull::CACHE_USE : AffineHull::CACHE_SKIP) | AffineHull::REDUNDANT_EQUATIONS_REMOVE;
  if (exactDelayed)
    options |= AffineHull::ORACLE_DELAYED;
  if (exactImmediate)
    options |= AffineHull::ORACLE_IMMEDIATE;
  if (exactNever)
    options |= AffineHull::ORACLE_NEVER;
  if (exactOnly)
    options |= AffineHull::ORACLE_ONLY;

  if (optimalFace)
  {
    soplex::DVectorRational objective;
    getSCIPObjective(scip, objective, true);

    OptimizationResult result;
    oracle->maximize(result, objective);

    if (result.isInfeasible())
    {
      std::cout << "Main polyhedron is empty.\n" << std::flush;
    }
    else if (result.isUnbounded())
    {
      std::cout << "Problem is unbounded.\n" << std::flush;
    }
    else
    {
      soplex::DSVectorRational faceNormal;
      faceNormal = objective;
      Face* face = new Face(n, soplex::LPRowRational(result.bestValue, faceNormal, result.bestValue));

      oracle->setFace(face);

      std::cout << "Computing affine hull of optimal face in " << n << "-dimensional space." << std::endl;

      UniqueRationalVectors points(n);
      UniqueRationalVectors rays(n);
      soplex::LPRowSetRational equations;

      AffineHull::Result hull;
      AffineHull::ProgressOutput hullOutput;
      hull.run(points, rays, equations, oracle, hullOutput, options);

      std::cout << "Its dimension is " << hull.dimension() << "." << std::endl;

      oracle->setFace(NULL);
      delete face;
    }
  }
  else
  {
    std::cout << "Computing affine hull of main polyhedron in " << n << "-dimensional space." << std::endl;

    UniqueRationalVectors points(n);
    UniqueRationalVectors rays(n);
    soplex::LPRowSetRational equations;
    mip.getConstraints(equations, false, true);
    mip.getFixedVariableEquations(equations);

    AffineHull::Result hull;
    AffineHull::ProgressOutput mainOutput(2);
    hull.run(points, rays, equations, oracle, mainOutput, options);

    std::cout << "Its dimension is " << hull.dimension() << "." << std::endl;

    if (improveEquations)
      manhattanNormImproveEquations(n, equations);

    if (showEquations)
    {
      std::cout << "\nList of equations:\n";
      oracle->printRows(std::cout, equations);
      std::cout << std::endl;
    }

    if (constraintFaces && hull.dimension() >= 0)
    {
      std::cout << "Computing the dimensions of its constraint-faces." << std::endl;

      soplex::LPRowSetRational inequalities;
      std::vector<std::string> names;
      mip.getConstraints(inequalities, true, false, &names);

      FilteredUniqueRationalVectors filteredPoints(points);
      FilteredUniqueRationalVectors filteredRays(rays);
      soplex::LPRowSetRational faceEquations;

      std::map<int, std::size_t> dimensionStatistics;
      dimensionStatistics[hull.dimension() - 1] = 0;
      for (int i = 0; i < inequalities.num(); ++i)
      {
        /// Create face.

        Face* face = new Face(n,
            soplex::LPRowRational(inequalities.rhs(i), inequalities.rowVector(i), inequalities.rhs(i)));
        oracle->setFace(face);

        /// Filter points and rays.

        soplex::DVectorRational faceNormal(n);
        faceNormal.clear();
        faceNormal.assign(face->normal());

        for (std::size_t p = points.first(); p < points.size(); p = points.next(p))
        {
          soplex::Rational activity = *points[p] * faceNormal;
          filteredPoints.set(p, activity == face->rhs());
        }
        for (std::size_t r = rays.first(); r < rays.size(); r = rays.next(r))
        {
          soplex::Rational activity = *rays[r] * faceNormal;
          filteredRays.set(r, activity == 0);
        }

        std::cout << "  Computing the dimension of face " << names[i] << " (" << (i + 1) << "/" << inequalities.num()
            << ") using " << filteredPoints.count() << " known points, " << filteredRays.count() << " rays and "
            << (equations.num() + 1) << " equations." << std::endl;

        AffineHull::Result faceHull;
        AffineHull::ProgressOutput faceOutput;

        faceEquations = equations;
        faceEquations.add(face->rhs(), face->normal(), face->rhs());
        faceHull.run(filteredPoints, filteredRays, faceEquations, oracle, faceOutput, options);

        std::cout << "  The dimension of face " << names[i] << " (" << (i + 1) << "/" << inequalities.num() << ") is "
            << faceHull.dimension() << "." << std::endl;

        std::map<int, std::size_t>::iterator iter = dimensionStatistics.find(faceHull.dimension());
        if (iter == dimensionStatistics.end())
          dimensionStatistics[faceHull.dimension()] = 1;
        else
          iter->second++;

        oracle->setFace(NULL);
        delete face;
      }
      std::cout << "Finished computation of constraint faces." << std::endl;
      for (std::map<int, std::size_t>::const_iterator iter = dimensionStatistics.begin();
          iter != dimensionStatistics.end(); ++iter)
      {
        std::cout << "Number of faces of dimension " << std::setw(3) << iter->first;
        if (iter->first == hull.dimension())
          std::cout << " (equations)";
        else if (iter->first == hull.dimension() - 1)
          std::cout << "    (facets)";
        else if (iter->first == 1)
          std::cout << "     (edges)";
        else if (iter->first == 0)
          std::cout << "  (vertices)";
        else if (iter->first == -1)
          std::cout << "     (empty)";
        else
          std::cout << "            ";
        std::cout << ": " << iter->second << "\n";
      }
      std::cout << std::flush;
    }
  }

  delete exactOracle;
  delete correctorOracle;
  delete scipOracle;

  SCIP_CALL_EXC(SCIPfree(&scip));

  soplex::Rational::freeListMem();
#endif

  return EXIT_SUCCESS;
}

