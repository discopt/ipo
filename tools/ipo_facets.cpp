#include <iostream>
#include <cmath>
#include <sstream>
#include <map>

#ifdef WITH_SCIP
#include <scip/scip.h>
#include <scip/scipdefplugins.h>
#include <scip/cons_linear.h>
#include <ipo/scip_exception.hpp>
#include <ipo/scip_oracles.h>
#endif

#include "ipo/affine_hull.h"
#include "ipo/facets.h"
#include "ipo/min_norm_2d.h"

using namespace ipo;

void printUsageAndExit(const char* program)
{
  std::cerr << "Usage: " << program << " [OPTIONS] MIP-MODEL\n";
  std::cerr << "Options:\n";
  std::cerr << std::flush;
  exit(EXIT_FAILURE);
}

int main(int argc, char** argv)
{
#ifdef WITH_SCIP

  {
    bool printCertificate = false;
    int fileArg = -1;
    for (int a = 1; a < argc; ++a)
    {
      if (fileArg < 0)
        fileArg = a;
      else
        printUsageAndExit(argv[0]);
    }
    if (fileArg < 0)
      printUsageAndExit(argv[0]);

    std::cout << "Computing facets of the mixed-integer hull of \"" << argv[fileArg] << "\".\n" << std::flush;

    SCIP* scip = NULL;
    SCIP_CALL_EXC(SCIPcreate(&scip));
    SCIP_CALL_EXC(SCIPincludeDefaultPlugins(scip));
    SCIP_CALL_EXC(SCIPsetIntParam(scip, "display/verblevel", 0));
    SCIP_CALL_EXC(SCIPreadProb(scip, argv[fileArg], NULL));
    SCIP_CALL_EXC(SCIPtransformProb(scip));

    Space  mipSpace;
    MixedIntegerProgram mip(mipSpace, scip);
    SCIP_CALL_EXC(SCIPfree(&scip));

    std::string name = argv[fileArg];
    SCIPOptimizationOracle* scipOracle = new SCIPOptimizationOracle(name, mip);
    MixedIntegerProgramCorrectorOracle* correctorOracle = new MixedIntegerProgramCorrectorOracle(name + "-corr", mip,
        scipOracle);
    ExactSCIPOptimizationOracle* exactOracle = new ExactSCIPOptimizationOracle(name + "-exact",
        "/home/xammy/software/exact-scip/scip-3.0.0-ex/bin/scip", mip, correctorOracle, 3600.0);
    OptimizationOracleBase* oracle = correctorOracle;
    std::size_t n = oracle->numVariables();

    /// Compute affine hull.

    UniqueRationalVectors points(n);
    UniqueRationalVectors rays(n);
    soplex::LPRowSetRational equations;
    AffineHull::ProgressOutput hullOutput;
    AffineHull::Result hull;
    hull.run(points, rays, equations, oracle, hullOutput, AffineHull::REDUNDANT_EQUATIONS_REMOVE);

    std::cout << "Dimension: " << hull.dimension() << std::endl;
    std::cout << "\n";

    manhattanNormImproveEquations(n, equations);

    std::cout << "\nList of equations:\n";
    oracle->printRows(std::cout, equations);
    std::cout << std::endl;

    /// Initialize cut LP.

    soplex::SoPlex spx;
    spx.setIntParam(soplex::SoPlex::SOLVEMODE, soplex::SoPlex::SOLVEMODE_RATIONAL);
    spx.setIntParam(soplex::SoPlex::SYNCMODE, soplex::SoPlex::SYNCMODE_AUTO);
    spx.setRealParam(soplex::SoPlex::FEASTOL, 0.0);
    spx.setBoolParam(soplex::SoPlex::RATREC, true);
    spx.setBoolParam(soplex::SoPlex::RATFAC, true);
    spx.setIntParam(soplex::SoPlex::OBJSENSE, soplex::SoPlex::OBJSENSE_MAXIMIZE);
    spx.setIntParam(soplex::SoPlex::VERBOSITY, soplex::SoPlex::VERBOSITY_ERROR);
    soplex::LPColSetRational cols;
    soplex::DSVectorRational vector;
    for (std::size_t v = 0; v < n; ++v)
      cols.add(mip.columns().maxObj(v), -soplex::infinity, vector, soplex::infinity);
    spx.addColsRational(cols);
    spx.addRowsRational(equations);

    /// Add inequalities from MIP.

    soplex::LPRowSetRational inequalities;
    mip.getConstraints(inequalities, true, false);
    spx.addRowsRational(inequalities);
    for (std::size_t v = 0; v < n; ++v)
      spx.changeBoundsRational(v, mip.columns().lower(v), mip.columns().upper(v));

    Separation::ProgressOutput separateOutput(2);
    Separation::Result separate(points, rays, hull.spanningDirections(), hull.spanningDirections(), hull.basicColumns(), oracle);

    soplex::DVectorRational denseSolution;
    soplex::DSVectorRational sparseSolution;
    denseSolution.reDim(n);
    while (true)
    {
      soplex::SPxSolver::Status status = spx.solve();
      if (status == soplex::SPxSolver::UNBOUNDED)
      {
        spx.getPrimalRayRational(denseSolution);
        sparseSolution.clear();
        sparseSolution = denseSolution;

        std::cout << "Relaxation LP is unbounded with ray ";
        oracle->printVector(std::cout, &sparseSolution);
        std::cout << std::endl;

        separate.separateRay(&sparseSolution, separateOutput);
        if (separate.violation() <= 0)
          break;
      }
      else if (status == soplex::SPxSolver::OPTIMAL)
      {
        spx.getPrimalRational(denseSolution);
        sparseSolution.clear();
        sparseSolution = denseSolution;

        std::cout << "Relaxation LP is bounded with optimum ";
        oracle->printVector(std::cout, &sparseSolution);
        std::cout << " of value " << spx.objValueRational() << "." << std::endl;

        separate.separatePoint(&sparseSolution, separateOutput);
        if (separate.violation() <= 0)
          break;
      }
      else
      {
        std::stringstream ss;
        ss << "Cut loop LP could not be solved to optimality. Status is " << status << ".";
        throw std::runtime_error(ss.str());
      }

      /// Obtain inequality and certificate.

      soplex::LPRowRational inequality;
      Separation::Certificate certificate;
      separate.inequality(inequality);
      separate.certificate(certificate);
      spx.addRowRational(inequality);
      if (separate.separatedFacet())
        std::cout << "Found a new facet: ";
      else if (separate.separatedEquation())
        std::cout << "Found a new equation: ";
      else
        std::cout << "Found a new face (neither facet nor equation!): ";

      manhattanNormImproveInequality(n, inequality, equations);

      oracle->printRow(std::cout, inequality);
      std::cout << ", certified by " << certificate.pointIndices.size() << " points and "
          << certificate.directionIndices.size() << " rays.\n" << std::endl;
      if (printCertificate)
      {
        for (std::size_t i = 0; i < certificate.pointIndices.size(); ++i)
        {
          std::cout << "  Certifying point: ";
          oracle->printVector(std::cout, points[certificate.pointIndices[i]]);
          std::cout << "\n";
        }
        for (std::size_t i = 0; i < certificate.directionIndices.size(); ++i)
        {
          std::cout << "  Certifying ray: ";
          oracle->printVector(std::cout, rays[certificate.directionIndices[i]]);
          std::cout << "\n";
        }
        std::cout << std::endl;
      }
    }

    separateOutput.printStatistics();

    if (exactOracle != NULL)
      delete exactOracle;
    delete correctorOracle;
    delete scipOracle;
  }

  soplex::Rational::freeListMem();

#endif

  return EXIT_SUCCESS;
}

