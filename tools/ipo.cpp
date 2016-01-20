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
#include "ipo/smallest_face.h"
#include "ipo/min_norm_2d.h"

using namespace ipo;

#ifdef WITH_SCIP
#define ORACLE_DEFAULT_SCIP
#define ORACLE_DEFAULT "scip"
#elif WITH_EXACT_SCIP
#define ORACLE_DEFAULT_EXACT_SCIP
#define ORACLE_DEFAULT "exactscip"
#else
#define ORACLE_DEFAULT_NONE
#define ORACLE_DEFAULT "";
#endif

void printUsageAndExit(const char* program)
{
  ///           01        11        21        31        41        51        61        71        81        91        <--End
  std::cerr << "Usage: " << program << " [OPTIONS] INSTANCE...\n\n";
  std::cerr << "Options for defining the base polyhedron P:\n";
  std::cerr << " --heuristic HEUR             Use the oracle HEUR as a heuristic in addition to the oracle.\n";
  std::cerr << " --oracle ORACLE              Use the oracle ORACLE to define the polyhedron P.\n";
  std::cerr << "                              See section on oracles below.\n";
  std::cerr << " --relaxation                 In case of a MIP, set all variable to continuous.\n";
  std::cerr << " --projection VARS            Consider the orthogonal projection onto a subset of variables.\n";
  std::cerr << "                              VARS is either a comma-separated list of variables or a regular\n";
  std::cerr << "                              expression enclosed with ^ and $. In this case it defines the\n";
  std::cerr << "                              projection onto all variables matching it.\n";
  std::cerr << " --restrict-face INEQ         Restrict the oracle to the specified face. See --face below.\n";
  std::cerr << "\n";
  std::cerr << "Oracles that can be used for ORACLE and HEUR (see IPO's build options):\n";
#ifdef WITH_SCIP
  std::cerr << " scip                         Use the MIP solver SCIP from the SCIP Optimization Suite\n";
  std::cerr << "                              (default oracle).\n";
#endif
#ifdef WITH_EXACT_SCIP
  std::cerr << " exactscip                    Use the exact MIP solver SCIP-ex from the SCIP Optimization Suite";
#ifdef ORACLE_DEFAULt_EXACT_SCIP
  std::cerr << "\n                              (default oracle)";
#endif
  std::cerr << ".\n";
#endif
  std::cerr << " PROGRAM                      For each oracle call, execute PROGRAM and communicate via stdin and\n";
  std::cerr << "                              stdout. See IPO's python directory for example implementations.\n";
#ifdef ORACLE_DEFAULT_NONE
  std::cerr << "           IPO is built without a default oracle (e.g., SCIP). See build options!\n";
#endif
  std::cerr << "\n";
  std::cerr << "Options for tasks (execution order is independent of argument order):\n";
  std::cerr << " --ambient-dimension          Print the ambient dimension.\n";
  std::cerr << " --variables                  Print the names of all variables in a comma-separated list.\n";
  std::cerr << " --instance-objective         Print the maximization objective from the given (MIP) instance.\n";
  std::cerr << " --maximize                   Maximize objectives specified by --objective arguments.\n";
  std::cerr << " --minimize                   Minimize objectives specified by --objective arguments.\n";
  std::cerr << " --dimension                  Compute the dimension.\n";
  std::cerr << " --equations                  Compute valid equations.\n";
  std::cerr << " --facet                      Attempt to separate every given point or direction using a\n";
  std::cerr << "                              facet-defining inequality.\n";
  std::cerr << " --smallest-face              Computes, for every given point in P, an objective vector that induces\n";
  std::cerr << "                              the smallest face containing the point.\n";
  std::cerr << " --facets                     Generate facets that are useful when maximizing the objectives\n";
  std::cerr << "                              specified by --objective or --objectives and the randomly generated\n";
  std::cerr << "                              ones (see --random).\n";
  std::cerr << "\n";
  std::cerr << "Options for input data (ordering w.r.t. tasks does not matter):\n";
  std::cerr << " --face INEQ                  Compute dimension/equations also for face induced by INEQ.\n";
  std::cerr << "                              INEQ is an inequality, potentially preceeded by a name with a colon,\n";
  std::cerr << "                              e.g., \"myface: x#1 + 3y#5 <= 7\". Usable multiple times.\n";
  std::cerr << " --faces FILE                 Each line of the file is considered as a single --face argument. If it\n";
  std::cerr << "                              does not match the pattern, the line is silently ignored, allowing it\n";
  std::cerr << "                              to be in LP format.\n";
  std::cerr << " --objective OBJ              Use the given objective vector for optimization/facet generation.\n";
  std::cerr << "                              OBJ is a sparse vector, e.g., \"x#1 + 3y#5\". Usable multiple times.\n";
  std::cerr << " --objectives FILE            Each line of the file is considered as a single --objective argument.\n";
  std::cerr << "                              If it does not match the pattern, the line is silently ignored.\n";
  std::cerr << " --random n                   Sample n objective vectors (uniformly at random from the sphere) to be\n";
  std::cerr << "                              used for facet generation.\n";
  std::cerr
      << " --point POINT                Use the given point for smallest-face computation and facet-separation.\n";
  std::cerr
      << "                              POINT is a sparse vector, e.g., \"(x#1, 3y#5)\". Usable multiple times.\n";
  std::cerr << " --direction DIR              Use the given unbounded direction for facet-separation.\n";
  std::cerr << "                              DIR is a sparse vector, e.g., \"(x#1, 3y#5)\". Usable multiple times.\n";
  std::cerr << "\n";
  std::cerr << "Further options (ordering does not matter):\n";
  std::cerr << " --cache on|off               Cache points and directions produced by IPO to speed up further\n";
  std::cerr << "                              computations. On by default.\n";
  std::cerr << " --readable on|off            Improve readability of inequalities/equations produced by IPO.\n";
  std::cerr << "                              On by default.\n";
  std::cerr << " --certificates on|off        Output, for every computed facet, points and directions that span it.\n";
  std::cerr << "                              Off by default.\n";
  std::cerr << " --trust-points on|off        If off (default), it is checked that points are in P when computing\n";
  std::cerr << "                              the smallest containing faces.\n";
  std::cerr << " --instance-objective on|off  Consider an objective from a (MIP) instance as if it was given as\n";
  std::cerr << "                              a --objective argument. On by default.\n";
  std::cerr << " --progress stdout|stderr|off Prints progress output to stdout or stderr. Off by default.\n";
  std::cerr << " --debug stdout|stderr|off    Prints debug output to stdout or stderr. Off by default.\n";
  std::cerr << std::flush;
  exit(EXIT_FAILURE);
}

int main(int argc, char** argv)
{
  std::string argOracle = ORACLE_DEFAULT
  ;
  std::string argHeuristic = "";
  bool argRelaxation = false;
  std::string argProjectionDescription = "";
  std::string argRestrictFaceDescription = "";
  bool optionCache = true;
  bool optionReadable = true;
  bool optionCertificates = false;
  bool optionTrustPoints = false;
  bool optionInstanceObjective = true;

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

    MixedIntegerProgram mip(scip);
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
    Separation::Result separate(points, rays, hull.spanningPoints(), hull.spanningRays(), hull.basicColumns(), oracle);

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
      << certificate.rayIndices.size() << " rays.\n" << std::endl;
      if (printCertificate)
      {
        for (std::size_t i = 0; i < certificate.pointIndices.size(); ++i)
        {
          std::cout << "  Certifying point: ";
          oracle->printVector(std::cout, points[certificate.pointIndices[i]]);
          std::cout << "\n";
        }
        for (std::size_t i = 0; i < certificate.rayIndices.size(); ++i)
        {
          std::cout << "  Certifying ray: ";
          oracle->printVector(std::cout, rays[certificate.rayIndices[i]]);
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

