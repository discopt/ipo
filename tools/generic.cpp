#include <iostream>
#include <cmath>
#include <sstream>
#include <map>

#include "scip_exception.hpp"
#include <scip/scip.h>
#include <scip/scipdefplugins.h>

#include "polycomb/scip_oracle.h"
#include "polycomb/scip_oracles.h"
#include "polycomb/shoot_unit.h"
#include "polycomb/shoot_random.h"
#include "polycomb/affine_hull.h"
#include "polycomb/facets.h"
#include "polycomb/reconstruct.h"

using namespace soplex;

int main(int argc, char** argv)
{
  if (argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " MIP-MODEL-FILE\n" << std::flush;
    return EXIT_FAILURE;
  }

  SCIP* scip = NULL;
  SCIP_CALL_EXC(SCIPcreate(&scip));
  SCIP_CALL_EXC(SCIPincludeDefaultPlugins(scip));
  SCIP_CALL_EXC(SCIPsetIntParam(scip, "display/verblevel", 0));
  SCIP_CALL_EXC(SCIPreadProb(scip, argv[1], NULL));
  SCIP_CALL_EXC(SCIPtransformProb(scip));

  polycomb::OldSCIPOPtimizationOracle* opt = new polycomb::OldSCIPOPtimizationOracle(
      "opt", scip, 1);

  polycomb::SCIPRelaxationSeparationOracle* sepaRelaxation =
      new polycomb::SCIPRelaxationSeparationOracle("sepa", scip);

  polycomb::BoxSeparationOracle* sepaBox = new polycomb::BoxSeparationOracle(
      "-1/2-box", opt, Rational(-1), Rational(2));

  std::size_t n = opt->numVariables();

  std::cout << "Oracles lives in " << n << "-dimensional space." << std::endl;

  soplex::LPRowSetRational equations;
  sepaRelaxation->populateRows(equations, false);

  std::cout << "Obtained " << equations.num()
      << " equations from initial description." << std::endl;

  polycomb::UniqueRationalVectors points(n);
  polycomb::UniqueRationalVectors rays(n);

  std::cout << "Shooting in unit directions" << std::flush;
  polycomb::shootUnitDirections(points, rays, opt);
  std::cout << "We know " << points.size() << " points and " << rays.size()
      << " rays." << std::endl;
//  points.printStatistics(std::cerr) << std::endl;

  LPRowSetRational newEquations;
  polycomb::VectorSubset spanningPoints, spanningRays;
  std::vector<bool> basicColumns;
  int dim = polycomb::computeAffineHull(opt, equations, points, rays,
      newEquations, 2, &spanningPoints, &spanningRays, &basicColumns);
  std::cout << "Mixed-integer hull has dimension " << dim << "." << std::endl;

  if (dim >= 0)
  {
    std::cout << "We know " << points.size() << " points and " << rays.size()
        << " rays." << std::endl;

    std::cout << "Given equations:\n";
    opt->printRows(std::cout, equations);
    polycomb::scaleVectorsIntegral(newEquations);
    std::cout << "New equations:\n";
    opt->printRows(std::cout, newEquations);

    polycomb::AugmentedSeparationOracle* sepaWithEquations =
        new polycomb::AugmentedSeparationOracle("box+eqns", sepaBox);
    sepaWithEquations->addRows(equations);
    sepaWithEquations->addRows(newEquations);

    DVectorRational originalObjective;
    polycomb::getSCIPObjective(scip, originalObjective);

    LPRowSetRational facets;
    std::vector<polycomb::FacetCertificate> certs;
    polycomb::computeFacets(opt, sepaWithEquations, points, rays,
        spanningPoints, spanningRays, basicColumns, originalObjective, facets,
        &certs);

    polycomb::scaleVectorsIntegral(facets);

    std::cout << "\nFacets:\n";
    for (int f = 0; f < facets.num(); ++f)
    {
      opt->printRow(std::cout, facets, f);
      std::cout << "\n  (certified by " << certs[f].pointIndices.size()
          << " points and " << certs[f].rayIndices.size() << " rays.)"
          << std::endl;
    }

    std::cout << "We know " << points.size() << " points, " << rays.size()
        << " rays and " << facets.num() << " facets." << std::endl;

    delete sepaWithEquations;
  }

  delete sepaBox;
  delete sepaRelaxation;
  delete opt;

  SCIP_CALL_EXC(SCIPfree(&scip));

  soplex::Rational::freeListMem();

  return EXIT_SUCCESS;
}

