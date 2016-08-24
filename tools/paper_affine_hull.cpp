
#include <ipo/common.h>

#ifdef WITH_SCIP
#ifdef NDEBUG
  #undef NDEBUG
  #include <scip/scip.h>
  #include <scip/scipdefplugins.h>
  #include <scip/cons_linear.h>
  #include <ipo/scip_exception.hpp>
  #include <ipo/scip_oracles.h>
  #define NDEBUG
#else
  #include <scip/scip.h>
  #include <scip/scipdefplugins.h>
  #include <scip/cons_linear.h>
  #include <ipo/scip_exception.hpp>
  #include <ipo/scip_oracles.h>
#endif
#endif

#include "ipo/scip_exception.hpp"
#include "ipo/affine_hull.h"
#include "ipo/scip_oracles.h"
#include "ipo/cache_oracle.h"

using namespace ipo;

int main(int argc, char** argv)
{
  SCIP* scip = NULL;
  SCIP_CALL_EXC(SCIPcreate(&scip));
  SCIP_CALL_EXC(SCIPincludeDefaultPlugins(scip));
  SCIP_CALL_EXC(SCIPreadProb(scip, argv[1], NULL));
  SCIP_CALL_EXC(SCIPtransformProb(scip));

  std::shared_ptr<SCIPOracle> scipOracle = std::make_shared<SCIPOracle>(argv[1], scip);
  std::shared_ptr<CacheOracle> cacheOracle = std::make_shared<CacheOracle>("cached " + std::string(argv[1]), scipOracle);

  SCIP_CALL_EXC(SCIPfree(&scip));

  std::vector<AffineHullHandler*> handlers;
  DebugAffineHullHandler handler(std::cout);
  handlers.push_back(&handler);

  std::vector<Vector> points, rays;
  std::vector<LinearConstraint> givenEquations, equations;
  affineHull(handlers, cacheOracle, givenEquations, points, rays, equations, 1);

  return 0;
}