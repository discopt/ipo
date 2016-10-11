
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
#include "ipo/facets.h"
#include "ipo/scip_oracles.h"
#include "ipo/cache_oracle.h"
#include "ipo/statistics_oracle.h"

using namespace ipo;
using namespace soplex;

int main(int argc, char** argv)
{
  SCIP* scip = NULL;
  SCIP_CALL_EXC(SCIPcreate(&scip));
  SCIP_CALL_EXC(SCIPincludeDefaultPlugins(scip));
  SCIP_CALL_EXC(SCIPsetIntParam(scip, "display/verblevel", 0));
  SCIP_CALL_EXC(SCIPreadProb(scip, argv[1], NULL));
  SCIP_CALL_EXC(SCIPtransformProb(scip));

  std::shared_ptr<SCIPOracle> scipOracle = std::make_shared<SCIPOracle>("SCIPOracle(" + std::string(argv[1]) + ")", scip);
  std::shared_ptr<StatisticsOracle> scipOracleStats = std::make_shared<StatisticsOracle>(scipOracle);

  std::shared_ptr<CacheOracle> cacheOracle = std::make_shared<CacheOracle>(scipOracleStats);
  std::shared_ptr<StatisticsOracle> cacheOracleStats = std::make_shared<StatisticsOracle>(cacheOracle);

  std::shared_ptr<OracleBase> oracle = cacheOracleStats;

  SCIP_CALL_EXC(SCIPfree(&scip));

  std::vector<AffineHullHandler*> affineHullHandlers;
  InnerDescription inner;
  AffineOuterDescription outer;
  affineHull(oracle, inner, outer, affineHullHandlers, 2, 0);
  std::cout << "Dimension: " << (long(inner.points.size() + inner.rays.size()) - 1) << std::endl;

  std::vector<FacetSeparationHandler*> facetSeparationHandlers;
  
  SoPlex spx;
  spx.setIntParam(SoPlex::SOLVEMODE, SoPlex::SOLVEMODE_RATIONAL);
  spx.setIntParam(SoPlex::SYNCMODE, SoPlex::SYNCMODE_AUTO);
  spx.setRealParam(SoPlex::FEASTOL, 0.0);
  spx.setBoolParam(SoPlex::RATREC, true);
  spx.setBoolParam(SoPlex::RATFAC, true);
  spx.setIntParam(SoPlex::OBJSENSE, SoPlex::OBJSENSE_MAXIMIZE);
  spx.setIntParam(SoPlex::VERBOSITY, SoPlex::VERBOSITY_ERROR);

  const MixedIntegerSet& mis = *scipOracle->mixedIntegerSet();
  LPColSetRational cols(mis.numVariables());
  DSVectorRational zero;
  for (std::size_t v = 0; v < oracle->space().dimension(); ++v)
  {
    cols.add(Rational(1), mis.variable(v).lowerBound, zero, mis.variable(v).upperBound);
  }
  spx.addColsRational(cols);
  addToLP(spx, mis.rowConstraints());

  DVectorRational solution(mis.numVariables());
  while (true)
  {
    std::cout << "Solving relaxation LP with " << spx.numRowsRational() << " rows." << std::endl;

    SPxSolver::Status status = spx.solve();
    if (status == SPxSolver::OPTIMAL)
    {
      spx.getPrimalRational(solution);
      ipo::Vector point = denseToVector(solution, false);
      InnerDescription certificate;
      LinearConstraint constraint;
      if (separatePoint(oracle, point, inner, facetSeparationHandlers, constraint, &certificate))
      {
        scaleIntegral(constraint);
        
        std::cout << "Separated point with ";
        oracle->space().printLinearConstraint(std::cout, constraint);
        std::cout << std::endl;
        
        addToLP(spx, constraint);
      }
      else
        break;
      
    }
    else
      assert(false);
  }
  
  std::cout << "Cut loop done." << std::endl;

  return 0;
}