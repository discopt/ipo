#include <gtest/gtest.h>
#include <ipo/oracles_scip.hpp>

#include <scip/scipdefplugins.h>

TEST(SCIP, Unbounded)
{
   SCIP* scip;
   SCIPcreate(&scip);
   SCIPincludeDefaultPlugins(scip);
   SCIPcreateProbBasic(scip, "unbounded");
   SCIP_VAR* x;
   SCIPcreateVarBasic(scip, &x, "x", 0.0, SCIPinfinity(scip), -1.0, SCIP_VARTYPE_INTEGER);
   SCIPaddVar(scip, x);

   SCIPreleaseVar(scip, &x);

   auto solver = new ipo::SCIPSolver(scip);
   auto oracle = solver->getOptimizationOracle();

   double obj[] = {1.0};
   ipo::OptimizationOracle::Query query;
   ipo::OptimizationOracle::Result result;
   oracle->maximize(obj, query, result);
   ASSERT_EQ(result.objectiveValues.size(), 1);
   ASSERT_FALSE(std::isfinite(result.objectiveValues[0]));
   ASSERT_FALSE(std::isfinite(result.dualBound));
}

TEST(SCIP, Infeasible)
{
   SCIP* scip;
   SCIPcreate(&scip);
   SCIPincludeDefaultPlugins(scip);
   SCIPcreateProbBasic(scip, "infeasible");
   SCIP_VAR* x;
   SCIPcreateVarBasic(scip, &x, "x", 0.0, SCIPinfinity(scip), -1.0, SCIP_VARTYPE_INTEGER);
   SCIPaddVar(scip, x);
   SCIP_VAR* y;
   SCIPcreateVarBasic(scip, &y, "y", 0.0, SCIPinfinity(scip), -1.0, SCIP_VARTYPE_INTEGER);
   SCIPaddVar(scip, y);
   SCIP_CONS* cons;
   SCIPcreateConsBasicLinear(scip, &cons, "cons", 0, nullptr, nullptr, 1.0, 1.0);
   SCIPaddCoefLinear(scip, cons, x, 2.0);
   SCIPaddCoefLinear(scip, cons, y, 2.0);
   SCIPaddCons(scip, cons);
   SCIPreleaseVar(scip, &x);
   SCIPreleaseVar(scip, &y);
   SCIPreleaseCons(scip, &cons);

   auto solver = new ipo::SCIPSolver(scip);
   auto oracle = solver->getOptimizationOracle();

   double obj[] = { 1.0, 2.0 };
   ipo::OptimizationOracle::Query query;
   ipo::OptimizationOracle::Result result;
   oracle->maximize(obj, query, result);
   ASSERT_EQ(result.objectiveValues.size(), 0);
   ASSERT_FALSE(std::isfinite(result.dualBound));
}


TEST(SCIP, Separate)
{
   SCIP* scip;
   SCIPcreate(&scip);
   SCIPincludeDefaultPlugins(scip);
   SCIPcreateProbBasic(scip, "cube");
   SCIP_VAR* x;
   SCIPcreateVarBasic(scip, &x, "x", -1.0, 1.0, 0.0, SCIP_VARTYPE_INTEGER);
   SCIPaddVar(scip, x);
   SCIP_VAR* y;
   SCIPcreateVarBasic(scip, &y, "y", -1.0, 1.0, 0.0, SCIP_VARTYPE_INTEGER);
   SCIPaddVar(scip, y);
   SCIP_CONS* cons;
   SCIPcreateConsBasicLinear(scip, &cons, "cons", 0, nullptr, nullptr, -1.0, 1.0);
   SCIPaddCoefLinear(scip, cons, x, 1.0);
   SCIPaddCoefLinear(scip, cons, y, 1.0);
   SCIPaddCons(scip, cons);
   SCIPreleaseVar(scip, &x);
   SCIPreleaseVar(scip, &y);
   SCIPreleaseCons(scip, &cons);

   auto solver = new ipo::SCIPSolver(scip);
   auto oracle = solver->getSeparationOracle();

   double vector[] = { 1.0, 1.0 };
   ipo::SeparationOracle::Query query;
   ipo::SeparationOracle::Result result;
   oracle->separate(vector, true, query, result);

   ASSERT_FALSE(result.nonzeroCoordinates.empty());
   
   result.reset();
   vector[0] = 0.5;
   vector[1] = 0.5;
   oracle->separate(vector, true, query, result);
   ASSERT_TRUE(result.nonzeroCoordinates.empty());
   oracle->separate(vector, false, query, result);
   ASSERT_FALSE(result.nonzeroCoordinates.empty());
}

