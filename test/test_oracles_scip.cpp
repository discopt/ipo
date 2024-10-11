#include <gtest/gtest.h>
#include <ipo/oracles_scip.hpp>

#include <scip/scipdefplugins.h>

TEST(SCIP, UnboundedDouble)
{
  // x >= 0
  // y >= 0
  // y == 2x + 1

  SCIP* scip;
  SCIPcreate(&scip);
  SCIPincludeDefaultPlugins(scip);
  SCIPcreateProbBasic(scip, "unbounded");
  SCIP_VAR* x;
  SCIP_VAR* y;
  SCIPcreateVarBasic(scip, &x, "x", 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_INTEGER);
  SCIPcreateVarBasic(scip, &y, "y", 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_INTEGER);
  SCIPaddVar(scip, x);
  SCIPaddVar(scip, y);

  SCIP_CONS* cons;
  SCIPcreateConsBasicLinear(scip, &cons, "xy", 0, 0, 0, 1.0, 1.0);
  SCIPaddCoefLinear(scip, cons, x, 2.0);
  SCIPaddCoefLinear(scip, cons, y, -1.0);
  SCIPaddCons(scip, cons);

  SCIPreleaseCons(scip, &cons);
  SCIPreleaseVar(scip, &x);
  SCIPreleaseVar(scip, &y);

  auto solver = std::make_shared<ipo::SCIPSolver>(std::move(scip));
  auto oracle = solver->getOptimizationOracle<double>();

  double obj[] = {1.0, 1.0};
  ipo::OptimizationOracle<double>::Query query;
  std::cout << "oracle->maximize()" << std::endl;
  auto response = oracle->maximize(obj, query);
  ASSERT_EQ(response.outcome, ipo::OptimizationOutcome::UNBOUNDED);
  ASSERT_FALSE(response.hasDualBound);
  ASSERT_FALSE(response.rays.empty());
  ASSERT_EQ(response.rays[0].vector->size(), 2);
  ASSERT_EQ(response.rays[0].vector->begin()->first, 0);
  ASSERT_EQ((response.rays[0].vector->begin()+1)->first, 1);
  ASSERT_NEAR(response.rays[0].vector->begin()->second / (response.rays[0].vector->begin()+1)->second,
    0.5, 1.0e-9);
}

TEST(SCIP, InfeasibleDouble)
{
  // x >= 0
  // y >= 0
  // 2x + 2y == 1

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

  auto solver = std::make_shared<ipo::SCIPSolver>(std::move(scip));
  auto oracle = solver->getOptimizationOracle<double>();

  double obj[] = { 1.0, 2.0 };
  ipo::OptimizationOracle<double>::Query query;
  auto response = oracle->maximize(obj, query);
  ASSERT_EQ(response.outcome, ipo::OptimizationOutcome::INFEASIBLE);
  ASSERT_FALSE(response.hasDualBound);
}

TEST(SCIP, SeparateDouble)
{
  // x in [-1.1,1.2]
  // y in [-1.3,1.4]
  // -1.5 <= x + y <= 1.5

  SCIP* scip;
  SCIPcreate(&scip);
  SCIPincludeDefaultPlugins(scip);
  SCIPcreateProbBasic(scip, "separate");
  SCIP_VAR* x;
  SCIPcreateVarBasic(scip, &x, "x", -1.1, 1.2, 0.0, SCIP_VARTYPE_INTEGER);
  SCIPaddVar(scip, x);
  SCIP_VAR* y;
  SCIPcreateVarBasic(scip, &y, "y", -1.3, 1.4, 0.0, SCIP_VARTYPE_INTEGER);
  SCIPaddVar(scip, y);
  SCIP_CONS* cons;
  SCIPcreateConsBasicLinear(scip, &cons, "cons", 0, nullptr, nullptr, -1.5, 1.5);
  SCIPaddCoefLinear(scip, cons, x, 1.0);
  SCIPaddCoefLinear(scip, cons, y, 1.0);
  SCIPaddCons(scip, cons);
  SCIPreleaseVar(scip, &x);
  SCIPreleaseVar(scip, &y);
  SCIPreleaseCons(scip, &cons);

  auto solver = std::make_shared<ipo::SCIPSolver>(std::move(scip));
  auto oracle = solver->getSeparationOracle<double>();

  double vector[] = { 1.0, 1.0 };
  ipo::SeparationQuery query;
  auto response = oracle->separate(vector, true, query);
  ASSERT_FALSE(response.constraints.empty());
  ASSERT_EQ(response.constraints[0].type(), ipo::ConstraintType::LESS_OR_EQUAL);
  ASSERT_NEAR(response.constraints[0].rhs(), 1.5, 1.0e-9);

  vector[0] = 0.5;
  vector[1] = 0.5;
  response = oracle->separate(vector, true, query);
  ASSERT_TRUE(response.constraints.empty());

  response = oracle->separate(vector, false, query);
  ASSERT_FALSE(response.constraints.empty());
  ASSERT_EQ(response.constraints[0].type(), ipo::ConstraintType::LESS_OR_EQUAL);
  ASSERT_TRUE(fabs(response.constraints[0].rhs() - 1.2) < 1.0e-9
    || fabs(response.constraints[0].rhs() - 1.4) < 1.0e-9
    || fabs(response.constraints[0].rhs() - 1.5) < 1.0e-9
  );

  vector[0] = 2.0;
  vector[1] = -2.0;
  response = oracle->separate(vector, true, query);
  ASSERT_FALSE(response.constraints.empty());
  ASSERT_EQ(response.constraints[0].type(), ipo::ConstraintType::LESS_OR_EQUAL);
  ASSERT_NEAR(response.constraints[0].rhs(), 1.2, 1.0e-9);
}

#if defined(IPO_RATIONAL_MIP_SCIP)

TEST(SCIP, UnboundedRational)
{
  // x >= 0
  // y >= 0
  // y >= 2x + 1

  SCIP* scip;
  SCIPcreate(&scip);
  SCIPincludeDefaultPlugins(scip);
  SCIPcreateProbBasic(scip, "unbounded");
  SCIP_VAR* x;
  SCIP_VAR* y;
  SCIPcreateVarBasic(scip, &x, "x", 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_INTEGER);
  SCIPcreateVarBasic(scip, &y, "y", 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_INTEGER);
  SCIPaddVar(scip, x);
  SCIPaddVar(scip, y);

  SCIP_CONS* cons;
  SCIPcreateConsBasicLinear(scip, &cons, "xy", 0, 0, 0, 1.0, 1.0);
  SCIPaddCoefLinear(scip, cons, x, 2.0);
  SCIPaddCoefLinear(scip, cons, y, -1.0);
  SCIPaddCons(scip, cons);

  SCIPreleaseCons(scip, &cons);
  SCIPreleaseVar(scip, &x);
  SCIPreleaseVar(scip, &y);

  auto solver = std::make_shared<ipo::SCIPSolver>(std::move(scip));
  auto oracle = solver->getOptimizationOracle<ipo::rational>();

  ipo::rational obj[] = {1.0, 1.0};
  ipo::OptimizationOracle<ipo::rational>::Query query;
  auto response = oracle->maximize(obj, query);
  ASSERT_EQ(response.outcome, ipo::OptimizationOutcome::UNBOUNDED);
  ASSERT_FALSE(response.hasDualBound);
  ASSERT_FALSE(response.rays.empty());
  ASSERT_EQ(response.rays[0].vector->size(), 2);
  ASSERT_EQ(response.rays[0].vector->begin()->first, 0);
  ASSERT_EQ((response.rays[0].vector->begin()+1)->first, 1);
  ASSERT_EQ(response.rays[0].vector->begin()->second / (response.rays[0].vector->begin()+1)->second,
    ipo::rational(1, 2));
}

TEST(SCIP, InfeasibleRational)
{
  // x >= 0
  // y >= 0
  // 2x + 2y == 1

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

  auto solver = std::make_shared<ipo::SCIPSolver>(std::move(scip));
  auto oracle = solver->getOptimizationOracle<ipo::rational>();

  ipo::rational obj[] = { 1, 2 };
  ipo::OptimizationOracle<ipo::rational>::Query query;
  auto response = oracle->maximize(obj, query);
  ASSERT_EQ(response.outcome, ipo::OptimizationOutcome::INFEASIBLE);
  ASSERT_FALSE(response.hasDualBound);
}

TEST(SCIP, SeparateRational)
{
  // x in [-1,2]
  // y in [-4,8]
  // -16 <= x + y <= 16

  SCIP* scip;
  SCIPcreate(&scip);
  SCIPincludeDefaultPlugins(scip);
  SCIPcreateProbBasic(scip, "separate");
  SCIP_VAR* x;
  SCIPcreateVarBasic(scip, &x, "x", -1, 2, 0.0, SCIP_VARTYPE_INTEGER);
  SCIPaddVar(scip, x);
  SCIP_VAR* y;
  SCIPcreateVarBasic(scip, &y, "y", -4, 8, 0.0, SCIP_VARTYPE_INTEGER);
  SCIPaddVar(scip, y);
  SCIP_CONS* cons;
  SCIPcreateConsBasicLinear(scip, &cons, "cons", 0, nullptr, nullptr, -16, 16);
  SCIPaddCoefLinear(scip, cons, x, 1.0);
  SCIPaddCoefLinear(scip, cons, y, 1.0);
  SCIPaddCons(scip, cons);
  SCIPreleaseVar(scip, &x);
  SCIPreleaseVar(scip, &y);
  SCIPreleaseCons(scip, &cons);

  auto solver = std::make_shared<ipo::SCIPSolver>(std::move(scip));
  auto oracle = solver->getSeparationOracle<ipo::rational>();

  ipo::rational vector[] = { 20, 20 };
  ipo::SeparationQuery query;
  auto response = oracle->separate(vector, true, query);
  ASSERT_FALSE(response.constraints.empty());
  ASSERT_EQ(response.constraints[0].type(), ipo::ConstraintType::LESS_OR_EQUAL);
  ASSERT_EQ(response.constraints[0].rhs(), ipo::rational(2));

  vector[0] = 0.5;
  vector[1] = 0.5;
  response = oracle->separate(vector, true, query);
  ASSERT_TRUE(response.constraints.empty());

  response = oracle->separate(vector, false, query);
  ASSERT_FALSE(response.constraints.empty());
  ASSERT_EQ(response.constraints[0].type(), ipo::ConstraintType::LESS_OR_EQUAL);
  ASSERT_EQ(response.constraints[0].rhs(), ipo::rational(2));

  vector[0] = 3.0;
  vector[1] = -3.0;
  response = oracle->separate(vector, true, query);
  ASSERT_FALSE(response.constraints.empty());
  ASSERT_EQ(response.constraints[0].type(), ipo::ConstraintType::LESS_OR_EQUAL);
  ASSERT_EQ(response.constraints[0].rhs(), ipo::rational(2));
}

#endif /* IPO_RATIONAL_MIP_SCIP */
