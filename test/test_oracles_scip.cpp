#include <gtest/gtest.h>
#include <ipo/oracles_scip.hpp>

#include <scip/scipdefplugins.h>

TEST(Oracles, SCIP)
{
  std::cout << "===== Oracles::SCIP::Double::Unbounded ===== " << std::endl;
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

    auto solver = std::make_shared<ipo::SCIPSolver>(scip);
    auto oracle = solver->getOptimizationOracleDouble();

    double obj[] = {1.0, 1.0};
    ipo::OptimizationOracle<double>::Query query;
    auto result = oracle->maximize(obj, query);
    ASSERT_TRUE(result.isUnbounded());
    ASSERT_TRUE(ipo::isPlusInfinity(result.dualBound));
    ASSERT_FALSE(result.rays.empty());
    ASSERT_EQ(result.rays[0].vector.size(), 2);
    ASSERT_EQ(result.rays[0].vector.begin()->first, 0);
    ASSERT_EQ((result.rays[0].vector.begin()+1)->first, 1);
    ASSERT_NEAR(result.rays[0].vector.begin()->second / (result.rays[0].vector.begin()+1)->second,
      0.5, 1.0e-9);
  }

  std::cout << "===== Oracles::SCIP::Double::Infeasible ===== " << std::endl;
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

    auto solver = std::make_shared<ipo::SCIPSolver>(scip);
    auto oracle = solver->getOptimizationOracleDouble();

    double obj[] = { 1.0, 2.0 };
    ipo::OptimizationOracle<double>::Query query;
    auto result = oracle->maximize(obj, query);
    ASSERT_FALSE(result.isUnbounded());
    ASSERT_FALSE(result.isFeasible());
    ASSERT_TRUE(ipo::isMinusInfinity(result.dualBound));
  }

  std::cout << "===== Oracles::SCIP::Double::Separate ===== " << std::endl;
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

    auto solver = std::make_shared<ipo::SCIPSolver>(scip);
    auto oracle = solver->getSeparationOracleDouble();

    double vector[] = { 1.0, 1.0 };
    ipo::SeparationOracle<double>::Query query;
    auto result = oracle->separate(vector, true, query);
    ASSERT_FALSE(result.constraints.empty());
    ASSERT_TRUE(ipo::isMinusInfinity(result.constraints[0].lhs()));
    ASSERT_NEAR(result.constraints[0].rhs(), 1.5, 1.0e-9);

    vector[0] = 0.5;
    vector[1] = 0.5;
    result = oracle->separate(vector, true, query);
    ASSERT_TRUE(result.constraints.empty());

    result = oracle->separate(vector, false, query);
    ASSERT_FALSE(result.constraints.empty());
    ASSERT_TRUE(ipo::isMinusInfinity(result.constraints[0].lhs()));
    ASSERT_TRUE(fabs(result.constraints[0].rhs() - 1.2) < 1.0e-9
      || fabs(result.constraints[0].rhs() - 1.4) < 1.0e-9
      || fabs(result.constraints[0].rhs() - 1.5) < 1.0e-9
    );

    vector[0] = 2.0;
    vector[1] = -2.0;
    result = oracle->separate(vector, true, query);
    ASSERT_FALSE(result.constraints.empty());
    ASSERT_TRUE(ipo::isMinusInfinity(result.constraints[0].lhs()));
    ASSERT_NEAR(result.constraints[0].rhs(), 1.2, 1.0e-9);
  }

#if defined(IPO_WITH_GMP) && defined(IPO_WITH_SOPLEX)

  std::cout << "===== Oracles::SCIP::Rational::Unbounded ===== " << std::endl;
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

    auto solver = std::make_shared<ipo::SCIPSolver>(scip);
    auto oracle = solver->getOptimizationOracleRational();

    ipo::rational obj[] = {1.0, 1.0};
    ipo::OptimizationOracle<ipo::rational>::Query query;
    auto result = oracle->maximize(obj, query);
    ASSERT_TRUE(result.isUnbounded());
    ASSERT_TRUE(ipo::isPlusInfinity(result.dualBound));
    ASSERT_FALSE(result.rays.empty());
    ASSERT_EQ(result.rays[0].vector.size(), 2);
    ASSERT_EQ(result.rays[0].vector.begin()->first, 0);
    ASSERT_EQ((result.rays[0].vector.begin()+1)->first, 1);
    ASSERT_EQ(result.rays[0].vector.begin()->second / (result.rays[0].vector.begin()+1)->second,
      ipo::rational(1, 2));
  }

  std::cout << "===== Oracles::SCIP::Rational::Infeasible ===== " << std::endl;
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

    auto solver = std::make_shared<ipo::SCIPSolver>(scip);
    auto oracle = solver->getOptimizationOracleRational();

    ipo::rational obj[] = { 1, 2 };
    ipo::OptimizationOracle<ipo::rational>::Query query;
    auto result = oracle->maximize(obj, query);
    ASSERT_FALSE(result.isUnbounded());
    ASSERT_FALSE(result.isFeasible());
    ASSERT_TRUE(ipo::isMinusInfinity(result.dualBound));
  }

  std::cout << "===== Oracles::SCIP::Rational::Separate ===== " << std::endl;
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

    auto solver = std::make_shared<ipo::SCIPSolver>(scip);
    auto oracle = solver->getSeparationOracleRational();

    ipo::rational vector[] = { 20, 20 };
    ipo::SeparationOracle<ipo::rational>::Query query;
    auto result = oracle->separate(vector, true, query);
    ASSERT_FALSE(result.constraints.empty());
    ASSERT_TRUE(ipo::isMinusInfinity(result.constraints[0].lhs()));
    ASSERT_EQ(result.constraints[0].rhs(), ipo::rational(2));

    vector[0] = 0.5;
    vector[1] = 0.5;
    result = oracle->separate(vector, true, query);
    ASSERT_TRUE(result.constraints.empty());

    result = oracle->separate(vector, false, query);
    ASSERT_FALSE(result.constraints.empty());
    ASSERT_TRUE(ipo::isMinusInfinity(result.constraints[0].lhs()));
    ASSERT_EQ(result.constraints[0].rhs(), ipo::rational(2));

    vector[0] = 3.0;
    vector[1] = -3.0;
      result = oracle->separate(vector, true, query);
    ASSERT_FALSE(result.constraints.empty());
    ASSERT_TRUE(ipo::isMinusInfinity(result.constraints[0].lhs()));
    ASSERT_EQ(result.constraints[0].rhs(), ipo::rational(2));
  }

#endif /* IPO_WITH_GMP && IPO_WITH_SOPLEX */

}
