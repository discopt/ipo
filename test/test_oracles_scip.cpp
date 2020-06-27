#include <gtest/gtest.h>
#include <ipo/oracles_scip.hpp>

#include <scip/scipdefplugins.h>

TEST(SCIP, Unbounded)
{
  std::cout << "===== SCIP.Unbounded ===== " << std::endl;
  
  SCIP* scip;
  SCIPcreate(&scip);
  SCIPincludeDefaultPlugins(scip);
  SCIPcreateProbBasic(scip, "unbounded");
  SCIP_VAR* x;
  SCIPcreateVarBasic(scip, &x, "x", 0.0, SCIPinfinity(scip), -1.0, SCIP_VARTYPE_INTEGER);
  SCIPaddVar(scip, x);

  SCIPreleaseVar(scip, &x);

  auto solver = std::make_shared<ipo::SCIPSolver>(scip);
  auto oracle = solver->getOptimizationOracle();

  double obj[] = {1.0};
  ipo::OptimizationOracle::Query query;
  ipo::OptimizationOracle::Result result;
  oracle->maximize(obj, query, result);
  ASSERT_TRUE(result.isUnbounded());
  ASSERT_TRUE(result.dualBound.isPlusInfinity());
}

TEST(SCIP, Infeasible)
{
  std::cout << "===== SCIP.Infeasible ===== " << std::endl;
  
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
  auto oracle = solver->getOptimizationOracle();

  double obj[] = { 1.0, 2.0 };
  ipo::OptimizationOracle::Query query;
  ipo::OptimizationOracle::Result result;
  oracle->maximize(obj, query, result);
  ASSERT_FALSE(result.isUnbounded());
  ASSERT_FALSE(result.isFeasible());
  ASSERT_TRUE(result.dualBound.isMinusInfinity());
}


TEST(SCIP, Separate)
{
  std::cout << "===== SCIP.Separate ===== " << std::endl;
  
  // x in [-1,1]
  // y in [-1,1]
  // -1 <= x + y <= 1
  

  SCIP* scip;
  SCIPcreate(&scip);
  SCIPincludeDefaultPlugins(scip);
  SCIPcreateProbBasic(scip, "simplex");
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

  auto solver = std::make_shared<ipo::SCIPSolver>(scip);
  auto oracle = solver->getSeparationOracle();

  double vector[] = { 1.0, 1.0 };
  ipo::SeparationOracle::Query query;
  ipo::SeparationOracle::Result result;
  oracle->separate(vector, true, query, result);
  ASSERT_FALSE(result.constraints.empty());

  result.reset();
  vector[0] = 0.5;
  vector[1] = 0.5;
  oracle->separate(vector, true, query, result);
  ASSERT_TRUE(result.constraints.empty());
  oracle->separate(vector, false, query, result);
  ASSERT_FALSE(result.constraints.empty());

  result.reset();
  vector[0] = 2.0;
  vector[1] = -2.0;
  oracle->separate(vector, true, query, result);
  ASSERT_FALSE(result.constraints.empty());
}

#if defined(IPO_WITH_GMP) && defined(IPO_WITH_SOPLEX)

TEST(SCIP, MakeRational)
{
  {
    std::cout << "===== SCIP.MakeRational.Bounded ===== " << std::endl;
    // max 3x + 4y + 2z
    // x in [0,1]
    // y in [0,1]
    // z in {0,1}
    // 7x+7y+7z <= 20
    
    SCIP* scip;
    SCIPcreate(&scip);
    SCIPincludeDefaultPlugins(scip);
    SCIPcreateProbBasic(scip, "Bounded");
    SCIP_VAR* x;
    SCIPcreateVarBasic(scip, &x, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS);
    SCIPaddVar(scip, x);
    SCIP_VAR* y;
    SCIPcreateVarBasic(scip, &y, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS);
    SCIPaddVar(scip, y);
    SCIP_VAR* z;
    SCIPcreateVarBasic(scip, &z, "z", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY);
    SCIPaddVar(scip, z);
    SCIP_CONS* cons;
    SCIPcreateConsBasicLinear(scip, &cons, "cons", 0, nullptr, nullptr, 0.0, 20.0);
    SCIPaddCoefLinear(scip, cons, x, 7.0);
    SCIPaddCoefLinear(scip, cons, y, 7.0);
    SCIPaddCoefLinear(scip, cons, z, 7.0);
    SCIPaddCons(scip, cons);
    SCIPreleaseVar(scip, &x);
    SCIPreleaseVar(scip, &y);
    SCIPreleaseVar(scip, &z);
    SCIPreleaseCons(scip, &cons);

    auto solver = std::make_shared<ipo::SCIPSolver>(scip);
    auto oracle = solver->getOptimizationOracle();
    const ipo::Space& space = *(solver->space());
    
    double vector[3];
    for (std::size_t i = 0; i < 3; ++i)
    {
      if (space[i] == "x")
        vector[i] = 3.0;
      else if (space[i] == "y")
        vector[i] = 4.0;
      else if (space[i] == "z")
        vector[i] = 2.0;
    }
    ipo::OptimizationOracle::Query query;
    ipo::OptimizationOracle::Result result;
    query.rational = true;
    oracle->maximize(vector, query, result);

    ASSERT_FALSE(result.points.empty());
    std::map<std::string, mpq_class> point;
    const ipo::Vector& best = result.points.front().vector;
    for (std::size_t i = 0; i < best.size(); ++i)
      point[space[best.coordinate(i)]] = best.rational(i);

    ASSERT_EQ(point["x"], mpq_class(6, 7));
    ASSERT_EQ(point["y"], mpq_class(1));
    ASSERT_EQ(point["z"], mpq_class(1));
  }

  {
    std::cout << "===== SCIP.MakeRational.Unbounded ===== " << std::endl;
    // max y
    // x >= 0
    // y >= 0
    // 3x - 4y = 1

    SCIP* scip;
    SCIPcreate(&scip);
    SCIPincludeDefaultPlugins(scip);
    SCIPcreateProbBasic(scip, "Unbounded");
    SCIP_VAR* x;
    SCIPcreateVarBasic(scip, &x, "x", 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS);
    SCIPaddVar(scip, x);
    SCIP_VAR* y;
    SCIPcreateVarBasic(scip, &y, "y", 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS);
    SCIPaddVar(scip, y);
    SCIP_CONS* cons;
    SCIPcreateConsBasicLinear(scip, &cons, "cons", 0, nullptr, nullptr, 1.0, 1.0);
    SCIPaddCoefLinear(scip, cons, x, 3.0);
    SCIPaddCoefLinear(scip, cons, y, -4.0);
    SCIPaddCons(scip, cons);
    SCIPreleaseVar(scip, &x);
    SCIPreleaseVar(scip, &y);
    SCIPreleaseCons(scip, &cons);

    auto solver = std::make_shared<ipo::SCIPSolver>(scip);
    auto oracle = solver->getOptimizationOracle();
    const ipo::Space& space = *(solver->space());

    double vector[2];
    for (std::size_t i = 0; i < 2; ++i)
    {
      if (space[i] == "x")
        vector[i] = 1.0;
      else if (space[i] == "y")
        vector[i] = 0.0;
    }
    ipo::OptimizationOracle::Query query;
    ipo::OptimizationOracle::Result result;
    query.rational = true;
    oracle->maximize(vector, query, result);

    ASSERT_FALSE(result.rays.empty());
    ASSERT_TRUE(result.isUnbounded());
    std::map<std::string, mpq_class> ray;
    ray[space[result.rays.begin()->vector.coordinate(0)]] = mpq_class(result.rays.begin()->vector.rational(0));
    ray[space[result.rays.begin()->vector.coordinate(1)]] = mpq_class(result.rays.begin()->vector.rational(1));
    ASSERT_EQ(ray["x"] / ray["y"], mpq_class(4, 3));
  }
}

#endif /* IPO_WITH_GMP && IPO_WITH_SOPLEX */
