#include <gtest/gtest.h>
#include <ipo/oracles_gurobi.hpp>

TEST(Gurobi, UnboundedDouble)
{
  // x >= 0
  // y >= 0
  // y == 2x + 1

  GRBenv* env = nullptr;
  GRBemptyenv(&env);
  GRBsetintparam(env, GRB_INT_PAR_OUTPUTFLAG, 0);
  GRBstartenv(env);
  GRBmodel* model = nullptr;
  double lbs[2] = { 0.0, 0.0 };
  char vtypes[2] = { 'I', 'I' };
  const char* names[2] = { "x", "y" };
  GRBnewmodel(env, &model, "UnboundedDouble", 2, nullptr, lbs, nullptr, vtypes, const_cast<char**>(&names[0]));
  int indices1[2] = { 0, 1 };
  double coefs1[2] = { 2.0, -1.0 };
  GRBaddconstr(model, 2, indices1, coefs1, '=', 1.0, nullptr);

  auto solver = std::make_shared<ipo::GurobiSolver>(std::move(model));
  auto oracle = solver->getOptimizationOracle<double>();

  double obj[] = { 1.0, 1.0 };
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

TEST(Gurobi, InfeasibleDouble)
{
  // x >= 0
  // y >= 0
  // 2x + 2y == 1

  GRBenv* env = nullptr;
  GRBemptyenv(&env);
  GRBsetintparam(env, GRB_INT_PAR_OUTPUTFLAG, 0);
  GRBstartenv(env);
  GRBmodel* model = nullptr;
  double lbs[2] = { 0.0, 0.0 };
  char vtypes[2] = { 'I', 'I' };
  const char* names[2] = { "x", "y" };
  GRBnewmodel(env, &model, "InfeasibleDouble", 2, nullptr, lbs, nullptr, vtypes, const_cast<char**>(&names[0]));
  int indices1[2] = { 0, 1 };
  double coefs1[2] = { 2.0, 2.0 };
  GRBaddconstr(model, 2, indices1, coefs1, '=', 1.0, nullptr);

  auto solver = std::make_shared<ipo::GurobiSolver>(std::move(model));
  auto oracle = solver->getOptimizationOracle<double>();

  double obj[] = { 1.0, 1.0 };
  ipo::OptimizationOracle<double>::Query query;
  std::cout << "oracle->maximize()" << std::endl;
  auto response = oracle->maximize(obj, query);
  ASSERT_EQ(response.outcome, ipo::OptimizationOutcome::INFEASIBLE);
  ASSERT_FALSE(response.hasDualBound);
}

TEST(Gurobi, SeparateDouble)
{
  // x in [-1.1,1.2]
  // y in [-1.3,1.4]
  // -1.5 <= x + y <= 1.5

  GRBenv* env = nullptr;
  GRBemptyenv(&env);
  GRBsetintparam(env, GRB_INT_PAR_OUTPUTFLAG, 0);
  GRBstartenv(env);
  GRBmodel* model = nullptr;
  double lbs[2] = { -1.1, -1.3 };
  double ubs[2] = {  1.2,  1.4 };
  char vtypes[2] = { 'I', 'I' };
  const char* names[2] = { "x", "y" };
  GRBnewmodel(env, &model, "SeparateDouble", 2, nullptr, nullptr, nullptr, vtypes, const_cast<char**>(&names[0]));
  GRBsetdblattrarray(model, GRB_DBL_ATTR_LB, 0, 2, lbs);
  GRBsetdblattrarray(model, GRB_DBL_ATTR_UB, 0, 2, ubs);
  int indices1[2] = { 0, 1 };
  double coefs1[2] = { 1.0, 1.0 };
  GRBaddconstr(model, 2, indices1, coefs1, '<', 1.5, nullptr);
  GRBaddconstr(model, 2, indices1, coefs1, '>', -1.5, nullptr);

  auto solver = std::make_shared<ipo::GurobiSolver>(std::move(model));
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


#if defined(IPO_RATIONAL_MIP_GUROBI)

TEST(Gurobi, UnboundedRational)
{
  // x >= 0
  // y >= 0
  // y == 2x + 1

  GRBenv* env = nullptr;
  GRBemptyenv(&env);
  GRBsetintparam(env, GRB_INT_PAR_OUTPUTFLAG, 0);
  GRBstartenv(env);
  GRBmodel* model = nullptr;
  double lbs[2] = { 0.0, 0.0 };
  char vtypes[2] = { 'I', 'I' };
  const char* names[2] = { "x", "y" };
  GRBnewmodel(env, &model, "UnboundedRational", 2, nullptr, lbs, nullptr, vtypes, const_cast<char**>(&names[0]));
  int indices1[2] = { 0, 1 };
  double coefs1[2] = { 2.0, -1.0 };
  GRBaddconstr(model, 2, indices1, coefs1, '=', 1.0, nullptr);

  auto solver = std::make_shared<ipo::GurobiSolver>(std::move(model));
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

TEST(Gurobi, InfeasibleRational)
{
  // x >= 0
  // y >= 0
  // 2x + 2y == 1

  GRBenv* env = nullptr;
  GRBemptyenv(&env);
  GRBsetintparam(env, GRB_INT_PAR_OUTPUTFLAG, 0);
  GRBstartenv(env);
  GRBmodel* model = nullptr;
  double lbs[2] = { 0.0, 0.0 };
  char vtypes[2] = { 'I', 'I' };
  const char* names[2] = { "x", "y" };
  GRBnewmodel(env, &model, "InfeasibleRational", 2, nullptr, lbs, nullptr, vtypes, const_cast<char**>(&names[0]));
  int indices1[2] = { 0, 1 };
  double coefs1[2] = { 2.0, 2.0 };
  GRBaddconstr(model, 2, indices1, coefs1, '=', 1.0, nullptr);

  auto solver = std::make_shared<ipo::GurobiSolver>(std::move(model));
  auto oracle = solver->getOptimizationOracle<ipo::rational>();

  ipo::rational obj[] = { 1, 2 };
  ipo::OptimizationOracle<ipo::rational>::Query query;
  auto response = oracle->maximize(obj, query);
  ASSERT_EQ(response.outcome, ipo::OptimizationOutcome::INFEASIBLE);
  ASSERT_FALSE(response.hasDualBound);
}

TEST(Gurobi, SeparateRational)
{
  // x in [-1,2]
  // y in [-4,8]
  // -1.5 <= x + y <= 1.5

  GRBenv* env = nullptr;
  GRBemptyenv(&env);
  GRBsetintparam(env, GRB_INT_PAR_OUTPUTFLAG, 0);
  GRBstartenv(env);
  GRBmodel* model = nullptr;
  double lbs[2] = { -1, -4 };
  double ubs[2] = { 2, 8 };
  char vtypes[2] = { 'I', 'I' };
  const char* names[2] = { "x", "y" };
  GRBnewmodel(env, &model, "SeparateRational", 2, nullptr, nullptr, nullptr, vtypes, const_cast<char**>(&names[0]));
  GRBsetdblattrarray(model, GRB_DBL_ATTR_LB, 0, 2, lbs);
  GRBsetdblattrarray(model, GRB_DBL_ATTR_UB, 0, 2, ubs);
  int indices1[2] = { 0, 1 };
  double coefs1[2] = { 1.0, 1.0 };
  GRBaddconstr(model, 2, indices1, coefs1, '<', 16, nullptr);
  GRBaddconstr(model, 2, indices1, coefs1, '>', -16, nullptr);

  auto solver = std::make_shared<ipo::GurobiSolver>(std::move(model));
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

#endif /* IPO_RATIONAL_MIP_GUROBI */
