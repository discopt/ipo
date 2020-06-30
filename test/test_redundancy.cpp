#include <gtest/gtest.h>
#include <../src/ipo/redundancy.hpp>

TEST(Redundancy, Real)
{
  auto red = ipo::EquationRedundancyCheck<double, ipo::RealIsZero>(3, ipo::RealIsZero(1.0e-9));
  ASSERT_EQ(red.rank(), 0);

  std::vector<double> vector(3);
  vector[0] = 1.0;
  vector[1] = 2.0;
  vector[2] = 3.0;
  auto eq1 = ipo::Constraint(1.0, ipo::Vector(vector), 1.0);
  ASSERT_EQ(red.test(eq1), ipo::EQUATION_INDEPENDENT);
  ASSERT_EQ(red.add(eq1), ipo::EQUATION_INDEPENDENT);
  ASSERT_EQ(red.rank(), 1);
  vector[0] = 1.0;
  vector[1] = 2.0;
  vector[2] = 4.0;
  auto eq2 = ipo::Constraint(1.0, ipo::Vector(vector), 1.0);
  ASSERT_EQ(red.test(eq2), ipo::EQUATION_INDEPENDENT);
  ASSERT_EQ(red.add(eq2), ipo::EQUATION_INDEPENDENT);
  ASSERT_EQ(red.rank(), 2);
  vector[0] = -3.0;
  vector[1] = -6.0;
  vector[2] = -10.0;
  auto eq3 = ipo::Constraint(-3.0, ipo::Vector(vector), -3.0);
  ASSERT_EQ(red.test(eq3), ipo::EQUATION_REDUNDANT);
  ASSERT_EQ(red.add(eq3), ipo::EQUATION_REDUNDANT);
  ASSERT_EQ(red.rank(), 2);
  vector[0] = -3.0;
  vector[1] = -6.0;
  vector[2] = -10.0;
  auto eq4 = ipo::Constraint(0.0, ipo::Vector(vector), 0.0);
  ASSERT_EQ(red.test(eq4), ipo::EQUATION_INFEASIBLE);
  ASSERT_EQ(red.add(eq4), ipo::EQUATION_INFEASIBLE);
  ASSERT_EQ(red.rank(), 2);
  vector[0] = -3.0;
  vector[1] = -5.0;
  vector[2] = -10.0;
  auto eq5 = ipo::Constraint(0.0, ipo::Vector(vector), 0.0);
  ASSERT_EQ(red.test(eq5), ipo::EQUATION_INDEPENDENT);
  ASSERT_EQ(red.add(eq5), ipo::EQUATION_INDEPENDENT);
  ASSERT_EQ(red.rank(), 3);
}


TEST(Redundancy, Rational)
{
  auto red = ipo::EquationRedundancyCheck<mpq_class, ipo::RationalIsZero>(3, ipo::RationalIsZero());
  ASSERT_EQ(red.rank(), 0);

  std::vector<mpq_class> vector(3);
  vector[0] = 1;
  vector[1] = 2;
  vector[2] = 3;
  auto eq1 = ipo::Constraint(mpq_class(1), ipo::Vector(vector), mpq_class(1));
  ASSERT_EQ(red.test(eq1), ipo::EQUATION_INDEPENDENT);
  ASSERT_EQ(red.add(eq1), ipo::EQUATION_INDEPENDENT);
  ASSERT_EQ(red.rank(), 1);
  vector[0] = 1;
  vector[1] = 2;
  vector[2] = 4;
  auto eq2 = ipo::Constraint(mpq_class(1), ipo::Vector(vector), mpq_class(1));
  ASSERT_EQ(red.test(eq2), ipo::EQUATION_INDEPENDENT);
  ASSERT_EQ(red.add(eq2), ipo::EQUATION_INDEPENDENT);
  ASSERT_EQ(red.rank(), 2);
  vector[0] = -3;
  vector[1] = -6;
  vector[2] = -10;
  auto eq3 = ipo::Constraint(mpq_class(-3), ipo::Vector(vector), mpq_class(-3));
  ASSERT_EQ(red.test(eq3), ipo::EQUATION_REDUNDANT);
  ASSERT_EQ(red.add(eq3), ipo::EQUATION_REDUNDANT);
  ASSERT_EQ(red.rank(), 2);
  vector[0] = -3;
  vector[1] = -6;
  vector[2] = -10;
  auto eq4 = ipo::Constraint(mpq_class(0), ipo::Vector(vector), mpq_class(0));
  ASSERT_EQ(red.test(eq4), ipo::EQUATION_INFEASIBLE);
  ASSERT_EQ(red.add(eq4), ipo::EQUATION_INFEASIBLE);
  ASSERT_EQ(red.rank(), 2);
  vector[0] = -3;
  vector[1] = -5;
  vector[2] = -10;
  auto eq5 = ipo::Constraint(mpq_class(0), ipo::Vector(vector), mpq_class(0));
  ASSERT_EQ(red.test(eq5), ipo::EQUATION_INDEPENDENT);
  ASSERT_EQ(red.add(eq5), ipo::EQUATION_INDEPENDENT);
  ASSERT_EQ(red.rank(), 3);
}
