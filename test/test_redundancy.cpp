#include <gtest/gtest.h>
#include <../src/ipo/redundancy.hpp>

TEST(LinearAlgebra, EquationRedundancyCheck)
{
  std::cout << "===== LinearAlgebra::EquationRedundancyCheck::Double ===== " << std::endl;
  {
    auto red = ipo::EquationRedundancyCheck<double, ipo::DoubleIsZero>(3, ipo::DoubleIsZero(1.0e-9));
    ASSERT_EQ(red.rank(), 0);
    sparse_vector<double> vector;

    vector.clear();
    vector.push_back(0, 1.0);
    vector.push_back(1, 2.0);
    vector.push_back(2, 3.0);
    auto eq1 = ipo::Constraint<double>(1.0, vector, 1.0);
    ASSERT_EQ(red.test(eq1), ipo::EQUATION_INDEPENDENT);
    ASSERT_EQ(red.add(eq1), ipo::EQUATION_INDEPENDENT);
    ASSERT_EQ(red.rank(), 1);
    
    vector.clear();
    vector.push_back(0, 1.0);
    vector.push_back(1, 2.0);
    vector.push_back(2, 4.0);
    auto eq2 = ipo::Constraint<double>(1.0, vector, 1.0);
    ASSERT_EQ(red.test(eq2), ipo::EQUATION_INDEPENDENT);
    ASSERT_EQ(red.add(eq2), ipo::EQUATION_INDEPENDENT);
    ASSERT_EQ(red.rank(), 2);

    vector.clear();
    vector.push_back(0, -3.0);
    vector.push_back(1, -6.0);
    vector.push_back(2, -10.0);
    auto eq3 = ipo::Constraint<double>(-3.0, vector, -3.0);
    ASSERT_EQ(red.test(eq3), ipo::EQUATION_REDUNDANT);
    ASSERT_EQ(red.add(eq3), ipo::EQUATION_REDUNDANT);
    ASSERT_EQ(red.rank(), 2);

    vector.clear();
    vector.push_back(0, -3.0);
    vector.push_back(1, -6.0);
    vector.push_back(2, -10.0);
    auto eq4 = ipo::Constraint<double>(0.0, vector, 0.0);
    ASSERT_EQ(red.test(eq4), ipo::EQUATION_INCONSISTENT);
    ASSERT_EQ(red.add(eq4), ipo::EQUATION_INCONSISTENT);
    ASSERT_EQ(red.rank(), 2);

    vector.clear();
    vector.push_back(0, -3.0);
    vector.push_back(1, -5.0);
    vector.push_back(2, -10.0);
    auto eq5 = ipo::Constraint<double>(0.0, vector, 0.0);
    ASSERT_EQ(red.test(eq5), ipo::EQUATION_INDEPENDENT);
    ASSERT_EQ(red.add(eq5), ipo::EQUATION_INDEPENDENT);
    ASSERT_EQ(red.rank(), 3);
  }

#if defined(IPO_WITH_GMP)
  std::cout << "===== LinearAlgebra::EquationRedundancyCheck::Rational ===== " << std::endl;
  {
    auto red = ipo::EquationRedundancyCheck<ipo::rational, ipo::RationalIsZero>(3, ipo::RationalIsZero());
    ASSERT_EQ(red.rank(), 0);
    sparse_vector<ipo::rational> vector;

    vector.clear();
    vector.push_back(0, 1);
    vector.push_back(1, 2);
    vector.push_back(2, 3);
    auto eq1 = ipo::Constraint<ipo::rational>(ipo::rational(1), vector, ipo::rational(1));
    ASSERT_EQ(red.test(eq1), ipo::EQUATION_INDEPENDENT);
    ASSERT_EQ(red.add(eq1), ipo::EQUATION_INDEPENDENT);
    ASSERT_EQ(red.rank(), 1);

    vector.clear();
    vector.push_back(0, 1);
    vector.push_back(1, 2);
    vector.push_back(2, 4);
    auto eq2 = ipo::Constraint<ipo::rational>(ipo::rational(1), vector, ipo::rational(1));
    ASSERT_EQ(red.test(eq2), ipo::EQUATION_INDEPENDENT);
    ASSERT_EQ(red.add(eq2), ipo::EQUATION_INDEPENDENT);
    ASSERT_EQ(red.rank(), 2);

    vector.clear();
    vector.push_back(0, -3);
    vector.push_back(1, -6);
    vector.push_back(2, -10);
    auto eq3 = ipo::Constraint<ipo::rational>(ipo::rational(-3), vector, ipo::rational(-3));
    ASSERT_EQ(red.test(eq3), ipo::EQUATION_REDUNDANT);
    ASSERT_EQ(red.add(eq3), ipo::EQUATION_REDUNDANT);
    ASSERT_EQ(red.rank(), 2);

    vector.clear();
    vector.push_back(0, -3);
    vector.push_back(1, -6);
    vector.push_back(2, -10);
    auto eq4 = ipo::Constraint<ipo::rational>(ipo::rational(0), vector, ipo::rational(0));
    ASSERT_EQ(red.test(eq4), ipo::EQUATION_INCONSISTENT);
    ASSERT_EQ(red.add(eq4), ipo::EQUATION_INCONSISTENT);
    ASSERT_EQ(red.rank(), 2);

    vector.clear();
    vector.push_back(0, -3);
    vector.push_back(1, -5);
    vector.push_back(2, -10);
    auto eq5 = ipo::Constraint<ipo::rational>(ipo::rational(0), vector, ipo::rational(0));
    ASSERT_EQ(red.test(eq5), ipo::EQUATION_INDEPENDENT);
    ASSERT_EQ(red.add(eq5), ipo::EQUATION_INDEPENDENT);
    ASSERT_EQ(red.rank(), 3);
  }

#endif /* IPO_WITH_GMP */

}
