#include <gtest/gtest.h>

#include <ipo/arithmetic.hpp>

#include <../src/ipo/redundancy.hpp>

TEST(LinearAlgebra, EquationRedundancyCheck)
{
  std::cout << "===== LinearAlgebra::EquationRedundancyCheck::Double ===== " << std::endl;
  {
    auto red = ipo::EquationRedundancyCheck<double>(3);
    ASSERT_EQ(red.rank(), 0);
    auto vector = std::make_shared<sparse_vector<double>>();
    vector->push_back(0, 1.0);
    vector->push_back(1, 2.0);
    vector->push_back(2, 3.0);
    auto eq1 = ipo::Constraint<double>(1.0, vector, 1.0);
    auto result = red.test(eq1);
    ASSERT_GT(result.maxViolation, 1.0e-12);
    ASSERT_TRUE(red.add(eq1, result.maxCoordinate, 1.0e-12));
    ASSERT_EQ(red.rank(), 1);
    
    vector = std::make_shared<sparse_vector<double>>();
    vector->push_back(0, 1.0);
    vector->push_back(1, 2.0);
    vector->push_back(2, 4.0);
    auto eq2 = ipo::Constraint<double>(1.0, vector, 1.0);
    result = red.test(eq2);
    ASSERT_GT(result.maxViolation, 1.0e-12);
    ASSERT_TRUE(red.add(eq2, result.maxCoordinate, 1.0e-12));
    ASSERT_EQ(red.rank(), 2);

    vector = std::make_shared<sparse_vector<double>>();
    vector->push_back(0, -3.0);
    vector->push_back(1, -6.0);
    vector->push_back(2, -10.0);
    auto eq3 = ipo::Constraint<double>(-3.0, vector, -3.0);
    result = red.test(eq3);
    ASSERT_LE(result.maxViolation, 1.0e-12);
    ASSERT_LE(fabs(result.rhs), 1.0e-12);
    ASSERT_EQ(red.rank(), 2);

    vector = std::make_shared<sparse_vector<double>>();
    vector->push_back(0, -3.0);
    vector->push_back(1, -6.0);
    vector->push_back(2, -10.0);
    auto eq4 = ipo::Constraint<double>(0.0, vector, 0.0);
    result = red.test(eq4);
    ASSERT_LE(result.maxViolation, 1.0e-12);
    ASSERT_GT(fabs(result.rhs), 1.0e-12);
    ASSERT_EQ(red.rank(), 2);

    vector = std::make_shared<sparse_vector<double>>();
    vector->push_back(0, -3.0);
    vector->push_back(1, -5.0);
    vector->push_back(2, -10.0);
    auto eq5 = ipo::Constraint<double>(0.0, vector, 0.0);
    result = red.test(eq5);
    ASSERT_GT(result.maxViolation, 1.0e-12);
    ASSERT_TRUE(red.add(eq5, result.maxCoordinate, 1.0e-12));
    ASSERT_EQ(red.rank(), 3);
  }

#if defined(IPO_RATIONAL)
  std::cout << "===== LinearAlgebra::EquationRedundancyCheck::Rational ===== " << std::endl;
  {
    auto red = ipo::EquationRedundancyCheck<ipo::rational>(3);
    ASSERT_EQ(red.rank(), 0);

    auto vector = std::make_shared<sparse_vector<ipo::rational>>();
    vector->push_back(0, 1);
    vector->push_back(1, 2);
    vector->push_back(2, 3);
    auto eq1 = ipo::Constraint<ipo::rational>(ipo::rational(1), vector, ipo::rational(1));
    auto result = red.test(eq1);
    ASSERT_GT(result.maxViolation, 0);
    ASSERT_TRUE(red.add(eq1, result.maxCoordinate, 0));
    ASSERT_EQ(red.rank(), 1);

    vector = std::make_shared<sparse_vector<ipo::rational>>();
    vector->push_back(0, 1);
    vector->push_back(1, 2);
    vector->push_back(2, 4);
    auto eq2 = ipo::Constraint<ipo::rational>(ipo::rational(1), vector, ipo::rational(1));
    result = red.test(eq2);
    ASSERT_GT(result.maxViolation, 0);
    ASSERT_TRUE(red.add(eq2, result.maxCoordinate, 0));
    ASSERT_EQ(red.rank(), 2);

    vector = std::make_shared<sparse_vector<ipo::rational>>();
    vector->push_back(0, -3);
    vector->push_back(1, -6);
    vector->push_back(2, -10);
    auto eq3 = ipo::Constraint<ipo::rational>(ipo::rational(-3), vector, ipo::rational(-3));
    result = red.test(eq3);
    ASSERT_EQ(result.maxViolation, 0);
    ASSERT_EQ(result.rhs, 0);
    ASSERT_EQ(red.rank(), 2);

    vector = std::make_shared<sparse_vector<ipo::rational>>();
    vector->push_back(0, -3);
    vector->push_back(1, -6);
    vector->push_back(2, -10);
    auto eq4 = ipo::Constraint<ipo::rational>(ipo::rational(0), vector, ipo::rational(0));
    result = red.test(eq4);
    ASSERT_EQ(result.maxViolation, 0);
    ASSERT_NE(result.rhs, 0);
    ASSERT_EQ(red.rank(), 2);

    vector = std::make_shared<sparse_vector<ipo::rational>>();
    vector->push_back(0, -3);
    vector->push_back(1, -5);
    vector->push_back(2, -10);
    auto eq5 = ipo::Constraint<ipo::rational>(ipo::rational(0), vector, ipo::rational(0));
    result = red.test(eq5);
    ASSERT_GT(result.maxViolation, 0);
    ASSERT_TRUE(red.add(eq5, result.maxCoordinate, 0));
    ASSERT_EQ(red.rank(), 3);
  }

#endif /* IPO_RATIONAL */

}
