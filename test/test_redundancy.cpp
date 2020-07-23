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
    ASSERT_EQ(red.test(eq1, 1.0e-12, 1.0e-12), ipo::EQUATION_INDEPENDENT);
    ASSERT_EQ(red.add(eq1, 1.0e-12, 1.0e-12, 1.0e-12), ipo::EQUATION_INDEPENDENT);
    ASSERT_EQ(red.rank(), 1);
    
    vector = std::make_shared<sparse_vector<double>>();
    vector->push_back(0, 1.0);
    vector->push_back(1, 2.0);
    vector->push_back(2, 4.0);
    auto eq2 = ipo::Constraint<double>(1.0, vector, 1.0);
    ASSERT_EQ(red.test(eq2, 1.0e-12, 1.0e-12), ipo::EQUATION_INDEPENDENT);
    ASSERT_EQ(red.add(eq2, 1.0e-12, 1.0e-12, 1.0e-12), ipo::EQUATION_INDEPENDENT);
    ASSERT_EQ(red.rank(), 2);

    vector = std::make_shared<sparse_vector<double>>();
    vector->push_back(0, -3.0);
    vector->push_back(1, -6.0);
    vector->push_back(2, -10.0);
    auto eq3 = ipo::Constraint<double>(-3.0, vector, -3.0);
    ASSERT_EQ(red.test(eq3, 1.0e-12, 1.0e-12), ipo::EQUATION_REDUNDANT);
    ASSERT_EQ(red.add(eq3, 1.0e-12, 1.0e-12, 1.0e-12), ipo::EQUATION_REDUNDANT);
    ASSERT_EQ(red.rank(), 2);

    vector = std::make_shared<sparse_vector<double>>();
    vector->push_back(0, -3.0);
    vector->push_back(1, -6.0);
    vector->push_back(2, -10.0);
    auto eq4 = ipo::Constraint<double>(0.0, vector, 0.0);
    ASSERT_EQ(red.test(eq4, 1.0e-12, 1.0e-12), ipo::EQUATION_INCONSISTENT);
    ASSERT_EQ(red.add(eq4, 1.0e-12, 1.0e-12, 1.0e-12), ipo::EQUATION_INCONSISTENT);
    ASSERT_EQ(red.rank(), 2);

    vector = std::make_shared<sparse_vector<double>>();
    vector->push_back(0, -3.0);
    vector->push_back(1, -5.0);
    vector->push_back(2, -10.0);
    auto eq5 = ipo::Constraint<double>(0.0, vector, 0.0);
    ASSERT_EQ(red.test(eq5, 1.0e-12, 1.0e-12), ipo::EQUATION_INDEPENDENT);
    ASSERT_EQ(red.add(eq5, 1.0e-12, 1.0e-12, 1.0e-12), ipo::EQUATION_INDEPENDENT);
    ASSERT_EQ(red.rank(), 3);
  }

#if defined(IPO_WITH_GMP)
  std::cout << "===== LinearAlgebra::EquationRedundancyCheck::Rational ===== " << std::endl;
  {
    auto red = ipo::EquationRedundancyCheck<mpq_class>(3);
    ASSERT_EQ(red.rank(), 0);

    auto vector = std::make_shared<sparse_vector<mpq_class>>();
    vector->push_back(0, 1);
    vector->push_back(1, 2);
    vector->push_back(2, 3);
    auto eq1 = ipo::Constraint<mpq_class>(mpq_class(1), vector, mpq_class(1));
    ASSERT_EQ(red.test(eq1, 0, 0), ipo::EQUATION_INDEPENDENT);
    ASSERT_EQ(red.add(eq1, 0, 0, 0), ipo::EQUATION_INDEPENDENT);
    ASSERT_EQ(red.rank(), 1);

    vector = std::make_shared<sparse_vector<mpq_class>>();
    vector->push_back(0, 1);
    vector->push_back(1, 2);
    vector->push_back(2, 4);
    auto eq2 = ipo::Constraint<mpq_class>(mpq_class(1), vector, mpq_class(1));
    ASSERT_EQ(red.test(eq2, 0, 0), ipo::EQUATION_INDEPENDENT);
    ASSERT_EQ(red.add(eq2, 0, 0, 0), ipo::EQUATION_INDEPENDENT);
    ASSERT_EQ(red.rank(), 2);

    vector = std::make_shared<sparse_vector<mpq_class>>();
    vector->push_back(0, -3);
    vector->push_back(1, -6);
    vector->push_back(2, -10);
    auto eq3 = ipo::Constraint<mpq_class>(mpq_class(-3), vector, mpq_class(-3));
    ASSERT_EQ(red.test(eq3, 0, 0), ipo::EQUATION_REDUNDANT);
    ASSERT_EQ(red.add(eq3, 0, 0, 0), ipo::EQUATION_REDUNDANT);
    ASSERT_EQ(red.rank(), 2);

    vector = std::make_shared<sparse_vector<mpq_class>>();
    vector->push_back(0, -3);
    vector->push_back(1, -6);
    vector->push_back(2, -10);
    auto eq4 = ipo::Constraint<mpq_class>(mpq_class(0), vector, mpq_class(0));
    ASSERT_EQ(red.test(eq4, 0, 0), ipo::EQUATION_INCONSISTENT);
    ASSERT_EQ(red.add(eq4, 0, 0, 0), ipo::EQUATION_INCONSISTENT);
    ASSERT_EQ(red.rank(), 2);

    vector = std::make_shared<sparse_vector<mpq_class>>();
    vector->push_back(0, -3);
    vector->push_back(1, -5);
    vector->push_back(2, -10);
    auto eq5 = ipo::Constraint<mpq_class>(mpq_class(0), vector, mpq_class(0));
    ASSERT_EQ(red.test(eq5, 0, 0), ipo::EQUATION_INDEPENDENT);
    ASSERT_EQ(red.add(eq5, 0, 0, 0), ipo::EQUATION_INDEPENDENT);
    ASSERT_EQ(red.rank(), 3);
  }

#endif /* IPO_WITH_GMP */

}
