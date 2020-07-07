#include <gtest/gtest.h>

#include <ipo/rational.hpp>

#include <../src/ipo/lu.hpp>

TEST(LinearAlgebra, LU)
{
  std::cout << "===== LinearAlgebra::IncrementalLUFactorization::Double ===== " << std::endl;
  {
    auto isZero = ipo::DoubleIsZero(1.0e-9);
    auto lu = ipo::IncrementalLUFactorization<double, ipo::DoubleIsZero>(isZero);

    /* Test matrix: 
    *  1 2 3
    *  4 5 6
    *  7 8 0
    */
    
    std::vector<double> row(3);
    std::vector<double> column(3);
    lu.extend(&row[0], &column[0], 1);
    row[0] = 4;
    column[0] = 2;
    lu.extend(&row[0], &column[0], 5);
    row[0] = 7;
    row[1] = 8;
    column[0] = 3;
    column[1] = 6;
    lu.extend(&row[0], &column[0], 0);

    column[0] = 1;
    column[1] = 0;
    column[2] = 0;
    lu.solveLeft(&column[0]);
    lu.solveUpper(&column[0]);
    ASSERT_NEAR(column[0], -16.0/9, 1.0e-9);
    ASSERT_NEAR(column[1], 14.0/9, 1.0e-9);
    ASSERT_NEAR(column[2], -1.0/9, 1.0e-9);
    column[0] = 0;
    column[1] = 1;
    column[2] = 0;
    lu.solveLeft(&column[0]);
    lu.solveUpper(&column[0]);
    ASSERT_NEAR(column[0], 8.0/9, 1.0e-9);
    ASSERT_NEAR(column[1], -7.0/9, 1.0e-9);
    ASSERT_NEAR(column[2], 2.0/9, 1.0e-9);
    column[0] = 0;
    column[1] = 0;
    column[2] = 1;
    lu.solveLeft(&column[0]);
    lu.solveUpper(&column[0]);
    ASSERT_NEAR(column[0], -1.0/9, 1.0e-9);
    ASSERT_NEAR(column[1], 2.0/9, 1.0e-9);
    ASSERT_NEAR(column[2], -1.0/9, 1.0e-9);
  }

#if defined(IPO_WITH_GMP)
  std::cout << "===== LinearAlgebra::IncrementalLUFactorization::Rational ===== " << std::endl;
  {
    auto lu = ipo::IncrementalLUFactorization<ipo::rational, ipo::RationalIsZero>(ipo::RationalIsZero());

    /* Test matrix: 
    *  1 2 3
    *  4 5 6
    *  7 8 0
    */
    
    std::vector<ipo::rational> row(3);
    std::vector<ipo::rational> column(3);
    lu.extend(&row[0], &column[0], 1);
    row[0] = 4;
    column[0] = 2;
    lu.extend(&row[0], &column[0], 5);
    row[0] = 7;
    row[1] = 8;
    column[0] = 3;
    column[1] = 6;
    lu.extend(&row[0], &column[0], 0);

    column[0] = 1;
    column[1] = 0;
    column[2] = 0;
    lu.solveLeft(&column[0]);
    lu.solveUpper(&column[0]);
    ASSERT_EQ(column[0], mpq_class(-16,9));
    ASSERT_EQ(column[1], mpq_class(14,9));
    ASSERT_EQ(column[2], mpq_class(-1,9));
    column[0] = 0;
    column[1] = 1;
    column[2] = 0;
    lu.solveLeft(&column[0]);
    lu.solveUpper(&column[0]);
    ASSERT_EQ(column[0], mpq_class(8,9));
    ASSERT_EQ(column[1], mpq_class(-7,9));
    ASSERT_EQ(column[2], mpq_class(2,9));
    column[0] = 0;
    column[1] = 0;
    column[2] = 1;
    lu.solveLeft(&column[0]);
    lu.solveUpper(&column[0]);
    ASSERT_EQ(column[0], mpq_class(-1,9));
    ASSERT_EQ(column[1], mpq_class(2,9));
    ASSERT_EQ(column[2], mpq_class(-1,9));
  }

#endif /* IPO_WITH_GMP */
}
