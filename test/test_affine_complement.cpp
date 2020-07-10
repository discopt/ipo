#include <gtest/gtest.h>

#include <ipo/rational.hpp>

#include <../src/ipo/affine_hull.cpp>

TEST(AffineHull, AffineComplement)
{
  std::cout << "===== AffineHull::AffineComplement::Double ===== " << std::endl;
  {
    auto ac = ipo::AffineComplement<double, ipo::DoubleIsZero>(4, ipo::DoubleIsZero(1.0e-9));
    ASSERT_EQ(ac.rank(), 0);
    sparse_vector<double> vector;

    vector.clear();
    vector.push_back(1, 1.0);
    ac.add(vector, 1.0, 1);
    ASSERT_EQ(ac.rank(), 1);

    vector.clear();
    vector.push_back(0, 1.0);
    vector.push_back(1, 1.0);
    ac.add(vector, 1.0, 0);
    ASSERT_EQ(ac.rank(), 2);

    vector.clear();
    vector.push_back(3, 1.0);
    ac.add(vector, 1.0, 4);
    ASSERT_EQ(ac.rank(), 3);
  }

#if defined(IPO_WITH_GMP)
  std::cout << "===== AffineHull::AffineComplement::Rational ===== " << std::endl;
  {
    auto ac = ipo::AffineComplement<ipo::rational, ipo::RationalIsZero>(4);
    ASSERT_EQ(ac.rank(), 0);
    sparse_vector<ipo::rational> vector;

    vector.push_back(1, ipo::rational(std::move(mpq_class(1))));
    ac.add(vector, 1, 1);
    ASSERT_EQ(ac.rank(), 1);

    vector.clear();
    vector.push_back(0, 1);
    vector.push_back(1, 1);
    ac.add(vector, 1, 0);
    ASSERT_EQ(ac.rank(), 2);

    vector.clear();
    vector.push_back(3, 1);
    ac.add(vector, 1, 4);
    ASSERT_EQ(ac.rank(), 3);
  }

#endif /* IPO_WITH_GMP */
}
