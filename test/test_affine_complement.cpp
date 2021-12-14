#include <gtest/gtest.h>

#include <ipo/arithmetic.hpp>

#include <../src/ipo/affine_hull.cpp>

TEST(AffineHull, AffineComplement)
{
  std::cout << "===== AffineHull::AffineComplement::Double ===== " << std::endl;
  {
    auto ac = ipo::AffineComplement<double>(4);
    ASSERT_EQ(ac.rank(), 0);
    sparse_vector<double> vector;

    vector.clear();
    vector.push_back(1, 1.0);
    ac.add(vector, 1.0, 1, 1.0e-12);
    ASSERT_EQ(ac.rank(), 1);

    vector.clear();
    vector.push_back(0, 1.0);
    vector.push_back(1, 1.0);
    ac.add(vector, 1.0, 0, 1.0e-12);
    ASSERT_EQ(ac.rank(), 2);

    vector.clear();
    vector.push_back(3, 1.0);
    ac.add(vector, 1.0, 4, 1.0e-12);
    ASSERT_EQ(ac.rank(), 3);
  }

#if defined(IPO_RATIONAL)
  std::cout << "===== AffineHull::AffineComplement::Rational ===== " << std::endl;
  {
    auto ac = ipo::AffineComplement<ipo::rational>(4);
    ASSERT_EQ(ac.rank(), 0);
    sparse_vector<ipo::rational> vector;

    vector.push_back(1, ipo::rational(std::move(ipo::rational(1))));
    ac.add(vector, 1, 1, 0.0);
    ASSERT_EQ(ac.rank(), 1);

    vector.clear();
    vector.push_back(0, 1);
    vector.push_back(1, 1);
    ac.add(vector, 1, 0, 0.0);
    ASSERT_EQ(ac.rank(), 2);

    vector.clear();
    vector.push_back(3, 1);
    ac.add(vector, 1, 4, 0.0);
    ASSERT_EQ(ac.rank(), 3);
  }

#endif /* IPO_RATIONAL */
}
