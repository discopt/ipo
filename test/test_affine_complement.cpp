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

#if defined(IPO_WITH_GMP)
  std::cout << "===== AffineHull::AffineComplement::Rational ===== " << std::endl;
  {
    auto ac = ipo::AffineComplement<mpq_class>(4);
    ASSERT_EQ(ac.rank(), 0);
    sparse_vector<mpq_class> vector;

    vector.push_back(1, mpq_class(std::move(mpq_class(1))));
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

#endif /* IPO_WITH_GMP */
}
