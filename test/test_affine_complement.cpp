#include <gtest/gtest.h>
#include <../src/ipo/affine_complement.hpp>

TEST(AffineComplement, Real)
{
  auto ac = ipo::AffineComplement<double, ipo::RealIsZero>(4, ipo::RealIsZero(1.0e-9));
  ASSERT_EQ(ac.rank(), 0);

  std::vector<double> vector(4);
  vector[0] = 0.0;
  vector[1] = 1.0;
  vector[2] = 0.0;
  vector[3] = 0.0;
  ac.add(ipo::Vector(vector), 1.0, 1);
  ASSERT_EQ(ac.rank(), 1);

  vector[0] = 1.0;
  vector[1] = 1.0;
  vector[2] = 0.0;
  vector[3] = 0.0;
  ac.add(ipo::Vector(vector), 1.0, 0);
  ASSERT_EQ(ac.rank(), 2);

  vector[0] = 0.0;
  vector[1] = 0.0;
  vector[2] = 0.0;
  vector[3] = 1.0;
  ac.add(ipo::Vector(vector), 1.0, 4);
  ASSERT_EQ(ac.rank(), 3);
}


TEST(AffineComplement, Rational)
{
  auto ac = ipo::AffineComplement<mpq_class, ipo::RationalIsZero>(4, ipo::RationalIsZero());
  ASSERT_EQ(ac.rank(), 0);

  std::vector<mpq_class> vector(4);
  vector[0] = 0;
  vector[1] = 1;
  vector[2] = 0;
  vector[3] = 0;
  auto row = ipo::Vector(vector);
  row.checkConsistency();
  ac.add(ipo::Vector(vector), 1, 1);
  ASSERT_EQ(ac.rank(), 1);

  vector[0] = 1;
  vector[1] = 1;
  vector[2] = 0;
  vector[3] = 0;
  ac.add(ipo::Vector(vector), 1, 0);
  ASSERT_EQ(ac.rank(), 2);

  vector[0] = 0;
  vector[1] = 0;
  vector[2] = 0;
  vector[3] = 1;
  ac.add(ipo::Vector(vector), 1, 4);
  ASSERT_EQ(ac.rank(), 3);
}
