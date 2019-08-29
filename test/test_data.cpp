#include <gtest/gtest.h>
#include <ipo/data.hpp>

TEST(Data, Value)
{
  auto value = ipo::Value(1.0);
  ASSERT_TRUE(value.isFinite());
  ASSERT_FALSE(value.isPlusInfinity());
  ASSERT_FALSE(value.isMinusInfinity());
  value = ipo::Value(std::numeric_limits<double>::infinity());
  ASSERT_FALSE(value.isFinite());
  ASSERT_TRUE(value.isPlusInfinity());
  ASSERT_FALSE(value.isMinusInfinity());
  value = ipo::Value(-std::numeric_limits<double>::infinity());
  ASSERT_FALSE(value.isFinite());
  ASSERT_FALSE(value.isPlusInfinity());
  ASSERT_TRUE(value.isMinusInfinity());

#if defined(IPO_WITH_GMP)
  value = ipo::Value(mpq_class(10,3));
  ASSERT_TRUE(value.isFinite());
  ASSERT_FALSE(value.isPlusInfinity());
  ASSERT_FALSE(value.isMinusInfinity());
#endif /* IPO_WITH_GMP */
}

TEST(Data, Vector)
{
  {
    auto vec = ipo::Vector();
    ASSERT_EQ(vec.usage(), 1);
    {
      auto copy = vec;
      ASSERT_EQ(vec.usage(), 2);
      ASSERT_EQ(copy.usage(), 2);

      vec = copy;
    }
    ASSERT_EQ(vec.usage(), 1);
    ASSERT_EQ(vec.size(), 0);

    vec = ipo::Vector(ipo::Vector());
    auto vec2 = ipo::Vector();

    std::vector<std::pair<std::size_t, double>> entries;
    entries.push_back(std::pair<std::size_t, double>(0, 1.0));
    entries.push_back(std::pair<std::size_t, double>(2, 3.0));
    vec = ipo::Vector(entries);
    ASSERT_EQ(vec.size(), 2);
    ASSERT_EQ(vec.coordinate(0), 0);
    ASSERT_EQ(vec.coordinate(1), 2);
    ASSERT_EQ(vec.real(0), 1.0);
    ASSERT_EQ(vec.real(1), 3.0);
  }

#if defined(IPO_WITH_GMP)
  {
    std::vector<std::pair<std::size_t, mpq_class>> entries;
    entries.push_back(std::pair<std::size_t, mpq_class>(0, mpq_class(1,2)));
    entries.push_back(std::pair<std::size_t, mpq_class>(2, mpq_class(-2,1)));
    auto vec = ipo::Vector(entries);
    ASSERT_TRUE(vec.isRational());
    ASSERT_EQ(vec.size(), 2);
    ASSERT_EQ(vec.coordinate(0), 0);
    ASSERT_EQ(vec.coordinate(1), 2);
    ASSERT_EQ(vec.rational(0), mpq_class(1,2));
    ASSERT_EQ(vec.rational(1), mpq_class(-2,1));
    ASSERT_EQ(vec.real(0), 0.5);
    ASSERT_EQ(vec.real(1), -2.0);
    ASSERT_EQ(entries[0].second, mpq_class(1,2));

    vec = ipo::Vector(entries);
    ASSERT_TRUE(vec.isRational());
    ASSERT_EQ(entries[0].second, mpq_class(1,2));
  }
#endif /* IPO_WITH_GMP */
}

TEST(Data, Product)
{
  {
    std::vector<std::pair<std::size_t, double>> a;
    std::vector<std::pair<std::size_t, double>> b;
    for (std::size_t i = 0; i < 100; ++i)
    {
      a.push_back(std::pair<std::size_t, double>(i, 0.1));
      b.push_back(std::pair<std::size_t, double>(i, 0.1));
    }

    ipo::Value product = ipo::Vector(a) * ipo::Vector(b);
    ASSERT_NE(product.real, 1.0);
    ASSERT_DOUBLE_EQ(product.real, 1.0);
  }

#if defined(IPO_WITH_GMP)
  {
    std::vector<std::pair<std::size_t, mpq_class>> a;
    std::vector<std::pair<std::size_t, mpq_class>> b;
    for (std::size_t i = 0; i < 100; ++i)
    {
      a.push_back(std::pair<std::size_t, mpq_class>(i, mpq_class(1,10)));
      b.push_back(std::pair<std::size_t, mpq_class>(i, mpq_class(1,10)));
    }

    ipo::Value product = ipo::Vector(a) * ipo::Vector(b);
    ASSERT_EQ(product.rational, 1);
    ASSERT_EQ(product.real, 1.0);
    ASSERT_DOUBLE_EQ(product.real, 1.0);
  }
#endif /* IPO_WITH_GMP */
}
