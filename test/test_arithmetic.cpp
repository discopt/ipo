#include <gtest/gtest.h>

#include <ipo/arithmetic.hpp>

#if defined(IPO_WITH_RATIONAL)

using namespace ipo;

TEST(Rational, Basics)
{
  extended_rational x(0);
  std::cout << x << std::endl;
  extended_rational y(3.14);
  std::cout << y << std::endl;

  extended_rational quiet_nan(std::numeric_limits<extended_rational>::quiet_NaN());
  std::cout << quiet_nan << std::endl;
  extended_rational signaling_nan(std::numeric_limits<extended_rational>::signaling_NaN());
  std::cout << signaling_nan << std::endl;
  extended_rational plus_infinity(std::numeric_limits<extended_rational>::infinity());
  std::cout << plus_infinity << std::endl;
  extended_rational minus_infinity(-std::numeric_limits<extended_rational>::infinity());
  std::cout << minus_infinity << std::endl;

  std::cout << (x + x) << "," << (x+y) << "," << (y+x) << "," << (x+quiet_nan) << "," << (quiet_nan+x) << ","
    << (x + plus_infinity) << "," << (plus_infinity + x) << "," << (x + minus_infinity) << "," << (minus_infinity + x)
    << "," << (plus_infinity + plus_infinity) << "," << (minus_infinity + minus_infinity) << ","
    << (plus_infinity + minus_infinity) << std::endl;

  std::cout << 1.0 / -0.0 << std::endl;
}

#endif /* IPO_WITH_RATIONAL */
