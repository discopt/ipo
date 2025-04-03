#include <gtest/gtest.h>

#include <ipo/lp.hpp>

#if defined(IPO_WITH_DOUBLE_LP)

TEST(LP, OptimizationDouble)
{
  ipo::LP<double> lp;
}

#endif /* IPO_WITH_DOUBLE_LP */

#if defined(IPO_WITH_RATIONAL_LP)

TEST(LP, OptimizationRational)
{
  ipo::LP<ipo::rational> lp;
}

#endif /* IPO_WITH_RATIONAL_LP */
