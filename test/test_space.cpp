#include <gtest/gtest.h>
#include <ipo/space.hpp>

TEST(Space, Basics)
{
   std::vector<std::string> vars;
   vars.push_back("x");
   vars.push_back("y");
   vars.push_back("z");

   ipo::Space space1(vars);
   ipo::Space space2;

   space2.addVariable("x");
   space2.addVariable("y");
   
   ASSERT_EQ(space2.dimension(), 2);
   ASSERT_NE(space1, space2);

   space2.addVariable("z");

   ASSERT_EQ(space1, space2);
   
   std::size_t coords[] = {1, 2};
   double vals[] = { -2.0, 1.0};
   ASSERT_EQ(space1.vectorToString(2, coords, vals), "(y=-2,z=1)");
   ASSERT_EQ(space1.linearFormToString(2, coords, vals), "-2y + z");

   ASSERT_EQ(space1.vectorToString(0, coords, vals), "()");
   ASSERT_EQ(space1.linearFormToString(0, coords, vals), "0");
}
