#include <gtest/gtest.h>
#include <ipo/space.hpp>

TEST(LinearAlgebra, Space)
{
  std::cout << "===== LinearAlgebra::Space ===== " << std::endl;
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

    sparse_vector<double> vector;
    ASSERT_EQ(space1.printVector(vector), "()");
    ASSERT_EQ(space1.printLinearForm(vector), "0");

    vector.push_back(1, -2.0);
    vector.push_back(2, 1.0);
    ASSERT_EQ(space1.printVector(vector), "(y=-2,z=1)");
    ASSERT_EQ(space1.printLinearForm(vector), "-2<y> + <z>");

    
  }
}
