#include <gtest/gtest.h>
#include <ipo/oracles_forest.hpp>
#include <ipo/oracles_projection.hpp>
#include <ipo/affine_hull.hpp>

TEST(Projection, OptimizationReal)
{
  typedef std::pair<std::size_t, std::size_t> E;

  size_t numNodes = 5;
  size_t numEdges = 10;
  E edges[] = {
    E(0,1), E(0,2), E(0,3), E(0,4), E(1,2), E(1,3), E(1,4), E(2,3), E(2,4), E(3,4)
  };

  auto sptreeOracle = std::make_shared<ipo::ForestRealOptimizationOracle>(numNodes, &edges[0], &edges[numEdges], true);
  auto projectedOracle = std::make_shared<ipo::ProjectionRealOptimizationOracle>(sptreeOracle);
  ASSERT_EQ(projectedOracle->addVariables("x_[1-9].*"), 6);
  auto projectedPolytope = std::make_shared<ipo::RealPolyhedron>(projectedOracle);

  auto result = ipo::affineHull(projectedPolytope);
  ASSERT_EQ(result.dimension, 6);
}

TEST(Projection, OptimizationRational)
{
  typedef std::pair<std::size_t, std::size_t> E;

  size_t numNodes = 5;
  size_t numEdges = 10;
  E edges[] = {
    E(0,1), E(0,2), E(0,3), E(0,4), E(1,2), E(1,3), E(1,4), E(2,3), E(2,4), E(3,4)
  };

  auto sptreeOracle = std::make_shared<ipo::ForestRationalOptimizationOracle>(numNodes, &edges[0], &edges[numEdges], true);
  auto projectedOracle = std::make_shared<ipo::ProjectionRationalOptimizationOracle>(sptreeOracle);
  ASSERT_EQ(projectedOracle->addVariables("x_[1-9].*"), 6);
  auto projectedPolytope = std::make_shared<ipo::RationalPolyhedron>(projectedOracle);

  auto result = ipo::affineHull(projectedPolytope);
  ASSERT_EQ(result.dimension, 6);
}
