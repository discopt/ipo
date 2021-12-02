#include <gtest/gtest.h>
#include <ipo/oracles_forest.hpp>

#include <ipo/affine_hull.hpp>

TEST(Forest, SpanningTreeReal)
{
  typedef std::pair<std::size_t, std::size_t> E;

  size_t numNodes = 7;
  size_t numEdges = 9;
  E edges[] = {
    E(0,1), E(0,2), E(0,3), E(0,4), E(1,2), E(2,3), E(3,4), E(4,1), E(5,6)
  };

  auto forestOracle = std::make_shared<ipo::ForestRealOptimizationOracle>(numNodes, &edges[0], &edges[numEdges], false);
  auto forestPolytope = std::make_shared<ipo::RealPolyhedron>(forestOracle);
  ipo::RealAffineHullResult result = ipo::affineHull(forestPolytope);
  ASSERT_EQ(result.dimension, 9);

  auto sptreeOracle = std::make_shared<ipo::ForestRealOptimizationOracle>(numNodes, &edges[0], &edges[numEdges], true);
  auto sptreePolytope = std::make_shared<ipo::RealPolyhedron>(sptreeOracle);
  result = ipo::affineHull(sptreePolytope);
  ASSERT_EQ(result.dimension, 7);
}

TEST(Forest, SpanningTreeRational)
{
  typedef std::pair<std::size_t, std::size_t> E;

  size_t numNodes = 7;
  size_t numEdges = 9;
  E edges[] = {
    E(0,1), E(0,2), E(0,3), E(0,4), E(1,2), E(2,3), E(3,4), E(4,1), E(5,6)
  };

  auto forestOracle = std::make_shared<ipo::ForestRationalOptimizationOracle>(numNodes, &edges[0], &edges[numEdges], false);
  auto forestPolytope = std::make_shared<ipo::RationalPolyhedron>(forestOracle);
  auto result = ipo::affineHull(forestPolytope);
  ASSERT_EQ(result.dimension, 9);

  auto sptreeOracle = std::make_shared<ipo::ForestRationalOptimizationOracle>(numNodes, &edges[0], &edges[numEdges], true);
  auto sptreePolytope = std::make_shared<ipo::RationalPolyhedron>(sptreeOracle);
  result = ipo::affineHull(sptreePolytope);
  ASSERT_EQ(result.dimension, 7);
}
