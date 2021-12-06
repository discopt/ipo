#include <gtest/gtest.h>
#include <ipo/oracles_forest.hpp>

#include <ipo/affine_hull.hpp>

TEST(Forest, SpanningTreeDouble)
{
  typedef std::pair<std::size_t, std::size_t> E;

  size_t numNodes = 7;
  size_t numEdges = 9;
  E edges[] = {
    E(0,1), E(0,2), E(0,3), E(0,4), E(1,2), E(2,3), E(3,4), E(4,1), E(5,6)
  };

  auto forestOracle = std::make_shared<ipo::ForestOptimizationOracle<double>>(numNodes, &edges[0], &edges[numEdges], false);
  auto forestPolytope = std::make_shared<ipo::Polyhedron<double>>(forestOracle);
  auto result = ipo::affineHull(forestPolytope);
  ASSERT_EQ(result.dimension, 9);

  auto sptreeOracle = std::make_shared<ipo::ForestOptimizationOracle<double>>(numNodes, &edges[0], &edges[numEdges], true);
  auto sptreePolytope = std::make_shared<ipo::Polyhedron<double>>(sptreeOracle);
  result = ipo::affineHull(sptreePolytope);
  ASSERT_EQ(result.dimension, 7);
}

#if defined(IPO_WITH_GMP)

TEST(Forest, SpanningTreeRational)
{
  typedef std::pair<std::size_t, std::size_t> E;

  size_t numNodes = 7;
  size_t numEdges = 9;
  E edges[] = {
    E(0,1), E(0,2), E(0,3), E(0,4), E(1,2), E(2,3), E(3,4), E(4,1), E(5,6)
  };

  auto forestOracle = std::make_shared<ipo::ForestOptimizationOracle<mpq_class>>(numNodes, &edges[0], &edges[numEdges], false);
  auto forestPolytope = std::make_shared<ipo::Polyhedron<mpq_class>>(forestOracle);
  auto result = ipo::affineHull(forestPolytope);
  ASSERT_EQ(result.dimension, 9);

  auto sptreeOracle = std::make_shared<ipo::ForestOptimizationOracle<mpq_class>>(numNodes, &edges[0], &edges[numEdges], true);
  auto sptreePolytope = std::make_shared<ipo::Polyhedron<mpq_class>>(sptreeOracle);
  result = ipo::affineHull(sptreePolytope);
  ASSERT_EQ(result.dimension, 7);
}

#endif /* IPO_WITH_GMP */
