#include <gtest/gtest.h>
#include <ipo/oracles_forest.hpp>
#include <ipo/dominant.hpp>
#include <ipo/affine_hull.hpp>

TEST(Dominant, OptimizationDouble)
{
  typedef std::pair<std::size_t, std::size_t> E;

  size_t numNodes = 5;
  size_t numEdges = 10;
  E edges[] = {
    E(0,1), E(0,2), E(0,3), E(0,4), E(1,2), E(1,3), E(1,4), E(2,3), E(2,4), E(3,4)
  };

  auto sptreeOracle = std::make_shared<ipo::ForestOptimizationOracle<double>>(numNodes, &edges[0], &edges[numEdges], true);
  auto dominantOracle = std::make_shared<ipo::DominantOptimizationOracle<double>>(sptreeOracle);
  auto poly = std::make_shared<ipo::Polyhedron<double>>(dominantOracle);

  auto result = ipo::affineHull(poly);
  ASSERT_EQ(result.dimension, 10);
}

#if defined(IPO_WITH_GMP)

TEST(Dominant, OptimizationRational)
{
  typedef std::pair<std::size_t, std::size_t> E;

  size_t numNodes = 5;
  size_t numEdges = 10;
  E edges[] = {
    E(0,1), E(0,2), E(0,3), E(0,4), E(1,2), E(1,3), E(1,4), E(2,3), E(2,4), E(3,4)
  };

  auto sptreeOracle = std::make_shared<ipo::ForestOptimizationOracle<mpq_class>>(numNodes, &edges[0], &edges[numEdges], true);
  auto dominantOracle = std::make_shared<ipo::DominantOptimizationOracle<mpq_class>>(sptreeOracle);
  auto poly = std::make_shared<ipo::Polyhedron<mpq_class>>(dominantOracle);

  auto result = ipo::affineHull(poly);
  ASSERT_EQ(result.dimension, 10);
}

#endif /* IPO_WITH_GMP */
