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
  auto projection = std::make_shared<ipo::Projection<double>>(sptreeOracle->space(), "x_[1-9].*");  
  auto projectedOracle = std::make_shared<ipo::ProjectionOptimizationOracle<double>>(sptreeOracle, projection);
  ASSERT_EQ(projectedOracle->space()->dimension(), 6);
  auto projectedPolytope = std::make_shared<ipo::Polyhedron<double>>(projectedOracle);

  auto result = ipo::affineHull(projectedPolytope);
  ASSERT_EQ(result.dimension, 6);
}

#if defined(IPO_WITH_GMP)

TEST(Projection, OptimizationRational)
{
  typedef std::pair<std::size_t, std::size_t> E;

  size_t numNodes = 5;
  size_t numEdges = 10;
  E edges[] = {
    E(0,1), E(0,2), E(0,3), E(0,4), E(1,2), E(1,3), E(1,4), E(2,3), E(2,4), E(3,4)
  };

  auto sptreeOracle = std::make_shared<ipo::ForestRationalOptimizationOracle>(numNodes, &edges[0], &edges[numEdges], true);
  auto projection = std::make_shared<ipo::Projection<mpq_class>>(sptreeOracle->space(), "x_[1-9].*");  
  auto projectedOracle = std::make_shared<ipo::ProjectionOptimizationOracle<mpq_class>>(sptreeOracle, projection);
  ASSERT_EQ(projectedOracle->space()->dimension(), 6);
  auto projectedPolytope = std::make_shared<ipo::Polyhedron<mpq_class>>(projectedOracle);

  auto result = ipo::affineHull(projectedPolytope);
  ASSERT_EQ(result.dimension, 6);
}

#endif /* IPO_WITH_GMP */
