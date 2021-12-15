#include <gtest/gtest.h>
#include <ipo/oracles_forest.hpp>
#include <ipo/oracles_polar.hpp>
#include <ipo/affine_hull.hpp>
#include <ipo/lp.hpp>

#if defined(IPO_RATIONAL_LP)

TEST(Polar, OptimizationRational)
{
  typedef std::pair<std::size_t, std::size_t> E;

  size_t numNodes = 4;
  size_t numEdges = 6;
  E edges[] = {
    E(0,1), E(0,2), E(0,3), E(1,2), E(1,3), E(2,3)
  };

  auto forestOracle = std::make_shared<ipo::ForestOptimizationOracle<ipo::rational>>(numNodes, &edges[0], &edges[numEdges], false);
  auto forestPolyhedron = std::make_shared<ipo::Polyhedron<ipo::rational>>(forestOracle);

  auto affineHull = ipo::affineHull(forestPolyhedron);
  ASSERT_EQ(affineHull.dimension, 6);

  ipo::LP<ipo::rational> lp;
  lp.setSense(ipo::LPSense::MAXIMIZE);
  int obj = 1;
  for (std::size_t c = 0; c < forestPolyhedron->space()->dimension(); ++c)
  {
    std::ostringstream name;
    name << "x#" << c;
    lp.addColumn(0, 1, ipo::rational(9+obj) / ipo::rational(10), name.str());
    ++obj;
  }
  lp.update();
  
  std::vector<int> nonzeroColumns;
  std::vector<ipo::rational> nonzeroCoefficients;

  auto polarOracle = std::make_shared<ipo::PolarSeparationOracle<ipo::rational>>(forestPolyhedron);
  polarOracle->setAffineHull(affineHull);

  while (true)
  {
    auto status = lp.solve();
    std::cout << "Main LP status: " << status << "." << std::endl;
    if (status == ipo::LPStatus::OPTIMAL)
    {
      std::cout << "The optimum is " << ipo::convertNumber<double>(lp.getObjectiveValue()) << "="
        << lp.getObjectiveValue() << "." << std::endl;
      assert(lp.hasPrimalSolution());
      std::vector<ipo::rational> solution = lp.getPrimalSolution();
      for (std::size_t c = 0; c < forestPolyhedron->space()->dimension(); ++c)
        std::cout << (c == 0 ? "Solution vector: " : ", ") << forestPolyhedron->space()->variable(c) << "=" << solution[c];
      std::cout << std::endl;

      auto sepaResponse = polarOracle->separate(&solution[0], true);
      std::cout << "Found " << sepaResponse.constraints.size() << " undominated constraints:" << std::endl;
      for (const auto& cons : sepaResponse.constraints)
      {
        std::cout << "Adding constraint ";
        forestPolyhedron->space()->printConstraint(std::cout, cons);
        std::cout << std::endl;
        nonzeroColumns.clear();
        nonzeroCoefficients.clear();
        for (const auto& iter : cons.vector())
        {
          nonzeroColumns.push_back(iter.first);
          nonzeroCoefficients.push_back(iter.second);
        }
        lp.addRow(cons.lhs(), nonzeroColumns.size(), &nonzeroColumns[0], &nonzeroCoefficients[0], cons.rhs());
      }
      if (sepaResponse.constraints.empty())
        break;
    }
    else
    {
      assert(false);
    }
  }
}

#endif /* IPO_RATIONAL_LP */
