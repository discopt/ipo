#include <ipo/polyhedron.hpp>

namespace ipo
{

  Polyhedron::Polyhedron(std::shared_ptr<OptimizationOracle> optimizationOracle)
    : _optimizationOracle(optimizationOracle)
  {

  }

  Polyhedron::~Polyhedron()
  {

  }

} /* namespace ipo */
