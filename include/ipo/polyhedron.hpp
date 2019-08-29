#pragma once

#include <ipo/config.hpp>
#include <ipo/export.hpp>
#include <ipo/space.hpp>
#include <ipo/oracles.hpp>

#include <memory>

namespace ipo
{

  class Polyhedron : public std::enable_shared_from_this<Polyhedron>
  {
  public:
    Polyhedron(std::shared_ptr<OptimizationOracle> optimizationOracle);

    ~Polyhedron();

    

  protected:
    std::shared_ptr<OptimizationOracle> _optimizationOracle;
  };

} /* namespace ipo */

