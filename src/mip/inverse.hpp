#pragma once

#include <iostream>
#include <sstream>
#include <random>

#include <ipo/oracles.hpp>
#include <ipo/affine_hull.hpp>
#include <ipo/lp.hpp>

namespace inverse
{

  template <typename Number>
  void run(std::shared_ptr<ipo::OptimizationOracle<Number>> oracle)
  {
    auto poly = std::make_shared<ipo::Polyhedron<Number>>(oracle);

    ipo::AffineHull<Number> affineHull;
    std::cerr << "Starting affine hull computation.\n" << std::flush;
    ipo::AffineHullQuery affQuery;
    // affQuery.timeLimit = timeLimit;
    affineHull = ipo::affineHull(poly, affQuery);
    std::cout << "Dimension: " << affineHull.dimension << " / " << poly->space()->dimension() << std::endl;

    assert(false);
  }

  int printUsage(const std::string& program)
  {
    std::cout << program << " [OPTIONS] FILE...\n";
    std::cout << "Solves inverse optimization problems on a polyhedron defined by FILE.\n";
    std::cout << "General options:\n";
    std::cout << " -h       Show this help and exit.\n";
    std::cout << "Oracle/polyhedron options:\n";
#if defined(IPO_RATIONAL)
    std::cout << " -x       Use exact arithmetic oracles instead of double precision.\n";
#endif /* IPO_RATIONAL */
    std::cout << std::flush;

    return EXIT_FAILURE;
  }

} /* namespace inverse */

