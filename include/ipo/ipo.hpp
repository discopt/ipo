#pragma once

#include <ipo/config.hpp>
#include <ipo/export.hpp>
#include <ipo/polyhedron.hpp>

namespace ipo
{
  typedef int AffineHullOptions;

  const int AFFINE_HULL_REAL = 0;
#if defined(IPO_WITH_GMP)
  const int AFFINE_HULL_RATIONAL = 1;
#endif

  int affineHull(std::shared_ptr<Polyhedron> polyhedron, std::vector<Vector>& innerPoints,
    std::vector<Vector>& innerRays, std::vector<Constraint>& outerEquations,
    const std::vector<Constraint>& knownEquations = std::vector<Constraint>(),
    double timeLimit = std::numeric_limits<double>::infinity(),
#if defined(IPO_WITH_GMP)
    AffineHullOptions options = AFFINE_HULL_RATIONAL
#else
    AffineHullOptions options = AFFINE_HULL_REAL
#endif /* IPO_WITH_GMP */
  );

} /* namespace ipo */
