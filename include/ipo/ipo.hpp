#pragma once

#include <ipo/config.hpp>
#include <ipo/export.hpp>
#include <ipo/polyhedron.hpp>
#include <ipo/rational.hpp>

namespace ipo
{

  IPO_EXPORT
  int affineHull(std::shared_ptr<Polyhedron<double, DoubleIsZero>> polyhedron,
    std::vector<sparse_vector<double>>& innerPoints,
    std::vector<sparse_vector<double>>& innerRays,
    std::vector<Constraint<double>>& outerEquations,
    const std::vector<Constraint<double>>& knownEquations = std::vector<Constraint<double>>(),
    double timeLimit = std::numeric_limits<double>::infinity());

#if defined(IPO_WITH_GMP)

  IPO_EXPORT
  int affineHull(std::shared_ptr<Polyhedron<rational, RationalIsZero>> polyhedron,
    std::vector<sparse_vector<rational>>& innerPoints,
    std::vector<sparse_vector<rational>>& innerRays,
    std::vector<Constraint<rational>>& outerEquations,
    const std::vector<Constraint<rational>>& knownEquations = std::vector<Constraint<rational>>(),
    double timeLimit = std::numeric_limits<double>::infinity());

#endif /* IPO_WITH_GMP */

} /* namespace ipo */
