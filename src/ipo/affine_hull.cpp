#include <ipo/ipo.hpp>

#include "redundancy.hpp"
#include "affine_complement.hpp"

#include <iostream>

namespace ipo
{

//   int affineHull(std::shared_ptr<Polyhedron> polyhedron, std::vector<Vector>& innerPoints,
//     std::vector<Vector>& innerRays, std::vector<Constraint>& outerEquations,
//     const std::vector<Constraint>& knownEquations, double timeLimit, AffineHullOptions options)
//   {
//     std::size_t n = polyhedron->space()->dimension();
// 
// #if defined(IPO_WITH_GMP)
//     if (options == AFFINE_HULL_RATIONAL)
//     {
//       auto redundancyCheck = EquationRedundancyCheck<mpq_class, ipo::RationalIsZero>(n);
//       for (auto equation : knownEquations)
//       {
//         if (redundancyCheck.add(equation) == EQUATION_INCONSISTENT)
//         {
//           // TODO: By inspecting the multipliers one could find a smaller inconsistent set.
// 
//           outerEquations.clear();
//           for (std::size_t e = 0; e < redundancyCheck.rank(); ++e)
//             outerEquations.push_back(redundancyCheck.getEquation(e));
//           outerEquations.push_back(equation);
//           return -1;
//         }
//       }
// 
//       std::cout << "Initial " << knownEquations.size() << " equations have rank " << redundancyCheck.rank() << "." << std::endl;
//     }
// #endif /* IPO_WITH_GMP */
//   }

  template <typename T, typename IsZero>
  int affineHullImplementation(std::shared_ptr<Polyhedron<T, IsZero>> polyhedron,
    std::vector<sparse_vector<T>>& innerPoints,
    std::vector<sparse_vector<T>>& innerRays,
    std::vector<Constraint<T>>& outerEquations,
    const std::vector<Constraint<T>>& knownEquations, double timeLimit, IsZero isZero)
  {
    auto redundancyCheck = EquationRedundancyCheck<T, IsZero>(polyhedron->space()->dimension(),
      isZero);
    for (auto equation : knownEquations)
    {
      if (redundancyCheck.add(equation) == EQUATION_INCONSISTENT)
      {
        // TODO: By inspecting the multipliers one could in principle find a smaller inconsistent set.

        outerEquations.clear();
        for (std::size_t e = 0; e < redundancyCheck.rank(); ++e)
          outerEquations.push_back(redundancyCheck.getEquation(e));
        outerEquations.push_back(equation);
        return -1;
      }
    }

    std::cout << "Initial " << knownEquations.size() << " equations have rank " << redundancyCheck.rank() << "." << std::endl;
  }
  
  int affineHull(std::shared_ptr<Polyhedron<double, DoubleIsZero>> polyhedron,
    std::vector<sparse_vector<double>>& innerPoints,
    std::vector<sparse_vector<double>>& innerRays,
    std::vector<Constraint<double>>& outerEquations,
    const std::vector<Constraint<double>>& knownEquations, double timeLimit)
  {
    return affineHullImplementation(polyhedron, innerPoints, innerRays, outerEquations,
      knownEquations, timeLimit, DoubleIsZero(1.0e-9));
  }

#if defined(IPO_WITH_GMP)

  int affineHull(std::shared_ptr<Polyhedron<rational, RationalIsZero>> polyhedron,
    std::vector<sparse_vector<rational>>& innerPoints,
    std::vector<sparse_vector<rational>>& innerRays,
    std::vector<Constraint<rational>>& outerEquations,
    const std::vector<Constraint<rational>>& knownEquations, double timeLimit)
  {
    return affineHullImplementation(polyhedron, innerPoints, innerRays, outerEquations,
      knownEquations, timeLimit, RationalIsZero());
  }

#endif /* IPO_WITH_GMP */

} /* namespace ipo */
