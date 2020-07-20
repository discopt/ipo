#include <ipo/constraint.hpp>

#include "reconstruct.hpp"

namespace ipo
{
  std::ostream& operator<<(std::ostream& stream, const Constraint<double>& constraint)
  {
    switch (constraint.type())
    {
    case EQUATION:
      return stream << constraint.vector() << " == " << constraint.rhs();
    case LESS_OR_EQUAL:
      return stream << constraint.vector() << " <= " << constraint.rhs();
    case GREATER_OR_EQUAL:
      return stream << constraint.vector() << " >= " << constraint.lhs();
    case RANGED:
      return stream << constraint.lhs() << " <= " << constraint.vector() << " <= "
        << constraint.rhs();
    default:
      return stream << "<Invalid constraint type>";
    }
  }

#if defined(IPO_WITH_GMP)

  std::ostream& operator<<(std::ostream& stream, const Constraint<rational>& constraint)
  {
    switch (constraint.type())
    {
    case EQUATION:
      return stream << constraint.vector() << " == " << constraint.rhs();
    case LESS_OR_EQUAL:
      return stream << constraint.vector() << " <= " << constraint.rhs();
    case GREATER_OR_EQUAL:
      return stream << constraint.vector() << " >= " << constraint.lhs();
    case RANGED:
      return stream << constraint.lhs() << " <= " << constraint.vector() << " <= "
        << constraint.rhs();
    default:
      return stream << "<Invalid constraint type>";
    }
  }

  Constraint<double> constraintToDouble(const Constraint<rational>& constraint)
  {
    auto vector = std::make_shared<sparse_vector<double>>();
    for (const auto& iter : constraint.vector())
      vector->push_back(iter.first, iter.second.approximation());
    return Constraint<double>(constraint.lhs().approximation(), vector,
      constraint.rhs().approximation(), constraint.type());
  }

  Constraint<rational> constraintToRational(const Constraint<double>& constraint)
  {
    auto vector = std::make_shared<sparse_vector<rational>>();
    for (const auto& iter : constraint.vector())
      vector->push_back(iter.first, rational(reconstruct(iter.second)));
    return Constraint<rational>(reconstruct(constraint.lhs()), vector, constraint.rhs(),
      constraint.type());
  }

#endif /* IPO_WITH_GMP */


}
