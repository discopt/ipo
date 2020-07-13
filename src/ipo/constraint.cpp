#include <ipo/constraint.hpp>

// #include <iostream> // TODO: Debug

#include "reconstruct.hpp"

namespace ipo
{
  std::ostream& operator<<(std::ostream& stream, const Constraint<double>& constraint)
  {
    if (constraint.lhs() == -std::numeric_limits<double>::infinity())
      return stream << constraint.vector() << " <= " << constraint.rhs();
    else if (constraint.rhs() == std::numeric_limits<double>::infinity())
      return stream << constraint.vector() << " >= " << constraint.lhs();
    else if (constraint.isEquation())
      return stream << constraint.vector() << " == " << constraint.rhs();
    else
      return stream << constraint.lhs() << " <= " << constraint.vector() << " <= " << constraint.rhs();
  }

#if defined(IPO_WITH_GMP)

  std::ostream& operator<<(std::ostream& stream, const Constraint<rational>& constraint)
  {
    if (constraint.lhs() == -std::numeric_limits<double>::infinity())
      return stream << constraint.vector() << " <= " << constraint.rhs();
    else if (constraint.rhs() == std::numeric_limits<double>::infinity())
      return stream << constraint.vector() << " >= " << constraint.lhs();
    else if (constraint.isEquation())
      return stream << constraint.vector() << " == " << constraint.rhs();
    else
      return stream << constraint.lhs() << " <= " << constraint.vector() << " <= " << constraint.rhs();
  }

  Constraint<double> constraintToDouble(const Constraint<rational>& constraint)
  {
    sparse_vector<double> vector;
    for (const auto& iter : constraint.vector())
      vector.push_back(iter.first, iter.second.approximation());
    return Constraint<double>(constraint.lhs().approximation(), std::move(vector),
      constraint.rhs().approximation());
  }

  Constraint<rational> constraintToRational(const Constraint<double>& constraint)
  {
    sparse_vector<rational> vector;
    for (const auto& iter : constraint.vector())
      vector.push_back(iter.first, rational(reconstruct(iter.second)));
    rational lhs, rhs;
    if (isMinusInfinity(constraint.lhs()))
      lhs = minusInfinity();
    else
      lhs = reconstruct(constraint.lhs());
    if (isPlusInfinity(constraint.rhs()))
      rhs = plusInfinity();
    else
      rhs = reconstruct(constraint.rhs());
    return Constraint<rational>(lhs, std::move(vector), rhs);
  }

#endif /* IPO_WITH_GMP */


}
