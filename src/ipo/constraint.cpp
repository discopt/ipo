#include <ipo/constraint.hpp>

#include <ipo/arithmetic.hpp>

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

  std::ostream& operator<<(std::ostream& stream, const Constraint<mpq_class>& constraint)
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

  Constraint<double> constraintToDouble(const Constraint<mpq_class>& constraint)
  {
    auto vector = std::make_shared<sparse_vector<double>>();
    for (const auto& iter : constraint.vector())
      vector->push_back(iter.first, iter.second.get_d());
    return Constraint<double>(constraint.lhs().get_d(), vector,
      constraint.rhs().get_d(), constraint.type());
  }

  Constraint<mpq_class> constraintToRational(const Constraint<double>& constraint)
  {
    auto vector = std::make_shared<sparse_vector<mpq_class>>();
    for (const auto& iter : constraint.vector())
      vector->push_back(iter.first, reconstruct(iter.second));
    return Constraint<mpq_class>(reconstruct(constraint.lhs()), vector, constraint.rhs(),
      constraint.type());
  }

#endif /* IPO_WITH_GMP */


}
