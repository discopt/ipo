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

  void scaleIntegral(Constraint<double>& constraint)
  {

  }

#if defined(IPO_WITH_GMP)

  void scaleIntegral(Constraint<mpq_class>& constraint)
  {
    int positiveNegative = 0;
    IntegralScaler scaler;
    for (const auto& iter : constraint.vector())
    {
      scaler(iter.second);
      if (sgn(iter.second) > 0)
        ++positiveNegative;
      else
        --positiveNegative;
    }

    if (scaler.factor() == 1 && positiveNegative >= 0)
      return;

    mpq_class factor = scaler.factor();
    if (positiveNegative < 0)
      factor *= -1;

    // Change types and swap lhs/rhs if we negate.
    if (sgn(factor) < 0 && constraint._type != EQUATION)
    {
      std::swap(constraint._lhs, constraint._rhs);
      if (constraint._type == LESS_OR_EQUAL)
        constraint._type = GREATER_OR_EQUAL;
      else if (constraint._type == GREATER_OR_EQUAL)
        constraint._type = LESS_OR_EQUAL;
    }

    // Scale numbers.
    constraint._lhs *= factor;
    constraint._rhs *= factor;
    sparse_vector<mpq_class> vector;
    for (const auto& iter : constraint.vector())
      vector.push_back(iter.first, iter.second * factor);
    constraint._vector = std::make_shared<sparse_vector<mpq_class>>(std::move(vector));
  }

#endif /* IPO_WITH_GMP */

}
