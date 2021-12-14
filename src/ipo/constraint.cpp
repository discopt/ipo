#include <ipo/constraint.hpp>

#include <ipo/arithmetic.hpp>

namespace ipo
{

#if defined(IPO_DOUBLE)

  std::ostream& operator<<(std::ostream& stream, const Constraint<double>& constraint)
  {
    switch (constraint.type())
    {
    case ConstraintType::EQUATION:
      return stream << constraint.vector() << " == " << constraint.rhs();
    case ConstraintType::LESS_OR_EQUAL:
      return stream << constraint.vector() << " <= " << constraint.rhs();
    case ConstraintType::GREATER_OR_EQUAL:
      return stream << constraint.vector() << " >= " << constraint.lhs();
    case ConstraintType::RANGED:
      return stream << constraint.lhs() << " <= " << constraint.vector() << " <= "
        << constraint.rhs();
    default:
      return stream << "<Invalid constraint type>";
    }
  }

#endif /* IPO_DOUBLE */

#if defined(IPO_RATIONAL)

  std::ostream& operator<<(std::ostream& stream, const Constraint<rational>& constraint)
  {
    switch (constraint.type())
    {
    case ConstraintType::EQUATION:
      return stream << constraint.vector() << " == " << constraint.rhs();
    case ConstraintType::LESS_OR_EQUAL:
      return stream << constraint.vector() << " <= " << constraint.rhs();
    case ConstraintType::GREATER_OR_EQUAL:
      return stream << constraint.vector() << " >= " << constraint.lhs();
    case ConstraintType::RANGED:
      return stream << constraint.lhs() << " <= " << constraint.vector() << " <= "
        << constraint.rhs();
    default:
      return stream << "<Invalid constraint type>";
    }
  }

#endif /* IPO_RATIONAL */

  void scaleIntegral(Constraint<double>& constraint)
  {

  }

#if defined(IPO_RATIONAL)

  void scaleIntegral(Constraint<rational>& constraint)
  {
    int positiveNegative = 0;
    IntegralScaler scaler;
    for (const auto& iter : constraint.vector())
    {
      scaler(iter.second);
      if (iter.second.sign() > 0)
        ++positiveNegative;
      else
        --positiveNegative;
    }

    if (scaler.factor() == 1 && positiveNegative >= 0)
      return;

    rational factor = scaler.factor();
    if (positiveNegative < 0)
      factor *= -1;

    // Change types and swap lhs/rhs if we negate.
    if (factor.sign() < 0 && constraint._type != ConstraintType::EQUATION)
    {
      std::swap(constraint._lhs, constraint._rhs);
      if (constraint._type == ConstraintType::LESS_OR_EQUAL)
        constraint._type = ConstraintType::GREATER_OR_EQUAL;
      else if (constraint._type == ConstraintType::GREATER_OR_EQUAL)
        constraint._type = ConstraintType::LESS_OR_EQUAL;
    }

    // Scale numbers.
    constraint._lhs *= factor;
    constraint._rhs *= factor;
    sparse_vector<rational> vector;
    for (const auto& iter : constraint.vector())
      vector.push_back(iter.first, iter.second * factor);
    constraint._vector = std::make_shared<sparse_vector<rational>>(std::move(vector));
  }

#endif /* IPO_RATIONAL */

}
