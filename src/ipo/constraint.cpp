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
    int positiveNegative = 0;
    double minAbsolute = std::numeric_limits<double>::infinity();
    for (const auto& iter : constraint.vector())
    {
      double x = iter.second;
      minAbsolute = std::min(minAbsolute, fabs(x));
       if (x > 0)
        ++positiveNegative;
      else
        --positiveNegative;
    }

    if (minAbsolute < 1.0e-9)
      minAbsolute = 1;
    if (minAbsolute == 1 && positiveNegative >= 0)
      return;

    double factor = 1.0 / minAbsolute;
    if (positiveNegative < 0)
      factor *= -1;

    // Change types and swap lhs/rhs if we negate.
    if (factor < 0 && constraint._type != ConstraintType::EQUATION)
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
    sparse_vector<double> vector;
    for (const auto& iter : constraint.vector())
      vector.push_back(iter.first, iter.second * factor);
    constraint._vector = std::make_shared<sparse_vector<double>>(std::move(vector));
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


  template <typename Number>
  bool ConstraintSet<Number>::exists(const Constraint<Number>& constraint, double hashProduct)
  {
    const double maxDistance = 1.0e-3;

    // Compute iterator after the last one whose distance is below maxDistance.

    auto beyondIter = _hashToIndex.lower_bound(hashProduct);
    auto firstIter = beyondIter;
    while (beyondIter != _hashToIndex.end())
    {
      assert(beyondIter->first >= hashProduct);

      if (beyondIter->first - hashProduct > maxDistance)
        break;
      ++beyondIter;
    }
    
    // Compute first iterator for which the distance is below maxDistance.

    while (firstIter != _hashToIndex.begin())
    {
      auto pred = firstIter;
      --pred;
      assert(hashProduct - pred->first >= 0);
      if (hashProduct - pred->first > maxDistance)
        break;
      firstIter = pred;
    }

    if (firstIter == beyondIter)
      return false;

//     auto bestIter = firstIter;
    double bestSquaredDist = std::numeric_limits<double>::infinity();
//     bool bestPositive = true; // whether constraint is close or -constraint is close.

    for (; firstIter != beyondIter; ++firstIter)
    {
      const Constraint<Number>& other = _data[firstIter->second].constraint;
      auto squaredDists = squaredEuclideanDistanceSigned(constraint.vector(), other.vector());
      double lhsDiff = (constraint.hasLhs() && other.hasLhs()) ? constraint.lhs() - other.lhs() : 0.0;
      double lhsSum = (constraint.hasLhs() && other.hasLhs()) ? constraint.lhs() + other.lhs() : 0.0;
      double rhsDiff = (constraint.hasRhs() && other.hasRhs()) ? constraint.rhs() - other.rhs() : 0.0;
      double rhsSum = (constraint.hasRhs() && other.hasRhs()) ? constraint.rhs() + other.rhs() : 0.0;
      squaredDists.first += lhsDiff * lhsDiff + rhsDiff * rhsDiff;
      squaredDists.second += lhsSum * lhsSum + rhsSum * rhsSum;
      if (squaredDists.first < bestSquaredDist)
      {
//         bestIter = firstIter;
        bestSquaredDist = squaredDists.first;
//         bestPositive = true;
      }
      if (squaredDists.second < bestSquaredDist)
      {
//         bestIter = firstIter;
        bestSquaredDist = squaredDists.second;
//         bestPositive = false;
      }
    }

    return bestSquaredDist < 1.0e-12;
  }

  template <>
  ConstraintSet<double>::ConstraintSet(std::size_t ambientDimension, const std::vector<Constraint<double>>& constraints)
    : _n(ambientDimension)
  {
    _hashVector = generateRandomVectorSphere(ambientDimension);
    for (const Constraint<double>& cons : constraints)
      this->push_back(cons);
  }

  template <>
  bool ConstraintSet<double>::exists(const Constraint<double>& constraint)
  {
    double product = 0.0;
    for (const auto& iter : constraint.vector())
      product += iter.second * _hashVector[iter.first];

    return exists(constraint, product);
  }

  template
  bool ConstraintSet<double>::exists(const Constraint<double>& constraint, double hashProduct);
  
#if defined(IPO_RATIONAL)
  template
  bool ConstraintSet<rational>::exists(const Constraint<rational>& constraint, double hashProduct);
#endif /* IPO_RATIONAL */
}
