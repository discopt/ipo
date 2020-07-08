#pragma once

#include <ipo/config.hpp>
#include <ipo/export.hpp>

#include <limits>
#include <ostream>
#include <ipo/sparse_vector.hpp>
#include <ipo/rational.hpp>

namespace ipo
{

  template <typename T>
  class Constraint
  {
  public:
    IPO_EXPORT
    Constraint(const T& lhs, const sparse_vector<T>& vector, const T& rhs)
      : _lhs(lhs), _vector(vector), _rhs(rhs)
    {
      assert(!isPlusInfinity(lhs));
      assert(!isMinusInfinity(rhs));
      assert(!(isMinusInfinity(lhs) && isPlusInfinity(rhs)));
#if !defined(NDEBUG)
      for (const auto& iter : vector)
        assert(iter.second != 0);
#endif
    }

    IPO_EXPORT
    Constraint(const T& lhs, sparse_vector<T>&& vector, const T& rhs)
      : _lhs(lhs), _vector(std::move(vector)), _rhs(rhs)
    {

    }

    IPO_EXPORT
    Constraint(T&& lhs, sparse_vector<T>&& vector, T&& rhs)
      : _lhs(std::move(lhs)), _vector(std::move(vector)), _rhs(std::move(rhs))
    {

    }

    IPO_EXPORT
    ~Constraint()
    {

    }

    IPO_EXPORT
    const T& lhs() const
    {
      return _lhs;
    }

    IPO_EXPORT
    const sparse_vector<T>& vector() const
    {
      return _vector;
    }

    IPO_EXPORT
    const T& rhs() const
    {
      return _rhs;
    }

    IPO_EXPORT
    bool isEquation() const
    {
      return _lhs == _rhs;
    }

    IPO_EXPORT
    bool isAlwaysSatisfied() const
    {
      return _vector.empty() && _lhs <= 0 && _rhs >= 0;
    }

    IPO_EXPORT
    bool isNeverSatisfied() const
    {
      return (_vector.empty() && (_lhs > 0 || _rhs < 0)) || (_lhs > _rhs);
    }

  protected:
    T _lhs;
    sparse_vector<T> _vector;
    T _rhs;
  };

  template <typename T>
  Constraint<T> alwaysSatisfiedConstraint()
  {
    return Constraint<T>(T(-std::numeric_limits<double>::infinity()), sparse_vector<T>(), T(1));
  }

  template <typename T>
  Constraint<T> neverSatisfiedConstraint()
  {
    return Constraint<T>(T(-std::numeric_limits<double>::infinity()), sparse_vector<T>(), T(-1));
  }

//   struct HashConstraint
//   {
//     std::size_t operator() (const Constraint& constraint) const
//     {
//       HashVector hash;
//       return hash(constraint.vector);
//     }
//   };  

  IPO_EXPORT
  std::ostream& operator<<(std::ostream& stream, const Constraint<double>& constraint);

#if defined(IPO_WITH_GMP)

  IPO_EXPORT
  std::ostream& operator<<(std::ostream& stream, const Constraint<rational>& constraint);

  IPO_EXPORT
  Constraint<double> constraintToDouble(const Constraint<rational>& constraint);

  IPO_EXPORT
  Constraint<rational> constraintToRational(const Constraint<double>& constraint);

#endif /* IPO_WITH_GMP */

}
