#pragma once

#include <ipo/config.hpp>
#include <ipo/export.hpp>

#include <ipo/sparse_vector.hpp>

#include <memory>
#include <limits>
#include <ostream>

#if defined(IPO_WITH_GMP)
#include <gmpxx.h>
#endif /* IPO_WITH_GMP */

namespace ipo
{
  enum ConstraintType
  {
    EQUATION,
    LESS_OR_EQUAL,
    GREATER_OR_EQUAL,
    RANGED
  };

  template <typename T>
  class Constraint
  {
  public:
    IPO_EXPORT
    Constraint(const T& lhs, std::shared_ptr<sparse_vector<T>> vector, const T& rhs,
      ConstraintType type)
      : _lhs(lhs), _vector(vector), _rhs(rhs), _type(type)
    {
      assert(_type != RANGED || lhs <= rhs);
      assert(_type != EQUATION || lhs == rhs);
    }

    IPO_EXPORT
    Constraint(T&& lhs, std::shared_ptr<sparse_vector<T>> vector, T&& rhs, ConstraintType type)
      : _lhs(std::move(lhs)), _vector(vector), _rhs(std::move(rhs)), _type(type)
    {
      assert(_type != RANGED || lhs <= rhs);
      assert(_type != EQUATION || lhs == rhs);
    }

    IPO_EXPORT
    Constraint(const T& lhs, std::shared_ptr<sparse_vector<T>> vector, const T& rhs)
      : _lhs(lhs), _vector(vector), _rhs(rhs), _type(lhs == rhs ? EQUATION : RANGED)
    {
      assert(lhs <= rhs);
    }

    IPO_EXPORT
    Constraint(T&& lhs, std::shared_ptr<sparse_vector<T>> vector, T&& rhs)
      : _lhs(std::move(lhs)), _vector(std::move(vector)), _rhs(std::move(rhs)),
      _type(lhs == rhs ? EQUATION : RANGED)
    {
      assert(lhs <= rhs);
    }

    IPO_EXPORT
    Constraint(const T& lhs, std::shared_ptr<sparse_vector<T>> vector)
      : _lhs(lhs), _vector(vector), _rhs(0), _type(GREATER_OR_EQUAL)
    {

    }

    IPO_EXPORT
    Constraint(T&& lhs, std::shared_ptr<sparse_vector<T>> vector)
      : _lhs(std::move(lhs)), _vector(std::move(vector)), _rhs(0),
      _type(GREATER_OR_EQUAL)
    {

    }

    IPO_EXPORT
    Constraint(std::shared_ptr<sparse_vector<T>> vector, const T& rhs)
      : _lhs(0), _vector(vector), _rhs(rhs), _type(LESS_OR_EQUAL)
    {

    }

    IPO_EXPORT
    Constraint(std::shared_ptr<sparse_vector<T>> vector, T&& rhs)
      : _lhs(0), _vector(std::move(vector)), _rhs(std::move(rhs)), _type(LESS_OR_EQUAL)
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
      return *_vector;
    }

    IPO_EXPORT
    std::shared_ptr<sparse_vector<T>> sharedVector()
    {
      return _vector;
    }

    IPO_EXPORT
    ConstraintType type() const
    {
      return _type;
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
    bool hasLhs() const
    {
      return _type != LESS_OR_EQUAL;
    }

    IPO_EXPORT
    bool hasRhs() const
    {
      return _type != GREATER_OR_EQUAL;
    }

    IPO_EXPORT
    bool isAlwaysSatisfied() const
    {
      return _vector->empty() && _lhs <= 0 && _rhs >= 0;
    }

    IPO_EXPORT
    bool isNeverSatisfied() const
    {
      return (_vector->empty() && (_lhs > 0 || _rhs < 0)) || (_lhs > _rhs);
    }

  protected:
    T _lhs;
    std::shared_ptr<sparse_vector<T>> _vector;
    T _rhs;
    ConstraintType _type;
  };

  template <typename T>
  Constraint<T> alwaysSatisfiedConstraint()
  {
    auto zero = std::make_shared<sparse_vector<T>>();
    return Constraint<T>(zero, T(1));
  }

  template <typename T>
  Constraint<T> neverSatisfiedConstraint()
  {
    auto zero = std::make_shared<sparse_vector<T>>();
    return Constraint<T>(zero, T(-1));
  }

  IPO_EXPORT
  std::ostream& operator<<(std::ostream& stream, const Constraint<double>& constraint);

#if defined(IPO_WITH_GMP)

  IPO_EXPORT
  std::ostream& operator<<(std::ostream& stream, const Constraint<mpq_class>& constraint);

  IPO_EXPORT
  Constraint<double> constraintToDouble(const Constraint<mpq_class>& constraint);

  IPO_EXPORT
  Constraint<mpq_class> constraintToRational(const Constraint<double>& constraint);

#endif /* IPO_WITH_GMP */

}
