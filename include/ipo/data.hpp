#pragma once

#include <ipo/config.hpp>
#include <ipo/export.hpp>

#include <cassert>
#include <vector>

#if defined(IPO_WITH_GMP)
#include <gmpxx.h>
#endif /* IPO_WITH_GMP */

namespace ipo
{
  struct Value
  {
#ifdef IPO_WITH_GMP
    mpq_class rational;
#endif /* IPO_WITH_GMP */
    double real;

    IPO_EXPORT
    Value(const Value& other)
      : real(other.real)
    {
#ifdef IPO_WITH_GMP
      rational = other.rational;
#endif /* IPO_WITH_GMP */
    }
    
    IPO_EXPORT
    Value(double x = 0.0)
      : real(x)
    {

    }

#ifdef IPO_WITH_GMP
    
    IPO_EXPORT
    Value(const mpq_class& x)
      : rational(x), real(x.get_d())
    {

    }

    IPO_EXPORT
    Value(mpq_class& x, bool swap = false)
      : real(x.get_d())
    {
      if (swap)
        rational.swap(x);
      else
        rational = x;
    }

    IPO_EXPORT
    Value(const mpq_t* x)
      : rational(*x), real(rational.get_d())
    {

    }

#endif /* IPO_WITH_GMP */

    IPO_EXPORT
    inline bool isPlusInfinity() const
    {
      return real == std::numeric_limits<double>::infinity();
    }

    IPO_EXPORT
    inline bool isMinusInfinity() const
    {
      return real == -std::numeric_limits<double>::infinity();
    }

    IPO_EXPORT
    inline bool isFinite() const
    {
      return real != std::numeric_limits<double>::infinity()
        && real != -std::numeric_limits<double>::infinity();
    }

    IPO_EXPORT
    inline bool operator<(const Value& other) const
    {
      if (isPlusInfinity() || other.isMinusInfinity())
        return false;
      if (isMinusInfinity() || other.isPlusInfinity())
        return true;
#if defined(IPO_WITH_GMP)
      if (real != other.real)
        return real < other.real;
      else
        return rational < other.rational;
#else
      return real < other.real;
#endif /* IPO_WITH_GMP */
    }

    IPO_EXPORT
    inline bool operator>(const Value& other) const
    {
      return other < *this;
    }

    IPO_EXPORT
    inline bool operator==(const Value& other) const
    {
#if defined(IPO_WITH_GMP)
      if (real == other.real)
        return rational < other.rational;
      else
        return real < other.real;
#else
      return real == other.real;
#endif /* IPO_WITH_GMP */ 
    }

    IPO_EXPORT
    inline bool operator!=(const Value& other) const
    {
      return !(*this == other);
    }

    IPO_EXPORT
    inline bool operator<=(const Value& other) const
    {
      return !(other < *this);
    }

    IPO_EXPORT
    inline bool operator>=(const Value& other) const
    {
      return !(*this < other);
    }
  };

  Value minusInfinity();

  Value plusInfinity();

  typedef unsigned int Index;

  struct Vector
  {
    IPO_EXPORT
    ~Vector();

    IPO_EXPORT
    Vector(const Vector& other);

    IPO_EXPORT
    Vector(Vector&& other)
    {
      _data = other._data;
      _data->usage++;
    }

    IPO_EXPORT
    Vector& operator=(const Vector& other);

    IPO_EXPORT
    Vector& operator=(Vector&& other)
    {
      _data = other._data;
      _data->usage++;
      return *this;
    }

    IPO_EXPORT
    inline bool operator==(const Vector& other) const
    {
      return _data == other._data;
    }

    IPO_EXPORT
    Vector();

    IPO_EXPORT
    Vector(const std::vector<std::size_t>& coordinates, const std::vector<double>& reals);

    IPO_EXPORT
    Vector(const std::vector<std::pair<std::size_t, double>>& nonzeros);

    IPO_EXPORT
    Vector(const std::vector<double>& entries);

    IPO_EXPORT
    std::size_t coordinate(std::size_t index) const;

    IPO_EXPORT
    double real(std::size_t index) const;

#if defined(IPO_WITH_GMP)
    IPO_EXPORT
    Vector(const std::vector<std::size_t>& coordinates, const std::vector<mpq_class>& rationals);

    IPO_EXPORT
    Vector(const std::vector<std::size_t>& coordinates, std::vector<mpq_class>& rationals,
      bool swap = false);

    IPO_EXPORT
    Vector(const std::vector<std::pair<std::size_t, mpq_class>>& nonzeros);

    IPO_EXPORT
    Vector(std::vector<std::pair<std::size_t, mpq_class>>& nonzeros, bool swap = false);

    IPO_EXPORT
    Vector(const std::vector<mpq_class>& entries);

    IPO_EXPORT
    Vector(std::vector<mpq_class>& entries, bool swap = false);

    IPO_EXPORT
    Vector(mpq_t* entries, std::size_t numEntries, bool swap = false);

    IPO_EXPORT
    const mpq_class& rational(std::size_t index) const;
#endif /* IPO_WITH_GMP */

    IPO_EXPORT
    inline Index size() const
    {
      assert(_data != nullptr);
      return _data->size;
    }

    IPO_EXPORT
    inline bool isRational() const
    {
      assert(_data != nullptr);
      return (_data->properties & 1);
    }

    IPO_EXPORT
    inline unsigned short usage() const
    {
      assert(_data != nullptr);
      return _data->usage;
    }

    IPO_EXPORT
    double squaredRealNorm() const;

#if defined(IPO_WITH_GMP)
    IPO_EXPORT
    mpq_class squaredRationalNorm() const;
#endif /* IPO_WITH_GMP */

    struct Header
    {
      unsigned short usage;
      unsigned short properties;
      Index size;
    };

  private:
    friend struct HashVector;

    Header* _data;
  };

  IPO_EXPORT
  Value operator*(const Vector& a, const Vector& b);
  
#if defined(IPO_WITH_GMP)

  IPO_EXPORT
  mpq_class operator*(const Vector& a, const mpq_class* b);
#endif /* IPO_WITH_GMP */  

  IPO_EXPORT
  double operator*(const Vector& a, const double* b);

  struct HashVector
  {
    std::size_t operator() (const Vector& vector) const
    {
      return std::hash<void*>()(vector._data);
    }
  };

  struct Constraint
  {
    Value lhs;
    Vector vector;
    Value rhs;

    IPO_EXPORT
    Constraint(const Value& l, const Vector& vec, const Value& r)
      : lhs(l), vector(vec), rhs(r)
    {
      assert(!lhs.isPlusInfinity());
      assert(!rhs.isMinusInfinity());
      assert(lhs.isFinite() || rhs.isFinite());
    }

    IPO_EXPORT
    inline bool isEquation() const
    {
      return lhs == rhs;
    }

    IPO_EXPORT
    inline bool isLessThanEqual() const
    {
      return lhs.isMinusInfinity();
    }

    IPO_EXPORT
    inline bool isGreaterThanEqual() const
    {
      return rhs.isPlusInfinity();
    }

    IPO_EXPORT
    inline bool isRanged() const
    {
      return lhs < rhs;
    }

    IPO_EXPORT
    inline bool operator==(const Constraint& other) const
    {
      return vector == other.vector && lhs == other.lhs && rhs == other.rhs;
    }

    IPO_EXPORT
    inline bool isAlwaysSatisfied() const
    {
      return vector.size() == 0 && lhs.real <= 0.0 && rhs.real >= 0.0;
    }

    IPO_EXPORT
    inline bool isNeverSatisfied() const
    {
      return (lhs.real > rhs.real) || (vector.size() == 0 && (lhs.real > 0.0 || rhs.real < 0.0));
    }

  };

  IPO_EXPORT
  Constraint alwaysSatisfiedConstraint();

  IPO_EXPORT
  Constraint neverSatisfiedConstraint();

  struct HashConstraint
  {
    std::size_t operator() (const Constraint& constraint) const
    {
      HashVector hash;
      return hash(constraint.vector);
    }
  };  

} /* namespace ipo */
