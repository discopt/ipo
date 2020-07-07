#pragma once

#include <ipo/config.hpp>
#include <ipo/export.hpp>

#if defined(IPO_WITH_GMP)

#include <gmpxx.h>
#include <cmath>
#include <cassert>
#include <ostream>

#endif /* IPO_WITH_GMP */

namespace ipo
{
  
#if defined(IPO_WITH_GMP)

  struct rational
  {
  public:
    IPO_EXPORT
    rational()
      : _exact(0), _approx(0)
    {

    }

    IPO_EXPORT
    rational(const rational& other)
      : _exact(other._exact), _approx(other._approx)
    {

    }

    IPO_EXPORT
    rational(double value)
      : _approx(value)
    {
      if (std::isfinite(value))
        _exact = value;
    }

    IPO_EXPORT
    rational(long numerator, long denominator)
      : _exact(numerator, denominator), _approx(_exact.get_d())
    {

    }

    IPO_EXPORT
    rational(rational&& other)
      : _approx(other._approx)
    {
      _exact.swap(other._exact);
    }

    IPO_EXPORT
    explicit operator double() const
    {
      return _approx;
    }

    IPO_EXPORT
    rational(const mpq_class& x)
      : _exact(x), _approx(x.get_d())
    {

    }

    IPO_EXPORT
    rational(mpq_class&& x)
      : _approx(x.get_d())
    {
      _exact.swap(x);
    }

    IPO_EXPORT
    rational(const mpq_t& x)
      : _exact(x), _approx(mpq_get_d(x))
    {

    }

    IPO_EXPORT
    ~rational()
    {

    }

    IPO_EXPORT
    inline rational& operator=(const rational& other)
    {
      _exact = other._exact;
      _approx = other._approx;
      return *this;
    }

    IPO_EXPORT
    inline rational& operator=(rational&& other)
    {
      _exact.swap(other._exact);
      _approx = other._approx;
      return *this;
    }

    IPO_EXPORT
    inline bool isFinite() const
    {
      return std::isfinite(_approx);
    }

    IPO_EXPORT
    inline bool isNAN() const
    {
      return std::isnan(_approx);
    }

    IPO_EXPORT
    inline bool isPlusInfinity() const
    {
      return _approx == std::numeric_limits<double>::infinity();
    }

    IPO_EXPORT inline bool isMinusInfinity() const
    {
      return _approx == -std::numeric_limits<double>::infinity();
    }

    IPO_EXPORT
    inline bool operator<(const rational& other) const
    {
      if (isPlusInfinity() || other.isMinusInfinity())
        return false;
      if (isMinusInfinity() || other.isPlusInfinity())
        return true;
      if (_approx != other._approx)
        return _approx < other._approx;
      else
        return _exact < other._exact;
    }

    IPO_EXPORT
    inline bool operator>(const rational& other) const
    {
      return other < *this;
    }

    IPO_EXPORT
    inline bool operator==(const rational& other) const
    {
      if (_approx == other._approx)
        return _exact == other._exact;
      else
        return false;
    }

    IPO_EXPORT
    inline bool operator!=(const rational& other) const
    {
      return !(*this == other);
    }

    IPO_EXPORT
    inline bool operator<=(const rational& other) const
    {
      return !(other < *this);
    }

    IPO_EXPORT
    inline bool operator>=(const rational& other) const
    {
      return !(*this < other);
    }

    IPO_EXPORT
    inline rational operator+(const rational& other) const
    {
      if (isFinite() && other.isFinite())
        return rational(this->_exact + other._exact);
      return rational(this->_approx + other._approx);
    }

    IPO_EXPORT
    inline rational operator-(const rational& other) const
    {
      if (isFinite() && other.isFinite())
        return rational(this->_exact - other._exact);
      return rational(this->_approx - other._approx);
    }

    IPO_EXPORT
    inline rational operator*(const rational& other) const
    {
      if (isFinite() && other.isFinite())
        return rational(this->_exact * other._exact);
      return rational(this->_approx * other._approx);
    }

    IPO_EXPORT
    inline rational operator/(const rational& other) const
    {
      if (isFinite() && other.isFinite())
        return rational(this->_exact / other._exact);
      return rational(this->_approx / other._approx);
    }

    IPO_EXPORT
    inline rational& operator+=(const rational& other)
    {
      this->_exact += other._exact;
      this->_approx = this->_exact.get_d();
      return *this;
    }

    IPO_EXPORT
    inline rational& operator-=(const rational& other)
    {
      this->_exact -= other._exact;
      this->_approx = this->_exact.get_d();
      return *this;
    }

    IPO_EXPORT
    inline rational& operator*=(const rational& other)
    {
      this->_exact *= other._exact;
      this->_approx = this->_exact.get_d();
      return *this;
    }

    IPO_EXPORT
    inline rational& operator*=(double other)
    {
      this->_exact *= other;
      this->_approx = this->_exact.get_d();
      return *this;
    }

    IPO_EXPORT
    inline rational& operator/=(const rational& other)
    {
      this->_exact /= other._exact;
      this->_approx = this->_exact.get_d();
      return *this;
    }

    IPO_EXPORT
    inline const mpq_class& get_mpq_class() const
    {
      return _exact;
    }

    IPO_EXPORT
    inline mpq_srcptr get_mpq_t() const
    {
      return _exact.get_mpq_t();
    }

    IPO_EXPORT
    inline double approximation() const
    {
      return _approx;
    }

  private:
    mpq_class _exact;
    double _approx;
  };

  rational minusInfinity();

  rational plusInfinity();

  std::ostream& operator<<(std::ostream& stream, const rational& x);

  IPO_EXPORT
  inline bool isPlusInfinity(const rational& x)
  {
    return x.isPlusInfinity();
  }

  IPO_EXPORT
  inline bool isMinusInfinity(const rational& x)
  {
    return x.isMinusInfinity();
  }

  IPO_EXPORT
  inline bool isNAN(const rational& x)
  {
    return x.isNAN();
  }

  IPO_EXPORT
  inline bool isFinite(const rational& x)
  {
    return x.isFinite();
  }

  struct RationalIsZero
  {
    bool operator()(const rational& x) const
    {
      return x == 0;
    }
  };
  
#endif /* IPO_WITH_GMP */  
  
  IPO_EXPORT
  inline bool isPlusInfinity(double x)
  {
    return x == std::numeric_limits<double>::infinity();
  }

  IPO_EXPORT
  inline bool isMinusInfinity(double x)
  {
    return x == -std::numeric_limits<double>::infinity();
  }

  IPO_EXPORT
  inline bool isNAN(double x)
  {
    return std::isnan(x);
  }

  IPO_EXPORT
  inline bool isFinite(double x)
  {
    return std::isfinite(x);
  }

  struct DoubleIsZero
  {
    double epsilon;

    DoubleIsZero(double eps)
      : epsilon(eps)
    {

    }

    bool operator()(double x) const
    {
      return fabs(x) < epsilon;
    }
  };
}
