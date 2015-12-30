#include "vector_2d.h"

#include <iostream>

namespace ipo {

  Integer2d::Integer2d() :
      x(0), y(0)
  {

  }

  Integer2d::Integer2d(const mpz_class& ix, const mpz_class& iy) :
      x(ix), y(iy)
  {

  }

  Integer2d::Integer2d(const Integer2d& other) :
      x(other.x), y(other.y)
  {

  }

  Integer2d::~Integer2d()
  {

  }

  mpz_class Integer2d::normalize()
  {
    if (x != 0 || y != 0)
    {
      mpz_class result = 1;
      mpz_gcd(result.get_mpz_t(), x.get_mpz_t(), y.get_mpz_t());
      x /= result;
      y /= result;
      return result;
    }
    else
      return 1;
  }

  mpz_class Integer2d::maximumNorm() const
  {
    mpz_class xAbs, yAbs;
    mpz_abs(xAbs.get_mpz_t(), x.get_mpz_t());
    mpz_abs(yAbs.get_mpz_t(), y.get_mpz_t());
    return xAbs >= yAbs ? xAbs : yAbs;
  }

  bool Integer2d::operator<(const Integer2d& other) const
  {
    if (x > 0)
    {
      if (other.x > 0)
      {
        mpz_class det = determinant(*this, other);
        if (det < 0)
          return true;
        else if (det > 0)
          return false;
        return x < other.x;
      }
      else if (other.x < 0)
        return true;
      else
        return other.y < 0;
    }
    else if (x < 0)
    {
      if (other.x >= 0)
        return false;
      else
      {
        mpz_class det = determinant(*this, other);
        if (det < 0)
          return false;
        else if (det > 0)
          return true;
        return x > other.x;
      }
    }
    else
    {
      if (y >= 0)
        return !(other.x == 0 && other.y >= 0 && other.y <= y);
      return other.x < 0 || (other.x == 0 && other.y <= y);
    }
  }

  Integer2d Integer2d::operator+(const Integer2d& rhs) const
  {
    return Integer2d(x + rhs.x, y + rhs.y);
  }

  Integer2d& Integer2d::operator+=(const Integer2d& rhs)
  {
    x += rhs.x;
    y += rhs.y;
    return *this;
  }

  Integer2d Integer2d::operator-(const Integer2d& rhs) const
  {
    return Integer2d(x - rhs.x, y - rhs.y);
  }

  Integer2d& Integer2d::operator-=(const Integer2d& rhs)
  {
    x -= rhs.x;
    y -= rhs.y;
    return *this;
  }

  mpz_class Integer2d::operator*(const Integer2d& rhs) const
  {
    return x * rhs.x + y * rhs.y;
  }

  Integer2d Integer2d::operator*(const mpz_class& rhs) const
  {
    return Integer2d(x * rhs, y * rhs);
  }

  Integer2d& Integer2d::operator*=(const mpz_class& rhs)
  {
    x *= rhs;
    y *= rhs;
    return *this;
  }

  std::ostream& operator<<(std::ostream& stream, const Integer2d& vector)
  {
    return stream << "(" << vector.x << "," << vector.y << ")";
  }

  mpz_class determinant(const Integer2d& left, const Integer2d& right)
  {
    return left.x * right.y - left.y * right.x;
  }

  Rational2d::Rational2d() :
      x(0), y(0)
  {
  }

  Rational2d::Rational2d(const mpq_class& ix, const mpq_class& iy) :
      x(ix), y(iy)
  {
  }
  Rational2d::Rational2d(const Rational2d& other) :
      x(other.x), y(other.y)
  {
  }

  Rational2d::~Rational2d()
  {
  }

  Integer2d Rational2d::normalized() const
  {
    Integer2d result(x.get_num() * y.get_den(), y.get_num() * x.get_den());
    result.normalize();
    return result;
  }

  mpq_class Rational2d::maximumNorm() const
  {
    mpq_class xAbs, yAbs;
    mpq_abs(xAbs.get_mpq_t(), x.get_mpq_t());
    mpq_abs(yAbs.get_mpq_t(), y.get_mpq_t());
    return xAbs >= yAbs ? xAbs : yAbs;
  }

  bool Rational2d::operator<(const Rational2d& other) const
  {
    if (x > 0)
    {
      if (other.x > 0)
      {
        mpq_class det = determinant(*this, other);
        if (det < 0)
          return true;
        else if (det > 0)
          return false;
        return x < other.x;
      }
      else if (other.x < 0)
        return true;
      else
        return other.y < 0;
    }
    else if (x < 0)
    {
      if (other.x >= 0)
        return false;
      else
      {
        mpq_class det = determinant(*this, other);
        if (det < 0)
          return false;
        else if (det > 0)
          return true;
        return x > other.x;
      }
    }
    else
    {
      if (y >= 0)
        return !(other.x == 0 && other.y >= 0 && other.y <= y);
      return other.x < 0 || (other.x == 0 && other.y <= y);
    }
  }

  Rational2d Rational2d::operator+(const Rational2d& rhs) const
  {
    return Rational2d(x + rhs.x, y + rhs.y);
  }

  Rational2d& Rational2d::operator+=(const Rational2d& rhs)
  {
    x += rhs.x;
    y += rhs.y;
    return *this;
  }

  Rational2d Rational2d::operator-(const Rational2d& rhs) const
  {
    return Rational2d(x - rhs.x, y - rhs.y);
  }

  Rational2d& Rational2d::operator-=(const Rational2d& rhs)
  {
    x -= rhs.x;
    y -= rhs.y;
    return *this;
  }

  mpq_class Rational2d::operator*(const Rational2d& rhs) const
  {
    return x * rhs.x + y * rhs.y;
  }

  Rational2d Rational2d::operator*(const mpq_class& rhs) const
  {
    return Rational2d(x * rhs, y * rhs);
  }

  Rational2d& Rational2d::operator*=(const mpq_class& rhs)
  {
    x *= rhs;
    y *= rhs;
    return *this;
  }

  std::ostream& operator<<(std::ostream& stream, const Rational2d& vector)
  {
    return stream << "(" << vector.x << "," << vector.y << ")";
  }

  mpq_class determinant(const Rational2d& left, const Rational2d& right)
  {
    return left.x * right.y - left.y * right.x;
  }

} /* namespace polycomb */
