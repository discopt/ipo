#ifndef IPO_VECTOR_2D_H_
#define IPO_VECTOR_2D_H_

#include "common.h"
#include "rational.h"

namespace ipo {

  class Integer2d
  {
  public:
    mpz_class x;
    mpz_class y;

  public:
    Integer2d();
    Integer2d(const mpz_class& x, const mpz_class& y);
    Integer2d(const Integer2d& other);
    virtual ~Integer2d();

    mpz_class normalize();
    mpz_class maximumNorm() const;
    bool operator<(const Integer2d& other) const;
    Integer2d operator+(const Integer2d& rhs) const;
    Integer2d& operator+=(const Integer2d& rhs);
    Integer2d operator-(const Integer2d& rhs) const;
    Integer2d& operator-=(const Integer2d& rhs);
    mpz_class operator*(const Integer2d& rhs) const;
    Integer2d operator*(const mpz_class& rhs) const;
    Integer2d& operator*=(const mpz_class& rhs);
  };

  std::ostream& operator<<(std::ostream& stream, const Integer2d& vector);
  mpz_class determinant(const Integer2d& left, const Integer2d& right);

  class Rational2d
  {
  public:
    mpq_class x;
    mpq_class y;

  public:
    Rational2d();
    Rational2d(const mpq_class& x, const mpq_class& y);
    Rational2d(const Rational2d& other);
    virtual ~Rational2d();

    Integer2d normalized() const;
    mpq_class maximumNorm() const;
    bool operator<(const Rational2d& other) const;
    Rational2d operator+(const Rational2d& rhs) const;
    Rational2d& operator+=(const Rational2d& rhs);
    Rational2d operator-(const Rational2d& rhs) const;
    Rational2d& operator-=(const Rational2d& rhs);
    mpq_class operator*(const Rational2d& rhs) const;
    Rational2d operator*(const mpq_class& rhs) const;
    Rational2d& operator*=(const mpq_class& rhs);
  };

  std::ostream& operator<<(std::ostream& stream, const Rational2d& vector);
  mpq_class determinant(const Rational2d& left, const Rational2d& right);

} /* namespace ipo */

#endif /* IPO_VECTOR_2D_H_ */
