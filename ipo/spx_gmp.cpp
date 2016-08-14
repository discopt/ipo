#include "spx_gmp.h"

namespace ipo {

  soplex::Rational mpq2rational(const mpq_class& x)
  {
    mpq_t y;
    mpq_init(y);
    mpq_set(y, x.get_mpq_t());
    soplex::Rational r(y);
    mpq_clear(y);
    return r;
  }

  soplex::Rational mpz2rational(const mpz_class& x)
  {
    mpq_class y = x;
    mpq_t z;
    mpq_init(z);
    mpq_set(z, y.get_mpq_t());
    soplex::Rational r(z);
    mpq_clear(z);
    return r;
  }

  mpq_class rational2mpq(const soplex::Rational& x)
  {
    mpq_class y(x.getMpqRef());
    return y;
  }

  mpz_class rational2mpzNum(const soplex::Rational& x)
  {
    mpq_class y = rational2mpq(x);
    return mpz_class(y.get_num());
  }

  mpz_class rational2mpzDen(const soplex::Rational& x)
  {
    mpq_class y = rational2mpq(x);
    return mpz_class(y.get_den());
  }

} /* namespace ipo */