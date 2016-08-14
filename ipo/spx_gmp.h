#ifndef IPO_SPX_GMP_H_
#define IPO_SPX_GMP_H_

#include "common.h"
#include "rational.h"

namespace ipo {

  soplex::Rational mpq2rational(const mpq_class& x);
  soplex::Rational mpz2rational(const mpz_class& x);
  mpq_class rational2mpq(const soplex::Rational& x);
  mpz_class rational2mpzNum(const soplex::Rational& x);
  mpz_class rational2mpzDen(const soplex::Rational& x);

} /* namespace ipo */

#endif /* IPO_SPX_GMP_H_ */
