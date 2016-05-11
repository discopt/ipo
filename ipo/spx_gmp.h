#ifndef IPO_SPX_GMP_H_
#define IPO_SPX_GMP_H_

#include "ipo.h"

soplex::Rational mpq2rational(const mpq_class& x);
soplex::Rational mpz2rational(const mpz_class& x);
mpq_class rational2mpq(const soplex::Rational& x);
mpz_class rational2mpzNum(const soplex::Rational& x);
mpz_class rational2mpzDen(const soplex::Rational& x);

void integralScaleScan(mpz_class& numScaler, mpz_class& denScaler, const soplex::Rational& x);
soplex::Rational integralScaleUpdate(const mpz_class& numScaler, const mpz_class& denScaler, const soplex::Rational& x);
void integralScaleUpdate(const mpz_class& numScaler, const mpz_class& denScaler, soplex::Rational& x);

#endif /* IPO_SPX_GMP_H_ */
