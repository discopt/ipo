#ifndef IPO_RATIONAL_H_
#define IPO_RATIONAL_H_

#include "common.h"

#include <string>

#include <gmpxx.h>
#define SOPLEX_WITH_GMP
#define SOPLEX_WITH_RATIONALPARAMS
#include <soplex.h>

namespace ipo {

  typedef soplex::Rational Rational;

    std::string rationalToString(const Rational& number);

  bool isIntegral(const Rational& number);

  class IntegralScaler
  {
  public:
    IntegralScaler();
    ~IntegralScaler();

    void operator()(const Rational& number);

    const Rational factor() const;

  private:
    mpz_class _numScaler;
    mpz_class _denScaler;
  };


} /* namespace ipo */

#endif /* IPO_RATIONAL_H_ */
