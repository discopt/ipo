#ifndef IPO_RATIONAL_H_
#define IPO_RATIONAL_H_

#include "common.h"

#include <gmpxx.h>
#define SOPLEX_WITH_GMP
#define SOPLEX_WITH_RATIONALPARAMS
#include <soplex.h>

namespace ipo {

  const double plusInfinity = soplex::infinity;
  const double minusInfinity = -soplex::infinity;

  typedef soplex::Rational Rational;

} /* namespace ipo */

#endif /* IPO_RATIONAL_H_ */
