#ifndef IPO_LINEAR_CONSTRAINT_IMPL_H_
#define IPO_LINEAR_CONSTRAINT_IMPL_H_

#include "linear_constraint-pub.h"

namespace ipo {

  LinearConstraint addScaled(char type, const LinearConstraint& a, int scaleA, const LinearConstraint& b, int scaleB);

} /* namespace ipo */

#endif /* IPO_LINEAR_CONSTRAINT_IMPL_H_ */