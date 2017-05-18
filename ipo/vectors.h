#ifndef IPO_VECTORS_IMPL_H_
#define IPO_VECTORS_IMPL_H_

#include <set>

#include "common.h"
#include "rational.h"

#include "vectors-pub.h"

namespace ipo {

  Vector addScaled(const Vector& a, int scaleA, const Vector& b, int scaleB);

  std::size_t differingIndex(const ReferenceCountedVector& a, const ReferenceCountedVector& b);

} /* namespace ipo */

#endif /* IPO_VECTORS_IMPL_H_ */
