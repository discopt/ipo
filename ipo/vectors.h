#ifndef IPO_VECTORS_IMPL_H_
#define IPO_VECTORS_IMPL_H_

#include <set>

#include "common.h"
#include "rational.h"

#include "vectors-pub.h"

namespace ipo {

  Vector addScaled(const Vector& a, int scaleA, const Vector& b, int scaleB);

  std::size_t differingIndex(const ReferenceCountedVector& a, const ReferenceCountedVector& b);

  void vectorToDense(const ReferenceCountedVector& source, soplex::VectorRational& target);
  void vectorToSparse(const ReferenceCountedVector& source, soplex::SVectorRational& target);
  MutableVector denseToVector(const soplex::VectorRational& source, bool saveMemory = true);
  MutableVector sparseToVector(const soplex::SVectorRational& source);

  Rational operator*(const soplex::VectorRational& a, const ReferenceCountedVector& b);
  Rational operator*(const ReferenceCountedVector& a, const soplex::VectorRational& b);
  soplex::VectorRational& operator+=(soplex::VectorRational& a, const ReferenceCountedVector& b);
  soplex::VectorRational& operator-=(soplex::VectorRational& a, const ReferenceCountedVector& b);

  MutableVector integralScaled(const Vector& vector, soplex::Rational* factor = NULL);
  void scaleIntegral(Vector& vector, soplex::Rational* factor = NULL);
  void scaleIntegral(std::vector<Vector>& vectors);

  void scaleIntegral(const soplex::VectorRational& vector, soplex::DVectorRational& scaled);

} /* namespace ipo */

#endif /* IPO_VECTORS_IMPL_H_ */
