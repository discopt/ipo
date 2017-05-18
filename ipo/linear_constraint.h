#ifndef IPO_LINEAR_CONSTRAINT_IMPL_H_
#define IPO_LINEAR_CONSTRAINT_IMPL_H_

#include "linear_constraint-pub.h"

#include "soplex_reproduce.h"

namespace ipo {

  LinearConstraint addScaled(char type, const LinearConstraint& a, int scaleA, const LinearConstraint& b, int scaleB);
  
  // DEBUG:
  
  void addToLP(soplex::ReproSoPlex& spx, const LinearConstraint& constraint);

  void addToLP(soplex::ReproSoPlex& spx, const std::vector<LinearConstraint>& constraints);


} /* namespace ipo */

#endif /* IPO_LINEAR_CONSTRAINT_IMPL_H_ */
