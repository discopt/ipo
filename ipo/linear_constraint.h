#ifndef IPO_LINEAR_CONSTRAINT_IMPL_H_
#define IPO_LINEAR_CONSTRAINT_IMPL_H_

#include "linear_constraint-pub.h"

namespace ipo {

  LinearConstraint addScaled(char type, const LinearConstraint& a, int scaleA, const LinearConstraint& b, int scaleB);

  void addToLP(soplex::SoPlex& spx, const LinearConstraint& constraint);

  void addToLP(soplex::SoPlex& spx, const std::vector<LinearConstraint>& constraints);
  
  LinearConstraint integralScaled(const LinearConstraint& constraint);
  void scaleIntegral(LinearConstraint& constraint);
  void scaleIntegral(std::vector<LinearConstraint>& constraints);


} /* namespace ipo */

#endif /* IPO_LINEAR_CONSTRAINT_IMPL_H_ */