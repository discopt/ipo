#ifndef IPO_MIN_NORM_2D_H_
#define IPO_MIN_NORM_2D_H_

#include "common.h"
#include "rational.h"
#include "linear_constraint.h"

namespace ipo {

  bool manhattanNormShortestCombination(std::size_t n, Vector& newTarget, const Vector& target, const Vector& source, 
    Rational& targetMultiplier, Rational& sourceMultiplier, Rational& norm);

  bool manhattanNormGreedyCombination(std::size_t n, Vector& newTarget, const Vector& target, const Vector& source, 
    Rational& targetMultiplier, Rational& sourceMultiplier, Rational& norm);

  void manhattanNormImproveEquations(std::size_t n, std::vector<LinearConstraint>& equations);

  void manhattanNormImproveInequality(std::size_t n, LinearConstraint& inequality,
      const std::vector<LinearConstraint>& equations);

} /* namespace ipo */

#endif /* IPO_MIN_NORM_2D_H_ */
