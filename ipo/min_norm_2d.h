#ifndef IPO_MIN_NORM_2D_H_
#define IPO_MIN_NORM_2D_H_

#include "ipo.h"

namespace ipo {

  bool manhattanNormShortestCombination(std::size_t n, soplex::DSVectorRational& newTarget,
      const soplex::SVectorRational& target, const soplex::SVectorRational& source, soplex::Rational& targetMultiplier,
      soplex::Rational& sourceMultiplier, soplex::Rational& norm);

  bool manhattanNormGreedyCombination(std::size_t n, soplex::DSVectorRational& newTarget,
      const soplex::SVectorRational& target, const soplex::SVectorRational& source, soplex::Rational& targetMultiplier,
      soplex::Rational& sourceMultiplier, soplex::Rational& norm);

  void manhattanNormImproveEquations(std::size_t n, soplex::LPRowSetRational& equations);

  void manhattanNormImproveInequality(std::size_t n, soplex::LPRowRational& inequality,
      const soplex::LPRowSetRational& equations);

} /* namespace ipo */

#endif /* IPO_MIN_NORM_2D_H_ */
