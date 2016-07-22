#ifndef IPO_RECONSTRUCT_H_
#define IPO_RECONSTRUCT_H_

#include <vector>

#include "common.h"

namespace ipo {

  /**
   * \brief Reconstructs a rational from a floating-point number.
   */

  void reconstruct(double x, mpq_class& approx, double maxError = 1.0e-6);

  /**
   * \brief Reconstructs a rational from a floating-point number.
   */

  void reconstruct(double x, soplex::Rational& approx, double maxError = 1.0e-6);

  /**
   * \brief Reconstructs a rational vector from from a vector of floating-point numbers.
   */

  void reconstruct(const std::vector<double>& vector, std::vector<mpq_class>& reconstructed, double maxError = 1.0e-6);

  /**
   * \brief Reconstructs a rational vector from from a vector of floating-point numbers.
   */

  void reconstruct(const std::vector<double>& vector, std::vector<soplex::Rational>& reconstructed, double maxError =
      1.0e-6);

} /* namespace ipo */

#endif /* IPO_RECONSTRUCT_H_ */
