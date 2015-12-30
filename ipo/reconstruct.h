#ifndef IPO_RECONSTRUCT_H_
#define IPO_RECONSTRUCT_H_

#include <vector>

#include <gmpxx.h>

#define SOPLEX_WITH_GMP
#include <rational.h>

namespace ipo {

  void reconstruct(double x, mpq_class& approx, double maxError = 1.0e-6);
  void reconstruct(double x, soplex::Rational& approx, double maxError = 1.0e-6);
  void reconstruct(const std::vector<double>& vector,
      std::vector<mpq_class>& reconstructed, double maxError = 1.0e-6);
  void reconstruct(const std::vector<double>& vector,
      std::vector<soplex::Rational>& reconstructed, double maxError = 1.0e-6);
  void scalePrimitiveIntegral(const std::vector<mpq_class>& input,
      std::vector<mpz_class>& output);
  void reconstructIntegralScaled(const std::vector<double>& vector,
      std::vector<mpz_class>& reconstructed);

}

#endif /* IPO_RECONSTRUCT_H_ */
