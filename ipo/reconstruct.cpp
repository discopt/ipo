#include "reconstruct.h"

#include <cmath>
#include <cassert>

namespace ipo {

  void reconstruct(double x, mpq_class& approx, double maxError)
  {
    const double input = x;

    /**
     * Reconstruction using continued fractions:
     *
     * x = floor(x) + x_1
     * (a_1,b_1,c_1,d_1) = (1,0,floor(x),1)
     *
     * Invariant: x = (d_i*x_i + c_i) / (b_i*x_i + a_i)
     *
     * Satisfied for i=1 because (1*x_1 + floor(x)) / (0*x_1 + 1) = x_1 + floor(x) = x
     *
     * Let z := floor(1/x_i) and x_{i+1} := 1/x_i - z. Then x_i = 1 / (z + x_{i+1}) holds.
     * (a,b,c,d)_{i+1} := ( (b+a*z), (a), (d+c*z), (c))
     *
     * Induction:
     *
     * x
     * = (d_i*x_i + c_i) / (b_i*x_i + a_i)
     * = (d_i + c_i*(z+x_{i+1}) / (b_i + a_i*(z+x_{i+1}))
     * = (c_i*x_{i+1} + (d_i + c_i*z)) / (a_i*x_{i+1} + (b_i + a_i*z))
     *     ^                 ^             ^                 ^
     *   d_{i+1}           c_{i+1}       b_{i+1}           a_{i+1}
     *
     * The final approximation is c_i / a_i via setting x_i = 0.
     **/

    double fx = floor(x);
    mpz_class a(1);
    mpz_class b(0);
    mpz_class c(mpq_class(round(fx)));
    mpz_class d(1);
    x -= fx;

    mpz_class z, bOld, dOld;
    while (true)
    {
      approx.get_num() = c;
      approx.get_den() = a;
      approx.canonicalize();
      double error = fabs(approx.get_d() - input);

      if (error < maxError)
        return;

      x = 1.0 / x;
      fx = floor(x);
      x -= fx;
      z = static_cast<long int>(round(fx));
      bOld = b;
      dOld = d;
      b = a;
      d = c;
      a = z * a + bOld;
      c = z * c + dOld;
    }
  }

  void reconstruct(double x, soplex::Rational& approx, double maxError)
  {
    mpq_class r;
    reconstruct(x, r, maxError);
    mpq_t q = { *r.get_mpq_t() };
    approx = soplex::Rational(q);
  }

  void reconstruct(const std::vector<double>& vector,
      std::vector<mpq_class>& reconstructed, double maxError)
  {
    reconstructed.resize(vector.size(), 0);
    for (std::size_t i = 0; i < vector.size(); ++i)
    {
      reconstruct(vector[i], reconstructed[i], maxError);
    }
  }

  void reconstruct(const std::vector<double>& vector,
      std::vector<soplex::Rational>& reconstructed, double maxError)
  {
    reconstructed.resize(vector.size(), 0);
    for (std::size_t i = 0; i < vector.size(); ++i)
    {
      reconstruct(vector[i], reconstructed[i], maxError);
    }
  }

  void scalePrimitiveIntegral(const std::vector<mpq_class>& input,
      std::vector<mpz_class>& output)
  {
    mpz_class gcd = 0;
    mpz_class lcm = 1;
    for (std::size_t i = 0; i < input.size(); ++i)
    {
      mpz_lcm(lcm.get_mpz_t(), lcm.get_mpz_t(), input[i].get_den_mpz_t());
      mpz_gcd(gcd.get_mpz_t(), gcd.get_mpz_t(), input[i].get_num_mpz_t());
    }
    output.resize(input.size());
    if (gcd > 0)
    {
      mpq_class scalar(lcm, gcd);
      for (std::size_t i = 0; i < input.size(); ++i)
      {
        mpq_class scaled = input[i] * scalar;
        scaled.canonicalize();
        assert(scaled.get_den() == 1);
        output[i] = scaled.get_num();
      }
    }
    else
    {
      for (std::size_t i = 0; i < input.size(); ++i)
        output[i] = 0;
    }

#ifdef IPO_DEBUG
    gcd = 0;
    for (std::size_t i = 0; i < output.size(); ++i)
      mpz_gcd(gcd.get_mpz_t(), gcd.get_mpz_t(), output[i].get_mpz_t());
    assert(gcd == 0 || gcd == 1);
#endif
  }

  void reconstructIntegralScaled(const std::vector<double>& vector,
      std::vector<mpz_class>& reconstructed)
  {
    std::vector<mpq_class> unscaled;
    reconstruct(vector, unscaled);
    scalePrimitiveIntegral(unscaled, reconstructed);
  }
}

