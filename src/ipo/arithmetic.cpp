#include <ipo/arithmetic.hpp>

#include <cmath>

namespace ipo
{
  DoubleIsZero::DoubleIsZero(double eps)
    : epsilon(eps)
  {

  }

  bool DoubleIsZero::operator()(double value) const
  {
    return fabs(value) < epsilon;
  }


  void mpq_reconstruct(mpq_t& result, double x, double maxError)
  {
    mpq_class approx = reconstruct(x, maxError);
    mpq_set(result, approx.get_mpq_t());
  }

  mpq_class reconstruct(double x, double maxError)
  {
    mpq_class approx;
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
        return approx;

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

    return approx;
  }
}
