#include <ipo/arithmetic.hpp>

#include <cmath>

#if defined(IPO_WITH_GMP)

void reconstructRational(mpq_ptr result, double x, double maxError)
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
    mpq_set_num(result, c.get_mpz_t());
    mpq_set_den(result, a.get_mpz_t());
    mpq_canonicalize(result);
    double error = fabs(mpq_get_d(result) - input);

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

mpq_class reconstructRational(double x, double maxError)
{
  mpq_class result;
  reconstructRational(result.get_mpq_t(), x, maxError);
  return result;
}

namespace ipo
{

  IntegralScaler::IntegralScaler()
    : _factor(1, 1)
  {

  }

  void IntegralScaler::operator()(const mpq_class& x)
  {
    mpz_gcd(_factor.get_den_mpz_t(), _factor.get_den_mpz_t(), x.get_num_mpz_t());
    mpz_lcm(_factor.get_num_mpz_t(), _factor.get_num_mpz_t(), x.get_den_mpz_t());
  }

  const mpq_class& IntegralScaler::factor() const
  {
    return _factor;
  }

}

#endif /* IPO_WITH_GMP */

