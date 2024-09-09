// #define IPO_DEBUG /* Uncomment to debug this file. */

#include <ipo/arithmetic.hpp>

#include <cmath>
#include <random>

#if defined(IPO_DEBUG)
#include <iostream>
#endif /* IPO_DEBUG */

namespace ipo
{

  double squaredEuclideanNorm(double* vector, std::size_t size)
  {
    double result = 0.0;
    for (std::size_t i = 0; i < size; ++i)
    {
      double x = vector[i];
      result += x * x;
    }
    return result;
  }

  double euclideanNorm(double* vector, std::size_t size)
  {
    return sqrt(squaredEuclideanNorm(vector, size));
  }

  double* generateRandomVectorSphere(std::size_t size)
  {
    double* result = new double[size];
    std::default_random_engine generator;
    std::normal_distribution<double> distribution;
    double squaredNorm = 0.0;
    while (squaredNorm < 1.0e-3)
    {
      for (std::size_t v = 0; v < size; ++v)
      {
        double x = distribution(generator);
        result[v] = x;
        squaredNorm += x*x;
      }
    }
    double norm = sqrt(squaredNorm);
    for (std::size_t v = 0; v < size; ++v)
      result[v] /= norm;

    return result;
  }

  double maxAbsoluteValue(const double* vector, std::size_t length)
  {
    double result = 0.0;
    const double* px = vector;
    while (length)
    {
      double x = fabs(*px);
      if (x > result)
        result = x;
      ++px;
      --length;
    }
    return result;
  }

#if defined(IPO_RATIONAL)

  void reconstructRational(mpq_ptr result, double x, double maxError)
  {
    const double input = x;
#if defined(IPO_DEBUG)
    std::cout << "reconstructRational(" << input << ")" << std::endl;
#endif /* IPO_DEBUG */

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
    mpz_t a,b,c,d;
    mpz_init_set_si(a, 1);
    mpz_init_set_si(b, 0);
    mpz_init_set_d(c, round(fx));
    mpz_init_set_si(d, 1);
    x -= fx;

    mpz_t z, bOld, dOld;
    mpz_init(z);
    mpz_init(bOld);
    mpz_init(dOld);
    while (true)
    {
      mpq_set_num(result, c);
      mpq_set_den(result, a);
      mpq_canonicalize(result);
      double error = fabs(mpq_get_d(result) - input);

#if defined(IPO_DEBUG)
      std::cout << "x = " << x << ", floorx = " << fx << ", round(floorx) = " << round(fx) << ", a = " << mpz_get_str(0, 10, a)
        << ", b = " << mpz_get_str(0, 10, b) << ", c = " << mpz_get_str(0, 10, c) << ", d = " << mpz_get_str(0, 10, d)
        << ", c/a = " << mpq_get_str(0, 10, result) << ", error = " << error << std::endl;
#endif /* IPO_DEBUG */

      if (error < maxError)
        return;

      x = 1.0 / x;
      fx = floor(x);
      x -= fx;
      mpz_set_si(z, static_cast<long int>(round(fx)));

#if defined(IPO_DEBUG)
      std::cout << "x = " << x << ", floorx = " << fx << ", round(floorx) = " << round(fx) << ", z = " << mpz_get_str(0, 10, z) << std::endl;
#endif /* IPO_DEBUG */

      mpz_set(bOld, b);
      mpz_set(dOld, d);
      mpz_set(b, a);
      mpz_set(d, c);
      mpz_mul(a, a, z);
      mpz_add(a, a, bOld);
      mpz_mul(c, c, z);
      mpz_add(c, c, dOld);
    }
  }

  rational reconstructRational(double x, double maxError)
  {
    rational result;
    reconstructRational(result.backend().data(), x, maxError);
    return result;
  }

  IntegralScaler::IntegralScaler()
    : _factor(1)
  {

  }

  void IntegralScaler::operator()(const rational& x)
  {
    if (_factor == 1)
      _factor = 1 / x;
    else
    {
      mpz_gcd(mpq_denref(_factor.backend().data()), mpq_denref(_factor.backend().data()), mpq_numref(x.backend().data()));
      mpz_lcm(mpq_numref(_factor.backend().data()), mpq_numref(_factor.backend().data()), mpq_denref(x.backend().data()));
    }
  }

  const rational& IntegralScaler::factor() const
  {
    return _factor;
  }

#endif /* IPO_RATIONAL */

} /* namespace ipo */

