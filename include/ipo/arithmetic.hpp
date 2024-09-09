#pragma once

#include <cstddef>

#include <ipo/config.hpp>
#include <ipo/export.hpp>

namespace ipo
{

  double squaredEuclideanNorm(double* vector, std::size_t size);

  double euclideanNorm(double* vector, std::size_t size);

  template <typename To, typename From>
  To convertNumber(const From& from)
  {
    return To::unimplemented;
  }

  template<>
  inline double convertNumber<double>(const double& from)
  {
    return from;
  }

  double* generateRandomVectorSphere(std::size_t size);

  double maxAbsoluteValue(const double* vector, std::size_t length);

} /* namespace ipo */

#if defined(IPO_RATIONAL)

#include <gmp.h>
#include <boost/multiprecision/number.hpp>
#include <boost/multiprecision/gmp.hpp>

namespace ipo
{
  using rational = boost::multiprecision::number<boost::multiprecision::gmp_rational, boost::multiprecision::et_off>;

  template<>
  inline rational convertNumber<rational>(const rational& x)
  {
    return x;
  }

  template<>
  inline double convertNumber<double>(const rational& x)
  {
    return x.convert_to<double>();
  }

  rational reconstructRational(double x, double maxError = 1.0e-9);

  void reconstructRational(mpq_ptr result, double x, double maxError = 1.0e-9 );

  template<>
  inline rational convertNumber<rational>(const double& x)
  {
    return reconstructRational(x);
  }

  class IntegralScaler
  {
  public:
    IntegralScaler();

    void operator()(const rational& x);

    const rational& factor() const;

  private:
    rational _factor;
  };

}

#endif /* IPO_RATIONAL */
