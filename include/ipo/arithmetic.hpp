#pragma once

#include <ipo/config.hpp>
#include <ipo/export.hpp>

#if defined(IPO_WITH_GMP)
#include <gmpxx.h>
#endif /* IPO_WITH_GMP */

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

#if defined(IPO_WITH_GMP)

template<>
inline mpq_class convertNumber<mpq_class>(const mpq_class& x)
{
  return x;
}

template<>
inline double convertNumber<double>(const mpq_class& x)
{
  return x.get_d();
}

IPO_EXPORT
mpq_class reconstructRational(double x, double maxError = 1.0e-15);

IPO_EXPORT
void reconstructRational(mpq_ptr result, double x, double maxError = 1.0e-15);

template<>
inline mpq_class convertNumber<mpq_class>(const double& x)
{
  return reconstructRational(x);
}

namespace ipo
{

  class IntegralScaler
  {
  public:
    IntegralScaler();

    void operator()(const mpq_class& x);

    const mpq_class& factor() const;

  private:
    mpq_class _factor;
  };

}

#endif /* IPO_WITH_GMP */

