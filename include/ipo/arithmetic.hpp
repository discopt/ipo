#pragma once

#include <ipo/config.hpp>
#include <ipo/export.hpp>

#if defined(IPO_WITH_GMP)
#include <gmpxx.h>
#endif /* IPO_WITH_GMP */

namespace ipo
{
  struct DoubleIsZero
  {
    double epsilon;

    IPO_EXPORT
    DoubleIsZero(double eps);

    IPO_EXPORT
    bool operator()(double value) const;
  };

  struct RationalIsZero
  {
    IPO_EXPORT
    RationalIsZero() = default;

    IPO_EXPORT
    bool operator()(const mpq_class& x) const
    {
      return x == 0;
    }
  };

  IPO_EXPORT
  inline
  double toDouble(double x)
  {
    return x;
  }

  IPO_EXPORT
  inline
  double toDouble(const mpq_class& x)
  {
    return x.get_d();
  }

  void mpq_reconstruct(mpq_t& result, double x, double maxError = 1.0e-12);

  mpq_class reconstruct(double x, double maxError = 1.0e-12);
    
}
