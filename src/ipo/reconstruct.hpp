#pragma once

#include <ipo/config.hpp>
#include <ipo/export.hpp>

#include <gmpxx.h>

namespace ipo
{

  void mpq_reconstruct(mpq_t& result, double x, double maxError = 1.0e-12);
  
  mpq_class reconstruct(double x, double maxError = 1.0e-12);
  
}
