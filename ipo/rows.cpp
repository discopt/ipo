#include "rows.h"

#include <gmpxx.h>

#include "unique_rational_vectors.h"
#include "min_norm_2d.h"

using namespace soplex;

namespace ipo {

  void scaleVectorIntegral(const soplex::VectorRational& input, soplex::DVectorRational& result)
  {
    result.reDim(input.dim());
    mpz_class numScaler = 0;
    mpz_class denScaler = 1;
    for (int i = 0; i < input.dim(); ++i)
      integralScaleScan(numScaler, denScaler, input[i]);

    /// Scale it.

    for (int i = 0; i < input.dim(); ++i)
      result[i] = integralScaleUpdate(numScaler, denScaler, input[i]);
  }

  void scaleRowIntegral(LPRowRational& row)
  {
    bool hasLhs = row.lhs() > -infinity;
    bool hasRhs = row.rhs() < infinity;

    /// Scan to find scaler.

    mpz_class numScaler = 0;
    mpz_class denScaler = 1;
    if (hasLhs)
      integralScaleScan(numScaler, denScaler, row.lhs());
    if (hasRhs)
      integralScaleScan(numScaler, denScaler, row.rhs());
    const SVectorRational& vector = row.rowVector();
    for (int p = vector.size() - 1; p >= 0; --p)
      integralScaleScan(numScaler, denScaler, vector.value(p));

    /// Scale the row.

    if (numScaler != 0 && (numScaler != 1 || denScaler != 1))
    {
      if (hasLhs)
        row.setLhs(integralScaleUpdate(numScaler, denScaler, row.lhs()));
      if (hasRhs)
        row.setRhs(integralScaleUpdate(numScaler, denScaler, row.rhs()));
      DSVectorRational vector;
      vector = row.rowVector();
      for (int p = vector.size() - 1; p >= 0; --p)
        integralScaleUpdate(numScaler, denScaler, vector.value(p));
      row.setRowVector(vector);
    }
  }

  void scaleRowsIntegral(LPRowSetRational& rows)
  {
    for (int i = 0; i < rows.num(); ++i)
    {
      bool hasLhs = rows.lhs(i) > -infinity;
      bool hasRhs = rows.rhs(i) < infinity;

      /// Scan to find scaler.

      mpz_class numScaler = 0;
      mpz_class denScaler = 1;
      if (hasLhs)
        integralScaleScan(numScaler, denScaler, rows.lhs(i));
      if (hasRhs)
        integralScaleScan(numScaler, denScaler, rows.rhs(i));
      const SVectorRational& vector = rows.rowVector(i);
      for (int p = vector.size() - 1; p >= 0; --p)
        integralScaleScan(numScaler, denScaler, vector.value(p));

      /// Scale the row.

      if (numScaler != 0 && (numScaler != 1 || denScaler != 1))
      {
        if (hasLhs)
          integralScaleUpdate(numScaler, denScaler, rows.lhs_w(i));
        if (hasRhs)
          integralScaleUpdate(numScaler, denScaler, rows.rhs_w(i));
        SVectorRational& vector = rows.rowVector_w(i);
        for (int p = vector.size() - 1; p >= 0; --p)
          integralScaleUpdate(numScaler, denScaler, vector.value(p));
      }
    }
  }

} /* namespace ipo */
