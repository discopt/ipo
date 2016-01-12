#ifndef IPO_ROWS_H_
#define IPO_ROWS_H_

#include "ipo.h"
#include "spx_gmp.h"

namespace ipo {

  void scaleVectorIntegral(const soplex::VectorRational& input, soplex::DVectorRational& result);
  void scaleRowIntegral(soplex::LPRowRational& row);
  void scaleRowsIntegral(soplex::LPRowSetRational& rows);

} /* namespace ipo */

#endif /* IPO_ROWS_H_ */
