#include <ipo/rational.hpp>

namespace ipo
{

//   rational minusInfinity()
//   {
//     return rational(-std::numeric_limits<double>::infinity());
//   }

//   rational plusInfinity()
//   {
//     return rational(std::numeric_limits<double>::infinity());
//   }

  std::ostream& operator<<(std::ostream& stream, const rational& x)
  {
//     if (x.isFinite())
    stream << x.get_mpq_class();
//     else if (x.isPlusInfinity())
//       stream << "inf";
//     else if (x.isMinusInfinity())
//       stream << "-inf";
//     else
//     {
//       assert(x.isNAN());
//       stream << "nan";
//     }
    return stream;
  }

}
