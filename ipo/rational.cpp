#include "rational.h"

namespace ipo {

  std::string rationalToString(const Rational& number)
  {
    std::stringstream stream;
    stream << number;
    return stream.str();
  }

  bool isIntegral(const Rational& number)
  {
    mpq_class y(number.getMpqRef());
    return y.get_den() == 1;
  }

  IntegralScaler::IntegralScaler()
    : _numScaler(0), _denScaler(1)
  {

  }

  IntegralScaler::~IntegralScaler()
  {

  }

  void IntegralScaler::operator()(const Rational& number)
  {
    const mpq_t& base = number.getMpqRef();
    mpz_gcd(_numScaler.get_mpz_t(), _numScaler.get_mpz_t(), mpq_numref(base));
    mpz_lcm(_denScaler.get_mpz_t(), _denScaler.get_mpz_t(), mpq_denref(base));
  }

  const Rational IntegralScaler::factor() const
  {
    mpq_t s;
    mpq_init(s);
    mpq_set_den(s, _numScaler.get_mpz_t());
    mpq_set_num(s, _denScaler.get_mpz_t());
    Rational result(s);
    mpq_clear(s);
    return result;
  }

} /* namespace ipo */
