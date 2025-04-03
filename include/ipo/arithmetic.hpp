#pragma once

#include <cstddef>
#include <sstream>
#include <limits>

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

  std::string formatNumberApprox(const double& x);

  double* generateRandomVectorSphere(std::size_t size);

  double maxAbsoluteValue(const double* vector, std::size_t length);

} /* namespace ipo */

#if defined(IPO_WITH_RATIONAL)

#include <gmp.h>
#include <boost/multiprecision/number.hpp>
#include <boost/multiprecision/gmp.hpp>


namespace ipo
{
  using rational = boost::multiprecision::number<boost::multiprecision::gmp_rational, boost::multiprecision::et_off>;

  class extended_rational
  {
   public:
    extended_rational()
      : _approximation(std::numeric_limits<double>::signaling_NaN())
    {

    }

    /**
     * \brief Copy constructor.
     */

    extended_rational(const extended_rational& other)
      : _rational(other._rational), _approximation(other._approximation)
    {

    }

    /**
     * \brief Move constructor.
     */

    extended_rational(extended_rational&& other)
      : _rational(std::move(other._rational)), _approximation(other._approximation)
    {

    }

    /**
     * \brief Constructor from rational.
     */

    extended_rational(const rational& other)
      : _rational(other), _approximation(0.0)
    {

    }

    /**
     * \brief Move constructor from rational.
     */

    extended_rational(rational&& other)
      : _rational(std::move(other)), _approximation(0.0)
    {

    }

    /**
     * \brief Constructor from a double.
     *
     * Checks if given number is finite.
     */

    extended_rational(double other)
      : _rational(std::isfinite(other) ? other : 0.0), _approximation(other)
    {

    }

    /**
     * \brief Constructor from an infinite or nan double.
     */

    extended_rational(const rational& rational, double other)
      : _rational(rational), _approximation(other)
    {

    }

    /**
     * \brief Assignment operator.
     */

    extended_rational& operator=(const extended_rational& other)
    {
      _rational = other._rational;
      _approximation = other._approximation;
      return *this;
    }

    /**
     * \brief Assignment from double.
     */

    extended_rational& operator=(double other)
    {
      _rational = std::isfinite(other) ? other : 0.0;
      _approximation = other;
      return *this;
    }

    /**
     * \brief Conversion to double.
     */

    explicit operator double() const
    {
      if (_approximation == 0.0)
        return static_cast<double>(_rational);
      else
        return _approximation;
    }

    /**
     * \brief Conversion to double.
     *
     * May update the internal double approximation.
     */

    explicit operator double()
    {
      if (_approximation == 0.0)
        _approximation = static_cast<double>(_rational);
      return _approximation;
    }

    /**
     * \brief Conversion to rational.
     */

    operator rational() const
    {
      if (std::isfinite(_approximation))
        return _rational;
      else
        throw std::runtime_error("not a rational number");
    }

    /**
     * \brief Addition with another \ref extended_rational.
     */

    extended_rational operator+(const extended_rational& other) const
    {
      if (std::isfinite(_approximation) && std::isfinite(other._approximation))
        return extended_rational(_rational + other._rational);
      else
      {
        assert(!std::isfinite(_approximation + other._approximation));
        return extended_rational(0, _approximation + other._approximation);
      }
    }

    /**
     * \brief Addition with a double.
     */

    extended_rational operator+(double other) const
    {
      if (std::isfinite(_approximation) && std::isfinite(other))
        return extended_rational(_rational + other);
      else
        return extended_rational(0, _approximation + other);
    }

    /**
     * \brief Addition with an int.
     */

    extended_rational operator+(int other) const
    {
      if (std::isfinite(_approximation))
        return extended_rational(_rational + other);
      else
        return extended_rational(0, _approximation + other);
    }

    /**
     * \brief Addition-assignment with another \ref extended_rational.
     */

    extended_rational& operator+=(const extended_rational& other)
    {
      if (std::isfinite(_approximation) && std::isfinite(other._approximation))
      {
        _rational += other._rational;
        _approximation = 0.0;
      }
      else
      {
        _rational = 0;
        _approximation += other._approximation;
        assert(!std::isfinite(_approximation));
      }
      return *this;
    }

    /**
     * \brief Addition-assignment with a double.
     */

    extended_rational& operator+=(double other)
    {
      if (std::isfinite(_approximation) && std::isfinite(other))
      {
        _rational += other;
        _approximation = 0.0;
      }
      else
      {
        _rational = 0;
        _approximation += other;
        assert(!std::isfinite(_approximation));
      }
      return *this;
    }

    /**
     * \brief Subtraction with another \ref extended_rational.
     */

    extended_rational operator-(const extended_rational& other) const
    {
      if (std::isfinite(_approximation) && std::isfinite(other._approximation))
        return extended_rational(_rational - other._rational);
      else
      {
        assert(!std::isfinite(_approximation - other._approximation));
        return extended_rational(0, _approximation - other._approximation);
      }
    }

    /**
     * \brief Subtraction with a double.
     */

    extended_rational operator-(double other) const
    {
      if (std::isfinite(_approximation) && std::isfinite(other))
        return extended_rational(_rational - other);
      else
        return extended_rational(0, _approximation - other);
    }

    /**
     * \brief Subtraction with an int.
     */

    extended_rational operator-(int other) const
    {
      if (std::isfinite(_approximation))
        return extended_rational(_rational - other);
      else
        return extended_rational(0, _approximation);
    }

    /**
     * \brief Unary minus.
     */

    extended_rational operator-() const
    {
      return extended_rational(-_rational, -_approximation);
    }

    /**
     * \brief Subtraction-assignment with another \ref extended_rational.
     */

    extended_rational& operator-=(const extended_rational& other)
    {
      if (std::isfinite(_approximation) && std::isfinite(other._approximation))
      {
        _rational -= other._rational;
        _approximation = 0.0;
      }
      else
      {
        _rational = 0;
        _approximation -= other._approximation;
        assert(!std::isfinite(_approximation));
      }
      return *this;
    }

    /**
     * \brief Subtraction-assignment with a double.
     */

    extended_rational& operator-=(double other)
    {
      if (std::isfinite(_approximation) && std::isfinite(other))
      {
        _rational += other;
        _approximation = 0.0;
      }
      else
      {
        _rational = 0;
        _approximation += other;
        assert(!std::isfinite(_approximation));
      }
      return *this;
    }

    /**
     * \brief Multiplication with another \ref extended_rational.
     */

    extended_rational operator*(const extended_rational& other) const
    {
      if (std::isfinite(_approximation) && std::isfinite(other._approximation))
        return extended_rational(_rational * other._rational);
      else
      {
        assert(!std::isfinite(_approximation * other._approximation));
        return extended_rational(0, _approximation * other._approximation);
      }
    }

    /**
     * \brief Multiplication with a double.
     */

    extended_rational operator*(double other) const
    {
      if (std::isfinite(_approximation) && std::isfinite(other))
        return extended_rational(_rational * other);
      else
        return extended_rational(0, _approximation * other);
    }

    /**
     * \brief Multiplication with an int.
     */

    extended_rational operator*(int other) const
    {
      if (std::isfinite(_approximation))
        return extended_rational(_rational * other);
      else
        return extended_rational(0, _approximation * other);
    }

    /**
     * \brief Multiplication-assignment with another \ref extended_rational.
     */

    extended_rational& operator*=(const extended_rational& other)
    {
      if (std::isfinite(_approximation) && std::isfinite(other._approximation))
      {
        _rational *= other._rational;
        _approximation = 0.0;
      }
      else
      {
        _rational = 0;
        _approximation *= other._approximation;
        assert(!std::isfinite(_approximation));
      }
      return *this;
    }

    /**
     * \brief Multiplication-assignment with a double.
     */

    extended_rational& operator*=(double other)
    {
      if (std::isfinite(_approximation) && std::isfinite(other))
      {
        _rational *= other;
        _approximation = 0.0;
      }
      else
      {
        _rational = 0;
        _approximation *= other;
        assert(!std::isfinite(_approximation));
      }
      return *this;
    }

    /**
     * \brief Division by another \ref extended_rational.
     */

    extended_rational operator/(const extended_rational& other) const
    {
      if (std::isfinite(_approximation) && std::isfinite(other._approximation))
      {
        if (other._rational != 0)
          return extended_rational(_rational / other._rational);
        else
          return extended_rational(0, INFINITY);
      }
      else
      {
        assert(!std::isfinite(_approximation / other._approximation));
        return extended_rational(0, _approximation / other._approximation);
      }
    }

    /**
     * \brief Division by a double.
     */

    extended_rational operator/(double other) const
    {
      if (std::isfinite(_approximation) && std::isfinite(other))
      {
        if (other != 0.0)
          return extended_rational(_rational / other);
        else
          return extended_rational(0, INFINITY);
      }
      else
        return extended_rational(0, _approximation / other);
    }

    /**
     * \brief Multiplication with an int.
     */

    extended_rational operator/(int other) const
    {
      if (std::isfinite(_approximation))
      {
        if (other != 0)
          return extended_rational(_rational / other);
        else
          return extended_rational(0, INFINITY);
      }
      else
        return extended_rational(0, _approximation / other);
    }

    /**
     * \brief Division-assignment by another \ref extended_rational.
     */

    extended_rational& operator/=(const extended_rational& other)
    {
      if (std::isfinite(_approximation) && std::isfinite(other._approximation))
      {
        if (other._rational != 0)
        {
          _rational /= other._rational;
          _approximation = 0.0;
        }
        else
        {
          _rational = 0;
          _approximation = INFINITY;
        }
      }
      else
      {
        _rational = 0;
        _approximation /= other._approximation;
        assert(!std::isfinite(_approximation));
      }
      return *this;
    }

    /**
     * \brief Division-assignment by a double.
     */

    extended_rational& operator/=(double other)
    {
      if (std::isfinite(_approximation) && std::isfinite(other))
      {
        if (other != 0.0)
        {
          _rational /= other;
          _approximation = 0.0;
        }
        else
        {
          _rational = 0;
          _approximation = INFINITY;
        }
      }
      else
      {
        _rational = 0;
        _approximation /= other;
        assert(!std::isfinite(_approximation));
      }
      return *this;
    }

    friend std::ostream& operator<<(std::ostream& stream, const extended_rational& enumber);

  private:


    /// Actual rational number.
    rational _rational;
    /// We use the approximation to indicate +/- infinity and NaN.
    double _approximation;
  };

  /**
   * \brief Write \ref extended_rational to \p stream.
   */

  inline
  std::ostream& operator<<(std::ostream& stream, const extended_rational& x)
  {
    if (std::isfinite(x._approximation))
      return stream << static_cast<double>(x._rational);
    else
      return stream << x._approximation;
  }

}

namespace std
{

  template <>
  struct numeric_limits<ipo::extended_rational>
  {
    static constexpr bool is_specialized = true;
    static constexpr bool is_signed = true;
    static constexpr bool is_integer = false;
    static constexpr bool is_exact = true;
    static constexpr bool has_infinity = true;
    static constexpr bool has_quiet_NaN = true;
    static constexpr bool has_signaling_NaN = true;

    // static constexpr std::float_round_style round_style = std::numeric_limits<int>::round_style;
    static constexpr bool is_iec559 = false;
    static constexpr bool is_bounded = false;
    static constexpr bool is_modulo = false;
    // static constexpr int digits = std::numeric_limits<float>::digits;
    // static constexpr int digits10 = std::numeric_limits<float>::digits10;
    // static constexpr int max_digits10 = std::numeric_limits<float>::max_digits10;
    // static constexpr int radix = std::numeric_limits<float>::radix;
    // static constexpr int min_exponent = std::numeric_limits<float>::min_exponent;
    // static constexpr int min_exponent10 = std::numeric_limits<float>::min_exponent10;
    // static constexpr int max_exponent = std::numeric_limits<float>::max_exponent;
    // static constexpr int max_exponent10 = std::numeric_limits<float>::max_exponent10;
    // static constexpr bool traps = std::numeric_limits<float>::traps;
    // static constexpr bool tinyness_before = std::numeric_limits<float>::tinyness_before;

    static ipo::extended_rational infinity()
    {
      return ipo::extended_rational(std::numeric_limits<double>::infinity());
    }

    static ipo::extended_rational quiet_NaN()
    {
      return ipo::extended_rational(std::numeric_limits<double>::quiet_NaN());
    }

    static ipo::extended_rational signaling_NaN()
    {
      return ipo::extended_rational(std::numeric_limits<double>::signaling_NaN());
    }
  };

} /* namespace std */

namespace ipo
{

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

  std::string formatNumberApprox(const rational& x);

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


#endif /* IPO_WITH_RATIONAL */
