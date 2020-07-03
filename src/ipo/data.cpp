#include <ipo/data.hpp>

namespace ipo
{
  
  Value minusInfinity()
  {
    return Value(-std::numeric_limits<double>::infinity());
  }

  Value plusInfinity()
  {
    return Value(std::numeric_limits<double>::infinity());
  }
  
  inline std::size_t* getFirstCoordinate(Vector::Header* data)
  {
    return (std::size_t*)((char*)(data) + sizeof(Vector::Header));
  }
  
  inline double* getFirstReal(Vector::Header* data)
  {
    return (double*)((char*)(data) + sizeof(Vector::Header)
      + sizeof(std::size_t) * data->size);
  }

  inline std::size_t calculateMemoryReal(std::size_t size)
  {
    return (std::size_t)( sizeof(Vector::Header)
      + (sizeof(std::size_t) + sizeof(double)) * size );
  }

#if defined(IPO_WITH_GMP)

  inline mpq_class* getFirstRational(Vector::Header* data)
  {
    return (mpq_class*)((char*)(data) + sizeof(Vector::Header)
      + (sizeof(std::size_t) + sizeof(double)) * data->size);
  }

  inline std::size_t calculateMemoryRational(std::size_t size)
  {
    return (std::size_t)( sizeof(Vector::Header)
      + (sizeof(std::size_t) + sizeof(double) + sizeof(mpq_class)) * size );
  }

#endif /* IPO_WITH_GMP */

  inline bool checkSorted(Vector::Header* data)
  {
    auto firstCoordinate = getFirstCoordinate(data);
    for (std::size_t i = 1; i < data->size; ++i)
    {
      if (firstCoordinate[i-1] >= firstCoordinate[i])
        return false;
    }
    return true;
  }

  Vector::~Vector()
  {
    assert(_data != nullptr);
    _data->usage--;
    if (_data->usage > 0)
      return;

    auto raw = (char*)(_data);
#if defined(IPO_WITH_GMP)
    if (isRational())
    {
      mpq_class* start = (mpq_class*)(raw 
        + sizeof(Header) + (sizeof(std::size_t) + sizeof(double)) * _data->size);
      for (std::size_t i = 0; i < _data->size; ++i)
      {
        start->~mpq_class();
        ++start;
      }
    }
#endif /* IPO_WITH_GMP */
    delete[] raw;
  }

  Vector::Vector(const Vector& other)
  {
    assert(other._data != nullptr);
    _data = other._data;
    _data->usage++;
  }

  Vector& Vector::operator=(const Vector& other)
  {
    assert(_data != nullptr);
    this->~Vector();

    assert(other._data != nullptr);
    _data = other._data;
    _data->usage++;
    return *this;
  }

  Vector::Vector()
  {
    _data = (Header*)(new char[calculateMemoryReal(0)]);
    _data->usage = 1;
    _data->properties = 1;
    _data->size = 0;
  }

  Vector::Vector(const std::vector<std::size_t>& coordinates, const std::vector<double>& reals)
  {
    assert(coordinates.size() == reals.size());
    _data = (Header*)(new char[calculateMemoryReal(coordinates.size())]);
    _data->usage = 1;
    _data->properties = 0;
    _data->size = coordinates.size();
    const auto firstCoordinate = getFirstCoordinate(_data);
    const auto firstReal = getFirstReal(_data);
    for (std::size_t i = 0; i < _data->size; ++i)
      firstCoordinate[i] = coordinates[i];
    for (std::size_t i = 0; i < _data->size; ++i)
      firstReal[i] = reals[i];

    assert(checkSorted(_data));
  }

  Vector::Vector(const std::vector<std::pair<std::size_t, double>>& nonzeros)
  {
    _data = (Header*)(new char[calculateMemoryReal(nonzeros.size())]);
    _data->usage = 1;
    _data->properties = 0;
    _data->size = nonzeros.size();
    const auto firstCoordinate = getFirstCoordinate(_data);
    const auto firstReal = getFirstReal(_data);
    for (std::size_t i = 0; i < _data->size; ++i)
    {
      firstCoordinate[i] = nonzeros[i].first;
      firstReal[i] = nonzeros[i].second;
    }

    assert(checkSorted(_data));
  }

  Vector::Vector(const std::vector<double>& entries)
  {
    std::size_t size = 0;
    for (std::size_t coordinate = 0; coordinate < entries.size(); ++coordinate)
    {
      if (entries[coordinate] != 0.0)
        ++size;
    }

    _data = (Header*)(new char[calculateMemoryReal(size)]);
    _data->usage = 1;
    _data->properties = 0;
    _data->size = size;
    const auto firstCoordinate = getFirstCoordinate(_data);
    const auto firstReal = getFirstReal(_data);
    std::size_t i = 0;
    for (std::size_t coordinate = 0; coordinate < entries.size(); ++coordinate)
    {
      if (entries[coordinate] != 0.0)
      {
        firstCoordinate[i] = coordinate;
        firstReal[i] = entries[coordinate];
        ++i;
      }
    }

    assert(checkSorted(_data));
  }

  std::size_t Vector::coordinate(std::size_t index) const
  {
    return getFirstCoordinate(_data)[index];
  }

  std::size_t Vector::findCoordinate(std::size_t coord) const
  {
    std::size_t left = 0;
    std::size_t right = this->size();
    while (left < right)
    {
      std::size_t mid = (left + right) / 2;
      std::size_t c = coordinate(mid);
      if (c == coord)
        return mid;
      else if (c < coord)
        left = mid + 1;
      else
        right = mid;
    }
    return std::numeric_limits<std::size_t>::max();
  }

  double Vector::real(std::size_t index) const
  {
    return getFirstReal(_data)[index];
  }

#if defined(IPO_WITH_GMP)

  Vector::Vector(const std::vector<std::size_t>& coordinates,
    const std::vector<mpq_class>& rationals)
  {
    assert(coordinates.size() == rationals.size());
    _data = (Header*)(new char[calculateMemoryRational(coordinates.size())]);
    _data->usage = 1;
    _data->properties = 1;
    _data->size = coordinates.size();
    auto firstCoordinate = getFirstCoordinate(_data);
    auto firstReal = getFirstReal(_data);
    auto firstRational = getFirstRational(_data);
    for (std::size_t i = 0; i < _data->size; ++i)
      firstCoordinate[i] = coordinates[i];
    for (std::size_t i = 0; i < _data->size; ++i)
    {
      mpq_class* x = new (&firstRational[i]) mpq_class(rationals[i]);
      firstReal[i] = x->get_d();
    }

    assert(checkSorted(_data));
  }

  Vector::Vector(const std::vector<std::size_t>& coordinates,
    std::vector<mpq_class>& rationals, bool swap)
  {
    assert(coordinates.size() == rationals.size());
    _data = (Header*)(new char[calculateMemoryRational(coordinates.size())]);
    _data->usage = 1;
    _data->properties = 1;
    _data->size = coordinates.size();
    auto firstCoordinate = getFirstCoordinate(_data);
    auto firstReal = getFirstReal(_data);
    auto firstRational = getFirstRational(_data);
    for (std::size_t i = 0; i < _data->size; ++i)
      firstCoordinate[i] = coordinates[i];
    for (std::size_t i = 0; i < _data->size; ++i)
    {
      mpq_class* x = new (&firstRational[i]) mpq_class();
      if (swap)
        x->swap(rationals[i]);
      else
        *x = rationals[i];
      firstReal[i] = x->get_d();
    }

    assert(checkSorted(_data));
  }

  Vector::Vector(const std::vector<std::pair<std::size_t, mpq_class>>& nonzeros)
  {
    _data = (Header*)(new char[calculateMemoryRational(nonzeros.size())]);
    _data->usage = 1;
    _data->properties = 1;
    _data->size = nonzeros.size();
    auto firstCoordinate = getFirstCoordinate(_data);
    auto firstReal = getFirstReal(_data);
    auto firstRational = getFirstRational(_data);
    for (std::size_t i = 0; i < _data->size; ++i)
    {
      firstCoordinate[i] = nonzeros[i].first;
      mpq_class* x = new (&firstRational[i]) mpq_class(nonzeros[i].second);
      firstReal[i] = x->get_d();
    }

    assert(checkSorted(_data));
  }

  Vector::Vector(std::vector<std::pair<std::size_t, mpq_class>>& nonzeros, bool swap)
  {
    _data = (Header*)(new char[calculateMemoryRational(nonzeros.size())]);
    _data->usage = 1;
    _data->properties = 1;
    _data->size = nonzeros.size();
    auto firstCoordinate = getFirstCoordinate(_data);
    auto firstReal = getFirstReal(_data);
    auto firstRational = getFirstRational(_data);
    for (std::size_t i = 0; i < _data->size; ++i)
    {
      firstCoordinate[i] = nonzeros[i].first;
      mpq_class* x = new (&firstRational[i]) mpq_class();
      if (swap)
        x->swap(nonzeros[i].second);
      else
        *x = nonzeros[i].second;
      firstReal[i] = x->get_d();
    }

    assert(checkSorted(_data));
  }

  Vector::Vector(const std::vector<mpq_class>& entries)
  {
    std::size_t size = 0;
    for (std::size_t i = 0; i < entries.size(); ++i)
    {
      if (entries[i] != 0)
        ++size;
    }

    _data = (Header*)(new char[calculateMemoryRational(size)]);
    _data->usage = 1;
    _data->properties = 1;
    _data->size = size;
    auto firstCoordinate = getFirstCoordinate(_data);
    auto firstReal = getFirstReal(_data);
    auto firstRational = getFirstRational(_data);
    std::size_t i = 0;
    for (std::size_t coordinate = 0; coordinate < entries.size(); ++coordinate)
    {
      if (entries[coordinate] != 0)
      {
        firstCoordinate[i] = coordinate;
        mpq_class* x = new (&firstRational[i]) mpq_class(entries[coordinate]);
        firstReal[i] = x->get_d();
        ++i;
      }
    }

    assert(checkSorted(_data));
  }

  Vector::Vector(mpq_t* entries, std::size_t numEntries, bool swap)
  {
    std::size_t size = 0;
    for (std::size_t coordinate = 0; coordinate < numEntries; ++coordinate)
    {
      if (mpq_sgn(entries[coordinate]) != 0)
        ++size;
    }

    _data = (Header*)(new char[calculateMemoryRational(size)]);
    _data->usage = 1;
    _data->properties = 1;
    _data->size = size;
    auto firstCoordinate = getFirstCoordinate(_data);
    auto firstReal = getFirstReal(_data);
    auto firstRational = getFirstRational(_data);
    std::size_t i = 0;
    for (std::size_t coordinate = 0; coordinate < numEntries; ++coordinate)
    {
      if (mpq_sgn(entries[coordinate]) != 0)
      {
        firstCoordinate[i] = coordinate;
        mpq_class* x = new (&firstRational[i]) mpq_class();
        if (swap)
          mpq_swap(x->get_mpq_t(), entries[coordinate]);
        else
          mpq_set(x->get_mpq_t(), entries[coordinate]);
        firstReal[i] = x->get_d();
        ++i;
      }
    }

    assert(checkSorted(_data));
  }

  Vector::Vector(std::vector<mpq_class>& entries, bool swap)
  {
    std::size_t size = 0;
    for (std::size_t i = 0; i < entries.size(); ++i)
    {
      if (entries[i] != 0)
        ++size;
    }

    _data = (Header*)(new char[calculateMemoryRational(size)]);
    _data->usage = 1;
    _data->properties = 1;
    _data->size = size;
    auto firstCoordinate = getFirstCoordinate(_data);
    auto firstReal = getFirstReal(_data);
    auto firstRational = getFirstRational(_data);
    std::size_t i = 0;
    for (std::size_t coordinate = 0; coordinate < entries.size(); ++coordinate)
    {
      if (entries[coordinate] != 0)
      {
        firstCoordinate[i] = coordinate;
        mpq_class* x = new (&firstRational[i]) mpq_class();
        if (swap)
          x->swap(entries[coordinate]);
        else
          *x = entries[coordinate];
        firstReal[i] = x->get_d();
        ++i;
      }
    }

    assert(checkSorted(_data));
  }

  const mpq_class& Vector::rational(std::size_t index) const
  {
    auto firstRational = getFirstRational(_data);
    return firstRational[index];
  }

#endif /* IPO_WITH_GMP */

  double Vector::squaredRealNorm() const
  {
    double result = 0.0;
    for (std::size_t p = 0; p < this->size(); ++p)
    {
      double x = this->real(p);
      result += x*x;
    }
    return result;
  }

#if defined(IPO_WITH_GMP)

  mpq_class Vector::squaredRationalNorm() const
  {
    mpq_class result = 0.0;
    for (std::size_t p = 0; p < this->size(); ++p)
    {
      const mpq_class& x = this->rational(p);
      result += x*x;
    }
    return result;
  }

#endif /* IPO_WITH_GMP */

  Value operator*(const Vector& a, const Vector& b)
  {
    Value result;
    if (a.size() == 0 || b.size() == 0)
      return result;

#if !defined(IPO_WITH_GMP)
    assert(!a.isRational() || !b.isRational());
#endif /* !IPO_WITH_GMP */
    
    if (a.isRational() && b.isRational())
    {
#if defined(IPO_WITH_GMP)
      std::size_t pa = 0;
      std::size_t pb = 0;
      std::size_t ia = a.coordinate(pa);
      std::size_t ib = b.coordinate(pb);
      while (true)
      {
        if (ia < ib)
        {
          ++pa;
          if (pa == a.size())
            break;
          ia = a.coordinate(pa);
        }
        else if (ia > ib)
        {
          ++pb;
          if (pb == b.size())
            break;
          ib = b.coordinate(pb);
        }
        else
        {
          result.rational += a.rational(pa) * b.rational(pb);
          ++pa;
          ++pb;
          if (pa == a.size() || pb == b.size())
            break;
          ia = a.coordinate(pa);
          ib = b.coordinate(pb);
        }
      }
      result.real = result.rational.get_d();
#endif /* IPO_WITH_GMP */
    }
    else
    {
      std::size_t pa = 0;
      std::size_t pb = 0;
      std::size_t ia = a.coordinate(pa);
      std::size_t ib = b.coordinate(pb);
      while (true)
      {
        if (ia < ib)
        {
          ++pa;
          if (pa == a.size())
            break;
          ia = a.coordinate(pa);
        }
        else if (ia > ib)
        {
          ++pb;
          if (pb == b.size())
            break;
          ib = b.coordinate(pb);
        }
        else
        {
          result.real += a.real(pa) * b.real(pb);
          ++pa;
          ++pb;
          if (pa == a.size() || pb == b.size())
            break;
          ia = a.coordinate(pa);
          ib = b.coordinate(pb);
        }
      }
#if defined(IPO_WITH_GMP)
      result.rational = result.real;
#endif /* IPO_WITH_GMP */
    }
    return result;
  }

#if defined(IPO_WITH_GMP)
  
  mpq_class operator*(const ipo::Vector& a, const mpq_class* b)
  {
    mpq_class product = 0;
    for (std::size_t i = 0; i < a.size(); ++i)
      product += a.real(i) * b[a.coordinate(i)];
    return product;
  }

#endif /* IPO_WITH_GMP */

  double operator*(const ipo::Vector& a, const double* b)
  {
    double product = 0;
    for (std::size_t i = 0; i < a.size(); ++i)
      product += a.real(i) * b[a.coordinate(i)];
    return product;
  }

  Constraint alwaysSatisfiedConstraint()
  {
    return Constraint(minusInfinity(), Vector(), Value());
  }

  Constraint neverSatisfiedConstraint()
  {
#if defined(IPO_WITH_GMP)
    return Constraint(minusInfinity(), Vector(), Value(mpq_class(-1)));
#else
    return Constraint(minusInfinity(), Vector(), Value(-1.0));
#endif /* IPO_WITH_GMP */
  }

} /* namespace ipo */
