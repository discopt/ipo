#include "vectors.h"

#include <algorithm>
#include <iostream>
#include <cmath>

namespace ipo {

  VectorData::Nonzero::Nonzero(std::size_t idx, const Rational& val)
    : index(idx), value(val)
  {
    approximation = (double)val;
  }

  VectorData::Nonzero::~Nonzero()
  {

  }

  bool VectorData::Nonzero::operator<(const VectorData::Nonzero& other) const
  {
    return index < other.index;
  }

  VectorData::VectorData()
    : _mutableUsage(0), _immutableUsage(0)
  {

  }

  VectorData::VectorData(std::size_t initialMemory)
    : _mutableUsage(0), _immutableUsage(0)
  {
    _nonzeros.reserve(initialMemory);
  }

  VectorData::VectorData(const VectorData& other)
    : _nonzeros(other._nonzeros), _mutableUsage(0), _immutableUsage(0)
  {

  }

  VectorData::~VectorData()
  {

  }

  VectorData& VectorData::operator=(const VectorData& other)
  {
    _nonzeros = other._nonzeros;
  }

  bool VectorData::operator==(const VectorData& other) const
  {
    if (_nonzeros.size() != other._nonzeros.size())
      return false;

    for (std::size_t i = 0; i < _nonzeros.size(); ++i)
    {
      if (_nonzeros[i].index != other._nonzeros[i].index)
        return false;
      if (_nonzeros[i].approximation != other._nonzeros[i].approximation)
        return false;
    }

    for (std::size_t i = 0; i < _nonzeros.size(); ++i)
    {
      if (_nonzeros[i].value != other._nonzeros[i].value)
        return false;
    }

    return true;
  }

  bool VectorData::operator<(const VectorData& other) const
  {
    if (_nonzeros.size() < other._nonzeros.size())
      return true;
    else if (_nonzeros.size() > other._nonzeros.size())
      return false;

    for (std::size_t i = 0; i < _nonzeros.size(); ++i)
    {
      if (_nonzeros[i].index < other._nonzeros[i].index)
        return true;
      if (_nonzeros[i].index > other._nonzeros[i].index)
        return false;
      if (_nonzeros[i].approximation < other._nonzeros[i].approximation)
        return true;
      if (_nonzeros[i].approximation > other._nonzeros[i].approximation)
        return false;
    }

    for (std::size_t i = 0; i < _nonzeros.size(); ++i)
    {
      if (_nonzeros[i].value < other._nonzeros[i].value)
        return true;
      if (_nonzeros[i].value > other._nonzeros[i].value)
        return false;
    }

    return false;
  }

  void VectorData::add(std::size_t index, const Rational& value)
  {
    assert(_immutableUsage == 0);
    assert(value != 0);
    _nonzeros.push_back(Nonzero(index, value));
  }

  void VectorData::sort()
  {
    assert(_immutableUsage == 0);
    std::sort(_nonzeros.begin(), _nonzeros.end());
  }

  bool VectorData::isSorted() const
  {
    if (_nonzeros.size() <= 1)
      return true;

    std::size_t max = _nonzeros.size() - 1;
    for (std::size_t i = 0; i < max; ++i)
    {
      if (_nonzeros[i].index >= _nonzeros[i+1].index)
        return false;
    }
    return true;
  }

  void VectorData::deleteIfUnused()
  {
    if (_mutableUsage == 0 && _immutableUsage == 0)
      delete this;
  }

  void VectorData::scale(int scalar)
  {
    assert(_immutableUsage == 0);
    if (scalar != 0)
    {
      for (std::size_t i = 0; i < _nonzeros.size(); ++i)
        _nonzeros[i].value *= scalar;
    }
    else
    {
      _nonzeros.clear();
    }
  }

  static VectorData* zeroVectorData = NULL;

  ReferenceCountedVector::ReferenceCountedVector()
  {
    if (!zeroVectorData)
    {
      zeroVectorData = new VectorData(0);
      zeroVectorData->_immutableUsage++;
    }
    _data = zeroVectorData;
  }

  MutableVector zeroVector()
  {
    return MutableVector(new VectorData());
  }

  MutableVector unitVector(std::size_t index)
  {
    VectorData* data = new VectorData(1);
    data->add(index, Rational(1));
    return MutableVector(data);
  }

  Rational operator*(const ReferenceCountedVector& a, const ReferenceCountedVector& b)
  {
    if (a.size() == 0 || b.size() == 0)
      return 0;

    Rational result = 0;
    std::size_t pa = 0;
    std::size_t pb = 0;
    std::size_t ia = a.index(pa);
    std::size_t ib = b.index(pb);
    while (true)
    {
      if (ia < ib)
      {
        ++pa;
        if (pa == a.size())
          break;
        ia = a.index(pa);
      }
      else if (ia > ib)
      {
        ++pb;
        if (pb == b.size())
          break;
        ib = b.index(pb);
      }
      else
      {
        result += a.value(pa) * b.value(pb);
        ++pa;
        ++pb;
        if (pa == a.size() || pb == b.size())
          break;
        ia = a.index(pa);
        ib = b.index(pb);
      }
    }
    return result;
  }

  Vector operator+(const Vector& a, const Vector& b)
  {
    return addScaled(a, 1, b, 1);
  }

  Vector operator-(const Vector& a, const Vector& b)
  {
    return addScaled(a, 1, b, -1);
  }

  Vector addScaled(const Vector& a, int scaleA, const Vector& b, int scaleB)
  {
    assert(scaleA != 0);
    assert(scaleB != 0);

    if (scaleA == 1 && b.size() == 0)
      return a;
    if (scaleB == 1 && a.size() == 0)
      return b;

    VectorData* data = new VectorData( a.size() + b.size());
    std::size_t pa = 0;
    std::size_t pb = 0;
    std::size_t ia, ib;
    while (true)
    {
      if (pa < a.size())
      {
        ia = a.index(pa);
        if (pb < b.size())
          ib = b.index(pb);
        else
          ib = std::numeric_limits<std::size_t>::max();
      }
      else
      {
        ia = std::numeric_limits<std::size_t>::max();
        if (pb < b.size())
          ib = b.index(pb);
        else
          break;
      }

      if (ia < ib)
      {
        data->add(ia, scaleA * a.value(pa));
        ++pa;
      }
      else if (ib < ia)
      {
        data->add(ib, scaleB * b.value(pa));
        ++pb;
      }
      else
      {
        Rational x = scaleA * a.value(pa) + scaleB * b.value(pb);
        if (x != 0)
          data->add(ia, x);
        ++pa;
        ++pb;
      }
    }
    return MutableVector(data);
  }

  std::size_t differingIndex(const ReferenceCountedVector& a, const ReferenceCountedVector& b)
  {
    if (a.size() == 0)
    {
      if (b.size() == 0)
        return std::numeric_limits<std::size_t>::max();
      else
        return b.index(0);
    }
    else if (b.size() == 0)
      return a.index(0);

    std::size_t pa = 0;
    std::size_t pb = 0;
    std::size_t ia = a.index(pa);
    std::size_t ib = b.index(pb);

    // Quick scan for comparing support structure and approximate values.

    while (true)
    {
      if (ia < ib)
        return ia;
      else if (ia > ib)
        return ib;
      else
      {
        if (a.approximation(pa) != b.approximation(pb))
          return ia;
        ++pa;
        ++pb;
        if (pa == a.size() && pb == b.size())
          break;
        ia = a.index(pa);
        ib = b.index(pb);

      }
    }

    // We now have same support, so we can compare exact values easier.

    for (std::size_t p = 0; p < a.size(); ++p)
    {
      if (a.value(p) != b.value(p))
        return a.index(p);
    }

    return std::numeric_limits<std::size_t>::max();
  }

  void vectorToDense(const ReferenceCountedVector& source, soplex::VectorRational& target)
  {
    assert(source.size() == 0 || source.index(source.size()-1) < target.dim());
    target.clear();
    for (std::size_t p = 0; p < source.size(); ++p)
      target[source.index(p)] = source.value(p);
  }

  void vectorToSparse(const ReferenceCountedVector& source, soplex::SVectorRational& target)
  {
    target.clear();
    for (std::size_t p = 0; p < source.size(); ++p)
      target.add(source.index(p), source.value(p));
  }

  MutableVector denseToVector(const soplex::VectorRational& source, bool saveMemory)
  {
    std::size_t mem;
    if (saveMemory)
    {
      mem = 0;
      for (std::size_t i = 0; i < source.dim(); ++i)
      {
        if (source[i] != 0)
          ++mem;
      }
    }
    else
      mem = source.dim();

    VectorData* data = new VectorData(mem);
    for (std::size_t i = 0; i < source.dim(); ++i)
    {
      if (source[i] != 0)
        data->add(i, source[i]);
    }
    return MutableVector(data);
  }

  MutableVector sparseToVector(const soplex::SVectorRational& source)
  {
    VectorData* data = new VectorData(source.size());
    for (int p = source.size() - 1; p >= 0; --p)
      data->add(source.index(p), source.value(p));
    data->sort();
    return MutableVector(data);
  }

  Rational operator*(const soplex::VectorRational& a, const ReferenceCountedVector& b)
  {
    assert(b.size() == 0 || b.index(b.size()-1) < a.dim());

    Rational result = 0;
    for (std::size_t p = 0; p < b.size(); ++p)
      result += a[b.index(p)] * b.value(p);
    return result;
  }

  Rational operator*(const ReferenceCountedVector& a, const soplex::VectorRational& b)
  {
    return b * a;
  }

  soplex::VectorRational& operator+=(soplex::VectorRational& a, const ReferenceCountedVector& b)
  {
    assert(b.size() == 0 || b.index(b.size() - 1) < a.dim());

    for (std::size_t p = 0; p < b.size(); ++p)
      a[b.index(p)] += b.value(p);

    return a;
  }

  soplex::VectorRational& operator-=(soplex::VectorRational& a, const ReferenceCountedVector& b)
  {
    assert(b.size() == 0 || b.index(b.size() - 1) < a.dim());

    for (std::size_t p = 0; p < b.size(); ++p)
      a[b.index(p)] -= b.value(p);

    return a;
  }

  MutableVector operator-(const ReferenceCountedVector& source)
  {
    VectorData* data = new VectorData(source.size());
    for (std::size_t p = 0; p < source.size(); ++p)
      data->add(source.index(p), -source.value(p));
    return MutableVector(data);
  }

  MutableVector integralScaled(const Vector& vector, Rational* factor)
  {
    // Compute scaling factor.

    IntegralScaler scaler;
    for (std::size_t p = 0; p < vector.size(); ++p)
      scaler(vector.value(p));
    const Rational theFactor = scaler.factor();
    if (factor != NULL)
      *factor = theFactor;

    // Scale it.

    VectorData* data = new VectorData(vector.size());
    for (std::size_t p = 0; p < vector.size(); ++p)
      data->add(vector.index(p), theFactor * vector.value(p));
    return MutableVector(data);
  }

  void scaleIntegral(Vector& vector, Rational* factor)
  {
    vector = integralScaled(vector, factor);
  }

  void scaleIntegral(std::vector<Vector>& vectors)
  {
    for (std::size_t i = 0; i < vectors.size(); ++i)
      scaleIntegral(vectors[i]);
  }

  void scaleIntegral(const soplex::VectorRational& vector, soplex::DVectorRational& scaled)
  {
    // Compute scaling factor.

    IntegralScaler scaler;
    for (std::size_t i = 0; i < vector.dim(); ++i)
      scaler(vector[i]);
    Rational factor = scaler.factor();

    // Scale it.

    scaled.reDim(vector.dim(), false);
    for (std::size_t i = 0; i < vector.dim(); ++i)
      scaled[i] = factor * vector[i];
  }

} /* namespace ipo */
