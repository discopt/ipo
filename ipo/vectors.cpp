#include "vectors.h"

#include <iostream>

namespace ipo {
  
  bool SparseVector::Nonzero::operator<(const SparseVector::Nonzero& other) const
  {
    return index < other.index;
  }
  
  SparseVector::SparseVector()
    : impl(NULL)
  {

  }

  SparseVector::SparseVector(std::size_t capacity)
  {
    impl = new Implementation;
    impl->usage = 1;
    impl->capacity = capacity;
    impl->size = 0;
    impl->nonzeros = new Nonzero[capacity];
  }
  
  SparseVector::SparseVector(const SparseVector& other)
  {
    if (other.impl != NULL)
    {
      impl = const_cast<Implementation*>(other.impl);
      ++impl->usage;
    }
    else
      impl = NULL;
  }

  SparseVector::~SparseVector()
  {
    if (impl != NULL)
    {
      --impl->usage;
      if (impl->usage == 0)
      {
        delete[] impl->nonzeros;
        delete impl;
      }
    }
  }
  
  SparseVector& SparseVector::operator=(const SparseVector& other)
  {
    if (impl != NULL)
    {
      --impl->usage;
      if (impl->usage == 0)
      {
        delete[] impl->nonzeros;
        delete impl;
      }
    }
    if (other.impl != NULL)
    {
      impl = const_cast<Implementation*>(other.impl);
      ++impl->usage;
    }
    else
      impl = NULL;
  }

  bool SparseVector::operator==(const SparseVector& other) const
  {
    if (impl == other.impl)
      return true;
    else if (impl->size != other.impl->size)
      return false;

    std::size_t s = impl->size;
    const Nonzero* a = impl->nonzeros;
    const Nonzero* b = other.impl->nonzeros;
    for (std::size_t p = 0; p < s; ++p)
    {
      if (a[p].index != b[p].index)
        return false;;
      if (a[p].approximation != b[p].approximation)
        return false;
    }
    for (std::size_t p = 0; p < s; ++p)
    {
      if (a[p].value != b[p].value)
        return false;
    }
    return true;
  }

  bool SparseVector::operator<(const SparseVector& other) const
  {
    if (impl == other.impl)
      return false;

    if (impl->size < other.impl->size)
      return true;
    else if (impl->size > other.impl->size)
      return false;

    std::size_t s = impl->size;
    const Nonzero* a = impl->nonzeros;
    const Nonzero* b = other.impl->nonzeros;
    for (std::size_t p = 0; p < s; ++p)
    {
      if (a[p].index < b[p].index)
        return true;
      if (a[p].index > b[p].index)
        return false;
      if (a[p].approximation < b[p].approximation)
        return true;
      if (a[p].approximation > b[p].approximation)
        return false;
    }
    for (std::size_t p = 0; p < s; ++p)
    {
      if (a[p].value < b[p].value)
        return true;
      if (a[p].value > b[p].value)
        return false;
    }
    return false;
  }

  void SparseVector::add(std::size_t index, const Rational& value)
  {
    assert(impl->size < impl->capacity);
    Nonzero& nz = impl->nonzeros[impl->size];
    nz.index = index;
    nz.value = value;
    nz.approximation = double(value);
    ++impl->size;
  }

  void SparseVector::sort()
  {
    std::sort(&impl->nonzeros[0], &impl->nonzeros[impl->size]);
  }

  bool SparseVector::isSorted() const
  {
    for (std::size_t p = 0; p + 1 < impl->size; ++p)
    {
      if (impl->nonzeros[p].index > impl->nonzeros[p+1].index)
        return false;
    }
    return true;
  }

  void SparseVector::swap(SparseVector& other)
  {
    std::swap(impl, other.impl);
  }

  SparseVector denseToSparseVector(const DenseVector& source)
  {
    std::size_t size = 0;
    for (std::size_t i = 0; i < source.dim(); ++i)
    {
      if (source[i] != 0)
        ++size;
    }
    SparseVector result(size);
    for (std::size_t i = 0; i < source.dim(); ++i)
    {
      if (source[i] != 0)
        result.add(i, source[i]);
    }
    return result;
  }

  void assign(DenseVectorApproximation& target, const DenseVector& source)
  {
    target = source;
  }

  void assign(DenseVector& target, const DenseVector& source)
  {
    target = source;
  }

  void assign(DenseVector& target, const SparseVector& source)
  {
    target.clear();
    assert(source.size() == 0 || source.index(source.size() - 1) < target.dim());
    for (std::size_t p = 0; p < source.size(); ++p)
      target[source.index(p)] = source.value(p);
  }

  void assign(SparseVector& target, const DenseVector& source)
  {
    std::size_t size = 0;
    for (std::size_t i = 0; i < source.dim(); ++i)
    {
      if (source[i] != 0)
        ++size;
    }
    target = SparseVector(size);
    for (std::size_t i = 0; i < source.dim(); ++i)
    {
      if (source[i] != 0)
        target.add(i, source[i]);
    }
  }

  void assign(SparseVector& target, const SparseVector& source)
  {
    target = source;
  }

  Rational scalarProduct(const DenseVector& a, const DenseVector& b)
  {
    assert(a.dim() == b.dim());
    return a * b;
  }

  Rational scalarProduct(const DenseVector& a, const SparseVector& b)
  {
    Rational result = 0;
    for (std::size_t p = 0; p < b.size(); ++p)
    {
      assert(b.index(p) < a.dim());
      result += a[b.index(p)] * b.value(p);
    }
    return result;
  }

  Rational scalarProduct(const SparseVector& a, const DenseVector& b)
  {
    return scalarProduct(b, a);
  }

  Rational scalarProduct(const SparseVector& a, const SparseVector& b)
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
  
  std::size_t differingIndex(const SparseVector& a, const SparseVector& b)
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
    while (true)
    {
      if (ia < ib)
        return ia;
      else if (ia > ib)
        return ib;
      else
      {
        if (a.value(pa) != b.value(pb))
          return ia;
        ++pa;
        ++pb;
        if (pa == a.size() && pb == b.size())
          return std::numeric_limits<std::size_t>::max();
        else if (pa == a.size())
          return b.index(pb);
        else if (pb == b.size())
          return a.index(pa);
        ia = a.index(pa);
        ib = b.index(pb);
      }
    }
  }

  void add(DenseVector& target, const SparseVector& source)
  {
    assert(source.size() == 0 || source.index(source.size()-1) < target.dim());
    for (std::size_t p = 0; p < source.size(); ++p)
      target[source.index(p)] += source.value(p);
  }

  void subtract(DenseVector& target, const SparseVector& source)
  {
    assert(source.size() == 0 || source.index(source.size()-1) < target.dim());
    for (std::size_t p = 0; p < source.size(); ++p)
      target[source.index(p)] -= source.value(p);
  }

  UniqueSparseVectors::Data::Data(double hsh, SparseVector& vec)
    : hash(hsh), vector(vec)
  {

  }

  UniqueSparseVectors::Data::Data(const UniqueSparseVectors::Data& other)
    : hash(other.hash), vector(const_cast<SparseVector&>(other.vector))
  {

  }


  UniqueSparseVectors::Data::~Data()
  {

  }

  bool UniqueSparseVectors::Data::operator<(const UniqueSparseVectors::Data& other) const
  {
    if (hash < other.hash)
      return true;
    else if (hash > other.hash)
      return false;
    return vector < other.vector;
  }

  UniqueSparseVectors::UniqueSparseVectors(std::size_t ambientDimension)
    : _hashVector(ambientDimension)
  {
    std::default_random_engine generator(0);
    std::normal_distribution<double> distribution;
    DenseVectorApproximation randomVector(ambientDimension);
    double norm = 0.0;
    while (norm == 0.0)
    {
      bool zero = false;
      for (std::size_t i = 0; i < ambientDimension; ++i)
      {
        double x = distribution(generator);
        if (fabs(x) < 1.0e-6)
        {
          zero = true;
          break;
        }
        randomVector[i] = x;
        norm += x*x;
      }
      if (zero)
        norm = 0.0;
    }

    norm = std::sqrt(norm);
    for (std::size_t i = 0; i < ambientDimension; ++i)
      _hashVector[i] = randomVector[i] / norm;
  }

  UniqueSparseVectors::~UniqueSparseVectors()
  {

  }
  
  void UniqueSparseVectors::clear()
  {
    _data.clear();
  }

  bool UniqueSparseVectors::insert(SparseVector& vector)
  {
    Data data(0.0, vector);
    for (std::size_t p = 0; p < vector.size(); ++p)
      data.hash += _hashVector[vector.index(p)] * vector.approximation(p);

    std::pair<std::set<Data>::iterator, bool> inserted = _data.insert(data);
    if (inserted.second)
    {
      return true;
    }
    else
    {
      vector = const_cast<SparseVector&>(inserted.first->vector);
      return false;
    }
  }
  
} /* namespace ipo */