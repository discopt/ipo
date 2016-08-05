#include "vectors.h"

#include <iostream>

namespace ipo {

  SparseVector::SparseVector(std::size_t capacity)
  {
    impl = new Implementation;
    impl->usage = 1;
    impl->capacity = capacity;
    impl->size = 0;
    impl->nonzeros = new Nonzero[capacity];
  }
  
  SparseVector::SparseVector(SparseVector& other)
  {
    impl = other.impl;
    ++impl->usage;
  }

  SparseVector::~SparseVector()
  {
    --impl->usage;
    if (impl->usage == 0)
    {
      delete[] impl->nonzeros;
      delete impl;
    }
  }
  
  SparseVector& SparseVector::operator=(SparseVector& other)
  {
    --impl->usage;
    if (impl->usage == 0)
    {
      delete[] impl->nonzeros;
      delete impl;
    }
    impl = other.impl;
    ++impl->usage;
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

  bool UniqueSparseVectors::insert(SparseVector& vector)
  {
    Data data(0.0, vector);
    for (std::size_t p = 0; p < vector.size(); ++p)
      data.hash += _hashVector[vector.index(p)] * vector.approximation(p);

    std::pair<std::set<Data>::iterator, bool> inserted = _data.insert(data);
    if (!inserted.second)
      vector = const_cast<SparseVector&>(inserted.first->vector);
  }
  
} /* namespace ipo */