#include "unique_vectors.h"

#include <vector>
#include <random>

namespace ipo {

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
