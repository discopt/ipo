#include "unique_vectors.h"

#include <vector>
#include <random>

namespace ipo {

  UniqueVectors::Data::Data(double hsh, Vector& vec)
    : hash(hsh), vector(vec)
  {

  }

  UniqueVectors::Data::Data(const UniqueVectors::Data& other)
    : hash(other.hash), vector(const_cast<Vector&>(other.vector))
  {

  }


  UniqueVectors::Data::~Data()
  {

  }

  bool UniqueVectors::Data::operator<(const UniqueVectors::Data& other) const
  {
    if (hash < other.hash)
      return true;
    else if (hash > other.hash)
      return false;
    return vector < other.vector;
  }

  UniqueVectors::UniqueVectors(std::size_t ambientDimension)
    : _hashVector(ambientDimension)
  {
    std::default_random_engine generator(0);
    std::normal_distribution<double> distribution;
    soplex::DVectorReal randomVector(ambientDimension);
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

  UniqueVectors::~UniqueVectors()
  {

  }

  void UniqueVectors::clear()
  {
    _data.clear();
  }

  bool UniqueVectors::insert(Vector& vector)
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
      vector = const_cast<Vector&>(inserted.first->vector);
      return false;
    }
  }

} /* namespace ipo */
