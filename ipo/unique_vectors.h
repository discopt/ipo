#ifndef IPO_UNIQUE_VECTORS_H_
#define IPO_UNIQUE_VECTORS_H_

#include <map>
#include <set>

#include "common.h"
#include "vectors.h"

namespace ipo {

  template <typename V>
  class VectorMap
  {
  public:
    typedef Vector Key;
    typedef V Value;

  protected:
    struct KeyData
    {
      double hash;
      Vector vector;

      KeyData(double hsh, Vector& vec)
        : hash(hsh), vector(vec)
      {

      }
      
      KeyData(const KeyData& other)
        : hash(other.hash), vector(const_cast<Vector&>(other.vector))
      {

      }

      ~KeyData()
      {

      }

      bool operator<(const KeyData& other) const
      {
        if (hash < other.hash)
          return true;
        else if (hash > other.hash)
          return false;
        return vector < other.vector;
      }
    };

  public:
    VectorMap(std::size_t ambientDimension)
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

    ~VectorMap()
    {

    }

    inline std::size_t size() const
    {
      return _data.size();
    }

    inline bool empty() const
    {
      return _data.empty();
    }

    void clear()
    {
      _data.clear();
    }

    bool insert(Vector& vector, const Value& value)
    {
      KeyData keyData(0.0, vector);
      for (std::size_t p = 0; p < vector.size(); ++p)
        keyData.hash += _hashVector[vector.index(p)] * vector.approximation(p);

      std::pair<typename std::map<KeyData, Value>::iterator, bool> inserted = _data.insert(std::make_pair(keyData, value));
      if (inserted.second)
        return true;
      else
      {
        vector = const_cast<Vector&>(inserted.first->first.vector);
        return false;
      }
    }

    Value& operator[](const Vector& vector)
    {
      Vector copy = vector;
      KeyData keyData(0.0, copy);
      for (std::size_t p = 0; p < vector.size(); ++p)
        keyData.hash += _hashVector[vector.index(p)] * vector.approximation(p);

      return _data.at(keyData);
    }

    const Value& operator[](const Vector& vector) const
    {
      KeyData keyData(0.0, vector);
      for (std::size_t p = 0; p < vector.size(); ++p)
        keyData.hash += _hashVector[vector.index(p)] * vector.approximation(p);

      return _data.at(keyData);
    }

  private:
    std::map<KeyData, Value> _data;
    soplex::DVectorReal _hashVector;
  };

  class UniqueVectors
  {
  protected:
    struct Data
    {
      double hash;
      Vector vector;

      Data(double hash, Vector& vector);
      Data(const Data& other);
      ~Data();

      bool operator<(const Data& other) const;
    };

  public:
    struct Iterator
    {
      Iterator(std::set<Data>::iterator iter)
        : _iter(iter)
      {

      }

      inline Iterator operator++()
      {
        ++_iter;
        return Iterator(_iter);
      }

      inline Vector& operator*()
      {
        return const_cast<Vector&>(_iter->vector);
      }

      inline bool operator==(const Iterator& other) const
      {
        return _iter == other._iter;
      }

      inline bool operator!=(const Iterator& other) const
      {
        return ! (*this == other);
      }

    private:
      std::set<Data>::iterator _iter;
    };

    UniqueVectors(std::size_t ambientDimension);

    ~UniqueVectors();

    inline std::size_t size() const
    {
      return _data.size();
    }

    void clear();

    bool insert(Vector& vector);

    inline Iterator begin()
    {
      return Iterator(_data.begin());
    }

    inline Iterator end()
    {
      return Iterator(_data.end());
    }

  private:
    std::set<Data> _data;
    soplex::DVectorReal _hashVector;
  };

} /* namespace ipo */

#endif /* IPO_UNIQUE_VECTORS_H_ */

