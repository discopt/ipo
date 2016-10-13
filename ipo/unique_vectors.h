#ifndef IPO_UNIQUE_VECTORS_H_
#define IPO_UNIQUE_VECTORS_H_

#include <set>

#include "common.h"
#include "vectors.h"

namespace ipo {

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

